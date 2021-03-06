from __future__ import division

import time
import os
import logging
import threading

from chimera.core.chimeraobject import ChimeraObject
from chimera.core.lock import lock
from chimera.core.event import event
from chimera.core.exceptions import ChimeraException, ClassLoaderException
from chimera.core.constants import SYSTEM_CONFIG_DIRECTORY

from chimera.interfaces.autoguider import Autoguider as IAutoguider
from chimera.interfaces.autoguider import GuiderStatus
from chimera.interfaces.autoguider import StarNotFoundException

from chimera.controllers.imageserver.imagerequest import ImageRequest
from chimera.controllers.imageserver.util         import getImageServer

from chimera.util.image import Image
from chimera.util.coord import Coord
from chimera.util.output import red, green

from astropy.io import fits
from astropysics.utils.alg import centroid

import numpy as np
import yaml

plot = True
try:
    import pylab as py
    py.ion()
except (ImportError, RuntimeError, ClassLoaderException):
    plot = False

class AutoGuider(ChimeraObject,IAutoguider):
    """
    Auto Guider
    ===========

    This controller is reponsible for proving auto-guiding capabilities
    to a telescope. It works by acquiring a star in a CCD and trying to
    guarantee that it remain in its initial position, offsetting the
    telescope accordingly. The operation is divided in two steps.

    1) Finding a guide star
    2) Guiding

    Finding a guide star consists of:

    1) Take an acquisition image
    2) Check if a predetermined guide-star is in the FoV or try to
        find a suitable bright star for guiding. Raise an exception if
        no guide star is found.
    3) Calculate the position of the star in the CCD (using centroid)

    Guiding consists of

    1) Take an image
    2) Calculate the centroid of the guide-star
    3) Offset the telescope
    4) Repeat until stopped or an error occur.

    """

    def __init__(self):
        ChimeraObject.__init__(self)

        self.imageRequest = None
        self.filter = None

        self.currentRun = None

        self._debugging = False
        self._debug_images = []
        self._debug_image = 0

        self._log_handler = None

        self._state = GuiderStatus.OFF

        self.abort = threading.Event()
        self.abort.clear()

    def getTel(self):
        return self.getManager().getProxy(self["telescope"])

    def getCam(self):
        return self.getManager().getProxy(self["camera"])

    def getFilter(self):
        return self.getManager().getProxy(self["filterwheel"])

    def getFocuser(self):
        return self.getManager().getProxy(self["focuser"])

    def getAutoFocus(self):
        return self.getManager().getProxy(self["autofocus"])

    def getPointVerify(self):
        return self.getManager().getProxy(self["point_verify"])

    def getSite(self):
        return self.getManager().getProxy(self["site"])

    def _getID(self):
        return "autoguider-%s" % time.strftime("%Y%m%d-%H%M%S")

    def _openLogger(self):

        if self._log_handler:
            self._closeLogger()

        self._log_handler = logging.FileHandler(os.path.join(SYSTEM_CONFIG_DIRECTORY, self.currentRun,
                                                             "autoguider.log"))
        self._log_handler.setFormatter(logging.Formatter(fmt="%(message)s"))
        self._log_handler.setLevel(logging.DEBUG)
        self.log.addHandler(self._log_handler)

    def _closeLogger(self):
        if self._log_handler:
            self.log.removeHandler(self._log_handler)
            self._log_handler.close()

    @lock
    def guide(self,position=None,exptime=None,interval=None,filter=None,binning=None,
              windowCCD=True,box=None,debug=False):
        """

        :param position: <type Position>. Ra/Dec of the guide star. Leave empty for automatic search.
        :param exptime: Exposure time of the acquisition images. This will set the frequency of updates as well.
        :param interval: Time to wait between guiding sequence. This will lower the frequency of updates.
        :param filter: Filter for the exposures.
        :param binning: Binning of the guide CCD
        :param windowCCD: Flag to window the acquisition images to the size of the guiding box. This may speed the
                          read-out of the CCD. But may not be a good idea in some cases.
        :param box: The size of the guiding box in pixels. If None will try to guess from the FWHM of the star
        :param debug: Run in debug mode?
        :return:
        """

        self._debugging = debug

        self.currentRun = self._getID()

        if not os.path.exists(os.path.join(SYSTEM_CONFIG_DIRECTORY, self.currentRun)):
            os.mkdir(os.path.join(SYSTEM_CONFIG_DIRECTORY, self.currentRun))

        self._openLogger()

        if debug:
            debug_file = open(os.path.join(debug, "autofocus.debug"), "r")
            debug_data = yaml.load(debug_file.read())

            start = debug_data["start"]
            end = debug_data["end"]
            step = debug_data["step"]

            debug_file.close()

        self.log.debug("Taking acquisition image.")

        self.imageRequest = ImageRequest()
        self.imageRequest["exptime"] = exptime or 1.
        self.imageRequest["frames"] = 1
        self.imageRequest["shutter"] = "OPEN"

        if filter:
            self.filter = filter
            self.log.debug("Using filter %s." % self.filter)
        else:
            self.filter = False
            self.log.debug("Using current filter.")

        if binning:
            self.imageRequest["binning"] = binning

        # No window in acquisition image

        if position:
            star_found = self._findGuideStar(self._takeImageAndResolveStars(),position)
        else:
            star_found = self._findBestStarToGuide(self._takeImageAndResolveStars())

        if not star_found:
            tries = 0

            while not star_found and tries < self["max_tries"]:
                star_found = self._findBestStarToGuide(self._takeImageAndResolveStars())
                tries += 1

            if not star_found:
                self._state = GuiderStatus.ERROR
                raise StarNotFoundException("Couldn't find a suitable star to guide on."
                                            "Giving up after %d tries." % tries)

        # Start guiding
        try:
            if windowCCD:
                # retake first image around guide star
                if not box:
                    box = np.int(np.ceil(star_found['FWHM_IMAGE']*3))
                self.imageRequest["window"] = "%i:%i,%i:%i"%(star_found['XWIN_IMAGE']-box/2,
                                                               star_found['XWIN_IMAGE']+box/2,
                                                               star_found['YWIN_IMAGE']-box/2,
                                                               star_found['YWIN_IMAGE']+box/2)
                # self.log.debug('Cutting CCD in subframe: %s'%(self.imageRequest['window']))
                # star_found['XWIN_IMAGE'],star_found['YWIN_IMAGE'] = self.getCentroid()
                #star_found = self._findBestStarToGuide(self._takeImageAndResolveStars())
                star_found['XWIN_IMAGE'] = box/2
                star_found['YWIN_IMAGE'] = box/2

            self.abort.clear()
            self.guideStart(star_found)

            msg=''
            self._state = GuiderStatus.GUIDING
            self.index = 1

            def process():
                nlost = 0

                while True:
                    offset = self.getOffset(star_found)
                    self.log.debug('Offset %f x %f'%(offset['N'],offset['E']))
                    if offset:
                        self.log.debug('Applying offset %s x %s'%(offset['N'],offset['E']))
                        self.applyOffset(offset)
                        self.offsetComplete(offset)
                    else:
                        nlost+=1

                    if self.abort.isSet():
                        self._state = GuiderStatus.ABORTED
                        msg = 'Sequence aborted.'
                        break

                    if nlost > self['nlost']:
                        self._state = GuiderStatus.ERROR
                        msg = 'Lost guide star. Stopping after %i tries.'%self['nlost']
                        break

                    if interval:
                        time.sleep(interval)

                self.guideStop(self._state,msg)

            t = threading.Thread(target=process)
            t.setDaemon(False)
            t.start()

        except Exception, e:
            self._state = GuiderStatus.ERROR
            self.guideStop(self._state,e)
            raise

        finally:
            # reset debug counter
            self._debug_image = 0

    def _takeImageAndResolveStars(self):

        frame = self._takeImage()
        stars = self._findStars(frame)

        return stars

    def _takeImage(self):

        if self._debugging:
            pass
            try:
                frame = self._debug_images[self._debug_image]
                self._debug_image += 1

                img = Image.fromFile(frame)
                srv = getImageServer(self.getManager())
                return srv.register(img)
            except IndexError:
                raise ChimeraException("Cannot find debug images")

        self.imageRequest["filename"] = os.path.join(SYSTEM_CONFIG_DIRECTORY, self.currentRun, "autoguider.fits")

        cam = self.getCam()

        if self.filter:
            filter = self.getFilter()
            filter.setFilter(self.filter)

        frame = cam.expose(self.imageRequest)

        if frame:
            return frame[0]
        else:
            raise Exception("Error taking image.")

    def _findStars(self, frame):

        config = {}
        config['PIXEL_SCALE'] = 0  # use WCS info
        config['BACK_TYPE'] = "AUTO"

        config['SATUR_LEVEL'] = self.getCam()["ccd_saturation_level"]

        # improve speed with higher threshold
        config['DETECT_THRESH'] = 3.0

        # no output, please
        config['VERBOSE_TYPE'] = "QUIET"

        # our "star" dict entry will contain all this members
        config['PARAMETERS_LIST'] = ["NUMBER",
                                     "XWIN_IMAGE", "YWIN_IMAGE",
                                     "FLUX_BEST", "FWHM_IMAGE",
                                     "FLAGS"]

        catalogName = os.path.splitext(frame.filename())[0] + ".catalog"
        configName = os.path.splitext(frame.filename())[0] + ".config"
        return frame.extract(config, saveCatalog=catalogName, saveConfig=configName)

    def _findGuideStar(self,catalog,position):
        return

    def _findBestStarToGuide(self,catalog):
        # simple plan: brighter star
        # FIXME: avoid "border" stars
        return self._findBrighterStar(catalog)

    def _findBrighterStar(self, catalog):

        fluxes = [star for star in catalog if star["FLAGS"] == 0]

        if not fluxes:  # empty catalog
            return False

        return max(fluxes, key=lambda star: star["FLUX_BEST"])


    def getOffset(self,position):

        ret = {'X':0.,'Y':0.,
               'N':0.,'E':0.,'Status':self.state()}

        try:
            frame = self._takeImage()
            fname = '/tmp/autoguider.fits'
            if os.path.exists(fname):
                os.remove(fname)
            frame.save(fname)
            img = fits.getdata(fname)
            # Extract some backgroud
            img -= (np.mean(img)*0.9)
            img[img < 0] = 0.

            pY,pX = centroid(img)
            ret['X'] = pX-position["XWIN_IMAGE"]
            ret['Y'] = pY-position["YWIN_IMAGE"]
            centerPos = frame.worldAt([position["XWIN_IMAGE"],
                                    position["YWIN_IMAGE"]])
            currPos = frame.worldAt([pX,
                                     pY])
            # offset = centerPos.angsep(currPos)
            ret['E'] = Coord.fromAS((centerPos.ra.AS - currPos.ra.AS) * np.cos(currPos.dec.R))
            ret['N'] = Coord.fromAS(centerPos.dec.AS - currPos.dec.AS)
        except:
            if self.abort.isSet():
                ret['Status'] = GuiderStatus.ABORTED
                return ret
            else:
                raise

        self.plot(frame,position,ret)

        return ret

    def plot(self,frame,position,offset):

        global plot

        if plot:
            fname = '/tmp/autoguider.fits'
            if os.path.exists(fname):
                os.remove(fname)
            frame.save(fname)
            img = fits.getdata(fname)
            hdr = fits.getheader(fname)
            py.figure(1)
            py.cla()
            py.imshow(img,origin='lower',interpolation='nearest',vmax=np.mean(img)*1.1,cmap=py.cm.gray)
            py.plot(position['XWIN_IMAGE'],position['YWIN_IMAGE'],'b+')
            py.plot(position['XWIN_IMAGE']+offset['X'],position['YWIN_IMAGE']+offset['Y'],'rx')
            py.plot(position['XWIN_IMAGE']-hdr['offsetY']+1,position['YWIN_IMAGE']-hdr['offsetX'],'r.')
            err = np.sqrt(((-hdr['offsetY']+1 - offset['X']) ** 2) + ((-hdr['offsetX'] - offset['Y']) ** 2))

            py.text(-3, -3, 'Offset: (%.2f/%.2f) px, Calculated: (%.2f/%.2f) px, Error: %3.2f' % (-hdr['offsetY']+1,
                                                                                                 -hdr['offsetX'],
                                                                                                 offset['X'],
                                                                                                 offset['Y'],
                                                                                                 err))
            py.savefig(os.path.join(SYSTEM_CONFIG_DIRECTORY, 'autoguider_%05i.png')%self.index)
            py.show()
            self.index += 1

    def applyOffset(self, offset):
        tel = self.getTel()
        tel.moveOffset(offset['N'],offset['E'])
        # return

    def stop(self):

        self.abort.set() # do this first!
        self.getCam().abortExposure(readout=False)
        self._state = GuiderStatus.OFF

    def state(self):
        return self._state