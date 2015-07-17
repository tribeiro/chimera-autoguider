from __future__ import division

"""
This scripts creates a video file simulating an autoguiding video. Image data is from DSS.
"""
import os
import urllib
import logging
import time
import shutil
from astropy.io import fits
from chimera.util.coord import Coord

from scipy.misc import imresize

import numpy as np
import matplotlib.pyplot as plt

from astropysics.utils.alg import centroid

from skimage.measure import label

log = logging.getLogger()
log.setLevel(logging.DEBUG)


def get_dss(ra, dec, width, height, px_size=None):
    url = "http://stdatu.stsci.edu/cgi-bin/dss_search?"
    query_args = {"r":
                      ra.strfcoord('%(h)02d:%(m)02d:%(s)04d'),
                  "d":
                      dec.strfcoord('%(d)02d:%(m)02d:%(s)04d', signed=True),
                  "f": "fits",
                  "e": "j2000",
                  "c": "gz",
                  "fov": "NONE"}

    # use POSS2-Red surbey ( -90 < d < -20 ) if below -25 deg declination, else use POSS1-Red (-30 < d < +90)
    # http://www-gsss.stsci.edu/SkySurveys/Surveys.htm
    if dec.D < -25:
        query_args["v"] = "poss2ukstu_red"
        query_args["h"] = height  # / 59.5  # ~1"/pix (~60 pix/arcmin) is the plate scale of DSS POSS2-Red
        query_args["w"] = width  # / 59.5
        img_scale = 60 / 59.5  # in arcsec / pix
    else:
        query_args["v"] = "poss1_red"
        query_args["h"] = height  # / 35.3 # 1.7"/pix (35.3 pix/arcmin) is the plate scale of DSS POSS1-Red
        query_args["w"] = width  # / 35.3
        img_scale = 60 / 35.3  # in arcsec / pix

    url += urllib.urlencode(query_args)

    log.debug("Attempting URL: " + url)
    try:
        t0 = time.time()
        dssfile = urllib.urlretrieve(url)[0]
        log.debug("download took: %.3f s" % (time.time() - t0))
        fitsfile = dssfile + ".fits.gz"
        shutil.copy(dssfile, fitsfile)
        hdulist = fits.open(fitsfile)
        pix = hdulist[0].data
        hdulist.close()
        os.remove(fitsfile)
    except Exception, e:
        log.warning("General error getting DSS image: " + str(e))

    if px_size:  # TODO: Check-me!
        # Re-scale image to the camera pixelsize
        imsize = pix.shape[0], pix.shape[1]
        imsize_resampled = int(round(imsize[0] * img_scale / px_size['x_ang'])), \
                           int(round(imsize[1] * img_scale / px_size['y_ang']))
        pix = imresize(pix, imsize_resampled)

        # imsize        img_scale
        # new_imsize    guider_scale

    return pix

# 1 - Telescope characteristics
telescope_platescale = 50  # arcsec/mm

# 2 - Autoguider center
# For an SBIG ST7XME the parameters of the offset between the CCD and the guider are:
guider_offset = {'x': 0, 'y': 6.27}  # in mm
guider_size = {'x': 4.86, 'y': 3.66, 'x_px': 657, 'y_px': 495}  # in mm / px and pixels
guider_pixelsize = {'x': 0.0074, 'y': 0.0074}  # in mm
guider_pixelsize.update({'x_ang': telescope_platescale * guider_pixelsize['x'],  # in mm
                         'y_ang': telescope_platescale * guider_pixelsize['y']})  # in arcsec / px

# For plotting purposes
ccd_size = {'x': 6.89, 'y': 4.59}  # in mm
ccd_pixelsize = {'x': 0.009, 'y': 0.009}  # in mm

# Get a DSS Image...
ra = Coord.fromDMS("12 53 37.08")
dec = Coord.fromDMS("-60 21 22.7")
img = get_dss(ra, dec, 60, 60)
# fits.writeto('test.fits', img)
#img = fits.getdata('test.fits')
img = np.array(np.copy(img), dtype=np.float32)
# Image sizes in pixels
imsize = int(round(guider_size['x'] / guider_pixelsize['x'])), int(round(guider_size['y'] / guider_pixelsize['y']))
# img = imresize(img, imsize)
# Image sizes in arcmins
aux_width = guider_size['x'] * telescope_platescale / 60
aux_height = guider_size['y'] * telescope_platescale / 60

# TODO: DELME
img_scale = 60 / 35.3  # in arcsec / pix
imsize = img.shape[0], img.shape[1]
imsize_resampled = int(round(imsize[0] * img_scale / guider_pixelsize['x_ang'])), \
                   int(round(imsize[1] * img_scale / guider_pixelsize['y_ang']))
pix = imresize(img, imsize_resampled)
# TODO: DELME END

px_border = 20  # Number of pixels to throw away due to border-effect
max_offset = 10  # Maximum offset between images, in pixels

ximg_corner = [int(round(pix.shape[0] / 2 - guider_size['y_px'] / 2)),
              int(round(pix.shape[1] / 2 - guider_size['x_px'] / 2))]
for i in range(300):
    offset = np.random.randint(-10, 10, size=2)
    img_corner = np.sum([ximg_corner, offset], axis=0)  # Add an offset to the corner pixel of the image
    img_cut_new = np.copy(
        pix[img_corner[0]:img_corner[0] + guider_size['y_px'], img_corner[1]:img_corner[1] + guider_size['x_px']])
    img_cut_new = img_cut_new[px_border:-px_border, px_border:-px_border]
    fig1 = plt.figure(1);
    plt.clf()
    plt.imshow(np.log10(img_cut_new), cmap=plt.cm.gray)
    plt.xlim(0, img_cut_new.shape[1])
    plt.ylim(img_cut_new.shape[0], 0)
    mask = img_cut_new < 1.5 * np.median(img_cut_new)  # Background / weak stars subtraction
    img_cut_new[mask] = 0
    segmap = label(img_cut_new > 0)
    fig2 = plt.figure(2);
    plt.clf()
    plt.imshow(segmap)
    plt.figure(1)
    cent_new = []
    for i_object in range(1, np.max(segmap) + 1):
        cent_new.append(centroid(img_cut_new * (segmap == i_object)))
        # plt.plot(cent_new[-1][1], cent_new[-1][0], '.', lw=2, color='black')
        # if i > 0:
        # plt.plot(cent[-1][1], cent[-1][0], '.', lw=2, color='black')
    cent_new = np.array(cent_new)

    if i > 0:
        print 'Applied offset in pixels:', offset
        s = min(len(cent_new), len(cent))
        if s > 0:
            aux_diff = cent[:s] - cent_new[:s]
            # Throw away the differences above the threshold
            # mask_cent_diff = np.sqrt(np.sum(aux_diff **2, axis=1)) <= max_offset
            mask_cent_diff = np.bitwise_and(max_offset <= cent[:, 1][:s], cent[:, 1][:s] <= guider_size['x_px'] - max_offset)
            mask_cent_diff *= np.bitwise_and(max_offset <= cent[:, 0][:s], cent[:, 0][:s] <= guider_size['y_px'] - max_offset)
            mask_cent_diff *= np.bitwise_and(max_offset <= cent_new[:, 1][:s], cent_new[:, 1][:s] <= guider_size['x_px'] - max_offset)
            mask_cent_diff *= np.bitwise_and(max_offset <= cent_new[:, 0][:s], cent_new[:, 0][:s] <= guider_size['y_px'] - max_offset)
            cent_diff = np.median(aux_diff[mask_cent_diff], axis=0)
            print mask_cent_diff
            print 'Centroid differences:', cent_diff
            err = np.sqrt(np.sum((cent_diff - offset) ** 2))
            print 'Error:', err
            color = 'white' if err == 0 else 'red'
            plt.text(5, 15, 'Offset: %s px, Calculated: %s px, Error: %3.2f' % (offset, cent_diff, err),
                     color=color)
            if color == 'red': print 'Attn: %i' % i
            for i_object in np.arange(s)[mask_cent_diff]:
                plt.plot(cent[i_object][1], cent[i_object][0], '.', lw=10, color='green')
                plt.plot(cent_new[i_object][1], cent_new[i_object][0], '.', lw=10, color='blue')
        else:
            plt.text(5, 15, 'Could not calculate offset!', color='red')

    cent = np.copy(cent_new)

    fig1.savefig('test/img_%04d.png' % i)

    # raw_input('Next...')
