from __future__ import division

"""
This scripts creates a video file simulating an autoguiding video. It takes observed images as input.
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
log.setLevel(logging.INFO)

def main(argv):

    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option('-f', '--filename',
                      help='Input list of fits files.'
                      , type='string')
    parser.add_option('-v', '--verbose',
                      help='Run in verbose mode.', action='store_true',
                      default=False)

    opt, args = parser.parse_args(argv)
    if opt.verbose:
        log.setLevel(logging.DEBUG)

    # 1 - Read list of images
    imglist = np.loadtxt(opt.filename,dtype='S')

    # 2 - First image is reference
    refimg = fits.open(imglist[0])

    # 2 - Telescope characteristics
    telescope_platescale = 133  # arcsec/mm

    # 2 - Autoguider center
    # For an SBIG ST7XME the parameters of the offset between the CCD and the guider are:
    guider_offset = {'x': 0, 'y': 0}  # in mm
    guider_size = {'x': refimg[0].header['CCDPXSZX']*1e-3*refimg[0].header['CCD_DIMX'],
                   'y': refimg[0].header['CCDPXSZY']*1e-3*refimg[0].header['CCD_DIMY'],
                   'x_px': refimg[0].header['CCD_DIMX'],
                   'y_px': refimg[0].header['CCD_DIMY']}  # in mm / px and pixels
    guider_pixelsize = {'x': refimg[0].header['CCDPXSZX']*1e-3,
                        'y': refimg[0].header['CCDPXSZY']*1e-3}  # in mm
    guider_pixelsize.update({'x_ang': telescope_platescale * guider_pixelsize['x'],
                             'y_ang': telescope_platescale * guider_pixelsize['y']})  # in arcsec / px

    # For plotting purposes
    ccd_size = {'x': guider_size['x'], 'y': guider_size['y']}  # in mm
    ccd_pixelsize = {'x': guider_pixelsize['x'], 'y': guider_pixelsize['y']}  # in mm

    # Get a DSS Image...
    img2 = np.array(refimg[0].data,copy=True)
    img = np.array(refimg[0].data,copy=True)
    # Image sizes in pixels
    imsize = int(round(guider_size['x'] / guider_pixelsize['x'])), int(round(guider_size['y'] / guider_pixelsize['y']))
    # img = imresize(img, imsize)
    # Image sizes in arcmins
    aux_width = guider_size['x'] * telescope_platescale / 60
    aux_height = guider_size['y'] * telescope_platescale / 60

    # TODO: DELME
    img_scale = 1.2  # in arcsec / pix
    imsize = img.shape[0], img.shape[1]
    imsize_resampled = int(round(imsize[0] * img_scale / guider_pixelsize['x_ang'])), \
                       int(round(imsize[1] * img_scale / guider_pixelsize['y_ang']))
    pix = imresize(img, imsize_resampled)
    # TODO: DELME END

    px_border = 20  # Number of pixels to throw away due to border-effect
    max_offset = 10  # Maximum offset between images, in pixels

    img_corner = [int(round(pix.shape[0] / 2 - guider_size['y_px'] / 2)),
                  int(round(pix.shape[1] / 2 - guider_size['x_px'] / 2))]

    mask = img < 1.1 * np.median(img)  # Background / weak stars subtraction
    #print img_new
    img[mask] = 0.
    segmap = label(img > 0)
    print segmap[segmap>0]
    # fig2 = plt.figure(2)
    # plt.clf()
    # plt.imshow(segmap)
    plt.figure(1)
    cent = []
    for i_object in range(1, np.max(segmap) + 1):
        mm = segmap == i_object
        if len(mm[mm]) > 1:
            cent.append(centroid(img * mm))
        # plt.plot(cent_new[-1][1], cent_new[-1][0], '.', lw=2, color='black')
        # if i > 0:
        # plt.plot(cent[-1][1], cent[-1][0], '.', lw=2, color='black')
    cent = np.array(cent)

    print img2.min(),img2.max(),img2.mean()

    plt.imshow(img2,cmap=plt.cm.gray,interpolation='nearest',origin='lower',
               vmin=img2.min(),
               vmax=img2.mean()*1.1)
    for i in range(len(cent)):
        plt.plot(cent[i][1],cent[i][0],marker='o',markerfacecolor='none',mec='r')
    plt.xlim(0, img.shape[1])
    plt.ylim(0, img.shape[0])

    # plt.show()
    # return 0

    for i in range(len(imglist)):
        #offset = np.random.randint(-10, 10, size=2)
        #img_corner = np.sum([img_corner, offset], axis=0)  # Add an offset to the corner pixel of the image
        #img_cut_new = np.copy(
        #    pix[img_corner[0]:img_corner[0] + guider_size['y_px'], img_corner[1]:img_corner[1] + guider_size['x_px']])
        img_cut_new = fits.getdata(imglist[i])
        img_cut_new = img_cut_new[px_border:-px_border, px_border:-px_border]
        fig1 = plt.figure(1)
        plt.clf()
        plt.imshow(img_cut_new, ,interpolation='nearest',origin='lower',
               vmin=img2.min(),
               vmax=img2.mean()*1.1)
        plt.xlim(0, img_cut_new.shape[1])
        plt.ylim(img_cut_new.shape[0], 0)
        mask = img_cut_new < 1.5 * np.median(img_cut_new)  # Background / weak stars subtraction
        img_cut_new[mask] = 0
        segmap = label(img_cut_new > 0)
        fig2 = plt.figure(2)
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
            #print 'Applied offset in pixels:', offset
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

if __name__ == '__main__':
    import sys
    main(sys.argv)