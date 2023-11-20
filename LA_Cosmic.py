# Copyright 2012 Bernd Husemann
#
#
#This file is part of PyCosmic.
#
#PyCosmic is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License  as published by
#the Free Software Foundation, either version 3 of the License, or
#any later version.
#
#PyCosmic is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with PyCosmic.  If not, see <http://www.gnu.org/licenses/>.


import sys, numpy
import astropy.io.fits as pyfits
from types import *
from LA_image import *

__version__ = "0.2"

def detCos(image,  out_clean,
           sigma_det=5, rlim=1.2, iteration=5, fwhm_gauss=2.0, replace_box=[5,5],
           replace_error=1e10, increase_radius=0, gain=1.0, rdnoise=10, verbose=True):
    """
           Detects and removes cosmics from astronomical images based on Laplacian edge
           detection scheme combined with a PSF convolution approach (Husemann  et al. in prep.).

           IMPORTANT:
           The image and the readout noise are assumed to be in units of electrons.
           The image also needs to be BIAS subtracted! The gain can be entered to convert the image from ADUs to electros, when this is down already set gain=1.0 as the default.

            Parameters
            --------------
            image: string
                    Name of the FITS file for which the comsics should be detected
	    out_mask: string
                    Name of the  FITS file with the bad pixel mask
            out_clean: string
                    Name of the  FITS file with the cleaned image
            rdnoise: float or string of header keyword
                    Value or FITS header keyword for the readout noise in electrons
            sigma_det: float, optional  with default: 5.0
                    Detection limit of edge pixel above the noise in (sigma units) to be detected as comiscs
            rlim: float, optional  with default: 1.2
                    Detection threshold between Laplacian edged and Gaussian smoothed image
            iter: integer, optional with default: 5
                    Number of iterations. Should be >1 to fully detect extended cosmics
            fwhm_gauss: float, optional with default: 2.0
                    FWHM of the Gaussian smoothing kernel in x and y direction on the CCD
            replace_box: array of two integers, optional with default: [5,5]
                    median box size in x and y to estimate replacement values from valid pixels
            replace_error: float, optional with default: 1e10
                    Error value for bad pixels in the comupted error image, will be ignored if empty
            increase_radius: integer, optional with default: 0
                    Increase the boundary of each detected cosmic ray pixel by the given number of pixels.
            verbose: bollean, optional  with default: True
                    Show information during the processing on the command line (0 - no, 1 - yes)


            References
            --------------
            B. Husemann et al. 2012  "", A&A, ??, ???

    """
    print('LACosmic started for image ', image)

    # convert all parameters to proper type
    iterations = iteration
    sigma = fwhm_gauss/2.354
    box_x = replace_box[0]
    box_y = replace_box[1]

    # load image from FITS file
    img = loadImage(image)
    try:
        gain=img.getHdrValue('gain')
    except KeyError:
        pass
    gain = float(gain)

    if gain!=1.0 and verbose==True:
      print ('Convert image from ADUs to electrons using a gain factor of %f' %(gain))

    img = img*gain
    #img.writeFitsData('test.fits')

    # create empty mask if no mask is present in original image
    if img._mask!=None:
        mask_orig=img.getMask()
    else:
        mask_orig=numpy.zeros(img.getDim(), dtype=bool)

    # create a new Image instance to store the initial data array
    img_original = Image(data=img.getData(), header=img.getHeader(), error = img.getError(),  mask=mask_orig)
    img.setData(mask=numpy.zeros(img.getDim(), dtype=bool))
    img.removeError()

    # estimate Poisson noise after roughly cleaning cosmics using a median filter
    try:
        rdnoise=float(img.getHdrValue('rdnoise'))
    except KeyError:
        rdnoise=float(rdnoise)
    if verbose==True:
        print ('A value of %f is used for the electron read-out noise.'%(rdnoise))


    # create empty mask
    select = numpy.zeros(img.getDim(),dtype=bool)

    # define Laplacian convolution kernal
    LA_kernel=numpy.array([[0,-1,0,],[-1,4,-1],[0,-1,0]])/4.0
    out=img

    # start iteration
    for i in range(iterations):
        if verbose:
            print ('Start iteration %i'%(i+1))
        # follow the LACosmic scheme to select pixel
        noise =out.medianImg((5, 5))
        select_neg2 = noise.getData()<=0
        noise.setData(data=0, select=select_neg2)
        noise=(noise+rdnoise**2).sqrt()
        result = []

        sub = out.subsampleImg(2) # subsample image
        conv= sub.convolveImg(LA_kernel) # convolve subsampled image with kernel
        select_neg = conv<0
        conv.setData(data=0, select=select_neg)  # replace all negative values with 0
        Lap = conv.rebin(2, 2) # rebin the data to original resolution
        S = Lap/(noise*2) # normalize Laplacian image by the noise
        S_prime = S-S.medianImg((5, 5)) # cleaning of the normalized Laplacian image
        fine=out.convolveGaussImg(sigma, sigma, mask=True) # convolve image with a 2D Gaussian

        fine_norm = out/fine
        select_neg = fine_norm<0
        fine_norm.setData(data=0, select=select_neg)
        sub_norm = fine_norm.subsampleImg(2) # subsample image
        Lap2 = (sub_norm).convolveImg(LA_kernel)
        Lap2 = Lap2.rebin(2, 2) # rebin the data to original resolution

        select = numpy.logical_or(numpy.logical_and((Lap2)>rlim, S_prime>sigma_det),  select)

        # print information on the screen if demanded
        if verbose:
            dim = img_original.getDim()
            det_pix = numpy.sum(select)
            print ('Total number of detected cosmics: %i out of %i pixels'%(numpy.sum(select), dim[0]*dim[1]))

        if i==iterations-1:
            img_original.setData(mask=True, select=select) # set the new mask
            if increase_radius>0:
                mask_img = Image(data=img_original._mask)
                mask_new=mask_img.convolveImg(kernel=numpy.ones((2*increase_radius+1, 2*increase_radius+1)))
                img_original.setData(mask=mask_new._data)
            out=img_original.replaceMaskMedian(box_x, box_y, replace_error=replace_error) # replace possible corrput pixel with zeros for final output
        else:
            out.setData(mask=True, select=select)# set the new mask
            out = out.replaceMaskMedian(box_x, box_y, replace_error=None)  # replace possible corrput pixel with zeros
    if verbose==True:
      print ('Cleaned image is stored in file: %s'%(out_clean))
      #print ('Cosmics mask is stored in file: %s'%(out_mask))
    out = out/gain
    out.writeFitsData(out_clean, extension_data=0, extension_mask=0, extension_error=0)
