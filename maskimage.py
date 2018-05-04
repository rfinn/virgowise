#!/usr/bin/env python

"""
GOAL:
- to mask extraneous sources in WISE image

PROCEDURE:
- run sextractor
- display image and segmentation image
- allow user to edit mask
- save mask

USEAGE:

from within ipython on laptop

%run ~/github/virgowise/maskimage.py --virgo '/Users/rfinn/github/Virgo/' 

REQUIRED PROGRAMS:
- sextractor

REQUIRED FILES:
- nsa fits catalog


NOTES:
"""

from astropy.io import fits
#import pyds9
import os
import numpy as np
from matplotlib import pyplot as plt
import glob

class maskimage():
    def __init__(self,imagename,sepath):
        self.image = imagename
        t = self.image.split('.fits')
        self.mask_image=t[0]+'-mask.fits'
        self.sepath = sepath
    def run_sextractor(self):

        ##############################################
        # RUN SEXTRACTOR AND MAKE OBJECT MASK
        ##############################################

        # run sextractor to generate a list of objects in the image and generate 'segmentation image'
        os.system('cp '+self.sepath+'/default.sex.wise .')
        os.system('cp '+self.sepath+'/default.param .')
        os.system('cp '+self.sepath+'/default.conv .')
        os.system('cp '+self.sepath+'/default.nnw .')
        os.system('sex %s -c default.sex.wise '%(self.image))
        # convert segmentation image to object mask by replacing the object ID of target galaxy with zeros
        #   parse sextractor output to get x,y coords of objects

        if os.path.exists(self.mask_image):
            os.remove(self.mask_image)
        os.rename('segmentation.fits',self.mask_image)

    def remove_central_object(self):
        sexout=fits.getdata('test.cat',2)
        sexnumber=sexout['NUMBER']
        xsex=sexout['X_IMAGE']
        ysex=sexout['Y_IMAGE']

        mask_image = fits.getdata(self.mask_image)
        xmax,ymax = mask_image.shape
        xmin=1
        ymin=1
        xcenter=1.*(xmax-xmin)/2.
        ycenter=1.*(ymax-ymin)/2.
        dist=np.sqrt((ycenter-ysex)**2+(xcenter-xsex)**2)

        #   find object ID
        objNumber=sexnumber[dist == min(dist)]
        objNumber=objNumber[0] # not sure why above line returns a list
        print 'object number = ',objNumber
        #   use iraf imexam to replace object ID values with zero
        self.replace_value(self.mask_image,objNumber)

    def replace_value(self,image,value):
        # sets value in image to zeros
        t = fits.open(image)
        t[0].data = t[0].data*(~(t[0].data == value))
        t.writeto(image,overwrite=True)

        
    def show_image_mask(self):
        a = fits.getdata(self.image)
        b = fits.getdata(self.mask_image)
        plt.figure(figsize = (8,4))
        plt.subplot(1,2,1)
        plt.imshow(a)
        plt.title(self.image)
        plt.subplot(1,2,2)
        plt.imshow(b)
        plt.title(self.mask_image)
        
    def edit_mask(self):
        ######################
        ### EDIT MASK
        ######################
        review_mask_flag = True
        while review_mask_flag:
            self.show_image_mask()
            flag=raw_input('edit the mask? \n')
            if flag.find('n') > -1:
                review_mask_flag = 0
            elif flag.find('y') > -1:
                editflag=int(raw_input('enter 0 to subtract an object, 1 to add an object\n'))
                if editflag == 0:
                    objid=float(raw_input('enter object id to remove from mask'))
                    self.replace_value(self.mask_image, int(objid))
                #elif editflag == 1:
                #    print 'entering imedit'
                #    a,b=runimedit(mask_image)
            elif flag.find('q') > -1:
                quitflag=1
                return quitflag

def mask_all(band):
    sepath = '/Users/rfinn/github/Virgo/sextractor/' 
    files = glob.glob('NSA*unwise-w'+str(band)+'-img-m.fits')
    for f in files:
        m = maskimage(f,sepath)
        m.run_sextractor()
        m.remove_central_object()
        m.edit_mask()
        
