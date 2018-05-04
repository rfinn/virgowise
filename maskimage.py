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

%run ~/github/virgowise/maskimage.py

sepath = '/Users/rfinn/github/Virgo/sextractor/'
nsapath = '/Users/rfinn/research/NSA/'
im = 'NSA-99882-unwise-w3-img-m.fits'
m = maskimage(im,sepath)
m.run_sextractor()
m.remove_central_object()
m.edit_mask()


See mask_all for how to run on multiple files

mask_all('NSA-*-unwise-w3-img-m.fits')

REQUIRED PROGRAMS:
- sextractor

REQUIRED FILES:
- sepath = path to sextractor files
- nsa fits catalog


NOTES:
"""

from astropy.io import fits
from astropy.wcs import WCS
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

    def remove_central_object(self,xobj=None,yobj=None):
        '''
        by default, this will remove the central object

        if object is  not centered (or to be sure), send in x,y of object
        '''
        
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

        if xobj == None:
            dist=np.sqrt((ycenter-ysex)**2+(xcenter-xsex)**2)
        else:
            dist=np.sqrt((yobj-ysex)**2+(xobj-xsex)**2)
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
        t.close()

        
    def show_image_mask(self,ra=None,dec=None):
        a = fits.open(self.image)
        wcs = WCS(a[0].header)
        b = fits.open(self.mask_image)
        plt.figure(figsize = (8,4))
        plt.subplots_adjust(wspace=.3,right=.95)
        ax = plt.subplot(1,2,1,projection=wcs)
        plt.imshow(a[0].data,origin='lower')
        if ra:
            ax.scatter(ra,dec,transform=ax.get_transform('fk5'),s=200,edgecolor='white',facecolor='none')
        plt.title(self.image)
        ax = plt.subplot(1,2,2,projection=wcs)
        plt.imshow(b[0].data,origin='lower')
        if ra:
            ax.scatter(ra,dec,transform=ax.get_transform('fk5'),s=200,edgecolor='white',facecolor='none')
        plt.title(self.mask_image)
        t = self.mask_image.split('.fits')
        plt.savefig(t[0]+'.png')
        a.close()
        b.close()
    def edit_mask(self,ra=None,dec=None):
        ######################
        ### EDIT MASK
        ######################
        review_mask_flag = True
        while review_mask_flag:
            plt.close()
            self.show_image_mask(ra=ra,dec=dec)

                
            flag=raw_input('edit the mask? \n')
            if flag.find('n') > -1:
                review_mask_flag = 0
            elif flag.find('y') > -1:
                editflag=int(raw_input('enter 0 to subtract an object, 1 to add an object\n'))
                if editflag == 0:
                    objid=float(raw_input('enter object id to remove from mask\n'))
                    self.replace_value(self.mask_image, int(objid))
                #elif editflag == 1:
                #    print 'entering imedit'
                #    a,b=runimedit(mask_image)
            elif flag.find('q') > -1:
                quitflag=1
                return quitflag
            
    def get_xy_image(self,ra,dec):
        h = fits.getheader(self.image)
        w = WCS(h)
        x, y = w.all_world2pix(ra,dec,1)
        return x,y

def mask_all(filestring,pause_flag=True,start_number = 0):
    sepath = '/Users/rfinn/github/virgowise/sextractor/'
    nsapath = '/Users/rfinn/research/NSA/'
    nsa = fits.getdata(nsapath+'nsa_v0_1_2.fits')
    nsadict=dict((a,b) for a,b in zip(nsa.NSAID,np.arange(len(nsa.NSAID)))) #useful for us can easily look up galaxy ID's       
    files = glob.glob(filestring)
    nfiles = len(files)
    for i in range(start_number, len(files)):
        f = files[i]
        print f

        # get nsaid from string
        # assumes this is second entry in line
        # obviously, this is very specific for WISE analysis!

        t = f.split('-')
        nsaid = int(t[1])
        ra = nsa.RA[nsadict[nsaid]]
        dec = nsa.DEC[nsadict[nsaid]]

        m = maskimage(f,sepath)
        x,y = m.get_xy_image(ra,dec)
        #

        m.run_sextractor()
        

        #x,y = get_xy_from_image(image,ra,dec)

        
        m.remove_central_object(xobj=x,yobj=y)
        m.edit_mask(ra=ra,dec=dec)

        print 'file ',i,' of ',nfiles
        if pause_flag:
            t = raw_input('hit any key to continue to next galaxy \n \t enter q to quit \n \t enter C to continue without pausing \n')
            if t.find('q') > -1:
                break
            elif t.find('C') > -1:
                pause_flag = False
        plt.close()
        i += 1
