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
import pyds9

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
        sexnumber=sexout['col1']
        xsex=sexout['col2']
        ysex=sexout['col3']
        xcenter=1.*(xmax-xmin)/2.
        ycenter=1.*(ymax-ymin)/2.
        dist=sqrt((ycenter-ysex)**2+(xcenter-xsex)**2)

        #   find object ID
        objIndex=where(dist == min(dist))
        objNumber=sexnumber[objIndex]
        objNumber=objNumber[0] # not sure why above line returns a list
        print 'object number = ',objNumber
        #   use iraf imexam to replace object ID values with zero
        iraf.imreplace(mask_image,value=0,lower=objNumber-.5,upper=objNumber+.5)

    def edit_mask(self):
        ######################
        ### EDIT MASK
        ######################
        d=pyds9.DS9()
        while review_mask_flag:
            d.set('frame delete all')
            s='file new '+self.image
            try:
                d.set(s)
                d.set('zoom to fit')
            except:
                print "couldn't access: ",s
            s='file new '+sex_image
            d.set(s)
            d.set('zoom to fit')
            s='file new '+mask_image
            d.set(s)
            d.set('zoom to fit')
            flag=raw_input('edit the mask? \n')
            if flag.find('n') > -1:
                review_mask_flag = 0
            elif flag.find('y') > -1:
                editflag=int(raw_input('enter 0 to subtract an object, 1 to add an object'))
                if editflag == 0:
                    objid=float(raw_input('enter object id to remove from mask'))
                    iraf.imreplace(mask_image,value=0,lower=objid-.5,upper=objid+.5)
                elif editflag == 1:
                    print 'entering imedit'
                    a,b=runimedit(mask_image)
            elif flag.find('q') > -1:
                quitflag=1
                return quitflag

