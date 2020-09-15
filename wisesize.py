
"""
GOAL:
- to retrieve unWISE images of galaxies given their ra, dec, image size
- run galfit

PROCEDURE:
- download unwise image for particular band as a tar file


USEAGE:

from within ipython on laptop


REQUIRED FILES:


NOTES:
* updating Aug 2020 to remove randomization of initial starting parameters
* just want to take in an RA and DEC (with option to provide and fix the BA and PA)
  - download WISE image
  - run galfit on W3 (and W4?)
  - return results


* down the road, need to handle galaxies that have multiple unwise images
  - combine the images
  - make a cutout from the combined image
"""


import os
#import pyds9
import numpy as np
import argparse
from astropy.io import fits
from astropy.visualization import simple_norm
from astropy import units as u
import wget
import tarfile
import glob
import gzip

import astropy.wcs as wcs

homedir = os.getenv("HOME")

os.sys.path.append(homedir+'/github/virgowise/')
import rungalfit as rg #This code has all the defined functions that I can use
os.sys.path.append(homedir+'/github/HalphaImaging/python3/')
import plot_cutouts_ha as cutouts #This code has all the defined functions that I can use




       
class galaxy():
   def __init__(self,ra,dec,size,name='galname',band='3'):
        '''
        galaxy for wise analysis

        params:
        -------
        * ra = ra of gal center in deg
        * dec = dec of gal center in deg
        * size = D25/2, or some other estimate of radius in arcsec

        optional params:
        ---------------
        * name = name of galaxy to use when saving images, default is galname
        * band = WISE band, default is '3'
          - 1 = W1
          - 2 = W2
          - 3 = W3
          - 4 = W4
        
        '''
        self.ra = ra#*u.deg
        self.dec = dec#*u.deg
        self.radius = size#*u.arcsec
        self.band = band
        self.galname = name
        self.image_rootname = self.galname+'-unwise-w'+str(self.band)
        self.image = self.image_rootname+'-img-m.fits'

        self.mask_image = self.galname+'-unwise-mask.fits'
        self.sigma_image = self.image_rootname+'-std-m.fits'
        self.invvar_image = self.image_rootname+'-invvar-m.fits'

        # remove the log file if it exists
        self.logfilename = self.galname+'-unwise-'+'w'+str(self.band)+'-log.txt'
        if os.path.exists(self.logfilename):
           os.remove(self.logfilename)

        # write the header line into the log file
        output=open(self.logfilename,'w')
        output.write('# xc xc_err yc yc_err mag mag_err re re_err nsersic nsrsic_err BA BA_err PA PA_err sky sky_err error chi2nu \n')
        # close log file
        output.close()

   def get_wise_image(self,makeplots=False):
        '''
        GOAL: Get the unWISE image from the unWISE catalog

        INPUT: nsaid used to grab unwise image information

        OUTPUT: Name of file to retrieve from

        '''


        baseurl = 'http://unwise.me/cutout_fits?version=allwise'
        imsize = self.radius*2

        imagenames,multiframe = cutouts.get_unwise_image(self.ra,self.dec,galid=self.galname,pixscale=1,imsize=self.radius*2,bands=self.band,makeplots=makeplots,subfolder=None)
        

        print(imagenames)
   def get_wise_image_old(self):
        '''
        GOAL: Get the unWISE image from the unWISE catalog

        INPUT: nsaid used to grab unwise image information

        OUTPUT: Name of file to retrieve from

        '''


        baseurl = 'http://unwise.me/cutout_fits?version=allwise'
        imsize = self.radius*2
        imurl = baseurl +'&ra=%.5f&dec=%.5f&size=%s&bands=%s'%(self.ra,self.dec,imsize,self.band)
        wisetar = wget.download(imurl)
        tartemp = tarfile.open(wisetar,mode='r:gz') #mode='r:gz'
        wnames = tartemp.getnames()
        self.wnames = wnames
        print(wnames)
        # check for multiple pointings - means galaxy is split between images
        self.multiframe = False
        if len(wnames) > 4:
            self.multiframe = True
            os.remove(wisetar)
            return

        wmembers = tartemp.getmembers()
        tartemp.extractall()
        for fname in wnames:
           t = fname.split('-')
           self.rename = self.galname+'-'+t[0]+'-'+t[2]+'-'+t[3]+'-'+t[4]

           #print self.rename
           if os.path.exists(self.rename): # this should only occur if multiple images are returned from wise
               os.remove(self.rename)
           os.rename(fname, self.rename)
           if self.rename.find('.gz') > -1:
              os.system('gunzip '+self.rename)
           if fname.find('img') > -1:
              self.inputimage = self.rename
        os.remove(wisetar)
        print('self.rename = ',self.rename)

        ##### DISPLAY IMAGE
        im = fits.getdata(self.rename)
        norm = simple_norm(im, stretch='asinh',percent=99)
        plt.imshow(im, norm=norm)
        plt.show()
   def set_image_names(self):
        '''
        GOAL:
        * Set the psf image name and some additional parameters for galfit 

        PARAMS:
        * none

        OUTPUT: PSF image and optional diffusion kernel, other parameters  

        '''
        #just using center til, doesn't matter usually
        self.psf_image = homedir+'/github/virgowise/wise_psfs/wise-w3-psf-wpro-09x09-05x05.fits' 
        self.psf_oversampling = 8
        #mask_image = 'testimage_mask.fits' no mask image 
        self.xminfit=0
        self.yminfit=0
        self.xmaxfit=100
        self.ymaxfit=100
        self.convolution_size=50

        ### NEED TO UPDATE MAGZP AND PSCALE FOR UNWISE
        self.magzp=22.5 #For WISE
        self.pscale=0.344 #For WISE
        
        #convflag=1 # apply psf convolution
        #constraintflag=1 # add a constraint file?
        #self.fitallflag=0
        self.ncomp=1

   def getpix(self):
        '''
        GOAL: Get pixel values that correspond to ra and dec

        INPUT: nsaid
        
        OUTPUT: xc, yc in pixels

        '''

        w=wcs.WCS(self.image)
        self.xc, self.yc = w.wcs_world2pix(self.ra, self.dec,1)


   def set_sersic_params(self):
        '''
        GOAL: Set random parameters for galfit

        INPUT: nsaid

        OUTPUT: 5 random parameters for nsersic, magnitude, effective radius, axis ratio, and position angle

        '''
        self.nsersic = 5.5*np.random.random()+.5
        self.mag =14*np.random.random()+2
        self.re = 60*np.random.random()
        self.BA = np.random.random()
        self.PA =181*np.random.random()-89.0

   def set_sersic_manual(self,n=2,m=7,re=5,BA=1,PA=0):
        '''
        GOAL: Set random parameters for galfit

        INPUT: nsaid

        OUTPUT: 5 random parameters for nsersic, magnitude, effective radius, axis ratio, and position angle

        '''
        self.nsersic = n
        self.mag =m
        self.re = re
        self.BA = BA
        self.PA = PA
                
   def initialize_galfit(self,convflag=True):
        '''
        GOAL: Preparing file to be run in galfit. Initialize galfit image parameters 

        INPUT: nsaid 

        OUTPUT: A definition of everything from galname to convflag, necessary for running galfit

        '''
        print('self.psfimage = ',self.psf_image)
        
        self.gal1 = rg.galfit(galname=self.image_rootname,image=self.image, mask_image = self.mask_image, sigma_image=self.sigma_image,psf_image=self.psf_image,psf_oversampling=self.psf_oversampling,xminfit=self.xminfit,yminfit=self.yminfit,xmaxfit=self.xmaxfit,ymaxfit=self.ymaxfit,convolution_size=self.convolution_size,magzp=self.magzp,pscale=self.pscale,ncomp=self.ncomp,convflag=convflag)
        
   def run_galfit_wise(self,fitBA=1,fitPA=1):
        '''
        GOAL: 
        * run galfit on one image

        optional params:
        ----------------
        * fitBA = set to 1 to let galfit fit the axis ratio BA of the galaxy
        * fitPA = set to 1 to let galfit fit the PA position angle of the galaxy

        OUTPUT: several output files

        '''
        #os.system('cp '+self.psf_image+' .')
        self.gal1.set_sersic_params(xobj=self.xc,yobj=self.yc,mag=self.mag,rad=self.re,nsersic=self.nsersic,BA=self.BA,PA=self.PA,fitmag=1,fitcenter=1,fitrad=1,fitBA=fitBA,fitPA=fitPA,fitn=1,first_time=0)
        self.gal1.set_sky(0)
        self.gal1.run_galfit()
   def get_galfit_results(self,printflag = False):
        '''
        GOAL: 
        -----
        * Grab results from galfit (xc, yc, mag, re, nsersic, BA, PA, sky, error, chi2nu) 
          and parse them into self.filename

        PARAMS:
        -------
        * none

        OPTIONAL PARAMS:
        ----------------
        * printflag = print fit results, default is False

        OUTPUT:
        -------
        * stores fit results in variables 

        '''

        self.filename = self.galname+'-unwise-'+'w'+str(self.band)+'-1Comp-galfit-out.fits'
        t = rg.parse_galfit_1comp(self.filename)
        if printflag:
            self.gal1.print_galfit_results(self.filename)
        
        header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_SKY','ERROR','CHI2NU']
        self.xc, self.xc_err = t[0]
        self.yc, self.yc_err = t[1]
        self.mag, self.mag_err = t[2]
        self.re, self.re_err = t[3]
        self.nsersic, self.nsersic_err = t[4]
        self.BA, self.BA_err = t[5]
        self.PA, self.PA_err = t[6]
        self.sky, self.sky_err = t[7]
        self.error = t[8]
        self.chi2nu = t[9]
        
   def write_results(self):
        '''
        GOAL: 
        * Put results from galfit into a logfile by appending values

        PARAMS:
        * none

        OUTPUT: 
        * logfile of parameter types and their associated outputs from galfit 

        '''
        self.get_galfit_results()

        output=open(self.logfilename,'a')
        # create string with best-fit parameters
        s = '%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f \n'%(self.xc,self.xc_err,self.yc,self.yc_err,self.mag,self.mag_err, self.re, self.re_err, self.nsersic, self.nsersic_err, self.BA, self.BA_err, self.PA, self.PA_err, self.sky, self.sky_err, self.error,self.chi2nu)
        output.write(s)
        output.close()
    

   def run_dmc(self, N=100,convflag=True):
        '''
        GOAL: 
        * Run galfit with monte carlo sampling to find all minima

        OPTIONAL PARAMS:
        * N = number of times to run galfit; default is 100
        * convflag = convolve galfit model with psf; default is True

        OUTPUT: 
        * galfit models for all possible minima

        '''
        #N is Number of random samples
        
        # set up arrays to store galfit output (e.g. xf, yf, rf, etc)
        # 7 output parameters, so Nx7 dimensional arrays
        X = np.empty((0,7))
        
        # uncertainty in fitted parameters
        dX = np.empty((0,7))
        C = np.empty((0)) #create list C (charge) with one vector (not array)

        # download the wise images if the user requests this
        if args.getwise:
            self.get_wise_image()

        # define image names
        self.set_image_names()

        # get the pixel coordinates of the galaxy
        # this uses the image header to translate RA and DEC into pixel coordinates
        self.getpix()

        # set up all of the inputs for galfit
        self.initialize_galfit(convflag=convflag)
        B = 1
        for i in range(N):
            E = 100000
            # this loop selects initial conditions that 
            while(np.random.random()>= np.exp(-B*E)): 
                nX = len(X) # number of unique fits that we already found
                D = np.zeros((1,7))             
                self.set_sersic_params() # select random initial conditions
                X0 = np.array([[self.xc,self.yc,self.mag,self.re,self.nsersic,self.BA,self.PA]])             
                for k in range(nX-1):
                    R = X0-X[k,:]
                    R7 = np.linalg.norm(R)**7
                    deltaD = C[k]*R/R7
                    D = D+deltaD
                E = np.linalg.norm(D)**2
            self.run_galfit_wise(fitBA=1,fitPA=1)
            self.get_galfit_results()
            # append best-fit values and errors to array         
            Qnew=np.array([[self.xc,self.yc,self.mag,self.re,self.nsersic,self.BA,self.PA]])
            dQnew = np.array([[self.xc_err,self.yc_err,self.mag_err,self.re_err,self.nsersic_err,self.BA_err,self.PA_err]])

            # initialize flag
            # the flag is 0 if it does not overlap the list of charges
            # the flag is 1 is it overlaps some other charge
            flag=0    
                    
            for j in range(nX): # loop over fitted variables (xc, yc, mag, re, nsersic, BA, PA)
                #compare the difference in fitted parameter with error in difference
                if (np.linalg.norm(Qnew-X[j,:])<= np.linalg.norm(dQnew+dX[j,:])): 
                    # overlaps
                    flag=1
                    C[j]=C[j]+1 # adding a charge
                    if np.linalg.norm(dQnew)<np.linalg.norm(dX[j,:]):
                        X[j,:]=Qnew
                        dX[j,:]=dQnew
                    break

            if flag==0: # if it does not overlap, then add it to the list
                if self.error==0:
                    X=np.append(X,Qnew,axis=0)
                    dX=np.append(dX,dQnew,axis=0)
                    C = np.append(C,[1],axis=0)
        return X    

   def run_simple(self, convflag=True):
        '''
        GOAL: 
        * Run galfit once 

        OPTIONAL PARAMS:
        * convflag = convolve galfit model with psf; default is True

        OUTPUT: 
        * galfit model

        '''
        # download the wise images if the user requests this
        self.get_wise_image()

        # define image names
        self.set_image_names()

        # get the pixel coordinates of the galaxy
        # this uses the image header to translate RA and DEC into pixel coordinates
        self.getpix()

        # set up all of the inputs for galfit
        self.initialize_galfit(convflag=convflag)
        self.set_sersic_params() # select random initial conditions
        self.run_galfit_wise(fitBA=1,fitPA=1)
        self.write_results()



    

    
