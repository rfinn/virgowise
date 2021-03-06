
"""
GOAL:
- to retrieve unWISE images of galaxies given their ra, dec, image size
- run galfit

PROCEDURE:
- download unwise image for particular band as a tar file


USEAGE:

from within ipython on laptop


REQUIRED FILES


NOTES:
* updating Aug 2020 to remove randomization of initial starting parameters
* just want to take in an RA and DEC (with option to provide and fix the BA and PA)
  - download WISE image
  - run galfit on W3 (and W4?)
  - return results


DONE:
* down the road, need to handle galaxies that have multiple unwise images
  - combine the images - DONE
  - make a cutout from the combined image
"""


import os
#import pyds9
import numpy as np
import argparse
from astropy.io import fits
from astropy.visualization import simple_norm
from astropy import units as u
from astropy.wcs import WCS
from matplotlib import pyplot as plt
from scipy.stats import scoreatpercentile
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



UNWISE_PIXSCALE = 2.75
LEGACY_PIXSCALE = 1
RSCALEFACTOR = 3
       
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

        # this function takes the wise image size in pixels
        imagenames,weightnames,multiframe = cutouts.get_unwise_image(self.ra,self.dec,galid=self.galname,pixscale=1,imsize=self.radius*RSCALEFACTOR/UNWISE_PIXSCALE,bands=self.band,makeplots=makeplots,subfolder=None)
        
        self.image = imagenames[0]
        
        self.sigma_image = weightnames[0]

        # read in image and get image size

        temp = fits.getdata(self.image)
        print(temp.shape)
        self.ximagesize,self.yimagesize = temp.shape
        
        #print(self.sigma_image)
        #print(imagenames)
   def get_legacy_images(self,makeplots=False,band='r'):
        '''
        GOAL: Get the legacy jpeg image 

        INPUT: nsaid used to grab unwise image information

        OUTPUT: Name of file to retrieve from

        '''


        self.legacy_filename_r,self.legacy_jpegname = cutouts.get_legacy_images(self.ra,self.dec,galid=self.galname,pixscale=1,imsize=self.radius*RSCALEFACTOR,bands=band,subfolder=None)
        

        # read in image and get image size
        
        #print(self.sigma_image)
        #print(imagenames)
        

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
        os.system('cp '+homedir+'/github/virgowise/wise_psfs/wise-w3-psf-wpro-09x09-05x05.fits .') 
        self.psf_image = 'wise-w3-psf-wpro-09x09-05x05.fits' 
        self.psf_oversampling = 8
        #mask_image = 'testimage_mask.fits' no mask image 
        self.xminfit=0
        self.yminfit=0
        self.xmaxfit=self.ximagesize
        self.ymaxfit=self.yimagesize
        self.convolution_size=self.ximagesize

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


   def set_sersic_params_random(self):
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

        self.filename = self.gal1.output_image
        #self.galname+'-unwise-'+'w'+str(self.band)+'-1Comp-galfit-out.fits'
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
        
   def write_results(self,printflag=False):
        '''
        GOAL: 
        * Put results from galfit into a logfile by appending values

        PARAMS:
        * none

        OUTPUT: 
        * logfile of parameter types and their associated outputs from galfit 

        '''
        self.get_galfit_results(printflag=printflag)

        output=open(self.logfilename,'a')
        # create string with best-fit parameters
        s = '%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f \n'%(self.xc,self.xc_err,self.yc,self.yc_err,self.mag,self.mag_err, self.re, self.re_err, self.nsersic, self.nsersic_err, self.BA, self.BA_err, self.PA, self.PA_err, self.sky, self.sky_err, self.error,self.chi2nu)
        output.write(s)
        output.close()
   def display_galfit_model(self,percentile1=.5,percentile2=99.5,p1residual=5,p2residual=99,cmap='viridis',zoom=None):
      '''
      ARGS:
      percentile1 = min percentile for stretch of image and model
      percentile2 = max percentile for stretch of image and model
      p1residual = min percentile for stretch of residual
      p2residual = max percentile for stretch of residual
      cmap = colormap, default is viridis
      '''
      # model name
      self.filename = self.galname+'-unwise-'+'w'+str(self.band)+'-1Comp-galfit-out.fits'
      pngname = self.galname+'-unwise-'+'w'+str(self.band)+'-1Comp-galfit-out.png'
      if self.convflag:
         pngname = self.galname+'-unwise-'+'w'+str(self.band)+'-1Comp-galfit-out-conv.png'
      image,h = fits.getdata(self.filename,1,header=True)
      model = fits.getdata(self.filename,2)
      residual = fits.getdata(self.filename,3)

      if zoom is not None:
         print("who's zoomin' who?")
         # display central region of image
         # figure out how to zoom

         # get image dimensions and center
         xmax,ymax = image.shape
         xcenter = int(xmax/2)
         ycenter = int(ymax/2)

         # calculate new size to display based on zoom factor
         new_xradius = int(xmax/2/(float(zoom)))
         new_yradius = int(ymax/2/(float(zoom)))

         # calculate pixels to keep based on zoom factor
         x1 = xcenter - new_xradius
         x2 = xcenter + new_xradius
         y1 = ycenter - new_yradius
         y2 = ycenter + new_yradius
         
         # check to make sure limits are not outsize image dimensions
         if (x1 < 1):
            x1 = 1
         if (y1 < 1):
            y1 = 1
         if (x2 > xmax):
            x2 = xmax
         if (y2 > ymax):
            y2 = ymax

         # cut images to new size
         image = image[x1:x2,y1:y2]
         model = model[x1:x2,y1:y2]
         residual = residual[x1:x2,y1:y2]         
         pass
      wcs = WCS(h)
      images = [image,model,residual]
      titles = ['image','model','residual']
      v1 = [scoreatpercentile(image,percentile1),
            scoreatpercentile(image,percentile1),
            scoreatpercentile(residual,p1residual)]
      v2 = [scoreatpercentile(image,percentile2),
            scoreatpercentile(image,percentile2),
            scoreatpercentile(residual,p2residual)]
      norms = [simple_norm(image,'asinh',max_percent=percentile2),
               simple_norm(image,'asinh',max_percent=percentile2),
               simple_norm(residual,'linear',max_percent=p2residual)]
               
      plt.figure(figsize=(14,6))
      plt.subplots_adjust(wspace=.0)
      for i,im in enumerate(images): 
         plt.subplot(1,3,i+1,projection=wcs)
         plt.imshow(im,origin='lower',cmap=cmap,vmin=v1[i],vmax=v2[i],norm=norms[i])
         plt.xlabel('RA')
         if i == 0:
            plt.ylabel('DEC')
         else:
            ax = plt.gca()
            ax.set_yticks([])
         plt.title(titles[i],fontsize=16)
      plt.savefig(pngname)
   def print_galfit_results(self):
      #self.filename = self.galname+'-unwise-'+'w'+str(self.band)+'-1Comp-galfit-out.fits'
      self.filename = self.gal1.output_image
      rg.print_galfit_results(self.filename)
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

   def run_simple(self, convflag=True,zoom=1,sersic_start=None):
        '''
        GOAL: 
        * Run galfit once 

        OPTIONAL PARAMS:
        * convflag = convolve galfit model with psf; default is True

        OUTPUT: 
        * galfit model

        '''
        self.convflag = convflag
        # download the wise images if the user requests this
        self.get_wise_image()

        # define image names
        self.set_image_names()

        # get the pixel coordinates of the galaxy
        # this uses the image header to translate RA and DEC into pixel coordinates
        self.getpix()

        # set up all of the inputs for galfit
        self.initialize_galfit(convflag=convflag)
        if sersic_start is not None:
           xc,yc,mag,re,n,BA,PA =  sersic_start
           self.set_sersic_manual(m=mag,re=re,n=n,BA=BA,PA=PA)
           # set fixed initial parameters
        else:
           self.set_sersic_manual()
      
        self.run_galfit_wise(fitBA=1,fitPA=1)
        self.write_results(printflag=True)
        #self.get_galfit_f
        self.display_galfit_model(zoom=zoom)


    

    
