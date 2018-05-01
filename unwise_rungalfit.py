
"""
GOAL:
- to retrieve unWISE images of galaxies given their NSA id

PROCEDURE:
- get list of NSAIDs
- read in NSA
- download unwise image for particular band as a tar file
- extract the images from the tar file
- rename the unwise images to something more user friendly
- delete the tar file

USEAGE:

from within ipython on laptop

%run ~/github/virgowise/using-rungalfit.py --nsapath '/Users/rfinn/github/Virgo/tables/' 

REQUIRED FILES:
- nsa fits catalog


NOTES:
"""


#Natasha goal for next week: TO SET UP A PARAMETER INPUT FILE AND TRY EVENTUALLY FOR GALFIT
import os
import pyds9
import numpy as np
import argparse
from astropy.io import fits
import wget
import tarfile
import glob
import gzip

#Need user to define galaxy image/sigma/psf path later on
parser = argparse.ArgumentParser(description ='Run galfit and store output with best fit parameters into a tar file')
parser.add_argument('--band', dest = 'band', default = '3', help = 'unWISE image band to download: 3=12um, 4=22um (can only do one at a time for now.')
parser.add_argument('--nsapath',dest = 'nsapath', default ='~/github/Virgo/tables/', help = 'Location of NSA tables fits')
parser.add_argument('--display',dest = 'display', default =True, help = 'display galfit results in ds9?  default = True')

args = parser.parse_args()

#os.sys.path.append('~/github/Virgo/programs/')
os.sys.path.append('~/github/virgowise/')
from rungalfit import * #This code has all the defined functions that I can use



class catalogs():
   def __init__(self,catalog_path):
       # read in nsa, wise, co catalogs
       self.nsatab = catalog_path + 'nsa.virgo.fits'
       self.wisetab = catalog_path + 'nsa_wise.virgo.fits'
       self.cotab = catalog_path + 'nsa_CO-Gianluca.virgo.fits'
       self.nsa = fits.getdata(self.nsatab)
       self.wise = fits.getdata(self.wisetab)
       self.co = fits.getdata(self.cotab)
       self.nsadict=dict((a,b) for a,b in zip(self.nsa.NSAID,np.arange(len(self.nsa.NSAID)))) #useful for us can easily look up galaxy ID's       
   def define_sample(self):
       #self.sampleflag = (self.wise.W3SNR>10) & (self.co.CO_DETECT==1)   
       self.w3_flag = self.wise.W3SNR>10 
       self.w4_flag = self.wise.W4SNR>5      
       self.co_flag = self.co.COdetected == '1'
       self.sampleflag = self.w3_flag & self.w4_flag & self.co_flag
       print 'number of galaxies in sample = ',sum(self.sampleflag)

       
class galaxy():
   def __init__(self,nsaid,band='3'):
        print 'hello galaxy NSA ',nsaid
        self.nsaid = nsaid
        self.band = band
        self.image_rootname = 'NSA-'+str(self.nsaid)+'-unwise-w'+str(self.band)
        self.image = self.image_rootname+'-img-m.fits'
        self.sigma_image = self.image_rootname+'-std-m.fits'
        self.invvar_image = self.image_rootname+'-invvar-m.fits'
   def get_wise_image(self):
        galindex = cats.nsadict[self.nsaid]
        baseurl = 'http://unwise.me/cutout_fits?version=allwise'
        imsize = '100'
        imurl = baseurl +'&ra=%.5f&dec=%.5f&size=%s&bands=%s'%(cats.nsa.RA[galindex],cats.nsa.DEC[galindex],imsize,self.band)
        wisetar = wget.download(imurl)
        tartemp = tarfile.open(wisetar,mode='r:gz') #mode='r:gz'
        wnames = tartemp.getnames()
        print wnames
        wmembers = tartemp.getmembers()
        tartemp.extractall()
        for fname in wnames:


           t = fname.split('-')
           
           self.rename = 'NSA-'+str(self.nsaid)+'-'+t[0]+'-'+t[2]+'-'+t[3]+'-'+t[4]

           print self.rename
           os.rename(fname, self.rename)
           if self.rename.find('.gz') > -1:
              os.system('gunzip '+self.rename)
           if fname.find('img') > -1:
              self.inputimage = self.rename
        os.remove(wisetar)
        print 'self.rename = ',self.rename
    
   def set_image_names(self,nsaid):
        self.psf_image = 'wise-w3-psf-wpro-09x09-05x05.fits' #just using center til, doesn't matter usually
        self.psf_oversampling = 80
        #mask_image = 'testimage_mask.fits' no mask image 
        self.xminfit=0
        self.yminfit=0
        self.xmaxfit=100
        self.ymaxfit=100
        self.convolution_size=100

        ### NEED TO UPDATE MAGZP AND PSCALE FOR UNWISE
        self.magzp=22.5 #For WISE
        self.pscale=0.344 #For WISE




        
        #convflag=1 # apply psf convolution
        #constraintflag=1 # add a constraint file?
        #self.fitallflag=0
        self.ncomp=1
   def set_sersic_params(self):
        # define first guess sersic parameters for galaxy 
        self.xc=50
        self.yc=50
        self.nsersic = 2
        self.mag = 7
        self.Re = 0.5
        self.BA = .4
        self.PA = 0
        
   def initialize_galfit(self,nsaid):
        self.gal1 = galfit(galname=self.image_rootname,image=self.image,sigma_image=self.sigma_image,psf_image=self.psf_image,psf_oversampling=self.psf_oversampling,xminfit=self.xminfit,yminfit=self.yminfit,xmaxfit=self.xmaxfit,ymaxfit=self.ymaxfit,convolution_size=self.convolution_size,magzp=self.magzp,pscale=self.pscale,ncomp=self.ncomp)
        
   def run_galfit_wise(self):
        self.gal1.set_sersic_params(xobj=self.xc,yobj=self.yc,mag=self.mag,rad=self.Re,nsersic=self.nsersic,BA=self.BA,PA=self.PA,fitmag=1,fitcenter=1,fitrad=1,fitBA=1,fitPA=1,fitn=1,first_time=0)
        self.gal1.set_sky(0)
        self.gal1.run_galfit()
        #self.gal1.display_results()
        #self.gal1.close_input_file()
        #self.gal1.print_params()
        #self.gal1.print_galfit_results()


############ MAIN PROGRAM ###############
if __name__ == "__main__":
    # READ IN CATALOGS
    cats = catalogs(args.nsapath)
    cats.define_sample()

    
    ngc_filament_ids = [56403,56409,56410,56411,56434,56455,56456,56462,56469,56482,56489,61690,61691,61692,67593,88353,90371,93403,94217,104307,104439,118647,119230,119289,120018,120053,142509,143682,143686,143827,143951,143986,162674,163136,163783,164224,164358]    
#119303 taking out bc bad image

    listname = ngc_filament_ids
    listname = cats.nsa.NSAID[cats.sampleflag]

    pause_flag = True
    for mynsaid in listname:
        mygals = galaxy(mynsaid,band=args.band)
        mygals.get_wise_image()
        mygals.set_image_names(mynsaid)
        mygals.set_sersic_params()
        mygals.initialize_galfit(mynsaid)
        mygals.run_galfit_wise()
        if pause_flag:
            t = raw_input('hit any key to continue to next galaxy \n \t enter q to quit \n \t enter C to continue without pausing \n')

            if t.find('q') > -1:
                break
            elif t.find('C') > -1:
                pause_flag = False
    
    #galaxy_index = np.arange(len(self.nsa.RA))[mygals.sampleflag]
    #mygals.get_wise_image(mygals.nsa.NSAID[galaxy_index[0]])

    
