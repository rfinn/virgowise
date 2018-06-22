
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

%run ~/github/virgowise/using-rungalfit.py --virgo '/Users/rfinn/github/Virgo/' 

For comparison with LCS galaxies

%run ~/Dropbox/pythonCode/LCSanalyzeblue.py

lcs_listname = s.s.NSAID[s.sampleflag]
band = '3'

%run /Users/rfinn/github/virgowise/unwise_rungalfit.py --nsapath '/Users/rfinn/research/NSA/' --nsafile 'full' --virgopath '/Users/rfinn/github/Virgo/'

t = process_list(lcs_listname,'3',convolution_flag = False)

REQUIRED FILES:
- nsa fits catalog


NOTES:
"""


#Natasha goal for next week: TO SET UP A PARAMETER INPUT FILE AND TRY EVENTUALLY FOR GALFIT
import os
#import pyds9
import numpy as np
import argparse
from astropy.io import fits
import wget
import tarfile
import glob
import gzip

import astropy.wcs as wcs

#Need user to define galaxy image/sigma/psf path later on
parser = argparse.ArgumentParser(description ='Run galfit and store output with best fit parameters into a tar file')
parser.add_argument('--band', dest = 'band', default = '3', help = 'unWISE image band to download: 3=12um, 4=22um (can only do one at a time for now.')
parser.add_argument('--virgopath',dest = 'virgopath', default ='/Users/rfinn/github/Virgo/', help = 'Location of NSA tables fits')
parser.add_argument('--nsapath',dest = 'nsapath', default ='/Users/rfinn/github/Virgo/tables/', help = 'Location of NSA tables fits')
parser.add_argument('--nsafile',dest = 'nsafile', default ='virgo', help = 'nsa file to use.  could be nsa.virgo.fits (default) or use full to get nsa_v0_1_2.fits.  make sure nsapath points to the right place.')
parser.add_argument('--display',dest = 'display', default =True, help = 'display galfit results in ds9?  default = True')
parser.add_argument('--getwise',dest = 'getwise', default =False, help = 'download wise images?  default is False.')



#os.sys.path.append('~/github/Virgo/programs/')
os.sys.path.append('/Users/rfinn/github/virgowise/')
from rungalfit import * #This code has all the defined functions that I can use



class catalogs():
   def __init__(self,catalog_path,virgoflag = False):
       # read in nsa, wise, co catalogs
       if virgoflag:
            self.nsatab = catalog_path + 'nsa.virgo.fits'
            self.wisetab = catalog_path + 'nsa_wise.virgo.fits'
            self.cotab = catalog_path + 'nsa_CO-Gianluca.virgo.fits'
            self.co = fits.getdata(self.cotab)
       else:
            self.nsatab = catalog_path + 'nsa_v0_1_2.fits'
            self.wisetab = catalog_path + 'nsa_v0_1_2_wise.fits'
       self.nsa = fits.getdata(self.nsatab)
       self.wise = fits.getdata(self.wisetab)

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
        self.mask_image = 'masks/'+self.image.split('.fits')[0]+'-mask.fits'
        self.sigma_image = self.image_rootname+'-std-m.fits'
        self.invvar_image = self.image_rootname+'-invvar-m.fits'

        # add code to remove the log file if it exists\
        logfilename = 'NSA-'+str(self.nsaid)+'-unwise-'+'w'+str(self.band)+'-log.txt'
        if os.path.exists(logfilename):
           os.remove(logfilename)

        # add code to write the header line into the log file
        output=open(logfilename,'w')
        output.write('# xc xc_err yc yc_err mag mag_err re re_err nsersic nsrsic_err BA BA_err PA PA_err sky sky_err error chi2nu \n')
        # close log file
        output.close

   def get_wise_image(self):
        galindex = cats.nsadict[self.nsaid]
        baseurl = 'http://unwise.me/cutout_fits?version=allwise'
        imsize = '100'
        imurl = baseurl +'&ra=%.5f&dec=%.5f&size=%s&bands=%s'%(cats.nsa.RA[galindex],cats.nsa.DEC[galindex],imsize,self.band)
        wisetar = wget.download(imurl)
        tartemp = tarfile.open(wisetar,mode='r:gz') #mode='r:gz'
        wnames = tartemp.getnames()
        print wnames
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
           self.rename = 'NSA-'+str(self.nsaid)+'-'+t[0]+'-'+t[2]+'-'+t[3]+'-'+t[4]

           print self.rename
           if os.path.exists(self.rename): # this should only occur if multiple images are returned from wise
               os.remove(self.rename)
           os.rename(fname, self.rename)
           if self.rename.find('.gz') > -1:
              os.system('gunzip '+self.rename)
           if fname.find('img') > -1:
              self.inputimage = self.rename
        os.remove(wisetar)
        print 'self.rename = ',self.rename
    
   def set_image_names(self):
        self.psf_image = 'wise-w3-psf-wpro-09x09-05x05.fits' #just using center til, doesn't matter usually
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
        galindex = cats.nsadict[self.nsaid]
        self.RA = cats.nsa.RA[galindex]
        self.DEC = cats.nsa.DEC[galindex]
        w=wcs.WCS(self.image)
        print self.RA, self.DEC
        self.xc, self.yc = w.wcs_world2pix(self.RA, self.DEC,1)
        print self.xc, self.yc

#   def set_sersic_params(self):
#        # define first guess sersic parameters for galaxy 
#        #self.xc=50
#        #self.yc=50
#        self.nsersic = 2
#        self.mag = 7
#        self.re = 5
#        self.BA = 1
#        self.PA = 0
#        self.BA = cats.nsa.SERSIC_BA[cats.nsadict[self.nsaid]]
#        self.PA = cats.nsa.SERSIC_PHI[cats.nsadict[self.nsaid]]
   def set_sersic_params(self):
        self.nsersic = 5.5*np.random.random()+.5
        self.mag =14*np.random.random()+2
        self.re = 60*np.random.random()
        self.BA = np.random.random()
        self.PA =181*np.random.random()-89.0
                
   def initialize_galfit(self,convflag=True):
        print 'self.psfimage = ',self.psf_image
        self.gal1 = galfit(galname=self.image_rootname,image=self.image, mask_image = self.mask_image, sigma_image=self.sigma_image,psf_image=self.psf_image,psf_oversampling=self.psf_oversampling,xminfit=self.xminfit,yminfit=self.yminfit,xmaxfit=self.xmaxfit,ymaxfit=self.ymaxfit,convolution_size=self.convolution_size,magzp=self.magzp,pscale=self.pscale,ncomp=self.ncomp,convflag=convflag)
        
   def run_galfit_wise(self,fitBA=1,fitPA=1):
        os.system('cp '+args.virgopath+'wisepsf/'+self.psf_image+' .')
        self.gal1.set_sersic_params(xobj=self.xc,yobj=self.yc,mag=self.mag,rad=self.re,nsersic=self.nsersic,BA=self.BA,PA=self.PA,fitmag=1,fitcenter=1,fitrad=1,fitBA=fitBA,fitPA=fitPA,fitn=1,first_time=0)
        self.gal1.set_sky(0)
        self.gal1.run_galfit()
        #self.gal1.display_results()
        #self.gal1.close_input_file()
        #self.gal1.print_params()
        #self.gal1.print_galfit_results()
   def get_galfit_results(self,printflag = False):

        self.filename = 'NSA-'+str(self.nsaid)+'-unwise-'+'w'+str(self.band)+'-1Comp-noconv-fitBAPA-galfit-out.fits'
        t = parse_galfit_1comp(self.filename)
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
        self.get_galfit_results()
        # Write a logfile 
        logfilename = 'NSA-'+str(self.nsaid)+'-unwise-'+'w'+str(self.band)+'-log.txt'
        output=open(logfilename,'a')
        # create string with best-fit parameters
        s = '%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f \n'%(self.xc,self.xc_err,self.yc,self.yc_err,self.mag,self.mag_err, self.re, self.re_err, self.nsersic, self.nsersic_err, self.BA, self.BA_err, self.PA, self.PA_err, self.sky, self.sky_err, self.error,self.chi2nu)
        output.write(s)
        output.close()

def process_list(listname,band,convolution_flag=True,getwise=True):

    pause_flag = True
    multiframe = np.zeros(len(listname),'bool')
    i=-1
    for mynsaid in listname:
        i += 1
        mygals = galaxy(mynsaid,band=band)
        if getwise:
            mygals.get_wise_image()
            multiframe[i] = mygals.multiframe
            if mygals.multiframe:
                print '\n NSA ',mynsaid,' is on multiple frames - skipping for now \n \ntarf'
                continue
        mygals.set_image_names()
        mygals.getpix()
        mygals.set_sersic_params()
        # set PA and BA to NSA values
        # fix these values
        mygals.initialize_galfit(convflag=0)
        mygals.run_galfit_wise(fitBA=1,fitPA=1)
        galfile = 'NSA-'+str(mygals.nsaid)+'-unwise-'+'w'+str(mygals.band)+'-1Comp-galfit-out.fits'
        altfilename = 'NSA-'+str(mygals.nsaid)+'-unwise-'+'w'+str(mygals.band)+'-1Comp-noconv-fitBAPA-galfit-out.fits'
        if os.path.exists(altfilename):
            os.remove(altfilename)
        os.rename(galfile,altfilename)
        mygals.run_galfit_wise(fitBA=0,fitPA=0)
        altfilename = 'NSA-'+str(mygals.nsaid)+'-unwise-'+'w'+str(mygals.band)+'-1Comp-noconv-galfit-out.fits'
        if os.path.exists(altfilename):
            os.remove(altfilename)
        os.rename(galfile,altfilename)
        # get output from no convolution - use mag and Re as input with convolution

        # skipping convolution for now
        # something is not set right
        # might be the oversampling number
        if convolution_flag:
            mygals.get_galfit_results()
            
            # keep PA and BA fixed to NSA values as we did with LCS
            mygals.PA = cats.nsa.SERSIC_PHI[cats.nsadict[mygals.nsaid]]
            mygals.BA = cats.nsa.SERSIC_BA[cats.nsadict[mygals.nsaid]]
            mygals.initialize_galfit(mynsaid)
            mygals.run_galfit_wise(fitBA=0,fitPA=0)
        if mygals.error>0: ###############Stops the program if galfit has a convergence issue
           break

        mygals.write_results()
        if pause_flag:
            t = raw_input('hit any key to continue to next galaxy \n \t enter q to quit \n \t enter C to continue without pausing \n')

            if t.find('q') > -1:
                break
            elif t.find('C') > -1:
                pause_flag = False
    return multiframe



############ MAIN PROGRAM ###############
if __name__ == "__main__":
    
    # READ IN CATALOGS
    args = parser.parse_args()
    if args.nsafile == 'virgo':
        cats = catalogs(args.nsapath,virgoflag=True)
    else:
        cats = catalogs(args.nsapath,virgoflag=False)
    if args.nsafile == 'virgo':
        cats.define_sample()
        listname = cats.nsa.NSAID[cats.sampleflag]

    #lcs = fits.getdata('/Users/rfinn/research/LocalClusters/NSAmastertables/LCS_all_size.fits')
    #ngc_filament_ids = [56403,56409,56410,56411,56434,56455,56456,56462,56469,56482,56489,61690,61691,61692,67593,88353,90371,93403,94217,104307,104439,118647,119230,119289,120018,120053,142509,143682,143686,143827,143951,143986,162674,163136,163783,164224,164358]    
#119303 taking out bc bad image

    #listname = ngc_filament_ids

    # call the line below to get started
    #process_list(listname,band=args.band)
    
    #galaxy_index = np.arange(len(self.nsa.RA))[mygals.sampleflag]
    #mygals.get_wise_image(mygals.nsa.NSAID[galaxy_index[0]])

    
