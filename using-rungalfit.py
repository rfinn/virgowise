

#  Goal  
# 
# The goal of this notebook is to show you how to use the galfit class in rungalfit.py
# 
# ### Required Modules ###
# 
# pyds9


#Natasha goal for next week: TO SET UP A PARAMETER INPUT FILE AND TRY EVENTUALLY FOR GALFIT
import os
import pyds9
import argparse
from astropy.io import fits


#Need user to define galaxy image/sigma/psf path later on
parser = argparse.ArgumentParser(description ='Run galfit and store output with best fit parameters into a tar file')
parser.add_argument('--l',dest = 'l', default =' /home/share/research/Virgo/galfitexample/WISEmodels/unwise-*p*-w3-*-m.fits', help = 'Locates list of images, sigma image, and psf image  of galaxy/galaxies path')
parser.add_argument('--t',dest = 't', default ='~/github/Virgo/tables/', help = 'Locates fits catalog in fits tables')


os.sys.path.append('/Users/rfinn/github/Virgo/programs/')#Dr.Finn local path

os.sys.path.append('~/github/Virgo/programs/')
from rungalfit import * #This code has all the defined functions that I can use
#Get catalog files
#os.system('cp ' + args.t + '/nsa.virgo.fits')
#os.system('cp ' + args.t + '/nsa_wise.virgo')
#os.system('cp ' + args.t + '/nsa_CO-HI.virgo.fits')

#import catalogs from tables folder in Virgo Github

'''
nsa = fits.getdata(args.t+'nsa.virgo.fits')
wise = fits.getdata(args.t+'nsa_wise.virgo.fits')
co = fits.getdata(args.t+'nsa_CO-HI.virgo.fits')



# define image properties
galname = 'test'
image = 'testimage.fits'
sigma_image = 'testimage-sigma.fits'
psf_image = 'testpsf.fits'
psf_oversampling = 1
mask_image = 'testimage_mask.fits'
xminfit=1
yminfit=1
xmaxfit=100
ymaxfit=100
convolution_size=100
magzp=20.
pscale=2.5
convflag=1 # apply psf convolution
constraintflag=1 # add a constraint file?
fitallflag=0
ncomp=1



# define first guess sersic parameters for galaxy 
xc=50
yc=50
nsersic = 2
mag = 7
Re = 10.
BA = .8
PA = 0


gal1 = galfit(galname=galname,image=image,sigma_image=sigma_image,psf_image=psf_image,psf_oversampling=psf_oversampling,mask_image=mask_image,xminfit=xminfit,yminfit=yminfit,xmaxfit=xmaxfit,ymaxfit=ymaxfit,convolution_size=convolution_size,magzp=magzp,pscale=pscale,convflag=convflag,constraintflag=constraintflag,fitallflag=fitallflag,ncomp=ncomp)


## Create an input file for galfit ##
#for args.l; want tofeed in args.l for list of input file?

sampleflag = (wise.W3SNR>10) & (co.CO_DETECT==1)
#Make a flag
if (wise.W3SNR>10)&(Co.CO_DETECT==1):
    images = sorted(glob.glob(args.l))#grabs each individual file with flag


# create output names using functions defined in rungalfit
gal1.create_output_names()
gal1.open_galfit_input()
gal1.write_image_params()
gal1.add_simple_sersic_object(1,'sersic',xc,yc,mag,Re,nsersic,BA,PA)
gal1.set_sky(0)
gal1.write_sky(2)
gal1.close_input_file()



'''


class galaxy():
   def __init__(self,catalog_path):
       # read in nsa, wise, co catalogs
       self.nsa = catalog_path+'nsa.virgo.fits'
       #self.wise
       #self.co
       self.nsadict=dict((a,b) for a,b in zip(self.nsa.NSAID,arange(len(self.nsa.NSAID))))
    def define_sample(self):
        self.sampleflag = (self.wise.W3SNR>10) & (self.co.CO_DETECT==1)
    def get_wise_image(self,nsaid):
        # download wise image
        # adapt your code to fill out this function
        # rename file as NSAID-W3-   get rid of RA and DEC
    def get_initial_from_sextractor(self):
        # we will do this together

        continue
    
    def set_image_names(self,nsaid):
        # fix this to match the actual filenames
        self.galname = nsaid+'-test'
        self.image = 'testimage.fits'
        sigma_image = 'testimage-sigma.fits'
        psf_image = 'testpsf.fits'
        psf_oversampling = 1
        mask_image = 'testimage_mask.fits'
        xminfit=1
        yminfit=1
        xmaxfit=100
        ymaxfit=100
        convolution_size=100
        magzp=20.
        pscale=2.5
        convflag=1 # apply psf convolution
        constraintflag=1 # add a constraint file?
        fitallflag=0
        ncomp=1
    def set_sersic_params(self):
        # define first guess sersic parameters for galaxy 
        xc=50
        yc=50
        nsersic = 2
        mag = 7
        Re = 10.
        BA = .8
        PA = 0
        
    def initialize_galfit(self):
        self.gal1 = galfit(galname=self.galname,image=self.image,sigma_image=sigma_image,psf_image=psf_image,psf_oversampling=psf_oversampling,mask_image=mask_image,xminfit=xminfit,yminfit=yminfit,xmaxfit=xmaxfit,ymaxfit=ymaxfit,convolution_size=convolution_size,magzp=magzp,pscale=pscale,convflag=convflag,constraintflag=constraintflag,fitallflag=fitallflag,ncomp=ncomp)
    def run_galfit(self):

        self.gal1.create_output_names()
        self.gal1.open_galfit_input()
        self.gal1.write_image_params()
        self.gal1.add_simple_sersic_object(1,'sersic',xc,yc,mag,Re,nsersic,BA,PA)
        self.gal1.set_sky(0)
        self.gal1.write_sky(2)
        self.gal1.close_input_file()


############ MAIN PROGRAM ###############
if __name__ == "__main__":
        
    mygals = galaxy(args.t)
    mygals.define_sample()

    galaxy_index = np.arange(len(self.nsa.RA))[mygals.sampleflag]
    mygals.get_wise_image(mygals.nsa.NSAID[galaxy_index[0]])

