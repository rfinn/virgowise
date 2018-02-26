

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
parser.add_argument('--t',dest = 't', default =' ~/github/Virgo/tables', help = 'Locates fits catalog in fits tables')


os.sys.path.append('/Users/rfinn/github/Virgo/programs/')#Dr.Finn local path

os.sys.path.append('~/github/Virgo/programs/')
from rungalfit import * #This code has all the defined functions that I can use
#Get catalog files
os.system('cp' + args.t + '/nsa.virgo.fits')
os.system('cp' + args.t + '/nsa_wise.virgo')
os.system('cp' + args.t + '/nsa_CO-HI.virgo.fits')

#import catalogs from tables folder in Virgo Github

Nsa = fits.getdata(nsa.virgo.fits)
Wise = fits.getdata(nsa_wise.virgo.fits)
Co = fits.getdata(nsa_CO-HI.virgo.fits)



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






