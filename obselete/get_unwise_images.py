#! /usr/bin/env python

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

%run ~/github/virgowise/get-unwise-images.py --nsapath '/Users/rfinn/github/Virgo/tables/' --band 3

REQUIRED FILES:
- nsa fits catalog


NOTES:
I want to download the unWISE images with the same NGC number as the Halpha images. I can do this by splitting up the radial output filenames and indexing just the nsatlas number of the galaxy,gonna need a way to convert that to the number in NED but for now: 

-1/28/18 I want to somehow download the tar files from unWISE but know what galaxy is which by outputting the nsatlas numbers which is what our Halpha images are in. Also commenting details of code.

5/1/18 - RAF Cleaning up code to remove remnants of Halpha.  Program now downloads images for a list of NSA ids.

"""

import os
import numpy as np
from astropy.io import fits
import wget
import tarfile
import argparse
import gzip


parser = argparse.ArgumentParser(description ='Download unWISE images for specific NSA galaxies in current directory')
parser.add_argument('--band', dest = 'band', default = '3', help = 'unWISE image band to download: 3=12um, 4=22um (can only do one at a time for now.')
parser.add_argument('--nsapath', dest = 'nsapath', default = '/Users/rfinn/github/Virgo/tables/', help = 'full path to Virgo NSA catalog (nsa.virgo.fits)')
args = parser.parse_args() #brings in these arguments above


###
### READ IN NSA CATALOG
###

nsatab = args.nsapath+'nsa.virgo.fits'
nsa = fits.getdata(nsatab)
nsadict = dict((a,b) for a,b in zip(nsa.NSAID,np.arange(len(nsa.NSAID))))

###
### GET NSA IDS FOR IMAGES TO DOWNLOAD
###

# hardcoding the list of NSA ids
# will make this smarter over the summer
# could run this on ALL galaxies in NSA virgo catalog!

#nsaid = ['56462','67595','164358','54578','61690','61692']
#nsaid = ['143986']

ngc_filament_ids = [56403,56409,56410,56411,56434,56455,56456,56462,56469,56482,56489,61690,61691,61692,67593,88353,90371,93403,94217,104307,104439,118647,119230,119289,120018,120053,142509,143682,143686,143827,143951,143986,162674,163136,163783,164224,164358]

nsaid = ngc_filament_ids


###
### START OF CODE TO DOWNLOAD IMAGES
###


baseurl = 'http://unwise.me/cutout_fits?version=allwise'  

# setting image size to 100 pixels for now
# may need to adjust this for larger galaxies
imsize = '100' # max size = 256 pixels



for id in nsaid:
    i = nsadict[int(id)]
    print 'NSA ID = ',id
    #print 'RA= ',nsa.RA[i]
    #print 'DEC= ',nsa.DEC[i]
    imurl = baseurl+'&ra=%.5f&dec=%.5f&size=%s&bands=%s'%(nsa.RA[i],nsa.DEC[i],imsize,args.band)
    print imurl
    wisetar = wget.download(imurl)  #trying to add in nsaid
    
    #print wisetar
    tartemp = tarfile.open(wisetar,mode='r:gz') #mode='r:gz'
    wnames = tartemp.getnames()
    wmembers = tartemp.getmembers()
    tartemp.extractall()
    for a in wnames:
        t = a.split('-')
        newname = 'NSA-'+str(id)+'-'+t[0]+'-'+t[2]+'-'+t[3]+'-'+t[4]
        os.rename(a, newname)
    # get rid of tar file
    os.remove(wisetar)

    





