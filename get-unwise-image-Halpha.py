#! /usr/bin/env python

"""
GOAL:
- to retrieve unWISE images of galaxies that we imaged in Halpha

PROCEDURE:
- query Halpha directory
- get list of NSAIDs
- read in NSA

USEAGE:

from within ipython on laptop

 %run ~/github/virgowise/get-unwise-image-Halpha.py --nsapath '/Users/rfinn/research/NSA/'





REQUIRED FILES:
- nsa fits catalog


NOTES:
I want to download the unWISE images with the same NGC number as the Halpha images. I can do this by splitting up the radial output filenames and indexing just the nsatlas number of the galaxy,gonna need a way to convert that to the number in NED but for now: 

-1/28/18 I want to somehow download the tar files from unWISE but know what galaxy is which by outputting the nsatlas numbers which is what our Halpha images are in. Also commenting details of code.


"""

import sys
import os
from astropy.io import fits
import numpy as np
from astroquery.sdss import SDSS
from astroquery.ned import Ned
from astropy import coordinates as coords
import astropy.units as u
from matplotlib import pyplot as plt
from astropy.io import fits
import wget
from scipy.stats import scoreatpercentile
import tarfile
import argparse
import glob
import gzip

parser = argparse.ArgumentParser(description ='Download only specific unWISE images that have the same ID number as the Halpha images')
parser.add_argument('--Halphaimagepath', dest = 'path', default = '/home/share/research/Virgo/reduced/NSAIDformatching/', help = 'full path to Halpha images')

parser.add_argument('--nsapath', dest = 'nsapath', default = '/home/share/catalogs/', help = 'full path to Halpha images')
 #ex. args.path; this is only for the specific directory in COMA to grab nsaid
parser.add_argument('--nsaid', dest = 'nsaid', default = '156774', help = 'Enter in nsaid')

args = parser.parse_args() #brings in these arguments above


#search_string = args.path + '/*-CS-snapshot.png' #allows me to index nsaid from Halpha images 
#print search_string
#input_images = glob.glob(search_string)
#nsaid=[]

#for f in input_images: #this physically index's the nsaid
#    t = f.split('/')
#    junk = t[-1].split('-')
#    nsaid.append(junk[2])
#    print nsaid

#nsaid = ['56462','67595','164358','54578','61690','61692']
nsaid = args.nsaid

#Now read in the nsa fits table to go the ra and dec


nsatab = args.nsapath+'nsa_v0_1_2.fits'
nsa = fits.getdata(nsatab)
nsadict = dict((a,b) for a,b in zip(nsa.NSAID,np.arange(len(nsa.NSAID))))

#baseurl = 'http://unwise.me/cutout_fits?version=allwise'
baseurl = 'http://unwise.me/cutout_fits?version=allwise'  
wisefilter1 = '3'
#wisefilter2 = '4'
imsize = '100' # max size = 256 pixels
bands='3'


for id in nsaid:
#for j in range(1):
    #id = nsaid[j]
    i = nsadict[int(id)]
    print i
    print 'RA= ',nsa.RA[i]
    print 'DEC= ',nsa.DEC[i]
    ra = nsa.RA[i]
    dec = nsa.DEC[i]
    imurl = baseurl+'&ra=%.5f&dec=%.5f&size=%s&bands=%s'%(ra,dec,imsize,bands)
    print imurl
    wisetar = wget.download(imurl)  #trying to add in nsaid
    
    
    tartemp = tarfile.open(wisetar,mode='r:gz') #mode='r:gz'
    #tartemp = tartemp.add(id)
    wnames = tartemp.getnames()
    #wnames.append(nsaid)
    wmembers = tartemp.getmembers()
    #path = '/home/share/research/Virgo/Galfit2018'
    tartemp.extractall()
    for a in wmembers:
        os.rename(a.name, str(id)+'-'+a.name)    

    #rename = []

    
    
    for i in wnames:
        split = i.split('-')
        print split
    
    #os.rename(old,new) path directories though :(
    #for i in wmembers:
    #    if wmembers[4]:
    #        print wmembers[4]
    #        #tartemp.extract(wmembers[4])
    #    if wmembers[8]:
    #        print wmembers[8]
    #        #tartemp.extract(wmembers[8])
    

'''
baseurl = 'http://unwise.me/cutout_fits?version=allwise'
wisefilter1 = '3' # 12 micron
wisefilter2 = '4' # 22 micron
imsize = '100' # max size = 256 pixels
bands='34'
ra = RA[id]
dec = DEC[id]
version=neo1ra=41&dec=10&size=100&bands=12
imurl = baseurl+'&ra=%.5f&dec=%.5f&size=%s&bands=%s'%(ra,dec,imsize,bands)


wisetar = wget.download(imurl)
tartemp = tarfile.open(wisetar,mode='r:gz') #mode='r:gz'


wnames = tartemp.getnames()
wmembers = tartemp.getmembers()

tartemp.extract(wmembers[0])
wisedat = fits.getdata(filelist[0])
'''



"""

#THIS IS EXACTLY DR.FINN'S CODE BUT HER'S IS IN .PYNB FORMAT SO I'M JUST GONNA COPY AND PASTE IT HERE :)))))))



galaxy_name = 'NGC5320' 
result_table = Ned.query_object(galaxy_name)

pos = coords.SkyCoord(ra=result_table['RA(deg)'][0]*u.deg,dec=result_table['DEC(deg)'][0]*u.deg, frame='icrs')



ra = pos.ra.deg
dec = pos.dec.deg
wisefilter = '3' # 12 micron
imsize = '100' # max size = 256 pixels
bands='34'
version=neo1&ra=41&dec=10&size=100&bands=12
imurl = baseurl+'&ra=%.5f&dec=%.5f&size=%s&bands=%s'%(ra,dec,imsize,bands)


#this will download a tar file
wisetar = wget.download(imurl)
tartemp = tarfile.open(wisetar,mode='r:gz') #mode='r:gz'


wnames = tartemp.getnames()
wmembers = tartemp.getmembers()

   #addnsaid=[]                                                                

   # input_files = glob.glob(wisetar)                                           
   # for i in input_files:                                                      
   #     n = i.split('.')                                                       
  #      addnsaid.append(id[3]) #I know I'm doing this part wrong, there's a ta\
r.add but it doesn't work for our type of file no matter what mod it's in. I wa\
nt to separate the file name, which the n does above, but now I want to add the\
 string, 'id' in, which is part of this loop, before the .tar.gz.

tartemp.extract(wmembers[0])
wisedat = fits.getdata(filelist[0])

"""
