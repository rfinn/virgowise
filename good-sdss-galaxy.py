from astroquery.sdss import SDSS
from astroquery.ned import Ned
from astropy import coordinates as coords
import astropy.units as u
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
import wget
from scipy.stats import scoreatpercentile
import os


galaxy_name = 'NGC5173'
result_table = Ned.query_object('NGC 5173')

#print result_table['RA(deg)'][0]


pos = coords.SkyCoord(ra=result_table['RA(deg)'][0]*u.deg,dec=result_table['DEC(deg)'][0]*u.deg, frame='icrs')
xid = SDSS.query_region(pos)
#print xid
sdsscoords = coords.SkyCoord(ra = xid['ra']*u.deg, dec=xid['dec']*u.deg,frame='icrs')

distance = pos.separation(sdsscoords)
match = (distance == min(distance))
print sdsscoords[match]

im = SDSS.get_images(matches = xid[match], band='g')

#Now retreieve image using wget:

run = (xid[match]['run'][0])
rerun = xid[match]['rerun'][0]
camcol = xid[match]['camcol'][0]
field = (xid[match]['field'][0])
band = 'r'

# taking a lot of this from core.py in astropy.sdss
baseurl = 'http://data.sdss3.org/sas'
IMAGING_URL_SUFFIX = ('{base}/dr{dr}/boss/photoObj/frames/'
                          '{rerun}/{run:04d}/{camcol}/'
                          'frame-{band}-{run:06d}-{camcol}-'
                          '{field:04d}.fits.bz2')
f = IMAGING_URL_SUFFIX.format(base=baseurl,run=run,dr='12',rerun=rerun,camcol=camcol,band=band,field=field)
print f

imfile = wget.download(f)
imdat = fits.getdata(imfile)
plt.imshow(imdat,origin='lower',vmin=scoreatpercentile(imdat,2),vmax=scoreatpercentile(imdat,98),cmap='gray')


#Get psf field file then run command in terminal to get image
# construct the url for the psField file
baseurl = 'http://data.sdss3.org/sas'
IMAGING_URL_SUFFIX = ('{base}/dr{dr}/boss/photo/redux/'
                          '{rerun}/{run:04d}/objcs/{camcol}/'
                          'psField-{run:06d}-{camcol}-'
                          '{field:04d}.fit')
psf = IMAGING_URL_SUFFIX.format(base=baseurl,run=run,dr='12',rerun=rerun,camcol=camcol,band=band,field=field)
print psf
#boss/photo/redux/%d/%d/objcs/%d/psField-%06d-%d-%04d.fit


# download the psField file
psfield_file = wget.download(psf)

# read in the psField
sdssbands = {'u':1,'g':2,'r':3,'i':4,'z':5}

row = 500
col = 500
band = 'r'

t=os.getcwd()
psf_image = t+'/'+galaxy_name+'-psf.fits'
print 'Run the following command from the linux prompt: \n'
print 'read_PSF '+psfield_file+' '+str(sdssbands[band]-1) +' '+str(row)+' '+str(col)+psf_image
#os.system('read_PSF '+psfield_file+' '+str(sdssbands[band]-1) +' '+str(row)+' '+str(col)+psf_image)










#BELOW IS THE OLD STUFF

#NOW GET PSF IMAGE SO THAT WE CAN USE THIS FOR GALFIT
# construct the url for the psField file
#baseurl = 'http://data.sdss3.org/sas'
#IMAGING_URL_SUFFIX = ('{base}/dr{dr}/boss/photo/redux/'
                         # '{rerun}/{run:04d}/objcs/{camcol}/'
                         # 'psField-{run:06d}-{camcol}-'
                         # '{field:04d}.fit')
#psf = IMAGING_URL_SUFFIX.format(base=baseurl,run=run,dr='12',rerun=rerun,camcol=camcol,band=band,field=field)
#print psf


# download the psField file
#psfield_file = wget.download(psf)

# read in the psField
#sdssbands = {'u':1,'g':2,'r':3,'i':4,'z':5}
# recall u,g,r,i,z == 0,1,2,3,4 so just add one to the index
#pstruct = fits.getdata(psfield_file,ext=sdssbands[band])





