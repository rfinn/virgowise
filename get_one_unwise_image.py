#THIS IS EXACTLY DR.FINN'S CODE BUT HER'S IS IN .PYNB FORMAT SO I'M JUST GONNA COPY AND PASTE IT HERE :)))))))

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


galaxy_name = 'NGC5012'
result_table = Ned.query_object(galaxy_name)

pos = coords.SkyCoord(ra=result_table['RA(deg)'][0]*u.deg,dec=result_table['DEC(deg)'][0]*u.deg, frame='icrs')


baseurl = 'http://unwise.me/cutout_fits?version=allwise'
ra = pos.ra.deg
dec = pos.dec.deg
wisefilter = '3' # 12 micron
imsize = '100' # max size = 256 pixels
bands='34'
#version=neo1&ra=41&dec=10&size=100&bands=12
imurl = baseurl+'&ra=%.5f&dec=%.5f&size=%s&bands=%s'%(ra,dec,imsize,bands)


# this will download a tar file
wisetar = wget.download(imurl)
tartemp = tarfile.open(wisetar,mode='r:gz') #mode='r:gz'


wnames = tartemp.getnames()
wmembers = tartemp.getmembers()

#tartemp.extract(wmembers[0])
#wisedat = fits.getdata(filelist[0])






