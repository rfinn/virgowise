#!/usr/bin/env python

import numpy as np
import numpy.ma as ma
import os
from astropy.io import fits
from astropy.io import ascii
from astropy import constants as c
from astropy import units as u
from astropy.table import Table, join, hstack, Column, MaskedColumn
from astropy.coordinates import SkyCoord
import sys
from astropy.cosmology import WMAP9 as cosmo

homedir = os.getenv('HOME')
sys.path.append(homedir+'/github/APPSS/')
from join_catalogs import make_new_cats, join_cats

zmin = 0.002
zmax = 0.025

a100 = fits.getdata(homedir+'/research/APPSS/tables/a100-sdss-wise.fits')
# cut on velocity
flag = a100['Vhelio'] < zmax*3e5
a100 = a100[flag]


### READ IN NSA
nsa = fits.getdata(homedir+'/research/NSA/nsa_v0_1_2.fits')

wise = fits.getdata(homedir+'/research/NSA/nsa_v0_1_2_wise.fits')
flag = nsa['ZDIST'] < zmax
nsa = Table(nsa[flag])
wise = Table(wise[flag])
wise = wise[wise.colnames[0:52]]
nsa = hstack([nsa,wise])
# this will have two RA, DEC columns
# rename first set to RA,DEC
nsa.rename_column('RA_1','RA')
nsa.rename_column('DEC_1','DEC')
voffset=300
a1002, a100_matchflag, nsa2, nsa_matchflag = make_new_cats(a100, nsa,RAkey1='RAdeg_Use',DECkey1='DECdeg_Use', velocity1=a100['Vhelio'], velocity2=nsa['ZDIST']*3e5, maxveloffset = voffset,maxoffset=15.)

# join a100-sdss and nsa into one table
joined_table = hstack([a1002,nsa2])

# print match statistics for full catalogs
print('AFTER MATCHING')
print('total number in A100 = ',sum(a100_matchflag))
print('total number in NSA = ',sum(nsa_matchflag))
print('number of unique galaxies = ',len(a1002))
print('number of matches between A100 and NSA = ',sum(a100_matchflag & nsa_matchflag))
print('number in A100 but not in NSA = ',sum(a100_matchflag & ~nsa_matchflag))
print('number in NSA but not in A100 = ',sum(~a100_matchflag & nsa_matchflag))

joined_table.write('a100-nsa-for-adap.fits',overwrite=True)
