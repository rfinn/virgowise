#!/usr/bin/env python

from matplotlib import pyplot as plt
from PIL import Image
from scipy.stats import scoreatpercentile


from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.visualization import simple_norm
from astropy import units as u
from astropy.nddata import Cutout2D


nuv = 'galex/VFID0566-NGC5989-NGC5989-nuv-158.fits'
leg = 'legacy/VFID0566-NGC5989-legacy-158.jpg'
legr = 'legacy/VFID0566-NGC5989-legacy-158-r.jpg'
w1 = 'unwise/VFID0566-NGC5989-unwise-2361p590-w1-img-m.fits'
w2 = 'unwise/VFID0566-NGC5989-unwise-2361p590-w2-img-m.fits'
w3 = 'unwise/VFID0566-NGC5989-unwise-2361p590-w3-img-m.fits'
w4 = 'unwise/VFID0566-NGC5989-unwise-2361p590-w4-img-m.fits'


