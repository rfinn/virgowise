#!/usr/bin/env python



# get r-band size measurements from NSA catalog

# get WISE 12mu size measurements

# get WISE 22um size measurement

# compare R(12) vs R(22)

# calculate the ratio of R(12)/R(r).  compare to LCS results

# calculate the ratio of R(22)/R(r).  compare to LCS results.

# are filament galaxies more like core or external samples?



from astropy.io import fits
import argparse
import numpy as np
#from virgowise.common import *
import rungalfit
import sys
import argparse

#Need user to define galaxy image/sigma/psf path later on
parser = argparse.ArgumentParser(description ='Read galfit output and store results in fits table')
#parser.add_argument('--band', dest = 'band', default = '3', help = 'unWISE image band to download: 3=12um, 4=22um (can only do one at a time for now.')
parser.add_argument('--nsapath',dest = 'nsapath', default ='~/github/Virgo/tables/', help = 'Location of NSA tables fits')

args = parser.parse_args()


os.sys.path.append('~/github/virgowise/')
from unwise-rungalfit import catalogs
cats = catalogs.catalogs(args.nsapath)
cats.define_sample()



ngc_filament_ids = [56403,56409,56410,56411,56434,56455,56456,56462,56469,56482,56489,61690,61691,61692,67593,88353,90371,93403,94217,104307,104439,118647,119230,119289,120018,120053,142509,143682,143686,143827,143951,143986,162674,163136,163783,164224,164358]

sample_ids = cats.nsa.NSAID[cats.sampleflag]

#removing the following galaxies
# 119303 - nothing there in unWISE image - could be coordinates or not detected

unwisepixelscale = 2.75 #arcsec/pixel, unWISE
# list of NGC filament galaxies (starting with those with Halpha)

def print_galfit_results(parse_output):            
    header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_SKY','ERROR','CHI2NU']
    t = parse_output
    for i in range(len(header_keywords)):
        try:
            print '%6s : %5.2f +/- %5.2f'%(header_keywords[i],t[i][0],t[i][1])
        except:
            print '%6s : %5.2f'%(header_keywords[i],t[i])


parser = argparse.ArgumentParser(description ='get galfit results for a sample of NGC galaxies')
parser.add_argument('--catalogpath', dest = 'catalogpath', default = '~/github/Virgo/tables/', help = 'path to NSA catalog (nsa.virgo.fits)')

args = parser.parse_args()


# read in nsa catalog

#print 'this is what I think the NSA table filename is'
#print args.catalogpath+'nsa.virgo.fits'

nsa = fits.getdata(args.catalogpath+'nsa.virgo.fits')
nsadict=dict((a,b) for a,b in zip(nsa.NSAID,np.arange(len(nsa.NSAID))))



class galaxy:
    def __init__(self, nsaid):
        self.nsaid = nsaid
        self.nsaindex = nsadict[int(nsaid)]
    def get_galfit_results(self,band='4',printflag = False):
        self.band = band
        #filename = 'NSA'+str(self.nsaid)-w+self.bands-1Comp-galfit-out.fits
        filename = 'NSA'+str(self.nsaid)+'-'+'w'+self.band+'-1Comp-galfit-out.fits'
        #filename = 'NSA'+str(self.nsaid)-1Comp-galfit-out.fits
        # extension 2 is the model
        

        t = rungalfit.parse_galfit_1comp(filename)
        if printflag:
            print_galfit_results(t)
        
        header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_SKY','ERROR','CHI2NU']
        self.xc, self.xc_err = t[0]
        self.yc, self.yc_err = t[1]
        self.mag, self.mag_err = t[2]
        self.re, self.re_err = t[3]
        self.sersicn, self.sersicn_err = t[4]
        self.BA, self.BA_err = t[5]
        self.PA, self.PA_err = t[6]
        self.sky, self.sky_err = t[7]
        self.error = t[8]
        self.chi2nu = t[9]
    def calc_sizeratio(self):
        self.sizeratio = self.re*unwisepixelscale/nsa.SERSIC_TH50[self.nsaindex]
        print 'NSAID %6i-w%1i:  R12, Re, R12/Re = %8.2f %8.2f %8.2f'%(int(self.nsaid),int(self.band),self.re*unwisepixelscale,nsa.SERSIC_TH50[self.nsaindex],self.sizeratio)


if __name__ == "__main__":

    nsa_re = np.zeros(len(sample_ids),'f')
    r12 = np.zeros(len(sample_ids),'f')
    r22 = np.zeros(len(sample_ids),'f')
    r12_err = np.zeros(len(sample_ids),'f')
    r22_err = np.zeros(len(sample_ids),'f')
    for i in range(len(sample_ids)):
        nsaid = sample_ids[i]
        g = galaxy(nsaid)

        #g.calc_sizeratio()
        nsa_re[i] = cats.nsa.SERSIC_TH50[g.nsaindex]
        g.get_galfit_results(band='3')
        r12[i] = g.re
        r12_err[i] = g.re_err
        try:
            g.get_galfit_results(band='4')
            r22[i] = g.re
            r22_err[i] = g.re_err
        except:
            print 'problem accessing band 4 data for NSA ',cats.nsa.NSAID[i]