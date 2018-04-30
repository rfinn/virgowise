#!/usr/bin/env python



# get r-band size measurements from NSA catalog

# get WISE 12mu size measurements

# get WISE 22um size measurement

# compare R(12) vs R(22)

# calculate the ratio of R(12)/R(r).  compare to LCS results

# calculate the ratio of R(22)/R(r).  compare to LCS results.

# are filament galaxies more like core or external samples?



from astropy.io import fits
import argsparse

from virgowise.common import *
import rungalfit

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
parser.add_argument('--catalogpath', dest = 'catalogpath', default = '/Users/rfinn/Research/NSA', help = 'path to NSA catalog (nsa_v0_1_2.fits)')

args = parser.parse_args()


# read in nsa catalog

nsa = fits.getdata(args.catalogpath+'nsa_v0_1_2.fits')
nsadict=dict((a,b) for a,b in zip(nsa.NSAID,arange(len(nsa.NSAID))))

for nsa in ngc_filament_ids:
    filename = 'NSA'+str(nsaid)+'-1Comp-galfit-out.fits' # extension 2 is the model

    t = rungalfit.parse_galfit_1comp(filename)
    print_galfit_results(t)


class galaxy:
    def __init__(self, nsaid):
        self.nsaid = nsaid
        self.nsaindex = nsadict(int(nsaid))
    def get_galfit_results(self):
        filename = 'NSA'+str(nsaid)+'-1Comp-galfit-out.fits' # extension 2 is the model

        t = rungalfit.parse_galfit_1comp(filename)
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


