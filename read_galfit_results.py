#!/usr/bin/env python

"""
GOAL:

- to read in results from galfit analysis of WISE 12 and 22um images

PROCEDURE:
- get r-band size measurements from NSA catalog
- get WISE 12mu size measurements
- get WISE 22um size measurement
- compare R(12) vs R(22)
- calculate the ratio of R(12)/R(r).  compare to LCS results
- calculate the ratio of R(22)/R(r).  compare to LCS results.
- are filament galaxies more like core or external samples?

USEAGE:

from within ipython on laptop

%run ~/github/virgowise/read_galfit_results.py --nsapath '/Users/rfinn/github/Virgo/tables/' 

then to generate plot of R22 vs R12:

plotresults() 

REQUIRED FILES:
- nsa fits catalog


NOTES:
"""




from astropy.io import fits
import argparse
import numpy as np
#from virgowise.common import *
import rungalfit
import sys
import os
import argparse

os.sys.path.append('/Users/rfinn/github/virgowise/')
from unwise_rungalfit import catalogs


mipspixelscale = 2.45
wisepixelscale = 2.75

#Need user to define galaxy image/sigma/psf path later on
parser = argparse.ArgumentParser(description ='Read galfit output and store results in fits table')
#parser.add_argument('--band', dest = 'band', default = '3', help = 'unWISE image band to download: 3=12um, 4=22um (can only do one at a time for now.')
parser.add_argument('--nsapath',dest = 'nsapath', default ='/Users/rfinn/github/Virgo/tables/', help = 'Location of NSA tables fits')
parser.add_argument('--makeplots',dest = 'makeplots', default =True, help = 'Plot results?  Default is True')
parser.add_argument('--band', dest = 'band', default = '3', help = 'unWISE image band to download: 3=12um, 4=22um (can only do one at a time for now.')
parser.add_argument('--virgopath',dest = 'virgopath', default ='/Users/rfinn/github/Virgo/', help = 'Location of NSA tables fits')
parser.add_argument('--nsafile',dest = 'nsafile', default ='virgo', help = 'nsa file to use.  could be nsa.virgo.fits (default) or use full to get nsa_v0_1_2.fits.  make sure nsapath points to the right place.')



args = parser.parse_args()

if args.makeplots:
    from matplotlib import pyplot as plt


# read in catalogs
cats = catalogs(args.nsapath,virgoflag=False)
#cats.define_sample()



ngc_filament_ids = [56403,56409,56410,56411,56434,56455,56456,56462,56469,56482,56489,61690,61691,61692,67593,88353,90371,93403,94217,104307,104439,118647,119230,119289,120018,120053,142509,143682,143686,143827,143951,143986,162674,163136,163783,164224,164358]


multiframe = np.array([ True, False, False, False, False, False, False, False, False,
       False, False, False, False, False, False, False,  True, False,
        True,  True, False, False, False, False, False, False, False,
       False, False, False, False, False, False, False, False, False,
       False,  True, False, False, False, False, False, False, False,
       False, False, False, False, False, False, False, False, False,
       False, False,  True, False, False, False, False, False, False,
       False, False, False, False,  True, False,  True, False, False,
       False, False, False, False, False, False,  True, False, False,
       False, False, False, False, False, False,  True, False, False,
       False, False,  True, False, False, False, False, False, False,
       False, False, False, False, False, False, False, False, False,
       False, False,  True, False, False,  True, False, False, False,
        True, False, False,  True, False, False, False,  True, False,
       False, False, False, False, False, False, False, False,  True,
       False, False, False, False,  True, False, False, False, False,
       False,  True,  True,  True,  True, False,  True,  True, False,
       False, False,  True, False, False, False, False, False, False,
       False, False, False,  True, False, False, False, False, False,
       False, False, False, False, False,  True, False,  True, False,
       False, False, False, False, False, False, False, False, False,
       False,  True,  True, False, False, False, False, False, False,
       False,  True, False,  True, False, False,  True, False, False,
       False, False, False, False, False,  True, False,  True, False,
       False, False, False, False, False, False,  True, False, False], dtype=bool)


lcsgals = np.array([ 70588,  70615,  70635,  70639,  70648,  70655,  70656,  70657,
        70658,  70696,  70703,  70704,  70706, 163566, 163584, 163589,
       170903, 171148,  43711,  43713,  43773,  43790,  43791,  43792,
        43796,  43817,  43830,  43836,  43851,  43872,  43873,  43875,
        69515,  69517,  69537,  69610,  69622,  69623,  69673, 146878,
       166866,  18152,  18161,  18232,  18234,  18235,  18236,  18266,
        18296,  18297,  18302, 165390,  68300,  68301,  68305,  68313,
        68330,  82064,  82068,  82076,  82168,  82170,  82172,  82182,
        82188,  82197,  82198,  82206,  82209,  82210,  98722, 166662,
       166699,  99058,  99060,  99478,  99535,  99538,  99540,  99644,
        99660,  99675,  99679,  99790,  99844,  99845,  99847,  99855,
        99862,  99865,  99882, 146606, 146613, 146636, 146659, 166658,
       166677, 166702, 113058, 113065, 113092, 113095, 113378, 113482,
       140153, 140156, 140160, 140187, 140197, 160491,  72561,  72623,
        72631,  72650,  72659,  72662,  72667,  72672,  72679,  72682,
        72690,  72706,  72733,  72738,  72749,  72750,  72767,  72778,
        72782,  72809,  72812,  79675,  79709,  79761,  79779, 146094,
       146115, 146119, 146121, 146127, 146130, 166165, 166167, 166181,
       166185,  72459,  72460,  72481,  72721,  79363,  79376,  79417,
        79551,  79558,  79561,  79564,  79565,  79591,  79605,  79607,
        79608,  79631,  79633,  79634,  79647,  79658,  79665,  79675,
        79680,  79706,  79709,  80763,  80769,  80873,  80882,  81018,
       145965, 145981, 145984, 145999, 146003, 146014, 166044, 166083,
       166090, 166096,  89063,  89078,  89101,  89108,  89129,  89134,
       103576, 103591, 103592, 103610, 103613, 103635, 103648, 103740,
       103742, 103791, 103816, 103832, 103889, 103899, 103903, 103917,
       103920, 103927, 103960, 103966, 104152, 104226, 104247, 104262,
       104285, 142655, 142668, 162769, 162792, 162798, 162814, 162838,
       162872], 'i')

#sample_ids = cats.nsa.NSAID[cats.sampleflag]

sample_ids = lcsgals[~multiframe]
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


def plotresults():
    plt.figure()
    plt.plot(r12,r22,'bo')
    plt.errorbar(r12,r22,xerr=r12_err,yerr=r22_err,fmt="None",ecolor='k')
    yline=np.arange(100)
    plt.plot(yline,yline,'k--')
    plt.xlabel('$R_{12} \ (pixels)$',fontsize=18)
    plt.ylabel('$R_{22} \ (pixels)$',fontsize=18)
    plt.title('Comparison of 12 and 22um Re from GALFIT')
    plt.gca().set_yscale('log')
    plt.gca().set_xscale('log')
    
    plt.axis([.2,80,.2,80])
    plt.subplots_adjust(bottom=.15)
    plt.savefig('r12_vs_r22.png')

def comparelcs():
    maxmag=10
    flag =  (m12 < maxmag) & (r12 < 20)
    plt.figure()
    plt.scatter(lcs_r24[flag],r12[flag],c=m12[flag],vmin=6,vmax=maxmag,s=80)#,cmap='jet')
    plt.colorbar()
    plt.errorbar(lcs_r24[flag],r12[flag],yerr=r12_err[flag],xerr=lcs_r24_err[flag],fmt="None",ecolor='k',alpha=0.5)
    xline = np.arange(0,30)
    plt.plot(xline,xline,'k--')
    plt.axis([-2,25,-2,25])
    plt.xlabel('Re(24) MIPS (arcsec)',fontsize=18)
    plt.ylabel('Re(12) WISE (arcsec)',fontsize=18)
    plt.savefig('MIPS24-WISE12.png')

# read in nsa catalog

#print 'this is what I think the NSA table filename is'
#print args.catalogpath+'nsa.virgo.fits'

nsa = fits.getdata(cats.nsatab)
nsadict=dict((a,b) for a,b in zip(nsa.NSAID,np.arange(len(nsa.NSAID))))

class lcs:
    def __init__(self,file):
        self.s = fits.getdata(infile)
        self.nsadict=dict((a,b) for a,b in zip(self.s.NSAID,np.arange(len(self.s.NSAID))))
        
class galaxy:
    def __init__(self, nsaid):
        self.nsaid = nsaid
        self.nsaindex = nsadict[int(nsaid)]
    def get_galfit_results(self,band='4',printflag = False):
        self.band = band
        #filename = 'NSA'+str(self.nsaid)-w+self.bands-1Comp-galfit-out.fits
        self.filename = 'NSA-'+str(self.nsaid)+'-unwise-'+'w'+self.band+'-1Comp-galfit-out.fits'
        #filename = 'NSA'+str(self.nsaid)-1Comp-galfit-out.fits
        # extension 2 is the model
        

        t = rungalfit.parse_galfit_1comp(self.filename)
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
    infile='/Users/rfinn/research/LocalClusters/NSAmastertables/LCS_all_size.fits'
    lcs = lcs(infile)
    nsa_re = np.zeros(len(sample_ids),'f')
    r12 = np.zeros(len(sample_ids),'f')
    r22 = np.zeros(len(sample_ids),'f')
    m12 = np.zeros(len(sample_ids),'f')
    r12_err = np.zeros(len(sample_ids),'f')
    r22_err = np.zeros(len(sample_ids),'f')
    lcs_r24 = np.zeros(len(sample_ids),'f')
    lcs_r24_err = np.zeros(len(sample_ids),'f')
    memb = np.zeros(len(sample_ids),'f')
    lcs_m24 = np.zeros(len(sample_ids),'f')
    for i in range(len(sample_ids)):
        lcs_r24[i] = lcs.s.fre1[lcs.nsadict[sample_ids[i]]]*mipspixelscale # in arcsec
        lcs_r24_err[i] = lcs.s.fre1err[lcs.nsadict[sample_ids[i]]]*mipspixelscale # in arcsec
        lcs_m24[i] = lcs.s.fmag1[lcs.nsadict[sample_ids[i]]]
        #memb[i] = lcs.membflag
        nsaid = sample_ids[i]
        g12 = galaxy(nsaid)

        #g.calc_sizeratio()
        nsa_re[i] = cats.nsa.SERSIC_TH50[g12.nsaindex]
        try:
            g12.get_galfit_results(band='3')
            r12[i] = g12.re*unwisepixelscale
            r12_err[i] = g12.re_err*unwisepixelscale
            m12[i] = g12.mag
        except IOError:
            print 'problem accessing band 3 data for NSA ',cats.nsa.NSAID[i]
            print 'filename = ',g12.filename
            print os.path.exists(g12.filename)

        g22 = galaxy(nsaid)
        try:
            g22.get_galfit_results(band='4')
            r22[i] = g22.re*unwisepixelscale
            r22_err[i] = g22.re_err*unwisepixelscale
        except IOError:
            print 'problem accessing band 4 data for NSA ',cats.nsa.NSAID[i]
            print 'filename = ',g22.filename
            print os.path.exists(g22.filename)
