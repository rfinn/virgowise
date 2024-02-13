#!/usr/bin/env python
import os
from matplotlib import pyplot as plt
import numpy as np

from astropy.table import Table
from astropy.io import fits
import warnings
warnings.filterwarnings("ignore")
homedir = os.getenv("HOME")
mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']
os.sys.path.append('/home/rfinn/github/LCS/python/Python3/')

plotdir = homedir+'/proposals/NASA2020/figures/'
plotdir = homedir+'/proposals/NSF2020/figures/'
from LCScommon import *
# using colors from matplotlib default color cycle
mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']

zmin = 0.0017
zmax = 0.025
sersicmin = 5.5
snrcut = 20

zmin = 0.002
zmax = 0.025
sersicmin = 5.5
snrcut = 20

class catalogs():
    def __init__(self,catalog_path,virgoflag = False):
        # read in nsa, wise, co catalogs
        if virgoflag:
            self.nsatab = catalog_path + 'nsa.virgo.fits'
            self.wisetab = catalog_path + 'nsa_wise.virgo.fits'
            self.cotab = catalog_path + 'nsa_CO-Gianluca.virgo.fits'
            self.co = fits.getdata(self.cotab)
        else:
            print('virgo is not true')
            self.nsatab = catalog_path + 'nsa_v0_1_2.fits'
            self.wisetab = catalog_path + 'nsa_v0_1_2_wise.fits'
            self.jmasstab = catalog_path + 'nsa_v1_2_fsps_v2.4_miles_chab_charlot_sfhgrid01.fits'
        self.nsa = fits.getdata(self.nsatab)
        self.wise = fits.getdata(self.wisetab)
        self.jmass = fits.getdata(self.jmasstab)
        self.a100nsa = Table.read(homedir+'/research/wisesize/a100-nsa-for-adap.fits')
        #zflag = (self.a100nsa['Vhelio'] < zmax*3e5) | (self.a100nsa['ZDIST'] < zmax)
        #self.a100nsa = self.a100nsa[zflag]
        self.nsadict=dict((a,b) for a,b in zip(self.nsa.NSAID,np.arange(len(self.nsa.NSAID)))) #useful for us can easily look up galaxy ID's       
        self.get_tempel()
        self.gmi = self.a100nsa['ABSMAG'][:,3] - self.a100nsa['ABSMAG'][:,5]
        absMag_i = self.a100nsa['ABSMAG'][:,5] - 5*np.log10(.7) # correct from HO=100 to 70
        logMstarTaylor=1.15+0.70*(self.gmi) -0.4*(absMag_i)
        #self.a100Flag = np.zeros(len(self.a100nsa),'bool')
        #self.nsaFlag = np.zeros(len(self.a100nsa),'bool')        
        #self.a100Flag = ~self.a100nsa['AGC'].mask
        #self.nsaFlag = ~self.a100nsa['NSAID'].mask
        self.a100Flag = self.a100nsa['AGC'] > 0
        self.nsaFlag = self.a100nsa['NSAID'] > 0        
        self.logMstar = self.a100nsa['logMstarTaylor']*self.a100Flag + logMstarTaylor*(~self.a100Flag & self.nsaFlag)
        self.color = self.a100nsa['gmi_corr']*self.a100Flag + self.gmi*(~self.a100Flag)
        self.calcagn()
        self.calcagn2()        
        self.define_sample()
        self.define_sample2()        
  
    def calcagn(self):
        self.AGNKAUFF= ((np.log10(self.nsa.O3FLUX/self.nsa.HBFLUX) > (.61/(np.log10(self.nsa.N2FLUX/self.nsa.HAFLUX)-.05)+1.3)) | (np.log10(self.nsa.N2FLUX/self.nsa.HAFLUX) > 0.))
#y=(.61/(x-.47)+1.19)
        self.AGNKEWLEY= ((np.log10(self.nsa.O3FLUX/self.nsa.HBFLUX) > (.61/(np.log10(self.nsa.N2FLUX/self.nsa.HAFLUX)-.47)+1.19)) | (np.log10(self.nsa.N2FLUX/self.nsa.HAFLUX) > 0.3))
        self.agnflag = self.AGNKAUFF
        #self.agnflag = self.AGNKEWLEY
    def calcagn2(self):
        self.AGNKAUFF2 = np.zeros(len(self.a100nsa),'bool')
        # this is true if the object in an AGN
        t= ((np.log10(self.a100nsa['O3FLUX']/self.a100nsa['HBFLUX']) > (.61/(np.log10(self.a100nsa['N2FLUX']/self.a100nsa['HAFLUX'])-.05)+1.3)) | (np.log10(self.a100nsa['N2FLUX']/self.a100nsa['HAFLUX']) > 0.))
        # in 2022, getting error about mask,
        # AttributeError: 'numpy.ndarray' object has no attribute 'mask'
        
        self.AGNKAUFF2 = t
        #self.AGNKAUFF2[~t.mask] = t[~t.mask].data
        self.t2 = t

        #y=(.61/(x-.47)+1.19)
        self.AGNKEWLEY2 = np.zeros(len(self.a100nsa),'bool')
        t= ((np.log10(self.a100nsa['O3FLUX']/self.a100nsa['HBFLUX']) > (.61/(np.log10(self.a100nsa['N2FLUX']/self.a100nsa['HAFLUX'])-.47)+1.19)) | (np.log10(self.a100nsa['N2FLUX']/self.a100nsa['HAFLUX']) > 0.3))
        self.t1 = t
        # in 2022, getting error about mask
        #self.AGNKEWLEY2[~t.mask] = t[~t.mask].data
        self.AGNKEWLEY2 = t        
        #self.agnflag2 = self.AGNKAUFF2
        self.agnflag2 = np.zeros(len(self.a100nsa),'bool')
        #self.agnflag2[~self.t1.mask] = self.t1[~self.t1.mask].data
        self.agnflag2 = self.t1
        
        # W1 - W2 > 0.8
        self.wiseagn = np.zeros(len(self.a100nsa),'bool')
        t = ((self.a100nsa['W1MPRO'] - self.a100nsa['W2MPRO']) > 0.8 )| \
            ((self.a100nsa['w1_mag'] - self.a100nsa['w2_mag']) > 0.8)
        
        self.wiseagn = t
        self.agn = self.agnflag2 | self.wiseagn
    def define_sample(self):
        self.w3_flag = self.wise.W3SNR>10 
        self.w4_flag = self.wise.W4SNR>5      
        #self.co_flag = self.co.COdetected == '1'
        #self.sampleflag = self.w3_flag & self.w4_flag & self.co_flag
        
        self.sampleflag = (self.wise.W3MPRO < 10.) & (self.nsa.SERSIC_TH50 > sersicmin) \
        & (self.nsa.ZDIST > zmin)& (self.nsa.ZDIST < zmax) & ~self.agnflag
        print('number of galaxies in sample = ',sum(self.sampleflag))
        self.allbutsizeflag = (self.wise.W3MPRO < 10.)  \
        & (self.nsa.ZDIST > zmin)& (self.nsa.ZDIST < zmax) & ~self.agnflag
        print('number of galaxies in sample = ',sum(self.sampleflag))
        #print('number of wise AGN = ',sum(self.wiseagn))
    def define_sample2(self):

        self.w3_flag = (self.a100nsa['w3_nanomaggies']*np.sqrt(self.a100nsa['w3_nanomaggies_ivar']) > snrcut) | (self.a100nsa['W3SNR'] > snrcut)
        #self.w3_flag = self.w3_snr > 10
        self.nuv_snr = self.a100nsa['NMGY'][:,1]*np.sqrt(self.a100nsa['NMGY_IVAR'][:,1])
        self.nuv_flag = (self.nuv_snr > snrcut)
        self.snr_flag = self.w3_flag | self.nuv_flag
                         
        self.w4_snr = (self.a100nsa['w4_nanomaggies']*np.sqrt(self.a100nsa['w4_nanomaggies_ivar']))
        self.w4_flag = self.w4_snr > snrcut
        self.sampleflag2 = self.snr_flag &\
             ((self.a100nsa['SERSIC_TH50'] > sersicmin) | (self.a100nsa['petroR50_r'] > sersicmin)) &\
             (((self.a100nsa['ZDIST'] > zmin)& (self.a100nsa['ZDIST'] < zmax)) | \
              ((self.a100nsa['Vhelio'] > zmin*3e5) & (self.a100nsa['Vhelio'] < zmax*3e5)))
        self.sizeflag2 = ((self.a100nsa['SERSIC_TH50'] > sersicmin) | (self.a100nsa['petroR50_r'] > sersicmin))        
        self.zflag2 = \
             (((self.a100nsa['ZDIST'] > zmin)& (self.a100nsa['ZDIST'] < zmax)) | \
              ((self.a100nsa['Vhelio'] > zmin*3e5) & (self.a100nsa['Vhelio'] < zmax*3e5))) 
             #& self.snr_flag 
        self.allbutsizeflag2 = \
             (((self.a100nsa['ZDIST'] > zmin)& (self.a100nsa['ZDIST'] < zmax)) | \
              ((self.a100nsa['Vhelio'] > zmin*3e5) & (self.a100nsa['Vhelio'] < zmax*3e5))) \
             & self.snr_flag 
        self.r50 = self.a100nsa['SERSIC_TH50']*self.nsaFlag + self.a100nsa['petroR50_r']*(~self.nsaFlag & self.a100Flag)
        self.r90 = 1.58*(self.a100nsa['PETROTH90']*self.nsaFlag + self.a100nsa['petroR90_r']*(~self.nsaFlag & self.a100Flag))
        self.D90 = 2*self.r90
        print('parent sample (zcut) = ',sum(self.allbutsizeflag2))
        print('resolved sample (not AGN) = ',sum(self.sampleflag2 & (~self.agnflag2)))
        print('resolved sample (not AGN) = %.3f'%(sum(self.sampleflag2 & (~self.agnflag2))/sum(self.allbutsizeflag2)))

        print('resolved sample (AGN) = ',sum(self.sampleflag2 & self.agnflag2))
        print('resolved sample (AGN) = %.3f'%(sum(self.sampleflag2 & (self.agnflag2))/sum(self.allbutsizeflag2)))
        print('resolved sample wise AGN = ',sum(self.sampleflag2 & self.wiseagn))
        print('resolved sample wise AGN = %.3f'%(sum(self.sampleflag2 & (self.wiseagn))/sum(self.allbutsizeflag2)))
        #self.sampleflag2 = self.sampleflag2 & ~self.agnflag2
    def write_sample_table(self):
        outfile = '/home/rfinn/research/wisesize/nsa-wisesize-catalog.fits'
        fits.writeto(outfile,self.nsa[self.sampleflag],overwrite=True)
        outfile = '/home/rfinn/research/wisesize/nsa-wisesize-catalogv2.fits'
        self.a100nsa.write(outfile,overwrite=True)
    def get_tempel(self):
        cat1 = homedir+'/research/wisesize/nsa_wisesize_catalog_Tempelfilaments.fits'
        tempel_filaments = Table.read(cat1)
        cat1 = homedir+'/research/wisesize/nsa_wisesize_catalog_Tempelgroups.fits'
        tempel_groups = Table.read(cat1)
        log_groups_mass_tempel = np.log10(tempel_groups['MNFW'])+12.
        ngal_with_env = (tempel_groups['MNFW']> 0) | (tempel_filaments['Separation'] > 0)
        print('number of galaxies with tempel env data = ',sum(ngal_with_env))
                                                         
        self.group_flag = (log_groups_mass_tempel>=13) & (log_groups_mass_tempel<14) #& (tempel_filaments['Ngal']<5)
        self.cluster_flag = (log_groups_mass_tempel>=14.) #& (log_groups_mass_tempel<15.5)

        self.field_flag = ((tempel_filaments['Separation']>0) & (tempel_filaments['Ngal']<2) & (tempel_groups['Separation']>0) & (np.isnan(log_groups_mass_tempel)) | (log_groups_mass_tempel<12))
        self.poorgroup_flag = (log_groups_mass_tempel>=12) & (log_groups_mass_tempel<13) & (tempel_filaments['Ngal']<5)

        self.filament_flag = (tempel_filaments['Ngal']>5) & (log_groups_mass_tempel<13.) #& (tempel_filaments['Ngal']<5)

        
    def plot_color_mstar(self):
        fsize = 18
        plt.figure(figsize=(12,4))
        plt.subplots_adjust(left=.1,bottom=.15,wspace=.2,right=.95)
        NUVr=self.nsa.ABSMAG[:,1] - self.nsa.ABSMAG[:,4]
        stellarmass = self.jmass.MSTAR_50
        plt.subplot(1,2,1)
        plt.hexbin(stellarmass[pflag],NUVr[pflag],cmap='gray_r',extent=(6,12.5,0,7))
        plt.axis([6,12.5,0,7])
        plt.xlabel(r'$ log_{10}\left( M_\star /M_\odot \right) $',fontsize=fsize)
        plt.ylabel(r'$ NUV - r $',fontsize=fsize)
        plt.title('This Proposal',fontsize=fsize)
        s = 'N = %i'%(sum(pflag))
        plt.text(6.1,6.1,s,fontsize=fsize-2)
        plt.subplot(1,2,2)
        plt.hexbin(stellarmass,NUVr,cmap='gray_r',extent=(6,12.5,0,7))
        plt.axis([6,12.5,0,7])
        plt.xlabel(r'$ log_{10}\left( M_\star /M_\odot \right) $',fontsize=fsize)
        plt.title('All NASA-Sloan Atlas',fontsize=fsize)
        allflag = (NUVr > 0) & (NUVr < 7) & (stellarmass > 6) & (stellarmass < 12.5) & ~self.agnflag
        s = 'N = %i'%(sum(allflag))
        plt.text(6.1,6.1,s,fontsize=fsize-2)
        #X,Y =  np.meshgrid(stellarmass[allflag],NUVr[allflag])
        plt.savefig(homedir+'proposals/NASA2018/sampleselection.pdf')
        #plt.ylabel(r'$ NUV - r $')
    def plot_color_mstar2(self,plotsingle=True):
        fsize = 18
        xmin,xmax=8,11.5
        ymin,ymax=0,1.7
        if plotsingle:
            #plt.figure(figsize=(12,4))
            plt.figure(figsize=(6,4))
            plt.subplots_adjust(left=.1,bottom=.15,wspace=.2,right=.95)            
        myflags = [self.sampleflag2,np.ones(len(self.sampleflag2),'bool')]# & self.agnflag2]

        NUVr=self.color
        stellarmass = self.logMstar
        allax=[]
        titles = ['This Proposal','A100+NSA (z<0.025)']
        titles = ['(a) Resolved Sample','A100+NSA']
        #for i in range(len(myflags)):
        for i in range(1):            
            #plt.subplot(1,2,i+1)
            plt.hexbin(stellarmass[ myflags[i]],NUVr[myflags[i]],cmap='gray_r',extent=(xmin,xmax,ymin,ymax),gridsize=50)
            #plt.plot(stellarmass[ myflags[0]],NUVr[myflags[0]],'k.',c=mycolors[1],alpha=.1,label='Resolved')            
            plt.axis([xmin,xmax,ymin,ymax])
            plt.xlabel(r'$ log_{10}\left( M_\star /M_\odot \right) $',fontsize=fsize)
            plt.ylabel(r'$ g - i $',fontsize=fsize)
            s = 'N = %i'%(sum(myflags[i]))
            plt.text(0.05,.9,s,transform=plt.gca().transAxes,fontsize=fsize-2)
            #X,Y =  np.meshgrid(stellarmass[allflag],NUVr[allflag])
            #plt.colorbar()
            allax.append(plt.gca())
            plt.title(titles[i],fontsize=fsize)
            xl = np.linspace(8.25,11,10)
            dmag = 0.15
            yl = .1*(xl-10)+1.15
            plt.plot(xl,yl,'r-')
            #plt.plot(xl,yl+dmag,'r--')
            #plt.plot(xl,yl-dmag,'r--')
            #plt.fill_between(xl,y1=yl+dmag,y2=yl-dmag,color='r',alpha=.08)
        #cb = plt.colorbar(fraction=.08,ax=allax)
        #cb.set_label('$N_{gal}$',fontsize=fsize)
        #plt.legend()
        if plotsingle:
            plt.savefig(homedir+'proposals/NASA2018/color-mass-adap202.pdf')
        #plt.ylabel(r'$ NUV - r $')
    def plot_positions_lcs(self):
        # plot positions
        plt.figure(figsize=(8,5))
        plt.subplots_adjust(bottom=.15,top=.95,left=.1,right=.95)
        plt.scatter(self.nsa.RA[self.sampleflag],self.nsa.DEC[self.sampleflag],s=10,alpha=.25,marker='.')#,label='This proposal')#,color='0.5',alpha='0.5')
        plt.legend(['This proposal'],loc='upper left')#,'LCS'])
        for i in range(len(clusternames)):
            s=' '+clusternames[i]
            text(clusterRA[clusternames[i]],clusterDec[clusternames[i]],s,fontsize=12,color='r')
            drawbox(cluster24Box[clusternames[i]],'r-')
        plt.xlim(100,275)
        plt.ylim(-5,70)
        plt.xlabel('$RA \ (deg) $',fontsize=16)
        plt.ylabel('$DEC \ (deg) $',fontsize=16)
        plt.savefig('/home/rfinn/proposals/NASA2020/WISEsize-LCS.pdf')
    def compare_size_resolution(self,plotsingle=True,nsaflag=False):
        if plotsingle:
            plt.figure()
        w3res = 6.5 # arcsec
        w4res = 12 # arcsec
        w3pixelscale = 2.75
        # plot histogram of SERSIC_TH90
        mybins = np.linspace(0,50,50)
        if nsaflag:
            flag1 = self.allbutsizeflag
            flag2 = self.sampleflag            
            cat = self.nsa
            rad = cat['SERSIC_TH50']
        else:
            flag1 = self.allbutsizeflag2
            flag2 = self.sampleflag2            
            cat = self.a100nsa
            rad = self.r50

        t = plt.hist(rad[flag1],bins=mybins,histtype='step',lw=2,label='A100+NSA')
        flag = self.sampleflag
        t = plt.hist(rad[flag2],bins=mybins,histtype='step',lw=2,label='GALFIT size',hatch='//')
        #plt.axvline(x=w3res,c='r',label='WISE 12um resolution',ls='--')
        plt.axvline(x=sersicmin,c='k',label='$W3 \ FWHM (6.5\'\')$',ls='--')
        plt.legend(fontsize=16)
        plt.xlim(-1,50)
        plt.xlabel('r-band Half-Light Radius (arcsec)',fontsize=16)
        plt.ylabel('Number of galaxies',fontsize=16)
        plt.title('(b) GALFIT modeling',fontsize=20)
        if plotsingle:
            plt.savefig(plotdir+'/rband-size-hist.png')
            plt.savefig(plotdir+'/rband-size-hist.pdf')        
    def compare_size_resolution2(self,plotsingle=True,nsaflag=False):
        '''histogram of petro r90 with wise resolution'''
        if plotsingle:
            plt.figure()
        w3res = 6.5# arcsec
        # for comparison, MIPS 24um resolution is 6"
        w4res = 12 # arcsec
        w3pixelscale = 2.75
        # plot histogram of SERSIC_TH90
        xmax=250
        nbins=50
        mybins = np.linspace(0,xmax,nbins)
        flag = self.allbutsizeflag

        if nsaflag:
            flag1 = self.allbutsizeflag
            flag2 = self.sampleflag            
            cat = self.nsa
            rad = cat['PETROTH90']
        else:
            flag1 = self.allbutsizeflag2
            flag2 = self.sampleflag2 & (self.r90 > 6)            
            cat = self.a100nsa
            rad = self.D90
        np_size_cut = 10*w3res
        t = plt.hist(rad[flag1],bins=mybins,histtype='step',lw=2,label='A100+NSA')
        t = plt.hist(rad[flag2],bins=mybins,histtype='step',lw=2,label='GALFIT Size',hatch='//')
        flag = flag2 & (self.D90 > np_size_cut)
        print('number in well-resolved sample = ',sum(flag))
        print('number in GALFIT sample = ',sum(self.sampleflag2))        
        t = plt.hist(rad[flag],bins=mybins,histtype='step',lw=2,label='Non-Param Size',hatch='\\\\')                
        #plt.axvline(x=w3res,c='r',label='WISE 12um resolution',ls='--')
        plt.axvline(x=w3res,c='k',label=r'$W3\ FWHM$',ls='--')
        #plt.axvline(x=np_size_cut,c='k',label=r'Non-param size',ls=':')        
        plt.legend(fontsize=16)
        plt.xlim(-1,xmax)
        plt.xlabel('r-band diameter (3.16xPetro R90) (arcsec)',fontsize=16)
        plt.ylabel('Number of galaxies',fontsize=16)
        plt.title('(c) Non-Param Profiles',fontsize=20)
        if plotsingle:
            plt.savefig(plotdir+'/rpetro-size-hist.png')
            plt.savefig(plotdir+'/rpetro-size-hist.pdf')
    def plot_sample(self):
        plt.figure(figsize=(16,4))
        plt.subplots_adjust(left=.05,wspace=.3,bottom=.2,right=.95)        
        plt.subplot(1,3,1)
        self.plot_color_mstar2(plotsingle=False)
        plt.subplot(1,3,2)
        self.compare_size_resolution(plotsingle=False,nsaflag=False)
        plt.subplot(1,3,3)
        self.compare_size_resolution2(plotsingle=False,nsaflag=False)
        plt.savefig(plotdir+'/sample-overview.png')
        plt.savefig(plotdir+'/sample-overview.pdf')        
    def plot_env(self):
        '''plot field, filament, groups/clusters from Tempel+2011 catalogs'''
        fig = plt.figure(figsize=(8,10))

        ax = fig.add_subplot(111)#, projection='lambert')



        nsubplot = 1
        plt.subplots_adjust(hspace=.25)
        env_flags = [self.field_flag,self.filament_flag,(self.cluster_flag )]
        labels=['$Field$','$Filament$','$Clusters$']
        for i in range(len(env_flags)):
            flag = env_flags[i]
            plt.subplot(3,1,nsubplot)
            plt.scatter(self.nsa['RA'][self.sampleflag], self.nsa['DEC'][self.sampleflag],label='_nolegend_',s=3,c='k',alpha=.05)            
            plt.scatter(self.nsa['RA'][self.sampleflag][flag], self.nsa['DEC'][self.sampleflag][flag],label=labels[i],s=5,c=mycolors[i],alpha=.8)
            # OR
            props = dict(boxstyle='square', facecolor='w', alpha=0.5)

            #plt.axis([-10,370,-95,95])
            plt.grid(True)
            plt.title(labels[i],fontsize=16)        
            #plt.legend(loc='lower center',markerscale=2,ncol=2,fontsize=12)
            plt.ylabel('$DEC \ (deg)$',fontsize=20)
            nsubplot += 1
            plt.xlim(105,270)
            plt.ylim(-5,70)
            plt.gca().invert_xaxis()
        # add groups and poor groups
        
        plt.scatter(self.nsa['RA'][self.sampleflag][self.group_flag], self.nsa['DEC'][self.sampleflag][self.group_flag],label='$Groups$',s=5,c=mycolors[i+1],alpha=.8)
        plt.scatter(self.nsa['RA'][self.sampleflag][self.poorgroup_flag], self.nsa['DEC'][self.sampleflag][self.poorgroup_flag],label='$Poor\ Groups$',s=5,c=mycolors[i+2],alpha=.8)
        plt.legend()
        plt.xlabel('$RA \ (deg)$',fontsize=20)                    
        plt.savefig(plotdir+'env-skyplot.pdf')
        plt.savefig(plotdir+'env-skyplot.png')
    def plot_env2(self):
        '''plot field, filament+groups/clusters from Tempel+2011 catalogs'''
        fig = plt.figure(figsize=(8,10))
        ax = fig.add_subplot(111)#, projection='lambert')



        nsubplot = 1
        plt.subplots_adjust(hspace=.25)
        env_flags = [self.field_flag,self.filament_flag]
        labels=['$Field$','$Filament$','$Clusters$']
        titles=['$Field \ & \ LCS \ Clusters$','$Filaments, \ Groups, \ & \ Clusters$','$Clusters$']        
        colors = [mycolors[0],'k']
        for i in range(len(env_flags)):
            flag = env_flags[i]
            plt.subplot(2,1,nsubplot)
            #plt.scatter(self.nsa['RA'][self.sampleflag], self.nsa['DEC'][self.sampleflag],label='_nolegend_',s=3,c='k',alpha=.05)            
            plt.scatter(self.nsa['RA'][self.sampleflag][flag], self.nsa['DEC'][self.sampleflag][flag],label=labels[i],s=5,c=colors[i],alpha=.8)
            # OR
            props = dict(boxstyle='square', facecolor='w', alpha=0.5)

            #plt.axis([-10,370,-95,95])
            plt.grid(True)
            plt.title(titles[i],fontsize=16)        
            #plt.legend(loc='lower center',markerscale=2,ncol=2,fontsize=12)
            plt.ylabel('$DEC \ (deg)$',fontsize=20)
            nsubplot += 1
            if i == 0:
                for j in range(len(clusternames)):
                    s=' '+clusternames[j]

                    plt.text(clusterRA[clusternames[j]],clusterDec[clusternames[j]],s,fontsize=12,color='r',bbox=dict(facecolor='white', alpha=0.5))
                    drawbox(cluster24Box[clusternames[j]],'r-')                    

            plt.xlim(105,270)
            plt.ylim(-5,70)
            plt.gca().invert_xaxis()
        # add groups and poor groups
        flag = self.group_flag
        plt.scatter(self.nsa['RA'][self.sampleflag][flag], self.nsa['DEC'][self.sampleflag][flag],label='$Groups$',s=5,c=mycolors[1],alpha=.8)
        flag = self.cluster_flag
        plt.scatter(self.nsa['RA'][self.sampleflag][flag], self.nsa['DEC'][self.sampleflag][flag],label='$Clusters$',s=5,c=mycolors[i+1],alpha=.8)
        plt.legend(fontsize=14,markerscale=3)
        plt.xlabel('$RA \ (deg)$',fontsize=20)
        plt.savefig(plotdir+'env-skyplot.pdf')
        plt.savefig(plotdir+'env-skyplot.png')

    def plot_env3(self):
        '''plot field, filament+groups/clusters from Tempel+2011 catalogs'''
        fig = plt.figure(figsize=(10,8))
        ax = fig.add_subplot(111)#, projection='lambert')



        nsubplot = 1
        plt.subplots_adjust(hspace=.25)
        env_flags = [self.field_flag,self.filament_flag]
        labels=['$Field$','$Filament$','$Clusters$']
        titles=['$Field \ & \ LCS \ Clusters$','$Filaments, \ Groups, \ & \ Clusters$','$Clusters$']        
        colors = [mycolors[0],'k']
        for i in range(len(env_flags)):
            flag = env_flags[i]
            plt.subplot(1,2,nsubplot)
            #plt.scatter(self.nsa['RA'][self.sampleflag], self.nsa['DEC'][self.sampleflag],label='_nolegend_',s=3,c='k',alpha=.05)            
            plt.scatter(self.nsa['RA'][self.sampleflag][flag], self.nsa['DEC'][self.sampleflag][flag],label=labels[i],s=5,c=colors[i],alpha=.8)
            # OR
            props = dict(boxstyle='square', facecolor='w', alpha=0.5)

            #plt.axis([-10,370,-95,95])
            plt.grid(True)
            plt.title(titles[i],fontsize=16)        
            #plt.legend(loc='lower center',markerscale=2,ncol=2,fontsize=12)
            plt.ylabel('$DEC \ (deg)$',fontsize=20)
            nsubplot += 1
            if i == 0:
                for j in range(len(clusternames)):
                    s=' '+clusternames[j]

                    plt.text(clusterRA[clusternames[j]],clusterDec[clusternames[j]],s,fontsize=12,color='r',bbox=dict(facecolor='white', alpha=0.5))
                    drawbox(cluster24Box[clusternames[j]],'r-')                    

            plt.xlim(105,270)
            plt.ylim(-5,70)
            plt.gca().invert_xaxis()
        # add groups and poor groups
        flag = self.group_flag
        plt.scatter(self.nsa['RA'][self.sampleflag][flag], self.nsa['DEC'][self.sampleflag][flag],label='$Groups$',s=5,c=mycolors[1],alpha=.8)
        flag = self.cluster_flag
        plt.scatter(self.nsa['RA'][self.sampleflag][flag], self.nsa['DEC'][self.sampleflag][flag],label='$Clusters$',s=5,c=mycolors[i+1],alpha=.8)
        plt.legend(fontsize=14,markerscale=3)
        plt.xlabel('$RA \ (deg)$',fontsize=20)
        plt.savefig(plotdir+'env-skyplot.pdf')
        plt.savefig(plotdir+'env-skyplot.png')
        
if __name__ == '__main__':
    cats = catalogs('/home/rfinn/research/NSA/',virgoflag=False)
