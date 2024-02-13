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
plotdir = homedir+'/proposals/NSF2022/figures/'
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

class plots():

    def plot_env22(self):
        '''plot field, filament+groups/clusters from Tempel+2011 catalogs'''
        fig = plt.figure(figsize=(14,5))
        plt.subplots_adjust(bottom=.15,left=.1,right=.95,top=.95)
        ax = fig.add_subplot(111)#, projection='lambert')



        nsubplot = 1
        plt.subplots_adjust(hspace=.25)
        env_flags = [self.field_flag,self.filament_flag]
        colors = [mycolors[0],'k']
        # plot full sample
        plt.scatter(self.nsa['RA'][self.sampleflag],self.nsa['DEC'][self.sampleflag],s=4,c=self.nsa['Z'][self.sampleflag],alpha=.7,label='This Study')
        # plot Virgo
        #plt.scatter(self.virgo['RA'],self.virgo['DEC'],s=1,c='c',alpha=.2,label='Virgo Pilot')
        for j in range(len(clusternames)):
            s=' '+clusternames[j]

            plt.text(clusterRA[clusternames[j]],clusterDec[clusternames[j]],s,fontsize=16,color='r',bbox=dict(facecolor='white',edgecolor='white', alpha=0.5))
            drawbox(cluster24Box[clusternames[j]],'r-')                    

        # plot LCS
            #plt.legend(loc='lower center',markerscale=2,ncol=2,fontsize=12)

        #plt.legend()
        plt.colorbar(label='redshift',fraction=.08)
        plt.xlim(105,270)
        plt.ylim(-5,70)
        plt.gca().invert_xaxis()
        # add groups and poor groups
        flag = self.group_flag
        #plt.scatter(self.nsa['RA'][self.sampleflag][flag], self.nsa['DEC'][self.sampleflag][flag],label='$Groups$',s=5,c=mycolors[1],alpha=.8)
        flag = self.cluster_flag
        #plt.scatter(self.nsa['RA'][self.sampleflag][flag], self.nsa['DEC'][self.sampleflag][flag],label='$Clusters$',s=5,c='k',alpha=.8)
        #plt.legend(fontsize=18,markerscale=3)
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)        
        plt.xlabel('$RA \ (deg)$',fontsize=26)
        plt.ylabel('$DEC \ (deg)$',fontsize=26)        
        plt.savefig(plotdir+'env-skyplot22.pdf')
        plt.savefig(plotdir+'env-skyplot22.png',dpi=300)


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
        plt.xlabel('$RA \ (deg)$',fontsize=2)
        plt.savefig(plotdir+'env-skyplot.pdf')
        plt.savefig(plotdir+'env-skyplot.png')


class catalogs(plots):
    def __init__(self,catalog_path,virgoflag = False):
        # read in nsa, wise, co catalogs
        if virgoflag:
            self.nsatab = catalog_path + 'nsa.virgo.fits'
            self.wisetab = catalog_path + 'nsa_wise.virgo.fits'
            self.cotab = catalog_path + 'nsa_CO-Gianluca.virgo.fits'
            self.co = fits.getdata(self.cotab)
        else:
            self.nsatab = catalog_path + 'nsa_v0_1_2.fits'
            self.wisetab = catalog_path + 'nsa_v0_1_2_wise.fits'
            self.jmasstab = catalog_path + 'nsa_v1_2_fsps_v2.4_miles_chab_charlot_sfhgrid01.fits'
        self.nsa = fits.getdata(self.nsatab)
        self.wise = fits.getdata(self.wisetab)
        self.jmass = fits.getdata(self.jmasstab)

        self.nsadict=dict((a,b) for a,b in zip(self.nsa.NSAID,np.arange(len(self.nsa.NSAID)))) #useful for us can easily look up galaxy ID's       
    def calcagn(self):
        self.AGNKAUFF= ((np.log10(self.nsa.O3FLUX/self.nsa.HBFLUX) > (.61/(np.log10(self.nsa.N2FLUX/self.nsa.HAFLUX)-.05)+1.3)) | (np.log10(self.nsa.N2FLUX/self.nsa.HAFLUX) > 0.))
#y=(.61/(x-.47)+1.19)
        self.AGNKEWLEY= ((np.log10(self.nsa.O3FLUX/self.nsa.HBFLUX) > (.61/(np.log10(self.nsa.N2FLUX/self.nsa.HAFLUX)-.47)+1.19)) | (np.log10(self.nsa.N2FLUX/self.nsa.HAFLUX) > 0.3))
        self.agnflag = self.AGNKAUFF
        self.agnflag = self.AGNKEWLEY
        
    def define_sample(self):
        self.w3_flag = self.wise.W3SNR>10 
        self.w4_flag = self.wise.W4SNR>5      
        #self.co_flag = self.co.COdetected == '1'
        #self.sampleflag = self.w3_flag & self.w4_flag & self.co_flag
        
        self.sampleflag = (self.wise.W3MPRO < 10.) & (self.nsa.PETROTH50 > 2.) \
        & (self.nsa.ZDIST > 0.005)& (self.nsa.ZDIST < .03) & ~self.agnflag
        print('number of galaxies in sample = ',sum(self.sampleflag))
        self.allbutsizeflag = (self.wise.W3MPRO < 10.)  \
        & (self.nsa.ZDIST > 0.005)& (self.nsa.ZDIST < .03) & ~self.agnflag
        print('number of galaxies in sample = ',sum(self.sampleflag))
    def write_sample_table(self):
        outfile = '/home/rfinn/research/wisesize/nsa-wisesize-catalog.fits'
        fits.writeto(outfile,self.nsa[self.sampleflag],overwrite=True)


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

    def get_virgo(self):
        cat = homedir+'/research/Virgo/tables-north/v2/vf_v2_main.fits'
        self.virgo = Table.read(cat)

if __name__ == '__main__':
    cats = catalogs('/home/rfinn/research/NSA/',virgoflag=False)
    cats.calcagn()
    cats.define_sample()
    cats.get_tempel()
    cats.get_virgo()
