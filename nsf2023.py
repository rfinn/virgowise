#!/usr/bin/env python
import os
from matplotlib import pyplot as plt
import numpy as np

from astropy.table import Table
from astropy.io import fits


import astropy.units as u
from astropy.coordinates import SkyCoord
from matplotlib.patches import Circle


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
sersicmin = 5.5 # arcsec
snrcut = 20


# for UAT Telescope Survey
zmin = 0.002
zmax = 0.025
sersicmin = 3

mindiameter = 10
snrcut = 20

filter_wavec = 656.3+np.array([0,4,8,12,16]) # nm
filter_width = 8.0 # nm
filter_wave_min = filter_wavec - filter_width/2
filter_wave_max = filter_wavec + filter_width/2

# redshift - delta lambda/lamba = lobs/lem - 1
#
filter_zcenter = filter_wavec/656.3 - 1
filter_zmin = filter_wave_min/656.3 - 1
filter_zmax = filter_wave_max/656.3 - 1

def get_filter(z):
    """
    takes in the redshift

    returns the zmin, zmax of the filter that the galaxy falls in
    """

    # find filter with center closest to the redshift

    dz = np.abs(z - filter_zcenter)
    bestFilter = dz == np.min(dz)

    # check that gal is within filter limits (just a double check...)
    if (z > filter_zmin[bestFilter]) & (z < filter_zmax[bestFilter]):
            return filter_zmin[bestFilter][0],filter_zmax[bestFilter][0]
    return -99, -99

    
def find_groups(ra,dec,maxoffset=0.4):
    ngroups=0
    # keep track on which are group centers
    groupCenterFlag = np.zeros(len(ra),'bool')
    
    maxgroupsize = 10
    cat1 = SkyCoord(ra*u.deg,dec*u.deg, frame='icrs')
    cat2 = SkyCoord(ra*u.deg,dec*u.deg, frame='icrs')
    
    # match and remove singles
    idxcat2all, idxcat1all, d2d, d3d = cat1.search_around_sky(cat2, maxoffset*u.deg)
    # remove self matches
    idxcat2 = idxcat2all[d2d > 0.]
    idxcat1 = idxcat1all[d2d > 0.]
    unique_values, counts = np.unique(idxcat2, return_counts=True)
    
    
    single_flag = np.ones(len(cat1),'bool')
    
    for id in unique_values[counts > 1]: 
        single_flag[id] = False
    singlemembs = np.arange(len(cat1))[single_flag]

    groupCenterFlag[single_flag] = np.ones(np.sum(single_flag),'bool')
    
    # track ra and dec of all pointings
    pointing_ra = []
    pointing_dec = []
 
    # get centers of single pointings
    #for index in singlemembs:
    #    pointing_ra.append(ra[index])
    #    pointing_dec.append(dec[index])
    singles_ra = []
    singles_dec = []
    for id in singlemembs:
        singles_ra.append(ra[id])
        singles_dec.append(dec[id])        
        
    #cat1 = cat1[~single_flag]
    #cat2 = cat2[~single_flag]

    #
    print()
    print("the number of single pointings = ",np.sum(single_flag))
    nsingles = np.sum(single_flag)
    print("number of remaining galaxies in groups = ",np.sum(counts > 1))
    print("testing: ",len(cat1),np.sum(single_flag)+ np.sum(counts > 1))
    
    #print(len(idxcat1),len(idxcat1all))

    nmembgroup = []
    groupmembs = []
    surviveFlag = ~groupCenterFlag
    print("number of galaxies surviving = ",np.sum(surviveFlag))
    remove_flag = np.zeros(len(cat1),'bool')
    nloop=0
    while maxgroupsize > 1:
        
        # find objects with match w/in max offset
        #print("matching catalogs...")
        #print()
        #print("number of galaxies surviving = ",np.sum(surviveFlag))
        #print(f"length of cat1[surviveFlag] = {len(cat1[surviveFlag])}")
        idxcat2all, idxcat1all, d2d, d3d = cat1[surviveFlag].search_around_sky(cat2[surviveFlag], maxoffset*u.deg)
        # remove self matches
        idxcat2 = idxcat2all[d2d > 0.]
        idxcat1 = idxcat1all[d2d > 0.]

        # translate this back to the indices of the uncut catalogs
        idxcat2 = np.arange(len(cat1))[surviveFlag][idxcat2]
        idxcat1 = np.arange(len(cat1))[surviveFlag][idxcat1]

        
        # find pointing with max number
        #print("finding unique ids ...")
        unique_values, counts = np.unique(idxcat2, return_counts=True)
        #print("number of unique groups = ",len(unique_values))
        maxgroupsize = np.max(counts)
        maxgroup_id = unique_values[counts == np.max(counts)]
        


        
        #print(f"maxgroup_id = {maxgroup_id}, max counts = {np.max(counts)}")
        if np.max(counts) == 1:
            #  * if max number is 1, then the remaining sources are single pointings
            ngroups += len(unique_values)

            
            return nsingles,ngroups,pointing_ra+singles_ra,pointing_dec+singles_dec,nmembgroup,groupmembs,singlemembs,groupCenterFlag,single_flag
    
        # save as group 1
        #print("appending ra and dec")
        #print(f"maxgroup_id = {maxgroup_id}")
        #print(maxgroup_id[0])
        nmembgroup.append(maxgroupsize)
        ngroups += 1
        # remove group 1 members from idxcat1 and idxcat 2
        #remove_flag = np.zeros(len(idxcat1),'bool')
        
        if len(maxgroup_id) > 1:
            membs = idxcat1[idxcat2 == maxgroup_id[0]]
            
            groupCenterFlag[int(maxgroup_id[0])] = True
            pointing_ra.append(cat2.ra[int(maxgroup_id[0])].value)
            pointing_dec.append(cat2.dec[int(maxgroup_id[0])].value)
            
        else:
            membs = idxcat1[idxcat2 == maxgroup_id]
            groupCenterFlag[int(maxgroup_id[0])] = True
            pointing_ra.append(cat2.ra[int(maxgroup_id[0])].value)
            pointing_dec.append(cat2.dec[int(maxgroup_id[0])].value)
            
        #print("number of membs = ",len(membs))
        #print("building remove_flag")
        #print("members = ",membs)
        #print(f"before eliminating, survive flag of membs = {surviveFlag[membs]}")        
        for m in membs:
            # remove members from idxcat1 and 2
            #groupmembs.append(cat1['NSAID'][m])
            #ragroup.append(cat1['RA'][m])
            #matchmemb_cat1 = cat1 == m
            remove_flag[m] = True
            surviveFlag[m] = False
            #matchmemb_cat2 = cat2 == m
            remove_flag[m] = True
        # keep non group members
        #print('number of galaxies to remove = ',np.sum(remove_flag))
        #print(len(idxcat2))
        #idxcat1 = idxcat1[~remove_flag]
        #idxcat2 = idxcat2[~remove_flag]
        #cat1 = cat1[~remove_flag]
        #cat2 = cat2[~remove_flag]
        #return idxcat1,idxcat2
        #print(len(idxcat2))
        #print(f"number of groups = {ngroups}; max group size = {maxgroupsize}",len(cat1))
        #print()
        #print(f"survive flag of membs = {surviveFlag[membs]}")
        #print()
        #if ngroups > 8:
        #    break
        
    return nsingles,ngroups,pointing_ra+singles_ra,pointing_dec+singles_dec,nmembgroup,groupmembs,singlemembs,groupCenterFlag,single_flag

def find_groups_3d(ra,dec,z,maxoffset=0.4,verbose=False,testing=False):
    ngroups=0
    # keep track on which are group centers
    groupCenterFlag = np.zeros(len(ra),'bool')
    
    maxgroupsize = 10
    c = SkyCoord(ra=1*u.deg, dec=2*u.deg, radial_velocity=10*u.km/u.s)

    # search around 3d will give separation in 3d coords, but that is not what we want
    # so sticking with 2d search, and then will figure out how to include redshift
    cat1 = SkyCoord(ra*u.deg,dec*u.deg, frame='icrs')
    cat2 = SkyCoord(ra*u.deg,dec*u.deg, frame='icrs')
    
    # match and remove singles
    idxcat2all, idxcat1all, d2d, d3d = cat1.search_around_sky(cat2, maxoffset*u.deg)
    # remove self matches
    idxcat2 = idxcat2all[d2d > 0.]
    idxcat1 = idxcat1all[d2d > 0.]
    unique_values, counts = np.unique(idxcat2, return_counts=True)
    
    
    single_flag = np.ones(len(cat1),'bool')
    
    for id in unique_values[counts > 1]: 
        single_flag[id] = False
    singlemembs = np.arange(len(cat1))[single_flag]

    groupCenterFlag[single_flag] = np.ones(np.sum(single_flag),'bool')
    
    # track ra and dec of all pointings
    pointing_ra = []
    pointing_dec = []
 
    # get centers of single pointings
    #for index in singlemembs:
    #    pointing_ra.append(ra[index])
    #    pointing_dec.append(dec[index])
    singles_ra = []
    singles_dec = []
    for id in singlemembs:
        singles_ra.append(ra[id])
        singles_dec.append(dec[id])        
        
    #cat1 = cat1[~single_flag]
    #cat2 = cat2[~single_flag]

    #
    nsingles = np.sum(single_flag)
    if verbose:
        print("the number of single pointings = ",np.sum(single_flag))
    
        print("number of remaining galaxies in groups = ",np.sum(counts > 1))
        print("testing: ",len(cat1),np.sum(single_flag)+ np.sum(counts > 1))
    
    #print(len(idxcat1),len(idxcat1all))

    nmembgroup = []
    groupmembs = []
    surviveFlag = ~groupCenterFlag
    if verbose:
        print("number of galaxies surviving = ",np.sum(surviveFlag))
    remove_flag = np.zeros(len(cat1),'bool')
    nloop=0
    while maxgroupsize > 1:
        
        # find objects with match w/in max offset
        #print("matching catalogs...")
        if verbose:
            print()
            print("number of galaxies surviving = ",np.sum(surviveFlag))
            print(f"length of cat1[surviveFlag] = {len(cat1[surviveFlag])}")
        idxcat2all, idxcat1all, d2d, d3d = cat1[surviveFlag].search_around_sky(cat2[surviveFlag], maxoffset*u.deg)
        # remove self matches
        #idxcat2 = idxcat2all[d2d > 0.]
        #idxcat1 = idxcat1all[d2d > 0.]

        # translate this back to the indices of the uncut catalogs
        idxcat2 = np.arange(len(cat1))[surviveFlag][idxcat2all]
        idxcat1 = np.arange(len(cat1))[surviveFlag][idxcat1all]
        
        # find pointing with max number
        #print("finding unique ids ...")
        unique_values, counts = np.unique(idxcat2, return_counts=True)
        #print("number of unique groups = ",len(unique_values))
        maxgroupsize = np.max(counts)
        maxgroup_id = unique_values[counts == np.max(counts)]

        if verbose:
            print(f"maxgroup_id = {maxgroup_id}, max counts = {np.max(counts)}")
            
        if np.max(counts) == 1:
            #  * if max number is 1, then the remaining sources are single pointings
            ngroups += len(unique_values)

            
            return nsingles,ngroups,pointing_ra+singles_ra,pointing_dec+singles_dec,nmembgroup,groupmembs,singlemembs,groupCenterFlag
    
        # save as group 1
        nmembgroup.append(maxgroupsize)
        ngroups += 1
        # remove group 1 members from idxcat1 and idxcat 2
        
        if len(maxgroup_id) > 1:
            membs = idxcat1[idxcat2 == maxgroup_id[0]]
            centerID = maxgroup_id[0]
            groupCenterFlag[int(maxgroup_id[0])] = True
            pointing_ra.append(cat2.ra[int(maxgroup_id[0])].value)
            pointing_dec.append(cat2.dec[int(maxgroup_id[0])].value)
            
        else:
            membs = idxcat1[idxcat2 == maxgroup_id[0]]
            groupCenterFlag[int(maxgroup_id[0])] = True
            pointing_ra.append(cat2.ra[int(maxgroup_id[0])].value)
            pointing_dec.append(cat2.dec[int(maxgroup_id[0])].value)
            centerID = maxgroup_id[0]

        # get redshift range of the filter for this galaxy

        zmin,zmax = get_filter(z[centerID])
        zflag = (z[membs] > zmin) & (z[membs] < zmax)
        if verbose:
            print(f"centerID={centerID}")
            print(f"membs={membs}")
            print(f"zmin,zmax={zmin:.4f},{zmax:.4f},zobj={z[centerID]:.4f},{np.sum(zflag)},{z[membs]}")
        
        #print("number of membs = ",len(membs))
        #print("building remove_flag")
        #print("members = ",membs)
        #print(f"before eliminating, survive flag of membs = {surviveFlag[membs]}")        
        for m in membs:
            # remove members from idxcat1 and 2
            #groupmembs.append(cat1['NSAID'][m])
            #ragroup.append(cat1['RA'][m])
            #matchmemb_cat1 = cat1 == m
            if (z[m] > zmin) & (z[m] < zmax):
                remove_flag[m] = True
                surviveFlag[m] = False
                #matchmemb_cat2 = cat2 == m
                remove_flag[m] = True
        if testing:
            if ngroups > 10:
                break
        
    return nsingles,ngroups,pointing_ra+singles_ra,pointing_dec+singles_dec,nmembgroup,groupmembs,singlemembs,groupCenterFlag

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

    def plot_color_mag(self):
        plt.figure()
        plt.scatter(self.Mabsr[self.sampleflag],self.gmr[self.sampleflag],'k.',alpha=.5)
        
        pass

    def plot_color_mstar2(self,plotsingle=True,sampleflag=None):
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
        titles = ['Primary Targets','A100+NSA (z<0.025)']
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

            # Plot red sequence
            xl = np.linspace(8.25,11,10)
            dmag = 0.15
            yl = .1*(xl-10)+1.15
            plt.plot(xl,yl,'r--',label='Red sequence')

            # Plot cut to use to separate blue galaxies
            plt.plot(xl,yl-.15,'b-',lw=2,label='color cut')
            
            #plt.plot(xl,yl+dmag,'r--')
            #plt.plot(xl,yl-dmag,'r--')
            #plt.fill_between(xl,y1=yl+dmag,y2=yl-dmag,color='r',alpha=.08)
        #cb = plt.colorbar(fraction=.08,ax=allax)
        #cb.set_label('$N_{gal}$',fontsize=fsize)
        plt.legend()
        if plotsingle:
            plt.savefig(homedir+'proposals/NASA2018/color-mass-adap202.pdf')
        #plt.ylabel(r'$ NUV - r $')
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
    def plot_color_mstar(self,plotsingle=True,sampleflag=None):
        fsize = 18
        xmin,xmax=7,11.5
        ymin,ymax=0,1.7
        if plotsingle:
            plt.figure(figsize=(12,4))
            #plt.figure(figsize=(10,4))
            plt.subplots_adjust(left=.1,bottom=.15,wspace=.2,right=.95)            
        myflags = [self.sampleflag2,self.zflag2]# & self.agnflag2]
        titles = ['(b) Primary Targets','(a) Parent Sample']        
        if sampleflag is not None:
            myflags.append(sampleflag)
            titles.append('(c) Full Sample')
        NUVr=self.color
        stellarmass = self.logMstar
        allax=[]

        #titles = ['(a) Resolved Sample','A100+NSA']
        nplot = 0
        indices = [1,0]
        if sampleflag is not None:
            indices.append(2)
        for i in indices:
        #for i in range(1):
            nplot += 1
            print(i,nplot)
            plt.subplot(1,3,nplot)
            if nplot < 3:
                plt.hexbin(stellarmass[ myflags[i]],NUVr[myflags[i]],cmap='gray_r',extent=(xmin,xmax,ymin,ymax),gridsize=50)
            #plt.plot(stellarmass[ myflags[0]],NUVr[myflags[0]],'k.',c=mycolors[1],alpha=.1,label='Resolved')            
            plt.axis([xmin,xmax,ymin,ymax])
            plt.xlabel(r'$ log_{10}\left( M_\star /M_\odot \right) $',fontsize=fsize)
            s = 'N = %i'%(sum(myflags[i]))
            plt.text(0.05,.9,s,transform=plt.gca().transAxes,fontsize=fsize-2)
            #X,Y =  np.meshgrid(stellarmass[allflag],NUVr[allflag])
            #plt.colorbar()
            allax.append(plt.gca())
            plt.title(titles[i],fontsize=fsize)

            # Plot red sequence
            xl = np.linspace(7.25,11,10)
            dmag = 0.15
            yl = .1*(xl-10)+1.15
            plt.plot(xl,yl,'r--',label='Red sequence')

            # Plot cut to use to separate blue galaxies
            plt.plot(xl,yl-.15,'b-',lw=2,label='color cut')

            if i == 1:
                plt.ylabel(r'$ g - i $',fontsize=fsize)
                plt.legend(loc='lower right')
            if i == 2:
                plt.plot(stellarmass[ myflags[i]],NUVr[myflags[i]],'b.',alpha=.1)
        if plotsingle:
            plt.savefig(homedir+'proposals/NSF2023-UAT-MRI/figures/color-mass-mri.pdf')
            plt.savefig(homedir+'proposals/NSF2023-UAT-MRI/figures/color-mass-mri.png')


 
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


        self.a100nsa = Table.read(homedir+'/research/wisesize/a100-nsa-for-adap.fits')
        #zflag = (self.a100nsa['Vhelio'] < zmax*3e5) | (self.a100nsa['ZDIST'] < zmax)
        #self.a100nsa = self.a100nsa[zflag]
        self.nsadict=dict((a,b) for a,b in zip(self.nsa.NSAID,np.arange(len(self.nsa.NSAID)))) #useful for us can easily look up galaxy ID's       
        self.get_tempel()
        self.gmr = self.a100nsa['ABSMAG'][:,3] - self.a100nsa['ABSMAG'][:,4]        
        self.gmi = self.a100nsa['ABSMAG'][:,3] - self.a100nsa['ABSMAG'][:,5]
        absMag_i = self.a100nsa['ABSMAG'][:,5] - 5*np.log10(.7) # correct from HO=100 to 70
        self.absMag_r = self.a100nsa['ABSMAG'][:,4] - 5*np.log10(.7) # correct from HO=100 to 70        
        logMstarTaylor=1.15+0.70*(self.gmi) -0.4*(absMag_i)
        #self.a100Flag = np.zeros(len(self.a100nsa),'bool')
        #self.nsaFlag = np.zeros(len(self.a100nsa),'bool')        
        #self.a100Flag = ~self.a100nsa['AGC'].mask
        #self.nsaFlag = ~self.a100nsa['NSAID'].mask
        self.a100Flag = self.a100nsa['AGC'] > 0
        self.nsaFlag = self.a100nsa['NSAID'] > 0        
        self.logMstar = self.a100nsa['logMstarTaylor']*self.a100Flag + logMstarTaylor*(~self.a100Flag & self.nsaFlag)
        self.color = self.a100nsa['gmi_corr']*self.a100Flag + self.gmi*(~self.a100Flag)
        #self.calcagn()
        self.calcagn2()        
        #self.define_sample()
        self.define_sample2()        
        self.get_super_coords()

    def get_super_coords(self):
        # create a super RA
        agcflag = self.a100nsa['W50'] > 0

        nsaflag = self.a100nsa['NSAID'] > 0

        self.superRA = np.zeros(len(self.a100nsa),'d')
        self.superDEC = np.zeros(len(self.a100nsa),'d')
        self.superZ = np.zeros(len(self.a100nsa),'d')                
        self.superRA[agcflag] = self.a100nsa['RAdeg_Use'][agcflag]
        self.superRA[nsaflag & ~agcflag] =  self.a100nsa['RA'][nsaflag & ~agcflag]
        
        self.superDEC[agcflag] =self.a100nsa['DECdeg_Use'][agcflag]
        self.superDEC[nsaflag & ~agcflag] =  self.a100nsa['DEC'][nsaflag & ~agcflag]
        

        self.superZ[agcflag] = + self.a100nsa['Vhelio'][agcflag]/3.e5
        self.superZ[nsaflag & ~agcflag] = self.a100nsa['Z'][nsaflag & ~agcflag]
        
        self.a100nsa.add_column(self.superRA,name='superRA')
        self.a100nsa.add_column(self.superDEC,name='superDEC')
        self.a100nsa.add_column(self.superZ,name='superZ')        
        
    def calcagn(self):
        self.AGNKAUFF= ((np.log10(self.nsa.O3FLUX/self.nsa.HBFLUX) > (.61/(np.log10(self.nsa.N2FLUX/self.nsa.HAFLUX)-.05)+1.3)) | (np.log10(self.nsa.N2FLUX/self.nsa.HAFLUX) > 0.))
#y=(.61/(x-.47)+1.19)
        self.AGNKEWLEY= ((np.log10(self.nsa.O3FLUX/self.nsa.HBFLUX) > (.61/(np.log10(self.nsa.N2FLUX/self.nsa.HAFLUX)-.47)+1.19)) | (np.log10(self.nsa.N2FLUX/self.nsa.HAFLUX) > 0.3))
        self.agnflag = self.AGNKAUFF
        self.agnflag = self.AGNKEWLEY

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
        

    def define_sample2(self):

        self.w3_flag = (self.a100nsa['w3_nanomaggies']*np.sqrt(self.a100nsa['w3_nanomaggies_ivar']) > snrcut) | (self.a100nsa['W3SNR'] > snrcut)
        #self.w3_flag = self.w3_snr > 10
        self.nuv_snr = self.a100nsa['NMGY'][:,1]*np.sqrt(self.a100nsa['NMGY_IVAR'][:,1])
        self.nuv_flag = (self.nuv_snr > snrcut)
        self.snr_flag = self.w3_flag | self.nuv_flag


        self.zflag2 = \
             (((self.a100nsa['ZDIST'] > zmin)& (self.a100nsa['ZDIST'] < zmax)) | \
              ((self.a100nsa['Vhelio'] > zmin*3e5) & (self.a100nsa['Vhelio'] < zmax*3e5))) 
        
        # Adding a color flag for the Halpha sample        
        self.color_flag = self.color < (0.1*(self.logMstar-10) + 1)

        self.mass_flag = self.logMstar > 8.
        
        self.w4_snr = (self.a100nsa['w4_nanomaggies']*np.sqrt(self.a100nsa['w4_nanomaggies_ivar']))
        self.w4_flag = self.w4_snr > snrcut

        self.sizeflag = ((self.a100nsa['SERSIC_TH50'] > sersicmin) | (self.a100nsa['petroR50_r'] > sersicmin))

        self.sizeflag2 = ((self.a100nsa['SERSIC_TH50'] > sersicmin) | (self.a100nsa['petroR50_r'] > sersicmin))        
             #& self.snr_flag 
        self.allbutsizeflag2 = self.zflag2 & self.snr_flag 
        self.r50 = self.a100nsa['SERSIC_TH50']*self.nsaFlag + self.a100nsa['petroR50_r']*(~self.nsaFlag & self.a100Flag)
        self.r90 = 1.58*(self.a100nsa['PETROTH90']*self.nsaFlag + self.a100nsa['petroR90_r']*(~self.nsaFlag & self.a100Flag))
        self.D90 = 2*self.r90

        self.diameter_flag = self.D90 > mindiameter

        # add a g-r cut

        self.sampleflag2 = self.zflag2 & self.mass_flag & self.color_flag & ~self.agnflag2 & self.sizeflag2 & self.diameter_flag
        print('parent sample (zcut ) = ',sum(self.zflag2))
        print('parent sample (zcut & mass ) = ',sum(self.zflag2 & self.mass_flag))        

        print("sample w/zcut,mass, color cut = ",sum(self.zflag2 & self.color_flag & self.mass_flag))
        print("sample w/zcut,mass, color cut, no AGN = ",sum(self.zflag2 & self.color_flag & self.mass_flag & (~self.agnflag2)))
        print("sample w/zcut,mass, color cut, no AGN,size = ",sum(self.zflag2 & self.color_flag & self.mass_flag & (~self.agnflag2) & self.sizeflag2))
        print("sample w/zcut,mass, color cut, no AGN,diam = ",sum(self.zflag2 & self.color_flag & self.mass_flag & (~self.agnflag2) & self.sizeflag2 & self.diameter_flag))          
        
        print('resolved sample (not AGN) = ',sum(self.sampleflag2 & (~self.agnflag2)))
        print('resolved sample (not AGN) = %.3f'%(sum(self.sampleflag2 & (~self.agnflag2))/sum(self.allbutsizeflag2)))

        #print('resolved sample (AGN) = ',sum(self.sampleflag2 & self.agnflag2))
        #print('resolved sample (AGN) = %.3f'%(sum(self.sampleflag2 & (self.agnflag2))/sum(self.allbutsizeflag2)))
        #print('resolved sample wise AGN = ',sum(self.sampleflag2 & self.wiseagn))
        #print('resolved sample wise AGN = %.3f'%(sum(self.sampleflag2 & (self.wiseagn))/sum(self.allbutsizeflag2)))
        #self.sampleflag2 = self.sampleflag2 & ~self.agnflag2
        
    def write_sample_table(self):
        outfile = '/home/rfinn/research/wisesize/nsa-wisesize-catalog-uat.fits'
        self.a100nsa[self.sampleflag2].write(outfile,overwrite=True)

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
    #cats.define_sample2()
    cats.get_tempel()
    cats.get_virgo()
