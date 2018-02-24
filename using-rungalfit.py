
# coding: utf-8

# # Goal # 
# 
# The goal of this notebook is to show you how to use the galfit class in rungalfit.py
# 
# ### Required Modules ###
# 
# pyds9

# In[1]:

import os
import pyds9


# In[2]:

os.sys.path.append('/Users/rfinn/github/Virgo/programs/')


# In[3]:

from rungalfit import *


# In[21]:

# define image properties
galname = 'test'
image = 'testimage.fits'
sigma_image = 'testimage-sigma.fits'
psf_image = 'testpsf.fits'
psf_oversampling = 1
mask_image = 'testimage_mask.fits'
xminfit=1
yminfit=1
xmaxfit=100
ymaxfit=100
convolution_size=100
magzp=20.
pscale=2.5
convflag=1 # apply psf convolution
constraintflag=1 # add a constraint file?
fitallflag=0
ncomp=1


# In[22]:

# define first guess sersic parameters for galaxy 
xc=50
yc=50
nsersic = 2
mag = 7
Re = 10.
BA = .8
PA = 0


# In[8]:

gal1 = galfit(galname=galname,image=image,sigma_image=sigma_image,psf_image=psf_image,psf_oversampling=psf_oversampling,mask_image=mask_image,xminfit=xminfit,yminfit=yminfit,xmaxfit=xmaxfit,ymaxfit=ymaxfit,convolution_size=convolution_size,magzp=magzp,pscale=pscale,convflag=convflag,constraintflag=constraintflag,fitallflag=fitallflag,ncomp=ncomp)


# ## Create an input file for galfit ##

# In[23]:

# create output names
gal1.create_output_names()
gal1.open_galfit_input()
gal1.write_image_params()
gal1.add_simple_sersic_object(1,'sersic',xc,yc,mag,Re,nsersic,BA,PA)
gal1.set_sky(0)
gal1.write_sky(2)
gal1.close_input_file()


# In[ ]:



