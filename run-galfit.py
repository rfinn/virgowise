'''
Natasha Collova

GOAL: This code will be used to run galfit with varying initial cond's and retrieve its output to be saved in a tar file for each galaxy; going to start with just 1 galaxy tho
'''
import glob
import os
from astropy.io import fits
from astropy.wcs import WCS
import argparse
import subprocess

parser = argparse.ArgumentParser(description ='Run galfit and store output with best fit parameters into a tar file')
parser.add_argument('--d',dest = 'd', default =' ~/usr/local/bin/galfit', help = 'Locates path of galfit program')
parser.add_argument('--galfit', dest = 'galfit', default = True, action = 'store_true', help = 'Run galfit to model galaxy/galaxies with varying initial conditions')
parser.add_argument('--store', dest = 'store', default = False, action = 'store_true', help = 'Store output from galfit as a table into a single tar file')


args = parser.parse_args()

#I need to learn how to write into an input file from python
#Locate input files of galaxies/create input files
galaxies = 

#Run galfit
if args.galfit:
    for f in galaxies:
        os.system
