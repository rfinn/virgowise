import os
import numpy as np
import argparse
from astropy.io import fits
import wget
import tarfile
import glob
import gzip

os.sys.path.append('/Users/rfinn/github/virgowise/')
from rungalfit import *
import g as g



#%run ~/github/virgowise/g.py --virgopath /home/rfinn/github/Virgo/ --nsapath /home/rfinn/github/Virgo/tables/ --getwise True


#def callg(X0):
#execfile('/home/obsastro10/github/virgowise/g.py')
#process_list(listname,band=3)
#    return Xnew
    

#def main():
#    Xmin = [.5,2.0,0,0,-89.0]
#    Xmax = [6.0,16.0,60.0,1,181.0]
#    MCsteps = 10
#    B = 1
#    Q = [[]]
    #for i in MCsteps:
    #E = 100000
    #    while(np.random.random()>= np.exp(-B*E)):
    #        nsersic = 5.5*np.random.random()+.5
    #        mag  = 14*np.random.random()+2
    #        re = 60*np.random.random()
    #        BA = np.random.random()
    #        PA =181*np.random.random()-89.0
    #        X0 = np.array(nsersic,mag,re,BA,PA)
    #        nQ = len(Q)
    #        D = np.array()
    #        for k in nQ-1:
    #            R = X0-Q[k,:]
    #            R5 = np.norm(R)**5
    #            deltaD = R/R5
    #            D = D+deltaD
    #        E = np.norm(D)**2
            
    #    Xnew = g(X0)
    #    Q = Q.append(Xnew)
 

#    X0=np.array([0.5,8.0,45.0,0.5,90])
#    Xnew = callg(X0)

X0=np.array([0.5,8.0,45.0,0.5,90])
#Xnew = callg(X0)
execfile('/home/obsastro10/github/virgowise/g.py')

#os.system('/home/obsastro10/github/virgowise/g.py')

Xnew=g.process_list(listname,band=3)


#if __name__ == "__main__":
#    main()
