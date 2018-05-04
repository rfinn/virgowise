# virgowise
code to analyze wise 12um images

(0) get-unwise-images.py : Download unwise images (you can also do
this in using-rungalfit.py)

from within ipython
%run ~/github/virgowise/get-unwise-images.py --nsapath '/Users/rfinn/github/Virgo/tables/' --band 3

(1) maskimage.py : create image masks using sextractor segmentation image

(2) using-rungalfit.py : run galfit on images

(3) read_galfit_results.py : reads galfit output fits files, makes a few plots
Still to do: to write code to extract fit parameters and store in a table




