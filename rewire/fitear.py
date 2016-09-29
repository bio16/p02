#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import funcs as ff
from h5py import File as h5
import argparse

#--- retrieve args
parser = argparse.ArgumentParser(
formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument(
'-inp', '--fname_inp',
type=str,
default='./LIT.h5',
help='input filename of the network/graph.',
)
pa = parser.parse_args()


fname_inp = pa.fname_inp #'./test.h5'
popt, pcov = ff.fit_ne_distribution(fname_inp, nbin=50)

f = h5(fname_inp,'r+')
f['fit/A'] = popt[0]
f['fit/mu'] = popt[1]
f['fit/sigma'] = popt[2]

#EOF
