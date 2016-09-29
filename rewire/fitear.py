#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import funcs as ff
from h5py import File as h5

fname_inp = './test.h5'
popt, pcov = ff.fit_ne_distribution(fname_inp, nbin=50)

f = h5(fname_inp,'r+')
f['fit/A'] = popt[0]
f['fit/mu'] = popt[1]
f['fit/sigma'] = popt[2]

#EOF
