#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
from pylab import figure, close, show, find
from h5py import File as h5
from numpy import (
exp, sqrt, power, nan, array, ones, zeros
)

M_PI = np.pi
M_E  = np.e


def gauss_normal(x):		# esto tiene integral =1.
	phi 	= exp(-.5 * x**2.) / sqrt(2.*M_PI)
	return phi

def gauss(x, A, mu, sig):	# esto tmb tiene integral =1.
	return A * (1./sig) * gauss_normal((x - mu) / sig)

def fit_ne_distribution(fname_inp, nbin=50):
    """
    - leer la coleccion de valores del nro de interacc esenciales
    - construir distribucion de dichos valores y fitear con gaussiana
    """
    from scipy.optimize import curve_fit
    # leemos la data de la simulacion anterior (la coleccion de valores
    # del nro de interacciones esenciales)
    with h5(fname_inp,'r') as f:
        ne = f['ne'][3:]  # obviamos sus primeros tres valores para 
                          # olvidarnos un poco de la configuracion original
                          # NOTE: la config original corresponde al `f['ne'][0]`.

    # hallamos su distribucion frecuentista
    # NOTE:
    #   - `hx_` son los bordes del dominio bineado
    #   - density=True es para q devuelva una distribucion con area=1.
    hc, hx_ = np.histogram(ne, bins=nbin, density=True) 
    hx = 0.5*(hx_[:-1] + hx_[1:])  # `hx` son los valores centrados de c/bin
    hx_mean = (hx*hc).sum()/hc.sum()

    # hacemos el ajuste de la distribucion frecuentista para hallar
    # su densidad de probabilidad
    popt, pcov = curve_fit(
        f = gauss, 
        xdata = hx,
        ydata = hc,
        p0 = [1., hx_mean, 30.], # semillas
    )
    return popt, pcov



#EOF
