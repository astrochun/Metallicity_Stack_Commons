import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii as asc
from astropy.table import vstack, hstack
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table, Column
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from os.path import exists
import numpy.ma as ma
from matplotlib.gridspec import GridSpec
from pylab import subplots_adjust
from astropy.convolution import Box1DKernel, convolve
from scipy.optimize import curve_fit
import scipy.integrate as integ
from mpl_toolkits.mplot3d import Axes3D
import sys

#For generalizing for several users
from getpass import getuser
from astropy import units as u
from chun_codes.cardelli import *

from . import k_dict

k_4363 = k_dict['OIII_4363']
k_5007 = k_dict['OIII_5007']

#Constants

a = 13205
b = 0.92506
c = 0.98062


def R_calculation(OIII4363, OIII5007, EBV):
  
    R_value = OIII4363/(OIII5007*(1+1/3.1))* 10**(0.4*EBV*(k_4363-k_5007))
    return R_value

def temp_calculation(R):
    #T_e = a(-log(R)-b)^(-c)
    T_e =  a*(-np.log10(R)-b)**(-1*c)     
    print(T_e)  
    return T_e

def metallicity_calculation(T_e, TWO_BETA, THREE_BETA):   #metallicity spelled wrong go back and change it if you have time
#(T_e,der_3727_HBETA, der_4959_HBETA, der_5007_HBETA, OIII5007, OIII4959, OIII4363, HBETA, OII3727, dustatt = False):
    #12 +log(O+/H) = log(OII/Hb) +5.961 +1.676/t_2 - 0.4logt_2 - 0.034t_2 + log(1+1.35x)
    #12 +log(O++/H) = log(OIII/Hb)+6.200+1.251/t_3 - 0.55log(t_3) - 0.014(t_3)
    #t_2 = 0.7*t_3 +0.17
    
    t_3 = T_e*1e-4
    t_2 = 0.7*t_3 +0.17
    x2 = 1e-4 * 1e3 * t_2**(-0.5)

    O_s_ion_log = np.log10(TWO_BETA) +5.961 +1.676/t_2 - 0.4*np.log10(t_2) - 0.034*t_2 + np.log10(1+1.35*x2)-12
    O_d_ion_log = np.log10(THREE_BETA)+6.200+1.251/t_3 - 0.55*np.log10(t_3) - 0.014*(t_3)-12

    O_s_ion = 10**(O_s_ion_log)
    O_d_ion = 10**(O_d_ion_log)
    com_O = O_s_ion + O_d_ion
    com_O_log = np.log10(com_O) +12

    return O_s_ion , O_d_ion, com_O_log, O_s_ion_log, O_d_ion_log
