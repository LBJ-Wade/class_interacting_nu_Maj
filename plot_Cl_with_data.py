
# coding: utf-8

# In[ ]:

# import necessary modules
#get_ipython().magic(u'matplotlib inline')
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.optimize import fsolve
import math
from matplotlib import rc
import matplotlib.patches as patches

from scipy.interpolate import interp1d
from matplotlib.ticker import FixedLocator
from math import floor
from mpl_toolkits.axes_grid1 import make_axes_locatable

# In[ ]:

# esthetic definitions for the plots

rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
# matplotlib.rc('font', **font)
matplotlib.mathtext.rcParams['legend.fontsize']='medium'
plt.rcParams["figure.figsize"] = [8.0,6.0]


l_TT_low,Dl_TT_low,err_TT_low= np.loadtxt("/Users/poulin/Dropbox/Labo/ProgrammeDarkAges/error_Planck/COM_PowerSpect_CMB-TT-loL-full_R2.02.txt",unpack=True,usecols=(0,1,2))
l_TE_low,Dl_TE_low,err_TE_low= np.loadtxt("/Users/poulin/Dropbox/Labo/ProgrammeDarkAges/error_Planck/COM_PowerSpect_CMB-TE-loL-full_R2.02.txt",unpack=True,usecols=(0,1,2))
l_TT_high,Dl_TT_high,err_TT_high= np.loadtxt("/Users/poulin/Dropbox/Labo/ProgrammeDarkAges/error_Planck/COM_PowerSpect_CMB-TT-hiL-binned_R2.02.txt",unpack=True,usecols=(0,3,4))
l_TE_high,Dl_TE_high,err_TE_high= np.loadtxt("/Users/poulin/Dropbox/Labo/ProgrammeDarkAges/error_Planck/COM_PowerSpect_CMB-TE-hiL-binned_R2.02.txt",unpack=True,usecols=(0,3,4))
l_EE_low,Dl_EE_low,err_EE_low= np.loadtxt("/Users/poulin/Dropbox/Labo/ProgrammeDarkAges/error_Planck/COM_PowerSpect_CMB-EE-loL-full_R2.02.txt",unpack=True,usecols=(0,1,2))
l_EE_high,Dl_EE_high,err_EE_high= np.loadtxt("/Users/poulin/Dropbox/Labo/ProgrammeDarkAges/error_Planck/COM_PowerSpect_CMB-EE-hiL-binned_R2.02.txt",unpack=True,usecols=(0,3,4))
lmin_phiphi,lmax_phiphi,cl_phiphi,err_phiphi= np.loadtxt("/Users/poulin/Dropbox/Labo/ProgrammeDarkAges/error_Planck/agressive_lensing.csv",unpack=True)

# In[ ]:
H0_values = [40,50,60,67.5,74,80,90,100]
##LCDM BESTFIT Planck 2018####
common_settings = {'output':'tCl,pCl,lCl,mPk',
                   'lensing':'yes',
                   'format':'camb',
                   'N_ncdm': 1,
                   'deg_ncdm': 1,
                   'm_ncdm': 0.06,
                   'T_ncdm': 0.71611,
                   'N_ur':  2.0328,
                   'P_k_max_1/Mpc':3.0,
                   'l_switch_limber':9,
                   'input_verbose':1,
                   'background_verbose':1,
                   'thermodynamics_verbose':1,
                   'perturbations_verbose':1,
                   'transfer_verbose':1,
                   'primordial_verbose':1,
                   'spectra_verbose':1,
                   'nonlinear_verbose':1,
                   'lensing_verbose':1,
                   'output_verbose':1}

M = Class()

T_cmb = 2.7225
factor1 = l_TT_high*(l_TT_high+1)/2./np.pi;
conversion1 = 1/(factor1*(T_cmb*1.e6)**2)
factor2 = l_TT_low*(l_TT_low+1)/2./np.pi;
conversion2 = 1/(factor2*(T_cmb*1.e6)**2)
plt.errorbar(l_TT_high, Dl_TT_high*conversion1, yerr=err_TT_high*conversion1, fmt='.')
plt.errorbar(l_TT_low, Dl_TT_low*conversion2, yerr=err_TT_low*conversion2, fmt='.')
for H0 in H0_values:
    M.set(common_settings)
    M.set({
    'omega_b':0.022383,
    'omega_cdm':0.12011,
    'ln10^{10}A_s':3.0448,
    'n_s':0.96605,
    'tau_reio':0.0543,
    })
    M.set({'H0':H0})
    M.compute()
    clM = M.lensed_cl(2500)
    ll_LCDM = clM['ell'][2:]
    clTT_LCDM = clM['tt'][2:]
    plt.plot(ll_LCDM,clTT_LCDM,lw=2,label=r'$H_0 = %.2f $km/s/Mpc'%(H0))
    plt.savefig('H0_%.2f.pdf'%(H0), bbox_inches='tight')
    plt.set_xlabel(r'$\ell$',fontsize=20)
    plt.set_ylabel(r'$C_\ell^\mathrm{TT}$',fontsize=20)
    plt.legend(frameon=False,prop={'size':12},loc='upper left',borderaxespad=0.)
    M.struct_cleanup()


# In[ ]:




# In[ ]:
