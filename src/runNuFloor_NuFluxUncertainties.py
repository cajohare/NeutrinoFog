#==============================================================================#
import sys
import os
sys.path.append('../src')
from numpy import *
from Params import *
from NeutrinoFuncs import BinnedNeutrinoRates
from WIMPFuncs import BinnedWIMPRate,MeanInverseSpeed_SHM,C_SI
from LabFuncs import FormFactorHelm
from Like import runDL_2D
#==============================================================================#
# Input
Nuc = Xe131
ex_min=1.0e-5
ex_max=1.0e4
sigma_min = 1e-50
sigma_max = 1e-43
#==============================================================================#
# Resolution
ne = 50
nm = 200
n_ex = 500
ns = 500
m_vals = logspace(log10(0.1),log10(1.0e4),nm)
#==============================================================================#
# Parameters
E_th = 1.0e-4
E_max = 200.0
m_vals = logspace(log10(0.1),log10(1.0e4),nm)
#==============================================================================#
# Make data
R_sig = BinnedWIMPRate(E_th,E_max,ne,m_vals,Nuc,C_SI,FormFactorHelm,MeanInverseSpeed_SHM)
#==============================================================================#
# Run limits
Flux_norm = NuFlux
#==============================================================================#
Flux_err = NuUnc*1
Flux_err[0:9] /= 10
R_nu = BinnedNeutrinoRates(E_th,E_max,ne,Nuc,Flux_norm)
runDL_2D('NuFloorXe_TenthSolarUncertainty',R_sig,R_nu,m_vals,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,Flux_norm,Flux_err,verbose=False)
#==============================================================================#
Flux_err = NuUnc*1
Flux_err[10] /= 10
R_nu = BinnedNeutrinoRates(E_th,E_max,ne,Nuc,Flux_norm)
runDL_2D('NuFloorXe_TenthAtmUncertainty',R_sig,R_nu,m_vals,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,Flux_norm,Flux_err,verbose=False)
#==============================================================================#
