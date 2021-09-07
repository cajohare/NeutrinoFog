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
ne = 50 # number of energy bins (anything >50 is accurate enough)
nm = 200 # number of mass points
n_ex = 500 # number of exposure points (I wouldn't go below 250)
ns = 500 # number of cross section points (I wouldn't go below 250)
ex_min = 1e-5 # minimum exposure
ex_max = 1e19 # maximum exposure
m_vals = logspace(log10(0.1),log10(1.0e4),nm) # mass points
#==============================================================================#
Flux_norm = NuFlux # See Params.py
Flux_err = NuUnc # See Params.py
E_th = 1.0e-4 # Threshold
E_max = 200.0 # Max recoil energy
sigma_min = 1e-50 # Minimum cross section to scan over
sigma_max = 1e-43 # Maximum cross section to scan over
inp = 'NuFloor'+sys.argv[1]+'_detailed_SI' # Filename to save data to
#==============================================================================#
# This is a bit ugly but works for now
if sys.argv[1]=='Xe':
    Nuc = Xe131
elif (sys.argv[1])=='Ar':
    Nuc = Ar40
elif (sys.argv[1])=='CaWO4':
    Nucs = [Ca40,W184,O16]
elif (sys.argv[1])=='F':
    Nuc = F19
elif (sys.argv[1])=='He':
    sigma_min = 1e-48
    sigma_max = 1e-41
    Nuc = He4
elif (sys.argv[1])=='Ge':
    Nuc = Ge74
elif (sys.argv[1])=='NaI':
    Nucs = [Na23,I127]
#==============================================================================#
if (sys.argv[1])=='CaWO4':
    f0 = 40/(40+184+4*16)
    f1 = 184/(40+184+4*16)
    f2 = 4*16/(40+184+4*16)
    R_sig = f0*BinnedWIMPRate(E_th,E_max,ne,m_vals,Nucs[0],C_SI,FormFactorHelm,MeanInverseSpeed_SHM)
    R_sig += f1*BinnedWIMPRate(E_th,E_max,ne,m_vals,Nucs[1],C_SI,FormFactorHelm,MeanInverseSpeed_SHM)
    R_sig += f2*BinnedWIMPRate(E_th,E_max,ne,m_vals,Nucs[2],C_SI,FormFactorHelm,MeanInverseSpeed_SHM)
    R_nu = f0*BinnedNeutrinoRates(E_th,E_max,ne,Nucs[0],Flux_norm)
    R_nu += f1*BinnedNeutrinoRates(E_th,E_max,ne,Nucs[1],Flux_norm)
    R_nu += f2*BinnedNeutrinoRates(E_th,E_max,ne,Nucs[2],Flux_norm)
elif sys.argv[1]=='NaI':
    f0 = 22/(22+127)
    f1 = 127/(22+127)
    R_sig = f0*BinnedWIMPRate(E_th,E_max,ne,m_vals,Nucs[0],C_SI,FormFactorHelm,MeanInverseSpeed_SHM)
    R_sig += f1*BinnedWIMPRate(E_th,E_max,ne,m_vals,Nucs[1],C_SI,FormFactorHelm,MeanInverseSpeed_SHM)
    R_nu = f0*BinnedNeutrinoRates(E_th,E_max,ne,Nucs[0],Flux_norm)
    R_nu += f1*BinnedNeutrinoRates(E_th,E_max,ne,Nucs[1],Flux_norm)
else:
    R_sig = BinnedWIMPRate(E_th,E_max,ne,m_vals,Nuc,C_SI,FormFactorHelm,MeanInverseSpeed_SHM)
    R_nu = BinnedNeutrinoRates(E_th,E_max,ne,Nuc,Flux_norm)
#==============================================================================#
runDL_2D(inp,R_sig,R_nu,m_vals,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,Flux_norm,Flux_err,verbose=False)
