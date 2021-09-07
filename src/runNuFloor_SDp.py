#==============================================================================#
import sys
import os
sys.path.append('../src')
from numpy import *
from Params import *
from NeutrinoFuncs import BinnedNeutrinoRates
from WIMPFuncs import BinnedWIMPRate,MeanInverseSpeed_SHM,C_SDp
from LabFuncs import FormFactorGaussian
from Like import runDL_2D
#==============================================================================#
ne = 50
nm = 200
n_ex = 500
ns = 500
ex_min = 1e-5
ex_max = 1e19
m_vals = logspace(log10(0.1),log10(1.0e4),nm)
#==============================================================================#
Flux_norm = NuFlux
Flux_err = NuUnc
E_th = 1.0e-4
E_max = 200.0
#==============================================================================#
if sys.argv[1]=='Xe':
    Nucs = [Xe131,Xe129]
    sigma_min = 1e-43
    sigma_max = 1e-35
elif sys.argv[1]=='Ge':
    Nucs = [Ge73]
    sigma_min = 1e-43
    sigma_max = 1e-35
elif sys.argv[1]=='F':
    Nucs = [F19]
    sigma_min = 1e-47
    sigma_max = 1e-40
elif sys.argv[1]=='NaI':
    Nucs = [Na23,I127]
    sigma_min = 1e-48
    sigma_max = 1e-38
elif sys.argv[1]=='Si':
    Nucs = [Si29]
    sigma_min = 1e-41
    sigma_max = 1e-33
#==============================================================================#
if sys.argv[1]=='NaI':
    f0 = 22/(22+127)
    f1 = 127/(22+127)
    R_sig = f0*Nucs[0].IsotopicFraction*BinnedWIMPRate(E_th,E_max,ne,m_vals,Nucs[0],C_SDp,FormFactorGaussian,MeanInverseSpeed_SHM)\
            + f1*Nucs[1].IsotopicFraction*BinnedWIMPRate(E_th,E_max,ne,m_vals,Nucs[1],C_SDp,FormFactorGaussian,MeanInverseSpeed_SHM)
    R_nu = f0*BinnedNeutrinoRates(E_th,E_max,ne,Nucs[0],Flux_norm)\
            +f1*BinnedNeutrinoRates(E_th,E_max,ne,Nucs[1],Flux_norm)
else:
    R_sig = 0
    for i in range(0,len(Nucs)):
        R_sig += Nucs[i].IsotopicFraction*BinnedWIMPRate(E_th,E_max,ne,m_vals,Nucs[i],C_SDp,FormFactorGaussian,MeanInverseSpeed_SHM)
    R_nu = BinnedNeutrinoRates(E_th,E_max,ne,Nucs[0],Flux_norm)
#==============================================================================#
runDL_2D('NuFloor'+sys.argv[1]+'_detailed_SDp',R_sig,R_nu,m_vals,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,Flux_norm,Flux_err,verbose=False)
