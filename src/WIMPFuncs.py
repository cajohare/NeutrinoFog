#===============================WIMPFuncs.py===================================#
# Created by Ciaran O'Hare 2020

# Contains all the functions for doing the WIMPy calculations

#==============================================================================#

import numpy as np
from numpy import pi, sqrt, exp, zeros, size, shape, array, linspace, logspace
from numpy import cos, sin, arctan2, arccos, trapz, ones, log10, ndim, meshgrid
from numpy import nan, isnan, column_stack, amin, amax, zeros_like
from numpy.linalg import norm
from scipy.special import erf
import Params
from Params import m_p_keV, c_km, seconds2year, m_p_kg, GeV_2_kg, c_cm, Jan1

#==============================================================================#
#-------------------- Energy-Time dependent recoil rate------------------------#

#---------------------------------- v_min -------------------------------------#
def MinimumWIMPSpeed(E_r,A,m_chi,delta=0):
    # E_r = recoil energy in keVr
    # A = nucleus mass number
    # m_chi = Wimp mass in GeV
    # delta = for inelastic scattering
    mu_p = 1.0e6*m_chi*m_p_keV/(1.0e6*m_chi + m_p_keV) # reduced proton mass
    m_N_keV = A*m_p_keV # nucleus mass in keV
    mu_N_keV = 1.0e6*m_chi*m_N_keV/(1.0e6*m_chi + m_N_keV) # reduced nucleus mass
    v_min = sqrt(1.0/(2*m_N_keV*E_r))*(m_N_keV*E_r/mu_N_keV + delta)*c_km
    return v_min

#---------------------------------- E_max -------------------------------------#
def MaxWIMPEnergy(A,m_chi,v_lab=245.6,v_esc=533.0):
    # A = nucleus mass number
    # v_lab = Lab velocity in km/s
    # m_chi = Wimp mass in GeV
    # v_esc = Escape speed in km/s
    m_N = m_p_keV*A
    mu_N = 1.0e6*m_N*m_chi/(1.0e6*m_chi+m_N)
    E_max_lim = 2.0*mu_N*mu_N*2.0*((v_esc+v_lab)/c_km)**2.0/m_N
    return E_max_lim

def MinimumWIMPSpeed(E_r,A,m_chi,delta=0):
    # E_r = recoil energy in keVr
    # A = nucleus mass number
    # m_chi = WIMP mass in GeV
    # delta = for inelastic scattering
    m_N_keV = A*m_p_keV # nucleus mass in keV
    mu_N_keV = 1.0e6*m_chi*m_N_keV/(1.0e6*m_chi + m_N_keV) # reduced nucleus mass
    v_min = sqrt(1.0/(2*m_N_keV*E_r))*(m_N_keV*E_r/mu_N_keV + delta)*c_km
    return v_min

#-------------------- Mean Inverse Speed (for Gaussian f(v)) --------------------------#
def MeanInverseSpeed_SHM(v_min,sig_v=167.0,v_esc=533.0,v_lab=245.6):
    N_esc = erf(v_esc/(sqrt(2.0)*sig_v))\
            -sqrt(2.0/pi)*(v_esc/sig_v)*exp(-v_esc**2.0/(2.0*sig_v**2.0))

    # Define:
    v_0 = sig_v*sqrt(2.0)
    x = v_min/v_0
    z = v_esc/v_0
    y = v_lab/v_0

    # Set up conditional terms
    g = zeros_like(v_min)
    g[(x<abs(y-z))&(z<y)] = (1.0/(v_0*y))
    g2 = (1.0/(2.0*N_esc*v_0*y))*(erf(x+y)-erf(x-y)-(4.0/sqrt(pi))*y*exp(-z**2))
    g3 = (1.0/(2.0*N_esc*v_0*y))*(erf(z)-erf(x-y)-(2.0/sqrt(pi))*(y+z-x)*exp(-z**2))

    # Apply conditions
    g[(x<abs(y-z))&(z>y)] = g2[(x<abs(y-z))&(z>y)]
    g[(abs(y-z)<x)&(x<(y+z))] = g3[(abs(y-z)<x)&(x<(y+z))]
    g[(y+z)<x] = 0.0

    return g

#--------------------Helm Form Factor-------------------------------------------#
def C_SI(Nuc,):
    A = Nuc.MassNumber
    return A**2

def C_SDp(Nuc):
    S_p = Nuc.ExpProtonSpin
    J = Nuc.NuclearSpin
    return (4/3)*((J+1)/J)*(S_p)**2

def C_SDn(Nuc):
    S_n = Nuc.ExpNeutronSpin
    J = Nuc.NuclearSpin
    return (4/3)*((J+1)/J)*(S_n)**2

def C_SDpn(Nuc):
    S_p = Nuc.ExpProtonSpin
    S_n = Nuc.ExpNeutronSpin
    J = Nuc.NuclearSpin
    return (4/3)*((J+1)/J)*(S_p**2+S_n**2)

def dRdE(E_r,m_chi,sigma_p,Nuc,NuclearEnhancementFactor,FormFactor,gvmin,rho_0=0.3):
    '''
    * Spin independent differentual recoil rate that takes in recoil energy in
    units of keVr and a proton cross section in units of cm^2 and outputs a rate
    in units of (ton year keVr)^-1

    * gvmin_function should be a function that takes in v_min in (km/s) and outputs
    g(v_min) in units of (km/s)^-1
    '''
    A = Nuc.MassNumber

    C = NuclearEnhancementFactor(Nuc)

    # DM-proton reduced mass (in units of keV)
    mu_p = 1.0e6*m_chi*m_p_keV/(1.0e6*m_chi + m_p_keV)

    # Rate constants (in units cm kg^-1 s^-2)
    R0 = (c_cm**2)*((rho_0*1.0e6*C*sigma_p)/(2*m_chi*GeV_2_kg*mu_p**2))

    # Mean inverse speed
    v_min = MinimumWIMPSpeed(E_r,A,m_chi)
    g = gvmin(v_min)/(1000.0*100.0) # convert to cm^-1 s

    # Compute rate = (Rate amplitude * gmin * form factor)
    FF = FormFactor(E_r,A)**2.0
    dR = R0*g*FF
    dR = dR*seconds2year*1000.0 # convert to (ton-year-keV)^-1
    return dR


def BinnedWIMPRate(E_th,E_max,ne,m_vals,Nuc,NuclearEnhancementFactor,FormFactor,gvmin,**kwargs):
    nm = size(m_vals)

    E_be = logspace(log10(E_th),log10(E_max),ne+1)
    R = zeros(shape=(nm,ne))
    for i in range(0,nm):
        E_r_max = MaxWIMPEnergy(Nuc.MassNumber,m_vals[i],**kwargs)
        Efine = logspace(log10(E_th),log10(E_r_max),1000)
        R_tot = trapz(dRdE(Efine,m_vals[i],1.0e-45,Nuc,NuclearEnhancementFactor,FormFactor,gvmin,**kwargs),Efine)

        dR = dRdE(E_be,m_vals[i],1.0e-45,Nuc,NuclearEnhancementFactor,FormFactor,gvmin,**kwargs)
        R[i,:] = 0.5*(E_be[1:]-E_be[0:-1])*(dR[1:]+dR[0:-1])
        R[i,:] = R_tot*R[i,:]/sum(R[i,:])
    return R
