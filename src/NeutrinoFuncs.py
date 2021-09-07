#===========================NeutrinoFuncs.py===================================#
# Created by Ciaran O'Hare 2021

# Contains functions for performing calculations of CEvNS

#==============================================================================#
from numpy import pi, sqrt, exp, zeros, size, shape, array, meshgrid, reshape
from numpy import trapz, interp, loadtxt, count_nonzero, flipud, arange, append
from numpy import loadtxt, digitize, log10, cos, sin, arccos, arctan2, count_nonzero
from numpy import logspace, linspace, ones, asarray, histogram2d,column_stack
from numpy import nan, isnan, amin, amax,argmin,argmax,cumsum, sum, around
from numpy import savetxt, histogram, minimum, tile, arcsin, squeeze
from numpy.linalg import norm
from Params import nufile_root, nufile_dir, nuname, n_Enu_vals, recoil_dir
from Params import mono, NuMaxEnergy, NuFlux, NuUnc, whichsolar, n_nu_tot
from Params import m_p_keV, c_km, seconds2year, m_p_keV,m_e,m_e_GeV
from Params import G_F_GeV, sinTheta_Wsq, N_A, Jan1, EarthRadius, eV2J
import LabFuncs
import Params

def MaxNuRecoilEnergies(Nuc): # Max recoil energies
    m_N = 0.93141941*(Nuc.MassNumber)*1.0e6
    E_r_max = 2*m_N*(1000.0*NuMaxEnergy)**2.0/(m_N+1000*NuMaxEnergy)**2.0
    return E_r_max

def GetNuFluxes(E_th,Nuc=Params.F19):
    # Reads each neutrino flux data file
    # the energies are stored in E_nu_all, fluxes in Flux_all

    # Figure out which backgrounds give recoils above E_th
    E_r_max = MaxNuRecoilEnergies(Nuc) # Max recoil energy for neutrino
    sel = range(1,n_nu_tot+1)*(E_r_max>E_th)
    sel = sel[sel!=0]-1
    n_nu = count_nonzero(E_r_max>E_th)
    E_nu_all = zeros(shape=(n_nu,n_Enu_vals))
    Flux_all = zeros(shape=(n_nu,n_Enu_vals))
    Flux_err = zeros(shape=(n_nu))
    Flux_norm = zeros(shape=(n_nu))
    Solar = zeros(n_nu,dtype=bool)
    Names = asarray([nuname[i] for i in sel])

    ii = 0
    for s in sel:
        if mono[s]:
            E_nu_all[ii,0] = NuMaxEnergy[s]
            Flux_all[ii,0] = NuFlux[s]
        else:
            data = loadtxt(nufile_dir+'normalised/'+nuname[s]+nufile_root,delimiter=',')
            E_nu_all[ii,:],Flux_all[ii,:] = data[:,0],data[:,1]
            Flux_all[ii,:] = Flux_all[ii,:]*NuFlux[s]

        Flux_norm[ii] = NuFlux[s]
        Flux_err[ii] = NuUnc[s] # Select rate normalisation uncertainties
        Solar[ii] = whichsolar[s]
        ii = ii+1
    return Names,Solar,E_nu_all,Flux_all,Flux_norm,Flux_err

def All_dRdE(E_r,t,Solar,E_nu_all,Flux_all,Nuc=Params.F19): # Time-Energy
    n_nu = shape(Flux_all)[0]
    ne = size(E_r)
    dR = zeros(shape=(n_nu,ne))
    for i in range(0,n_nu):
        dR[i,:] = dRdE_nu(E_r,t,\
                            Solar[i],E_nu_all[i,:],Flux_all[i,:],Nuc=Nuc)
    return dR

# CEvNS event rate vs recoil energy E_r
def dRdE(E_r,t,sol,E_nu,Flux,Nuc=Params.F19):
    N = Nuc.NumberOfNeutrons
    Z = Nuc.NumberOfProtons
    Q_W = 1.0*N-(1-4.0*sinTheta_Wsq)*Z # weak nuclear hypercharge
    m_N_GeV = 0.93141941*(N+Z) # nucleus mass in GeV
    m_N_keV = m_N_GeV*1.0e6 # nucleus mass in keV

    dRdE = zeros(shape=shape(E_r))
    FF = LabFuncs.FormFactorHelm(E_r,N+Z)**2.0
    ne = size(E_r)

    if Flux[1]>0.0:
        for i in range(0,ne):
            diff_sigma = (G_F_GeV**2.0 /(4.0*pi))*(Q_W**2.0)*m_N_GeV*(1.0 \
                        -(m_N_keV*E_r[i])/(2.0*(E_nu*1000.0)**2.0))*\
                        (0.197e-13)**2.0*(1.0e-6)*1000.0/(1.0*N+1.0*Z)*(N_A)
            diff_sigma[diff_sigma<0.0] = 0.0
            dRdE[i] = trapz(diff_sigma*Flux*FF[i],x=E_nu)
    else:
        for i in range(0,ne):
            diff_sigma = (G_F_GeV**2.0 /(4.0*pi))*(Q_W**2.0)*m_N_GeV*(1.0 \
                        -(m_N_keV*E_r[i])/(2.0*(E_nu[0]*1000.0)**2.0))*\
                        (0.197e-13)**2.0*(1.0e-6)*1000.0/(1.0*N+1.0*Z)*(N_A)
            if diff_sigma>0: # for monochromatic nu's
                dRdE[i] = diff_sigma*Flux[0]*E_nu[0]*FF[i]

    if sol:
        fMod = LabFuncs.EarthSunDistanceMod(t)
    else:
        fMod = 1.0

    # Convert into /ton/year/keV
    dRdE = fMod*dRdE*1000*seconds2year
    return dRdE

def R_Indiv(s,E_th,E_max,Nuc=Params.F19,f_eff=None,f_eres=None):
    nfine = 1000
    Efine = logspace(-3.0,log10(200.0),nfine)
    dR = dRdE_Indiv(s,Efine,Nuc=Nuc,f_eff=f_eff,f_eres=f_eres)
    mask = (Efine<E_max)&(Efine>E_th)
    R = trapz(dR[mask],Efine[mask])
    return R

def dRdE_Indiv(s,E_r,Nuc=Params.F19,f_eff=None,f_eres=None):
    # Load nu flux
    nfine = 1000
    Efine = logspace(-3.0,log10(200.0),nfine)

    data = loadtxt(nufile_dir+'normalised/'+nuname[s]+nufile_root,delimiter=',')
    E_nu = data[:,0]
    Flux = NuFlux[s]*data[:,1]

    sol = whichsolar[s]
    dR = dRdE(Efine,Jan1*ones(shape=nfine),sol,E_nu,Flux,Nuc=Nuc)

    if f_eff != None: dR *= LabFuncs.f_eff(Efine) # Correct for efficiency

    if f_eres != None: dR = LabFuncs.SmearE(Efine,dR,f_eres(Efine)) # Smear by energy resolution

    dR = interp(E_r,Efine,dR)
    return dR

# Specific examples, for quick access
def R_hep(E_th,E_max,Nuc=Params.F19,f_eff=None,f_eres=None):
    return R_Indiv(2,E_th,E_max,Nuc=Nuc,f_eff=f_eff,f_eres=f_eres)

def R_8B(E_th,E_max,Nuc=Params.F19,f_eff=None,f_eres=None):
    return R_Indiv(5,E_th,E_max,Nuc=Nuc,f_eff=f_eff,f_eres=f_eres)

def R_AtmNu(E_th,E_max,Nuc=Params.F19,f_eff=None,f_eres=None):
    return R_Indiv(10,E_th,E_max,Nuc=Nuc,f_eff=f_eff,f_eres=f_eres)

def R_DSNB(E_th,E_max,Nuc=Params.F19,f_eff=None,f_eres=None):
    return R_Indiv(9,E_th,E_max,Nuc=Nuc,f_eff=f_eff,f_eres=f_eres)

def dRdE_hep(E_r,Nuc=Params.F19,f_eff=None,f_eres=None):
    return dRdE_Indiv(2,E_r,Nuc=Nuc,f_eff=f_eff,f_eres=f_eres)

def dRdE_8B(E_r,Nuc=Params.F19,f_eff=None,f_eres=None):
    return dRdE_Indiv(5,E_r,Nuc=Nuc,f_eff=f_eff,f_eres=f_eres)

def dRdE_AtmNu(E_r,Nuc=Params.F19,f_eff=None,f_eres=None):
    return dRdE_Indiv(10,E_r,Nuc=Nuc,f_eff=f_eff,f_eres=f_eres)

def dRdE_DSNB(E_r,Nuc=Params.F19,f_eff=None,f_eres=None):
    return dRdE_Indiv(9,E_r,Nuc=Nuc,f_eff=f_eff,f_eres=f_eres)

def dRdE_TotalSolar(E_r,Nuc=Params.F19,f_eff=None,f_eres=None):
    dR = 0
    for s in range(0,9):
        dR += dRdE_Indiv(s,E_r,Nuc=Nuc,f_eff=f_eff,f_eres=f_eres)
    return dR

def BinnedNeutrinoRates(E_th,E_max,ne,Nuc,Flux_norm=NuFlux):
    E_be = logspace(log10(E_th),log10(E_max),ne+1)
    E_r_max = MaxNuRecoilEnergies(Nuc)
    # Load neutrinos
    Names,Solar,E_nu_all,Flux_all,Flux_norm_default,errs = GetNuFluxes(E_th)
    n_nu = shape(E_nu_all)[0]
    R = zeros((n_nu,ne))
    for i in range(0,n_nu):
        Efine = logspace(log10(E_th),log10(E_r_max[i]),1000)
        R_tot = trapz(dRdE(Efine,0,Solar[i],E_nu_all[i,:],Flux_all[i,:],Nuc),Efine)

        dR = dRdE(E_be,0,Solar[i],E_nu_all[i,:],Flux_all[i,:],Nuc)
        R[i,:] = 0.5*(E_be[1:]-E_be[0:-1])*(dR[1:]+dR[0:-1])
        R[i,:] = (Flux_norm[i]/Flux_norm_default[i])*R_tot*R[i,:]/sum(R[i,:])
    return R

#===============================Reactor neutrinos data=========================#
import pandas as pd
def ReactorFlux(E_nu,Loc):
    # from https://arxiv.org/abs/1101.2663 TABLE VI.
    U235_c = array([3.217,-3.111,1.395,-0.3690,0.04445,-0.002053])
    U238_c = array([0.4833,0.1927,-0.1283,-0.006762,0.00233,-0.0001536])
    P239_c = array([6.413,-7.432,3.535,-0.8820,0.1025,-0.00455])
    P241_c = array([3.251,-3.204,1.428,-0.3675,0.04245,-0.001896])

    def Sk(coeff,E_nu):
        Sk = exp(coeff[0]
                 +coeff[1]*E_nu**1.0
                 +coeff[2]*E_nu**2.0
                 +coeff[3]*E_nu**3.0
                 +coeff[4]*E_nu**4.0
                 +coeff[5]*E_nu**5.0)
        return Sk

    # Get powers of nearby reactors and coordinates
    df = pd.read_excel(nufile_dir+'reactor/Reactors_2019.xlsx',header=None)
    Powers = df.loc[:,6].values*1.0
    Coords = df.loc[:,2:3].values
    eff1 = sum(df.loc[:,7:].values,1)/(100*12.0)
    Powers *= eff1

    # Calculate distance
    phi1 = Loc.Latitude*pi/180.0
    phi2 = Coords[:,0]*pi/180.0
    dphi = phi1-phi2
    dlambda = (Loc.Longitude-Coords[:,1])*pi/180.0
    a = sin(dphi/2)**2.0 + cos(phi1)*cos(phi2)*sin(dlambda/2)**2.0
    bearing = 2*arctan2(sqrt(a),sqrt(1-a))
    Distances = 2*EarthRadius*sin(bearing/2.0)

    # Fission fraction and average energies from https://arxiv.org/abs/1101.2663
    # k = 0,3 for U235,U238,U239,P241
    fk = array([0.58,0.07,0.3,0.05])
    Ek = array([202.36,205.99,211.12,214.26])
    phi = fk[0]*Sk(U235_c,E_nu)
    phi += fk[1]*Sk(U238_c,E_nu)
    phi += fk[2]*Sk(P239_c,E_nu)
    phi += fk[3]*Sk(P241_c,E_nu)

    # Calculate flux by multiplying spectrum
    # by Emission rate*power/(4*pi*distance^2)
    # Where Emission rate = (6.0*0.75/sum(fk*Ek)) accounting for operating efficiency
    Phi = 0.0
    #eff = 0.75
    for i in range(0,size(Distances)):
        Phi += phi*(Powers[i]/eV2J)*(6.0/sum(fk*Ek))\
                /(4*pi*(Distances[i]*100)**2.0)
    return Phi
