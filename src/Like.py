#====================================Like.py===================================#
# Created by Ciaran O'Hare 2021

# Contains functions for interfacing with the fortran code in src/like
# the fortran likelihood code needs to be compiled first by running the make
# file in src/like
#==============================================================================#


from __future__ import print_function
from numpy import pi, sqrt, exp, zeros, size, shape, array, append, flipud, gradient
from numpy import trapz, interp, loadtxt, log10, log, savetxt, vstack, transpose
from numpy import ravel,tile,mean,inf,nan,amin,amax
from scipy.ndimage.filters import gaussian_filter1d
from scipy.integrate import cumtrapz
from numpy.linalg import norm
from scipy.special import gammaln
from Params import *
import LabFuncs
import NeutrinoFuncs
import WIMPFuncs
import shlex
import subprocess
import pprint

def Floor_2D(data,filt=True,filt_width=3,Ex_crit=1e10):
    sig = data[1:,0]
    m = data[0,1:]
    n = size(m)
    ns = size(sig)
    Ex = flipud(transpose(data[1:,1:].T))
    Ex[Ex>Ex_crit] = nan
    Exmin = amin(Ex[Ex>0])
    Ex[Ex==0] = Exmin
    DY = zeros(shape=shape(Ex))
    for j in range(0,n):
        y = log10(Ex[:,j])
        if filt:
            y = gaussian_filter1d(gaussian_filter1d(y,sigma=3),filt_width)
            dy = gradient(y,log10(sig[2])-log10(sig[1]))
            dy = gaussian_filter1d(dy,filt_width)
        else:
            dy = gradient(y,log10(sig[2])-log10(sig[1]))

        DY[:,j] = dy
    NUFLOOR = zeros(shape=n)
    #for j in range(0,n):
    #    DY[:,j] = gaussian_filter1d(DY[:,j],filt_width)
    for j in range(0,n):
        for i in range(0,ns):
            if DY[ns-1-i,j]<=-2.0:
                i0 = ns-1-i
                i1 = i0+10
                NUFLOOR[j] = 10.0**interp(-2,DY[i0:i1+1,j],log10(sig[i0:i1+1]))
                DY[ns-1-i:-1,j] = nan
                break
    DY = -DY
    DY[DY<2] = 2
    return m,sig,NUFLOOR,DY

def NuFloor_1event(mvals,Nuc,nths=100):
    # Load neutrino fluxes
    Names,solar,E_nu_all,Flux_all,Flux_norm,Flux_err = NeutrinoFuncs.GetNuFluxes(0.0)
    n_nu = shape(Flux_all)[0]
    E_ths = logspace(log10(0.0001),log10(100.0),nths)
    t = 0

    R = zeros(shape=nths)
    for i in range(0,n_nu):
        R = R+NeutrinoFuncs.dRdE(E_ths,t,solar[i],E_nu_all[i,:],Flux_all[i,:],Nuc)
    cumR = flipud(cumtrapz(flipud(E_ths),flipud(R)))
    cumR = append(cumR,cumR[-1])
    Exposures = 1.0/cumR

    nm = size(mvals)

    DL = zeros(shape=(nm,nths))
    for j in range(0,nths-10):
        Evals = logspace(log10(E_ths[j]),log10(1000.0),200)
        for i in range(0,nm):
            m = mvals[i]
            Nw = Exposures[j]*trapz(WIMPFuncs.dRdE(Evals,m,1.0e-45,Nuc,\
                    WIMPFuncs.C_SI,LabFuncs.FormFactorHelm,WIMPFuncs.MeanInverseSpeed_SHM),Evals)
            if Nw>0:
                DL[i,j] = 2.3*1.0e-45/Nw
    DL[DL<0] = inf
    DL[DL==0] = inf

    nu1 = amin(DL,1)
    return nu1

#==============================================================================#
# Both of these functions save WIMP/neutrino data in a format that can be then
# read by the fortran code
def SaveWIMPData(inp,R_sig,m_vals):
    nTot_bins = shape(R_sig)[1]
    nm = shape(R_sig)[0]
    hdr1 = str(nm)+' '+str(nTot_bins)
    dat1 = zeros(shape=(nm,nTot_bins+1))
    dat1[:,1:] = R_sig
    dat1[:,0] = m_vals
    savetxt(recoil_dir+'RD_sig_'+inp+'.txt',dat1,header=hdr1)
    return

def SaveNuData(inp,R_nu,Flux_norm,Flux_err):
    nTot_bins = shape(R_nu)[1]
    n_nu = shape(R_nu)[0]
    hdr2 = str(n_nu)+' '+str(nTot_bins)
    dat2 = zeros(shape=(n_nu,nTot_bins+2))
    dat2[:,2:] = R_nu
    dat2[:,0] = Flux_norm
    dat2[:,1] = Flux_err
    savetxt(recoil_dir+'RD_bg_'+inp+'.txt',dat2,header=hdr2)
    return
#==============================================================================#
# These are functions that call the compiled fortran code from python.
def runDL_fort(inp,ex_min=1.0e-1,ex_max=1.0e7,n_ex=9,\
                  verbose=False):
    savetxt(recoil_dir+'Ex_'+inp+'.txt',array([[ex_min],[ex_max],[n_ex]]))
    command = "../src/like/./runDL "+inp
    if verbose:
        command += " 1"

    process = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE)
    while True:
        output = process.stdout.readline()
        if process.poll()==0:
            break
        if output:
            print(output.strip().decode("utf-8"))
    rc = process.poll()
    return rc

def runDL_2D_fort(inp,sigma_min=1e-50,sigma_max=1e-41,ns=200,ex_min=1.0e-1,ex_max=1.0e7,n_ex=100,\
                  verbose=False):
    savetxt(recoil_dir+'Ex_'+inp+'.txt',array([[ex_min],[ex_max],[n_ex]]))
    savetxt(recoil_dir+'Sig_'+inp+'.txt',array([[sigma_min],[sigma_max],[ns]]))

    command = "../src/like/./runDL_2D "+inp
    if verbose:
        command += " 1"

    process = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE)
    while True:
        output = process.stdout.readline()
        if process.poll()==0:
            break
        if output:
            print(output.strip().decode("utf-8"))
    rc = process.poll()
    return rc

def runDL(inp,R_sig,R_nu,m_vals,ex_min,ex_max,n_ex,Flux_norm,Flux_err,verbose=True):
    SaveWIMPData(inp,R_sig,m_vals)
    SaveNuData(inp,R_nu,Flux_norm,Flux_err)
    rc = runDL_fort(inp,ex_min=ex_min,ex_max=ex_max,n_ex=n_ex,verbose=verbose)
    return

def runDL_2D(inp,R_sig,R_nu,m_vals,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,Flux_norm,Flux_err,verbose=True):
    SaveWIMPData(inp,R_sig,m_vals)
    SaveNuData(inp,R_nu,Flux_norm,Flux_err)
    rc = runDL_2D_fort(inp,ex_min=ex_min,ex_max=ex_max,n_ex=n_ex,sigma_min=sigma_min,sigma_max=sigma_max,ns=ns,verbose=verbose)
    return








def lnPF(Nob,Nex): # SUM OF LOG(POISSON PDF)
    # in principle there should be a log gamma here
    # (or a factorial if using real data)
    # but it always cancels in the likelihood ratio
    # so it's commented out for speed
    L = sum(Nob*log(Nex) - Nex) #- gammaln(Nob+1.0))
    return L

def lnChi2(Nob,Nex): # SUM OF LOG(POISSON PDF)
    L = -0.5*sum((Nob-Nex)**2.0/Nex)
    return L

def lnGF(x,mu,sig): # SUM OF LOG(GAUSSIAN PDF)
    L = sum(-1.0*log(sig)-0.5*log(2.0*pi)-(x-mu)**2.0/(2.0*sig**2.0))
    return L


# def DL_gradient(sig,Ex_vals,sig_f,filt=True,filt_width=3):
#     y = log10(sig)
#     yc = (y[1:]+y[0:-1])/2
#     dEx = log10(Ex_vals[1:])-log10(Ex_vals[0:-1])
#     dsig = y[1:]-y[0:-1]
#     if filt:
#         dy = gaussian_filter1d(dEx/dsig,filt_width)
#     else:
#         dy = dEx/dsig
#     dy_f = interp(log10(sig_f),flipud(yc),flipud(dy))
#     dy_f[sig_f<amin(10.0**y)] = -2.5
#     dy_f[sig_f>amax(10.0**y)] = nan
#     return dy_f
#
# def MakeNuFloor_2D(data,filt=True,filt_width=2,ns=400,sigma_min=1e-50,sigma_max=1e-41):
#     sig_f = logspace(log10(sigma_min),log10(sigma_max),ns)
#     sig = data[1:,1:]
#     sig[sig==0] = sigma_max
#     m = data[0,1:]
#     nm = size(m)
#     Ex_vals = data[1:,0]
#
#     dy = zeros((ns,nm))
#
#     for i in range(0,nm):
#         if filt:
#             sig_i = 10.0**gaussian_filter1d(log10(sig[:,i]),filt_width)
#         else:
#             sig_i = sig[:,i]
#         sig_i[0] = sig[0,i]
#         sig_i[-1] = sig[-1,i]
#         dy[:,i] = DL_gradient(sig_i,Ex_vals,sig_f,filt=filt,filt_width=filt_width)
#     return m,sig_f,dy
