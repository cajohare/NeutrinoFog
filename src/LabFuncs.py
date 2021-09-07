#================================LabFuncs.py===================================#
# Created by Ciaran O'Hare 2021

# Description:
# Contains an assortment of functions that are all related to the 'Lab' somehow
# e.g. the nuclear form factor, lab velocity etc.


# Contains:
#####
# Efficiency functions for Ar and Xe
# Energy resolution functions for Ar and Xe
# FormFactorHelm: Only Form factor being used atm
#####

##### Resolutions
# Smear: Applies angular resolution to a recoil map as a function of direction
# SmearE: Applies energy resolution to a recoil spectrum as a function of energy
#####


##### Lab velocity
# LabVelocity: Full lab velocity in (N,W,Z) with Earth rotation
# LabVelocitySimple: Simplified Lab velocity in galactic coordinates
# JulianDay: JulianDay at dd-mm-yyyy hh:hh
# EarthVelocity: Earth velocity to second order in eccentricity
# EarthVector: Earth radius vector to second order in eccentricity
#####

##### Solar direction:
# EarthSunDistance: Distance between Earth and Sun as a function of time
# SolarDirection: Direction of the sun at a given time
#####

##### Co-ordinate transformations
# eqt2lab: Equatorial system to laboratory system
# gal2eqt: Galactic system to equatorial system
# gal2lab: Galactic system to lab system
#####
#==============================================================================#

import numpy as np
from numpy import cos, sin, pi, floor, exp, sqrt, size, zeros, shape, arccos
from numpy import array, trapz, arctan2, sign, histogram2d
from numpy import random, percentile, loadtxt, savetxt, interp, flipud
from scipy.stats import truncexpon
from scipy.special import erf
import Params
from Params import Jan1,AstronomicalUnit,EarthRadius,Msun,bigG

#==============================Form Factors====================================#
def FormFactorHelm(E_r,A):
    q = sqrt(2*A*931.5*1000*E_r)*1.0e-12/1.97e-7
    c1 = 1.23*A**(1.0/3.0)-0.6
    s = 0.9
    R_1 = sqrt(c1**2 + (7.0/3.0)*pi**2.0*(0.52**2.0) - 5*s**2.0)
    F = (3*(sin(q*R_1) - q*R_1*cos(q*R_1))*exp(-q*q*s*s/2.0)/(q*R_1)**3)
    F[q==0.0] = 1.0
    return F

def FormFactorGaussian(E_r,A):
    q = sqrt(2*A*931.5*1000*E_r)*1.0e-12/1.97e-7 # q = sqrt(2 m_N E_r)
    R = 0.92*A**(1/3)+2.68-0.78*sqrt((A**(1/3) - 3.8)**2 + 0.2)
    F = exp(-q*R/2)
    return F

#==============================Lab Velocity====================================#
# Peculiar velocity
v_pec = Params.SHMpp.PeculiarVelocity

# Earth orbital params
vv_earthrev = 29.79
eccentricity = 0.016722
eccentricity_deg = 0.9574
orb_long_ecliptic = 13.0+1.0
lat_ecl_gal = np.array([-5.5303,59.575,29.812])
long_ecl_gal = np.array([266.141,-13.3485,179.3212])
e1 = array([0.9941,0.1088,0.0042])
e2 = array([-0.0504,0.4946,-0.8677])
w_p = 2*pi/365 # orbital freq.
t1 = 79
ve = 29.79 # Earth's revolution
vrot = 0.47 # Earth's rotation

#------------------------------------------------------------------------------#
# Simple LabVelocity outputs in Galactic coordinates
def LabVelocitySimple(day,v_LSR=233.0):
    # day measured from Jan1
    vsun = array([0.0,v_LSR,0.0])+v_pec
    v_lab = vsun + EarthVelocity(day)
    return v_lab

# Only use the longer LabVelocity if the transformation into the Lab coordinate
# system is needed.
def LabVelocity(JD, Loc=Params.GranSasso, v_LSR=233.0):
    lat = Loc.Latitude
    lon = Loc.Longitude

    # Convert day into phase of Earth rotation t_lab
    UT = 24*(JD+0.5-floor(JD+0.5)) #Universal time
    MJD = JD - 2400000.5 #Modified Julian Day
    T_0 = (floor(MJD)-55197.5)/36525.0
    t_GAST = (101.0308 + 36000.770*T_0 + 15.04107*UT)/15.0
    t_lab = t_GAST + lon/15
    t_lab = 15*t_lab #Lab time in degrees


    # Galactic (LSR) Rotation
    vtemp = np.array([0.0,v_LSR,0.0])
    v_galrot = gal2lab(vtemp,t_lab, lat) #transform to lab co-ords

    # Peculiar solar Motion
    vtemp1 = v_pec
    v_solar = gal2lab(vtemp1,t_lab, lat) # transform to lab co-ords

    #Earth's revolution (first calculate in galactic frame then transform)
    e = eccentricity
    lambda_0 = orb_long_ecliptic
    L = 281.0298 + 36000.77*T_0 + 0.04107*UT
    g = 357.9258 + 35999.05*T_0 + 0.04107*UT
    lambda_sun = L + (1.915 - 0.0048*T_0)*sin(g*pi/180.0)\
         + 0.020*sin(2*g*pi/180.0)
    beta = lat_ecl_gal
    lambda_i = long_ecl_gal
    v_earthrev1 = vv_earthrev*(1-e*sin(pi/180.0*(lambda_sun-lambda_0)))*\
         (cos(beta*pi/180.0)*sin(pi/180.0*(lambda_sun-lambda_i)))
    v_earthrev = gal2lab(v_earthrev1,t_lab, lat) #transform to lab co-ords

    # Earth's rotation (already in lab co-ords)
    v_earthrot = 0.465102*cos(lat*pi/180)*np.array([0.0,-1.0,0.0])

    # Add them all together (delete as needed)
    v_lab = np.array([0.,0.,0.])
    v_lab += v_earthrot
    v_lab += v_earthrev
    v_lab += v_solar
    v_lab += v_galrot

    return v_lab

def JulianDay(month, day, year, hour): # Calculates time in JD for a given date
    year_r = year+4800-floor((14-month)/12.0)
    month_r = month+12*floor((14-month)/12.0)-3
    JulianDay = day + floor((153*month_r+2)/5.0) + 365*year_r\
                + floor(year_r/4.0) - floor(year_r/100.0)\
                + floor(year_r/400.0) - 32045 + (hour-12.0)/24.0
    return JulianDay

def EarthVelocity(day):
    # Second order in eccentricity
    # day measured from Jan1
    lambda_p = 102.93*pi/180.0
    th = w_p*(day-t1)
    v_E = cos(th)*(e1-2*eccentricity*sin(lambda_p)*e2) \
          +sin(th)*(e2+2*eccentricity*sin(lambda_p)*e1) \
          -eccentricity*(cos(2*th)*(cos(lambda_p)*e1-sin(lambda_p)*e2) \
          +sin(2*th)*(sin(lambda_p)*e1+cos(lambda_p)*e2))
    return vv_earthrev*v_E

def EarthVector(day):
    # Earth's orbital radius vectors
    # day measured from Jan1
    # Second order in Earth's eccentricity
    a_earth = AstronomicalUnit/1.0e3
    tp = 3
    lamb_p = 102*pi/180
    g = w_p*(day-tp)
    nu = g + 2.*eccentricity*sin(g)*(5.0/4.0)+eccentricity**2.0*sin(2*g)
    r = a_earth*(1-eccentricity**2.0)/(1+eccentricity*cos(nu))
    r_earth = r*(-sin(lamb_p+nu)*e1 + cos(lamb_p+nu)*e2)
    return r_earth

#==========================Solar direction=====================================#
def EarthSunDistance(JD): # Earth-sun distance at Julian Day (JD)
    D = JD-2451545.0
    g = 357.529 + 0.98560028*D
    g = g*pi/180.0
    r_es = 1.00014 - 0.01671*cos(g) - 0.00014*cos(2*g)
    r_es = r_es*AstronomicalUnit
    return r_es

#------------------------------------------------------------------------------#
def SolarDirection(JD,Loc=Params.GranSasso): # Solar direction in lab coords at Julian Day (JD)

    lat = Loc.Latitude
    lon = Loc.Longitude

    # Compute RA and dec of Sun
    #JD = day+Jan1
    n = JD - 2451545.0
    Omega = 2.1429-0.0010394594*n
    L = 4.8950630 + 0.017202791698*n
    g = 6.2400600 + 0.0172019699*n
    ll = L+0.03341607*sin(g) + 0.00034894*sin(2*g)\
        - 0.0001134 - 0.0000203*sin(Omega)
    ep = 0.4090928 - 6.214e-9*n + 0.0000396*cos(Omega)
    ra = np.arctan2((cos(ep)*sin(ll)),cos(ll)) # Right ascension of Sun
    dec = np.arcsin(sin(ep)*sin(ll)) # Declination of sun

    # Solar vector
    x_sun1 = np.array([0.,0.,0.])
    x_sun1[0] = cos(dec)*cos(ra)
    x_sun1[1] = cos(dec)*sin(ra)
    x_sun1[2] = sin(dec)

    # Lab time conversion
    UT = 24*(JD+0.5-floor(JD+0.5))
    MJD = JD - 2400000.5
    T_0 = (floor(MJD)-55197.5)/36525.0
    t_GAST = (101.0308 + 36000.770*T_0 + 15.04107*UT)/15.0
    t_lab = t_GAST + lon/15.0
    t_lab = 15*t_lab # DEGREES

    # Convert vector from equatorial system into lab system
    x_sun = eqt2lab(x_sun1,t_lab,lat)
    return x_sun

def EarthSunDistanceMod(JD):
    # Solar neutrinos:
    # Flux is scaled by 1/EarthSunDistance^2 but since Flux is already averaged
    # We need to also divide by Integral(1/R^2) over one year
    # Integral_inv_EarthSun_sq is defined in params.f95
    Integral_inv_EarthSun_sq = 4.468864372000642e-23 # integral(1/R^2) over 1 year
    f = (1.0/Integral_inv_EarthSun_sq)*(1.0/EarthSunDistance(JD)**2.0)
    return f

#------------------------------------------------------------------------------#


#==============================================================================#
#---------------------------Coordinate trans.----------------------------------#
def eqt2lab(vp,t_lab,lat): # Equatorial (x_e,y_e,z_e) to Laboratory (N,W,Z)
    t = t_lab*pi/180.0
    latr = lat*pi/180.0
    v = vp*0.0
    v[0] = -cos(t)*sin(latr)*vp[0] - sin(t)*sin(latr)*vp[1] + cos(latr)*vp[2]
    v[1] = sin(t)*vp[0] - cos(t)*vp[1]
    v[2] = cos(t)*cos(latr)*vp[0] + cos(latr)*sin(t)*vp[1] + sin(latr)*vp[2]
    return v

def gal2eqt(vp): # Galactic (x_g,y_g,z_g) to Equatorial (x_e,y_e,z_e)
    v = 0.0*vp
    v[0] = -0.06699*vp[0] + 0.4927*vp[1] - 0.8676*vp[2]
    v[1] = -0.8728*vp[0]  - 0.4503*vp[1] - 0.1884*vp[2]
    v[2] = -0.4835*vp[0]  + 0.7446*vp[1] + 0.4602*vp[2]
    return v

def gal2lab(v,t_lab, lat): # Galactic (x_g,y_g,z_g) to Laboratory (N,W,Z)
    vp = gal2eqt(v)
    return eqt2lab(vp, t_lab, lat)
#==============================================================================#
