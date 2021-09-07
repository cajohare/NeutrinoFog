module params
  implicit none

!================================param.f95=============================================!
! contains all parameters that need to be defined globally, seems like there is a lot  !
! but it just makes everything a lot nice to look at and less confusing this way.	   !
!																					   !
!																					   !
! Look at this as a dictionary if you want to know what any parameter is for           !
!======================================================================================!

character(len=100),parameter :: homedir='/Users/ciaranohare/Work/NeutrinoFog'
character(len=100),parameter :: mylimitsdir = '/Users/ciaranohare/Work/NeutrinoFog/data/WIMPlimits/mylimits'



!------------------------------Nelder Mead stuff---------------------------------------!
integer,parameter :: MAXFUNEVALS = 20000 ! Maximum function evaluations
integer,parameter :: IPRINT = -1 ! Print results from minimisaition
integer,parameter :: NLOOP = 2 ! Number of iterations before looping
integer,parameter :: IQUAD = 0 ! Can't remeber what this does
double precision,parameter :: SIMP = 0.1 ! Nor this
double precision,parameter :: STOPCR0 = 1.0d-10 ! Accuracy of max likelihood



!------------------------------Dark matter stuff---------------------------------------!
double precision :: sigma_p ! WIMP cross section
double precision :: m_chi ! WIMP mass
double precision :: rho_0 ! Local DM density
double precision :: v_esc ! Escape speed
double precision :: v_LSR ! Local standard of rest
double precision :: sig_v ! Halo dispersion
double precision,dimension(3),parameter :: v_pec = (/11.1d0,12.24d0,7.25d0/) ! Solar vel.
double precision,dimension(:,:),allocatable :: v_lab_all ! Time dependent lab velocity
!---------------------------------------------------------------------------------------!






!--------------------------------------Backgrounds---------------------------------------!
integer :: n_bg ! Number of backgrounds
integer :: n_nu_tot
integer,parameter :: n_nu_tot_solar = 9
integer,dimension(:),allocatable :: solar
double precision,dimension(:),allocatable :: R_bg ! Background rate normalisations
double precision,dimension(:),allocatable :: R_bg_err ! Background rate uncertainties
double precision,dimension(:),allocatable :: NuUnc ! All neutrino uncertainties
double precision,dimension(:),allocatable  :: NuFlux ! All neutrino fluxes
double precision,dimension(:,:),allocatable :: Flux_all ! All neutrino fluxes
double precision,dimension(:,:),allocatable :: E_nu_all ! All neutrino energies
! Neutrino max energies (useful for figuring out which ones need to be considere)
double precision,dimension(15),parameter :: NuMaxEnergy = (/0.42341d0,1.44d0,18.765d0,0.3843d0,0.8613d0,16.34d0,1.193d0,1.7285d0,1.7365d0,91.201d0,981.75d0,4.54d0,2.33d0,1.3572d0,1.1418d1/)
! Neutrino fluxes and uncertainties from Vinyoles et al Barcelona GS98 SSM
double precision,dimension(9),parameter :: NuFlux_B17GS98 = (/5.98d10,1.44d8,7.98d3,4.93d8,4.50d9,5.45d6,2.78d8,2.05d8,5.29d6/)
double precision,dimension(9),parameter :: NuUnc_B17GS98 = (/0.006d0, 0.01d0, 0.3d0,0.06d0, 0.06d0, 0.12d0, 0.15d0 ,0.17d0 ,0.2d0/)
! Neutrino fluxes and uncertainties from Vinyoles et al Barcelona AGSS09met SSM
double precision,dimension(9),parameter :: NuFlux_B17AGSS = (/6.03d10,1.46d8,8.25d3,4.5d8,4.50d9,4.5d6,2.04d8,1.44d8,3.26d6/)
double precision,dimension(9),parameter :: NuUnc_B17AGSS = (/0.005d0, 0.009d0, 0.3d0,0.06d0, 0.06d0, 0.12d0, 0.14d0 ,0.16d0 ,0.18d0/)
double precision :: B8Flux_Bergstrom = 5.16d6
double precision :: B8Unc_Bergstrom = 0.02d0

double precision :: NuFlux_Atm = 10.54d0
double precision :: NuUnc_Atm = 0.25d0

double precision :: NuFlux_DSNB = 85.7d0
double precision :: NuUnc_DSNB = 0.5d0

double precision,dimension(3),parameter :: NuFlux_Geo = (/3808776.91874,3352686.94783,21639789.2056/) ! U,Th,K
double precision,dimension(3),parameter :: NuUnc_Geo = (/0.2,0.257,0.168/) ! U,Th,K
double precision :: NuFlux_Reactor = 208537.673299
double precision :: NuUnc_Reactor = 0.08
!---------------------------------------------------------------------------------------!






!------------------------------Experimental stuff---------------------------------------!
double precision :: Exposure ! Exposure in ton-years
double precision :: E_th ! Min analysis energy
double precision :: E_max ! Max analysis energy
double precision,dimension(1000) :: E_vals ! Energy values for detector performance
double precision,dimension(1000,2) :: energyres_data ! energy resolution (He,F)
double precision,dimension(1000,2) :: efficiency_data ! efficiency (He,F)
integer :: nucleus(2) ! target nucleus = (Number of neutrons, number of protons)
logical :: energyres_on ! Switch for turning efficiency on
logical :: efficiency_on ! Switch for turning efficiency on
!---------------------------------------------------------------------------------------!





!------------------------------Binning parameters---------------------------------------!
integer :: nE_bins ! Number of bins in Energy (from E_th to E_max)
integer :: nTot_bins ! Number of bins used for likelihood
integer :: nT_bins ! Number of bins in Time (from Jan 1 to Dec 31)
double precision,dimension(:),allocatable :: E_bin_centers ! Center of energy bins
double precision,dimension(:),allocatable :: E_bin_edges ! Edges of energy bins
double precision,dimension(:),allocatable :: T_bin_centers ! Center of time bins
!------------------------------------------------------------------------------------------!





!--------------------------Signal and background models----------------------------------!
double precision,dimension(:),allocatable :: RD_wimp ! Indiv WIMP Signal model
double precision,dimension(:,:),allocatable :: RD_bg,RD_sig ! Background/Signal model
double precision,dimension(:),allocatable :: N_obs ! N of observed events (for likelihood)
!------------------------------------------------------------------------------------------!




!--------------------------------------General constants----------------------------------!
double precision,parameter :: pi = 3.141592653589793
double precision,parameter :: G_F_GeV=1.16637d-5 ! GeV**-2 ! Fermi constan in GeV
double precision,parameter :: G_F_keV=G_F_GeV/(1.0d6)**2.0d0 ! Fermi constant in keV
double precision,parameter :: N_A=6.02214d23 ! Avocado's constant
double precision,parameter :: sinTheta_Wsq=0.2387d0 ! sin(Theta_W) weak mixing angle
double precision,parameter :: AstronomicalUnit = 1.49597892d11 ! Astronomical Unit
double precision,parameter :: EarthRadius = 6371.01d0*1000.0d0 ! Earth Radius
double precision,parameter :: Integral_inv_EarthSun_sq = 4.468864372000642d-23 ! 1/r^2 int
double precision,parameter :: Jan1 = 2458849.5! Julian Day of January 1 2020
double precision,parameter :: vv_earthrev = 29.79d0 ! Earth revolution speed
double precision,parameter :: eccentricity = 0.016722d0 ! eccentricity of earth orbit
double precision,parameter :: eccentricity_deg = 0.9574d0 ! eccentricity in degrees
double precision,parameter :: orb_long_ecliptic = 13.0d0+1.0d0 !
double precision,dimension(3),parameter :: lat_ecl_gal = (/-5.5303d0,59.575d0,29.812d0/)
double precision,dimension(3),parameter :: long_ecl_gal =(/266.141d0,-13.3485d0,179.3212d0/)
!------------------------------------------------------------------------------------------!




!--------------------------------------Various Locations-----------------------------------!
double precision, parameter :: lat = 54.5591d0  ! latitude of expt in degrees
double precision, parameter :: long = 0.8310d0 ! longitude of expt in degrees
double precision,dimension(2),parameter :: Boulby = (/54.5591d0,0.8310d0/)
double precision,dimension(2),parameter :: GranSasso = (/2.4691, 13.5654/)
double precision,dimension(2),parameter :: Kamioka = (/36.2381, 137.1863/)
double precision,dimension(2),parameter :: SNOlab = (/46.4719, -81.1868/)
double precision,dimension(2),parameter :: Stawell = (/-37.0576, 142.7754/)
double precision,dimension(2),parameter :: Oahu = (/21.4389, -158.0001/)
!------------------------------------------------------------------------------------------!




!--------------------------------------Various Nuclei--------------------------------------!
! nucleus = (Number of neutrons, number of protons)
integer,dimension(2),parameter :: Helium = (/2,2/)
integer,dimension(2),parameter :: Fluorine = (/10,9/)
integer,dimension(2),parameter :: Xenon131 = (/77,54/)
integer,dimension(2),parameter :: Xenon129 = (/74,54/)
integer,dimension(2),parameter :: Germanium74 = (/42,32/)
integer,dimension(2),parameter :: Argon40 = (/22,18/)
integer,dimension(2),parameter :: Silicon = (/14,14/)
integer,dimension(2),parameter :: Calcium = (/20,20/)
integer,dimension(2),parameter :: Iodine = (/73,53/)
integer,dimension(2),parameter :: Carbon12 = (/6,6/)
!------------------------------------------------------------------------------------------!




!-------------------------------- Needed for code timing ----------------------------------!
integer :: mytime(3)
real :: clock_start,clock_stop
!------------------------------------------------------------------------------------------!




end module params
