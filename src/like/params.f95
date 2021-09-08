module params
  implicit none

!================================param.f95=============================================!
! contains all parameters that need to be defined globally, seems like there is a lot  !
! but it just makes everything a lot nice to look at and less confusing this way.	   !
!																					   !
!																					   !
! Look at this as a dictionary if you want to know what any parameter is for           !
!======================================================================================!

character(len=100),parameter :: homedir='../'
character(len=100),parameter :: mylimitsdir = '../data/WIMPlimits/mylimits'


!------------------------------Nelder Mead stuff---------------------------------------!
integer,parameter :: MAXFUNEVALS = 20000 ! Maximum function evaluations
integer,parameter :: IPRINT = -1 ! Print results from minimisaition
integer,parameter :: NLOOP = 2 ! Number of iterations before looping
integer,parameter :: IQUAD = 0
double precision,parameter :: SIMP = 0.1 ! Nor this
double precision,parameter :: STOPCR0 = 1.0d-10 ! Accuracy of max likelihood

!------------------------------Dark matter stuff---------------------------------------!
double precision :: sigma_p ! WIMP cross section
double precision :: m_chi ! WIMP mass

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

!------------------------------Experimental stuff---------------------------------------!
double precision :: Exposure ! Exposure in ton-years

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


!-------------------------------- Needed for code timing ----------------------------------!
integer :: mytime(3)
real :: clock_start,clock_stop
!------------------------------------------------------------------------------------------!




end module params
