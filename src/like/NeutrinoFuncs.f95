module NeutrinoFuncs
  use params
  use LabFuncs
  use util
  implicit none

contains

!================================NeutrinoFuncs.f95==========================================!
!
! Contents:
! 1. BackgroundRecoilDistribution: Generate Background model RD_bg = (Num of bg events)/Rate
!
! 2. Neutrino event rates and fluxes
! GetNuFluxes = Determines which neutrino backgrounds are needed and loads them in
!
!===========================================================================================!


!===========================================================================================!
! Non-neutrino backgrounds would be added here if they were needed
! Directional detector effects would also be applied here
subroutine BackgroundRecoilDistribution
	double precision :: RD(nTot_bins,n_bg)
	integer :: i,i1,i2,ii,j,s,k

	RD_bg = 0.0d0
	RD = 0.0d0

	! Load all Neutrino RDs
	call NeutrinoRD(nTot_bins,RD)

	! Get total rates and rescale RD_bg by them
	do s = 1,n_bg
		R_bg(s) = sum(RD(:,s))
		RD(:,s) = RD(:,s)/R_bg(s)
	end do
  RD_bg = RD
	! Multiply whole thing by Exposure so RD = Num events/R_bg
	RD_bg = RD_bg*Exposure

end subroutine BackgroundRecoilDistribution
!---------------------------------------------------------------------------------------------!




!===============================Neutrino Recoil distributions=================================!
subroutine NeutrinoRD(n1,RD) ! Generates an RD for all neutrinos
	double precision :: RD(n1,n_bg),fE_r1,fE_r2,E_r1,E_r2,fEmin_r1,fEmin_r2,dpix
	double precision :: dRdE(nE_bins),eff(nE_bins+1),eff_HT(nE_bins+1),R_tot,Flux_t
	integer :: i,j,si,n1,i1,i2,ii,k,n_nu
	RD = 0.0D0
	n_nu = n_bg

	! Load efficiency curve
	eff = Efficiency(E_bin_edges,nE_bins+1)

  do si = 1,n_nu
    ! Integration to get dR/dE -> R
    E_r1 = E_bin_edges(1)
    fE_r1 =  NeutrinoRecoilEnergySpectrum(E_r1,E_nu_all(:,si),Flux_all(:,si))
    do j = 1,nE_bins
      E_r2 = E_bin_edges(j+1)
	    fE_r2 =  NeutrinoRecoilEnergySpectrum(E_r2,E_nu_all(:,si),Flux_all(:,si))
      dRdE(j) = (E_r2-E_r1)*(fE_r1*eff(j) + fE_r2*eff(j+1))/2.0d0
      E_r1 = E_r2
      fE_r1 = fE_r2
    end do
    RD(1:nE_bins,si) = dRdE

    if (nT_bins.gt.0) then
      ! Correct Solar neutrinos for annual modulation
      ii = 1
      do i = 1,nT_bins
      	i1 = ii
      	i2 = i1+nE_bins-1
      	if (solar(si).eq.1) then
      		! Solar neutrinos:
      		! Flux is scaled by 1/EarthSunDistance^2 but since Flux is already averaged
      		! We need to also divide by Integral(1/R^2) over one year
      		! Integral_inv_EarthSun_sq is defined in params.f95
      		RD(i1:i2,si) = (dRdE/Integral_inv_EarthSun_sq)&
      				*(1.0d0/EarthSunDistance(T_bin_centers(i))**2.0d0)
      	else
      		! dsnb and atm neutrinos:
      		! Currently am not incorporating any modulation since the
      		! event rate is nonexistent anyway
      		RD(i1:i2,si) = dRdE
      	end if
      	ii = i2+1
      end do
      RD(:,si) = RD(:,si)/(1.0d0*nT_bins)
    end if


  end do

end subroutine NeutrinoRD
!---------------------------------------------------------------------------------------------!







!========================================Neutrino data========================================!
subroutine GetNuFluxes
	! Reads each neutrino flux data file
	! Each flux file has 1000 rows apart from monochromatic ones
	! the energies are stored in E_nu_all, fluxes in Flux_all
	integer :: i,nvals,ii,j
	integer :: sel(n_nu_tot),s
	double precision :: E_r_max
	nvals = 1000 ! All neutrino backgrounds saved with 1000 entries

	! ORDER OF NEUTRINOS
	! 1. pp
	! 2. pep
	! 3. hep
	! 4. 7Be [343 keV]
	! 5. 7Be [843 keV]
	! 6. 8B
	! 7. 13 N
	! 8. 15 O
	! 9. 17 F
	! 10. DSNB
	! 11. Atm
  ! 12. Geo-U
  ! 13. Geo-Th
  ! 14. Geo-K
  ! 15. Reactor

	! Monochromatic neutrinos (2, 4, 5) have a negative value for E_nu which is
	! Used to tell the rate formula to use monochromatic result

	! Figure out which backgrounds give recoils above E_th
	sel = 0
	n_bg = 0
	do i = 1,n_nu_tot
		E_r_max = MaxNuRecoilEnergy(i) ! Max recoil energy for neutrino
		if (E_r_max.gt.E_th) then
			sel(i) = i ! Set sel if neutrino needs to be loaded
			n_bg = n_bg+1
		end if
		write(*,*) E_r_max,E_th,sel(i),i
	end do


	! Load in all backgrounds with sel != 0
	allocate(E_nu_all(1000,n_bg))
	allocate(Flux_all(1000,n_bg))
	allocate(R_bg_err(n_bg))
	allocate(R_bg(n_bg))
  allocate(solar(n_bg))
  solar = 0
	E_nu_all = 0.0d0
	Flux_all = 0.0d0
	ii = 1
	do i = 1,n_nu_tot
		s = sel(i)
		if (s.gt.0) then
        if (s.le.n_nu_tot_solar) then
          solar(i) = 1
        end if
		    if (s.eq.1) then
		       open(unit=10,file='../../data/neutrinos/pp-1000.txt')
		    else if (s.eq.2) then
		       open(unit=10,file='../../data/neutrinos/pep-1000.txt')
		    else if (s.eq.3) then
		       open(unit=10,file='../../data/neutrinos/hep-1000.txt')
		    else if (s.eq.4) then
		       open(unit=10,file='../../data/neutrinos/7Be1-1000.txt')
		    else if (s.eq.5) then
		       open(unit=10,file='../../data/neutrinos/7Be2-1000.txt')
		    else if (s.eq.6) then
		       open(unit=10,file='../../data/neutrinos/8B-1000.txt')
		    else if (s.eq.7) then
		       open(unit=10,file='../../data/neutrinos/13N-1000.txt')
		    else if (s.eq.8) then
		       open(unit=10,file='../../data/neutrinos/15O-1000.txt')
		    else if (s.eq.9) then
		       open(unit=10,file='../../data/neutrinos/17F-1000.txt')
		    else if (s.eq.10) then
		       open(unit=10,file='../../data/neutrinos/DSNB-1000.txt')
		    else if (s.eq.11) then
		       open(unit=10,file='../../data/neutrinos/Atm-1000.txt')
        else if (s.eq.12) then
          open(unit=10,file='../../data/neutrinos/GeoU-1000.txt')
        else if (s.eq.13) then
          open(unit=10,file='../../data/neutrinos/GeoTh-1000.txt')
        else if (s.eq.14) then
          open(unit=10,file='../../data/neutrinos/GeoK-1000.txt')
        else if (s.eq.15) then
          open(unit=10,file='../../data/neutrinos/Reactor-1000.txt')
		    end if
		    do j = 1,nvals
		      read(10,*) E_nu_all(j,ii),Flux_all(j,ii)
		    end do
		    close(10)
		    Flux_all(:,ii) = NuFlux(i)*Flux_all(:,ii) ! Select rate normalisations
			R_bg_err(ii) = NuUnc(i) ! Select rate normalisation uncertainties
			ii = ii+1
		end if
	end do

end subroutine GetNuFluxes

 !---------------------------------------------------------------------------------------------!
function MaxNuRecoilEnergy(s) result(E_r_max) ! Max recoil energy for neutrino number s
	integer :: s
	double precision :: E_r_max,m_N_keV
	m_N_keV = 0.93141941*(nucleus(1)+nucleus(2))*1.0e6
	E_r_max = 2*m_N_keV*(1000.0*NuMaxEnergy(s))**2.0&
	   /(m_N_keV+1000*NuMaxEnergy(s))**2.0
end function

!===================================nu spectra=================================!
function NeutrinoRecoilEnergySpectrum(E_r,E_nu,Flux) result(dRdE)
	integer :: N,Z,nn,i
	double precision :: Q_W,E_r,dRdE,m_N_keV,m_N_GeV
	double precision, dimension(1000) :: E_nu,Flux,diff_sigma,dR
	nn = 1000
	N = nucleus(1) ! number of neutrons
	Z = nucleus(2) ! number of protons
	Q_W = N-(1-4.0d0*sintheta_Wsq)*Z ! weak nuclear hypercharge
	m_N_GeV = 0.93141941d0*(N+Z) ! nucleus mass in GeV
	m_N_keV = m_N_GeV*1.0d6 ! nucleus mass in keV

	! differential cross section
	diff_sigma = (G_F_GeV**2.0d0/(4.0d0*pi))*(Q_W**2.0d0)*&
	     m_N_GeV*(1.0d0-(m_N_keV*E_r)/(2.0d0*(E_nu*1000.0d0)**2.0d0))*&
	     (0.197d-13)**2.0d0*(1.0d-6)*1000.0d0/(1.0d0*N+1.0d0*Z)*(N_A)*&
	     FormFactorHelm(E_r,N+Z)**2.0d0

	! diff_sigma goes negative for kinematically forbidden energies so just set them 0 and sum
	where (diff_sigma.lt.0.0d0)
	   diff_sigma=0.0d0
	end where

	! Event rate is integral(cross section x flux)
	if (Flux(2).gt.0.0) then
		dR = diff_sigma*Flux
	 	dRdE = sum(0.5d0*(E_nu(2:nn)-E_nu(1:nn-1))*(dR(2:nn)+dR(1:nn-1)))
	else
		dRdE = diff_sigma(1)*Flux(1)*E_nu(1) ! for monochromatic nu's
	end if

	! Convert into /ton/year/keV
	dRdE = dRdE*(365.0d0*3600.0d0*24d0)*(1000.0)
end function NeutrinoRecoilEnergySpectrum

end module NeutrinoFuncs
