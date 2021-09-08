module like
  use params
  use util
  implicit none

contains

!====================================like.f95==================================!
! Module for doing all of the Likelihood analysis
!
! Contents:
!
! 1. Compute limits for input experiment
! DisocveryLimit calculates a single DL for one Exposure
! DiscoveryLimit_vs_MassCrossSection calculates the required exposure for a given grid of masses and cross sections

! 2. Likelihoods
! lnPF = sum over log(Poisson pdf)
! lnGF = sum over log(Gaussian pdf)
! llhood1 = signal+background likelihood
! llhood0 = background likelihood
! llhood1_bigN = signal+background likelihood for when N is very large
! llhood0_bigN = background likelihood for when N is very large
!==============================================================================!





!=====================================DISCOVERY LIMIT==========================!
subroutine DiscoveryLimit(m_vals,nm,sigma_min,sigma_max,ns,verbose,DL)
	! * Calculates discovery limits (DL) from m_chi = m_min to m_max
	! * Parameters of likelihood are X = (sigma_p,R_bg)
	! * Assumes signal model RD_wimp and background model RD_bg
	! * RD_wimp = (Number of WIMP events in each bin)/sigma_p
	! * RD_bg = (Number of BG events each bin)/R_bg
	! * Parameters enter likelihood function to rescale RD_wimp and RD_bg
	! * Each mass scans cross sections from sigma_min to sigma_min
	! * Done this way because the minimisation is fastest for low event numbers
	! * Uses Asimov asymptotic result to get median limit for 3sigma detection
	integer :: nf,nm,i,ns,im,j,ifault0,si,ii,verbose
	double precision :: sigma_p_vals(ns),DL(nm),m_vals(nm),sigma_min,sigma_max
	double precision :: x_in0(n_bg),x_in1(n_bg+1),step0(n_bg),N_exp(nTot_bins),N_exp_bg(nTot_bins)
	double precision :: D01,L1,L0,var(2),N_tot_bg,D_prev,s_prev

	! Mass and cross section discretisation
	sigma_p_vals = logspace(sigma_min,sigma_max,ns)
	DL = 0.0d0

	! GENERATE BACKGROUND DATA
	N_exp_bg = 0.0d0
	do si = 1,n_bg
		N_exp_bg = N_exp_bg + RD_bg(:,si)
	end do
	N_tot_bg = sum(N_exp_bg)
	! MASS SCAN:
	do im = 1,nm
		m_chi = m_vals(im)
		RD_wimp = RD_sig(im,:)/1.0d-45	! Call WIMP recoil distribution for each new mass
		! CROSS SECTION SCAN
    X_in0 = R_bg

		do j = 1,ns
			sigma_p = sigma_p_vals(j)
			if (sum(RD_wimp*sigma_p).gt.0.1d0) then	! Generally need >0.5 events to see DM
        N_exp = N_exp_bg + RD_wimp*sigma_p
        N_obs = N_exp  ! Observed=Expected for Asimov data
        X_in1= (/log10(sigma_p),R_bg/)
        step0 = R_bg_err*X_in0

          if (sum(RD_wimp*sigma_p).gt.1.0d3) then
            call llhood1_bigN(X_in1,L1) ! Asimov data maximises likelihood at correct value
            call llhood0_bigN(X_in0,L0)
            call MINIM(X_in0,step0,n_bg,L0,MAXFUNEVALS,IPRINT,STOPCR0,NLOOP,IQUAD,SIMP,VAR,llhood0_bigN,IFAULT0)
          else
            call llhood1(X_in1,L1) ! Asimov data maximises likelihood at correct value
            call llhood0(X_in0,L0)
            call MINIM(X_in0,step0,n_bg,L0,MAXFUNEVALS,IPRINT,STOPCR0,NLOOP,IQUAD,SIMP,VAR,llhood0,IFAULT0)
          end if

				! Test statistic
				D01 = -2.0*(L1-L0)
        !write(*,*) j,sigma_p,D01
				if (D01.ge.9.0d0) then ! Median 3sigma detection -> D = 9
					! Do interpolation to find discovery limit cross section
					DL(im) = 10.0d0**(interp1D((/D_prev,D01/),(/log10(s_prev),log10(sigma_p)/),2,9.0d0))
					exit
				end if
				s_prev = sigma_p ! Reset for interpolation
				D_prev = D01
			end if
		end do
    if (j.eq.1) then
      DL(im) = sigma_min
    end if
    if (verbose.eq.1) then
		    write(*,*) 'm = ',m_chi,'DL = ',DL(im),'Sig:',sum(RD_wimp*sigma_p),'BG:',N_tot_bg
    end if
		!stop
	end do
end subroutine


subroutine DiscoveryLimit_vs_MassCrossSection(m_vals,nm,sigma_p_vals,ns,ex_vals,n_ex,DL)
  integer :: i,j,k,k1,nm,nf,ns,n_ex,si,ifault0,ii
  double precision :: x_in0(n_bg),x_in1(n_bg+1),step0(n_bg),N_exp(nTot_bins),N_exp_bg(nTot_bins)
	double precision :: m_min,m_max,ex_min,ex_max,m,sigma_min,sigma_max
  double precision :: DL(ns,nm),m_vals(nm),ex_vals(n_ex),sigma_p_vals(ns)
  double precision :: N_tot_bg,D01,L0,L1,D_prev,ex_prev,var(2)

  ! mvals = array of masses to scan over
  ! sigma_p_vals = array of cross sections to scan over
  ! ex_vals = array of exposures to try
  ! DL = discovery limit for the 2D grid of masses and cross sections (output)

  ! This routine scans over masses and cross sections to find the exposure needed for the median experiment to observe that cross section
  ! Some parts of the mass, cross section grid will require insane exposures, so those parts are sped up by minimising a chi2 rather than a poisson Likelihood
  ! The initial guess for the minimisation is X_in0 which is set to R_bg initially (background amplitudes), but then is held over through each iteration for speed

  DL = 0.0d0
  ! Mass scan
  do i = 1,nm
    k1 = 1
    Exposure = 1.0
    m_chi = m_vals(i)
		RD_wimp = RD_sig(i,:)/1.0d-45	! Call WIMP recoil distribution for each new mass (RD_sig is always defined for sigma = 10^-45 cm^2)
    X_in0 = R_bg ! Initial guess set of parameters
    if (sum(RD_wimp).gt.0.0) then
      ! cross section scan
      do j=1,ns
        sigma_p = sigma_p_vals(ns+1-j) ! Go backwards through the cross sections since the larger ones require smaller exposures
        ex_prev = ex_min
        do k = k1,n_ex
          Exposure = ex_vals(k)
          RD_wimp = RD_wimp*Exposure ! Always scale signal by exposure
          RD_bg = RD_bg*Exposure ! Always scale background by exposure

          ! Background events:
          N_exp_bg = 0.0d0
        	do si = 1,n_bg
        		N_exp_bg = N_exp_bg + RD_bg(:,si)
        	end do
        	N_tot_bg = sum(N_exp_bg)


          if (sum(RD_wimp*sigma_p).gt.0.5d0) then	! Generally need >0.5 events to see DM
            N_exp = N_exp_bg + RD_wimp*sigma_p
            N_obs = N_exp  ! Observed=Expected for Asimov data
            X_in1= (/log10(sigma_p),R_bg/) ! signal+background model parameters
            step0 = R_bg_err*X_in0 ! initial step size

              if (sum(RD_wimp*sigma_p).gt.1.0d3) then
                call llhood1_bigN(X_in1,L1) ! Asimov data maximises likelihood at correct value
                call llhood0_bigN(X_in0,L0)
                call MINIM(X_in0,step0,n_bg,L0,MAXFUNEVALS,IPRINT,STOPCR0,NLOOP,IQUAD,SIMP,VAR,llhood0_bigN,IFAULT0)
              else
                call llhood1(X_in1,L1) ! Asimov data maximises likelihood at correct value
                call llhood0(X_in0,L0)
                call MINIM(X_in0,step0,n_bg,L0,MAXFUNEVALS,IPRINT,STOPCR0,NLOOP,IQUAD,SIMP,VAR,llhood0,IFAULT0)
              end if

            D01 = -2.0*(L1-L0) ! Test statistic
            if (D01.ge.9.0d0) then ! Median 3sigma detection -> D = 9
              ! Do interpolation to find discovery limit cross section
              if (k.eq.1) then
                DL(j,i) = ex_min
                k1 = 1
              else
                DL(j,i) = 10.0d0**(interp1D((/D_prev,D01/),(/log10(ex_prev),log10(Exposure)/),2,9.0d0)) ! interpolate to find where DL crosses 9
                k1 = k-1
              end if
              write(*,*) i,'of',nm,'|| m = ',m_chi,'|| sigma = ',sigma_p,'|| Ex = ',DL(j,i),'|| N_obs = ',sum(N_obs),L0,L1
              RD_bg = RD_bg/Exposure
              RD_wimp = RD_wimp/Exposure
              exit
            end if
            ex_prev = Exposure ! Reset for interpolation
            D_prev = D01
            !write(*,*) Exposure,D01
          end if
          RD_bg = RD_bg/Exposure
          RD_wimp = RD_wimp/Exposure
        end do

        if ((k-1).eq.n_ex) then
          DL(j:ns,i) = ex_max
          exit
        end if
      end do
    end if
  end do
end subroutine




!===================================USEFUL SUMS================================!
function lnPF(nbins,Nob,Nex)! SUM OF LOG(POISSON PDF)
	! Uses log-gamma function to generalise the factorial
	! to non-integer number of observed events (i.e. Asimov data)
	double precision :: lnPF,Nex(nbins),Nob(nbins)
	integer :: ii,nbins
	lnPF = 0.0d0
	do ii = 1,nbins
		lnPF = lnPF + Nob(ii)*log(Nex(ii)) - Nex(ii) - lgamma(Nob(ii)+1.0)
	end do
end function

function lnGF(Rob,Rex,Rer) ! SUM OF LOG(GAUSSIAN PDF)
	! for x=Rob, mu = Rex, sig = Rer
	double precision :: Rob(:),Rex(:),Rer(:),lnGF
	lnGF = sum(-1.0d0*log(Rer)-0.5d0*log(2.0d0*pi)&
	       	 -(Rob-Rex)**2.0d0/(2.0d0*Rer**2.0d0))
end function




!===================================LIKELIHOODS================================!
!---------------------------------SIGNAL+BACKGROUND----------------------------!
 subroutine llhood1(X,  LL)
    ! input: Parameters X = (log10(sigma_p),background rates(1:n_bg))
	! output: -1*LogLikelihood = LL
	double precision :: X(n_bg+1),LL,N_exp0(nTot_bins),N_exp1(nTot_bins)
	integer :: i
	! Background events
	N_exp0 = 0.0d0 ! Expected number of events
	do i = 1,n_bg
		N_exp0 = N_exp0 + X(i+1)*RD_bg(:,i)/R_bg(i) ! Sum over backgrounds
	end do

	! Signal events
	N_exp1 = N_exp0 + RD_wimp*(10.0d0**X(1)) ! Add signal events sig


	! LL = log(Poiss. for N_obs events) + log(Gauss. for R_bg normalisations)
	LL = -1.0*(lnPF(nTot_bins,N_obs,N_exp1)+lnGF(X(2:n_bg+1),R_bg,R_bg_err*R_bg))

 end subroutine llhood1

!---------------------------------SIGNAL+BACKGROUND----------------------------!
subroutine llhood0(X,  LL)
    ! input: Parameters X = background rates(1:n_bg)
	! output: -1*LogLikelihood = LL
    double precision :: X(n_bg),LL,N_exp0(nTot_bins)
	integer :: i
	N_exp0 = 0.0d0 ! Expected number of events
	do i = 1,n_bg
		N_exp0 = N_exp0 + X(i)*RD_bg(:,i)/R_bg(i) ! Sum over backgrouds
	end do

	! LL = log(Poiss. for N_obs events) + log(Gauss. for R_bg normalisations)
	LL = -1.0*(lnPF(nTot_bins,N_obs,N_exp0)+lnGF(X,R_bg,R_bg_err*R_bg))
  end subroutine llhood0




!---------------------------------Large N limit----------------------------!
subroutine llhood1_bigN(X,  LL)
   ! input: Parameters X = (log10(sigma_p),background rates(1:n_bg))
 ! output: -1*LogLikelihood = LL
 double precision :: X(n_bg+1),LL,N_exp0(nTot_bins),N_exp1(nTot_bins)
 integer :: i
 ! Background events
 N_exp0 = 0.0d0 ! Expected number of events
 do i = 1,n_bg
   N_exp0 = N_exp0 + X(i+1)*RD_bg(:,i)/R_bg(i) ! Sum over backgrounds
 end do

 ! Signal events
 N_exp1 = N_exp0 + RD_wimp*(10.0d0**X(1)) ! Add signal events sig

 LL = 0.5*sum((N_obs-N_exp1)**2.0/N_exp1)
 LL = LL-lnGF(X(2:n_bg+1),R_bg,R_bg_err*R_bg)

end subroutine llhood1_bigN

subroutine llhood0_bigN(X,  LL)
   ! input: Parameters X = background rates(1:n_bg)
 ! output: -1*LogLikelihood = LL
   double precision :: X(n_bg),LL,N_exp0(nTot_bins)
 integer :: i
 N_exp0 = 0.0d0 ! Expected number of events
 do i = 1,n_bg
   N_exp0 = N_exp0 + X(i)*RD_bg(:,i)/R_bg(i) ! Sum over backgrouds
 end do

 LL = -1.0*(-0.5*sum((N_obs-N_exp0)**2.0/N_exp0)+lnGF(X,R_bg,R_bg_err*R_bg))
end subroutine llhood0_bigN













! !==================================Generate limits=============================!
! subroutine Limit_vs_Mass(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DL)
! 	! Limits from m_chi = m_min to m_max, with nm values
! 	! Cross section scan from sigma_min to sigma_max with ns values
! 	! DL = Discovery Limit
! 	double precision :: m_min,m_max,sigma_min,sigma_max,m_vals(nm),DL(nm)
! 	integer :: i,nm,nf,ns
!   write(*,*) 'Nucleus = ',nucleus,'Exposure = ',Exposure,'ton years'
!   write(*,*) '[E_th,E_max] = ',E_th,E_max,'keV'
!   write(*,*) '----------------------------------------------------'
!   call GetNuFluxes ! Load Neutrinos
!   call PreAllocate ! Allocate data size (readout dependent)
!   call BackgroundRecoilDistribution ! Load Background	model
!   call SHM ! Load halo model
!   call DiscoveryLimit(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DL)
!   call UnAllocate ! Reset
! end subroutine
!
! subroutine Limit_vs_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL,Nsig,Nbg)
! 	! Limits for single m_chi vs exposure from ex_min to ex_max
! 	! Cross section scan from sigma_min to sigma_max with ns values
! 	! DL = Discovery Limit
! 	double precision :: sigma_min1,ex_min,ex_max,m,sigma_min,sigma_max,m_vals(1)
!   double precision :: DL(n_ex),ex_vals(n_ex),Nsig(n_ex),Nbg(n_ex)
! 	integer :: j,i,nf,ns,n_ex,si,ns_reduced
!   write(*,*) 'Nucleus = ',nucleus,'Exposure = ',Exposure,'ton years'
!   write(*,*) '[E_th,E_max] = ',E_th,E_max,'keV'
!   write(*,*) '----------------------------------------------------'
!   Exposure = 1.0
!   call GetNuFluxes ! Load Neutrinos
!   call PreAllocate ! Allocate data size (readout dependent)
!   call BackgroundRecoilDistribution ! Load Background	model
!   call SHM ! Load halo model
!   ex_vals = logspace(ex_min,ex_max,n_ex)
!   sigma_min1 = sigma_min
!   do i=1,n_ex
!     j = n_ex+1-i
!     Exposure = ex_vals(j)
!     RD_bg = RD_bg*Exposure
!     ns_reduced = floor(ns*(log10(sigma_max)-log10(sigma_min1))/(log10(sigma_max)-log10(sigma_min1)))
!     call DiscoveryLimit(m,m,1,sigma_min1,sigma_max,ns_reduced,	m_vals,DL(j))
!     sigma_min1 = DL(j)
!     ! numbers of events
!     Nsig(j) = sum(RD_wimp*DL(j))
!     Nbg(j) = 0.0d0
!     do si = 1,n_bg
!       Nbg(j) = Nbg(j) + sum(R_bg(si)*RD_bg(:,si))
!     end do
!
!     RD_bg = RD_bg/Exposure
!   end do
!   call UnAllocate ! Reset
! end subroutine
!
!
! subroutine Limit_vs_MassCrossSection(m_min,m_max,nm,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,DL)
!   integer :: i,j,k,k1,nm,nf,ns,n_ex,si,ifault0,ii
! 	double precision :: m_min,m_max,ex_min,ex_max,m,sigma_min,sigma_max
!   double precision :: DL(nm,ns),m_vals(nm),ex_vals(n_ex),sigma_p_vals(ns)
!   double precision :: N_tot_bg,D01,L0,L1,D_prev,ex_prev,var(2)
!   double precision,dimension(:),allocatable :: x_in0,x_in1,step0, N_exp,N_exp_bg
!   write(*,*) 'Nucleus = ',nucleus,'Exposure = ',ex_min,ex_max
!   write(*,*) '[E_th,E_max] = ',E_th,E_max,'keV'
!   write(*,*) '----------------------------------------------------'
!
!   Exposure = 1.0
!   call GetNuFluxes ! Load Neutrinos
!   call PreAllocate ! Allocate data size (readout dependent)
!   allocate(N_exp_bg(nTot_bins))
!   allocate(N_exp(nTot_bins))
!   allocate(x_in0(n_bg))
!   allocate(step0(n_bg))
!   allocate(x_in1(n_bg+1))
!   call BackgroundRecoilDistribution ! Load Background	model
!   call SHM ! Load halo model
!   ex_vals = logspace(ex_min,ex_max,n_ex)
!   m_vals = logspace(m_min,m_max,nm)
!   sigma_p_vals = logspace(sigma_min,sigma_max,ns)
!   DL = 0.0d0
!
!   do i = 1,nm
!     k1 = 1
!     m_chi = m_vals(i)
!     Exposure = 1.0
!     call WIMPRecoilDistribution	! Call WIMP recoil distribution for each new mass
!
!     if (sum(RD_wimp).gt.0.0) then
!       do j=1,ns
!         sigma_p = sigma_p_vals(ns+1-j)
!
!         do k = k1,n_ex
!           Exposure = ex_vals(k)
!           RD_wimp = RD_wimp*Exposure
!           RD_bg = RD_bg*Exposure
!           N_exp_bg = 0.0d0
!           do si = 1,n_bg
!             N_exp_bg = N_exp_bg + R_bg(si)*RD_bg(:,si)
!           end do
!           N_tot_bg = sum(N_exp_bg)
!
!           if (sum(RD_wimp*sigma_p).gt.0.5d0) then	! Generally need >0.5 events to see DM
!             N_exp = N_exp_bg + RD_wimp*sigma_p
!             N_obs = N_exp  ! Observed=Expected for Asimov data
!             X_in1= (/log10(sigma_p),R_bg/)
!             X_in0 = R_bg
!             step0 = R_bg_err*R_bg
!
!             if (sum(RD_wimp*sigma_p).gt.1.0d3) then
!              call llhood1_bigN(X_in1,L1) ! Asimov data maximises likelihood at correct value
!              call llhood0_bigN(X_in0,L0)
!              call MINIM(X_in0,step0,n_bg,L0,MAXFUNEVALS,IPRINT,STOPCR0,NLOOP,IQUAD,SIMP,VAR,llhood0_bigN,IFAULT0)
!             else
!               call llhood1(X_in1,L1) ! Asimov data maximises likelihood at correct value
!               call llhood0(X_in0,L0)
!               call MINIM(X_in0,step0,n_bg,L0,MAXFUNEVALS,IPRINT,STOPCR0,NLOOP,IQUAD,SIMP,VAR,llhood0,IFAULT0)
!             end if
!
!             D01 = -2.0*(L1-L0)
!             if (D01.ge.9.0d0) then ! Median 3sigma detection -> D = 9
!               ! Do interpolation to find discovery limit cross section
!               if (k.eq.1) then
!                 DL(i,j) = ex_min
!                 k1 = 1
!               else
!                 DL(i,j) = 10.0d0**(interp1D((/D_prev,D01/),(/log10(ex_prev),log10(Exposure)/),2,9.0d0))
!                 k1 = k-1
!               end if
!               write(*,*) 'm = ',m_chi,'|| sigma = ',sigma_p,'|| DL = ',DL(i,j),'|| N_obs = ',sum(N_obs),L0,L1
!               RD_bg = RD_bg/Exposure
!               RD_wimp = RD_wimp/Exposure
!               exit
!             end if
!             ex_prev = Exposure ! Reset for interpolation
!             D_prev = D01
!             !write(*,*) Exposure,D01
!           end if
!           RD_bg = RD_bg/Exposure
!           RD_wimp = RD_wimp/Exposure
!         end do
!
!         if ((k-1).eq.n_ex) then
!           DL(i,j:ns) = ex_max
!           exit
!         end if
!       end do
!     end if
!   end do
!   call UnAllocate ! Reset
! end subroutine


end module like
