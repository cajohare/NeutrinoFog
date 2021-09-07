program runDL
	use params
	use util
	use like
	implicit none

	integer :: i,ns,nm,n_ex,cnt,verbose
	double precision :: m_min,m_max,sigma_min,sigma_max,ex_min,ex_max,junk_db
	double precision,dimension(:),allocatable :: m_vals,ex_vals,DL_ex
	double precision,dimension(:),allocatable :: DL,sigma_p_vals
	character(len=10) :: junk1,junk2,junk3
	character(len=200) :: filename1,inp
  CHARACTER(len=1) :: arg


	cnt = COMMAND_ARGUMENT_COUNT()
	if(cnt.eq.0) then
  	write(*,*)'ERROR, AT LEAST ONE COMMAND-LINE ARGUMENTS REQUIRED, STOPPING'
  	stop
	end if

  CALL get_command_argument(1, inp)
	if (cnt.eq.2) then
		CALL get_command_argument(2, arg)
		read(arg,*) verbose
	else
		verbose = 0
	end if

	if ((verbose.ne.0).and.(verbose.ne.1)) then
		write(*,*)'ERROR, VERBOSE MUST BE 0 OR 1, STOPPING'
		stop
	end if

	! Seed rand and random
	call itime(mytime)
	call cpu_time(clock_start)

	write(*,*) 'Starting'

	! Resolution of floor
	sigma_min = 1.0d-51
	sigma_max = 1.0d-41
	ns = 200

	! Read in signal distribution
	open(unit=111,file=trim(homedir)//'/data/recoils/RD_sig_'//trim(inp)//'.txt')
	read(111,*) junk1,nm,nTot_bins
	allocate(RD_sig(nm,nTot_bins))
	allocate(RD_wimp(nTot_bins))
	allocate(m_vals(nm))
	allocate(DL(nm))
	do i = 1,nm
		read(111,*) m_vals(i),RD_sig(i,:)
	end do
	close(111)

	! Read in background distribution
	open(unit=222,file=trim(homedir)//'/data/recoils/RD_bg_'//trim(inp)//'.txt')
	read(222,*) junk1,n_bg,nTot_bins
	allocate(RD_bg(nTot_bins,n_bg))
	allocate(R_bg(n_bg))
	allocate(R_bg_err(n_bg))
	do i = 1,n_bg
		read(222,*) R_bg(i),R_bg_err(i),RD_bg(:,i)
	end do
	close(222)



	! write neutrino floor
	filename1 = trim(mylimitsdir)//'/DL'//trim(inp)//'.txt'
 	open(unit=1000,file=trim(filename1))
 	write(1000,*) 0.0,m_vals

	! Exposure values
	open(unit=333,file=trim(homedir)//'/data/recoils/Ex_'//trim(inp)//'.txt')
	read(333,*) ex_min,ex_max,junk_db
	n_ex = int(junk_db)
	allocate(ex_vals(n_ex))
	ex_vals = logspace(ex_min,ex_max,n_ex)
	close(333)

	write(*,*) '====NuFloor===='
	write(*,*) 'Reading '//trim(inp)
	write(*,*) 'nm =',nm
	write(*,*) 'nTot_bins =',nTot_bins
	write(*,*) 'n_nu =',n_bg
	write(*,*) 'Exposure = ',ex_min,'to',ex_max
	do i = 1,n_ex
		Exposure = ex_vals(i)
		RD_sig = RD_sig*Exposure
		RD_bg = RD_bg*Exposure
		call DiscoveryLimit(m_vals,nm,sigma_min,sigma_max,ns,verbose,DL)
		write(*,*) i,'of',n_ex,'Exposure = ',Exposure,'ton-year||',DL(nm)
		write(1000,*) Exposure,DL
		RD_sig = RD_sig/Exposure
		RD_bg = RD_bg/Exposure
	end do
	close(1000)

	call cpu_time(clock_stop); write(*,*) 'Time elapsed = ',clock_stop-clock_start
	write(*,*) filename1
	write(*,*) '=====DONE====='
	stop
end program runDL
