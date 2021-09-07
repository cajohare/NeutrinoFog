module LabFuncs
	use params
	use util

	implicit none

contains


!================================LabFuncs.f95==================================!
! Contains all functions needed to compute lab dependent and detector dependent effects

! Contents:

! 1. Nuclear Form Factors
! FormFactorHelm: Helm form factor, currently the only one I can be bothered with

! 2. Detector performance
! EnergyResolution: charge detection resolution sig_E [keV]
! Efficiency: charge detection efficiency eff = (0 -> 1)
! LoadDetector: loads readout detector performance data files

! 3. Lab velocity
! LabVelocitySimple: Simple lab velocity in Galactic system without Earth rotation
! LabVelocity: Full lab velocity in (N,W,Z) with Earth rotation
! JulianDay: JulianDay at dd-mm-yyyy hh:hh

! 4. Solar
! EarthSunDistance: Distance between Earth and Sun as a function of time

! 5. Co-ordinate transformations
! gal2eqt: Galactic system to equatorial system
! gal2lab: Galactic system to lab system
! Rgal2lab: Full galactic system to lab system
! lab2gal: Lab system to galactic system
!==============================================================================!






!====================================Form Factors==============================!
!----------------------------------------HELM----------------------------------!
function FormFactorHelm(E_r,A)
	integer :: A
	double precision :: E_r,FormFactorHelm,q,c1,s,R_1
	q = sqrt(2*A*931.5*1000*E_r)*1.0d-12/1.97d-7

	c1 = 1.23d0*A**(1.0d0/3.0d0)-0.6d0
	s = 0.9d0
	R_1 = sqrt(c1**2 + (7.0d0/3.0d0)*pi**2.0d0*(0.52d0**2.0d0) - 5*s**2.0)

	FormFactorHelm = (3*(sin(q*R_1) - q*R_1*cos(q*R_1))*exp(-q*q*s*s/2.0)/(q*R_1)**3)
	if (q.eq.0.0d0) then
	FormFactorHelm = 1.0d0
	end if
end function FormFactorHelm
!------------------------------------------------------------------------------!



!==============================Detector Performanc=============================!
!----------------------------Load Detector Resolution curve--------------------!
function EnergyResolution(E_r,ni) result(sig_E)
	integer :: i,nuc,ni
	double precision :: E_r(ni),sig_E(ni)
	if (energyres_on) then
		do i = 1,ni
			sig_E(i) = E_vals(i)*interp1D(E_vals,energyres_data(:,1),1000,E_r(i))
		end do
	else
		sig_E = 1.0d0
	end if
end function

!---------------------------------Load  Efficiency curve-----------------------!
function Efficiency(E_r,ni) result(eff)
	integer :: i,nuc,ni
	double precision :: E_r(ni),eff(ni),s1,Ec,EC2
	if (efficiency_on) then
		do i = 1,ni
			eff(i) = interp1D(E_vals,efficiency_data(:,1),1000,E_r(i))
		end do
	else
		eff = 1.0d0
	end if
end function

!----------------------------Load Everything-----------------------------------!
subroutine LoadDetector(ro, DetectorName)
	integer :: ro,i
	character(len=100) :: DetectorName
	double precision :: Ei

	! Set readout name
	if (ro.eq.0) then
		DetectorName = 'Ideal'
	elseif (ro.eq.1) then
		DetectorName = 'Pixel'
	end if

	! Load all the data files needed
	open(unit=1300,file='../../data/detectors/energyres/'//trim(DetectorName)//'-EnergyRes.txt')
	open(unit=1301,file='../../data/detectors/efficiency/'//trim(DetectorName)//'-Efficiency.txt')

	! Allocate all data to the right arrays
	do i = 1,1000
		read(1300,*) E_vals(i),energyres_data(i,:)
		read(1301,*) Ei,efficiency_data(i,:)
	end do

	! Remeber to close the file or you will let "him" out
	close(1300)
	close(1301)

end subroutine
!------------------------------------------------------------------------------!



!=================================Lab Velocity=================================!
!------------------------------------------------------------------------------!
function LabVelocitySimple(day) result(vlab)
	! Simple version in Galactic system without Earth ro.
	double precision :: day,e1(3),e2(3),t1,w,vlab(3)
	w = 2*pi/365.0d0
	t1 = 79.0d0
	e1 = (/0.9941d0,0.1088d0,0.0042d0/)
	e2 = (/-0.0504d0,0.4946d0,-0.8677d0/)
	vlab = vv_earthrev*(cos(w*(day-t1))*e1 + sin(w*(day-t1))*e2) ! Earth rev.
	vlab = vlab + v_pec ! Add peculiar velocity
	vlab(2) = vlab(2) + v_LSR ! add LSR velocity
end function LabVelocitySimple

!------------------------------------------------------------------------------!
function LabVelocity(day) ! More complex version in Lab frame system
	double precision :: day,JD,LabVelocity(3),UT,MJD,T_0,t_GAST,t_lab
	double precision :: vv_galrot,v_galrot(3),v_solar(3),v_earthrot(3)
	double precision :: e,lambda_0,L,g,lambda_sun,beta(3),lambda_i(3),v_earthrev(3)

	! Convert day into phase of Earth rotation t_lab
	JD = day+Jan1
	UT = 24*(JD+0.5-floor(JD+0.5)) ! Universal time
	MJD = JD - 2400000.5 ! Modified Julian Day
	T_0 = (floor(MJD)-55197.5)/36525.0
	t_GAST = (101.0308 + 36000.770*T_0 + 15.04107*UT)/15.0
	t_lab = t_GAST + long/15
	t_lab = 15*t_lab ! Lab time in degrees

	! Galactic (LSR) Rotation
	vv_galrot = v_LSR
	v_galrot = (/0.0d0,vv_galrot,0.0d0/)
	call gal2lab(v_galrot,t_lab) ! transform to lab co-ords

	! Peculiar solar Motion
	v_solar = v_pec
	call gal2lab(v_solar,t_lab) ! transform to lab co-ords

	! Earth's revolution (first calculate in galactic frame then transform)
	e = eccentricity
	lambda_0 = orb_long_ecliptic
	L = 281.0298 + 36000.77*T_0 + 0.04107*UT
	g = 357.9258 + 35999.05*T_0 + 0.04107*UT
	lambda_sun = L + (1.915 - 0.0048*T_0)*sin(g*pi/180.0)&
				+ 0.020*sin(2*g*pi/180.0)
	beta = lat_ecl_gal
	lambda_i = long_ecl_gal
	v_earthrev = vv_earthrev*(1-e*sin(pi/180.0*(lambda_sun-lambda_0)))*&
				(cos(beta*pi/180.0)*sin(pi/180.0*(lambda_sun-lambda_i)))
	call gal2lab(v_earthrev,t_lab) ! transform to lab co-ords

	! Earth's rotation
	v_earthrot = 0.465102*cos(lat*pi/180)*(/0.0,-1.0,0.0/) ! already in lab co-ords

	! Total
	LabVelocity = v_earthrot+v_earthrev+v_solar+v_galrot

end function LabVelocity

!------------------------------------------------------------------------------!
function JulianDay(month,day,year,hour) ! calculates Julian day from input date
	integer :: month,day,year,year_r,month_r
	double precision :: hour,JulianDay
	year_r = year+4800-floor((14-month)/12.0)
	month_r = month+12*floor((14-month)/12.0)-3
	JulianDay = day + floor((153*month_r+2)/5.0) + 365*year_r &
				+ floor(year_r/4.0) -floor(year_r/100.0) &
				+ floor(year_r/400.0) - 32045 + (hour-12.0)/24.0
end function JulianDay
!------------------------------------------------------------------------------!





!===================================Solar direction============================!
!------------------------------------------------------------------------------!
function EarthSunDistance(day) result(r_es) ! Earth-sun distance at time = day
	double precision :: day,r_es,g,D
	double precision :: JD
	JD = day+Jan1
	D = JD-2451545.0
	g = 357.529 + 0.98560028*D
	g = g*pi/180
	r_es = 1.00014 - 0.01671*cos(g) - 0.00014*cos(2*g)
	r_es = r_es*AstronomicalUnit
end function  EarthSunDistance

!==============================================================================!







!=================================Co-ord transformations=======================!
!------------------------------------------------------------------------------!
subroutine eqt2lab(v,t_lab) ! Equatorial (x_e,y_e,z_e) to Laboratory (N,W,Z)
	double precision:: v(3),t_lab,t,vp(3),latr
	t = t_lab*pi/180.0
	latr = lat*pi/180.0
	vp = v
	v(1) = -cos(t)*sin(latr)*vp(1) - sin(t)*sin(latr)*vp(2) + cos(latr)*vp(3)
	v(2) = sin(t)*vp(1) - cos(t)*vp(2)
	v(3) = cos(t)*cos(latr)*vp(1) + cos(latr)*sin(t)*vp(2) + sin(latr)*vp(3)
end subroutine eqt2lab

!------------------------------------------------------------------------------!
subroutine gal2eqt(v) ! Galactic (x_g,y_g,z_g) to Equatorial (x_e,y_e,z_e)
	double precision :: v(3),vp(3)
	vp = v
	v(1) = -0.06699*vp(1) + 0.4927*vp(2) - 0.8676*vp(3)
	v(2) = -0.8728*vp(1) -0.4503*vp(2) -0.1884*vp(3)
	v(3) = -0.4835*vp(1) + 0.7446*vp(2) + 0.4602*vp(3)
end subroutine gal2eqt

!------------------------------------------------------------------------------!
subroutine gal2lab(v,t_lab) ! Galactic (x_g,y_g,z_g) to Laboratory (N,W,Z)
	double precision :: v(3),t_lab
	call gal2eqt(v)
	call eqt2lab(v,t_lab)
end subroutine gal2lab

!------------------------------------------------------------------------------!
function Rgal2lab(v,day) result(vrot)
	! Full Galactic (x_g,y_g,z_g) to Laboratory (N,W,Z)
	double precision :: vrot(3),v(3),t_lab,JD,day
	double precision :: UT,MJD,T_0,t_GAST
	JD = Jan1 + day
	! Lab time conversion
	UT = 24*(JD+0.5-floor(JD+0.5))
	MJD = JD - 2400000.5
	T_0 = (floor(MJD)-55197.5)/36525.0
	t_GAST = (101.0308 + 36000.770*T_0 + 15.04107*UT)/15.0
	t_lab = t_GAST + long/15
	t_lab = 15*t_lab ! DEGREES

	vrot = v
	call gal2eqt(vrot)
	call eqt2lab(vrot,t_lab)
end function

!------------------------------------Inverse of above--------------------------!
subroutine lab2gal(v,day)
	double precision :: day,JD,v(3) ! input
	double precision :: UT,MJD,T_0,t_GAST,t_lab,vp(3),latr
	double precision :: M(3,3),l2e(3,3),e2g(3,3),e2l(3,3),g2e(3,3)
	JD = day+Jan1

	! Lab time conversion
	UT = 24*(JD+0.5-floor(JD+0.5)) ! Universal time
	MJD = JD - 2400000.5 ! Modified Julian Day
	T_0 = (floor(MJD)-55197.5)/36525.0
	t_GAST = (101.0308 + 36000.770*T_0 + 15.04107*UT)/15.0
	t_lab = t_GAST + long/15
	t_lab = 15*t_lab ! Lab time in degrees

	t_lab = t_lab*pi/180.0
	latr = lat*pi/180.0

	g2e(1,:) = (/-0.06699d0, 0.4927d0, -0.8676d0/)
	g2e(2,:) = (/-0.8728d0,-0.4503d0,-0.1884d0/)
	g2e(3,:) = (/-0.4835d0,0.7446d0,0.4602d0/)
	call inverse(g2e,e2g,3)

	e2l(1,:) = (/-sin(latr)*cos(t_lab), -sin(latr)*sin(t_lab), cos(latr)/)
	e2l(2,:) = (/sin(t_lab), -cos(t_lab),0.0d0/)
	e2l(3,:) = (/cos(latr)*cos(t_lab), cos(latr)*sin(t_lab), sin(latr)/)
	call inverse(e2l,l2e,3)
	M =  matmul(e2g,l2e)
	vp = v
	v(1) = M(1,1)*vp(1) + M(1,2)*vp(2) + M(1,3)*vp(3)
	v(2) = M(2,1)*vp(1) + M(2,2)*vp(2) + M(2,3)*vp(3)
	v(3) = M(3,1)*vp(1) + M(3,2)*vp(2) + M(3,3)*vp(3)
end subroutine lab2gal
!------------------------------------------------------------------------------!




end module LabFuncs
