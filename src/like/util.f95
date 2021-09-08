module util
	use params
	implicit none

contains


!================================util.f95=============================================!
! Just contains a bunch of boring but useful/essential stuff
! linspace: linearly spaced array of n values from x_lower to x_upper
! logspace: logarithmically spaced array of n values from x_lower to x_upper
! interp1D: interpolates set of points defined by (x,y) at a different set of x1 points
! erf: error function
! erfi: complex error function
! inverse: matrix inverse
!=====================================================================================!



!================================USEFUL FUNCTIONS=====================================!
!------------------------------Linearly Spaced Array----------------------------------!
function linspace(x_lower,x_upper,n) result(V)
	integer :: n,i
	double precision :: x_lower,x_upper,V(n),w
	w = (x_upper-x_lower)/(1.0*n-1.0d0)
	V = (/(i,i=0,n-1)/)*w + x_lower
end function

!------------------------Logarithmically Spaced Array---------------------------------!
function logspace(x_lower,x_upper,n) result(V)
	double precision :: V(n),x_upper,x_lower,wid
	integer :: n,i
	if (n>1) then
		wid = (log10(x_upper)-log10(x_lower))/(1.0d0*n-1.0d0)
		do i = 1,n
			V(i) = x_lower*10.0d0**((i-1)*wid)
		end do
	else
		V(1) = x_lower
	end if
end function

!------------------------------Linear Interpolation--------------------------------!
! TO DO: Merge these two functions for fuck's sake
subroutine interp1(E_full,dRdE_full_s,nbins_full,E_rec,   y)
	integer :: nbins_full
	double precision,dimension(nbins_full) :: E_full,dRdE_full_s
	double precision :: x,xa,xb,ya,yb,y,E_rec
	integer :: ia,ib,iaa(1)
	iaa = minloc(abs(E_full(:)-E_rec))
	write(*,*) E_rec,iaa
	ia = iaa(1)
	if (abs(E_full(ia+1)-E_rec).gt.abs(E_full(ia-1)-E_rec)) then
	ia = ia-1
	end if
	ib = ia+1
	xa = E_full(ia)
	xb = E_full(ib)
	ya = dRdE_full_s(ia)
	yb = dRdE_full_s(ib)
	y = ya + (yb-ya)*(E_rec-xa)/(xb-xa)
	!write(*,*) E_rec,xa,xb,ya,yb,y
end subroutine interp1
function interp1D(x,y,n,x1) result(y1)
	integer :: n,n1,i,i1,i2
	double precision :: x(n),y(n),x1,y1,xi
	xi = x1
	i1 = sum(minloc(abs(x-xi)))
	if (xi.lt.x(i1)) then
		i1 = i1-1
		i2 = i1+1
		y1 = (y(i1)*(x(i2)-xi)+y(i2)*(xi-x(i1)))/(x(i2)-x(i1))
	elseif (xi.gt.x(i1)) then
		i2 = i1+1
		y1 = (y(i1)*(x(i2)-xi)+y(i2)*(xi-x(i1)))/(x(i2)-x(i1))
	else
		y1 = y(i1)
	end if
end function






!================================OTHER  FUNCTIONS=====================================!
!----------------------------- error function-----------------------------------------!
function erf(x)
	double precision :: dumerfc, x, erf
	double precision :: t, z
	z = abs(x)
	t = 1.0d0 / ( 1.0d0 + 0.5d0 * z )
	dumerfc =       t * exp(-z * z - 1.26551223d0 + t *	    &
	   ( 1.00002368d0 + t * ( 0.37409196d0 + t *		&
	   ( 0.09678418d0 + t * (-0.18628806d0 + t *		&
	   ( 0.27886807d0 + t * (-1.13520398d0 + t *		&
	   ( 1.48851587d0 + t * (-0.82215223d0 + t * 0.17087277d0 )))))))))
	if ( x.lt.0.0d0 ) dumerfc = 2.0d0 - dumerfc
	erf = 1.0d0 - dumerfc
end function erf
!-----------------------complex error function----------------------------------------!
! Copied from some library. It does seem to work
! just don't look at if for too long
function erfi(x)
	double precision :: erfi,x
	erfi = (2.0d0/sqrt(pi))*exp(x**2.0d0)*Daw(x)
end function

function DAW(XX) RESULT(fn_val)
	INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
	REAL (dp), INTENT(IN)  :: xx
	REAL (dp)              :: fn_val
	INTEGER    :: i
	REAL (dp)  :: frac, sump, sumq, w2, x, y

	REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
	                         six25 = 6.25_dp, one225 = 12.25_dp, two5 = 25.0_dp
	REAL (dp), PARAMETER  :: XSMALL = 1.05D-08, XLARGE = 9.49D+07,   &
	                         XMAX = 2.24D+307
	REAL (dp), PARAMETER  :: P1(10) = (/  &
	        -2.69020398788704782410D-12, 4.18572065374337710778D-10,  &
	        -1.34848304455939419963D-08, 9.28264872583444852976D-07,  &
	        -1.23877783329049120592D-05, 4.07205792429155826266D-04,  &
	        -2.84388121441008500446D-03, 4.70139022887204722217D-02,  &
	        -1.38868086253931995101D-01, 1.00000000000000000004D+00 /)
	REAL (dp), PARAMETER  :: Q1(10) = (/  &
	         1.71257170854690554214D-10, 1.19266846372297253797D-08,  &
	         4.32287827678631772231D-07, 1.03867633767414421898D-05,  &
	         1.78910965284246249340D-04, 2.26061077235076703171D-03,  &
	         2.07422774641447644725D-02, 1.32212955897210128811D-01,  &
	         5.27798580412734677256D-01, 1.00000000000000000000D+00 /)
	REAL (dp), PARAMETER  :: P2(10) = (/  &
	        -1.70953804700855494930D+00,-3.79258977271042880786D+01,  &
	         2.61935631268825992835D+01, 1.25808703738951251885D+01,  &
	        -2.27571829525075891337D+01, 4.56604250725163310122D+00,  &
	        -7.33080089896402870750D+00, 4.65842087940015295573D+01,  &
	        -1.73717177843672791149D+01, 5.00260183622027967838D-01 /)
	REAL (dp), PARAMETER  :: Q2(9) = (/  &
	         1.82180093313514478378D+00, 1.10067081034515532891D+03,  &
	        -7.08465686676573000364D+00, 4.53642111102577727153D+02,  &
	         4.06209742218935689922D+01, 3.02890110610122663923D+02,  &
	         1.70641269745236227356D+02, 9.51190923960381458747D+02,  &
	         2.06522691539642105009D-01 /)
	REAL (dp), PARAMETER  :: P3(10) = (/  &
	        -4.55169503255094815112D+00,-1.86647123338493852582D+01,  &
	        -7.36315669126830526754D+00,-6.68407240337696756838D+01,  &
	         4.84507265081491452130D+01, 2.69790586735467649969D+01,  &
	        -3.35044149820592449072D+01, 7.50964459838919612289D+00,  &
	        -1.48432341823343965307D+00, 4.99999810924858824981D-01 /)
	REAL (dp), PARAMETER  :: Q3(9) = (/  &
	         4.47820908025971749852D+01, 9.98607198039452081913D+01,  &
	         1.40238373126149385228D+01, 3.48817758822286353588D+03,  &
	        -9.18871385293215873406D+00, 1.24018500009917163023D+03,  &
	        -6.88024952504512254535D+01,-2.31251575385145143070D+00,  &
	         2.50041492369922381761D-01 /)
	REAL (dp), PARAMETER  :: P4(10) = (/  &
	        -8.11753647558432685797D+00,-3.84043882477454453430D+01,  &
	        -2.23787669028751886675D+01,-2.88301992467056105854D+01,  &
	        -5.99085540418222002197D+00,-1.13867365736066102577D+01,  &
	        -6.52828727526980741590D+00,-4.50002293000355585708D+00,  &
	        -2.50000000088955834952D+00, 5.00000000000000488400D-01 /)
	REAL (dp), PARAMETER  :: Q4(9) = (/  &
	         2.69382300417238816428D+02, 5.04198958742465752861D+01,  &
	         6.11539671480115846173D+01, 2.08210246935564547889D+02,  &
	         1.97325365692316183531D+01,-1.22097010558934838708D+01,  &
	        -6.99732735041547247161D+00,-2.49999970104184464568D+00,  &
	         7.49999999999027092188D-01 /)
	x = xx
	IF (ABS(x) > xlarge) THEN
	  IF (ABS(x) <= xmax) THEN
	    fn_val = half / x
	  ELSE
	    fn_val = zero
	  END IF
	ELSE IF (ABS(x) < xsmall) THEN
	  fn_val = x
	ELSE
	  y = x * x
	  IF (y < six25) THEN
	    sump = p1(1)
	    sumq = q1(1)
	    DO  i = 2, 10
	      sump = sump * y + p1(i)
	      sumq = sumq * y + q1(i)
	    END DO
	    fn_val = x * sump / sumq
	  ELSE IF (y < one225) THEN
	    frac = zero
	    DO  i = 1, 9
	      frac = q2(i) / (p2(i)+y+frac)
	    END DO
	    fn_val = (p2(10)+frac) / x
	  ELSE IF (y < two5) THEN
	    frac = zero
	    DO  i = 1, 9
	      frac = q3(i) / (p3(i)+y+frac)
	    END DO
	    fn_val = (p3(10)+frac) / x
	  ELSE
	    w2 = one / x / x
	    frac = zero
	    DO  i = 1, 9
	      frac = q4(i) / (p4(i)+y+frac)
	    END DO
	    frac = p4(10) + frac
	    fn_val = (half + half*w2*frac) / x
	  END IF
	END IF
	RETURN
end function daw



!===================================Matrix Inverse=====================================!
subroutine inverse(a,c,n)
!============================================================
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed
! during the calculation
!===========================================================
	implicit none
	integer n
	double precision a(n,n), c(n,n)
	double precision L(n,n), U(n,n), b(n), d(n), x(n)
	double precision coeff
	integer i, j, k

	! step 0: initialization for matrices L and U and b
	! Fortran 90/95 aloows such operations on matrices
	L=0.0
	U=0.0
	b=0.0

	! step 1: forward elimination
	do k=1, n-1
	  do i=k+1,n
	  coeff=a(i,k)/a(k,k)
	  L(i,k) = coeff
	  do j=k+1,n
	    a(i,j) = a(i,j)-coeff*a(k,j)
	  end do
	  end do
	end do

	! Step 2: prepare L and U matrices
	! L matrix is a matrix of the elimination coefficient
	! + the diagonal elements are 1.0
	do i=1,n
	  L(i,i) = 1.0
	end do
	! U matrix is the upper triangular part of A
	do j=1,n
	  do i=1,j
	U(i,j) = a(i,j)
	  end do
	end do

	! Step 3: compute columns of the inverse matrix C
	do k=1,n
	  b(k)=1.0
	  d(1) = b(1)
	! Step 3a: Solve Ld=b using the forward substitution
	  do i=2,n
	d(i)=b(i)
	do j=1,i-1
	  d(i) = d(i) - L(i,j)*d(j)
	end do
	  end do
	! Step 3b: Solve Ux=d using the back substitution
	  x(n)=d(n)/U(n,n)
	  do i = n-1,1,-1
	x(i) = d(i)
	do j=n,i+1,-1
	  x(i)=x(i)-U(i,j)*x(j)
	end do
	x(i) = x(i)/u(i,i)
	  end do
	! Step 3c: fill the solutions x(n) into column k of C
	  do i=1,n
	c(i,k) = x(i)
	  end do
	  b(k)=0.0
	end do
end subroutine inverse


!-------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------!

end module util
