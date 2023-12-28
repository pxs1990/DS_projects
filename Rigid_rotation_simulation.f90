program cre_data
    implicit none

    integer, parameter :: ni = 50, nj = 50, nk = 50 ! Data dimensions
	real, parameter :: pi = 4.0*atan(1.0)
	
	! Data Generation
	real, parameter :: xi = 1.0, gma = 1.0
	real, parameter :: re = 1.0
	real :: nu, radi
	
	integer :: i,j,k
	real :: h
	real, dimension(3,3) :: asp, sp
	real, dimension(0:ni,0:nj,0:nk) :: u,v,w,x,y,z
	real ::  R
	
	! pds declarations
    real, dimension(0:ni,0:nj,0:nk) :: u_x, v_x, w_x, u_y, v_y, w_y, u_z, v_z, w_z, Q, delta

	! For calc eigenvalues
	real, parameter :: z0(3) = (/0.0, 0.0, 1.0/)
	real, dimension(3) :: vr, vg
	real, dimension(3,3) :: a, tt, qqq 
	real, dimension(0:ni,0:nj,0:nk) :: liutexEigr
	real :: aa, bb, cc, discr, qq, rr, aaaa, bbbb
	real :: del1, del2, del3, temp, beta
	real :: eig3r
	complex :: eig1c, eig2c
	
	! Liutex 
	real, dimension(ni,nj,nk) :: liutex_x, liutex_y, liutex_z,  liutex_ci
	real, dimension(ni,nj,nk) :: liutex_mag, luitex_mag_x, liutex_mag_y, liutex_mag_z
	real, dimension(ni,nj,nk) :: isLocExtrOmg, omgGradXLiutex
	real :: rm, ljmtmp, ljmtmp2
	
	! Creating Subsurface document
	open(21, file='rigid_data.dat', form='formatted', action='write')
	write(21,*) 'variables = "x","y","z","u","v","w","Q","delta"'
	write(21,*) 'zone f=point, I=',ni+1,', J=',nj+1,', K=',nk+1

	h = 1.0 ! Step size
	
    do k=0,nk
        do j=0,nj
            do i=0,ni
				x(i,j,k) = h*(i - ni/2.0)
				y(i,j,k) = h*(j - nj/2.0)
				z(i,j,k) = h*(k - nk/2.0)
				
				u(i,j,k) = y(i,j,k)
				v(i,j,k) = -x(i,j,k)
				w(i,j,k) = 0.0
				
				! Partials
				u_x(i,j,k) = 0.0
				u_y(i,j,k) = 1.0
				u_z(i,j,k) = 0.0
				v_x(i,j,k) = -1.0
				v_y(i,j,k) = 0.0
				v_z(i,j,k) = 0.0
				w_x(i,j,k) = 0.0
				w_y(i,j,k) = 0.0
				w_z(i,j,k) = 0.0

				
				! A matrix
				a(1,1) = u_x(i,j,k)
				a(1,2) = u_y(i,j,k)
				a(1,3) = u_z(i,j,k)
				a(2,1) = v_x(i,j,k)
				a(2,2) = v_y(i,j,k)
				a(2,3) = v_z(i,j,k)
				a(3,1) = w_x(i,j,k)
				a(3,2) = w_y(i,j,k)
				a(3,3) = w_z(i,j,k)
				

				! Frobenius norm of symmetric(sp) and antisymmetric parts(asp)
				asp = 0.5*(a-transpose(a))
				sp = 0.5*(a+transpose(a))
				
				! Calculating Q
				Q(i,j,k) = 0.5*(norm2(asp)**2 - norm2(sp)**2)

				! Calculating R for delta
				R = -u_x(i,j,k)*(v_y(i,j,k)*w_z(i,j,k) - w_y(i,j,k)*v_z(i,j,k)) +		&
					u_y(i,j,k)*(v_x(i,j,k)*w_z(i,j,k) - w_x(i,j,k)*v_z(i,j,k)) -		&
					u_z(i,j,k)*(v_x(i,j,k)*w_y(i,j,k) - w_x(i,j,k)*v_y(i,j,k))
				! Calculating delta
				delta(i,j,k) = (Q(i,j,k)/3)**3 + (R/2)**2
				
				!==================================================================================
				! Calculate the eigenvalues and eigenvector of delta, Q, 
				!CHARARISTISTIC EQUATION:=	coef_1 + coef_2*lambda + coef_3*lambda**2 - lambda**3
				!==================================================================================
				
				! coef_1 = u_z*v_x*w_y + u_y*v_z*w_z - u_z*v_y*w_z + u_x*v_y*w_z - u_x*v_z*w_y - u_y*v_x*w_z
				! coef_2 = u_y*v_x - u_x*w_z - u_x*v_y + u_z*w_z + v_z*w_y - v_y*w_z 
				! coef_3 = v_y + w_z + u_x 
				
				
				
				! Writing x, y, z, u, v, w, Q, delta to file
				write(21,100) x(i,j,k),y(i,j,k),z(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),Q(i,j,k),delta(i,j,k)
				
				! Partial derivatives
				write(22,101) u_x(i,j,k),u_y(i,j,k),u_z(i,j,k),v_x(i,j,k),v_y(i,j,k),v_z(i,j,k),w_x(i,j,k),w_y(i,j,k),w_z(i,j,k)
            end do
        end do
    end do

	100 format(8f19.10)
	101 format(9f13.13)
	
end program


function char_eq(coeff_1, coeff_2, coeff3)
	
	implicit none

end function char_eq


! calculate the rortex
! a velocity gradient tensor, vor: vorticity,vr: rotational axis, rorMag:rotational strength
subroutine cal_rortex(a, vor, vr, rorMag)
	implicit none

	real, intent(in) :: a(3,3)
	real, intent(in) :: vor(3)
	real, intent(out) :: vr(3), rorMag

	real:: aa,bb,cc,delta,rr,aaaa,bbbb,qq,delta1,delta2,delta3,temp
	complex:: eig1c,eig2c,eig3r
	real :: tt(3,3)
	! Cubic Formula
	! Reference: Numerical Recipes in FORTRAN 77, Second Edition
	! 5.6 Quadratic and Cubic Equations
	! Page 179
	!---------------------------------------------------------------------

	! cubic equation
	! x**3 + aa * x**2 + bb * x + cc = 0

	! coefficients of characteristic equation

	aa = -(a(1,1)+a(2,2)+a(3,3))

	tt = matmul(a,a)

	bb = -0.5*(tt(1,1)+tt(2,2)+tt(3,3)-(a(1,1)+a(2,2)+a(3,3))**2)

	cc = -(a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) 	&
		-a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))		&
		+a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1)))

	! discriminant of characteristic equation
	delta = 18*aa*bb*cc-4*aa**3*cc+aa**2*bb**2-4*bb**3-27*cc**2

	qq = (aa**2-3*bb)/9.0
	rr = (2*aa**3-9*aa*bb+27*cc)/54.0

	! delta = rr**2 - qq**3
	! alleviate round error
	delta = -delta/108
	vr=0.0
	rorMag=0.0

	if(delta > 0.0) then ! one real root and two complex conjugate roots

		aaaa = -sign(1.0, rr)*(abs(rr)+sqrt(delta))**(1.0/3.0)

		if(aaaa == 0.0) then
			bbbb = 0.0
		else
			bbbb = qq/aaaa
		end if

		eig1c = cmplx(-0.5*(aaaa+bbbb)-aa/3.0, 0.5*sqrt(3.0)*(aaaa-bbbb))
		eig2c = cmplx(real(eig1c), -aimag(eig1c)) !original wrong here
		eig3r = aaaa+bbbb-aa/3.0

		! real right eigenvector

		delta1 = (a(1,1)-eig3r)*(a(2,2)-eig3r) - a(2,1)*a(1,2)
		delta2 = (a(2,2)-eig3r)*(a(3,3)-eig3r) - a(2,3)*a(3,2)
		delta3 = (a(1,1)-eig3r)*(a(3,3)-eig3r) - a(1,3)*a(3,1)

		if(delta1 == 0.0 .and. delta2 == 0.0 .and. delta3 == 0.0) then
			write(*,*) 'ERROR: delta1 = delta2 = delta3 = 0.0'
			write(*,*) a(1,1)-eig3r,  a(1,2),       a(1,3)
			write(*,*) a(2,1),        a(2,2)-eig3r, a(2,3)
			write(*,*) a(3,1),        a(3,2),       a(3,3)-eig3r
			!  write(*,*) i, j, k
			stop
		end if

		if(abs(delta1) >= abs(delta2) .and. abs(delta1) >= abs(delta3)) then

		  vr(1) = (-(a(2,2)-eig3r)*a(1,3) + a(1,2)*a(2,3))/delta1
		  vr(2) = (a(2,1)*a(1,3) - (a(1,1)-eig3r)*a(2,3))/delta1
		  vr(3) = 1.0

		else if(abs(delta2) >= abs(delta1) .and. abs(delta2) >= abs(delta3)) then

		  vr(1) = 1.0
		  vr(2) = (-(a(3,3)-eig3r)*a(2,1) + a(2,3)*a(3,1))/delta2
		  vr(3) = (a(3,2)*a(2,1) - (a(2,2)-eig3r)*a(3,1))/delta2

		else if(abs(delta3) >= abs(delta1) .and. abs(delta3) >= abs(delta2)) then

		   vr(1) = (-(a(3,3)-eig3r)*a(1,2) + a(1,3)*a(3,2))/delta3
		   vr(2) = 1.0
		   vr(3) = (a(3,1)*a(1,2) - (a(1,1)-eig3r)*a(3,2))/delta3

		else

			write(*,*) 'ERROR: '
			write(*,*) delta1, delta2, delta3
			stop

		end if

		temp = sqrt(vr(1)**2+vr(2)**2+vr(3)**2)

		vr(1) = vr(1)/temp
		vr(2) = vr(2)/temp
		vr(3) = vr(3)/temp
		temp=dot_product(vor,vr)
		vr(1)=sign(1.0,temp)*vr(1)
		vr(2)=sign(1.0,temp)*vr(2)
		vr(3)=sign(1.0,temp)*vr(3)

		rorMag=(temp - sqrt(temp**2-4*aimag(eig2c)*aimag(eig2c)))
	end if

end subroutine cal_rortex

subroutine rotation(u, v, r)
!-------------------------------------------------------------------------------
! calculate rotation matrix r which rotates unit vector u to unit vector v
!-------------------------------------------------------------------------------

	implicit none

	real, intent(in) :: u(3)
	real, intent(in) :: v(3)
	real, intent(out) :: r(3, 3)

	real :: a(3)
	real :: aa
	real :: t
	real :: alpha
	real :: c, s

	real, parameter :: eps = 1.0e-10

	! a = u x v
	a(1) = u(2)*v(3)-u(3)*v(2)
	a(2) = u(3)*v(1)-u(1)*v(3)
	a(3) = u(1)*v(2)-u(2)*v(1)

	! norm
	aa = sqrt(a(1)**2+a(2)**2+a(3)**2)

	if(aa < eps) then

		r = 0.0
		r(1,1) = 1.0
		r(2,2) = 1.0
		r(3,3) = 1.0

	else

		a = a/aa
		t = u(1)*v(1)+u(2)*v(2)+u(3)*v(3)
		if(t > 1.0) t = 1.0
		if(t < -1.0) t = -1.0
		alpha = acos(t)

		c = cos(alpha)
		s = sin(alpha)

		r(1,1) = a(1)**2*(1-c)+c
		r(1,2) = a(1)*a(2)*(1-c)-a(3)*s
		r(1,3) = a(1)*a(3)*(1-c)+a(2)*s

		r(2,1) = a(2)*a(1)*(1-c)+a(3)*s
		r(2,2) = a(2)**2*(1-c)+c
		r(2,3) = a(2)*a(3)*(1-c)-a(1)*s

		r(3,1) = a(3)*a(1)*(1-c)-a(2)*s
		r(3,2) = a(3)*a(2)*(1-c)+a(1)*s
		r(3,3) = a(3)**2*(1-c)+c

	end if

end subroutine rotation