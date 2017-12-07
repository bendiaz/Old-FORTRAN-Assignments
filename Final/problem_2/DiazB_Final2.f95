module confusion
	use numtype
	implicit none
	real(dp), parameter :: hbar = 1._dp , omega = 1._dp, m = 1._dp, ombee = 2._dp
	integer, parameter :: n = 5
	

end module confusion

program DiazB_MT2

	use integr
	use confusion
	implicit none
	real(dp) :: x
	integer :: i, j, k
	real(dp) :: psiduck, psiduck2
	real(dp), dimension(0:n,0:n) :: Energy, Vmatrix ,identity, Umatrix, Tmatrix, H, &
	pert, NewHam, oldham, hermitian
	real(dp), dimension(n+1, n+1) :: Ham, Cmat
	complex(dp),dimension(n+1,n+1)::expmat, EigenMat, Xmat
	real(dp), dimension(1:maxint) :: xint, wint
	real(dp) :: a, b, res1, resV, resU, orthogonal, t
    integer :: nint, ifail, info
	integer, parameter :: lwork = 16*n
	real(dp) :: work(lwork)
	real(dp), dimension(1:n+1) :: w
	
	
	t = 5._dp
    nint = 300
    b=10._dp
    a = -10._dp

	do i = 0, n 
		do j =0,n 
			if (i-j==0) then
			 	Energy(i,j) = (float(i)+0.5_dp)*hbar*omega

			else
				Energy(i,j) = 0
			end if
		end do
	end do

	do i = 0, n 
		do j =0, n 
			if (abs(i-j) == 5) then
			 	pert(j,i) = hbar*ombee /2
			else 		  
				pert(j,i) = 0 
			end if
		end do
	end do

	call d01bcf (0,a,b,0._dp,0._dp, nint, wint, xint, ifail)
		!call cc11(nint, a,b, xint, wint)

	do i =0, n
		do j =0, n 
			resV = 0._dp
			resU = 0._dp
			res1 = 0._dp

			

				 do k =0, nint
				 	resV = resV+wint(k)*psi2(xint(k),i) * well(xint(k))*psi2(xint(k),j) !ombee
				 	resU = resU+wint(k)*psi(xint(k),i) * potential(xint(k))*psi(xint(k),j) ! omega
					res1 = res1+wint(k)*psi(xint(k),i) *1._dp*psi(xint(k),j)
					

				end do
				Umatrix(i,j) = resU
				identity(i,j) = res1
				Vmatrix(i,j) = resV

			end do
		end do

 				Tmatrix(0:n, 0:n) = Energy(0:n,0:n) - Umatrix(0:n,0:n)
				H(0:n, 0:n) = Tmatrix(0:n, 0:n) + Vmatrix(0:n,0:n)
				NewHam(0:n,0:n) = pert(0:n,0:n)+H(0:n,0:n)
				Ham(1:n+1, 1:n+1) = NewHam(0:n,0:n)
				

		! print*, 'Pert:'
		! do i = 0,n
		! 	print'(''(''31f9.5'')'')', pert(i, 0:n)
		! end do


		! print*, 'I:'
		! do i = 0,n
		! 	print'(''(''31f9.5'')'')', identity(i, 0:n)
		! end do

		! print*, 'U:'
		! do i = 0,n
		! 	print'(''(''31(f9.5, 3x)'')'')', Umatrix(i, 0:n)
		! end do

		! print*, 'Energy:'
		! do i = 0,n
		! 	print'(''(''31(f9.5, 3x)'')'')', Energy(i, 0:n)
		! end do

		! print*, 'T:'
		! do i = 0,n
		! 	print'(''(''31(f9.5, 3x)'')'')', Tmatrix(i, 0:n)
		! end do

		! print*, 'Hamiltonian:'
		! do i = 1,n+1
		! 	print'(''(''31(f9.5, 3x)'')'')', Ham(i, 1:n+1)
		! end do

 
 		call dsyev('v','u',n+1,Ham,n+1, w, work,lwork,info)
	 	print*, 'Eigenvalues of H' 
	 	!print*, info
		 print*, 'Omega =' , omega


		 do i = 1, n+1
		 	! if (w(i)<0) then
		print*, 'Eigenenergy:'	
		 print *,i, w(i)
		   print*, 'Eigenstate:'
		   print'(''(''31(f9.5, 3x)'')'')', Ham(i, 1:n+1)
		 	! end if
		end do

		 hermitian(0:n,0:n) = (transpose(Ham(1:n+1,1:n+1))) 

		! print*, 'Is it hermitian?' 
		! 	print'(''(''31(f9.5, x)'')'')', hermitian


			print *, 'completeness'

		Cmat(1:n+1,1:n+1) = matmul(ham(1:n+1,1:n+1),(transpose(ham(1:n+1,1:n+1))))


	
		do i = 1, n
			print'(''(''31(f9.5, 3x)'')'')', Cmat(i,1:n+1)
		end do

		do i = 1,n+1
				EigenMat(i,1:n+1) = exp(-iic*w(i)*t) * Cmat(i,1:n+1)
		
		end do



			  Xmat(1:n+1,1:n+1) =matmul(ham(1:n+1,1:n+1), EigenMat(1:n+1,1:n+1) )
			  expmat(1:n+1,1:n+1) = matmul(xmat(1:n+1,1:n+1), transpose(ham(1:n+1,1:n+1)))

		 print *, 'Exp (-iH(t))'
		 	do i = 1, n+1
		 		print'(10("(",f12.5,",",f12.5,")"),2x)', expmat(i,1:n+1)
			end do


		
			





contains 
	function psi(x,n) 
		
		implicit none 
		integer :: n, i
		real(dp) :: x, r, psi
		real(dp), dimension(0:n) :: psiduck

		r = sqrt(m*omega/hbar)*x

		psiduck(0) = (m*omega/(hbar*pi))**(.25_dp) * exp(-m* omega* x**2/(2*hbar))

		do i = 0, n-1
			psiduck(i+1) = 2*r/sqrt(2*float(i+1))*psiduck(i) - sqrt(float(i)/float(i+1)) * psiduck(i-1)

		end do
		psi = psiduck(n)
			
	end function psi


	function psi2(x,n) 
		
		implicit none 
		integer :: n, i
		real(dp) :: x, r, psi2
		real(dp), dimension(0:n) :: psiduck2

		r = sqrt(m*ombee/hbar)*x

		psiduck2(0) = (m*ombee/(hbar*pi))**(.25_dp) * exp(-m* ombee* x**2/(2*hbar))

		do i = 0, n-1
			psiduck2(i+1) = 2*r/sqrt(2*float(i+1))*psiduck2(i) - sqrt(float(i)/float(i+1)) * psiduck2(i-1)

		end do
		psi2 = psiduck2(n)
			
	end function psi2


	 function potential(x) result(U)
	 	use confusion
		
		
		real(dp):: U
		real(dp) ::  x

		U = 0.5_dp * omega**2 * hbar * x**2

	end function potential


	function well(x) result(V)
	 	use confusion
		
		
		real(dp):: V
		real(dp) ::  x

		V = 0.5_dp * ombee**2 * hbar * x**2

	end function well

	end program DiazB_MT2
