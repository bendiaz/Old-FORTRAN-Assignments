program DiazB_A5
	use numtype
	implicit none
	real(dp) :: m
	integer :: i, j, k, ndim 
	integer, parameter :: n=1000, lwork = 16*n
	real(dp), dimension(n,n) :: x, x2, p, p2, ham, e, a
	integer :: info
	real(dp) :: work(lwork)
	real(dp), dimension(n) :: w
	real(dp), dimension(4):: omsq
	

	m =1 
	
	
	do i = 1, n !top to bottom index
		do j =1, n !left to right
			if((i - j) == 1) then
				x(j,i) = sqrt(float(j))
				p(j,i) = -sqrt(float(j))
			else if ((j - i) == 1) then
				p(j,i) = sqrt(float(i)) 
				x(j,i) = sqrt(float(i))
			else
				x(i,j) = 0
				p(i,j) = 0
			end if
		end do
	end do
	
	!Eigenvalues from 1 multiplied by sqrt(omega)
	omsq = (/0.5, 1.0, 2.0, 3.0/)


	! print*, "Matrix X"
	
	! do i = 1, n
	! 	print'(''(''100f10.5'')'')', x(i,1:n)
	! end do

	! print*, "Matrix P"
	
	! do i = 1, n
	! 	print'(''(''100f10.5'')'')', p(i,1:n)
	! end do

	! print*, "Matrix H"
	
	! do i = 1, n
	! 	print'(''(''100f10.5'')'')', ham(i,1:n)
	! end do

	 
	do k = 1, size(omsq)
		ham(1:n,1:n)= -1/(4._dp*m) * matmul(p(1:n,1:n),p(1:n,1:n)) &
		+((m*omsq(k))/4._dp)* matmul(x(1:n,1:n),x(1:n,1:n))
		 a(1:n, 1:n) = ham(1:n,1:n)

		 info = 0

	 	call dsyev('n','u',n,a,n, w, work,lwork,info)
	 	print*, 'Eigenvalues of H for Omega^2 =', omsq(k)
	 	!print*, info
		 print '(100(f10.5, 3x))', w(1:n)
		 print *, 
		 	do j = 0 , n-1
	 			write(100+k,*) j, w(j+1)
				write(200+k,*) (j+1/2._dp)*sqrt(omsq(k))
	 	

	 		end do
	 end do

	 
	 	


end program DiazB_A5