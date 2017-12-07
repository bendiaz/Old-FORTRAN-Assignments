program DiazB_A2

	use numtype
	implicit none
	
	real(dp) :: x, omega, hbar, m, Energy
	integer :: i, j, n, k
	real(dp), dimension(0:10) :: psi
	
	
	hbar = 1
	omega = 2
	m = 1
	
	
	do j=0, 10
	
		do i = -5000, 5000
		
			x = float(i)*.001
		
			psi(0)  = (m*omega/hbar/pi)**(.25) * exp(-m* omega* x**2/2*hbar)
			
			do n= 0, 10
			
				psi(n+1) = (2*sqrt(m* omega / hbar)*x)/((2.0*(n+1.0))**(.5))* psi(n) - (float(n)/(n+1.0))**(.5) * psi(n-1)
				
				Energy = hbar*omega*(float(j)+.5)
				
			end do
			
			
		
			write(10+j,*) x, psi(j)+Energy
		
		end do
		
	end do
		
end program DiazB_A2
