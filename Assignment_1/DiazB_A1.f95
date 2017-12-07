program DiazB_A1
	use numtype
 
	implicit none
 
	real(dp) :: nbe, B 
	real(dp), parameter :: T1=6500
	real(dp), parameter :: T2=300
	real(dp), parameter :: T3 = 2.3
 	real(dp), parameter :: h= 6.62E-34_dp
	real(dp), parameter :: c= 3.00E8_dp
	real(dp), parameter :: kb= 1.38E-23_dp
	real(dp) :: v
	real(dp) :: q = 2* (kb **3) / (c*h)**2
	real(dp) :: Eprime 
	
	
	integer:: i, n, j
 
	
	
	
	do i =1, 10000
		
		Eprime = float(i)*0.001
			
		B = q* T1**3 * Eprime**3 / (exp(Eprime)-1)
			
		v = Eprime*kb*T1/h
			
		write(1,*) Eprime, B
		
	end do
	
	
	do i =1, 10000
		
		Eprime = float(i)*0.001
			
		B = q* T2**3 * Eprime**3 / (exp(Eprime)-1)
			
		v = Eprime*kb*T2/h
			
		write(2,*) Eprime, B
		
	end do
	
	
	
		do i =1, 10000
		
		Eprime = float(i)*0.001
			
		B = q* T3**3 * Eprime**3 / (exp(Eprime)-1)
			
		v = Eprime*kb*T3/h
			
		write(3,*) Eprime, B
		
	end do
		
	

end program 
