program DiazB_MT1p1

	use numType
	implicit none
	real(dp) :: x, xmin, xmax, dx
	integer :: i, imax
    real(dp) :: F,  a
    
	
	xmin = -3._dp
	xmax = 3._dp
	imax = 50._dp
	dx =(xmax - xmin)/imax
	a = 0.5_dp
	x = 0._dp
	
	
	
	do i = 0, imax
	
		x = xmin + i*dx
			if (x==0) F = 1._dp	
			if (x > 0) F = inc_gam(x**2)/sqrt(pi)
			if(x < 0) F = 2- (inc_gam(x**2))/sqrt(pi)
				
			write(1, *) x, F				
				
	end do
	
	contains
	
		function inc_gam(z) result(ss)
		
			implicit none
			real(dp) :: z, ss
			integer :: nn, n
			
			nn = 10000
			ss = 0._dp
			
			do n = nn, 1, -1
			
				ss = z + ((n-a)/(1+(n/ss)))
				
				  
	    		end do
	    	
	    			ss =exp(-z)*(z**a)*(1/ss)
			
			
				
		end function inc_gam
	
	
		
end program DiazB_MT1p1