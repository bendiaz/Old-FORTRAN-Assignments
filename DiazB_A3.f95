program DiazB_A3

	use numType
	implicit none
	real(dp) :: x, xmin, xmax, dx
	integer :: mu, i, imax
    real(dp) :: F,  sigma, z
    	
	xmin = -3
	xmax = 8
	imax = 50
	dx =(xmax - xmin)/imax
	mu = 2
	sigma =1.5
	
	do i = 0, imax
	
		x = xmin + i*dx
		z = (x- mu) /(sigma* (2**(.5)))
		F = .5_dp*(1+erf_cf(z))
		
		write(11, *) x, F
	
	end do
	
	contains
	
		function w_cf(z) result(ss)
		
			implicit none
			complex(dp) :: z, ss
			integer :: nn, n
			
			nn = 10000
			ss = 1
			
			do n = nn, 1, -1
			
				ss = n/2._dp/(z-ss)
					  
			
			end do
			
			 if (z == 0) then
	    
	    	ss =1._dp
	    	
	    		else 
	    		ss =iic/sqrt(pi)*1/(z-ss)
	    		
	    	end if
		
		end function w_cf
		
		function erfc_cf(x) result(ss)
		
			implicit none
			real(dp) :: x, ss, xx
			
			xx = abs(x)
			
			ss = exp(-xx**2) * w_cf(iic*xx)			
			
			if(x < 0) ss = 2 - ss
		
		end function erfc_cf
		
		function erf_cf(x) result(ss)
		
			implicit none
			real(dp) :: x, ss
			
			ss = 1 - erfc_cf(x)
			
		end function erf_cf

end program DiazB_A3