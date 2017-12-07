module setup
	use numtype
	implicit none
	
	integer, parameter :: n_eq = 3
	real(dp):: q, a, omega_d

end module setup

program DiazB_A4
	use setup
	implicit none
	real(dp), dimension(n_eq) :: y
	real(dp) :: tmax, amin, amax, da, dt, t 
	integer :: n, nmax
	
	q = 2._dp
	a = 4._dp
	omega_d = 2._dp/3
	
	nmax = 600
	amin = .9_dp 
	amax = 1.6_dp 
	da = (amax - amin) / nmax
	
	tmax = 600*2*pi
	
	t = 0._dp
	dt = 0.01_dp
	
	do n = 0, nmax
	
		a =amin+ n*da
		
		t= 0._dp
		
		y(1) = 0._dp !omega
		y(2) = 0.2 !theta
		y(3) = 0._dp !phi - driving phase
	
		do while (t<tmax) 
	
			
			if (t > 9*tmax/10) then
		
				if(mod(y(3), 2*pi) <= dt)  write(1, *) a, y(1)
				
			end if
				call rk4step(t, y, dt)
				
		end do
	end do
print* , "done"
end program DiazB_A4