module setup
	use numtype
	implicit none
	
	integer, parameter :: n_eq = 3
	real(dp):: b, a, c

end module setup

program DiazB_Final1
	use setup
	implicit none
	real(dp), dimension(n_eq) :: y, f
	real(dp) :: tmax,  dt, t 
	integer :: n, nmax
	
	b = 0.2_dp
	a = 0.2_dp
	c = 5.7_dp

		
	y(1) = -5._dp !x(t)
	y(2) = 0._dp !y(t)
	y(3) = 0._dp !z(t)
	
	
	tmax = 100
	
	t = 0._dp
	dt = 0.001_dp
	
	do while (t< tmax)
	
		call rk4step(t, y, dt)
		 write(1, *) t, y(1)
		 write(2,*) t, y(2)
		 write(3,*) t, y(3)
		 write(4,*) y(1),y(2), y(3)
		 t = t + dt
	end do
print* , "done"
end program DiazB_Final1