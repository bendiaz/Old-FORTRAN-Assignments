module setup
	use numtype
	implicit none
	
	integer, parameter :: n_eq = 12
	real(dp), parameter :: grav = 6.673e-11_dp, mass_sun = 1.9891e+30_dp, &
							mass_earth = 5.9736e+24_dp

end module setup


program DiazB_MT1p2

	use setup
	implicit none
	real(dp) :: t, dt
	real(dp), dimension(n_eq) :: y
	real(dp) :: tmax
	
	
	tmax = 60*60*24*365
	
	t = 0._dp
	dt = 60*60*24
	
	
	y(1) = 1.496e+11_dp !x, y, z
	y(2) = 0._dp 
	y(3) = 0._dp
	
	y(4) = 0._dp	!v
	y(5) = 29.783e+3_dp 
	y(6) = 0._dp
	
	!Moon position from Sun
	y(7) = 3.844e+8_dp + 1.496e+11_dp    !x, y, z
	y(8) = 0._dp
	y(9) = 0._dp
	!Moon relative velocity
	y(10) = 0._dp				!v
	y(11) = 2*pi*3.844e+8_dp/(27*24*60*60) + 29.783e+3_dp   
	y(12) = 0._dp
	
	do while(t<tmax)
	!month
		if (t < 60*60*24*30) then
			write (1, *) y(7) , y(8)
			write(2,*) y(1), y(2)
		end if
	
	
		!full year here
		write(3,*) y(1), y(2)
		write(4,*) y(7), y(8)
		
		call rkf45step(t, y, dt)
		
		
	end do	
	print*, 'Done'
end program DiazB_MT1p2