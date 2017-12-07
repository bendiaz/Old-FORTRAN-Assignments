subroutine rk4step(t, y, dt)


	
	use setup
	implicit none
	
	real(dp), intent(inout) :: t ! these variables can either be the input or output; can be changed later
	real(dp), intent(in) :: dt ! dt can only be an input
	real(dp), dimension(n_eq), intent(inout) :: y ! y becomes a vector of two
	real(dp), dimension(n_eq) :: k1, k2, k3, k4, dy
	
	
	! this part is all the math that calls the physics from the sub-subroutine
	call deriv(t, y, dt, k1) ! calling at certain t's and y's to get k1 with dt "step size" (h = dt)
	call deriv(t + dt/2, y + k1/2, dt, k2) ! k2 from notes and so on
	call deriv(t + dt/2, y + k2/2, dt, k3)
	call deriv(t + dt, y + k3, dt, k4)
	
	dy = (k1 + 2*k2 + 2*k3 + k4)/6
	
	t = t + dt ! new t
	y = y + dy ! new y
	
	
	contains
	
		
		subroutine deriv(t, y, dt, k) ! the physics of how to actually do a derivative
		
			
			
			implicit none
			
			real(dp), intent(in) :: t, dt ! for the derivative, we want t and y to only be inputs
			real(dp), dimension(n_eq), intent(in) :: y
			real(dp), dimension(n_eq), intent(out) :: k
			real(dp), dimension(n_eq) :: f
			real(dp):: b, a, c

			b = 0.2_dp
			a = 0.2_dp
			c = 5.7_dp


			f(1) = -(y(2)+y(3))
			f(2) = y(1)+ a*y(2)
			f(3) = b + y(3)*(y(1)-c)
			
			
			k = f*dt
			
		
		
		end subroutine deriv
	


end subroutine rk4step