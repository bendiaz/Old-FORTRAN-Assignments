program contourintegral

    use numtype
    implicit none
    real(dp):: weight(100), abscis(100), r0, t
    integer :: nint, ifail, i
    complex(dp) :: zz, z, residue
    complex(dp), external :: ff
    real(dp):: z0


    nint = 20
    ifail = 0
    
    call d01bcf(0, 0._dp, 2*pi,0._dp, 0._dp, nint, weight, abscis, ifail)

    z0 = 0 !point at which we integrate around 
    r0 = 0.01


    residue = 0._dp
    do  i= 1, nint

        t = abscis(i)
        z = r0*exp((0._dp, 1._dp)*t)
        zz = z0 + z 
        residue = residue + (0.5_dp * ff(zz) * z * weight(i))
        

    end do

    print'(''(''2(f9.5)'')'')', residue


end program contourintegral


function ff(z)
    use numtype
    implicit none
    complex(dp) :: z, ff

    ff = 1._dp/tan(pi*z)
    
end function ff