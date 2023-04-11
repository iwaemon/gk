module midpoint_integral
    implicit none
contains

subroutine integral_midpoint(a,b,n,func,sum)
    implicit none
    interface
        function func(x)
            double precision,intent(in) :: x
            doubleprecision :: func
        end function func
    end interface     
    doubleprecision,intent(in) :: a,b
    integer,intent(in) :: n
    doubleprecision,intent(out) :: sum
    integer :: i
    doubleprecision :: x,h
    h = dble((b-a)/n)
    sum = 0d0
    do i = 1,n
        x = a + (dble(i)-0.5d0)*h
        sum = sum + func(x)
    end do
    sum = sum*h
end subroutine integral_midpoint

subroutine integral_midpoint_complex(a,b,n,func,sum)
    implicit none
    interface
        function func(x)
            double precision,intent(in) :: x
            complex(kind(0d0)) :: func
        end function func
    end interface     
    doubleprecision,intent(in) :: a,b
    integer,intent(in) :: n
    complex(kind(0d0)),intent(out) :: sum
    integer :: i
    doubleprecision :: x,h
    h = dble((b-a)/n)
    sum = 0d0
    do i = 1,n
        x = a + (dble(i)-0.5d0)*h
        sum = sum + func(x)
    end do
    sum = sum*h
end subroutine integral_midpoint_complex

end module midpoint_integral