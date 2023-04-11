module mod_test_gk
use gk
use midpoint_integral
implicit none
contains

function test_func(x, a, b) result(retval)
    implicit none
    DOUBLE PRECISION :: x
    DOUBLE PRECISION :: a, b
    DOUBLE PRECISION :: retval
    retval = exp(-x*x) + a + b
end function test_func

subroutine integrate_func_midpoint(sum)
    implicit none
    DOUBLE PRECISION, INTENT(OUT) :: sum
    DOUBLE PRECISION :: a, b
    DOUBLE PRECISION :: x_min, x_max
    INTEGER :: n
    a = 0d0
    b = 0d0
    x_min = -10d0
    x_max = 10d0
    n = 100
    call integral_midpoint(x_min, x_max, n, wrapper, sum)
    contains
    function wrapper(x) result(retval)
        implicit none
        DOUBLE PRECISION,INTENT(IN) :: x
        DOUBLE PRECISION :: retval
        retval = test_func(x, a, b)
    end function wrapper
end subroutine integrate_func_midpoint

subroutine integrate_func_midpoint_write(sum, file_name, file_unit)
    implicit none
    DOUBLE PRECISION, INTENT(OUT) :: sum
    CHARACTER(100), INTENT(IN) :: file_name
    INTEGER, INTENT(IN) :: file_unit
    DOUBLE PRECISION :: a, b
    DOUBLE PRECISION :: x_min, x_max
    INTEGER :: n
    a = 0d0
    b = 0d0
    x_min = -10d0
    x_max = 10d0
    n = 100
    open(file_unit, file= file_name, status= "replace")
    call integral_midpoint(x_min, x_max, n, wrapper, sum)
    close(file_unit)
    contains
    function wrapper(x) result(retval)
        implicit none
        DOUBLE PRECISION,INTENT(IN) :: x
        DOUBLE PRECISION :: retval
        retval = test_func(x, a, b)
        write(file_unit,*) x, retval
    end function wrapper
end subroutine integrate_func_midpoint_write

subroutine integrate_func_gk(sum)
    implicit none
    DOUBLE PRECISION, INTENT(OUT) :: sum
    DOUBLE PRECISION :: a, b
    DOUBLE PRECISION :: x_min, x_max
    DOUBLE PRECISION :: eps
    INTEGER :: ier
    a = 0d0
    b = 0d0
    x_min = -10d0
    x_max = 10d0
    eps = 1d-4
    call dqag_k2(wrapper, x_min, x_max, eps, sum, ier)
    contains
    function wrapper(x) result(retval)
        implicit none
        DOUBLE PRECISION,INTENT(IN) :: x
        DOUBLE PRECISION :: retval
        retval = test_func(x, a, b)
    end function wrapper
end subroutine integrate_func_gk

subroutine integrate_func_gk_write(sum, file_name, file_unit)
    implicit none
    DOUBLE PRECISION, INTENT(OUT) :: sum
    CHARACTER(100), INTENT(IN) :: file_name
    INTEGER, INTENT(IN) :: file_unit
    DOUBLE PRECISION :: a, b
    DOUBLE PRECISION :: x_min, x_max
    DOUBLE PRECISION :: eps
    INTEGER :: ier
    a = 0d0
    b = 0d0
    x_min = -10d0
    x_max = 10d0
    eps = 1d-4
    open(file_unit, file= file_name, status= "replace")
    call dqag_k2(wrapper, x_min, x_max, eps, sum, ier)
    close(file_unit)
    contains
    function wrapper(x) result(retval)
        implicit none
        DOUBLE PRECISION,INTENT(IN) :: x
        DOUBLE PRECISION :: retval
        retval = test_func(x, a, b)
        write(file_unit,*) x, retval
    end function wrapper
end subroutine integrate_func_gk_write

function func_osc(z) result(retval)
    implicit none
    complex(kind(0d0)), intent(in) :: z
    complex(kind(0d0)) :: retval
    retval = 
end function func_osc

end module mod_test_gk

program main
use mod_test_gk
implicit none
DOUBLE PRECISION, PARAMETER :: pi = 4d0*atan(1d0)
DOUBLE PRECISION :: sum_mid, sum_gk
CHARACTER(100) :: file_name_mid, file_name_gk
file_name_mid = "mid.dat"
file_name_gk = "gk.dat"


call integrate_func_midpoint(sum_mid)
call integrate_func_gk(sum_gk)

call integrate_func_midpoint_write(sum_mid, file_name_mid ,10)
call integrate_func_gk_write(sum_gk, file_name_gk, 11)

write(*,*) sum_mid, sum_gk, sqrt(pi)

end program main