! Author: Arpan Purkait
! Date: 2025-06-07
program main
    implicit none
    integer, parameter :: n = 1000
    real :: x(n), y(n)
    integer :: i
    open(unit=1, file='output.txt')
    call euler(dydx, 0.0, 1.0, 2.0, n, x, y)
    do i = 1, n
    write(1,*) x(i), y(i)
    enddo
    close(1)
contains

    real function dydx(x, y)
        real, intent(in) :: x, y
        dydx = (2 - 10 * y) * 10  
    end function dydx

end program main

subroutine euler(f, x0, y0, xEnd, n, x, y)
    interface
        real function f(x,y) 
            real, intent(in) :: x,y
        end function f
    end interface

    real :: x0, y0,xEnd
    integer, intent(in) :: n
    real, intent(out) :: x(n), y(n)
    integer :: i
    real :: h

    h = (abs(xEnd - x0)/n)
    x(1) = x0
    y(1) = y0

    do i = 1, n-1
        x(i + 1) = x(i) + h
        y(i + 1) = y(i) + h * f(x(i), y(i))
    end do

end subroutine euler

subroutine rk2(f, x0, y0, xEnd, n, x, y)
    interface
        real, function f(x,y)
            real, intent(in) :: x, y
        end function f 
    end interface
    real, intent(in) :: x0, y0, xEnd
    integer, intent(in) :: n
    real, intent(out) :: x(n), y(n)
    integer :: i
    real :: h, k1, k2
    h = (abs(xEnd - x0) / n)
    x(1) = x0 
    y(1) = y0
    do i = 1, n-1 
        k1 = h * f(x(i), y(i))
        k2 = h * f(x(i) + h, y(i) + k1)
        x(i + 1) = x(i) + h
        y(i + 1) = y(i) + (k1 + k2) / 2.0
    end do
end subroutine rk2

subroutine rk4(f, x0, y0, xEnd, n, x, y)
    interface
        real function f(x,y)
            real, intent(in) :: x, y
        end function f 
    end interface
    real, intent(in) :: x0, y0, xEnd
    integer, intent(in) :: n
    real, intent(out) :: x(n), y(n) 
    integer :: i 
    real :: h, k1, k2, k3, k4
    h = (abs(xEnd - x0) / n)
    x(1) = x0
    y(1) = y0
    do i = 1, n-1
        k1 = h * f(x(i), y(i))
        k2 = h * f(x(i) + h / 2.0, y(i) + k1 / 2.0)
        k3 = h * f(x(i) + h / 2.0, y(i) + k2 / 2.0)
        k4 = h * f(x(i) + h, y(i) + k3)
        x(i + 1) = x(i) + h
        y(i + 1) = y(i) + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    end do
end subroutine rk4