! Author: Arpan Purkait
! Date: 2025-06-07
program main
    integer, parameter :: n = 11
    real :: x(n), y(n)
    real, external :: trapz, simpson, trapz2
    x = (/0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 1.9, 2.0/)
    y = (/0.72, 0.78, 0.92, 0.96, 0.98, 1.00, 1.2, 1.22, 1.24, 1.26, 1.27/)

    print *, 'Integral =', trapz(fi, 0.00, 1.00, 0.001)
    print *, 'Integral =', simpson(fi, 0.00, 1.00, 0.001)
    print *, 'Integral =', trapz2(n, x, y)

contains
    real function fi(x)
        real, intent(in) :: x
        fi = x
    end function fi
end program main

real function trapz(f, a, b, h) result(s)
    interface
        real function f(x)
            real, intent(in) :: x
        end function f
    end interface
    real, intent(in) :: a, b, h
    integer :: n, i
    s = 0.0
    n = int(abs(b - a) / h)
    s = f(a) + f(b)
    do i = 1, n - 1
        s = s + 2.0 * f(a + i * h)
    end do
    s = s * h / 2.0
end function trapz

real function simpson(f,a,b,h) result(s)
    interface
    real function f(x)
        real, intent(in) :: x
        end function f 
    end interface
    real, intent(in) :: a, b, h
    integer :: n, i
    s = 0.0
    n = int(abs(b - a) / h)
    s = f(a) + f(b)
    do i = 1, n-1, 2
        s = s + 4.0 * f(a + i * h)
    end do
    do i = 2, n-2, 2
        s = s + 2.0 * f(a + i * h)
    end do
    s = s * h / 3.0
end function simpson

real function trapz2(n,x,y) result(s)
    integer :: n, i
    real :: x(n), y(n)
    s = 0.0
    do i = 1, n - 1
        s = s + (((y(i) + y(i + 1)) * (x(i + 1) - x(i))) / 2.0)
    end do
end function trapz2