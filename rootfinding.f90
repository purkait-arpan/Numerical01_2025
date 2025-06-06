program main
    implicit none
    real, external :: bisection, newton_raphson, secant
    write(*,*) bisection(fi,0.0,7.0,0.0001)
    write(*,*) newton_raphson(fi, dfi, 2.0, 0.0001)
    write(*,*) secant(fi, 0.0, 7.0, 0.0001)

contains
    real function fi(x)
        real, intent(in) :: x
        fi = (x**2 - 4.0)
    end function fi
    real function dfi(x)
        real, intent(in) :: x
        dfi = 2.0 * x
    end function dfi

end program main

real function bisection(f,a,b,tol) result(root)
    interface 
        real function f(x) 
            real, intent(in) :: x
        end function f 
    end interface
    real, intent(in) :: a, b, tol
    real :: p, q, er
    real :: c
    p = a; q = b; er = tol
    if (abs(f(p)) < tol) then
        root = p
        return
    end if
    if (abs(f(q)) < tol) then
        root = q
        return
    end if
    if (f(p) * f(q) > 0.0) then
        print *, "No root found in the interval [", p, ",", q, "]"
        return
    end if
    10  c = (p + q) / 2.0
    if (abs(f(c)) < tol) goto 100
    if (f(c) * f(p) < 0.0) then
        q = c
    else
        p = c
    end if
    goto 10
    100 root = c
end function bisection

real function newton_raphson(f,df,a,tol) result(root)
    interface 
        real function f(x)
            real, intent(in) :: x
        end function f 
        real function df(x)
            real, intent(in) :: x
        end function df
    end interface
    real, intent(in) :: a, tol
    real :: x, er
    x = a; er = tol
    if (abs(f(x)) < tol) then
        root = x
        return
    end if
    10 x = x - f(x) / df(x)
    if (abs(f(x)) < tol) goto 100
    if (abs(f(x)) > tol) goto 10
    100 root = x
end function newton_raphson

real function secant(f,a,b,tol) result(root)
    interface 
        real function f(x)
            real, intent(in) :: x
        end function f 
    end interface
    real, intent(in) :: a, b, tol
    real :: x0, x1, x_new
    x0 = a
    x1 = b
    if (abs(f(x0)) < tol) then
        root = x0
        return
    end if
    if (abs(f(x1)) < tol) then
        root = x1
        return
    end if
    10 x_new = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
    if (abs(f(x_new)) < tol) goto 100
    x0 = x1; x1 = x_new
    goto 10
    100 root = x_new
end function secant