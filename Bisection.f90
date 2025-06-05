program bisection
    real :: a, b, c, x
    f(x) = (x**2.0 - 4.0)
    a = 1.0
    b = 3.5
    if ((f(a)*f(b)) > 0) then
    print*, "not brackted"
    stop
    end if
    if (f(a) == 0.0) then
    print*, "root at a"
    stop
    end if
    if (f(b) == 0.0) then
    print*, "root at b"
    stop
    end if
10  c = (a + b) / 2.0
    if (abs(f(c)) > 1.0e-6) then
        if (f(a)*f(c) < 0) then 
        b = c
        else
        a = c 
        endif
        goto 10
    endif
    print *, "root at c = ", c
end program bisection