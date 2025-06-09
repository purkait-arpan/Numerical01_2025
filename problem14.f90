program cr
real, external :: trapz2
real :: x(100),y(100)
call rk2(f,0.0,10.0,5.0,100,x,y)
print*, trapz2(100,x,y)
contains
real function f(x)
real, intent(in) :: x
f = (10.0 - x / 0.0001) / 10.0
end function f
end program cr
