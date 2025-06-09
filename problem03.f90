program franhofer
real, external :: bisecection
print *, bisection(f,6.28,9.41,0.001)
contains
real function f(x)
real, intent(in) :: x
f = x*cos(x) - sin(x)
end function f
end program franhofer
