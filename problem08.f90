program inertia
write(*,'(F10.8)')trapz(f,-0.2,0.2,0.0001)
contains
real function f(x)
real, intent(in) :: x
real :: p, a, b
a = 0.2; b = 0.1
p = 1.0/(acos(-1.0)*a*b)
f = p * 2.0 * sqrt(b**2*(1.0-(x**2/a**2)))*x**2
end function f
end program inertia
