program pendulum
real, external :: trapz
print*, trapz(f,0.0,acos(-1.0)/2.0,0.0001)
contains
real function f(x)
real, intent(in) :: x
f = 4.0*sqrt(1.0/9.8) / sqrt(1.0 - (sin(acos(-1.0)*5.0/360.0)*sin(x))**2)
!f = 2.0*sqrt(2.0/9.8)/sqrt(cos(x)-cos(acos(-1.0)*5.0/360.0))
end function f
end program pendulum
