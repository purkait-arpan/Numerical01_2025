program tf
real, external :: secant
print*, secant(f,4.1,4.2,0.001)
contains
real function f(x)
real, intent(in) ::x
f = 25.0*(32.0+100.0*sin(acos(-1.0)/4.0)/5.0)*(1.0 - exp(-1.0*x/5.0)) - 5.0*32.0*x
end function f
end program tf

