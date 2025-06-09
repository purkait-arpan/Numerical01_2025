program realgas
real, external :: secant
print*, secant(f,1.4,1.2,0.0001)
contains
real function f(x)
real, intent(in) :: x
f = x**3 - (0.0023 + (0.003686*273))*x**2 + 0.00874*x - 0.00874*0.0023
end function f
end program realgas
