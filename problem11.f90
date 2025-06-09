program newton2ndL
integer, parameter :: n = 10
real :: x(n), y(n)
call euler(f,0.0,0.0,5.5,n,x,y)
do i = 1, n
print*, x(i), y(i)
enddo 
contains
real function f(x,y)
real, intent(in) :: x,y
f = 9.8
end function f
end program newton2ndL
