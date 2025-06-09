program newton2ndL
integer, parameter :: n = 100
real :: x(n), y(n)
call meuler(f,0.0,0.0,5.5,n,x,y)
do i = 1, n
write(1,*) x(i), y(i)
enddo 
call system('gnuplot -p -e "p ''fort.1'' w p"')
contains
real function f(x,y)
real, intent(in) :: x,y
f = 9.8 - 2.0*y**2
end function f
end program newton2ndL

! Modified Euler Method
subroutine meuler(f,x0,y0,xend,n,x,y)
interface
real function f(x,y)
real, intent(in) :: x,y
end function f
end interface
integer, intent(in)  :: n
real,    intent(in)  :: x0, y0, xend
real,    intent(out) :: x(n), y(n)
real                 :: y_(n), h
h = abs(xend - x0) / n
x(1) = x0
y(1) = y0
y_(1) = y0
do i = 1, n-1
x(i+1) = x(i) + h
y_(i+1) = y(i) + h * f(x(i),y(i))
y(i+1) = y(i) + h * 0.5 * (f(x(i),y(i)) + f(x(i+1),y_(i+1)))
enddo
end subroutine meuler
