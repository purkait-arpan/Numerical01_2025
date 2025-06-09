program ho
integer, parameter :: n = 1000
real :: t(n),x(n),y(n)
call euler1(f1,f2,0.0,1.0,1.0,10.0,n,t,x,y)
do i = 1, n
write(1,*) t(i),x(i),y(i)
enddo
call system ('gnuplot -p -e "p ''fort.1'' u 1:2 w p, ''fort.1'' u 2:3 w p" ')
contains
real function f1(t,x,y)
real, intent (in) :: t, x, y
f1 = y
end function f1
real function f2(t,x,y)
real, intent (in) :: t, x, y
f2 = -1.0*y - x
end function f2
end program ho



subroutine euler1(f1,f2,t0,x0,y0,tend,n,t,x,y)
interface
real function f1(t,x,y)
real, intent (in) :: t, x, y
end function f1
real function f2(t,x,y)
real, intent (in) :: t, x, y
end function f2
end interface
integer, intent(in) :: n
real,    intent(in) :: t0,x0,y0,tend
real,    intent(out) :: t(n), x(n), y(n)
integer :: i
real    :: h
h = (abs(tend - t0)) / n
t(1) = t0
x(1) = x0
y(1) = y0
do i = 1, n-1
t(i+1) = t(i) + h
x(i+1) = x(i) + h * f1(t(i),x(i),y(i))
y(i+1) = y(i) + h * f2(t(i),x(i),y(i))
enddo 
end subroutine euler1
