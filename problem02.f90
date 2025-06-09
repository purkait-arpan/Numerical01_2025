program archimedes
real, external :: bisection
real :: p(2), r(5)
r = (/0.5,1.0,1.5,2.0,2.5/)
p = (/0.2,0.7/)
open(1,file='archimedes.dat',status='unknown')
do i = 1, 2
do j = 1, 5
write(1,*) bisection(f,r(j),p(i),0.0,(2.0*r(j)),0.001), r(j), p(i)
enddo
enddo

contains
real function f(h,r,p)
real, intent(in) :: h, r, p
f = h**3 - 3.0*r*h**2 + 4.0*p*r**3
end function  f
end program archimedes

real function bisection(f,r,p,a,b,tol) result(root)
interface 
real function f(h,r,p)
real, intent(in) :: h,r,p
end function f
end interface
real, intent(in) :: r,p,a,b,tol
real :: x0, x1, x_new
x0 = a; x1 = b
if (abs(f(x0,r,p)) < tol) then
root = x0
return
endif 
if (abs(f(x1,r,p)) < tol) then
root = x1
return
endif
if (f(x1,r,p)*f(x0,r,p) > 0.0) then
print*, "Root not Bracketed"
return
end if
10 x_new = (x0+x1) / 2.0
if(abs(f(x_new,r,p)) < tol) goto 100
if(f(x0,r,p)*f(x_new,r,p) < 0.0) then
x1 = x_new
else
x0 = x_new
endif
goto 10
100 root = x_new
end function bisection
