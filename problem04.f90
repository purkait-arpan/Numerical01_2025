program kepler
real, external:: newton_raphson
integer :: i
real :: M,si,r,theta
open(1,file='orbit.dat',status='unknown')

do i = 1, 365
M = (2.0*acos(-1.0)*i)/365.2564
si = newton_raphson(f,df,M,M,0.001)
r = (1.0 - 0.01671022*cos(si))
theta = 2.0*atan((sqrt((1+0.01671022)*(1-0.01671022)))*(tan(si/2.0)))
write(1,*) r, theta
enddo
i = 215
M = (2.0*acos(-1.0)*i)/365.2564
si = newton_raphson(f,df,M,M,0.001)
r = (1.0 - 0.01671022*cos(si))
theta = 2.0*atan((sqrt((1+0.01671022)*(1-0.01671022)))*(tan(si/2.0)))
write(2,*) r, theta


contains
real function f(M,x)
real, intent(in) :: M,x
f = M - 1.0*x + 0.01671022*cos(x)
end function f
real function df(x)
real, intent(in) :: x
df = -1.0 + 0.01671022*cos(x)
end function df
end program kepler

function newton_raphson(f,df,M,a,tol) result(root)
    interface 
        real function f(M,x)
            real, intent(in) :: M,x
        end function f 
        real function df(x)
            real, intent(in) :: x
        end function df
    end interface
    real, intent(in) :: M, a, tol
    real :: x, roo
    x = a
    if (abs(f(M,x)) < tol) then
        root = x
        return
    end if
    10 x = x - f(M,x) / df(x)
    if (abs(f(M,x)) > tol) goto 10
    root = x
end function newton_raphson
