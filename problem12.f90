program newton2ndL
integer, parameter :: n = 100
real :: x(n), y(n)
call euler(f,0.0,0.0,5.5,n,x,y)
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
