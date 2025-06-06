program main
integer, parameter :: n = 7
integer :: i
real:: x(n), y(n)
real, external :: forward, backward, lagrange

x = (/ 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0 /)
y = (/ 0.0061, 0.0123, 0.0234, 0.0424, 0.0738, 0.1124, 0.1992 /)

write(*,'(F4.1,F10.6)') (x(i), y(i) , i = 1, n)

write(*,'(F10.6)') forward(5.0, n, x, y)  
write(*,'(F10.6)') backward(45.0, n, x, y) 
write(*,'(F10.6)') lagrange(5.0, n, x, y) 
write(*,'(F10.6)') lagrange(45.0, n, x, y) 

end program main

function forward(xp,n,x,y) result(yp)
integer :: n, i, j, k
real :: x(n), y(n), d_table(n,n)
real :: xp,yp, h, u, mult
d_table = 0.0
d_table(:,1) = y
h = x(2) - x(1) 
k = 1
do i = 1, n
  if ((xp - x(i)) > 0.0) k = i
enddo
yp = y(k)
u = (xp - x(k)) / h
mult = 1.0

do j = 2, n
  do i = 1, n-j+1
    d_table(i,j) = d_table(i+1,j-1) - d_table(i, j-1)
  enddo
enddo

do j = 2, n
  mult = (mult * (u - real(j-2))) / real(j-1)
  yp = yp + mult * d_table(k,j)
enddo
end function forward

function backward(xp,n,x,y) result(yp)
integer :: n, i, j, k
real :: x(n), y(n), d_table(n,n)
real :: xp,yp, h, u, mult
d_table = 0.0
d_table(:,1) = y
h = x(2) - x(1)
k = 1
do i = 1, n
    if (xp >= x(i)) then
        k = i + 1
    endif
enddo
u = (xp - x(k)) / h
yp = y(k)
mult = 1.0  
do j = 2, n
    do i = n, j, -1
        d_table(i,j) = d_table(i,j-1) - d_table(i-1,j-1)
    enddo
enddo
do i = 2, n
    mult = mult * (u + real(i - 2)) / real(i-1)
    yp = yp + mult * d_table(k, i)
enddo
end function backward


function lagrange(xp,n,x,y) result(yp)
implicit none
integer :: n, i, j
real :: x(n), y(n), xp,yp,co
yp = 0.00
do i = 1, n
  co = y(i)
  do j = 1, n
    if (j /= i) then
      co = co * (xp - x(j)) / (x(i) - x(j))
    endif
  enddo
  yp = yp + co
enddo
end function lagrange













