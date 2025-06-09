real function forward(xp,n,x,y) result(yp)
integer, intent(in) :: n
integer :: i, j, k
real, intent(in) :: x(n), y(n)
real, intent(in) :: xp
real :: d_table(n,n)
real :: mult, h, u
d_table = 0.0; d_table(:,1) = y
h = x(2) - x(1)
do i = 1, n
if (xp >= x(i)) k = i
enddo
u = (xp - x(k)) / h
do i = 2, n
do j = 1, n-i+1
d_table(j,i) = d_table(j+1,i-1) - d_table(j,i-1)
enddo
enddo
mult = 1.0
yp = y(k)
do i = 2, n
mult = mult * (u - real(i-2)) / real(i-1)
yp = yp + mult * d_table(k,i)
enddo
end function forward

real function backward(xp,n,x,y) result(yp)
integer, intent(in) :: n
integer :: i,j,k
real, intent(in) :: x(n), y(n)
real, intent(in) :: xp
real :: d_table(n,n)
real :: mult, h, u
d_table = 0.0; d_table(:,1) = y
h = x(2) - x(1)
do i = 1, n
if (xp >= x(i)) k = i+1
enddo
u = (xp - x(k)) / h
do i = 2, n
do j = n, i, -1
d_table(j,i) = d_table(j,i-1) - d_table(j-1,i-1)
enddo
enddo
mult = 1.0
yp = y(k)
do i = 2, n
mult = mult * (u + real(i-2)) / real(i-1)
yp = yp + mult*d_table(k,i)
enddo
end function backward

real function lagrange(xp,n,x,y) result(yp)
integer, intent(in) :: n
integer :: i, j
real, intent(in) :: x(n), y(n)
real :: co
yp = 0.0
do i = 1,n
co = y(i)
do j = 1,n
if (i /= j) co = co * (xp - x(j)) / (x(i) - x(j))
enddo
yp = yp + co 
enddo
end function lagrange
