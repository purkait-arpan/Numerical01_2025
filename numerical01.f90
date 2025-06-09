! Arpan Purkait
real function bisection(f,a,b,tol) result(root)
    interface 
        real function f(x) 
            real, intent(in) :: x
        end function f 
    end interface
    real, intent(in) :: a, b, tol
    real :: x0, x1
    real :: c
    x0 = a
    x1 = b
    if (abs(f(x0)) < tol) then
        root = x0
        return
    end if
    if (abs(f(x1)) < tol) then
        root = x1
        return
    end if
    if (f(x0) * f(x1) > 0.0) then
        print *, "No root found in the interval [", a, ",", b, "]"
        return
    end if
    10  c = (x0 + x1) / 2.0
    if (abs(f(c)) < tol) goto 100
    if (f(c) * f(x0) < 0.0) then
        x1 = c
    else
        x0 = c
    end if
    goto 10
    100 root = c
    end function bisection

real function newton_raphson(f,df,a,tol) result(root)
    interface 
        real function f(x)
            real, intent(in) :: x
        end function f 
        real function df(x)
            real, intent(in) :: x
        end function df
    end interface
    real, intent(in) :: a, tol
    real :: x
    x = a
    if (abs(f(x)) < tol) then
        root = x
        return
    end if
    10 x = x - f(x) / df(x)
    if (abs(f(x)) > tol) goto 10
    root = x
    end function newton_raphson

real function secant(f,a,b,tol) result(root)
    interface 
        real function f(x)
            real, intent(in) :: x
        end function f 
    end interface
    real, intent(in) :: a, b, tol
    real :: x0, x1, x_new
    x0 = a
    x1 = b
    if (abs(f(x0)) < tol) then
        root = x0
        return
    end if
    if (abs(f(x1)) < tol) then
        root = x1
        return
    end if
    10 x_new = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
    x0 = x1; x1 = x_new
    if (abs(f(x_new)) > tol) goto 10
    root = x_new
    end function secant

real function forward(xp,n,x,y) result(yp)
  integer, intent(in) :: n
  real, intent(in) :: x(n), y(n)
  real, intent(in) :: xp
  real :: d_table(n,n)
  integer :: i, j, k
  real :: h, u, mult
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

real function backward(xp,n,x,y) result(yp)
  integer, intent(in) :: n
  real, intent(in) :: x(n), y(n)
  real, intent(in) :: xp
  real :: d_table(n,n)
  integer :: i, j, k
  real :: h, u, mult
  d_table = 0.0
  d_table(:,1) = y
  h = x(2) - x(1)
  k = 1
  do i = 1, n
      if (xp >= x(i)) k = i + 1
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

real function lagrange(xp,n,x,y) result(yp)
  integer, intent(in) :: n
  real, intent(in) :: x(n), y(n)
  real, intent(in) :: xp
  integer :: i, j
  real :: co
  yp = 0.00
  do i = 1, n
    co = y(i)
    do j = 1, n
      if (j /= i) co = co * (xp - x(j)) / (x(i) - x(j))
    enddo
    yp = yp + co
  enddo
    end function lagrange

real function trapz(f, a, b, h) result(s)
    interface
        real function f(x)
            real, intent(in) :: x
        end function f
    end interface
    real, intent(in) :: a, b, h
    integer :: n, i
    s = 0.0
    n = int(abs(b - a) / h)
    s = f(a) + f(b)
    do i = 1, n - 1
        s = s + 2.0 * f(a + i * h)
    end do
    s = s * h / 2.0
    end function trapz

real function simpson(f,a,b,h) result(s)
    interface
    real function f(x)
        real, intent(in) :: x
        end function f 
    end interface
    real, intent(in) :: a, b, h
    integer :: n, i
    s = 0.0
    n = int(abs(b - a) / h)
    s = f(a) + f(b)
    do i = 1, n-1, 2
        s = s + 4.0 * f(a + i * h)
    end do
    do i = 2, n-2, 2
        s = s + 2.0 * f(a + i * h)
    end do
    s = s * h / 3.0
    end function simpson

real function trapz2(n,x,y) result(s)
    integer, intent(in) :: n
    real, intent(in) :: x(n), y(n)
    integer :: i
    s = 0.0
    do i = 1, n - 1
        s = s + (((y(i) + y(i + 1)) * (x(i + 1) - x(i))) / 2.0)
    end do
    end function trapz2

subroutine euler(f, x0, y0, xEnd, n, x, y)
    interface
        real function f(x,y) 
            real, intent(in) :: x,y
        end function f
    end interface

    real, intent(in) :: x0, y0,xEnd
    integer, intent(in) :: n
    real, intent(out) :: x(n), y(n)
    integer :: i
    real :: h

    h = (abs(xEnd - x0)/n)
    x(1) = x0
    y(1) = y0

    do i = 1, n - 1
        x(i + 1) = x(i) + h
        y(i + 1) = y(i) + h * f(x(i), y(i))
    end do
    end subroutine euler

subroutine rk2(f, x0, y0, xEnd, n, x, y)
    interface
        real function f(x,y)
            real, intent(in) :: x, y
        end function f 
    end interface
    real, intent(in) :: x0, y0, xEnd
    integer, intent(in) :: n
    real, intent(out) :: x(n), y(n)
    integer :: i
    real :: h, k1, k2
    h = (abs(xEnd - x0) / n)
    x(1) = x0 
    y(1) = y0
    do i = 1, n-1 
        k1 = h * f(x(i), y(i))
        k2 = h * f(x(i) + h, y(i) + k1)
        x(i + 1) = x(i) + h
        y(i + 1) = y(i) + (k1 + k2) / 2.0
    end do
    end subroutine rk2

subroutine rk4(f, x0, y0, xEnd, n, x, y)
    interface
        real function f(x,y)
            real, intent(in) :: x, y
        end function f 
    end interface
    real, intent(in) :: x0, y0, xEnd
    integer, intent(in) :: n
    real, intent(out) :: x(n), y(n) 
    integer :: i 
    real :: h, k1, k2, k3, k4
    h = (abs(xEnd - x0) / n)
    x(1) = x0
    y(1) = y0
    do i = 1, n-1
        k1 = h * f(x(i), y(i))
        k2 = h * f(x(i) + h / 2.0, y(i) + k1 / 2.0)
        k3 = h * f(x(i) + h / 2.0, y(i) + k2 / 2.0)
        k4 = h * f(x(i) + h, y(i) + k3)
        x(i + 1) = x(i) + h
        y(i + 1) = y(i) + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    end do
    end subroutine rk4
