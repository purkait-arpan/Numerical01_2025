program main
  real, external :: bisection
  print *, "sqrt of 2 is (calculated) :", bisection(f,1.0,2.0,0.0001)
  !print *, "sqrt of 2 is (calculated) :", newton_raphson(f,df,2.0,0.0001)
  print *, "Actual value :", sqrt(2.0) 

contains
    real function f(x)
      real, intent(in) :: x
      f = (x**2 - 2)
    end function f
    real function df(x)
      real, intent(in) :: x
      df = 2.0*x
    end function df
end program main
