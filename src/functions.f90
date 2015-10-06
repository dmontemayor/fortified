module functions
  use type_kinds
  implicit none
  private

  public:: linear, logistic, gaussian,mexhat

contains
  function linear(x)
    real(double),intent(in)::x(:)
    real(double)::linear(size(x))
    linear=x
  end function linear
!---------
  function logistic(x)
    real(double),intent(in)::x(:)
    real(double)::logistic(size(x))
    logistic=1/(1+exp(-x))
  end function logistic
!---------
  function gaussian(x)
    real(double),intent(in)::x(:)
    real(double)::gaussian(size(x))
    gaussian=exp(-0.5_double*x**2)
  end function gaussian
!---------
  function mexhat(x)
    real(double),intent(in)::x(:)
    real(double)::mexhat(size(x))
    mexhat=gaussian(x)
    mexhat=mexhat*(1-x*x)
  end function mexhat

end module functions

