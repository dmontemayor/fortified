module functions
  use type_kinds
  implicit none
  private

  public:: logistic, dlogistic, bernoulli
  public:: gaussian, dgaussian, mexhat
  public:: dtanh

contains
!---------
  function bernoulli(x)
    real(double),intent(in)::x(:)
    real(double)::bernoulli(size(x)),r(size(x))
    bernoulli=1/(1+exp(-x))
    call random_number(r)
    where(r.LT.bernoulli)
       bernoulli=1._double
    elsewhere
       bernoulli=0._double
    end where
  end function bernoulli
!---------
  function logistic(x)
    real(double),intent(in)::x(:)
    real(double)::logistic(size(x))
    logistic=1/(1+exp(-x))
  end function logistic
!---------
  function dlogistic(x)
    real(double),intent(in)::x(:)
    real(double)::dlogistic(size(x))
    dlogistic=logistic(x)
    dlogistic=dlogistic*(1.0_double-dlogistic)
  end function dlogistic
!---------
  function dtanh(x)
    real(double),intent(in)::x(:)
    real(double)::dtanh(size(x))
    dtanh=tanh(x)
    dtanh=1.0_double-dtanh*dtanh
  end function dtanh
!---------
  function gaussian(x)
    real(double),intent(in)::x(:)
    real(double)::gaussian(size(x))
    gaussian=exp(-0.5_double*x**2)
  end function gaussian
!---------
  function dgaussian(x)
    real(double),intent(in)::x(:)
    real(double)::dgaussian(size(x))
    dgaussian=-x*gaussian(x)
  end function dgaussian
!---------
  function mexhat(x)
    real(double),intent(in)::x(:)
    real(double)::mexhat(size(x))
    mexhat=gaussian(x)
    mexhat=mexhat*(1-x*x)
  end function mexhat

end module functions

