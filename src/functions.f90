module functions
  use type_kinds
  use math
  implicit none
  private

  public:: softplus,logistic,dlogistic
  public:: bernoulli,poisson
  public:: gaussian,mexhat
  public:: softmax,dsoftmax
  public:: dtanh

contains
  !---------
  function softplus(x)
    real(double),intent(in)::x(:)
    real(double)::softplus(size(x))
    softplus=log(1+exp(x))
  end function softplus
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
    dlogistic=dlogistic*(1.0_double-dlogistic)
  end function dlogistic
  !---------
  function bernoulli(x)
    real(double),intent(in)::x(:)
    real(double)::bernoulli(size(x)),r(size(x))
    bernoulli=logistic(x)
    where(r.LT.bernoulli)
       bernoulli=1._double
    elsewhere
       bernoulli=0._double
    end where
  end function bernoulli
  !---------
  function poisson(x)
    real(double),intent(in)::x(:)
    real(double)::poisson(size(x)),r(size(x))
    poisson=exp(-x)
    call random_number(r)
    where(r.LT.poisson)
       poisson=1._double
    elsewhere
       poisson=0._double
    end where
  end function poisson
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
    mexhat=gaussian(x)*(1-x*x)
  end function mexhat
  !---------
  function softmax(x)
    real(double),intent(in)::x(:)
    real(double)::softmax(size(x))
    integer::i
    softmax=exp(x)
    softmax=softmax/sum(softmax)
  end function softmax
  !---------
  function dsoftmax(x)
    real(double),intent(in)::x(:)
    real(double)::dsoftmax(size(x),size(x)),s(size(x))
    integer::i,j
    s=softmax(x)
    do i=1,size(x)
       do j=i,size(x)
          dsoftmax(i,j)=s(i)*(krondelta(i,j)-s(j))
       end do
    end do
  end function dsoftmax
  !---------
  function dtanh(x)
    real(double),intent(in)::x(:)
    real(double)::dtanh(size(x))
    integer::i
    dtanh=tanh(x)
    dtanh=1.0_double-dtanh*dtanh
  end function dtanh

end module functions

