subroutine corrcoefm(n,m,X,R)
! This subroutine is used to calculate the Pearson correlation coefficient matrix of an input matrix.
! X(n,m) is the input matrix. R(m,m) is the output Pearson correlation coefficient matrix.
! n is the number of observations, m is the number of variables.
use paramod
implicit none

  integer(kind=INP),intent(in)       :: n,m
  real(kind=RLP),intent(in)          :: X(n,m)
  real(kind=RLP),intent(out)         :: R(m,m)
  integer(kind=INP)                  :: i
  real(kind=RLP)                     :: Xm(m),Xrm(n,m),Xdev(m,1)

  ! calculate the mean values of each variables
  Xm=sum(X,1)/n

  ! remove the mean values for each varivables
  forall(i=1:n) Xrm(i,:)=X(i,:)-Xm

  ! calculate the deviations of each variables, note this is not standered deviations
  Xdev(:,1)=sqrt(sum(Xrm*Xrm,1))

  !calculate the Pearson correlation coefficient matrix
  R=matmul(transpose(Xrm),Xrm)/matmul(Xdev,transpose(Xdev))

end subroutine corrcoefm