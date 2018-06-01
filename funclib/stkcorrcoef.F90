subroutine stkcorrcoef(ntwd,nre,data,stkval)
! This subroutine is used to calculate the stacked Pearson correlation coefficient of an input matrix.
! data(ntwd,nre) is the input matrix. RCC(nre,nre) is the Pearson correlation coefficient matrix.
! ntwd is the number of time points (observations), nre is the number of stations (variables).
! stkval is the final output stacking correlation coefficient.
use paramod
implicit none

  integer(kind=INP),intent(in)       :: ntwd,nre
  real(kind=RLP),intent(in)          :: data(ntwd,nre)
  real(kind=RLP),intent(out)         :: stkval
  integer(kind=INP)                  :: ii,jj
  real(kind=8)                       :: RCC(nre,nre),Xm(nre),Xrm(ntwd,nre),Xdev(nre,1)

  ! calculate the mean values of each variables
  Xm=sum(data,1)/ntwd

  ! remove the mean values for each varivables
  forall(ii=1:ntwd) Xrm(ii,:)=data(ii,:)-Xm

  ! calculate the deviations of each variables, note this is not standered deviations
  Xdev(:,1)=sqrt(sum(Xrm*Xrm,1))

  ! calculate the Pearson correlation coefficient matrix
  RCC=matmul(transpose(Xrm),Xrm)/matmul(Xdev,transpose(Xdev))

  ! calculate the stacking correlation coefficient
  forall (ii=1:nre,jj=1:nre,ii>=jj) RCC(ii,jj)=0
  stkval=SUM(ABS(RCC))

RETURN
end subroutine stkcorrcoef