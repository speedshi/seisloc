subroutine stkcorrcoef(ntwd,nre,data,ncoe,stkval)
! This subroutine is used to calculate the stacked Pearson correlation coefficient of an input matrix.
! data(ntwd,nre) is the input matrix. RCC(nre,nre) is the Pearson correlation coefficient matrix.
! ntwd is the number of time points (observations), nre is the number of stations (variables).
! ncoe is the number of caculated correlation coefficients of the corresponding station pairs or groups.
! stkval is the final output normalized stacking correlation coefficient.
use paramod
implicit none

  integer(kind=INP),intent(in)       :: ntwd,nre
  real(kind=RLP),intent(in)          :: ncoe,data(ntwd,nre)
  real(kind=RLP),intent(out)         :: stkval
  integer(kind=INP)                  :: ii,jj,nps
  real(kind=RLP)                     :: RCC(nre,nre),Xm(nre),Xrm(ntwd,nre),Xdev(nre,1)

  ! calculate the mean values of each variables
  Xm=sum(data,1)/ntwd

  ! remove the mean values for each varivables
  forall(ii=1:ntwd) Xrm(ii,:)=data(ii,:)-Xm

  ! calculate the deviations of each variables, note this is not standered deviations
  Xdev(:,1)=sqrt(sum(Xrm*Xrm,1))

  ! calculate the Pearson correlation coefficient matrix
  RCC=matmul(transpose(Xrm),Xrm)/matmul(Xdev,transpose(Xdev))

  ! set the low triangle part to be 0
  forall (ii=1:nre,jj=1:nre,ii>=jj) RCC(ii,jj)=0

  ! set the infinite values to NAN; all values should be no larger than 1, otherwise are infinite values
  where(ABS(RCC)>1) RCC=NAN

  ! count the number of NAN values in the correlation matrix
  nps=COUNT(ISNAN(RCC))

  ! set the NAN to 0
  where(ISNAN(RCC)) RCC=0

  ! calculate the normalized stacking correlation coefficient
  stkval=SUM(ABS(RCC))/REAL(ncoe-nps)

RETURN
end subroutine stkcorrcoef


subroutine stkcorrcoef1(ntwd,nre,data,RCC,Xm,Xdev,ncoe,iii,jjj,nps,stkval)
! This subroutine is used to calculate the stacked Pearson correlation coefficient of an input matrix.
! data(ntwd,nre) is the input matrix. RCC(nre,nre) is the Pearson correlation coefficient matrix.
! ntwd is the number of time points (observations), nre is the number of stations (variables).
! ncoe is the number of caculated correlation coefficients of the corresponding station pairs or groups.
! stkval is the final output normalized stacking correlation coefficient.
use paramod
implicit none

  integer(kind=INP),intent(in)       :: ntwd,nre
  real(kind=RLP),intent(in)          :: ncoe
  real(kind=RLP),intent(inout)       :: data(ntwd,nre),RCC(nre,nre),Xm(nre),Xdev(nre)
  real(kind=RLP),intent(out)         :: stkval
  integer(kind=INP),intent(inout)    :: iii,jjj,nps

  ! calculate the mean values of each variables
  Xm=sum(data,1)/ntwd

  ! remove the mean values for each varivables
  forall(iii=1:ntwd) data(iii,:)=data(iii,:)-Xm

  ! calculate the deviations of each variables, note this is not standered deviations
  Xdev=sqrt(sum(data*data,1))

  ! calculate the Pearson correlation coefficient matrix, only for up triangle part
  RCC=0
  forall (iii=1:nre,jjj=1:nre,iii<jjj) RCC(iii,jjj)=sum(data(:,iii)*data(:,jjj))/(Xdev(iii)*Xdev(jjj))

  ! set the infinite values to NAN; all values should be no larger than 1, otherwise are infinite values
  where(ABS(RCC)>1) RCC=NAN

  ! count the number of NAN values in the correlation matrix
  nps=COUNT(ISNAN(RCC))

  ! set the NAN to 0
  where(ISNAN(RCC)) RCC=0

  ! calculate the normalized stacking correlation coefficient
  stkval=SUM(ABS(RCC))/REAL(ncoe-nps)

RETURN
end subroutine stkcorrcoef1
