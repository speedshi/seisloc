subroutine combnum(nre,mcmdim,ncoe)
! This subroutine is used to calculate the combination number: C_nre^mcmdim
use paramod
implicit none

  integer(kind=INP),intent(in)  :: nre,mcmdim
  real(kind=RLP),intent(out)    :: ncoe
  integer(kind=8)               :: id,ztemp,mtemp

  ztemp=1
  do id=0,mcmdim-1
    ztemp=ztemp*(nre-id)
  enddo
  mtemp=1
  do id=1,mcmdim
    mtemp=mtemp*id
  enddo
  ncoe=REAL(ztemp/mtemp)

RETURN
end subroutine combnum