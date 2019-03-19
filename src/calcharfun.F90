subroutine calcharfun(wfdata,dimloadd,cfuntp)
! This subroutine is used to calculate the characteristic function of the input data.
! Note the input data 'wfdata' is a matrix, ntsload*nre, row->time samples, colume->different stations.
use paramod
implicit none
  integer(kind=INP),intent(in) :: cfuntp,dimloadd(2)
  real(kind=RLP),intent(out)   :: wfdata(dimloadd(1),dimloadd(2))

  if (cfuntp==1) then
    ! envelope
    CALL gdtmansg(dimloadd(1),dimloadd(2),wfdata)
  elseif (cfuntp==2) then
    ! absolute value
    wfdata=ABS(wfdata)
  elseif (cfuntp==3) then
    ! non-negative value
    where(wfdata<0) wfdata=0
  elseif (cfuntp==4) then
    ! square value
    wfdata=wfdata**2
  else
    ! wrong input
    write(*,*) 'Wrong input parameter for cfuntp! Program stop!'
    STOP
  endif

RETURN
end subroutine calcharfun