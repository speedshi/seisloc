subroutine evprocesst(nsemax,nloadd,seiseloc,nsevt,nenpro,spaclim,timelim,nnse,nseall,eventsloc)
! This subroutine is used to re-process some seismic events.
! Because the dataset is devided into many data segment,
! the seismic evnets at the end of each time segment need to
! be reprocessed to assure they conform space and time limit.
! Input:------------------------------------------------------------------------------------------
! 'nsemax': row number of the matrix 'seissloc', also the maximum number of identified events for different data segments.
! 'nloadd': the number of data segments, also the number of loading times.
! 'seiseloc': information of identified seismic events for different data segment, Time-X-Y-Z-Migration_value, 5nsemax*nloadd.
! 'nsevt': the total number of identified seismic events for different data segment, nloadd*1.
! 'nenpro': the total number of identified seismic events that do not need re-processing for different data segment, nloadd*1.
! 'spaclim': the spatial limit to restrict identifying too close seismic events in space.
! 'timelim': the time limit to restrict identifying too close seismic events in time.
! 'nseall': total number of all potential seismic events
! Output:-----------------------------------------------------------------------------------------
! 'nnse': total number of the final identified seismic events.
! 'eventloc': information of the final identified seismic events, Time-X-Y-Z-Migration_value, 5nseall*nloadd.

use paramod
implicit none

  integer(kind=INP),intent(in)    :: nsemax,nloadd,nseall
  integer(kind=INP),intent(in)    :: nsevt(nloadd),nenpro(nloadd)
  real(kind=RLP),intent(inout)    :: seiseloc(5*nsemax,nloadd)
  real(kind=RLP),intent(in)       :: spaclim,timelim
  integer(kind=INP),intent(out)   :: nnse
  real(kind=RLP),intent(out)      :: eventsloc(5*nseall)
  integer(kind=INP)               :: iload,ie,ii
  real(kind=RLP)                  :: etime,ctime,pepos(3),dist,emv

  do iload=1,nloadd-1
    if (nsevt(iload+1)>0) then
      ! the next data segment has effective seismic events
      do ie=nenpro(iload)+1,nsevt(iload)
        etime=seiseloc(5*ie-4,iload)
        ctime=seiseloc(1,iload+1)
        pepos=seiseloc((5*ie-3):(5*ie-1),iload)
        emv=seiseloc(5*ie,iload)
        ii=1
        do while(ctime-etime<=timelim .AND. ii<=nsevt(iload+1))
          if (seiseloc(5*ii,iload+1)>0) then
            dist=SQRT(SUM((pepos-seiseloc((5*ii-3):(5*ii-1),iload+1))**2))
            if (dist<=spaclim) then
              if (emv>seiseloc(5*ii,iload+1)) then
                seiseloc(5*ii,iload+1)=-1
              else
                seiseloc(5*ie,iload)=-1
                EXIT
              endif
            endif
          endif
          ii=ii+1
          ctime=seiseloc(5*ii-4,iload+1)
        enddo
      enddo
    endif
  enddo

  ! form the final event location results
  nnse=0
  do iload=1,nloadd
    do ie=1,nsevt(iload)
      if (seiseloc(5*ie,iload)>0) then
        nnse=nnse+1
        eventsloc((5*nnse-4):(5*nnse))=seiseloc((5*ie-4):(5*ie),iload)
      endif
    enddo
  enddo

RETURN
end subroutine evprocesst