subroutine eventidf_stm(nssot,nntld,event_sp,event_mv,npsit,timelim,spaclim,dt0,pst0,nnpse,seisinfo,nsevt,nenpro)
! This subroutine is used to identify the seismic events for continuous searching origin times.
! The parameter 'timelim' indicates the close located seismic events must have a time interval larger than 'timelim'.
! Input:------------------------------------------------------------------------------------
! 'nssot': the maximum number of seismic events we can identified for one searching origin time.
! 'nntld': the number of searching origin times.
! 'event_sp': locations (X-Y-Z) of potential seismic events, X-Y-Z, 3nssot*nntld.
! 'event_mv': migration values of potential seismic events, nssot*nntld.
! 'npsit': the number of actual identified potential seismic events for different searching origin times, nntld*1.
! 'timelim': the time limit to restrict identifying too close seismic events in time.
! 'spaclim': the spatial limit to restrict identifying too close seismic events in space.
! 'dt0': the time interval of searching origin times.
! 'pst0': searching origin times, nntld*1.
! 'nnpse': the total number of the input potential seismic events.
! Output:-----------------------------------------------------------------------------------
! 'seisinfo': information of identified seismic events, Time-X-Y-Z-Migration_value, 5nnpse*1.
! 'nsevt': the total number of identified seismic events.
! 'nenpro': the total number of identified seismic events that do not need re-processing using results of next data segment.
use paramod
implicit none

  integer(kind=INP),intent(in)    :: nssot,nntld,nnpse
  integer(kind=INP),intent(in)    :: npsit(nntld)
  real(kind=RLP),intent(in)       :: pst0(nntld),timelim,spaclim,dt0
  real(kind=RLP),intent(inout)    :: event_sp(3*nssot,nntld),event_mv(nssot,nntld)
  real(kind=RLP),intent(out)      :: seisinfo(5*nnpse)
  integer(KIND=INP),intent(out)   :: nsevt,nenpro
  integer(KIND=INP)               :: nnt,it,id,in,inid,iflag
  real(KIND=RLP)                  :: pepos(3),psmv,dist

  nnt=INT(timelim/dt0)

  if (nnt>nntld) then
    ! process some extreme conditions
    nnt=nntld
  endif

  if (nnt==0) then
    ! all the input seismic events fulfill the time limit requirement
    nsevt=nnpse
    nenpro=nsevt
    forall (it=1:nntld)
      forall (id=1:npsit(it))
        seisinfo((SUM(npsit(1:(it-1)))+id)*5-4)=pst0(it) ! origin time
        seisinfo(((SUM(npsit(1:(it-1)))+id)*5-3):((SUM(npsit(1:(it-1)))+id)*5-1))=event_sp((id*3-2):(id*3),it) ! X-Y-Z position
        seisinfo((SUM(npsit(1:(it-1)))+id)*5)=event_mv(id,it) ! migration value
      end forall
    end forall
    RETURN
  else
    nsevt=0
    do it=1,nntld-nnt
      do id=1,npsit(it)
        if (event_mv(id,it)>0) then
          ! effective seismic events
          iflag=1
          pepos=event_sp((id*3-2):(id*3),it) ! position
          psmv=event_mv(id,it) ! migration value
          INNERT: do in=it+1,it+nnt
                    do inid=1,npsit(in)
                      if (event_mv(inid,in)>0) then
                        ! the compared seismic event is effective
                        dist=SQRT(SUM((pepos-event_sp((inid*3-2):(inid*3),in))**2))
                        if (dist<=spaclim) then
                        ! two events are close and within the space limit
                          if (psmv>event_mv(inid,in)) then
                            event_mv(inid,in)=-1
                          else
                            iflag=0
                            EXIT INNERT
                          endif
                        endif
                      endif
                    enddo
          enddo INNERT
        else
          iflag=0
        endif
        if (iflag==1) then
          ! find and add a seismic event
          nsevt=nsevt+1 ! number of seismic events
          seisinfo(nsevt*5-4)=pst0(it) ! origin time
          seisinfo((nsevt*5-3):(nsevt*5-1))=pepos ! position of seismic event
          seisinfo(nsevt*5)=psmv ! migration value
        endif
      enddo
    enddo

    nenpro=nsevt

    ! the last few searching origin times
    do it=nntld-nnt+1,nntld
      do id=1,npsit(it)
        if (event_mv(id,it)>0) then
          ! effective seismic events
          iflag=1
          pepos=event_sp((id*3-2):(id*3),it) ! position
          psmv=event_mv(id,it) ! migration value
          INNERT2: do in=it+1,nntld
                    do inid=1,npsit(in)
                      if (event_mv(inid,in)>0) then
                        ! the compared seismic event is effective
                        dist=SQRT(SUM((pepos-event_sp((inid*3-2):(inid*3),in))**2))
                        if (dist<=spaclim) then
                        ! two events are close and within the space limit
                          if (psmv>event_mv(inid,in)) then
                            event_mv(inid,in)=-1
                          else
                            iflag=0
                            EXIT INNERT2
                          endif
                        endif
                      endif
                    enddo
          enddo INNERT2
        else
          iflag=0
        endif
        if (iflag==1) then
          ! find and add a seismic event
          nsevt=nsevt+1 ! number of seismic events
          seisinfo(nsevt*5-4)=pst0(it) ! origin time
          seisinfo((nsevt*5-3):(nsevt*5-1))=pepos ! position of seismic event
          seisinfo(nsevt*5)=psmv ! migration value
        endif
      enddo
    enddo

  endif

RETURN
end subroutine eventidf_stm


subroutine eventidf_stm2(nssot,nntld,event_sp,event_mv,npsit,timelim,spaclim,dt0,pst0,nnpse,seisinfo,nsevt,nenpro)
! This subroutine is used to identify the seismic events for continuous searching origin times.
! Different with 'eventidf_stm', this subroutine first sort the events according to the migration value and then pick the suitable events.
! The parameter 'timelim' indicates the close located seismic events must have a time interval larger than 'timelim'.
! Input:------------------------------------------------------------------------------------
! 'nssot': the maximum number of seismic events we can identified for one searching origin time.
! 'nntld': the number of searching origin times.
! 'event_sp': locations (X-Y-Z) of potential seismic events, X-Y-Z, 3nssot*nntld, here we receive it as an one-dimension vector, (3*nssot*nntld)*1.
! 'event_mv': migration values of potential seismic events, nssot*nntld, here we receive it as an one-dimension vector, (nssot*nntld)*1.
! 'npsit': the number of actual identified potential seismic events for different searching origin times, nntld*1.
! 'timelim': the time limit to restrict identifying too close seismic events in time.
! 'spaclim': the spatial limit to restrict identifying too close seismic events in space.
! 'dt0': the time interval of searching origin times.
! 'pst0': searching origin times, nntld*1.
! 'nnpse': the total number of the input potential seismic events.
! Output:-----------------------------------------------------------------------------------
! 'seisinfo': information of identified seismic events, Time-X-Y-Z-Migration_value, 5nnpse*1.
! 'nsevt': the total number of identified seismic events.
! 'nenpro': the total number of identified seismic events that do not need re-processing using results of next data segment.
use paramod
implicit none

  integer(kind=INP),intent(in)    :: nssot,nntld,nnpse
  integer(kind=INP),intent(in)    :: npsit(nntld)
  real(kind=RLP),intent(in)       :: pst0(nntld),timelim,spaclim,dt0
  real(kind=RLP),intent(inout)    :: event_sp(3*nssot*nntld),event_mv(nssot*nntld)
  real(kind=RLP),intent(out)      :: seisinfo(5*nnpse)
  integer(KIND=INP),intent(out)   :: nsevt,nenpro
  integer(KIND=INP)               :: nne,indx(nssot*nntld),it,id,iflag,neffc
  real(KIND=RLP)                  :: event_tm(nssot*nntld),temp_tm(nnpse),temp_sp(3*nnpse),temp_mv(nnpse),petm,pepos(3),dist,tnit

  nne=nssot*nntld
  forall (it=1:nntld) event_tm(((it-1)*nssot+1):(it*nssot))=pst0(it)

  temp_tm=pst0(nntld)+9999

  ! sort the migration values in ascending order for all origin times, note is in ascending order.
  CALL hpsort_eps_epw(nne,event_mv,indx,1E-5)

  neffc=COUNT(event_mv>0) ! number of effective potential seismic events (migration_value > 0)
  
  nsevt=1
  temp_tm(1)=event_tm(indx(nne))
  temp_sp(1:3)=event_sp((3*indx(nne)-2):(3*indx(nne)))
  temp_mv(1)=event_mv(nne)
  do it=nne-1,nne-neffc+1,-1
    iflag=1
    petm=event_tm(indx(it))
    pepos=event_sp((3*indx(it)-2):(3*indx(it)))

    do id=1,nsevt
      dist=SQRT(SUM((pepos-temp_sp((3*id-2):(3*id)))**2))
      tnit=ABS(petm-temp_tm(id))
      if (dist<=spaclim .AND. tnit<= timelim) then
        ! within the space and time limit
        iflag=0
        EXIT
      endif
    enddo
    
    if (iflag==1) then
      nsevt=nsevt+1
      temp_tm(nsevt)=petm
      temp_sp((3*nsevt-2):(3*nsevt))=pepos
      temp_mv(nsevt)=event_mv(it)
    endif
  enddo

  ! sort the identified events according to their origin times
  CALL hpsort_eps_epw(nsevt,temp_tm,indx,1E-5)

  nenpro=COUNT(temp_tm<=(pst0(nntld)-timelim))
  seisinfo=0

  forall (it=1:nsevt)
    seisinfo(5*it-4)=temp_tm(it)
    seisinfo((5*it-3):(5*it-1))=temp_sp((3*indx(it)-2):(3*indx(it)))
    seisinfo(5*it)=temp_mv(indx(it))
  end forall

RETURN
end subroutine eventidf_stm2