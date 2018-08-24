subroutine eventidf_3ds(nsr,migvol_3d,soupos,vthrd,spaclim,nssot,event_sp,event_mv,npsit)
! This subroutine is used to identify seismic events in the 3D migration volume (3D space domain at the same origin time).
! The parameter 'spaclim' indicates the seismic events must separate at a distance larger than 'spaclim'.
! The parameter 'nssot' limits the maximum number of seismic events we can identify in the 3D volume.
! Input:---------------------------------------------------------------------------------
! 'nsr': totoal number of searching source points.
! 'migvol_3d': migration volume of all source points for a particular searching origin time, nsr*1.
! 'soupos': spatial positions of all source points, X-Y-Z (m), nsr*3.
! 'vthrd': used to set threshold value (veset) for pick seismic events in the migration volume, migration_value>=veset are viewed as potential seismic events.
!  If 0<vthrd<1, veset=vthrd; if vthrd<=0, veset=mean+3*std; if vthrd>=1, veset=mean+vthrd*std.
! 'spaclim': spatial limit to restrict identifying too close seismic events in space.
! 'nssot': the maximum number of seismic events we can identified in the migration volume.
! Output:--------------------------------------------------------------------------------
! 'npsit': the number of identified seismic events for this migration volume, npsit<=nssot.
! 'event_sp': locations (X-Y-Z) of potential seismic events, X-Y-Z, 3nssot*1.
! 'event_mv': migration values of potential seismic events, nssot*1.
use paramod
implicit none

  integer(kind=INP),intent(in)  :: nsr,nssot
  real(kind=RLP),intent(in)     :: vthrd,spaclim,migvol_3d(nsr),soupos(nsr,3)
  integer(kind=INP),intent(out) :: npsit
  real(kind=RLP),intent(out)    :: event_sp(3*nssot),event_mv(nssot)
  integer(kind=INP)             :: nseis,id,ii,pflag
  real(kind=RLP)                :: pepos(3),dist,mmean,mstd,veset
  integer(kind=INP),allocatable :: indx(:)
  real(kind=RLP),allocatable    :: mvseis(:),spseis(:,:)

  ! check and set the threshold value
  if (vthrd>0 .OR. vthrd<1) then
    ! user defined threshold value, between 0 and 1
    veset=vthrd
  elseif (vthrd<=0) then
    ! adaptive threshold, use predefined default outlier indication
    mmean=SUM(migvol_3d)/nsr
    mstd=SQRT(SUM((migvol_3d-mmean)**2)/(nsr-1.0))
    veset=mmean+3.0*mstd
  else
    ! adaptive threshold, use user defined outlier indication
    mmean=SUM(migvol_3d)/nsr
    mstd=SQRT(SUM((migvol_3d-mmean)**2)/(nsr-1.0))
    veset=mmean+vthrd*mstd
  endif

  npsit=0
  event_sp=0
  event_mv=0

  ! caculate the number of potential seismic events above the pre-set threshold.
  nseis=COUNT(migvol_3d>=veset)
  if (nseis==0) then
    ! migration values of all grid positions are below threshold, no seismic event.
    RETURN
  else
    if (nssot==1) then
      ! only one event is dentified for this migration volume
      allocate(indx(1))
      indx=MAXLOC(migvol_3d)
      npsit=1
      event_sp(1:3)=soupos(indx(1),:)
      event_mv(1)=migvol_3d(indx(1))
      deallocate(indx)
    else
      ! need to identify more than one event
      allocate(mvseis(nseis),spseis(nseis,3),indx(nseis))
      ii=0
      do id=1,nsr
        if (migvol_3d(id)>=veset) then
          ii=ii+1
          mvseis(ii)=migvol_3d(id)
          spseis(ii,:)=soupos(id,:)
        endif
      enddo

      ! sort the migration values in ascending order, note is in ascending order.
      CALL hpsort_eps_epw(nseis,mvseis,indx,1E-5)
    
      ! the point with maximum migration value certainly should be an seismic event
      event_sp(1:3)=spseis(indx(nseis),:)
      event_mv(1)=mvseis(nseis)
      npsit=1

      ! identify seismic events that satisfy the pre-set space limit
      do ii=nseis-1,1,-1
        pepos=spseis(indx(ii),:)
        pflag=1
        do id=1,npsit
          dist=SQRT(SUM((pepos-event_sp((3*id-2):(3*id)))**2))
          if (dist<=spaclim) then
            ! within the distance limit, not a qualified seismic event.
            pflag=0
            EXIT
          endif
        enddo

        if (pflag==1) then
          ! fullfill the distance limit, add a new seismic event
          npsit=npsit+1
          event_sp((3*npsit-2):(3*npsit))=pepos
          event_mv(npsit)=mvseis(ii)
          if (npsit==nssot) then
            ! reach the maximum number of events, return to the main program
            RETURN
          endif
        endif
      enddo
      deallocate(mvseis,spseis,indx)
    endif

  endif

RETURN
end subroutine eventidf_3ds


subroutine hpsort_eps_epw(n,ra,ind,eps)
  !---------------------------------------------------------------------
  ! sort an array ra(1:n) into ascending order using heapsort algorithm,
  ! and considering two elements being equal if their values differ
  ! for less than "eps".
  ! n is input, ra is replaced on output by its sorted rearrangement.
  ! create an index table (ind) by making an exchange in the index array
  ! whenever an exchange is made on the sorted data array (ra).
  ! in case of equal values in the data array (ra) the values in the
  ! index array (ind) are used to order the entries.
  ! if on input ind(1)  = 0 then indices are initialized in the routine,
  ! if on input ind(1) != 0 then indices are assumed to have been
  !                initialized before entering the routine and these
  !                indices are carried around during the sorting process
  !
  ! no work space needed !
  ! free us from machine-dependent sorting-routines !
  !
  ! adapted from Numerical Recipes pg. 329 (new edition)
  !
  use paramod
  implicit none
  !-input/output variables
  integer(kind=INP),intent(in)   :: n  
  real(kind=RLP),intent(in)      :: eps
  real(kind=RLP),intent(inout)   :: ra(n)
  integer(kind=INP),intent(out)  :: ind(n)
  !-local variables
  integer(kind=INP) :: i,ir,j,l,iind
  real(kind=RLP)    :: rra
!
  ! initialize index array
  forall(i=1:n)  ind(i)=i
  
  ! nothing to order
  IF (n.lt.2) return  
  ! initialize indices for hiring and retirement-promotion phase
  l = n / 2 + 1  

  ir = n  

  sorting: do 
  
    ! still in hiring phase
    IF ( l .gt. 1 ) then  
       l    = l - 1  
       rra  = ra (l)  
       iind = ind (l)  
       ! in retirement-promotion phase.
    ELSE  
       ! clear a space at the end of the array
       rra  = ra (ir)  
       !
       iind = ind (ir)  
       ! retire the top of the heap into it
       ra (ir) = ra (1)  
       !
       ind (ir) = ind (1)  
       ! decrease the size of the corporation
       ir = ir - 1  
       ! done with the last promotion
       IF ( ir .eq. 1 ) then  
          ! the least competent worker at all !
          ra (1)  = rra  
          !
          ind (1) = iind  
          exit sorting  
       ENDIF
    ENDIF
    ! wheter in hiring or promotion phase, we
    i = l  
    ! set up to place rra in its proper level
    j = l + l  
    !
    DO while ( j .le. ir )  
       IF ( j .lt. ir ) then  
          ! compare to better underling
          IF ( hslt( ra (j),  ra (j + 1) ) ) then  
             j = j + 1  
          !else if ( .not. hslt( ra (j+1),  ra (j) ) ) then
             ! this means ra(j) == ra(j+1) within tolerance
           !  if (ind (j) .lt.ind (j + 1) ) j = j + 1
          ENDIF
       ENDIF
       ! demote rra
       IF ( hslt( rra, ra (j) ) ) then  
          ra (i) = ra (j)  
          ind (i) = ind (j)  
          i = j  
          j = j + j  
       !else if ( .not. hslt ( ra(j) , rra ) ) then
          !this means rra == ra(j) within tolerance
          ! demote rra
         ! if (iind.lt.ind (j) ) then
         !    ra (i) = ra (j)
         !    ind (i) = ind (j)
         !    i = j
         !    j = j + j
         ! else
             ! set j to terminate do-while loop
         !    j = ir + 1
         ! endif
          ! this is the right place for rra
       ELSE
          ! set j to terminate do-while loop
          j = ir + 1  
       ENDIF
    ENDDO
    ra (i) = rra  
    ind (i) = iind  

  END DO sorting    
contains 

  !  internal function 
  !  compare two real number and return the result

  logical function hslt(a,b)
    REAL(kind=RLP) :: a,b
    IF( abs(a-b) <  eps ) then
      hslt = .false.
    ELSE
      hslt = ( a < b )
    end if
  end function hslt

  !
end subroutine hpsort_eps_epw