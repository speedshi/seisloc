program migrationloc
! This program is used to perform waveform migration to locate the source position and origin time.
! All parameters use SI unit, i.e.: length in meter, time in second.
use paramod
implicit none

  integer(kind=INP) :: migtp,phasetp,cfuntp,nre,nsr,ntdata,ntpwd,ntswd,ntwd,nssot,nt,ntsload,ntmod,nloadd,nst0
  integer(kind=INP) :: iload,idxt,id,it,ir,ntoverlap,mcmdim,migplan,nntld,nnpse,nsemax,nseall,nnse
  integer(kind=INP) :: time_array_s(8), time_array_e(8)
  character(100)    :: dfname,opname,infname,outfname
  real(kind=RLP)    :: dt,tdatal,tpwind,tswind,dt0,vthrd,ncoe,stkz,stkn,stke,spaclim,timelim,twallc
  integer(kind=INP),allocatable             :: tvtn(:),tvpn(:),tvsn(:),s0idsg(:),npsit(:),nsevt(:),nenpro(:)
  real(kind=RLP),allocatable,dimension(:,:) :: soupos,travelt,travelp,travels,zsdatain,nsdatain,esdatain,zwfdata,nwfdata,ewfdata,wfex,zwfex,nwfex,ewfex,event_sp,event_mv,seiseloc
  real(kind=RLP),allocatable,dimension(:)   :: st0,btime,migvol_3d,pst0,seisinfo,eventsloc

  ! start timing. Wall-clock time, program must run in the same month
  call date_and_time(values=time_array_s)

  ! set or create the folder name for input and output
  infname="./data/"
  outfname="./results/"
  CALL system("mkdir "//TRIM(outfname))

  ! read in migration parameters
  open(unit=100,file=TRIM(infname)//'migpara.dat')
  read(100,*) migtp
  read(100,*) phasetp
  read(100,*) cfuntp
  read(100,*) nre
  read(100,*) nsr
  read(100,*) dfname
  read(100,*) dt
  read(100,*) tdatal
  read(100,*) tpwind
  read(100,*) tswind
  read(100,*) dt0
  read(100,*) vthrd
  read(100,*) mcmdim
  read(100,*) spaclim
  read(100,*) timelim
  read(100,*) nssot
  close(unit=100)
 
  ! determine the number of time points for the whole seismic data
  ntdata=NINT(tdatal/dt)+1
  ! determine the number of time points of the P- and S-phase time window used in the migration process. 
  ! For MCM is the coherency analysis time window length. 
  ! For traditional migration, stacking is performed using data starting from the calculated arrival times and within this time window (not centred at the arrival times).
  ! ntwd describes how many time points are used in the migration or coherency calculation.
  ntpwd=NINT(tpwind/dt)+1
  ntswd=NINT(tswind/dt)+1

  ! load the potential source position file
  allocate(soupos(nsr,3))
  open(unit=10,file=TRIM(infname)//'soupos.dat',form='unformatted',status='old',access='stream',action='read')
  read(10) soupos
  close(unit=10)

  ! allocate memory, load the traveltime table and decide the migration plan
  if (phasetp==0) then
    ! only P-phase is used in the migration
    ntwd=ntpwd
    allocate(travelt(nsr,nre),tvtn(nre),wfex(ntwd,nre))
    ! load traveltime table for P-phase
    open(unit=11,file=TRIM(infname)//'travelp.dat',form='unformatted',status='old',access='stream',action='read')
    read(11) travelt
    close(unit=11)
    ! determine the overlap time points for loading the data
    ntoverlap=CEILING(MAXVAL(travelt)/dt)+ntwd+1
    if (migtp==0) then
      ! using P-phase and MCM
      migplan=1
      ! calculate the number of effective station groups used in the migration
      CALL combnum(nre,mcmdim,ncoe)
    else
      ! using P-phase and conventional DSI
      migplan=2
      ! calculate the number of effective stations used in the migration
      ncoe=REAL(nre)
    endif
  elseif (phasetp==1) then
    ! only S-phase is used in the migration
    ntwd=ntswd
    allocate(travelt(nsr,nre),tvtn(nre),wfex(ntwd,nre))
    ! load traveltime table for S-phase
    open(unit=12,file=TRIM(infname)//'travels.dat',form='unformatted',status='old',access='stream',action='read')
    read(12) travelt
    close(unit=12)
    ! determine the overlap time points for loading the data
    ntoverlap=CEILING(MAXVAL(travelt)/dt)+ntwd+1
    if (migtp==0) then
      ! using S-phase and MCM
      migplan=1
      ! calculate the number of effective station groups used in the migration
      CALL combnum(nre,mcmdim,ncoe)
    else
      ! using S-phase and conventional DSI
      migplan=2
      ! calculate the number of effective stations used in the migration
      ncoe=REAL(nre)
    endif
  else
    ! both P- and S-phases are used in the migration
    allocate(travelp(nsr,nre),travels(nsr,nre),tvpn(nre),tvsn(nre),zwfex(ntpwd,nre),nwfex(ntswd,nre),ewfex(ntswd,nre))
    ! load the travetime tables of P- and S-phases
    open(unit=11,file=TRIM(infname)//'travelp.dat',form='unformatted',status='old',access='stream',action='read')
    read(11) travelp
    close(unit=11)
    open(unit=12,file=TRIM(infname)//'travels.dat',form='unformatted',status='old',access='stream',action='read')
    read(12) travels
    close(unit=12)
    ntoverlap=MAX(CEILING(MAXVAL(travelp)/dt)+ntpwd,CEILING(MAXVAL(travels)/dt)+ntswd)+1
    if (migtp==0) then
      ! using both P- and S-phases and MCM, three component data
      migplan=3
      ! calculate the number of effective station groups used in the migration
      CALL combnum(nre,mcmdim,ncoe)
    else
      ! using both P- and S-phases and conventional DSI
      migplan=4
      ! calculate the number of effective stations used in the migration
      ncoe=REAL(nre)
    endif
  endif
  
  ! check the maximum number of potential seismic events for a single origin time.
  if (nssot<1) then
    nssot=100
  endif
  ! set the number of time points for loading in the data per time. For very large dataset, we need to load and process segment data many times.
  ! For 250 station, single precision, nt=10^6, the memory cost of the data is about 1G.
  nt=1000000
  ! determine the number of time points that really loaded each time, the loading data is overlaped to assure continuously origin time searching
  ntsload=nt+ntoverlap
  ! determine the number of searching origin times for the whole dataset
  nst0=INT((tdatal-ntoverlap*dt)/dt0)+1

  ! allocate matrixs
  allocate(st0(nst0),zwfdata(ntsload,nre),nwfdata(ntsload,nre),ewfdata(ntsload,nre),migvol_3d(nsr))

  ! determine the searching origin times
  forall (idxt=1:nst0) st0(idxt)=(idxt-1)*dt0

  ! determine the times for loading the whole dataset
  ntmod=MOD(ntdata,nt)
  if (ntmod>ntoverlap) then
    nloadd=INT(REAL(ntdata)/REAL(nt))+1
  elseif (ntmod>0) then
    nloadd=INT(REAL(ntdata)/REAL(nt))
  else
    nloadd=NINT(REAL(ntdata)/REAL(nt))
  endif

  ! calculate the starting time of each data segment (btime)
  ! calculate the starting index of origin times for each data segment (s0idsg)
  allocate(s0idsg(nloadd+1),btime(nloadd),nsevt(nloadd),nenpro(nloadd))
  btime(1)=0
  s0idsg(1)=1
  s0idsg(nloadd+1)=nst0+1
  if (nloadd>1) then
    forall (iload=2:nloadd)
      btime(iload)=(iload-1)*nt*dt
      s0idsg(iload)=COUNT(st0<btime(iload))+1
    endforall
  endif

  ! read input data, the size of the whole dataset is nre*ntdata. Note the storage sequence of the data should first be 'station', then the 'time'.
  ! steam I/O
  open(unit=13,file=TRIM(infname)//TRIM(dfname)//'.z',form='unformatted',status='old',access='stream',action='read')
  open(unit=14,file=TRIM(infname)//TRIM(dfname)//'.n',form='unformatted',status='old',access='stream',action='read')
  open(unit=15,file=TRIM(infname)//TRIM(dfname)//'.e',form='unformatted',status='old',access='stream',action='read')

  ! staring migration and location according to different migration plan
    ! migration plan: single phase and MCM
    do iload=1,nloadd
      if (iload<nloadd) then
        allocate(zsdatain(nre,ntsload),nsdatain(nre,ntsload),esdatain(nre,ntsload))
      else
        ! different memory allocation for the 'iload=nloadd' case
        if (ntmod>0) then
          deallocate(zwfdata,nwfdata,ewfdata)
          allocate(zsdatain(nre,ntmod),nsdatain(nre,ntmod),esdatain(nre,ntmod),zwfdata(ntmod,nre),nwfdata(ntmod,nre),ewfdata(ntmod,nre))
        else
          deallocate(zwfdata,nwfdata,ewfdata)
          allocate(zsdatain(nre,nt),nsdatain(nre,nt),esdatain(nre,nt),zwfdata(nt,nre),nwfdata(nt,nre),ewfdata(nt,nre))
        endif
      endif
      ! read the current data segment, note we set the position of the loading data segment
      read(13,POS=RLP*(iload-1)*nt+1) zsdatain
      read(14,POS=RLP*(iload-1)*nt+1) nsdatain
      read(15,POS=RLP*(iload-1)*nt+1) esdatain
      ! transpose the matrix, change the storage sequence to accelerate the calculation
      zwfdata=TRANSPOSE(zsdatain)
      nwfdata=TRANSPOSE(nsdatain)
      ewfdata=TRANSPOSE(esdatain)
      deallocate(zsdatain,nsdatain,esdatain)

      ! allocate memory
      nntld=s0idsg(iload+1)-s0idsg(iload) ! number of searching origin times for the current loading dataset
      allocate(pst0(nntld),npsit(nntld),event_sp(3*nssot,nntld),event_mv(nssot,nntld))
      ! the searching origin times for the current loading dataset
      pst0=st0(s0idsg(iload):(s0idsg(iload+1)-1))

      ! migration and imaging
      ! both P- and S-phases
      !$OMP PARALLEL DO PRIVATE(it,id,ir,tvpn,tvsn,zwfex,nwfex,ewfex,stkz,stkn,stke,migvol_3d) SCHEDULE(DYNAMIC)
      do it=1,nntld
        do id=1,nsr
          ! calculate the starting time points for migration
          tvpn=NINT((travelp(id,:)+pst0(it)-btime(iload))/dt)+1
          tvsn=NINT((travels(id,:)+pst0(it)-btime(iload))/dt)+1
          ! extract the waveforms
          forall (ir=1:nre)
            zwfex(:,ir)=zwfdata(tvpn(ir):(tvpn(ir)+ntpwd-1),ir)
            nwfex(:,ir)=nwfdata(tvsn(ir):(tvsn(ir)+ntswd-1),ir)
            ewfex(:,ir)=ewfdata(tvsn(ir):(tvsn(ir)+ntswd-1),ir)
          endforall
          ! calculate the stacking value of this imaging point
          CALL stkcorrcoef(ntpwd,nre,zwfex,ncoe,stkz)
          CALL stkcorrcoef(ntswd,nre,nwfex,ncoe,stkn)
          CALL stkcorrcoef(ntswd,nre,ewfex,ncoe,stke)
          migvol_3d(id)=0.6*stkz+0.2*(stkn+stke)
        enddo
        ! find potential source locations in 3D space domain for a particular origin time
        CALL eventidf_3ds(nsr,migvol_3d,soupos,vthrd,spaclim,nssot,event_sp(1,it),event_mv(1,it),npsit(it))
      enddo
      !$OMP END PARALLEL DO

      ! calculate total number of potential seismic events for all searching origin times
      nnpse=SUM(npsit)
      if (nnpse>0) then
        ! allocate memory for saving information of real seismic events
        allocate(seisinfo(5*nnpse))
        ! find locations of seismic events that fullfill the time limit
        CALL eventidf_stm(nssot,nntld,event_sp,event_mv,npsit,timelim,spaclim,dt0,pst0,nnpse,seisinfo,nsevt(iload),nenpro(iload))
        write(opname,*) iload
        open(unit=66,file=TRIM(outfname)//TRIM(ADJUSTL(opname)),form='unformatted',status='replace',access='stream',action='write')
        write(66) seisinfo(1:(5*nsevt(iload)))
        close(unit=66)
        deallocate(seisinfo)
      else
        ! no seismic events in this loading time period, no output
        nsevt(iload)=0
        nenpro(iload)=0
      endif
      deallocate(pst0,npsit,event_sp,event_mv)
    enddo

  close(unit=13)
  close(unit=14)
  close(unit=15)

  ! deallocate memory
  deallocate(zwfdata,nwfdata,ewfdata,soupos,migvol_3d,st0,s0idsg,btime)
  if (phasetp==0 .OR. phasetp==1) then
    ! only single phase is used in the migration
    deallocate(travelt,tvtn,wfex)
  else
    ! both P- and S-phases are used in the migration
    deallocate(travelp,travels,tvpn,tvsn,zwfex,nwfex,ewfex)
  endif

  ! load results of all data segments and combine them to form the final results
  nsemax=MAXVAL(nsevt)
  allocate(seiseloc(5*nsemax,nloadd))
  seiseloc=0
  do iload=1,nloadd
    if (nsevt(iload)>0) then
      write(opname,*) iload
      open(unit=666,file=TRIM(outfname)//TRIM(ADJUSTL(opname)),form='unformatted',status='old',access='stream',action='read')
      read(666) seiseloc(1:(5*nsevt(iload)),iload)
      close(unit=666)
    endif
  enddo

  !re-process some seismic events
  nseall=SUM(nsevt)
  allocate(eventsloc(5*nseall))
  if (nloadd>1) then
    ! more than one data segment
    CALL evprocesst(nsemax,nloadd,seiseloc,nsevt,nenpro,spaclim,timelim,nnse,nseall,eventsloc)
  else
    ! only one data segment
    nnse=nseall
    eventsloc=seiseloc(:,1)
  endif

  ! output the final results
  open(unit=100,file=TRIM(outfname)//'event_location.dat',form='formatted',status='replace',action='write')
  write(unit=100,fmt="(I8)") nnse
  do id=1,nnse
    write(unit=100,fmt="(5(2XF18.6))") eventsloc((5*id-4):(5*id))
  enddo
  close(unit=100)

  deallocate(nsevt,nenpro,seiseloc,eventsloc)

  ! end timing
  call date_and_time(values=time_array_e)
  twallc=(time_array_e(3)-time_array_s(3))*24*3600+(time_array_e(5)-time_array_s(5))*3600+(time_array_e(6)-time_array_s(6))*60+(time_array_e(7)-time_array_s(7))+0.001*(time_array_e(8)-time_array_s(8))
  write(*,*) "Program running time:", twallc, "second."

end program migrationloc
