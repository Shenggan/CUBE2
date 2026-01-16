module pencil_fft
  use omp_lib 
  use parameters
  implicit none
  save

  integer(4),parameter :: NULL=0
  integer(8) planx,plany,planz,iplanx,iplany,iplanz
  real        rho1(nw,nw,nw)[nn,nn,*], rho1_copy(nw,nw,nw)
  complex     rho1k(nw/2,nw,nw)
  real        rxyz(nw_global+2  ,nw,npen)
  complex     cxyz(nw_global/2+1,nw,npen)
  !complex     cyyyxz(npen,nn,nn,nw/2+1,npen)
  complex     cyyxz(nw,     nn,nw/2+1,npen)[nn,nn,*]
  complex     czzzxy(npen,nn,nn,nw/2+1,npen)[nn,nn,*]
  !equivalence(cyyyxz,cyyxz,rho1,rho1k)
  !equivalence(cyyyxz,cyyxz)
  equivalence(rxyz,cxyz)
  equivalence(rho1_copy,rho1k)

  contains

  subroutine pencil_fft_forward
    use omp_lib
    implicit none
    save
    !sync all
    call system_clock(ttt1,t_rate)
    call system_clock(tt1,t_rate)
    call c2x
    sync all; call system_clock(tt2,t_rate); if (head) print*,'c2x',real(tt2-tt1)/t_rate

    call system_clock(tt1,t_rate)    
    call sfftw_execute(planx)
    sync all; call system_clock(tt2,t_rate); if (head) print*,'xtran',real(tt2-tt1)/t_rate

    call system_clock(tt1,t_rate)
    call x2y
    sync all; call system_clock(tt2,t_rate); if (head) print*,'x2y',real(tt2-tt1)/t_rate

    call system_clock(tt1,t_rate)
    call sfftw_execute(plany)
    sync all; call system_clock(tt2,t_rate); if (head) print*,'ytran',real(tt2-tt1)/t_rate

    call system_clock(tt1,t_rate)
    call y2z
    sync all; call system_clock(tt2,t_rate); if (head) print*,'y2z',real(tt2-tt1)/t_rate

    call system_clock(tt1,t_rate)
    call sfftw_execute(planz)
    sync all; call system_clock(tt2,t_rate); if (head) print*,'ztran',real(tt2-tt1)/t_rate

    call system_clock(tt1,t_rate)
    call z2y
    sync all; call system_clock(tt2,t_rate); if (head) print*,'z2y',real(tt2-tt1)/t_rate
    
    call system_clock(tt1,t_rate)
    call y2x
    sync all; call system_clock(tt2,t_rate); if (head) print*,'y2x',real(tt2-tt1)/t_rate
    sync all
    call system_clock(ttt2,t_rate); if (head) print*,'pen_forward time',real(ttt2-ttt1)/t_rate
  endsubroutine

  subroutine pencil_fft_backward
    use omp_lib
    implicit none
    save
    sync all
    call x2y
    call y2z
    call sfftw_execute(iplanz)
    call z2y
    call sfftw_execute(iplany)
    call y2x
    call sfftw_execute(iplanx)
    call x2c
    rho1=rho1/nw_global/nw_global/nw_global
    sync all
  endsubroutine

  subroutine c2x
    use omp_lib
    implicit none
    save
    integer(8) i1,islab
    complex,allocatable :: ctransfer(:,:,:)[:,:,:]
    !complex,allocatable :: rho1k_send(:,:,:)[:,:,:]

    if (.true.) then
    call system_clock(t1,t_rate)
    allocate(ctransfer(nw/2,nw,nn)[nn,nn,*])
    rho1_copy=rho1
    do islab=1,npen ! loop over cells in z, extract slabs
      ctransfer=rho1k(:,:,islab::npen) ! nn slabs of rho1k copied to ctransfer1
      sync all
      do i1=1,nn ! loop over parts in x, get slabs from each y node
        cxyz(nw*(i1-1)/2+1:nw*i1/2,:,islab)=ctransfer(:,:,m2)[i1,m1,m3]
      enddo
      sync all
    enddo
    deallocate(ctransfer)
    sync all; call system_clock(t2,t_rate); if (head) print*,'c2x 1',real(t2-t1)/t_rate
    endif

    if (.false.) then
    call system_clock(t1,t_rate)
    !allocate(rho1k_send(nw/2,nw,nw)[nn,nn,*])
    !rho1k_send=rho1k
    sync all
    do i1=1,nn
      !rho1k_send=rho1k(:,:,npen*(m2-1)+1:npen*m2); sync all
      rxyz(nw*(i1-1)+1:nw*i1,:,:)=rho1(:,:,npen*(m2-1)+1:npen*m2)[i1,m1,m3]
    enddo
    !deallocate(rho1k_send)
    sync all; call system_clock(t2,t_rate); if (head) print*,'c2x 2',real(t2-t1)/t_rate
    endif
  endsubroutine

  subroutine x2y
    use omp_lib
    implicit none
    save
    integer(8) i1,islab
    complex,allocatable :: ctransfer(:,:,:)[:,:,:],cyxzx(:,:,:,:)[:,:,:]
    
    if (.true.) then
    call system_clock(t1,t_rate)
    allocate(ctransfer(nw,nw/2+1,nn)[nn,nn,*])
    do islab=1,npen ! loop over z
      do i1=1,nn ! loop over squares in x direction
        ctransfer(:,:,i1)=transpose(cxyz(nw/2*(i1-1)+1:nw/2*i1+1,:,islab))
      enddo
      sync all
      do i1=1,nn
        cyyxz(:,i1,:,islab)=ctransfer(:,:,m1)[i1,m2,m3]
      enddo
      sync all
    enddo
    deallocate(ctransfer)
    sync all; call system_clock(t2,t_rate); if (head) print*,'x2y 1',real(t2-t1)/t_rate
    endif

    if (.false.) then
    call system_clock(t1,t_rate)
    allocate(cyxzx(nw,nw/2+1,npen,nn)[nn,nn,*])
    !!$omp paralleldo default(shared) schedule(dynamic) private(islab,i1)
    do islab=1,npen
      do i1=1,nn
        cyxzx(:,:,islab,i1)=transpose(cxyz(nw/2*(i1-1)+1:nw/2*i1+1,:,islab))
      enddo
    enddo
    !!$omp endparalleldo
    sync all; call system_clock(t2,t_rate); if (head) print*,'x2y 2',real(t2-t1)/t_rate

    call system_clock(t1,t_rate)
    do i1=1,nn
      cyyxz(:,i1,:,:)=cyxzx(:,:,:,m1)[i1,m2,m3]
    enddo
    deallocate(cyxzx)
    sync all; call system_clock(t2,t_rate); if (head) print*,'x2y 3',real(t2-t1)/t_rate
    endif
  endsubroutine

  subroutine y2z
    use omp_lib
    implicit none
    save
    integer(8) i1,i2,islab
    complex,allocatable :: ctransfer(:,:,:,:)[:,:,:]!,cyyyxz_t(:,:,:,:,:)[:,:,:]
    allocate(ctransfer(npen,npen,nn,nn)[nn,nn,*])
    !allocate(cyyyxz_t(npen,nn,nn,nw/2+1,npen)[nn,nn,*])
    
    if (.true.) then
    call system_clock(t1,t_rate)
    do islab=1,nw/2+1 ! loop over slices in x direction
      do i2=1,nn
      do i1=1,nn
        ctransfer(:,:,i1,i2)=transpose(cyyxz(npen*(i1-1)+1:npen*i1,i2,islab,:))
      enddo
      enddo
      sync all
      do i2=1,nn
      do i1=1,nn
        czzzxy(:,i1,i2,islab,:)=ctransfer(:,:,m2,m3)[m1,i1,i2]
      enddo
      enddo
      sync all
    enddo
    sync all; call system_clock(t2,t_rate); if (head) print*,'y2z 1',real(t2-t1)/t_rate
    endif
    
    if (.false.) then
    call system_clock(t1,t_rate)
    !!$omp paralleldo default(shared) schedule(dynamic) private(islab,i2,i1)
    do islab=1,nw/2+1
      do i2=1,nn
      do i1=1,nn
        !cyyyxz_t(:,i1,i2,islab,:)=transpose(cyyyxz(:,i1,i2,islab,:))
        cyyxz(npen*(i1-1)+1:npen*i1,i2,islab,:)=transpose(cyyxz(npen*(i1-1)+1:npen*i1,i2,islab,:))
      enddo
      enddo
    enddo
    !!$omp endparalleldo
    sync all; call system_clock(t2,t_rate); if (head) print*,'y2z 2',real(t2-t1)/t_rate
    
    call system_clock(t1,t_rate)
    do i2=1,nn
    do i1=1,nn
      !czzzxy(:,i1,i2,:,:)=cyyyxz_t(:,m2,m3,:,:)[m1,i1,i2]
      czzzxy(:,i1,i2,:,:)=cyyxz(npen*(m2-1)+1:npen*m2,m3,:,:)[m1,i1,i2]
    enddo
    enddo
    sync all; call system_clock(t2,t_rate); if (head) print*,'y2z 3',real(t2-t1)/t_rate
    endif
    deallocate(ctransfer)
    !deallocate(cyyyxz_t)
  endsubroutine

  subroutine z2y
    implicit none
    save
    integer(8) i1,i2,islab
    complex,allocatable :: ctransfer(:,:,:,:)[:,:,:]
    if (.true.) then
    allocate(ctransfer(npen,npen,nn,nn)[nn,nn,*])
    do islab=1,nw/2+1 ! loop over slices in x direction
      do i2=1,nn
      do i1=1,nn
        ctransfer(:,:,i1,i2)=transpose(czzzxy(:,i1,i2,islab,:))
      enddo
      enddo
      sync all
      do i2=1,nn
      do i1=1,nn
        !cyyyxz(:,i1,i2,islab,:)=ctransfer(:,:,m2,m3)[m1,i1,i2]
        cyyxz(npen*(i1-1)+1:npen*i1,i2,islab,:)=ctransfer(:,:,m2,m3)[m1,i1,i2]
      enddo
      enddo
      sync all
    enddo
    deallocate(ctransfer)
    endif

    if (.false.) then
    do islab=1,nw/2+1
    do i2=1,nn
    do i1=1,nn
      czzzxy(:,i1,i2,islab,:)=transpose(czzzxy(:,i1,i2,islab,:))
    enddo
    enddo
    enddo
    sync all
    do i2=1,nn
    do i1=1,nn
      cyyxz(npen*(i1-1)+1:npen*i1,i2,:,:)=czzzxy(:,m2,m3,:,:)[m1,i1,i2]
    enddo
    enddo
    endif
    sync all
  endsubroutine

  subroutine y2x
    implicit none
    save
    integer(8) i1,islab
    complex,allocatable :: ctransfer(:,:,:)[:,:,:],cxyyz(:,:,:,:)[:,:,:]
    if (.true.) then
    allocate(ctransfer(nw/2+1,nw,nn)[nn,nn,*])
    do islab=1,npen ! loop over z
      do i1=1,nn ! loop over squares in x direction
        ctransfer(:,:,i1)=transpose(cyyxz(:,i1,:,islab))
      enddo
      sync all
      do i1=1,nn
        cxyz(nw/2*(i1-1)+1:nw/2*i1+1,:,islab)=ctransfer(:,:,m1)[i1,m2,m3]
      enddo
      sync all
    enddo
    deallocate(ctransfer)
    endif

    if (.false.) then
    allocate(cxyyz(nw/2+1,nn,nw,npen)[nn,nn,*])
    do islab=1,npen
      do i1=1,nn
        cxyyz(:,i1,:,islab)=transpose(cyyxz(:,i1,:,islab))
      enddo
    enddo
    sync all
    do i1=1,nn
      cxyz(nw/2*(i1-1)+1:nw/2*i1+1,:,:)=cxyyz(:,m1,:,:)[i1,m2,m3]
    enddo
    sync all
    deallocate(cxyyz)
    endif
  endsubroutine

  subroutine x2c
    implicit none
    save
    integer(8) i1,islab
    real,allocatable :: rtransfer(:,:,:)[:,:,:], rxyz_send(:,:,:)[:,:,:]
    if (.true.) then
    allocate(rtransfer(nw,nw,nn)[nn,nn,*])
    do islab=1,npen
      do i1=1,nn
        rtransfer(:,:,i1)=rxyz(nw*(i1-1)+1:nw*i1,:,islab)
      enddo
      sync all
      do i1=1,nn
        rho1(:,:,islab+(i1-1)*npen)=rtransfer(:,:,m1)[m2,i1,m3]
      enddo
      sync all
    enddo
    deallocate(rtransfer)
    endif

    if (.false.) then
      allocate(rxyz_send(nw_global+2  ,nw,npen)[nn,nn,*])
      rxyz_send=rxyz; sync all
      do i1=1,nn
        rho1(:,:,(i1-1)*npen+1:i1*npen)=rxyz_send(nw*(m1-1)+1:nw*m1,:,:)[m2,i1,m3]
      enddo
      deallocate(rxyz_send)
    endif
    sync all
  endsubroutine

  subroutine create_penfft_plan
    implicit none
    save
    include 'fftw3.f'
    call sfftw_plan_many_dft_r2c(planx,1,nw*nn,nw*npen,cxyz,NULL,1,nw*nn+2,cxyz,NULL,1,nw*nn/2+1,FFTW_MEASURE)
    call sfftw_plan_many_dft_c2r(iplanx,1,nw*nn,nw*npen,cxyz,NULL,1,nw*nn/2+1,cxyz,NULL,1,nw*nn+2,FFTW_MEASURE)
    call sfftw_plan_many_dft(plany,1,nw*nn,(nw/2+1)*npen,cyyxz,NULL,1,nw*nn,cyyxz,NULL,1,nw*nn,FFTW_FORWARD,FFTW_MEASURE)
    call sfftw_plan_many_dft(iplany,1,nw*nn,(nw/2+1)*npen,cyyxz,NULL,1,nw*nn,cyyxz,NULL,1,nw*nn,FFTW_BACKWARD,FFTW_MEASURE)
    call sfftw_plan_many_dft(planz,1,nw*nn,(nw/2+1)*npen,czzzxy,NULL,1,nw*nn,czzzxy,NULL,1,nw*nn,FFTW_FORWARD,FFTW_MEASURE)
    call sfftw_plan_many_dft(iplanz,1,nw*nn,(nw/2+1)*npen,czzzxy,NULL,1,nw*nn,czzzxy,NULL,1,nw*nn,FFTW_BACKWARD,FFTW_MEASURE)
    sync all
  endsubroutine

  subroutine destroy_penfft_plan
    implicit none
    save
    include 'fftw3.f'
    call sfftw_destroy_plan(planx)
    call sfftw_destroy_plan(iplanx)
    call sfftw_destroy_plan(plany)
    call sfftw_destroy_plan(iplany)
    call sfftw_destroy_plan(planz)
    call sfftw_destroy_plan(iplanz)
    sync all
  endsubroutine
endmodule
