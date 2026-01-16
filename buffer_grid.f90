subroutine buffer_grid
  use omp_lib
  use variables
  implicit none
  save
  call tic(3)
  call buffer_np
  !call buffer_vc
  call redistribute_cdm
  call toc(3)
endsubroutine

subroutine buffer_np
  use omp_lib
  use variables
  implicit none
  save
  call tic(7)
  if (head) then
    print*, 'buffer_np'
  endif
  sync all
  !!$omp workshare
  rhoc(:0,:,:,1,:,:)=rhoc(nt-ncb+1:nt,:,:,nnt,:,:)[inx,icy,icz]
  rhoc(:0,:,:,2:,:,:)=rhoc(nt-ncb+1:nt,:,:,:nnt-1,:,:)
  rhoc(nt+1:,:,:,nnt,:,:)=rhoc(1:ncb,:,:,1,:,:)[ipx,icy,icz]
  rhoc(nt+1:,:,:,:nnt-1,:,:)=rhoc(1:ncb,:,:,2:,:,:)
  !!$omp endworkshare
  sync all
  !!$omp workshare
  rhoc(:,:0,:,:,1,:)=rhoc(:,nt-ncb+1:nt,:,:,nnt,:)[icx,iny,icz]
  rhoc(:,:0,:,:,2:,:)=rhoc(:,nt-ncb+1:nt,:,:,1:nnt-1,:)
  rhoc(:,nt+1:,:,:,nnt,:)=rhoc(:,1:ncb,:,:,1,:)[icx,ipy,icz]
  rhoc(:,nt+1:,:,:,:nnt-1,:)=rhoc(:,1:ncb,:,1:,2:,:)
  !!$omp endworkshare
  sync all
  !!$omp workshare
  rhoc(:,:,:0,:,:,1)=rhoc(:,:,nt-ncb+1:nt,:,:,nnt)[icx,icy,inz]
  rhoc(:,:,:0,:,:,2:)=rhoc(:,:,nt-ncb+1:nt,:,:,:nnt-1)
  rhoc(:,:,nt+1:,:,:,nnt)=rhoc(:,:,1:ncb,:,:,1)[icx,icy,ipz]
  rhoc(:,:,nt+1:,:,:,:nnt-1)=rhoc(:,:,1:ncb,:,:,2:)
  !!$omp endworkshare
  sync all
  call toc(7)
endsubroutine

#ifdef old
subroutine buffer_vc
  use variables
  implicit none
  save
  call tic(8)
  if (head) print*, 'buffer_vc'
  sync all
  vfield(:,:0,:,:,1,:,:)=vfield(:,nt-ncb+1:nt,:,:,nnt,:,:)[inx,icy,icz]
  vfield(:,:0,:,:,2:,:,:)=vfield(:,nt-ncb+1:nt,:,:,:nnt-1,:,:)
  sync all
  vfield(:,nt+1:,:,:,nnt,:,:)=vfield(:,1:ncb,:,:,1,:,:)[ipx,icy,icz]
  vfield(:,nt+1:,:,:,:nnt-1,:,:)=vfield(:,1:ncb,:,:,2:,:,:)
  sync all
  vfield(:,:,:0,:,:,1,:)=vfield(:,:,nt-ncb+1:nt,:,:,nnt,:)[icx,iny,icz]
  vfield(:,:,:0,:,:,2:,:)=vfield(:,:,nt-ncb+1:nt,:,:,1:nnt-1,:)
  sync all
  vfield(:,:,nt+1:,:,:,nnt,:)=vfield(:,:,1:ncb,:,:,1,:)[icx,ipy,icz]
  vfield(:,:,nt+1:,:,:,:nnt-1,:)=vfield(:,:,1:ncb,:,1:,2:,:)
  sync all
  vfield(:,:,:,:0,:,:,1)=vfield(:,:,:,nt-ncb+1:nt,:,:,nnt)[icx,icy,inz]
  vfield(:,:,:,:0,:,:,2:)=vfield(:,:,:,nt-ncb+1:nt,:,:,:nnt-1)
  sync all
  vfield(:,:,:,nt+1:,:,:,nnt)=vfield(:,:,:,1:ncb,:,:,1)[icx,icy,ipz]
  vfield(:,:,:,nt+1:,:,:,:nnt-1)=vfield(:,:,:,1:ncb,:,:,2:)
  sync all
  call toc(8)
endsubroutine
#endif

subroutine redistribute_cdm
  use omp_lib
  use variables
  implicit none
  save
  integer(8),parameter :: unit8=1
  integer(8) nshift,checkxp0,checkxp1,checkpid0,checkpid1
  integer i,iy,iz
  call tic(9)
  if (head) then
    print*,''
    print*, 'redistribute_cdm'
    call system_clock(t1,t_rate)
  endif
  ! check
  overhead_image=sum(rhoc*unit8)/real(np_image_max,8)
  sync all
  if (head) then
    do i=2,nn**3
      overhead_image=max(overhead_image,overhead_image[i])
    enddo
  endif
  sync all
  overhead_image=overhead_image[1]
  sync all
  !print*,sum(rhoc*unit8),np_image_max,sum(rhoc*unit8)/real(np_image_max,8)
  if (head) then
    print*, '  image overhead',overhead_image*100,'% full'
    print*, '  consumed image_buffer =',overhead_image*image_buffer,'/',image_buffer
    tcat(51,istep)=overhead_image*image_buffer
  endif
  sync all

  ! shift to right
  !checkxp0=sum(xp*unit8)
  !checkpid0=sum(pid)
  nshift=np_image_max-sim%nplocal
  !$omp parallelsections default(shared)
  !$omp section
    xp(:,nshift+1:np_image_max)=xp(:,1:sim%nplocal)
    xp(:,1:nshift)=0
  !$omp section
    vp(:,nshift+1:np_image_max)=vp(:,1:sim%nplocal)
    vp(:,1:nshift)=0
# ifdef PID
    !$omp section
    pid(nshift+1:np_image_max)=pid(1:sim%nplocal)
    pid(1:nshift)=0
# endif
  !$omp endparallelsections
  
  ! shift back
  call spine_image(rhoc,idx_b_l,idx_b_r,ppl0,pplr,pprl,ppr0,ppl,ppr)
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
    do iz=1,nt
    do iy=1,nt
      xp(:,ppl0(iy,iz,itx,ity,itz)+1:ppr0(iy,iz,itx,ity,itz))=xp(:,nshift+ppl(iy,iz,itx,ity,itz)+1:nshift+ppr(iy,iz,itx,ity,itz))
      vp(:,ppl0(iy,iz,itx,ity,itz)+1:ppr0(iy,iz,itx,ity,itz))=vp(:,nshift+ppl(iy,iz,itx,ity,itz)+1:nshift+ppr(iy,iz,itx,ity,itz))
#     ifdef PID
        pid(ppl0(iy,iz,itx,ity,itz)+1:ppr0(iy,iz,itx,ity,itz))=pid(nshift+ppl(iy,iz,itx,ity,itz)+1:nshift+ppr(iy,iz,itx,ity,itz))
#     endif
    enddo
    enddo
  enddo
  enddo
  enddo
  xp(:,idx_b_r(nt+ncb,nt+ncb,nnt,nnt,nnt):)=0
  pid(idx_b_r(nt+ncb,nt+ncb,nnt,nnt,nnt):)=0
  
  if (head) then
    call system_clock(t2,t_rate)
    print*, '  elapsed time =',real(t2-t1)/t_rate,'secs'
    print*, ''
  endif
  sync all
  call toc(9)
endsubroutine
