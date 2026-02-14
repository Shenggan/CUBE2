subroutine kick
  use omp_lib
  use variables
  use pencil_fft
  use cubefft
  use iso_c_binding
  implicit none

  interface
      subroutine c_pp_force_kernel(isort, ires, ixyz3, apm3, ratio_sf, &
                                   rhoc, idx_b_r, xp, vp, &
                                   mass_p_cdm, a_mid, dt, &
                                   f2max_t, vmax_t, ptotal, ttotal)  &
                                   bind(C, name='c_pp_force_kernel')
          use iso_c_binding
          integer(c_int), dimension(*), intent(in) :: isort, ires, ixyz3
          real(c_float), dimension(7), intent(in) :: apm3
          integer(c_int), dimension(7), intent(in) :: ratio_sf
          integer(c_int), dimension(*), intent(in) :: rhoc
          integer(c_int64_t), dimension(*), intent(in) :: idx_b_r
          integer(c_int16_t), dimension(*), intent(in) :: xp
          real(c_float), dimension(*), intent(inout) :: vp
          real(c_float), value :: mass_p_cdm, a_mid, dt
          real(c_float), intent(out) :: f2max_t
          real(c_float), dimension(3), intent(out) :: vmax_t
          real(c_double), dimension(*), intent(out) :: ptotal, ttotal
      end subroutine c_pp_force_kernel
  end interface

  logical,parameter :: PM3=.true.
  logical,parameter :: PP=.true.
  integer inode,itile,np_phy,nl(3),nh(3),i1,i2,it,i3,isort(ns3),ires(ns3),rcp,npgrid,nptile
  integer,allocatable :: ll(:),hoc(:,:,:)
  integer(8),allocatable :: ip_local(:)
  real,allocatable :: rho_th(:,:,:),phi_th(:,:,:),force_th(:,:,:,:),xf(:,:),vf(:,:),af(:,:)
  real(8),allocatable :: ptotal(:),ttotal(:)
  real(8) pt
  real npc_max(ns3),overden_phy
  
  call tic(15)
  
  if (head) print*,'PM1' ! =====================================================
  call tic(11)
  iapm=1; pm%pm_layer=1; pm%iapm=iapm; pm%nwork=nt; pm%nstart(:)=1; pm%nend(:)=nc
  pm%tile1(:)=1; pm%tile2(:)=nnt; pm%nloop=1; pm%nex=2; pm%nphy=nt; pm%gridsize=ratio_cs
  pm%tile_shift=1; pm%utile_shift=0; pm%m1=1; pm%m2=nt; pm%nforce=nc; pm%m1phi(:)=-2; pm%m2phi(:)=nc+3
  pm%sigv1=sigma_vi; pm%sigv2=sigma_vi; pm%nc1(:)=1; pm%nc2(:)=nt; testrho=0; pm%f2max=0
  call particle_mesh(rho1)
  sim%dt_pm1=0.5*sqrt( 1. / (sqrt(pm%f2max)*a_mid*G_grid) )
  call toc(11)
  if (head) then
    print*,'  pm%f2max =',pm%f2max; print*,'    sum of PM1 density =',testrho
    print*,'  real time =',tcat(11,istep),'secs'; print*,''
  endif
  sync all

  if (head) print*,'PM2' ! =====================================================
  call tic(12)
  vmax=0; f2max=0; pt=0
  allocate(ptotal(nnt**3),ttotal(nnt**3)); ptotal=0
  !$omp paralleldo num_threads(nteam) schedule(dynamic,1) default(shared) &
  !$omp& private(itile,iteam,pm,rho_th,force_th,vmax_team,f2max_team) reduction(max:vmax,f2max) reduction(+:pt)
  do itile=1,nnt**3
    iteam=omp_get_thread_num()+1; iapm=1; pm%pm_layer=2; pm%iapm=iapm; pm%nwork=ngt; pm%nstart(:)=1-ngb
    pm%nend(:)=ngp+ngb; pm%nend(1)=pm%nend(1)+2; pm%tile1(:)=ixyz2(:,itile); pm%tile2(:)=ixyz2(:,itile)
    pm%nloop=ncb; pm%nex=ngb+1; pm%nphy=ngp; pm%gridsize=1; pm%tile_shift=0; pm%utile_shift=0; 
    pm%m1=1-ngb; pm%m2=ngp+ngb; pm%nforce=ngp; pm%m1phi(:)=1-ngb; pm%m2phi(:)=ngp+ngb; pm%m2phi(1)=pm%m2phi(1)+2
    pm%sigv1=sigma_vi; pm%sigv2=sigma_vi_new; pm%nc1(:)=1; pm%nc2(:)=nt
    allocate(rho_th(1-pm%nex:pm%nphy+pm%nex,1-pm%nex:pm%nphy+pm%nex,1-pm%nex:pm%nphy+pm%nex))
    allocate(force_th(ndim,0:pm%nforce+1,0:pm%nforce+1,0:pm%nforce+1))
    call pm_thread(rho2(:,:,:,iteam),rho_th,force_th,pm,itile,iteam,f2max_team,vmax_team)
    deallocate(rho_th,force_th)
    f2max=max(f2max,f2max_team); vmax=max(vmax,vmax_team)
  enddo
  !$omp endparalleldo
  sigma_vi=sigma_vi_new;
  sim%vsim2phys=(1.5/sim%a)*box*h0*100.*sqrt(omega_m)/ng_global
  sim%dt_pm2=0.5*sqrt( 1. / (sqrt(f2max)*a_mid*G_grid) )
  call toc(12)
  if (head) then
    print*,'  sum of PM2 density',sum(ptotal),pt; print*,'  vmax =',vmax; print*,'  f2max =',f2max
    !print*,'  team time =',ttotal
    print*,'  total team time =',sum(ttotal),'secs'; print*,'  real time =',tcat(12,istep),'secs'; print*,''
  endif
  deallocate(ptotal,ttotal)
  sync all
  
  if (PM3 .or. PP) then  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    napm=0
    do itile=1,ns3 ! count clustering
      nl=(ixyz3(1:3,itile)-1)*ntt+1; nh=ixyz3(1:3,itile)*ntt
      i1=ixyz3(4,itile); i2=ixyz3(5,itile); i3=ixyz3(6,itile); 
      np_phy=sum(rhoc(nl(1):nh(1),nl(2):nh(2),nl(3):nh(3),i1,i2,i3))
      overden_phy=(real(np_phy)/(ntt*np_nc)**3)-1
      npc_max(itile)=maxval(rhoc(nl(1):nh(1),nl(2):nh(2),nl(3):nh(3),i1,i2,i3))/64.
      if (npc_max(itile)<24) then
        ires(itile)=2
      elseif (npc_max(itile)<114) then
        ires(itile)=3
      elseif (npc_max(itile)<320) then
        ires(itile)=4
      elseif (npc_max(itile)<1350) then
        ires(itile)=5
      !elseif (npc_max(itile)<2000) then
      !  ires(itile)=5
      else
        ires(itile)=6
      endif
      !ires(itile)=4
      napm(ires(itile))=napm(ires(itile))+1
    enddo
    call indexed_sort(ns3,-npc_max,isort)
    !print*,npc_max
    !print*,isort
    !print*,npc_max(isort)
  endif  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sync all

  if (PM3) then ! =====================================================
    if (head) print*, 'PM3'
    call tic(13)
    allocate(ptotal(ns3),ttotal(ns3)); ptotal=0
    pm%f2max=0; f2_max_pp=0; vmax=0; f2max=0; pt=0; 
    !$omp paralleldo num_threads(nteam) schedule(dynamic,1)&
    !$omp& default(shared) private(it,iteam,itile,iapm,nl,nh,i1,i2,i3,pm,rho_th,force_th,vmax_team,f2max_team)&
    !$omp& reduction(max:vmax,f2max) reduction(+:pt)
    do it=1,ns3
      iteam=omp_get_thread_num()+1; itile=isort(it); iapm=ires(itile); nl=(ixyz3(1:3,itile)-1)*ntt+1; nh=ixyz3(1:3,itile)*ntt
      i1=ixyz3(4,itile); i2=ixyz3(5,itile); i3=ixyz3(6,itile); pm%pm_layer=3; pm%iapm=iapm; pm%nwork=nft(iapm)
      pm%nstart(:)=1-nfb(iapm); pm%nend(:)=nfp(iapm)+nfb(iapm); pm%nend(1)=pm%nend(1)+2
      pm%tile1(:)=ixyz3(4:6,itile); pm%tile2(:)=ixyz3(4:6,itile); pm%nloop=1; pm%nex=nfb(iapm)+1; pm%nphy=nfp(iapm); 
      pm%gridsize=1./ratio_sf(iapm); pm%tile_shift=0; pm%utile_shift=ixyz3(1:3,itile)-1; pm%m1=1-nfb(iapm); 
      pm%m2=nfp(iapm)+nfb(iapm);pm%nforce=nfp(iapm); pm%m1phi(:)=1-nfb(iapm); pm%m2phi(:)=nfp(iapm)+nfb(iapm); 
      pm%m2phi(1)=pm%m2phi(1)+2; pm%sigv1=sigma_vi; pm%sigv2=sigma_vi
      pm%nc1(:)=(ixyz3(1:3,itile)-1)*ntt+1; pm%nc2(:)=ixyz3(1:3,itile)*ntt
      allocate(rho_th(1-pm%nex:pm%nphy+pm%nex,1-pm%nex:pm%nphy+pm%nex,1-pm%nex:pm%nphy+pm%nex))
      allocate(force_th(ndim,0:pm%nforce+1,0:pm%nforce+1,0:pm%nforce+1))
      selectcase(iapm)
      case(2)
        call pm_thread(rho3_2(:,:,:,iteam),rho_th,force_th,pm,itile,iteam,f2max_team,vmax_team)
      case(3)
        call pm_thread(rho3_4(:,:,:,iteam),rho_th,force_th,pm,itile,iteam,f2max_team,vmax_team)
      case(4)
        call pm_thread(rho3_6(:,:,:,iteam),rho_th,force_th,pm,itile,iteam,f2max_team,vmax_team)
      case(5)
        call pm_thread(rho3_8(:,:,:,iteam),rho_th,force_th,pm,itile,iteam,f2max_team,vmax_team)
      case(6)
        call pm_thread(rho3_12(:,:,:,iteam),rho_th,force_th,pm,itile,iteam,f2max_team,vmax_team)
      endselect
      deallocate(rho_th,force_th)
      f2max=max(f2max,f2max_team); vmax=max(vmax,vmax_team)
    enddo
    !$omp endparalleldo
    sim%dt_pm3=0.5*sqrt( 1. / (sqrt(f2max)*a_mid*G_grid) )
    call toc(13)
    if (head) then
      print*,'  sum of PM3 density',this_image(),sum(ptotal),pt; print*,'  vmax =',vmax; print*,'  f2max =',f2max
      !print*,'  team time =',ttotal
      print*,'  total team time =',sum(ttotal),'secs'; print*,'  real time =',tcat(13,istep),'secs'; print*,''
    endif
    deallocate(ptotal,ttotal)
  endif !PM3
  sync all

  if (PP) then ! =====================================================
    if (head) print*,'PP'
    call tic(14)
    vmax=0; f2max=0; pt=0
    
    allocate(ptotal(nnt*nnt*nnt),ttotal(nnt*nnt*nnt))
    call c_pp_force_kernel(isort, ires, ixyz3, apm3, ratio_sf, &
                           rhoc, idx_b_r, xp, vp, &
                           sim%mass_p_cdm, a_mid, dt, &
                           f2max, vmax, ptotal, ttotal)

    ! allocate(ptotal(ns3),ttotal(ns3))
    ! !$omp paralleldo num_threads(nteam) schedule(dynamic,1)&
    ! !$omp& default(shared) private(it,iteam,itile,iapm,nl,nh,i1,i2,i3,pm,rho_th,force_th,vmax_team,f2max_team)&
    ! !$omp& private(rcp,npgrid,nptile,ll,hoc,ip_local,xf,vf,af)&
    ! !$omp& reduction(max:vmax,f2max) reduction(+:pt)
    ! do it=1,ns3
    !   iteam=omp_get_thread_num()+1; itile=isort(it); iapm=ires(itile);i1=ixyz3(4,itile); i2=ixyz3(5,itile); i3=ixyz3(6,itile);
    !   pm%pm_layer=3; pm%iapm=iapm; pm%nwork=nft(iapm); pm%nstart(:)=1-nfb(iapm); pm%nend(:)=nfp(iapm)+nfb(iapm); pm%nend(1)=pm%nend(1)+2
    !   pm%tile1(:)=ixyz3(4:6,itile); pm%tile2(:)=ixyz3(4:6,itile);pm%nloop=1; pm%nex=nfb(iapm)+1; pm%nphy=nfp(iapm); 
    !   pm%gridsize=1./ratio_sf(iapm); pm%tile_shift=0; pm%utile_shift=ixyz3(1:3,itile)-1; pm%m1=1-nfb(iapm); pm%m2=nfp(iapm)+nfb(iapm);pm%nforce=nfp(iapm); pm%m1phi(:)=1-nfb(iapm); pm%m2phi(:)=nfp(iapm)+nfb(iapm); pm%m2phi(1)=pm%m2phi(1)+2
    !   pm%sigv1=sigma_vi; pm%sigv2=sigma_vi; pm%nc1(:)=(ixyz3(1:3,itile)-1)*ntt+1; pm%nc2(:)=ixyz3(1:3,itile)*ntt
    !   nl=pm%nc1; nh=pm%nc2(:); rcp=ratio_sf(iapm); npgrid=ntt*rcp
    !   nptile=sum(rhoc(nl(1)-1:nh(1)+1,nl(2)-1:nh(2)+1,nl(3)-1:nh(3)+1,pm%tile1(1),pm%tile1(2),pm%tile1(3)))
    !   allocate(ll(nptile),hoc(1-rcp:npgrid+rcp,1-rcp:npgrid+rcp,1-rcp:npgrid+rcp))
    !   allocate(ip_local(nptile),xf(3,nptile),vf(3,nptile),af(3,nptile))
    !   call PP_Force(ll,hoc,ip_local,xf,vf,af,pm,itile,iteam,f2max_team,vmax_team)
    !   deallocate(ll,hoc,ip_local,xf,vf,af)
    !   f2max=max(f2max,f2max_team); vmax=max(vmax,vmax_team)
    ! enddo
    ! !$omp endparalleldo
    sim%dt_pp=0.5*sqrt( 1. / (sqrt(f2max)*a_mid*G_grid) )  
    call toc(14)
    if (head) then 
      print*,'  sum of PP count',sum(ptotal); print*,'  vmax =',vmax; print*,'  f2max =',f2max
      !print*,'  team time =',ttotal
      print*,'  total team time =',sum(ttotal),'secs'; print*,'  real time =',tcat(14,istep),'secs'; print*,''
    endif
    deallocate(ptotal,ttotal)
  endif ! PP
  sync all

  sim%dt_vmax=vbuf*14/maxval(vmax)
  sim%vz_max=vmax(3)
  if (head) then ! gather info from all nodes
    do inode=2,nn**3
      napm(:)=napm(:)+napm(:)[inode]
      sim%dt_pm1=min(sim%dt_pm1,sim[inode]%dt_pm1)
      sim%dt_pm2=min(sim%dt_pm2,sim[inode]%dt_pm2)
      sim%dt_pm3=min(sim%dt_pm1,sim[inode]%dt_pm3)
      sim%dt_pp=min(sim%dt_pp,sim[inode]%dt_pp)
      sim%dt_vmax=min(sim%dt_vmax,sim[inode]%dt_vmax)
    enddo
  endif
  sync all ! broadcast
  napm(:)=napm(:)[1]
  sim%dt_pm1=sim[1]%dt_pm1
  sim%dt_pm2=sim[1]%dt_pm2
  sim%dt_pm3=sim[1]%dt_pm3
  sim%dt_pp=sim[1]%dt_pp
  sim%dt_vmax=sim[1]%dt_vmax
  sync all

  call toc(15)
  if (head) then 
    tcat(61:67,istep)=real(napm,kind=kind(tcat))
    print*,'    napm =',napm
    print*,'  real time =',tcat(15,istep),'secs'
  endif

  contains




  subroutine PP_Force(ll,hoc,ip_local,xf,vf,af,pm,itile,iteam,f2max_t,vmax_t)
    use omp_lib
    implicit none
    type(type_pm) pm
    integer(4) i_neighbor,nptile,itile,iteam,iapm,rcp,npgrid,nl(3),nh(3),ncount,tteam1,tteam2,tteam_rate
    integer i,j,k,ilayer,np,l,idx(3),ip1,ip2,ipenz,ipeny,ipenx,nz,ny,nx,iz,iy,ix,n_iter,i_iter
    integer ll(:),hoc(1-ratio_sf(pm%iapm):,1-ratio_sf(pm%iapm):,1-ratio_sf(pm%iapm):)!,ip_local(:)
    integer(8) ip,nzero,ip_local(:)
    real xvec(3),rvec(3),rmag,fpp(3),pp_range,xf(:,:),vf(:,:),af(:,:),f2max_t,vmax_t(:)

    call system_clock(tteam1,tteam_rate)
    iapm=pm%iapm
    pp_range=apm3(iapm); rcp=ratio_sf(iapm); npgrid=ntt*rcp
    hoc=0; ll=0; ip1=0; af=0; f2max_t=0; vmax_t=0
    do k=pm%nc1(3)-1,pm%nc2(3)+1
    do j=pm%nc1(2)-1,pm%nc2(2)+1
    do i=pm%nc1(1)-1,pm%nc2(1)+1
      np=rhoc(i,j,k,pm%tile1(1),pm%tile1(2),pm%tile1(3))
      nzero=idx_b_r(j,k,pm%tile1(1),pm%tile1(2),pm%tile1(3))-sum(rhoc(i:,j,k,pm%tile1(1),pm%tile1(2),pm%tile1(3)))
      do l=1,np
        ip=nzero+l; ip1=ip1+1; ip_local(ip1)=ip
        xvec=[i,j,k]-1+(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
        xf(:,ip1)=xvec*ratio_cs
        !vf(:,ip1)=tan(pi*real(vp(:,ip))/real(nvbin-1))/(sqrt(pi/2)/(sigma_vi*vrel_boost))
        vf(:,ip1)=vp(:,ip)
        idx=floor(rcp*(xvec-pm%utile_shift*ntt))+1
        ll(ip1)=hoc(idx(1),idx(2),idx(3))
        hoc(idx(1),idx(2),idx(3))=ip1
      enddo
    enddo
    enddo
    enddo
    
    !self pp
    !print*,'self pp'
    !$omp paralleldo default(shared) num_threads(nnest) schedule(dynamic)&
    !$omp& private(k,j,i,ip1,ip2,rvec,rmag,fpp)
    do k=1,npgrid
    do j=1,npgrid
    do i=1,npgrid
      ip1=hoc(i,j,k)
      do while(ip1/=0) ! particle A
        ip2=ll(ip1)
        do while (ip2/=0) ! particle B in the same cell
          rvec=xf(:,ip2)-xf(:,ip1); rmag=norm2(rvec)
          if (0.<rmag .and. rmag<pp_range) then
            fpp=rvec/rmag*(F_ra(rmag,app)-F_ra(rmag,pp_range))
            af(:,ip1)=af(:,ip1)+fpp
            af(:,ip2)=af(:,ip2)-fpp
          endif
          ip2=ll(ip2)
        enddo ! particle B
        ip1=ll(ip1)
      enddo ! particle A
    enddo
    enddo
    enddo
    !$omp endparalleldo

    !print*,'neigh pp'
    do ipenz=0,1
    do ipeny=0,2 ! divide all pencils into 6 bundles 
    do ipenx=0,2
      nz=(npgrid+1-ipenz)/2+1
      ny=(npgrid+1-ipeny)/3+1
      nx=(npgrid+1-ipenx)/3+1
      n_iter=nz*ny*nx ! number of pencils
      !$omp paralleldo default(shared) num_threads(nnest) schedule(dynamic)&
      !$omp& private(i_iter,iz,iy,ix,k,j,i,ip1,ip2,rvec,rmag,fpp,i_neighbor,idx)
      do i_iter=0,n_iter-1 ! (loop over all [k,j] pencils of the bundle)
        !iz=i_iter/ny
        !iy=mod(i_iter,ny)
        iz=i_iter/(ny*nx)
        iy=(i_iter-ny*nx*iz)/nx
        ix=mod(i_iter,nx)
        k=ipenz+2*iz
        j=ipeny+3*iy
        i=ipenx+3*ix
        !do i=0,npgrid+1
          ip1=hoc(i,j,k)
          do while(ip1/=0) ! particle A
            do i_neighbor=1,n_neighbor ! neighbor cells
              idx=[i,j,k]+ijk(:,i_neighbor)
              if (minval(idx)<0 .or. maxval(idx)>npgrid+1) cycle
              ip2=hoc(idx(1),idx(2),idx(3))
              do while (ip2/=0) ! particle B
                rvec=xf(:,ip2)-xf(:,ip1); rmag=norm2(rvec)
                if (rmag<pp_range) then
                  fpp=rvec/rmag*(F_ra(rmag,app)-F_ra(rmag,pp_range))
                  af(:,ip1)=af(:,ip1)+fpp
                  af(:,ip2)=af(:,ip2)-fpp
                endif
                ip2=ll(ip2)
              enddo
            enddo
            ip1=ll(ip1)
          enddo
        !enddo ! i
      enddo ! i_iter (loop over k,j)
      !$omp endparalleldo
    enddo ! ipenx
    enddo ! ipeny
    enddo ! ipenz

    !print*,'update v'
    ncount=0
    af=af*sim%mass_p_cdm*sim%mass_p_cdm
    do k=1,npgrid
    do j=1,npgrid
    do i=1,npgrid
      ip1=hoc(i,j,k)
      do while (ip1/=0)
        ncount=ncount+1
        f2max_t=max(f2max_t,sum(af(:,ip1)**2))
        ip=ip_local(ip1)
        !vreal=vf(:,ip1)+af(:,ip1)*a_mid*dt/6/pi
        vp(:,ip)=vp(:,ip)+af(:,ip1)*a_mid*dt/6/pi
        vmax_t=max(vmax_t,abs(vp(:,ip)))
        !vp(:,ip)=nint(real(nvbin-1)*atan(sqrt(pi/2)/(sigma_vi*vrel_boost)*vreal)/pi,kind=izipv)
        ip1=ll(ip1)
      enddo
    enddo
    enddo
    enddo
    ptotal(itile)=ncount
    call system_clock(tteam2,tteam_rate)
    ttotal(itile)=real(tteam2-tteam1)/tteam_rate
  endsubroutine PP_Force



  !! call pm_thread(rho2(:,:,:,iteam),rho_th,force_th,pm,itile,iteam)
  subroutine pm_thread(density,rho_th,force_th,pm,itile,iteam,f2max_t,vmax_t)
    use omp_lib
    implicit none
    integer,parameter :: nlayer=3 ! thread save for TSC intepolation
    type(type_pm) pm
    real density(pm%nstart(1):,pm%nstart(2):,pm%nstart(3):)
    !complex density_k(pm%nwork/2+1,pm%nwork,pm%nwork)
    real rho_th(1-pm%nex:,1-pm%nex:,1-pm%nex:),force_th(1:,0:,0:,0:)
    integer(4) itile,iteam,iapm,i_0,j_0,k_0,i_n(4),j_n(4),k_n(4),tteam1,tteam2,tteam_rate
    integer i,j,k,l,np,idx0(ndim,p+1),ilayer,i1,i2,i3,ix,iy,iz,n1(ndim),n2(ndim)
    integer(8) ip,nzero
    real(8) xpos(ndim)
    real dx(ndim,p+1),l3(ndim),dv(ndim),vreal(ndim),f2max_t,vmax_t(:)
    !equivalence(density,density_k)

    call system_clock(tteam1,tteam_rate)
    iapm=pm%iapm; density=0
    do iz=pm%tile1(3),pm%tile2(3) ! interpolate particles to density field
    do iy=pm%tile1(2),pm%tile2(2)
    do ix=pm%tile1(1),pm%tile2(1)
      !print*,ix,iy,iz
      rho_th=0
      do ilayer=0,nlayer-1
        !print*,ilayer
        !$omp paralleldo default(shared) num_threads(nnest) schedule(dynamic,1)&
        !$omp& private(k,j,i,np,nzero,l,ip,xpos,idx0,l3,dx,i3,i2,i1)
        do k=pm%nc1(3)-pm%nloop+ilayer,pm%nc2(3)+pm%nloop,nlayer
        !do k=pm%nc1(3)-pm%nloop,pm%nc2(3)+pm%nloop
        do j=pm%nc1(2)-pm%nloop,pm%nc2(2)+pm%nloop
        do i=pm%nc1(1)-pm%nloop,pm%nc2(1)+pm%nloop
          np=rhoc(i,j,k,ix,iy,iz)
          nzero=idx_b_r(j,k,ix,iy,iz)-sum(rhoc(i:,j,k,ix,iy,iz))
          do l=1,np
            ip=nzero+l
            xpos=([i,j,k]-1)+(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
            xpos=xpos*ratio_cs/pm%gridsize-pm%utile_shift*nfp(iapm)
            idx0(:,2)=floor(xpos)+1
            idx0(:,1)=idx0(:,2)-1
            idx0(:,3)=idx0(:,2)+1
            l3=xpos-floor(xpos)
            dx(:,1)=(1-l3)**2/2
            dx(:,3)=l3**2/2
            dx(:,2)=1-dx(:,1)-dx(:,3)
            do i3=1,p+1
            do i2=1,p+1
            do i1=1,p+1
              rho_th(idx0(1,i1),idx0(2,i2),idx0(3,i3))=rho_th(idx0(1,i1),idx0(2,i2),idx0(3,i3))+dx(1,i1)*dx(2,i2)*dx(3,i3)
            enddo
            enddo
            enddo
          enddo
        enddo
        enddo
        enddo
        !$omp endparalleldo
      enddo ! ilayer
      !print*,'    from rho_ex',pm%m1,pm%m2
      n1(:)=pm%nstart+pm%tile_shift*([ix,iy,iz]-1)*pm%nphy
      n2(:)=n1(:)+pm%nwork-1
      !print*,'     to density (',n1(1),':',n2(1),',',n1(2),':',n2(2),',',n1(3),':',n2(3),')'
      density(n1(1):n2(1),n1(2):n2(2),n1(3):n2(3))=rho_th(pm%m1:pm%m2,pm%m1:pm%m2,pm%m1:pm%m2)*sim%mass_p_cdm
      ptotal(itile)=ptotal(itile)+sum(1d0*rho_th(1:pm%nphy,1:pm%nphy,1:pm%nphy))
      pt=pt+sum(1d0*rho_th(1:pm%nphy,1:pm%nphy,1:pm%nphy))
    enddo
    enddo
    enddo

    selectcase(pm%pm_layer)
    case(0)
    case(1)
    case(2)
      call sfftw_execute(plan2(iteam))
      rho2k(:,:,:,iteam)=rho2k(:,:,:,iteam)*Gk2
      call sfftw_execute(iplan2(iteam))
    case(3)
      call sfftw_execute(plan3(iteam,iapm))
      selectcase(iapm)
      case(2)
        rho3k_2(:,:,:,iteam)=rho3k_2(:,:,:,iteam)*Gk3_2
      case(3)
        rho3k_4(:,:,:,iteam)=rho3k_4(:,:,:,iteam)*Gk3_4
      case(4)
        rho3k_6(:,:,:,iteam)=rho3k_6(:,:,:,iteam)*Gk3_6
      case(5)
        rho3k_8(:,:,:,iteam)=rho3k_8(:,:,:,iteam)*Gk3_8
      case(6)
        rho3k_12(:,:,:,iteam)=rho3k_12(:,:,:,iteam)*Gk3_12
      endselect
      call sfftw_execute(iplan3(iteam,iapm))
    endselect
    call potential_to_force(force_th,density,pm)
    f2max_t=maxval(sum(force_th(:,1:pm%nforce,1:pm%nforce,1:pm%nforce)**2,1))
    vmax_t=0
    
    do iz=pm%tile1(3),pm%tile2(3) ! update velocity
    do iy=pm%tile1(2),pm%tile2(2)
    do ix=pm%tile1(1),pm%tile2(1)
      !$omp paralleldo num_threads(nnest) default(shared) schedule(dynamic)&
      !$omp& private(k,j,i,np,nzero,l,ip,xpos,idx0,l3,dx,vreal,dv,i3,i2,i1)&
      !$omp& reduction(max:vmax_t)
      do k=pm%nc1(3),pm%nc2(3)
      do j=pm%nc1(2),pm%nc2(2)
      do i=pm%nc1(1),pm%nc2(1)
        np=rhoc(i,j,k,ix,iy,iz)
        nzero=idx_b_r(j,k,ix,iy,iz)-sum(rhoc(i:,j,k,ix,iy,iz))
        do l=1,np ! loop over particle
          ip=nzero+l
          xpos=pm%tile_shift*([ix,iy,iz]-1)*nt+([i,j,k]-1)+(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
          xpos=xpos*ratio_cs/pm%gridsize-pm%utile_shift*nfp(iapm)
          idx0(:,2)=floor(xpos)+1
          idx0(:,1)=idx0(:,2)-1
          idx0(:,3)=idx0(:,2)+1
          l3=xpos-floor(xpos)
          dx(:,1)=(1-l3)**2/2
          dx(:,3)=l3**2/2
          dx(:,2)=1-dx(:,1)-dx(:,3)
          !vreal=tan(pi*real(vp(:,ip))/real(nvbin-1))/(sqrt(pi/2)/(pm%sigv1*vrel_boost))
          dv=0
          do i3=1,p+1
          do i2=1,p+1
          do i1=1,p+1
            dv=dv+force_th(:,idx0(1,i1),idx0(2,i2),idx0(3,i3))*dx(1,i1)*dx(2,i2)*dx(3,i3)
          enddo
          enddo
          enddo
          !vreal=vreal+dv*a_mid*dt/6/pi
          vp(:,ip)=vp(:,ip)+dv*a_mid*dt/6/pi
          !vmax_t=max(vmax_t,abs(vreal+vfield(:,i,j,k,ix,iy,iz)))
          vmax_t=max(vmax_t,abs(vp(:,ip)))
          !vp(:,ip)=nint(real(nvbin-1)*atan(sqrt(pi/2)/(pm%sigv2*vrel_boost)*vreal)/pi,kind=izipv)
        enddo
      enddo
      enddo
      enddo
      !$omp endparalleldo
    enddo
    enddo
    enddo
    call system_clock(tteam2,tteam_rate)
    ttotal(itile)=real(tteam2-tteam1)/tteam_rate
  endsubroutine pm_thread

  subroutine potential_to_force(force,phi,pm)
    use omp_lib
    implicit none
    type(type_pm) pm
    integer i_0,j_0,k_0,i_n(4),j_n(4),k_n(4)
    real force(1:,0:,0:,0:)
    real phi(pm%m1phi(1):,pm%m1phi(2):,pm%m1phi(3):)
    phi=phi/(real(pm%nwork)**3) ! dividing factor in cubefft
    !$omp paralleldo num_threads(nnest) default(shared) schedule(dynamic)&
    !$omp& private(k_0,j_0,i_0,i_n,j_n,k_n)
    do k_0=0,pm%nforce+1
    do j_0=0,pm%nforce+1
    do i_0=0,pm%nforce+1
      i_n=i_0+[-2,-1,1,2]
      j_n=j_0+[-2,-1,1,2]
      k_n=k_0+[-2,-1,1,2]
      force(1,i_0,j_0,k_0)=sum(phi(i_n,j_0,k_0)*weight)
      force(2,i_0,j_0,k_0)=sum(phi(i_0,j_n,k_0)*weight)
      force(3,i_0,j_0,k_0)=sum(phi(i_0,j_0,k_n)*weight)
    enddo
    enddo
    enddo
    !$omp endparalleldo
  endsubroutine

  subroutine indexed_sort(N,ARRIN,INDX)
    implicit none
    integer(4) N ! number of halos to sort
    integer(4) IR
    integer(4) INDX(:),INDXT,J,L,I
    real ARRIN(N),Q

    DO 11 J=1,N
      INDX(J)=J
    11    CONTINUE
    L=N/2+1
    IR=N
    10    CONTINUE
      IF(L.GT.1)THEN
        L=L-1
        INDXT=INDX(L)
        Q=ARRIN(INDXT)
      ELSE
        INDXT=INDX(IR)
        Q=ARRIN(INDXT)
        INDX(IR)=INDX(1)
        IR=IR-1
        IF(IR.EQ.1)THEN
          INDX(1)=INDXT
          RETURN
        ENDIF
      ENDIF
      I=L
      J=L+L
    20      IF(J.LE.IR)THEN
        IF(J.LT.IR)THEN
          IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
        ENDIF
        IF(Q.LT.ARRIN(INDX(J)))THEN
          INDX(I)=INDX(J)
          I=J
          J=J+J
        ELSE
          J=IR+1
        ENDIF
      GO TO 20
      ENDIF
      INDX(I)=INDXT
    GO TO 10
  endsubroutine indexed_sort
endsubroutine kick
