module variables
  use omp_lib
  use parameters
  implicit none
  

  integer(8),parameter :: np_image=(nc*np_nc)**3*merge(2,1,body_centered_cubic) ! average number of particles per image
  integer(8),parameter :: np_image_max=np_image*(nte*1./nt)**3*image_buffer
  integer(8),parameter :: np_tile_max=np_image/nnt**3*(nte*1./nt)**3*tile_buffer
  real,parameter :: dt_max=1
  real,parameter :: G_grid=1.0/6.0/pi

  integer itx,ity,itz,ifx,ify,ifz,iapm,napm(7)[*]
  integer(8) plan1,plan2(nteam),plan3(nteam,2:7),iplan1,iplan2(nteam),iplan3(nteam,2:7),plan0,iplan0 ! FFT plans
  integer(4) ijk(3,n_neighbor),iteam,ipm2,ipm3,ixyz2(3,nnt**3),ixyz3(6,(nnt*nns)**3)
  integer(OMP_integer_kind) ith,nth

  real dt[*],dt_old[*],dt_mid[*],dt_e,da[*],a_mid[*],f2_max_pp[*]
  real vmax(3),vmax_team(3),f2max,f2max_team,overhead_tile[*],overhead_image[*],sigma_vi,sigma_vi_new,svz(500,2),svr(100,2)
  real(8) testrho,std_vsim_c[*],std_vsim_res[*],std_vsim[*]
  
  real,allocatable :: Gk1(:,:,:),Gk2(:,:,:),Gk3_2(:,:,:),Gk3_4(:,:,:)
  real,allocatable :: Gk3_6(:,:,:),Gk3_8(:,:,:),Gk3_12(:,:,:),Gk3_16(:,:,:),Gk0(:,:,:) ! Green's functions
  integer(izipx),allocatable :: xp(:,:)[:,:,:],xp_new(:,:)
  real(4),allocatable :: vp(:,:)[:,:,:],vp_new(:,:)
  integer(4),allocatable :: rhoc(:,:,:,:,:,:)[:,:,:]
  !real(4),allocatable :: vfield(:,:,:,:,:,:,:)[:,:,:]
#ifdef PID
    integer(izipi),allocatable :: pid(:)[:,:,:],pid_new(:)
#endif

  integer,parameter :: ns3=(nnt*nns)**3
  real(4) speed(4,7,ns3)[*],speed_global(4,7,ns3,nn**3)

  real rho0(ng+2,ng,ng)
  complex rho0k(ng/2+1,ng,ng)
  equivalence(rho0,rho0k)
  real rho2(1-ngb:ngp+ngb+2,1-ngb:ngp+ngb,1-ngb:ngp+ngb,nteam)
  complex rho2k(ngt/2+1,ngt,ngt,nteam)
  equivalence(rho2,rho2k)

  real rho3_2(1-nfb(2):nfp(2)+nfb(2)+2,1-nfb(2):nfp(2)+nfb(2),1-nfb(2):nfp(2)+nfb(2),nteam)
  real rho3_4(1-nfb(3):nfp(3)+nfb(3)+2,1-nfb(3):nfp(3)+nfb(3),1-nfb(3):nfp(3)+nfb(3),nteam)
  real rho3_6(1-nfb(4):nfp(4)+nfb(4)+2,1-nfb(4):nfp(4)+nfb(4),1-nfb(4):nfp(4)+nfb(4),nteam)
  real rho3_8(1-nfb(5):nfp(5)+nfb(5)+2,1-nfb(5):nfp(5)+nfb(5),1-nfb(5):nfp(5)+nfb(5),nteam)
  real rho3_12(1-nfb(6):nfp(6)+nfb(6)+2,1-nfb(6):nfp(6)+nfb(6),1-nfb(6):nfp(6)+nfb(6),nteam)
  !real rho3_16(1-nfb(7):nfp(7)+nfb(7)+2,1-nfb(7):nfp(7)+nfb(7),1-nfb(7):nfp(7)+nfb(7))
  complex rho3k_2(nft(2)/2+1,nft(2),nft(2),nteam)
  complex rho3k_4(nft(3)/2+1,nft(3),nft(3),nteam)
  complex rho3k_6(nft(4)/2+1,nft(4),nft(4),nteam)
  complex rho3k_8(nft(5)/2+1,nft(5),nft(5),nteam)
  complex rho3k_12(nft(6)/2+1,nft(6),nft(6),nteam)
  !complex rho3k_16(nft(7)/2+1,nft(7),nft(7))
  equivalence(rho3_2,rho3k_2)
  equivalence(rho3_4,rho3k_4)
  equivalence(rho3_6,rho3k_6)
  equivalence(rho3_8,rho3k_8)
  equivalence(rho3_12,rho3k_12)
  !equivalence(rho3_16,rho3k_16)

  integer(8),dimension(1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt),codimension[nn,nn,*] :: idx_b_l,idx_b_r
  integer(8),dimension(nt,nt,nnt,nnt,nnt),codimension[nn,nn,*] :: ppl0,pplr,pprl,ppr0,ppl,ppr
  
  type type_pm
    integer pm_layer ! 1=coarse, 2=standard, 3=fine
    integer iapm
    integer nwork ! grid number to assign density. = nt,ngt,nft(iapm)
    integer nstart(ndim),nend(ndim) ! density index range; for cubefft, do FFT padding 
    integer tile1(ndim),tile2(ndim) ! tiles to work on
    integer nloop ! additional layer of coarse cells to loop over particles
    integer nex ! additional layers of grids for density assignment
    integer nphy ! physical density assignment per tile. =nt,ngp,nfp(iapm)
    real gridsize ! determines gravity strength
    real f2max ! get maximum force squared
    real sigv1,sigv2 ! in/out-come velocity conversion coefficient
    integer tile_shift,utile_shift(ndim) ! tile and subtile offset
    integer m1,m2 ! density assignment per tile. = [1,nt],[1-ngb,ngp+ngb],[1-nfb(iapm),nfp(iapm)+nfb(iapm)]
    integer nforce ! grid number to update velocity. = nc,ngp,nfp(iapm)
    integer m1phi(ndim),m2phi(ndim) ! potential index range. [-2,nc+3],
    integer nc1(ndim),nc2(ndim) ! coarse grid iteration
  endtype
  type(type_pm) pm

contains
  subroutine buffer_potential(pot,n)
    implicit none
    integer(8) n
    real pot(-2:n+3,-2:n+3,-2:n+3)[nn,nn,*]
    sync all
    pot(:0,:,:)=pot(n-2:n,:,:)[inx,icy,icz]
    pot(n+1:,:,:)=pot(1:3,:,:)[ipx,icy,icz]
    sync all
    pot(:,:0,:)=pot(:,n-2:n,:)[icx,iny,icz]
    pot(:,n+1:,:)=pot(:,1:3,:)[icx,ipy,icz]
    sync all
    pot(:,:,:0)=pot(:,:,n-2:n)[icx,icy,inz]
    pot(:,:,n+1:)=pot(:,:,1:3)[icx,icy,ipz]
    sync all
  endsubroutine

  subroutine spine_tile(rhoce,idx_ex_r,pp_l,pp_r,ppe_l,ppe_r)
    !! make a particle index (cumulative sumation) on tile
    !! used in update_particle, initial_conditions
    !! input:
    !! rhoce -- particle number density on tile, with 2x buffer depth
    !! output:
    !! idx_ex_r -- last extended index on extended right boundary
    !! ppe_r -- last extended index on physical right boundary
    !! pp_l -- first physical index on physical left boundary
    !! pp_r -- last physical index on physical right boundary
    !! ppe_l -- first extended index on physical left boundary
    use omp_lib
    implicit none
    integer(4),intent(in) :: rhoce(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
    integer(8),dimension(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb),intent(out) :: idx_ex_r
    integer(8),dimension(nt,nt),intent(out) :: pp_l,pp_r,ppe_l,ppe_r
    integer(8) nsum,np_phy
    integer igz,igy
    ! spine := yz-plane to record cumulative sum
    nsum=0
    do igz=1-2*ncb,nt+2*ncb
    do igy=1-2*ncb,nt+2*ncb
      nsum=nsum+sum(rhoce(:,igy,igz))
      idx_ex_r(igy,igz)=nsum ! right index
    enddo
    enddo
    !$omp paralleldo default(shared) private(igz,igy)
    do igz=1,nt
    do igy=1,nt
      ppe_r(igy,igz)=idx_ex_r(igy,igz)-sum(rhoce(nt+1:,igy,igz))
    enddo
    enddo
    !$omp endparalleldo
    nsum=0
    do igz=1,nt
    do igy=1,nt
      pp_l(igy,igz)=nsum+1
      np_phy=sum(rhoce(1:nt,igy,igz))
      nsum=nsum+np_phy
      pp_r(igy,igz)=nsum
      ppe_l(igy,igz)=ppe_r(igy,igz)-np_phy+1
    enddo
    enddo
  endsubroutine

  subroutine spine_image(rhoc,idx_b_l,idx_b_r,ppe_l0,ppe_lr,ppe_rl,ppe_r0,ppl,ppr)
    !! make a particle index (cumulative sumation) on image
    !! used in buffer_grid
    !! input:
    !! rhoc -- particle number density on image, with 1x buffer depth
    !! output:
    !! idx_b_l -- zeroth extended index on extended left boundary
    !! idx_b_r -- last extended index on extended right boundary
    !! ppl -- zeroth physical index on physical left boundary
    !! ppr -- last physical index on physical right boundary
    !! ppe_r0 -- last extended index on physical right boundary
    !! ppe_l0 -- zeroth extended index on physical left boundary
    !! ppe_rl -- last extended index on inner right boundary
    !! ppe_lr -- zeroth extended index on inner left boundary
    implicit none
    integer(4),intent(in) :: rhoc(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)
    integer(8),dimension(1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt),intent(out) :: idx_b_l,idx_b_r
    integer(8),dimension(nt,nt,nnt,nnt,nnt),intent(out) :: ppe_l0,ppe_lr,ppe_rl,ppe_r0,ppl,ppr
    integer(8) nsum,nsum_p,np_phy,ihz,ihy,ihx,igz,igy,ctile_mass(nnt,nnt,nnt),ctile_mass_p(nnt,nnt,nnt)
    ! spine := yz-plane to record cumulative sum
    nsum=0;nsum_p=0
    do ihz=1,nnt ! sum cumulative tile mass first
    do ihy=1,nnt
    do ihx=1,nnt
      ctile_mass(ihx,ihy,ihz)=nsum
      ctile_mass_p(ihx,ihy,ihz)=nsum_p
      nsum=nsum+sum(rhoc(:,:,:,ihx,ihy,ihz))
      nsum_p=nsum_p+sum(rhoc(1:nt,1:nt,1:nt,ihx,ihy,ihz))
    enddo
    enddo
    enddo
    !$omp paralleldo default(shared) private(ihz,ihy,ihx,nsum,igz,igy)
    do ihz=1,nnt ! calculate extended spine cumulative index on both sides
    do ihy=1,nnt
    do ihx=1,nnt
      nsum=ctile_mass(ihx,ihy,ihz)
      do igz=1-ncb,nt+ncb
      do igy=1-ncb,nt+ncb
        idx_b_l(igy,igz,ihx,ihy,ihz)=nsum
        nsum=nsum+sum(rhoc(:,igy,igz,ihx,ihy,ihz))
        idx_b_r(igy,igz,ihx,ihy,ihz)=nsum
      enddo
      enddo
    enddo
    enddo
    enddo
    !$omp endparalleldo
    !$omp paralleldo default(shared) private(ihz,ihy,ihx,nsum_p,igz,igy,np_phy)
    do ihz=1,nnt ! calculate physical spine
    do ihy=1,nnt
    do ihx=1,nnt
      nsum_p=ctile_mass_p(ihx,ihy,ihz)
      do igz=1,nt
      do igy=1,nt
        ppl(igy,igz,ihx,ihy,ihz)=nsum_p
        np_phy=sum(rhoc(1:nt,igy,igz,ihx,ihy,ihz))
        nsum_p=nsum_p+np_phy
        ppr(igy,igz,ihx,ihy,ihz)=nsum_p

        ppe_r0(igy,igz,ihx,ihy,ihz)=idx_b_r(igy,igz,ihx,ihy,ihz)-sum(rhoc(nt+1:,igy,igz,ihx,ihy,ihz))
        ppe_rl(igy,igz,ihx,ihy,ihz)= ppe_r0(igy,igz,ihx,ihy,ihz)-sum(rhoc(nt-ncb+1:nt,igy,igz,ihx,ihy,ihz))

        ppe_l0(igy,igz,ihx,ihy,ihz)=idx_b_l(igy,igz,ihx,ihy,ihz)+sum(rhoc(:0,igy,igz,ihx,ihy,ihz))
        ppe_lr(igy,igz,ihx,ihy,ihz)=ppe_l0(igy,igz,ihx,ihy,ihz) +sum(rhoc(1:ncb,igy,igz,ihx,ihy,ihz))
      enddo
      enddo
    enddo
    enddo
    enddo
    !$omp endparalleldo
  endsubroutine

  real function interp_sigmav(aa,rr)
    implicit none
    integer(8) ii,i1,i2
    real aa,rr,term_z,term_r
    i1=1
    i2=500
    do while (i2-i1>1)
      ii=(i1+i2)/2
      if (aa>svz(ii,1)) then
        i1=ii
      else
        i2=ii
      endif
    enddo
    term_z=svz(i1,2)+(svz(i2,2)-svz(i1,2))*(aa-svz(i1,1))/(svz(i2,1)-svz(i1,1))
    i1=1
    i2=100
    do while (i2-i1>1)
      ii=(i1+i2)/2
      if (rr>svz(ii,1)) then
        i1=ii
      else
        i2=ii
      endif
    enddo
    term_r=svr(i1,2)+(svr(i2,2)-svr(i1,2))*(rr-svr(i1,1))/(svr(i2,1)-svr(i1,1))
    interp_sigmav=term_z*term_r
    print*,term_z,term_r
  endfunction

endmodule
