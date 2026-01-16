!! read existing random seed to generate random numbers
!#define READ_SEED
#define USE_PKIC

!! read same random noise to generate same initial conditions
!#define READ_NOISE
! #define READ_RHO
!#define macbook
!#define ELUCID_big

program initial_conditions
  use omp_lib
  use pencil_fft
  use iso_fortran_env, only : int64
  implicit none

  logical,parameter :: correct_kernel=.false.
  logical,parameter :: write_potential=.true.
  logical,parameter :: write_phik=.true.
  logical,parameter :: norm_as_sigma8=.true.
  integer,parameter :: nk=219 ! transfer function length

  !! buffer
  integer,parameter :: n_buffer=1 ! depth in coarse cells
  integer,parameter :: nb=16 ! replace ngb
  real tf(7,nk) ! transfer function array

  integer(8) i,j,k,ip,l,nzero,kg,jg,ig,ii,jj,kk,itx,ity,itz,idx,imove,g(3),iq(3)
  integer(4) seedsize,ilayer,nlayer,istat,icore
  real kr,kx,ky,kz,kmax,temp_r,temp_theta,pow,phi8,temp8[*]
  real(8) v8,norm_As,xq(3),gradphi(3),vreal(3)
  integer(int64) time64

  integer(4),allocatable :: iseed(:),rhoce(:,:,:)
  real,allocatable :: phi(:,:,:)[:,:,:],rseed_all(:,:)!,vfield(:,:,:,:)
  complex,allocatable :: delta_k(:,:,:)

  !! zip format arrays
  integer(8),parameter :: npt=nt*np_nc ! np (number of particle) / tile / dimension (dim)
  integer(8),parameter :: npb=n_buffer*np_nc ! np / buffer depth
  integer(8),parameter :: npmax=2*(npt+2*npb)**3 ! maximum np in memory
  integer(1),allocatable :: rholocal(:,:,:)
  integer(8) idx_ex_r(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb) ! index_extended_rhs
  integer(8),dimension(nt,nt) :: pp_l,pp_r,ppe_l,ppe_r ! for particle index

  integer(izipx),allocatable :: xp(:,:)
  real(4),allocatable :: vp(:,:)
#ifdef PID
    integer(8),allocatable :: pid(:)
#endif
  real grad_max(3)[*],vmax(3),vf ! max of gradient of phi; max velocity; velocity factor
  real diff_x(2,2,2),diff_y(2,2,2),diff_z(2,2,2)
  real(4) svz(500,2),svr(100,2) ! for velocity conversion between integers and reals
  real(8) sigma_vc,sigma_vf ! coarse- and fine- grid based velocity dispersion
  real(8) std_vsim_c,std_vsim_res,std_vsim ! velocity dispersion measured in simulation

  call system_clock(ttt1,t_rate)
  if (head .or. image==nn**3) print*, this_image(),num_images()
  call omp_set_num_threads(ncore) !! set number of cores for OpenMP
  call geometry !! initialize some basic variables and parameters
  call system('mkdir -p '//opath//'image'//image2str(image)) !! create output directory

  sim%cur_checkpoint=1 !! set current checkpoint to 1
  open(16,file='./z_checkpoint.txt',status='old') !! open redshift list to do checkpoint
  read(16,fmt='(f8.4)') z_checkpoint(sim%cur_checkpoint)
  close(16)

  if (head) then
    print*, ''
    print*, 'CUBE IC run on',int(nn**3,2),' images x',int(ncore,1),'cores'
    print*, 'IC resolution: ngic =', int(ngic*nn,2)
    print*, 'CUBE resolution: ng =', int(ng*nn,2)
    print*, 'at redshift=',z_checkpoint(sim%cur_checkpoint)
    if (body_centered_cubic) then
      print*, 'To genterate npglobal = 2 x', int(np_nc*ng/ratio_cs*nn,2),'^3'
    else
      print*, 'To genterate npglobal =', int(np_nc*ng/ratio_cs*nn,2),'^3'
    endif
    print*, 'Box size', box
    print*, 'body_centered_cubic =',body_centered_cubic
    print*, 'output: ', opath
    print*, '-----------------------------------------'
    call system('mkdir -p '//opath//'code')
    call system('cp ./utilities/*.f90 '//opath//'code/')
    call system('cp ./*.f90 '//opath//'code/')
    call system('cp ./z_*.txt '//opath//'code/')
    print*, 'export PYTHONUNBUFFERED=1'
    call system('export PYTHONUNBUFFERED=1')
    print*, 'python ./ic.py'
    call system('python ./ic.py')
  endif

  allocate(phi(-nb:ngic+nb+1,-nb:ngic+nb+1,-nb:ngic+nb+1)[nn,nn,*])
  !allocate(vfield(3,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb))
  allocate(rhoce(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb))
  allocate(rholocal(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb))
  allocate(xp(3,npmax),vp(3,npmax))
  allocate(delta_k(nyquist+1,nw,npen))
#ifdef PID
  allocate(pid(npmax))
#endif
  sync all
  !! initialize the simulation information in 'sim' structure
  sim%nplocal=0 ! local (on this image) CDM particle number
  sim%a=1./(1+z_checkpoint(sim%cur_checkpoint)) ! scale factor of the universe
  sim%t=0 ! time
  sim%tau=0 ! conformal time

  sim%timestep=0 ! number of timesteps
  sim%dt_pm1=9999;  sim%dt_pm2=9999;  sim%dt_pm3=9999;  sim%dt_pp=9999;  sim%dt_vmax=9999

  sim%cur_checkpoint=1 ! current checkpoint
  sim%box=box ! box size in Mpc/h
  sim%image=image ! image number
  sim%nn=nn ! image per (/) dim
  sim%nnt=nnt ! number of tiles / image / dim
  sim%nt=nt ! number of coarse grid / tile / dim
  sim%ncb=ncb ! buffer depth in coarse cell; always 6
  sim%izipx=izipx ! integer format for CDM position
  sim%izipv=izipv ! integer format for CDM velocity

  sim%h0=h0 ! Hubble
  sim%omega_m=omega_m
  sim%omega_l=omega_l
  sim%s8=s8 ! sigma_8
  sim%vsim2phys=(150./sim%a)*box*sqrt(omega_m)/ng_global ! velocity unit
  sim%z_i=z_checkpoint(sim%cur_checkpoint) ! initial redshift
  sync all
  phi=0
  tf=0
  sync all

  if (head) print*,''
  if (head) print*,'Creating FFT plans'
  call system_clock(t1,t_rate)
  call create_penfft_plan
  call system_clock(t2,t_rate)
  if (head) print*, '  elapsed time =',real(t2-t1)/t_rate,'secs';
  sync all

# ifdef READ_RHO
   if (head) print*,''
   if (head) print*,'Reading existed delta_L'
   open(11,file=output_dir()//'delta_L'//output_suffix(),access='stream')
   read(11) rho1
   close(11)
   rho1=rho1*Dgrow(sim%a)
   sync all
   if (head .or. image==nn**3) print*,'  delta_L',rho1(1:4,1,1)
   if (head .or. image==nn**3) print*,'  rms of delta',sqrt(sum(rho1**2*1.d0)/nw_global/nw_global/nw_global)

   call pencil_fft_forward ! Fourier transform ! transforms rho1 to cxyz
   if (head) cxyz(1,1,1)=0 ! zero frequency
   delta_k=cxyz
   sync all
# else

#ifdef USE_PKIC
   if(head) print*,'    read Pk_ic from '//opath//'IC/Pcb_ic.txt'
   open(11,file=opath//'IC/Pcb_ic.txt',form='formatted')
   read(11,*) tf(:2,:)
   close(11)
   
   tf(2,:)= sqrt( tf(2,:) * (( real(nw_global)/real(box) )**3))
#else
   if(head) print*,'    read TF_ic from IC/IC_rescaled_transfer_z'//trim(adjustl(str_z))//'.txt'
   open(11,file='IC/IC_rescaled_transfer_z'//trim(adjustl(str_z))//'.txt',form='formatted')
   read(11,*) tf
   close(11)
   print*,tf(:2,:)
   tf(2,:)= sqrt(tf(1,:)**(n_s) * (( real(nw_global)/real(box) )**3) ) * tf(2,:)
  ! ! transferfnc --------------------------------------
  ! if (head) print*,''
  ! if (head) print*,'Transfer function'
  ! call system_clock(t1,t_rate)
  ! open(11,file='./tf_wmap9/camb_82031666_transfer_out_z0.dat.txt',form='formatted')
  ! read(11,*) tf
  ! close(11)
  ! tf(2,:)= tf(1,:)**(3+n_s) * tf(2,:)**2 / (2*pi**2)
  if (norm_as_sigma8) then
    if (head) print*,'  scale amplitude as sigma_8'
    if (head) print*,'  sigma_8 =',s8
    v8=0; tf(4,1)=tf(1,2)/2 ! calculate v8
    do k=2,nk-1
      tf(4,k)=(tf(1,k+1)-tf(1,k-1))/2
    enddo
    tf(4,nk)=tf(1,nk)-tf(1,nk-1)    
    kmax=2*pi*sqrt(3.)*nyquist/box
    do k=1,nk
      if (tf(1,k)>kmax) exit
      v8=v8+tf(2,k)*tophat(tf(1,k)*8)**2*tf(4,k)/tf(1,k)
    enddo
    if (head) print*, '  v8, (s8^2/v8) =', v8, s8**2/v8
    tf(2,:)=tf(2,:)*(s8**2/v8)*Dgrow(sim%a)**2
  else
    if (head) print*,'  scale as A_s'
    !norm_As=2.*pi**2*h0**4*(h0/0.05)**(n_s-1)
    norm_As=3.878
    if (head) print*,'  norm_As =',norm_As
    if (head) print*,'  A_s =',A_s
    tf(2,:)=A_s*norm_As*tf(2,:)*Dgrow(sim%a)**2
  endif
#endif
  sync all



   ! noisemap -------------------------------------
   if (head) print*,''
   if (head) print*,'Generating random noise'
# ifdef READ_NOISE
   open(11,file=output_dir()//'noise'//output_suffix(),access='stream')
   read(11) rho1
   close(11)
   if (head .or. this_image()==num_images()) print*, '  READ IN NOISE MAP:', rho1(1,1,1), rho1(ngic,ngic,ngic)
# else
   call random_seed(size=seedsize) ! generate random seed
   if (head) print*,'  min seedsize =', seedsize
   seedsize=max(seedsize,12)
   allocate(iseed(seedsize))
   allocate(rseed_all(seedsize,nn**3))
#ifdef READ_SEED
   if (head) print*, '  Copy and read seeds from ../confings/'
   open(11,file=output_dir()//'seed'//output_suffix(),status='old',access='stream')
   read(11) iseed
   close(11)
   ! Input iseed
   call random_seed(put=iseed)
   if (head) print*, 'iseed', iseed
#else
   ! Generate at least 12 seeds according to system clock
   call system_clock(time64)
   do i = 1, seedsize
      iseed(i) = lcg(time64) + image*137
      !print*,'time64,iseed(',int(i,1),')=',time64,iseed(i)
   enddo
   ! Input iseed to system
   call random_seed(put=iseed) ! generate random seed
   ! Write iseed into file
   open(11,file=output_dir()//'seed'//output_suffix(),status='replace',access='stream')
   write(11) iseed
   close(11)
#endif
   deallocate(iseed)
   deallocate(rseed_all)
   sync all
   call random_number(rho1) ! generate random numbers ! rho1 is defined in pencil_fft ! tophat [0,1) distribution
   open(11,file=output_dir()//'noise'//output_suffix(),status='replace',access='stream')
   write(11) rho1 ! write random numbers into a file
   close(11)
   if (head .or. this_image()==num_images()) print*, '  noise',int(image,1),rho1(1:2,1,1)
# endif
   sync all

  ! Box-Muller transform ----------------------------------------------
  ! convert to standard Normal (Gaussian) distribution
  if (head) print*,'  Box-Muller transform'
  sync all
  !$omp paralleldo&
  !$omp& default(shared) &
  !$omp& private(k,j,i,temp_theta,temp_r)
  do k=1,ngic
  do j=1,ngic
  do i=1,ngic,2
    temp_theta=2*pi*rho1(i,j,k)
    temp_r=sqrt(-2*log(1-rho1(i+1,j,k)))
    rho1(i,j,k)=temp_r*cos(temp_theta)
    rho1(i+1,j,k)=temp_r*sin(temp_theta)
  enddo
  enddo
  enddo
  !$omp endparalleldo
  sync all
  call system_clock(t2,t_rate)
  sync all
  if (head) print*, '  elapsed time =',real(t2-t1)/t_rate,'secs';
  sync all
  ! delta_field ----------------------------------------------------
  ! apply transfer function in Fourier space (Wiener filter) and get delta_L
  if (head) print*, ''
  if (head) print*, 'delta field'
  call system_clock(t1,t_rate)
  if (head) print*, '  ftran'
  !sync all
  !print*, this_image(), rho1(1:4,1,1)
  sync all
  call pencil_fft_forward ! Fourier transform ! transforms rho1 to cxyz
  sync all

  if (head .or. image==nn**3) then
    print*, cxyz(1:6,1,1)
  endif
  sync all;
  !print*,'ifft'
  !sync all
  !call pencil_fft_backward
  !sync all
  !print*, this_image(), rho1(1:4,1,1)
  !sync all

  if (head) print*, '  Wiener filter'
  !$omp paralleldo&
  !$omp& default(shared) &
  !$omp& private(k,j,i,kg,jg,ig,kz,ky,kx,kr,pow)
  do k=1,npen
  do j=1,ngic
  do i=1,nyquist+1
    ! global grid in Fourier space for i,j,k
    kg=(nn*(icz-1)+icy-1)*npen+k
    jg=(icx-1)*ngic+j
    ig=i
    kz=mod(kg+nyquist-1,nw_global)-nyquist
    ky=mod(jg+nyquist-1,nw_global)-nyquist
    kx=ig-1
    kr=2*real(pi,kind=4)*sqrt(kx**2+ky**2+kz**2)/box ! kr is |k_n|

    cxyz(i,j,k) = cxyz(i,j,k) * interp_tf(kr,1,2)

  enddo
  enddo
  enddo
  !$omp endparalleldo
  if (head) cxyz(1,1,1)=0 ! zero frequency
  delta_k=cxyz ! backup (Fourier) k-space delta_L
  sync all

  if (head) print*,'  btran'
  call pencil_fft_backward ! inverse Fourier transform ! transforms cxyz to rho1
  if (head .or. image==nn**3) print*,'  delta_L',rho1(1:4,1,1)
  if (head .or. image==nn**3) print*,'  rms of delta',sqrt(sum(rho1**2*1.d0)/nw_global/nw_global/nw_global)
  sync all

  if (head) print*,'  write delta_L into file'
  if (head) print*,'  growth factor Dgrow(',sim%a,') =',Dgrow(sim%a) ! growth factor

  open(11,file=output_dir()//'delta_L'//output_suffix(),status='replace',access='stream')
  do i=1,ngic
    write(11) rho1(:,:,i)/Dgrow(sim%a) ! write layer by layer to avoid bug
  enddo
  close(11); sync all

  call system_clock(t2,t_rate)
  if (head) print*, '  elapsed time =',real(t2-t1)/t_rate,'secs';
#endif

  ! Potential field ----------------------------------------------------
  ! solve Poisson equation in Fourier space to get the gravitational potential
  if (head) print*, ''
  if (head) print*, 'Potential field'
  call system_clock(t1,t_rate)
  !$omp paralleldo&
  !$omp& default(shared) &
  !$omp& private(k,j,i,kg,jg,ig,kz,ky,kx,kr)
  do k=1,npen ! construct the potential kernel -1/k
  do j=1,ngic
  do i=1,nyquist+1
    kg=(nn*(icz-1)+icy-1)*npen+k
    jg=(icx-1)*ngic+j
    ig=i
    kz=mod(kg+nyquist-1,nw_global)-nyquist
    ky=mod(jg+nyquist-1,nw_global)-nyquist
    kx=ig-1
    kz=2*sin(pi*kz/nw_global)
    ky=2*sin(pi*ky/nw_global)
    kx=2*sin(pi*kx/nw_global)
    kr=kx**2+ky**2+kz**2
    kr=max(kr,1.0/nw_global**2) ! avoid kr being 0
    cxyz(i,j,k)=-4*pi/kr ! raw potential kernel
  enddo
  enddo
  enddo
  !$omp endparalleldo
  if (head) cxyz(1,1,1)=0 ! zero frequency
  sync all

  if (correct_kernel) then ! correct the potential kernel
    if (head) print*, '  correct kernel'
    call pencil_fft_backward
    temp8=0
    if (image==1) temp8=temp8+rho1(9,1,1)+rho1(1,9,1)+rho1(1,1,9)
    sync all
    if (icx==nn .and. icy==1 .and. icz==1) temp8=temp8+rho1(ngic-7,1,1)
    sync all
    if (icx==1 .and. icy==nn .and. icz==1) temp8=temp8+rho1(1,ngic-7,1)
    sync all
    if (icx==1 .and. icy==1 .and. icz==nn) temp8=temp8+rho1(1,1,ngic-7)
    sync all
    if (head) then
      do i=2,nn**3
        temp8=temp8+temp8[i]
      enddo
    endif
    sync all
    phi8=temp8[1]/6
    if (head) print*,'  phi8 =',phi8
    sync all
    if (head) print*, '  construct Ewald potential kernel in real space'
    !$omp paralleldo default(shared) private(k,j,i,kg,jg,ig,kx,ky,kz,kr)
    do k=1,ngic
    do j=1,ngic
    do i=1,ngic
      kg=k+ngic*(icz-1)
      jg=j+ngic*(icy-1)
      ig=i+ngic*(icx-1)
      kx=mod(kg+nyquist-1,nw_global)-nyquist
      ky=mod(jg+nyquist-1,nw_global)-nyquist
      kz=mod(ig+nyquist-1,nw_global)-nyquist
      kr=sqrt(kx**2+ky**2+kz**2)
      if (kr>8) then
        rho1(i,j,k)=rho1(i,j,k)-(phi8+1/8.)
      elseif (kr>0) then
        rho1(i,j,k)=-1/kr
      elseif (kr==0) then
        rho1(i,j,k)=-2.5
      endif
    enddo
    enddo
    enddo
    !$omp endparalleldo
    sync all
    call pencil_fft_forward
    ! now the potential kernel is in Fourier space
  endif

  cxyz=real(cxyz)*delta_k ! phi(k) = kernel(k) * delta_L(k)
  if (write_phik) then
    if (head) print*, '  write phi1 into file'
    open(11,file=output_name('phik'),status='replace',access='stream')
    write(11) cxyz
    close(11)
  endif
  call pencil_fft_backward ! real space phi is stored in rho1

  phi(1:ngic,1:ngic,1:ngic)=rho1 ! phi1
  if (head .or. image==nn**3) print*,'  phi',phi(1:4,1,1); sync all
  if (write_potential) then
    if (head) print*, '  write phi1 into file'
    open(11,file=output_name('phi1'),status='replace',access='stream')
    write(11) rho1
    close(11)
  endif

  sync all
  ! buffer phi ---------------------------------------------------
  if (head) print*, '  buffer phi'
  phi(:0,:,:)=phi(ngic-nb:ngic,:,:)[inx,icy,icz]; sync all
  phi(ngic+1:,:,:)=phi(1:nb+1,:,:)[ipx,icy,icz]; sync all
  phi(:,:0,:)=phi(:,ngic-nb:ngic,:)[icx,iny,icz]; sync all
  phi(:,ngic+1:,:)=phi(:,1:nb+1,:)[icx,ipy,icz]; sync all
  phi(:,:,:0)=phi(:,:,ngic-nb:ngic)[icx,icy,inz]; sync all
  phi(:,:,ngic+1:)=phi(:,:,1:nb+1)[icx,icy,ipz]; sync all

  if (head) print*, '  destroying FFT plans'
  call destroy_penfft_plan; sync all
  call system_clock(t2,t_rate)
  if (head) print*, '  elapsed time =',real(t2-t1)/t_rate,'secs'; sync all

  ! zip checkpoints ------------------------------------------------
  if (head) print*, ''
  if (head) print*, 'zip checkpoints'
  vf=vfactor(sim%a) ! velocity factor
  if (head) print*, '  vf =',vf
  !! maximum gradient of phi
  grad_max(1)=maxval(abs(phi(-nb:ngic+nb-1,:,:)-phi(-nb+2:ngic+nb+1,:,:)))
  grad_max(2)=maxval(abs(phi(:,-nb:ngic+nb-1,:)-phi(:,-nb+2:ngic+nb+1,:)))
  grad_max(3)=maxval(abs(phi(:,:,-nb:ngic+nb-1)-phi(:,:,-nb+2:ngic+nb+1)))
  sync all
  if (head) then
    do i=2,nn**3 ! reduce
      grad_max(:)=max(grad_max(:),grad_max(:)[i])
    enddo
  endif
  sync all
  grad_max(:)=grad_max(:)[1] ! broadcast
  vmax=grad_max/(8*nic*pi)*vf ! maximum velocity
  sim%dt_vmax=vbuf*16./maxval(abs(vmax)) ! constrain dt by maximum velocity
  sim%vz_max=vmax(3)
  nlayer=2*ceiling(grad_max(3)*np_nc/(8*nic*pi*ratio_cs))+3 ! for OpenMP
  if (head) then
    print*, '  grad_max',grad_max
    print*, '  max dsp',grad_max/(8*nic*pi)
    print*, '  vmax',vmax
    print*, '  vz_max',sim%vz_max
    if (maxval(grad_max)/(8*nic*pi)>=nb) then
      print*, '  particle dsp > buffer' ! particle might move beyond buffer depth
      print*, maxval(grad_max)/(8*nic*pi),nb
      stop
    endif
    print*,'  Thread save nlayer =',nlayer
    print*,''
  endif
  sync all

  ! velocity conversion as a function of redshift and scale
  open(11,file='./velocity_conversion/sigmav_z.bin',access='stream')
  read(11) svz
  close(11)
  open(11,file='./velocity_conversion/sigmav_r.bin',access='stream')
  read(11) svr
  close(11)

  sigma_vf=interp_sigmav(sim%a,box/ng_global) ! sigma(v) on scale of fine grid, in km/s
  sigma_vc=interp_sigmav(sim%a,box/nc_global) ! sigma(v) on scale of coarse grid, in km/s
  sim%sigma_vres=sqrt(sigma_vf**2-sigma_vc**2) ! sigma(v) residual, in km/s
  sim%sigma_vi=sim%sigma_vres/sim%vsim2phys/sqrt(3.) ! sigma(v_i) residual, in sim unit
  sync all
  if (head) then
    print*, ''
    print*, 'Read velocity dispersion prediction'
    print*,'sigma_vf(a=',sim%a,', r=',box/ng_global,'Mpc/h)=',real(sigma_vf,4),'km/s'
    print*,'sigma_vc(a=',sim%a,', r=',box/nc_global,'Mpc/h)=',real(sigma_vc,4),'km/s'
    print*,'sigma_vres=',real(sim%sigma_vres,4),'km/s'
    print*,'sigma_vi =',real(sim%sigma_vi,4),'(simulation unit)'
  endif

  sync all
  ! create particles (no communication) ----------------------------
  if (head) print*,''
  if (head) print*, 'Create particles'
  diff_x(1,:,:)=-0.5; diff_x(2,:,:)=0.5
  diff_y(:,1,:)=-0.5; diff_y(:,2,:)=0.5
  diff_z(:,:,1)=-0.5; diff_z(:,:,2)=0.5
  call system_clock(t1,t_rate)
  open(11,file=output_name('xp'),status='replace',access='stream') ! position list
  open(12,file=output_name('vp'),status='replace',access='stream') ! velocitie list
  open(13,file=output_name('np'),status='replace',access='stream') ! particle number list (integer density on coarse grid)
  !open(14,file=output_name('vc'),status='replace',access='stream') ! velocity_coarse-grid list
#ifdef PID
  if (head) print*, '  also create PID'
  open(15,file=output_name('id'),status='replace',access='stream') ! particle ID list
#endif
  if (body_centered_cubic .and. ratio_cs/np_nc/2==0) stop 'ratio_cs/np_nc/2 = 0, unsuitable for body centered cubic'
  !vfield=0
  std_vsim_c=0; std_vsim_res=0; std_vsim=0; !np_prev=0
  do itz=1,nnt ! loop tile by tile
  do ity=1,nnt
  do itx=1,nnt
    !print*,'tile',itx,ity,itz
    rhoce=0; rholocal=0
    ! first loop to create integer number density field "rhoce", and coarse velocity field "vfield"
    do ilayer=0,nlayer-1
      !$omp paralleldo default(shared) num_threads(ncore) schedule(dynamic,1) private(k,j,i,imove,kk,jj,ii,xq,gradphi,g,vreal)
      do k=1-npb+ilayer,npt+npb,nlayer
      do j=1-npb,npt+npb
      do i=1-npb,npt+npb
      do imove=0,merge(1,0,body_centered_cubic)
        kk=ngp*nic*(itz-1)+(ratio_cs*nic/np_nc)*(k-1)+1+imove*(ratio_cs*nic/np_nc/2)
        jj=ngp*nic*(ity-1)+(ratio_cs*nic/np_nc)*(j-1)+1+imove*(ratio_cs*nic/np_nc/2)
        ii=ngp*nic*(itx-1)+(ratio_cs*nic/np_nc)*(i-1)+1+imove*(ratio_cs*nic/np_nc/2)
        xq=([i,j,k]-1d0)/np_nc + (0.5d0/nic+imove*(ratio_cs/np_nc/2))/ratio_cs ! Lagrangian position q, in coarse grid
        !gradphi(1)=phi(ii+1,jj,kk)-phi(ii-1,jj,kk)
        !gradphi(2)=phi(ii,jj+1,kk)-phi(ii,jj-1,kk)
        !gradphi(3)=phi(ii,jj,kk+1)-phi(ii,jj,kk-1)
        
        !xq=([i,j,k]-1d0)/np_nc + (imove*(ratio_cs/np_nc/2))/ratio_cs
        gradphi(1)=sum(phi(ii-1:ii,jj-1:jj,kk-1:kk)*diff_x)
        gradphi(2)=sum(phi(ii-1:ii,jj-1:jj,kk-1:kk)*diff_y)
        gradphi(3)=sum(phi(ii-1:ii,jj-1:jj,kk-1:kk)*diff_z)
        
        g=ceiling(xq-gradphi/(8*nic*pi*ratio_cs)) ! g(1:3) is the grid index
        rhoce(g(1),g(2),g(3))=rhoce(g(1),g(2),g(3))+1 ! increment in an extended density field
        !print*, xq,ii,jj,kk; pause
        !vreal=-gradphi/(8*nic*pi)*vf
        !vfield(:,g(1),g(2),g(3))=vfield(:,g(1),g(2),g(3))+vreal ! record vfield according to real particles
      enddo ! imove
      enddo
      enddo
      enddo !k
      !$omp endparalleldo
    enddo ! ilayer
    !vfield(1,:,:,:)=vfield(1,:,:,:)/merge(1,rhoce,rhoce==0)
    !vfield(2,:,:,:)=vfield(2,:,:,:)/merge(1,rhoce,rhoce==0)
    !vfield(3,:,:,:)=vfield(3,:,:,:)/merge(1,rhoce,rhoce==0)
    !print*,'a', sum(rhoce),npb,npt
    call spine_tile(rhoce,idx_ex_r,pp_l,pp_r,ppe_l,ppe_r) ! spine structure to speed up indexing
    !print*,'b'
    ! second loop to create particles, according to rhoce
    ! also assign vp according to vfield
    do ilayer=0,nlayer-1
      !$omp paralleldo default(shared) num_threads(ncore) schedule(dynamic,1)&
      !$omp& private(k,j,i,imove,kk,jj,ii,xq,gradphi,g,idx,vreal,iq)
      do k=1-npb+ilayer,npt+npb,nlayer
      do j=1-npb,npt+npb
      do i=1-npb,npt+npb
      do imove=0,merge(1,0,body_centered_cubic)
        kk=ngp*nic*(itz-1)+(ratio_cs*nic/np_nc)*(k-1)+1+imove*(ratio_cs*nic/np_nc/2)
        jj=ngp*nic*(ity-1)+(ratio_cs*nic/np_nc)*(j-1)+1+imove*(ratio_cs*nic/np_nc/2)
        ii=ngp*nic*(itx-1)+(ratio_cs*nic/np_nc)*(i-1)+1+imove*(ratio_cs*nic/np_nc/2)
        xq=([i,j,k]-1d0)/np_nc + (0.5d0/nic+imove*(ratio_cs/np_nc/2))/ratio_cs
        !gradphi(1)=phi(ii+1,jj,kk)-phi(ii-1,jj,kk)
        !gradphi(2)=phi(ii,jj+1,kk)-phi(ii,jj-1,kk)
        !gradphi(3)=phi(ii,jj,kk+1)-phi(ii,jj,kk-1)
        
        !xq=([i,j,k]-1d0)/np_nc + (imove*(ratio_cs/np_nc/2))/ratio_cs
        gradphi(1)=sum(phi(ii-1:ii,jj-1:jj,kk-1:kk)*diff_x)
        gradphi(2)=sum(phi(ii-1:ii,jj-1:jj,kk-1:kk)*diff_y)
        gradphi(3)=sum(phi(ii-1:ii,jj-1:jj,kk-1:kk)*diff_z)
        
        g=ceiling(xq-gradphi/(8*nic*pi*ratio_cs))
        !print*, g
        rholocal(g(1),g(2),g(3))=rholocal(g(1),g(2),g(3))+1 ! number of particles have been written in the current coarse cell
        idx=idx_ex_r(g(2),g(3))-sum(rhoce(g(1):,g(2),g(3)))+rholocal(g(1),g(2),g(3))
        xp(:,idx)=floor((xq-gradphi/(8*nic*pi*ratio_cs))/x_resolution,kind=8)
        
        vp(:,idx)=-gradphi/(8*pi*nic)*vf
        
        !vreal=-gradphi/(8*pi*nic)*vf
        !vreal=vreal-vfield(:,g(1),g(2),g(3)) ! save relative velocity
        !vp(:,idx)=nint(real(nvbin-1)*atan(sqrt(pi/2)/(sim%sigma_vi*vrel_boost)*vreal)/pi,kind=izipv)
#       ifdef PID
          !! particle IDs are unique according to their Lagrangian positions, not necessarily consecutive
          iq = ([icx,icy,icz]-1)*ng + ([itx,ity,itz]-1)*ngp + (ratio_cs/np_nc)*([i,j,k]-1)+imove
          iq = modulo(iq,ng_global)
          pid(idx)=iq(1)+ng_global*iq(2)+ng_global**2*iq(3)+1
#       endif
      enddo ! imove
      enddo
      enddo
      enddo ! k
      !$omp endparalleldo
    enddo ! ilayer

    do k=1,nt ! delete buffer particles
    do j=1,nt
      xp(:,pp_l(j,k):pp_r(j,k))=xp(:,ppe_l(j,k):ppe_r(j,k))
      vp(:,pp_l(j,k):pp_r(j,k))=vp(:,ppe_l(j,k):ppe_r(j,k))
#     ifdef PID
        pid(pp_l(j,k):pp_r(j,k))=pid(ppe_l(j,k):ppe_r(j,k))
#     endif
    enddo
    enddo

    ! velocity analysis
    !$omp paralleldo&
    !$omp& default(shared) &
    !$omp& private(k,j,i,nzero,l,ip,vreal)&
    !$omp& reduction(+:std_vsim_c,std_vsim_res,std_vsim)
    do k=1,nt
    do j=1,nt
    do i=1,nt
      nzero=pp_r(j,k)-sum(rhoce(i:nt,j,k))
      !std_vsim_c=std_vsim_c+sum(vfield(:,i,j,k)**2)
      do l=1,rhoce(i,j,k)
        ip=nzero+l
        !vreal=tan(pi*real(vp(:,ip))/real(nvbin-1))/(sqrt(pi/2)/(sim%sigma_vi*vrel_boost))
        !std_vsim_res=std_vsim_res+sum(vreal**2)
        !vreal=vreal+vfield(:,i,j,k)
        !std_vsim=std_vsim+sum(vreal**2)
      enddo
    enddo
    enddo
    enddo
    !$omp endparalleldo

    write(11) xp(:,1:pp_r(nt,nt))
    write(12) vp(:,1:pp_r(nt,nt))
    write(13) rhoce(1:nt,1:nt,1:nt)
    !write(14) vfield(:,1:nt,1:nt,1:nt)
#   ifdef PID
      write(15) pid(1:pp_r(nt,nt))
#   endif
    sim%nplocal=sim%nplocal+pp_r(nt,nt)
  enddo
  enddo
  enddo ! end of tile loop

  close(11)
  close(12)
  close(13)
  !close(14)
# ifdef PID
    close(15)
# endif
  call system_clock(t2,t_rate)
  sync all

  if (head) print*, '  elapsed time =',real(t2-t1)/t_rate,'secs';

  if (head) then
    !print*,''
    !print*,'Velocity analysis on head node'
    !std_vsim_res=sqrt(std_vsim_res/sim%nplocal)
    !std_vsim_c=sqrt(std_vsim_c/nc/nc/nc)
    !std_vsim=sqrt(std_vsim/sim%nplocal)
    !print*,'  std_vsim         ',real(std_vsim*sim%vsim2phys,4),'km/s'
    !print*,'  std_vsim_c       ',real(std_vsim_c*sim%vsim2phys,4),'km/s'
    !print*,'  std_vsim_res     ',real(std_vsim_res*sim%vsim2phys,4),'km/s'
    !print*,'  std_vi (sim unit)',real(std_vsim_res/sqrt(3.),4),'(simulation unit)'
    !print*,''
  endif
  sync all

  if (head .or. this_image()==num_images()) print*,'image',image,', nplocal',sim%nplocal
  sim%npglobal=0
  if (head) then
    do i=1,nn**3
      sim%npglobal=sim%npglobal+sim[i]%nplocal
    enddo
  endif
  sync all
  sim%npglobal=sim[1]%npglobal
  if (head) then
    print*, 'npglobal =',sim%npglobal
    if (sim%npglobal/=merge(2,1,body_centered_cubic)*(np_nc*nc*nn)**3) then
      print*, 'warning: incorrect npglobal'
    endif
  endif
  sim%mass_p_cdm=real(ng_global**3,kind=8)/sim%npglobal ! particle mass in unit of fine cell
  if (head) print*,'sim%mass_p_cdm =',sim%mass_p_cdm
  sync all

  call print_header(sim)

  sync all
  open(10,file=output_name('info'),status='replace',access='stream')
  write(10) sim
  close(10)
  sync all

  call system_clock(ttt2,t_rate)
  if (head) print*, 'total elapsed time =',real(ttt2-ttt1)/t_rate,'secs';
  if (head) print*, 'initial condition done'

  contains

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
  endfunction

  real function interp_tf(kr,ix,iy) ! interpolation in log space
    implicit none
    integer(4) ix,iy
    integer(8) ii,i1,i2
    real kr,xx,yy,x1,x2,y1,y2
    i1=1
    i2=nk
    do while (i2-i1>1)
      ii=(i1+i2)/2
      if (kr>tf(ix,ii)) then
        i1=ii
      else
        i2=ii
      endif
    enddo
    x1=log(tf(ix,i1))
    y1=log(tf(iy,i1))
    x2=log(tf(ix,i2))
    y2=log(tf(iy,i2))
    xx=log(kr)
    yy=y1+(y2-y1)*(xx-x1)/(x2-x1)
    interp_tf=exp(yy)
  endfunction


  function tophat(x)
    implicit none
    real :: x,tophat
    !if (x/=0) then
    !  tophat=3*(sin(x)-cos(x)*x)/x**3
    !else
    !  tophat=1
    !endif
    tophat=merge(1d0,3*(sin(x*1d0)-cos(x*1d0)*x*1d0)/(x*1d0)**3,x==0)
  endfunction tophat

  function Dgrow(a) ! growth function
    implicit none
    real, parameter :: om=omega_m
    real, parameter :: ol=omega_l
    real a
    real Dgrow
    real g,ga,hsq,oma,ola
    hsq=om/a**3+(1-om-ol)/a**2+ol
    oma=om/(a**3*hsq)
    ola=ol/hsq
    g=2.5*om/(om**(4./7)-ol+(1+om/2)*(1+ol/70))
    ga=2.5*oma/(oma**(4./7)-ola+(1+oma/2)*(1+ola/70))
    Dgrow=a*ga/g
  end function Dgrow

  real function vfactor(a)
    implicit none
    real :: a
    real :: H,km,lm
    lm=omega_l/omega_m
    km=(1-omega_m-omega_l)/omega_m
    H=2/(3*sqrt(a**3))*sqrt(1+a*km+a**3*lm)
    vfactor=a**2*H
  endfunction

  function lcg(s) ! Linear congruential generator
    implicit none
    integer(4) :: lcg
    integer(int64) :: s
    if (s == 0) then
       s = 104729
    else
       s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  endfunction

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
end
