subroutine particle_initialization
  use parameters
#ifdef PID
  use variables, only: xp,xp_new,vp,vp_new,pid,pid_new
  use variables, only: rhoc,sigma_vi,np_image_max,np_tile_max
#else
  use variables, only: xp,xp_new,vp,vp_new
  use variables, only: rhoc,sigma_vi,np_image_max,np_tile_max
#endif
  implicit none
  
  character(100) fn10,fn11,fn12,fn13,fn14,fn15
  call tic(41)
  if (head) then
    print*, ''
    print*, 'particle_initialization'
    print*, '  at redshift',z_checkpoint(sim%cur_checkpoint)
  endif
  fn10=output_name('info')
  fn11=output_name('xp')
  fn12=output_name('vp')
  fn13=output_name('np')
  fn14=output_name('vc')
  fn15=output_name('id')

  open(10,file=fn10,status='old',access='stream')
  read(10) sim
  close(10)
  sim%a=1./(1+z_checkpoint(sim%cur_checkpoint))
  sigma_vi=sim%sigma_vi

  if (sim%izipx/=izipx .or. sim%izipv/=izipv) then
    print*, '  zip format incompatable'
    close(12)
    stop
  endif

  allocate(xp(3,np_image_max)[nn,nn,*],xp_new(3,np_tile_max))
  allocate(vp(3,np_image_max)[nn,nn,*],vp_new(3,np_tile_max))
  allocate(rhoc(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)[nn,nn,*])
  !allocate(vfield(3,1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)[nn,nn,*])
#ifdef PID
  allocate(pid(np_image_max)[nn,nn,*],pid_new(np_tile_max))
#endif
  !$omp parallelsections default(shared)
  !$omp section
    !omp workshare
    open(11,file=fn11,status='old',access='stream'); read(11) xp(:,:sim%nplocal); close(11)
    !omp endworkshare
  !$omp section
    open(12,file=fn12,status='old',access='stream'); read(12) vp(:,:sim%nplocal); close(12)
  !$omp section
    open(13,file=fn13,status='old',access='stream'); read(13) rhoc(1:nt,1:nt,1:nt,:,:,:); close(13)
!!  !$omp section
!!    open(14,file=fn14,status='old',access='stream'); read(14) vfield(:,1:nt,1:nt,1:nt,:,:,:); close(14)
# ifdef PID
  !$omp section
    open(15,file=fn15,status='old',access='stream'); read(15) pid(:sim%nplocal); close(15)
    !print*, '  image',this_image(),' PID range: ',minval(pid(:sim%nplocal)),maxval(pid(:sim%nplocal))
    if (minval(pid(:sim%nplocal))<1) then
      print*, 'warning: pid are not all positive'
    endif
# endif
  if (this_image()==1 .or. this_image()==num_images()) print*,'  from image',this_image(),' read',sim%nplocal,' particles'
  !$omp endparallelsections
  call toc(41)
  if (head) then
    print*,'  npglobal    =', sim%npglobal
    print*,'  mass_p_cdm =', sim%mass_p_cdm
    print*,'  vsim2phys =',sim%vsim2phys, ' (km/s)/(1.0)'
    print*,'  sigma_vi    =',sigma_vi,'(simulation unit)'
    print*,'  elapsed time =',tcat(6,0),'secs'
    print*,''
  endif
  sync all
endsubroutine
