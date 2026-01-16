module parameters
  implicit none
  
  ! output directory
  character(*),parameter :: opath='/mnt/18T/cube/usr_out/cbh/1000_128/'

  ! zip parameters
  integer,parameter :: ndim=3
  integer,parameter :: izipx=2 ! size to store xp as
  integer,parameter :: izipv=2 ! size to store vp as
  integer, parameter :: izipi = 8 ! if pids are on, size to store as

  ! cell resolution parameters
  integer,parameter :: nn=1
  real,parameter :: box=100
  integer,parameter :: ncore=32
  integer,parameter :: nteam=8
  integer,parameter :: nnest=4
  integer(8),parameter :: ng=512

  integer,parameter :: ratio_cs=4
  integer,dimension(7),parameter :: ratio_sf=[1,2,4,6,8,12,16]
  integer,dimension(7),parameter :: ratio_cf=ratio_cs*ratio_sf ! 8,16,24,32
  integer,parameter :: nnt=4 ! number of tiles /image/dim
  integer,parameter :: nns=4 ! number of subtiles /tile/dim

  ! particle resolution parameters
  integer(8),parameter :: np_nc=ratio_cs
  integer,parameter :: nic=1   ! refined resolution for IC
  integer,parameter :: ngic=ng*nic
  logical,parameter :: body_centered_cubic=.false.

  ! force resolution parameters
  real,parameter :: apm1c=3.5            ! PM1's softening
  real,parameter :: apm1=apm1c*ratio_cs  ! PM2's range = 14
  real,parameter :: apm2=3.5             ! PM2's softening = 3.5
  real,dimension(7),parameter :: apm2f=apm2*ratio_sf  ! PM3's range =(3.5),7,14,21,28
  real,parameter :: apm3f=3.5            ! PM3's softening = 3.5
  real,dimension(7),parameter :: apm3=3.5/ratio_sf    ! PP's range
  real,dimension(7),parameter :: appr=apm3
  real,parameter :: dt_refine=1 ! rescale dt
  real,parameter :: app=0.06 ! target=0.05 !apm3f/ratio_sf ! 0.875 pp softening: apm3f/ratio_sf
  integer,parameter :: nrange=1 ! 1 actual PP computing range
  integer,parameter :: n_neighbor=((1+2*nrange)**3-1)/2 ! 13

  ! Green's function
  integer,parameter :: n_int=2
  integer,parameter :: p=2 ! order of interpolation. 1=CIC, 2=TSC
  real,parameter :: alpha=4./3
  real,dimension(4),parameter :: weight = [(alpha-1)/4,-alpha/2,alpha/2,(1-alpha)/4]

  ! derived parameters
  integer(8),parameter :: nc=ng/ratio_cs ! 192 coarse cells /image/dim, must be integer
  integer(8),parameter :: nt=nc/nnt ! 96
  integer,parameter :: ntt=nt/nns ! 24
  integer(8),parameter :: ngp=ng/nnt ! 384 physical tile
  integer,parameter :: ngb=16 ! tile buffer, floor(apm1)+2
  integer(8),parameter :: ngt=ngp+2*ngb ! 416
  integer(8),parameter :: ng_global=ng*nn
  integer(8),parameter :: nc_global=nc*nn

  integer(8),dimension(7),parameter :: nfp=(ngp/nns)*ratio_sf ! 192*[1,2,4,6,8,12,16]
  !integer(8),parameter :: nf=nfp
  integer,dimension(7),parameter :: nfb=ratio_cf
  integer(8),dimension(7),parameter :: nft=nfp+2*nfb
  
  ! ngrid /image/dim for pencil-fft
# ifdef FFTFINE
    integer(8),parameter :: nw=ngic ! standard grid fft, for IC, dsp, convert_zip_to_xv
# else
    integer(8),parameter :: nw=nc ! coarse grid fft, for main code
# endif
  integer(8),parameter :: npen=nw/nn ! ng /dim in shorter side of the pencil, for pencil decomposition
  integer(8),parameter :: nw_global=nw*nn
  integer(8),parameter :: nyquist=nw_global/2

  integer(8),parameter :: ncb=ngb/ratio_cs ! nc in buffer /dim, single side; 6 by default
  integer(8),parameter :: nte=nt+2*ncb ! extended nc

  
  real,parameter :: image_buffer=2.0
  real,parameter :: tile_buffer=3.0
  real,parameter :: vbuf=0.9

  ! derived parameters
  integer(8),parameter :: nvbin=int(2,8)**(8*izipv)
  integer(8),parameter :: ishift=-(int(2,8)**(izipx*8-1))
  real(8),parameter :: rshift=0.5-ishift

  real(8),parameter :: pi=4*atan(1d0)

  ! cosmological parameters
  real,parameter :: s8=0.821131408 ! use -Dsigma_8 in initial_conditions
  integer,parameter :: zdim=2 ! the dimension being the redshift direction

  ! background parameters
  real, parameter :: h0=0.70

  real, parameter :: omega_cdm=0.2327 ! cdm energy
  real, parameter :: omega_bar=0.0463 ! baryon energy, goes into cdm
  real, parameter :: omega_mhd=0.0 ! mhd energy, evolved separately
  real, parameter :: omega_m=omega_cdm+omega_bar+omega_mhd ! total matter

  !real, parameter :: omega_l = 1-omega_m-omega_r
  real, parameter :: omega_l=1-omega_m
  real, parameter :: wde=-1 ! de equation of state

  ! initial conditions
  real,parameter :: f_nl=0
  real,parameter :: g_nl=0
  real,parameter :: n_s=0.972
  real,parameter :: A_s=2.6025e-09
  real,parameter :: k_o=0.05/h0

  integer(8),parameter :: istep_max=250000
  real,parameter :: ra_max=0.1
  real(8),parameter :: v_resolution=2.1/(int(2,8)**(izipv*8))
  real(8),parameter :: x_resolution=1.0/(int(2,8)**(izipx*8))
  !real(8),parameter :: vdisp_boost=1.0
  real(8),parameter :: vrel_boost=2.5
  real,parameter :: b_link=0.20 ! linking length for FoF
  real,parameter :: np_halo_min=10 ! minimum number of particles to be a halo

  !! MPI image variables !!
  integer(8) image,rank,icx,icy,icz,inx,iny,inz,ipx,ipy,ipz
  integer(8) m1,m2,m3,m
  logical head
  ! checkpoint variables
  integer(8),parameter :: nmax_redshift=1000
  integer(8) n_checkpoint[*],n_halofind[*]
  real z_checkpoint(nmax_redshift)[*],z_halofind(nmax_redshift)[*]
  logical checkpoint_step[*],halofind_step[*],final_step[*]

  ! timing variables
  integer(4) istep,t1,t2,tt1,tt2,ttt1,ttt2,t_start,t_end,t_rate,tictoc(2,0:100),tall_start,tall_end,tall_rate
  real(4) tcat(100,0:10000)

  type sim_header
    integer(8) nplocal,npglobal
    integer(8) izipx,izipv
    integer(8) image
    integer(8) nn,nnt,nt,ncell,ncb
    integer(8) timestep
    integer(8) cur_checkpoint
    integer(8) cur_halofind
    real a, t, tau
    real dt_pp, dt_pm3, dt_pm2, dt_pm1, dt_vmax
    real mass_p_cdm
    real box
    real h0
    real omega_m
    real omega_l
    real s8
    real vsim2phys
    real sigma_vres
    real sigma_vi
    real z_i
    real vz_max
  endtype
  type(sim_header) sim[*]


  type type_halo_catalog_header
    integer nhalo_tot,nhalo,ninfo
    real linking_parameter
  endtype
  type type_halo_catalog_array
    real hmass ! number of particles ! 1:1
    real xv(6)
  endtype
  integer,parameter :: ninfo=7 ! number of real numbers per halo in the halo catalog
  character(4) b_link_string


  contains
    subroutine print_header(s)
      type(sim_header),intent(in) :: s
      if (this_image()==1) then
      print*,'-------------------------------- CUBE info --------------------------------'
      print*,'| np local/global =',s%nplocal,s%npglobal
      print*,'| a,t,tau         =',s%a,s%t,s%tau
      print*,'| timestep        =',s%timestep
      print*,'| dt PM123,PP     =',s%dt_pm1,s%dt_pm2,s%dt_pm3,s%dt_pp
      print*,'| dt v            =',s%dt_vmax
      print*,'| cur_checkpoint  =',int(s%cur_checkpoint,2)
      print*,'| cur_halofind    =',int(s%cur_halofind,2)
      print*,'| mass_p          =',s%mass_p_cdm
      print*,'| box             =',s%box, 'Mpc/h'
      print*,'| image           =',s%image
      print*,'| nn              =',s%nn
      print*,'| nnt             =',s%nnt
      print*,'| ncb             =',s%ncb
      print*,'| izip x,v        =',int(s%izipx,1),int(s%izipv,1)
      print*,'| h_0             =',s%h0,'*100 km/s/Mpc'
      print*,'| omega_m         =',s%omega_m
      print*,'| omega_l         =',s%omega_l
      print*,'| sigma_8         =',s%s8
      print*,'| vsim2phys       =',s%vsim2phys, '(km/s)/(1.0)'
      print*,'| sigma_vres      =',s%sigma_vres,'(km/s)'
      print*,'| sigma_vi        =',s%sigma_vi,'(simulation unit)'
      print*,'| z_i             =',s%z_i
      print*,'| vz_max          =',s%vz_max
      print*,'------------------------------------------------------------------------------'
      endif
      sync all
    endsubroutine

    include 'basic_functions.f08'
endmodule
