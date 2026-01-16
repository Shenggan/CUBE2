#define write_gadget_snapshot

program convert_format
  use parameters
  implicit none
  integer(8) i,j,k,l,nplocal,npglobal,itx,ity,itz,nlast,ip,np,idx1(3),idx2(3)
  integer cur_checkpoint,itile,iq
  real mass_p,pos1(3),dx1(3),dx2(3)
  integer(izipx),allocatable :: xp(:,:)
  integer(izipv),allocatable :: vp(:,:)
  integer(izipi),allocatable :: pid(:)
  integer(4),allocatable :: rhoc(:,:,:,:,:,:)
  real(4),allocatable :: xv(:,:),vc(:,:,:,:,:,:,:)
  character(20) str_z,str_i
  character(6) str_snap,str_rank
#ifdef write_gadget_snapshot
  type gadget_snapshot_header
    integer(4) Npart(6)
    real(8) Massarr(6),Time,Redshift
    integer(4) Flagsfr,Flagfbk,Nall(6),FlagCooling,NumFiles
    real(8) BoxSize,Omega0,OmegaLambda,HubbleParam
    integer(4) FlagAge,FlagMetals,NallHW(6),flag_entr_ics,unused(15)
  endtype
  type(gadget_snapshot_header) io_header
#endif

  call geometry
  if (head) then
    print*,''; print*, 'Convert to gadget format. Checkpoint at:'
    open(16,file='z_checkpoint.txt',status='old')
    do i=1,nmax_redshift
      read(16,end=71,fmt='(f8.4)') z_checkpoint(i); print*, i, '   z=',z_checkpoint(i)
    enddo
    71 n_checkpoint=i-1; close(16); print*,''
  endif
  allocate(rhoc(nt,nt,nt,nnt,nnt,nnt),vc(3,nt,nt,nt,nnt,nnt,nnt)); sync all
  n_checkpoint=n_checkpoint[1]; z_checkpoint(:)=z_checkpoint(:)[1]

  do cur_checkpoint =n_checkpoint,n_checkpoint
    write(str_snap,'(i6)') cur_checkpoint+100000
    !write(str_rank,'(i6)') rank
    !if (head) print*, '  will write into file:',&
    !  output_dir()//'snapshot_'//str_snap(4:6)//'.'//trim(adjustl(str_rank))
    sim%cur_checkpoint=cur_checkpoint
    if (head) print*, 'Start analyzing redshift ',z2str(z_checkpoint(cur_checkpoint))
    if (head) print*, output_name('info')
    open(11,file=output_name('info'),access='stream'); read(11) sim; close(11)
    if (sim%izipx/=izipx .or. sim%izipv/=izipv) stop 'zip format incompatable'
    mass_p=1.
    nplocal=sim%nplocal; npglobal=sim%npglobal
    if (head) then
      print*, 'mass_p =',mass_p
      print*, 'nplocal =',nplocal
      print*, 'npglobal =',npglobal
      print*, 'sigma_vi =',sim%sigma_vi
      print*, 'vsim2phys =',sim%vsim2phys
    endif
    allocate(xp(3,nplocal),vp(3,nplocal),xv(6,nplocal),pid(nplocal))
    open(11,file=output_name('xp'),access='stream'); read(11) xp; close(11)
    open(11,file=output_name('vp'),access='stream'); read(11) vp; close(11)
    open(11,file=output_name('np'),access='stream'); read(11) rhoc; close(11)
    open(11,file=output_name('vc'),access='stream'); read(11) vc; close(11)
    open(11,file=output_name('id'),access='stream'); read(11) pid; close(11)
    if (head) print*,'check PID range:',minval(pid),maxval(pid)
    !deallocate(xp,vp,xv,pid); cycle
    nlast=0

    itile=0

    do itz=1,nnt
    do ity=1,nnt
    do itx=1,nnt

      itile=itile+1
      write(str_rank,'(i6)') itile-1
      print*, '  will write into file:',&
        output_dir()//'snapshot_'//str_snap(4:6)//'.'//trim(adjustl(str_rank))

      iq=0
      do k=1,nt
      do j=1,nt
      do i=1,nt
        np=rhoc(i,j,k,itx,ity,itz)
        do l=1,np
          ip=nlast+l; iq=iq+1
          pos1=nt*([itx,ity,itz]-1)+ ([i,j,k]-1) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
          xv(1:3,iq)=pos1*real(ng)/real(nc)+([m1,m2,m3]-1)*ng
          xv(4:6,iq)=vc(:,i,j,k,itx,ity,itz)+tan((pi*real(vp(:,ip)))/real(nvbin-1))/(sqrt(pi/2)/(sim%sigma_vi*vrel_boost))*sim%vsim2phys
        enddo
        nlast=nlast+np
      enddo
      enddo
      enddo

      if (head) print*,'write gadget snapshot'
      if (head) print*,output_dir()//'snapshot_'//str_snap(4:6)//'.'//trim(adjustl(str_rank))
      io_header%Npart=[0,int(iq,kind=4),0,0,0,0]
      io_header%Massarr=[0.,27.75*omega_m*(box**3)/npglobal,0.,0.,0.,0.] 
      io_header%Time=sim%a; io_header%Redshift=1./sim%a-1
      io_header%FlagSfr=0; io_header%Flagfbk=0
      io_header%Nall=[0,int(npglobal,kind=4),0,0,0,0]
      io_header%FlagCooling=0
      io_header%NumFiles=nn**3; io_header%BoxSize=box*1000
      io_header%Omega0=omega_m; io_header%OmegaLambda=omega_l
      io_header%HubbleParam=h0
      io_header%FlagAge=0; io_header%FlagMetals=0
      io_header%NallHW=0; io_header%flag_entr_ics=0; io_header%unused=0
      xv(1:3,1:iq)=xv(1:3,:)*box*1000/real(ng) ! positions in kpc/h, real(4)
      xv(4:6,1:iq)=xv(4:6,:)/sqrt(sim%a) ! velocities in km/sec/sqrt(a), real(4)
      open(11,file=output_dir()//'snapshot_'//str_snap(4:6)//'.'//trim(adjustl(str_rank)),status='replace',form='unformatted')
        write(11) io_header ! 1st block: 256 byte header
        write(11) xv(1:3,1:iq) ! 2nd block: positions
        write(11) xv(4:6,1:iq) ! 3rd block: velocities
        write(11) int(pid(nlast-iq+1:nlast),kind=8) ! 4th block: PIDs, in integer(4)
      close(11)

      sync all
      if (head) then 
        print*,'io_header'; print*,io_header
        print*,'position '; print*,xv(1:3,1),xv(1:3,nplocal)
        print*,'velocity '; print*,xv(4:6,1),xv(4:6,nplocal)
        print*,'PID      '; print*,pid(1),pid(nplocal)
        print*,'write_gadget_snapshot done.'
      endif
    enddo
    enddo
    enddo
    deallocate(xv,pid,xp,vp)
  enddo
  deallocate(rhoc,vc)
end
