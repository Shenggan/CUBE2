!#define gadget

program CUBE_FoF
  use variables
  implicit none
  integer i,j,k,l,cur_checkpoint,itile,np,nptile,nfof,n_refine,n1,n2,idx(3),nh_tile,nh[*]
  integer iq1,iq2,iq3,i_neighbor,jq(3)
  integer(8) nlast,nzero,ip,jp,np_iso,np_mem,np_head,nlist
  integer(4),allocatable :: np_halo_all(:),np_halo(:)
  integer(izipi),allocatable :: hoc(:,:,:),ll(:),llgp(:),hcgp(:),ecgp(:),iph_halo_all(:),iph_halo(:)
  real rp2,rsq,xf_hoc(3)
  real,allocatable :: xv(:,:),dxv(:,:),xv_mean(:,:)
  type(type_halo_catalog_header) halo_header
  type(type_halo_catalog_array),allocatable :: hcat(:)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call geometry
  if (head) then
    open(16,file='z_checkpoint.txt',status='old')
    do i=1,nmax_redshift-1
      read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
    enddo
    71 n_checkpoint=i-1
    close(16)
    if (n_checkpoint==0) stop 'z_checkpoint.txt empty'
  endif
  sync all
  n_checkpoint=n_checkpoint[1]
  z_checkpoint(:)=z_checkpoint(:)[1]

  if (head) print*,'  initialize FoF cell neighbors'
  l=0
  do j=-1,1
  do i=-1,1
    l=l+1; ijk(:,l)=[i,j,-1]
  enddo
  enddo
  do i=-1,1
    l=l+1; ijk(:,l)=[i,-1,0]
  enddo
  l=l+1; ijk(:,l)=[-1,0,0]
  !do i=1,13
  !  print*,ijk(:,i)
  !enddo

  if (head) print*,'  initialize tile index'
  ipm2=0;
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
    ipm2=ipm2+1
    ixyz2(:,ipm2)=[itx,ity,itz]
  enddo
  enddo
  enddo

  n_refine=1
  nfof=nt*n_refine*ratio_cs
  n1=1-ncb*n_refine*ratio_cs; n2=(nt+ncb)*n_refine*ratio_cs
  rp2=b_link**2

  do cur_checkpoint=100,100
    sim%cur_checkpoint=cur_checkpoint
    if (head) print*, ''
    if (head) print*, 'FoF at redshift ',z2str(z_checkpoint(cur_checkpoint))
    if (head) print*, '  read checkpoint header',output_name('info')
    call particle_initialization
    deallocate(xp_new,vp_new,pid_new)
    call buffer_grid
    call buffer_x
    call buffer_v

    halo_header%nhalo_tot=0; nh=0
    halo_header%nhalo=0; halo_header%ninfo=ninfo; halo_header%linking_parameter=b_link
    if (head) print*, output_name('halo')
    open(21,file=output_name('halo'),status='replace',access='stream')
    write(21) halo_header
    
    do itile=1,nnt**3 ! work on each tile
      nlast=0
      if (head) print*, ixyz2(:,itile)
      nptile=sum(rhoc(:,:,:,ixyz2(1,itile),ixyz2(2,itile),ixyz2(3,itile)))
      !print*,nfof,n1,n2,nptile
      allocate(xv(6,nptile))
      do k=1-ncb,nt+ncb
      do j=1-ncb,nt+ncb
      do i=1-ncb,nt+ncb
        np=rhoc(i,j,k,ixyz2(1,itile),ixyz2(2,itile),ixyz2(3,itile))
        nzero=idx_b_r(j,k,ixyz2(1,itile),ixyz2(2,itile),ixyz2(3,itile)) &
              -sum(rhoc(i:,j,k,ixyz2(1,itile),ixyz2(2,itile),ixyz2(3,itile)))
        do l=1,np
          ip=nzero+l; jp=nlast+l
          xv(1:3,jp)=ratio_cs*([i,j,k]-1+(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution) ! coarse grid
          xv(4:6,jp)=vp(:,ip)
        enddo
        nlast=nlast+np
      enddo
      enddo
      enddo
      !print*,minval(xv),maxval(xv),nlast;

      ! create hoc ll
      allocate(hoc(n1:n2,n1:n2,n1:n2),ll(nptile),llgp(nptile),hcgp(nptile),ecgp(nptile))
      hoc=0; ll=0
      do ip=1,nptile
        idx=ceiling(xv(1:3,ip)*n_refine) ! index of the grid
        if (minval(idx)<n1 .or. maxval(idx)>n2) then
          print*, 'idx out of range'
          print*, idx
          stop
        endif
        ll(ip)=hoc(idx(1),idx(2),idx(3)) ! linked list 
        hoc(idx(1),idx(2),idx(3))=ip ! head of chain
        hcgp(ip)=ip ! initialize hcgp(ip)=ip for isolated particles
      enddo
      llgp=0; ecgp=0; ! initialize group link list
      if (head) print*, 'percolation'

      ! loop over fof cells
      do iq3=n1,n2
      do iq2=n1,n2
      do iq1=n1,n2
        ip=hoc(iq1,iq2,iq3)
        do while (ip/=0)
          jp=ll(ip)
          do while (jp/=0)
            rsq=sum((xv(1:3,ip)-xv(1:3,jp))**2)
            if (rsq<=rp2) call merge_chain(ip,jp)
            jp=ll(jp)
          enddo
          do i_neighbor=1,13
            jq=[iq1,iq2,iq3]+ijk(:,i_neighbor)
            if (minval(jq)<n1 .or. maxval(jq)>n2) cycle
            jp=hoc(jq(1),jq(2),jq(3))
            do while (jp/=0)
              rsq=sum((xv(1:3,ip)-xv(1:3,jp))**2)
              if (rsq<=rp2) call merge_chain(ip,jp)
              jp=ll(jp)
            enddo
          enddo
          ip=ll(ip) ! find next particle in the chain
        enddo
      enddo
      enddo
      enddo

      !print*, 'count'
      ! count isolated, member and leader particles
      np_iso=0; np_mem=0; np_head=0;
      do i=1,nptile
        if (hcgp(i)==i) then
          np_iso=np_iso+1
        elseif (hcgp(i)==0) then
          np_mem=np_mem+1
        else
          np_head=np_head+1
        endif
      enddo
      !print*,'N_iso,mem,head =',np_iso,np_mem,np_head

      !print*,'fof particles'
      nlist=np_head
      allocate(xv_mean(6,nlist),iph_halo_all(nlist),dxv(6,nlist),np_halo_all(nlist))

      nh_tile=0;
      do ip=1,nptile
        if (hcgp(ip)/=ip .and. hcgp(ip)/=0) then
          np=0; jp=hcgp(ip) ! hoc of the group
          do while (jp/=0) ! loop over group members
            np=np+1; dxv(:,np)=xv(:,jp);
            jp=llgp(jp) ! next particle in chain
          enddo
          if (np>=np_halo_min .and. minval(sum(dxv(1:3,:np),2)/np)>=0 .and. & 
                                    maxval(sum(dxv(1:3,:np),2)/np)<nt*ratio_cs) then
            nh_tile=nh_tile+1
            np_halo_all(nh_tile)=np
            iph_halo_all(nh_tile)=hcgp(ip) ! hoc of the halo "ip-header"
            xv_mean(:,nh_tile)=sum(dxv(:,:np),2)/np+(ixyz2(:,itile)-1)*ngp+([icx,icy,icz]-1)*ng
          endif
        endif
      enddo
      !print*,'nh_tile',nh_tile
      !print*, minval(xv_mean(:,:nhalo)),maxval(xv_mean(:,:nhalo))

      ! transfer data to smaller arrays
      allocate(hcat(nh_tile),iph_halo(nh_tile)); iph_halo=iph_halo_all(:nh_tile); hcat%hmass=np_halo_all(:nh_tile)
      do i=1,6
        hcat%xv(i)=xv_mean(i,1:nh_tile)
      enddo
      nh=nh+nh_tile
      write(21) hcat
      deallocate(np_halo_all,xv_mean,dxv,iph_halo_all,iph_halo)
      deallocate(xv,hoc,ll,llgp,hcgp,ecgp,hcat)
    enddo ! itile

    halo_header%nhalo=nh; sync all
    if (head) then
      do i=2,nn**3
        nh=nh+nh[i]
      enddo
    endif
    sync all; nh=nh[1]; halo_header%nhalo_tot=nh; 
    rewind(21); write(21) halo_header; close(21)
    sync all
    deallocate(xp,vp,rhoc,pid)
    if (head) print*,'halo_header',halo_header
  enddo ! cur_checkpoint



  contains

  subroutine merge_chain(ii,jj)
    ! llgp is a linked list: llgp(ip) means ip->llgp
    ! ip1->ip2->...->ipn
    ! ip1 is the head of chain (hoc)
    ! ipn is the end of chain (eoc)
    integer(8) ii,jj,ihead,jhead,iend,jend,ipart
    jend=merge(jj,ecgp(jj),ecgp(jj)==0)
    iend=merge(ii,ecgp(ii),ecgp(ii)==0)
    if (iend==jend) return ! same chain
    ihead=max(hcgp(ii),hcgp(iend))
    jhead=max(hcgp(jj),hcgp(jend))
    ipart=jhead ! change eoc of the chain-j
    do while (ipart/=0)
      ecgp(ipart)=iend ! set chain-j's eoc to be iend
      ipart=llgp(ipart) ! next particle
    enddo
    llgp(jend)=ihead ! link j group to i group
    ecgp(ii)=iend
    hcgp(iend)=jhead ! change hoc
    hcgp(jend)=0 ! set jend as a member
  endsubroutine
  
end
