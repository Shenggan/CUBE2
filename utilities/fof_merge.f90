program CUBE_FoF_merge
  use variables
  implicit none
  integer i,cur_checkpoint
  type(type_halo_catalog_header) halo_header
  type(type_halo_catalog_array),allocatable :: hcat(:),hcat_max(:)[:]

  open(16,file='z_checkpoint.txt',status='old')
  do i=1,nmax_redshift-1
    read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
  enddo
  71 n_checkpoint=i-1
  close(16)
  if (n_checkpoint==0) stop 'z_checkpoint.txt empty'

  do cur_checkpoint=79,79
    sim%cur_checkpoint=cur_checkpoint
    print*, 'FoF_merge at redshift: ',z2str(z_checkpoint(cur_checkpoint))
    image=1; open(20,file=output_name('halo_all'),status='replace',access='stream')
    do i=1,nn**3
      image=i
      open(21,file=output_name('halo'),status='old',access='stream');
        read(21) halo_header
        if (i==1) write(20) halo_header
        allocate(hcat(halo_header%nhalo))
        read(21) hcat
        write(20) hcat
        deallocate(hcat)
      close(21)
    enddo
    close(20)
    print*,'wrote:',output_name('halo_all')
  enddo

end