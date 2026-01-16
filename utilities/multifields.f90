! split or merge fields
program multifields
  use parameters
  implicit none
  integer ix,iy,iz,idx(3),flag
  real field_local(ng,ng,ng)[nn,nn,*],field_global(ng_global,ng_global,ng_global)
  call geometry
  sync all

  flag=2 ! 1=split, 2=merge

  if (flag==1) then
    if (head) print*, 'split field'
    !! split field into nn*3 sub fields
    ! set opath to be the output
    open(11,file='../output/universe37/image1/noise_1.bin',access='stream')
    read(11) field_global
    close(11)
    idx=[icx,icy,icz]*ng ! upper index of each dim
    field_local=field_global(idx(1)-ng+1:idx(1),idx(2)-ng+1:idx(2),idx(3)-ng+1:idx(3))
    open(11,file=output_dir()//'noise'//output_suffix(),access='stream')
    write(11) field_local
    close(11)
  elseif (flag==2) then
    if (head) print*, 'merge fields'
    !! merge nn*3 sub-fields into one field
    ! set opath to be the input
    open(11,file=output_dir()//'0.000_delta_c'//output_suffix(),access='stream')
    read(11) field_local
    close(11)
    sync all
    if (head) then
      do iz=1,nn
      do iy=1,nn
      do ix=1,nn
        idx=[ix,iy,iz]*ng ! upper index of each dim
        field_global(idx(1)-ng+1:idx(1),idx(2)-ng+1:idx(2),idx(3)-ng+1:idx(3))=field_local(:,:,:)[ix,iy,iz]
      enddo
      enddo
      enddo
      open(11,file='../output/universe37/image1/0.000_delta_c_merge_1.bin',access='stream')
      write(11) field_global
      close(11)
    endif
    sync all
  endif
end
