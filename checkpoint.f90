subroutine checkpoint
  use omp_lib
  use variables
  implicit none
  save

  character(100) fn10,fn11,fn12,fn13,fn14,fn15

  if (head) print*, 'checkpoint'

  fn10=output_name('info')
  fn11=output_name('xp')
  fn12=output_name('vp')
  fn13=output_name('np')
  fn14=output_name('vc')
  fn15=output_name('id')

  sim%vsim2phys=(1.5/sim%a)*box*100.*sqrt(omega_m)/ng_global
  sim%sigma_vi=sigma_vi
  !$omp parallelsections default(shared)
  !$omp section
  open(10,file=fn10,status='replace',access='stream'); write(10) sim; close(10)
  !$omp section
  open(11,file=fn11,status='replace',access='stream'); write(11) xp(:,:sim%nplocal); close(11)
  !$omp section
  open(12,file=fn12,status='replace',access='stream'); write(12) vp(:,:sim%nplocal); close(12)
  !$omp section
  open(13,file=fn13,status='replace',access='stream'); write(13) rhoc(1:nt,1:nt,1:nt,:,:,:); close(13)
!!  !$omp section
!!  open(14,file=fn14,status='replace',access='stream'); write(14) vfield(:,1:nt,1:nt,1:nt,:,:,:); close(14)
# ifdef PID
    !$omp section
    open(15,file=fn15,status='replace',access='stream'); write(15) pid(:sim%nplocal); close(15)
    if (head) print*, 'check PID range: ',minval(pid(:sim%nplocal)),maxval(pid(:sim%nplocal))
    if (minval(pid(:sim%nplocal))<1) then
      print*, 'pid are not all positive'
    endif
# endif
  !$omp endparallelsections
  if (this_image()==1 .or. this_image()==num_images()) print*,'  image',this_image(),'wrote',sim%nplocal,'CDM particles'
  sync all
endsubroutine
