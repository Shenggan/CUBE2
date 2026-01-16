!#define macbook
!  if defined macbook, then omp-fftw is disabled.
module cubefft
  contains
  subroutine create_cubefft_plan
    use,intrinsic :: ISO_C_BINDING
    use variables
    implicit none
    
    !include 'fftw3.f'
    include 'fftw3.f03'
    
#ifndef macbook
    integer istat,nthreads
    call sfftw_init_threads(istat)
    if(head) print*, '    sfftw_init_threads status',istat
    nthreads=omp_get_max_threads()
    if(head) print*, '    omp_get_max_threads() =',nthreads
    call sfftw_plan_with_nthreads(nnest)
#endif
    do iteam=1,nteam
      call sfftw_plan_dft_r2c_3d(plan2(iteam),ngt,ngt,ngt,rho2(:,:,:,iteam),rho2(:,:,:,iteam),FFTW_MEASURE)
      call sfftw_plan_dft_c2r_3d(iplan2(iteam),ngt,ngt,ngt,rho2(:,:,:,iteam),rho2(:,:,:,iteam),FFTW_MEASURE)
      call sfftw_plan_dft_r2c_3d( plan3(iteam,2),nft(2),nft(2),nft(2),rho3_2(:,:,:,iteam),rho3_2(:,:,:,iteam),FFTW_MEASURE)
      call sfftw_plan_dft_c2r_3d(iplan3(iteam,2),nft(2),nft(2),nft(2),rho3_2(:,:,:,iteam),rho3_2(:,:,:,iteam),FFTW_MEASURE)
      call sfftw_plan_dft_r2c_3d( plan3(iteam,3),nft(3),nft(3),nft(3),rho3_4(:,:,:,iteam),rho3_4(:,:,:,iteam),FFTW_MEASURE)
      call sfftw_plan_dft_c2r_3d(iplan3(iteam,3),nft(3),nft(3),nft(3),rho3_4(:,:,:,iteam),rho3_4(:,:,:,iteam),FFTW_MEASURE)
      call sfftw_plan_dft_r2c_3d( plan3(iteam,4),nft(4),nft(4),nft(4),rho3_6(:,:,:,iteam),rho3_6(:,:,:,iteam),FFTW_MEASURE)
      call sfftw_plan_dft_c2r_3d(iplan3(iteam,4),nft(4),nft(4),nft(4),rho3_6(:,:,:,iteam),rho3_6(:,:,:,iteam),FFTW_MEASURE)
      call sfftw_plan_dft_r2c_3d( plan3(iteam,5),nft(5),nft(5),nft(5),rho3_8(:,:,:,iteam),rho3_8(:,:,:,iteam),FFTW_MEASURE)
      call sfftw_plan_dft_c2r_3d(iplan3(iteam,5),nft(5),nft(5),nft(5),rho3_8(:,:,:,iteam),rho3_8(:,:,:,iteam),FFTW_MEASURE)
      call sfftw_plan_dft_r2c_3d( plan3(iteam,6),nft(6),nft(6),nft(6),rho3_12(:,:,:,iteam),rho3_12(:,:,:,iteam),FFTW_MEASURE)
      call sfftw_plan_dft_c2r_3d(iplan3(iteam,6),nft(6),nft(6),nft(6),rho3_12(:,:,:,iteam),rho3_12(:,:,:,iteam),FFTW_MEASURE)
      !call sfftw_plan_dft_r2c_3d( plan3(iteam,7),nft(7),nft(7),nft(7),rho3_16(:,:,:,iteam),rho3_16(:,:,:,iteam),FFTW_MEASURE)
      !call sfftw_plan_dft_c2r_3d(iplan3(iteam,7),nft(7),nft(7),nft(7),rho3_16(:,:,:,iteam),rho3_16(:,:,:,iteam),FFTW_MEASURE)
    enddo
    
    !call sfftw_plan_dft_r2c_3d(plan0,ng,ng,ng,rho0,rho0,FFTW_MEASURE)
    !call sfftw_plan_dft_c2r_3d(iplan0,ng,ng,ng,rho0,rho0,FFTW_MEASURE)
  endsubroutine create_cubefft_plan

  subroutine destroy_cubefft_plan
    use,intrinsic :: ISO_C_BINDING
    use variables
    implicit none
    !include 'fftw3.f'
    include 'fftw3.f03'
    call sfftw_destroy_plan(plan2)
    call sfftw_destroy_plan(iplan2)
    do iapm=2,6 ! resolution
    do iteam=1,nteam 
      call sfftw_destroy_plan( plan3(iteam,iapm))
      call sfftw_destroy_plan(iplan3(iteam,iapm))
    enddo
    enddo
    !call sfftw_destroy_plan(plan0)
    !call sfftw_destroy_plan(iplan0)
    !call fftw_cleanup_threads()
  endsubroutine destroy_cubefft_plan

endmodule
