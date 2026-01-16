subroutine Green_3D(Gk,nglobal,n1,n2,n3,kshift,apm_soft,apm_range,gridsize)
  use omp_lib
  use parameters
  implicit none
  integer(8) nglobal,n1,n2,n3
  real apm_soft,apm_range,gridsize,Gk(n1,n2,n3)
  integer k1,k2,k3,i1,i2,i3,kshift
  real(8) fftspace(nglobal),kvec(ndim),u2n_2,Dvec(ndim),u2n,ur_cn(ndim),kvec_n(ndim)
  real(8) Uk2_n,kmag_n,y,S_soft,S_range,Rvec_n(ndim)
  if (head) then
    call system_clock(t1,t_rate)
    print*,'  Green_3D'
    print*,'    mesh =',int(nglobal,kind=2)
    print*,'    n1,n2,n3 =',int(n1,kind=2),int(n2,kind=2),int(n3,kind=2)
    print*,'    apm_soft, apm_range =',apm_soft,apm_range
    print*,'    gridsize =',gridsize
    print*,'    n_int =',int(n_int,kind=2)
    print*,'    kshift =',kshift
  endif
  sync all

  fftspace=2*pi*(1./nglobal)*(mod([(k1,k1=1,nglobal)]+nglobal/2-1,nglobal)-nglobal/2)
  !$omp paralleldo default(shared) schedule(dynamic)&
  !$omp& private(k3,k2,k1,kvec,u2n_2,Dvec,u2n,ur_cn,i3,i2,i1,kvec_n,kmag_n,Uk2_n)&
  !$omp& private(y,S_soft,S_range,Rvec_n)
  do k3=1,n3
  do k2=1,n2
  do k1=1,n1
    kvec(1)=fftspace(k1)
    kvec(2)=fftspace(k2+kshift*n2*(icx-1))
    kvec(3)=fftspace(k3+kshift*n3*((icz-1)*nn+icy-1))
    !if (k1*k2*k3==1 .and. kshift==1) then
    !  print*, icx,icy,icz,k1,k2+kshift*n2*(icx-1), k3+kshift*n3*((icz-1)*nn+icy-1), kvec
    !endif
    u2n_2=product(1-sin(kvec/2)**2+(2./15)*sin(kvec/2)**4)
    Dvec=alpha*sin(kvec)+(1-alpha)*sin(2*kvec)/2
    !u2n=0
    ur_cn=0
    do i3=-n_int,n_int
    do i2=-n_int,n_int
    do i1=-n_int,n_int
      kvec_n=kvec+2*pi*[i1,i2,i3]
      kmag_n=norm2(kvec_n)
      Uk2_n=product(merge(1d0,sin(kvec_n/2)/(kvec_n/2),kvec_n==0))**(2*p+2)
      !u2n=u2n+Uk2_n
      y=apm_soft*kmag_n/2
      if (kmag_n==0) then
        S_soft=1
      else
        S_soft=12*(y**-4)*(2-2*cos(y)-y*sin(y))
      endif
      if (apm_range==0) then
        S_range=0
      elseif (kmag_n==0) then
        S_range=1
      else
        y=apm_range*kmag_n/2
        S_range=12*(y**-4)*(2-2*cos(y)-y*sin(y))
      endif
      Rvec_n=kvec_n*(S_soft**2-S_range**2)/kmag_n**2
      ur_cn=ur_cn+Uk2_n*Rvec_n
    enddo
    enddo
    enddo
    u2n=product(1-sin(kvec/2)**2+(2./15)*sin(kvec/2)**4);
    Gk(k1,k2,k3)=4*pi*sum(Dvec*ur_cn)/sum(Dvec**2)/u2n**2/gridsize**(ndim-1)
  enddo
  enddo
  enddo
  !$omp endparalleldo
  Gk(::nglobal/2,::nglobal/2,::nglobal/2)=0;
  sync all

  if (head) then
    call system_clock(t2,t_rate)
    print*,'    elapsed time =',real(t2-t1)/t_rate,'secs';
    print*,''
  endif
  sync all
endsubroutine

subroutine Green_3D_iso(Gk,nglobal,n1,n2,n3,apm_soft,apm_range,gridsize)
  use omp_lib
  use parameters
  implicit none
  integer(8) nglobal,n1,n2,n3
  real apm_soft,apm_range,gridsize
  real Gk(n1,n2,n3)
  integer k1,k2,k3,i1,i2,i3
  real(8) fftspace(nglobal),kvec(ndim),u2n_2,Dvec(ndim),u2n,ur_cn(ndim),kvec_n(ndim)
  real(8) Uk2_n,kmag_n,y,S_soft,S_range,Rvec_n(ndim)

  if (head) then
    call system_clock(t1,t_rate)
    print*,'  Initialize Green''s function'
    print*,'    mesh =',int(nglobal,kind=2)
    print*,'    n1,n2,n3 =',int(n1,kind=2),int(n2,kind=2),int(n3,kind=2)
    print*,'    apm_soft, apm_range =',apm_soft,apm_range
    print*,'    gridsize =',gridsize
    print*,'    n_int =',int(n_int,kind=2)
  endif

  fftspace=2*pi*(1./nglobal)*(mod([(k1,k1=1,nglobal)]+nglobal/2-1,nglobal)-nglobal/2)
  !$omp paralleldo default(shared) schedule(dynamic)&
  !$omp& private(k3,k2,k1,kvec,u2n_2,Dvec,u2n,ur_cn,i3,i2,i1,kvec_n,kmag_n,Uk2_n)&
  !$omp& private(y,S_soft,S_range,Rvec_n)
  do k3=1,n3
  do k2=1,n2
  do k1=1,n1
    kvec=[fftspace(k1),fftspace(k2),fftspace(k3)]
    u2n_2=product(1-sin(kvec/2)**2+(2./15)*sin(kvec/2)**4)
    Dvec=alpha*sin(kvec)+(1-alpha)*sin(2*kvec)/2
    !u2n=0
    ur_cn=0
    do i3=-n_int,n_int
    do i2=-n_int,n_int
    do i1=-n_int,n_int
      kvec_n=kvec+2*pi*[i1,i2,i3]
      kmag_n=norm2(kvec_n)
      Uk2_n=product(merge(1d0,sin(kvec_n/2)/(kvec_n/2),kvec_n==0))**(2*p+2)
      !u2n=u2n+Uk2_n
      y=apm_soft*kmag_n/2
      if (kmag_n==0) then
        S_soft=1
      else
        S_soft=12*(y**-4)*(2-2*cos(y)-y*sin(y))
      endif
      if (apm_range==0) then
        S_range=0
      elseif (kmag_n==0) then
        S_range=1
      else
        y=apm_range*kmag_n/2
        S_range=12*(y**-4)*(2-2*cos(y)-y*sin(y))
      endif
      Rvec_n=kvec_n*(S_soft**2-S_range**2)/kmag_n**2
      Rvec_n=Rvec_n*(kmag_n*nglobal/2-sin(kmag_n*nglobal/2))/(kmag_n*nglobal/2);
      ur_cn=ur_cn+Uk2_n*Rvec_n
    enddo
    enddo
    enddo
    u2n=product(1-sin(kvec/2)**2+(2./15)*sin(kvec/2)**4);
    Gk(k1,k2,k3)=4*pi*sum(Dvec*ur_cn)/sum(Dvec**2)/u2n**2/gridsize**(ndim-1)
  enddo
  enddo
  enddo
  !$omp endparalleldo
  Gk(::nglobal/2,::nglobal/2,::nglobal/2)=0;
  sync all

  if (head) then
    call system_clock(t2,t_rate)
    print*,'    elapsed time =',real(t2-t1)/t_rate,'secs';
    print*,''
  endif
  sync all
endsubroutine
