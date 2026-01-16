subroutine timestep
  use variables
  implicit none
  save
  integer ntry
  real ra,da_1,da_2,a_next,z_next

  dt_old=dt
  sim%timestep=sim%timestep+1
  call tic(1)
  if (head) then
#ifndef SPEEDTEST
    open(101,file=output_name('tcat'),status='replace',access='stream')
    write(101) int(sim%timestep,kind=4), tcat(:,:sim%timestep)
    close(101)
#endif
    print*, ''
    print*, '-------------------------------------------------------'
    print*, 'timestep    :',sim%timestep
    dt_e=dt_max
    ntry=0
    do
      ntry=ntry+1
      call expansion(sim%a,dt_e,da_1,da_2)
      da=da_1+da_2
      ra=da/(sim%a+da)
      if (ra>ra_max) then
        dt_e=dt_e*(ra_max/ra)
      else
        exit
      endif
      if (ntry>10) exit
    enddo
    dt=min(dt_e,sim%dt_pm1,sim%dt_pm2,sim%dt_pm3,dt_refine*sim%dt_pp,sim%dt_vmax)
    call expansion(sim%a,dt,da_1,da_2)
    da=da_1+da_2
    ! check if checkpointing is needed
    checkpoint_step=.false.
    halofind_step=.false.
#   ifdef HALOFIND
      z_next=max(z_checkpoint(sim%cur_checkpoint),z_halofind(sim%cur_halofind))
#   else
      z_next=z_checkpoint(sim%cur_checkpoint)
#   endif
    a_next=1.0/(1+z_next)
    if (da>=a_next-sim%a) then
      if (z_next==z_checkpoint(sim%cur_checkpoint)) then
        checkpoint_step=.true.
        if (sim%cur_checkpoint==n_checkpoint) final_step=.true.
      endif
#   ifdef HALOFIND
      if (z_next==z_halofind(sim%cur_halofind)) halofind_step=.true.
#   endif
      do while (abs((sim%a+da)/a_next-1)>=1e-6)
        dt=dt*(a_next-sim%a)/da
        call expansion(sim%a,dt,da_1,da_2)
        da=da_1+da_2
        print*, 'a+da, dt, z+dz, err_a', sim%a+da, dt, 1.0/(sim%a+da)-1.0, (sim%a+da)/a_next-1
      enddo
    endif

    ra=da/(sim%a+da)
    a_mid=sim%a+(da/2)

    tcat(41,istep)=sim%a
    tcat(42,istep)=a_mid
    tcat(43,istep)=sim%a+da

    !print*, 'tau         :',tau,tau+dt
    print*, 'z         :',1.0/sim%a-1.0,1.0/(sim%a+da)-1.0
    print*, 'a         :',sim%a,a_mid,sim%a+da
    print*, 'expansion :',ra
    print*, 'dt        :',dt
    print*, 'dt_e      :',dt_e
    print*, 'dt_pm1    :',sim%dt_pm1
    print*, 'dt_pm2    :',sim%dt_pm2
    print*, 'dt_pm3    :',sim%dt_pm3
    print*, 'dt_pp     :',sim%dt_pp
    print*, 'dt_vmax   :',sim%dt_vmax
    print*, ''
    sim%tau=sim%tau+dt
    sim%t=sim%t+dt
    sim%a=sim%a+da
  endif
  sync all

  !a=a[1]
  a_mid=a_mid[1]
  dt=dt[1]
  checkpoint_step=checkpoint_step[1]
#ifdef HALOFIND
  halofind_step=halofind_step[1]
#endif
  final_step=final_step[1]
  call toc(1)
  sync all
endsubroutine timestep

#ifdef new_expansion
subroutine expansion(a0,dt0,da1,da2)
  !! Expansion subroutine :: Hy Trac -- trac@cita.utoronto.ca
  !! Added Equation of State for Dark Energy :: Pat McDonald -- pmcdonal@cita.utoronto.ca
  use variables
  implicit none
  save
  real(4) :: a0,dt0,dt_x,da1,da2
  real(8) :: a_x,adot,addot,atdot,arkm,am1rrm,a3rlm,omHsq
  real(8), parameter :: e = 2.718281828459046
  !! Expand Friedman equation to third order and integrate
  dt_x=dt0/2
  omHsq=4.0/9.0

  a_x=a0
  a3rlm=a_x**(-3*wde)*omega_l/omega_m
  arkm=a_x*(1.0-omega_m-omega_l)/omega_m
  am1rrm=a_x**(-1.)*omega_r/omega_m

  adot=sqrt(omHsq*a_x**3*(1.0+arkm+a3rlm+am1rrm))
  addot=a_x**2*omHsq*(1.5+2.0*arkm+1.5*(1.0-wde)*a3rlm+am1rrm)
  atdot=a_x*adot*omHsq*(3.0+6.0*arkm+1.5*(2.0-3.0*wde)*(1.0-wde)*a3rlm+am1rrm)

  da1=adot*dt_x+(addot*dt_x**2)/2.0+(atdot*dt_x**3)/6.0

  a_x=a0+da1

  a3rlm=a_x**(-3*wde)*omega_l/omega_m
  arkm=a_x*(1.0-omega_m-omega_l)/omega_m
  am1rrm=a_x**(-1.)*omega_r/omega_m

  adot=sqrt(omHsq*a_x**3*(1.0+arkm+a3rlm+am1rrm))
  addot=a_x**2*omHsq*(1.5+2.0*arkm+1.5*(1.0-wde)*a3rlm+am1rrm)
  atdot=a_x*adot*omHsq*(3.0+6.0*arkm+1.5*(2.0-3.0*wde)*(1.0-wde)*a3rlm+am1rrm)

  da2=adot*dt_x+(addot*dt_x**2)/2.0+(atdot*dt_x**3)/6.0

endsubroutine expansion
#else
subroutine expansion(a0,dt0,da1,da2)
  !! Expansion subroutine :: Hy Trac -- trac@cita.utoronto.ca
  !! Added Equation of State for Dark Energy :: Pat McDonald -- pmcdonal@cita.utoronto.ca
  use variables
  implicit none
  save
  real(4) :: a0,dt0,dt_x,da1,da2
  real(8) :: a_x,adot,addot,atdot,arkm,a3rlm,omHsq
  real(8), parameter :: e = 2.718281828459046
  !! Expand Friedman equation to third order and integrate
  dt_x=dt0/2
  a_x=a0
  omHsq=4.0/9.0
  a3rlm=a_x**(-3*wde)*omega_l/omega_m
  !a3rlm=a_x**(-3*wde - 3*w_a)*(omega_l/omega_m)*e**(3*w_a*(a_x - 1))
  arkm=a_x*(1.0-omega_m-omega_l)/omega_m

  adot=sqrt(omHsq*a_x**3*(1.0+arkm+a3rlm))
  !    addot=a_x**2*omHsq*(1.5+2.0*arkm+3.0*a3rlm)
  addot=a_x**2*omHsq*(1.5+2.0*arkm+1.5*(1.0-wde)*a3rlm)
  !    addot=a_x**2*omHsq*(1.5+2.0*arkm+1.5*(1.0-wde + w_a*(a_x - 1))*a3rlm)
  !    atdot=a_x*adot*omHsq*(3.0+6.0*arkm+15.0*a3rlm)
  atdot=a_x*adot*omHsq*(3.0+6.0*arkm+1.5*(2.0-3.0*wde)*(1.0-wde)*a3rlm)
  !    atdot=a_x*adot*omHsq*(3.0+6.0*arkm+1.5*(3*w_a**2*a_x**2 + 6*w_a*(1-wde-w_a)+ (2.0-3.0*(wde+w_a))*(1.0-(wde+w_a)))*a3rlm)

  da1=adot*dt_x+(addot*dt_x**2)/2.0+(atdot*dt_x**3)/6.0

  a_x=a0+da1
  omHsq=4.0/9.0
  !    a3rlm=a_x**3*omega_l/omega_m
  a3rlm=a_x**(-3*wde)*omega_l/omega_m
  !    a3rlm=a_x**(-3*wde - 3*w_a)*(omega_l/omega_m)*e**(3*w_a*(a_x - 1))
  arkm=a_x*(1.0-omega_m-omega_l)/omega_m

  adot=sqrt(omHsq*a_x**3*(1.0+arkm+a3rlm))
  !    addot=a_x**2*omHsq*(1.5+2.0*arkm+3.0*a3rlm)
  addot=a_x**2*omHsq*(1.5+2.0*arkm+1.5*(1.0-wde)*a3rlm)
  !    addot=a_x**2*omHsq*(1.5+2.0*arkm+1.5*(1.0-wde + w_a*(a_x - 1))*a3rlm)
  !    atdot=a_x*adot*omHsq*(3.0+6.0*arkm+15.0*a3rlm)
  atdot=a_x*adot*omHsq*(3.0+6.0*arkm+1.5*(2.0-3.0*wde)*(1.0-wde)*a3rlm)
  !    atdot=a_x*adot*omHsq*(3.0+6.0*arkm+1.5*(3*w_a**2*a_x**2 + 6*w_a*(1-wde-w_a)+ (2.0-3.0*(wde+w_a))*(1.0-(wde+w_a)))*a3rlm)

  da2=adot*dt_x+(addot*dt_x**2)/2.0+(atdot*dt_x**3)/6.0

endsubroutine expansion
#endif
