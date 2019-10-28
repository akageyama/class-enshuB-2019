!=========================================================
!
!  plot-two-balls-move-bounce-gravity/src/main.f90
!    * With gravity
!
!=========================================================

program main
  use const_m
  use particles_m
  use vect2d_m
  implicit none
  integer(DI) :: loop
  integer(DI), parameter :: LOOP_MAX = 500
  real(DR) :: dt = 1.e-3_DR
  real(DR) :: time = 0.0_DR

  call const__print

  call particles__set
  call particles__print

  do loop = 0, LOOP_MAX
    call particles__move(time, dt)
    time = time + dt
    if ( mod(loop, 10)==0 ) then
      call particles__print_energy(time)
      call particles__save
    end if
  end do
end program main
