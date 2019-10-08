!=========================================================
!
!  main.f90
!    * Shows the basic usage of Fortran modules.
!
!=========================================================

program main
  use const_m
  use particles_m
  use vect2d_m
  implicit none
  integer(DI) :: loop
  integer(DI), parameter :: LOOP_MAX = 1000
  real(DR) :: dt = 5.e-2_DR

  type(particles_t) :: particles

  call const__print

  call particles__set
  call particles__print

  do loop = 0, LOOP_MAX
    call particles__move(dt)
    if ( mod(loop, LOOP_MAX/100)==0 ) then
      call particles__save
    end if
  end do
end program main
