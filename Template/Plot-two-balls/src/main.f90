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

  type(particles_t) :: particles

  call const__print

  call particles__set
  call particles__print

  call particles__save
end program main
