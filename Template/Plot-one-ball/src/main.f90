!=========================================================
!
!  main.f90
!    * Shows the basic usage of Fortran modules.
!
!=========================================================

program main
  use const_m
  use particle_m
  use vect2d_m
  implicit none

  call const__print

  call particle__set(0.3_DR, 0.7_DR)
  call particle__print
  call particle__save
end program main
