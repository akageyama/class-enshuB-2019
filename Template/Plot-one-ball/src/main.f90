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

  type(particle_t) :: p

  call const__print

  call p%set(0.3_DR, 0.7_DR)
  call p%print
  call p%save
end program main
