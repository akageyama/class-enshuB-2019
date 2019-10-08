!=========================================================
!
!  main.f90
!    * Shows the basic usage of Fortran modules.
!
!=========================================================

program main
  use const_m
  implicit none

  call const__print
  print *, 'SI = ', SI
end program
