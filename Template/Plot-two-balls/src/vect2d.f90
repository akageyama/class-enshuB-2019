!=========================================================
!
!  vect2d.f90
!    * 2-dimensional vector
!
!=========================================================

module vect2d_m
  use const_m
  implicit none
  private

  type, public :: vect2d_t
    real(DR) :: x, y
  end type vect2d_t

end module vect2d_m
