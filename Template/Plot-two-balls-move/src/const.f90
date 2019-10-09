!=========================================================
!
!  const.f90
!    * Mathematical, physical, and Fortran constants.
!
!=========================================================

module const_m
  implicit none
  private
  public :: const__print
  public :: SI, DI, SR, DR, PI, TWOPI

  !<< Fortran constants >>!
  integer, parameter :: SI = selected_int_kind(6)
  integer, parameter :: DI = selected_int_kind(15)
  integer, parameter :: SR = selected_real_kind(6)
  integer, parameter :: DR = selected_real_kind(15)

  !<< Mathematical constants >>!
  real(DR), parameter ::    PI = 3.1415926535897932_DR
  real(DR), parameter :: TWOPI = PI*2


contains


  subroutine print_DR(string,d)
    character(len=*), intent(in) :: string
    real(DR)   , intent(in) :: d
    print *, string, d
  end subroutine print_DR


  subroutine print_SI(string,i)
    character(len=*), intent(in) :: string
    integer(SI)   , intent(in) :: i

    print *, string, i
  end subroutine print_SI


!  private
!=================
!  public

  subroutine const__print
    print *, '==================================='
    call print_SI('const SI', SI)
    call print_SI('const DI', DI)
    call print_SI('const SR', SR)
    call print_SI('const DR', DR)
    call print_DR('const PI', PI)
    call print_DR('const TWOPI', TWOPI)
    print *, '==================================='
  end subroutine const__print

end module const_m
