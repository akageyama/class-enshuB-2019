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
  public :: G_ACCEL

  !<< Fortran constants >>!
  integer, parameter :: SI = selected_int_kind(6)
  integer, parameter :: DI = selected_int_kind(15)
  integer, parameter :: SR = selected_real_kind(6)
  integer, parameter :: DR = selected_real_kind(15)

  !<< Mathematical constants >>!
  real(DR), parameter ::    PI = 3.1415926535897932_DR
  real(DR), parameter :: TWOPI = PI*2

  !<< Physical constants >>!
  real(DR) :: G_ACCEL = 9.80665_DR
                        ! [m/s^2] Gravity acceleration

  integer(SI), parameter :: SOUT = 6


contains


  subroutine print_DR(string,d)
    character(len=*), intent(in) :: string
    real(DR)   , intent(in) :: d

    !integer(SI), parameter :: TOTAL_LENGTH = 60
    !character(len=22) :: string_for_value
    !character(len=TOTAL_LENGTH) :: line

    !line = repeat('.',TOTAL_LENGTH)

    !write(string_for_value,'(1pe22.15)') d
    !line(1:len_trim(string)) = trim(string)
    !line(TOTAL_LENGTH-22:TOTAL_LENGTH) = string_for_value
    !write(SOUT,*) line
    print *, string // ': ', d
  end subroutine print_DR


  subroutine print_SI(string,i)
    character(len=*), intent(in) :: string
    integer(SI), intent(in) :: i

    ! integer(SI), parameter :: TOTAL_LENGTH = 60
    ! character(len=12) :: string_for_i01
    ! character(len=TOTAL_LENGTH) :: line

    ! line = repeat('.', TOTAL_LENGTH)

    ! write(string_for_i01,'(a1,i0)') ' ', i ! put a space in front of i
    ! line(1:len_trim(string)) = trim(string)
    ! line(TOTAL_LENGTH-len_trim(string_for_i01):TOTAL_LENGTH)    &
    !    = trim(string_for_i01)
    ! write(SOUT,*) line
    print *, string // ': ' , i
  end subroutine print_SI


!  private
!=================
!  public

  subroutine const__print
    call print_SI('const SI', SI)
    call print_SI('const DI', DI)
    call print_SI('const SR', SR)
    call print_SI('const DR', DR)
    call print_DR('const PI', PI)
    call print_DR('const TWOPI', TWOPI)
  end subroutine const__print

end module const_m
