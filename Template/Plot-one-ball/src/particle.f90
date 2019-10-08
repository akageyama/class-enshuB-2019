!=========================================================
!
!  particle.f90
!    * A particle
!
!=========================================================

module particle_m
  use const_m
  use vect2d_m
  implicit none
  private
  public :: particle_t
  public :: particle__set,  &
            particle__print,  &
            particle__save

  type :: particle_t
    type(vect2d_t) :: pos
  end type particle_t

  type(particle_t) :: particle

contains

  subroutine particle__print
    print *, '(x,y)=(',  &
                     particle%pos%x,  &
                     ', ',  &
                     particle%pos%y, &
                   ')'
  end subroutine particle__print

  subroutine particle__save
    logical, save :: first_time = .true.
    integer(SI) :: file_num = 10
    if ( first_time ) then
      open(file_num, file='./particle.pos.data')
      first_time = .false.
    end if
    write(file_num,*) particle%pos%x, particle%pos%y
  end subroutine particle__save

  subroutine particle__set(x, y)
    real(DR), intent(in) :: x, y

    particle%pos%x = x
    particle%pos%y = y
  end subroutine particle__set

end module particle_m
