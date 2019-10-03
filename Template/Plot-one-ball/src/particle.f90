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

  type, public :: particle_t
    type(vect2d_t) :: pos
  contains
    procedure :: set => particle__set
    procedure :: print => particle__print
    procedure :: save => particle__save
  end type particle_t

contains

  subroutine particle__print(self)
    class(particle_t), intent(in) :: self
    print *, '(x,y)=(',  &
                     self%pos%x,  &
                     ', ',  &
                     self%pos%y, &
                   ')'
  end subroutine particle__print

  subroutine particle__save(self)
    class(particle_t), intent(in) :: self
    logical, save :: first_time = .true.
    integer(SI) :: file_num = 10
    if ( first_time ) then
      open(file_num, file='particle.pos.data')
      first_time = .false.
    end if
    write(file_num,*) self%pos%x, self%pos%y
  end subroutine particle__save

  subroutine particle__set(self, x, y)
    class(particle_t), intent(out) :: self
    real(DR), intent(in) :: x, y

    self%pos%x = x
    self%pos%y = y
  end subroutine particle__set

end module particle_m
