!=========================================================
!
!  particles.f90
!    * Multiple particles
!
!=========================================================

module particles_m
  use const_m
  use vect2d_m
  implicit none
  private
  public :: particles__print,  &
            particles__save,  &
            particles__set

  integer(SI), parameter :: N_PARTICLES = 2

  type, public :: particles_t
    type(vect2d_t), dimension(N_PARTICLES) :: pos
  end type particles_t

  type(particles_t) :: particles

contains

  subroutine particles__print
    integer(SI) :: n

    do n = 1 , N_PARTICLES
      print *, n, ': (px,py)=(',  &
                              particles%pos(n)%x,  &
                              ', ',  &
                              particles%pos(n)%y, &
                            ')'
    end do
  end subroutine particles__print

  subroutine particles__save
    integer(SI) :: n
    integer(SI) :: file_num = 10
    integer(SI), save :: counter = 0

    open(file_num, file='particles.pos.data')
      do n = 1 , N_PARTICLES
        write(file_num,*) particles%pos(n)%x, particles%pos(n)%y
      end do
    close(file_num)

    counter = counter + 1
  end subroutine particles__save

  subroutine particles__set
    real(DR) :: rand(N_PARTICLES)

    call random_number(rand)
    particles%pos(:)%x = rand(:)

    call random_number(rand)
    particles%pos(:)%y = rand(:)
  end subroutine particles__set

end module particles_m
