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
  public :: particles__move,  &
            particles__print,  &
            particles__save,  &
            particles__set

  integer(SI), parameter :: N_PARTICLES = 2

  type, public :: particles_t
    type(vect2d_t), dimension(N_PARTICLES) :: pos
    type(vect2d_t), dimension(N_PARTICLES) :: vel
  end type particles_t

  type(particles_t) :: particles

contains

  function int_to_str5(i) result(str)
    integer(SI), intent(in) :: i
    character(len=5) :: str
    !===============
    !  Convert an integer into 5 characters.
    !          i =   10 --> str="00010"
    !          i =  123 --> str="00123"
    !          i = -123 --> str="XXXXX"
    !===============
    if ( i < 0 .or. i > 99999 ) then
      str = 'XXXXX'
    else
      write(str, '(i5.5)') i
    end if
  end function int_to_str5

  subroutine particles__move(dt)
    real(DR), intent(in) :: dt

    particles%pos(:)%x = particles%pos(:)%x + particles%vel(:)%x * dt
    particles%pos(:)%y = particles%pos(:)%y + particles%vel(:)%y * dt

  end subroutine particles__move

  subroutine particles__print
    integer(SI) :: n

    do n = 1 , N_PARTICLES
      print *, n, ': (px,py)=(',  &
                              particles%pos(n)%x,  &
                              ', ',  &
                              particles%pos(n)%y, &
                            ')'
      print *, n, ': (vx,vy)=(',  &
                              particles%vel(n)%x,  &
                              ', ',  &
                              particles%vel(n)%y, &
                            ')'
    end do
  end subroutine particles__print

  subroutine particles__save
    integer(SI) :: n
    integer(SI) :: file_num = 10
    integer(SI), save :: counter = 0

    open(file_num, file='particles.pos.data.'//int_to_str5(counter))
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
    particles%vel(:)%x = 0.01_DR

    call random_number(rand)
    particles%pos(:)%y = rand(:)
    particles%vel(:)%y = -0.02_DR
  end subroutine particles__set

end module particles_m
