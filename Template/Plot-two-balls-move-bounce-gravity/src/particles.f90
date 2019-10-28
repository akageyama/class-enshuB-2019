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
  public :: particles
  public :: particles__move, &
            particles__print, &
            particles__print_energy, &
            particles__save, &
            particles__set

  integer(SI), parameter :: N_PARTICLES = 2

  real(DR), parameter :: MASS    = 1.e-3_DR  ! [Kg]

  real(DR), parameter :: RADIUS = 0.02_DR
  real(DR), parameter :: XMIN = 0.0_DR
  real(DR), parameter :: XMAX = 1.0_DR
  real(DR), parameter :: YMIN = 0.0_DR
  real(DR), parameter :: YMAX = 1.0_DR
  real(DR), parameter :: XMIN2 = XMIN + RADIUS
  real(DR), parameter :: XMAX2 = XMAX - RADIUS
  real(DR), parameter :: YMIN2 = YMIN + RADIUS
  real(DR), parameter :: YMAX2 = YMAX - RADIUS

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

    integer(SI) :: n

    do n = 1, N_PARTICLES
      particles%pos(n)%x = particles%pos(n)%x + particles%vel(n)%x * dt
      particles%pos(n)%y = particles%pos(n)%y + particles%vel(n)%y * dt
      particles%vel(n)%y = particles%vel(n)%y - G_ACCEL * dt

      if ( particles%pos(n)%x < XMIN2 .or. particles%pos(n)%x > XMAX2 ) then
        particles%vel(n)%x = -particles%vel(n)%x
      end if
      if ( particles%pos(n)%y < YMIN2 .or. particles%pos(n)%y > YMAX2 ) then
        particles%vel(n)%y = -particles%vel(n)%y
      end if
    end do

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


  subroutine particles__print_energy(t)
    real(DR), intent(in) :: t
    integer(SI) :: n
    real(DR) :: vel_sq
    real(DR) :: energy_velocity, energy_potential, sum

    sum = 0.0_DR

    do n = 1, N_PARTICLES
      vel_sq = particles%vel(n)%x**2 + particles%vel(n)%y**2
      energy_velocity = 0.5_DR * MASS * vel_sq
      energy_potential = MASS * G_ACCEL * particles%pos(n)%y
      sum = sum + ( energy_velocity + energy_potential )
    end do

    print *, 't, Total energy = ', t, sum
  end subroutine particles__print_energy


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
    particles%vel(:)%x = 0.1_DR

    call random_number(rand)
    particles%pos(:)%y = rand(:)
    particles%vel(:)%y = -0.2_DR
  end subroutine particles__set

end module particles_m
