!=========================================================
!
!  particles.f90
!    * To solve equation of motion by Runge-Kutta method
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
  real(DR), parameter :: MASS = 1.e-3_DR ! [Kg]

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

  interface runge_kutta
     module procedure rk_step_01_02_03, &
                      rk_step_04
  end interface


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


  subroutine bounce(ps)
    type(particles_t), intent(inout) :: ps
    integer(SI) :: n

    do n = 1, N_PARTICLES
      if ( ps%pos(n)%x < XMIN2 .and. ps%vel(n)%x < 0.0_DR ) then
        ps%vel(n)%x = -ps%vel(n)%x
      end if
      if ( ps%pos(n)%x > XMAX2 .and. ps%vel(n)%x > 0.0_DR ) then
        ps%vel(n)%x = -ps%vel(n)%x
      end if
      if ( ps%pos(n)%y < YMIN2 .and. ps%vel(n)%y < 0.0_DR ) then
        ps%vel(n)%y = -ps%vel(n)%y
      end if
      if ( ps%pos(n)%y > YMAX2 .and. ps%vel(n)%y > 0.0_DR ) then
        ps%vel(n)%y = -ps%vel(n)%y
      end if
    end do
  end subroutine bounce


  subroutine equation_of_motion( ps, t, dt, dps)
    type(particles_t), intent(in) :: ps
    real(DR), intent(in) :: t, dt
    type(particles_t), intent(out) :: dps

    ! (d/dt)(pos%x) = vel%x  ==> d(pos%x) = (vel%x)*dt
    ! (d/dt)(pos%y) = vel%y  ==> d(pos%y) = (vel%y)*dt
    ! (d/dt)(vel%x) = 0.0    ==> d(vel%x) =  0.0
    ! (d/dt)(vel%y) = -G     ==> d(vel%y) = -G*dt

    integer(SI) :: n

    do n = 1, N_PARTICLES
      dps%pos(n)%x = ps%vel(n)%x * dt
      dps%pos(n)%y = ps%vel(n)%y * dt
      dps%vel(n)%x = 0.0_DR
      dps%vel(n)%y = - G_ACCEL * dt
    end do
  end subroutine equation_of_motion


  subroutine rk_step_01_02_03( step, ps, dps, rk )
    character(len=3), intent(in) :: step
    type(particles_t), intent(in) :: ps, dps
    type(particles_t), intent(out) :: rk
!
!   * 1st to 3rd steps in 4-stage-Runge-Kutta method.
!
    integer(SI) :: n

    if ( step=='1st' .or. step=='2nd' ) then
      do n = 1, N_PARTICLES
        rk%pos(n)%x = ps%pos(n)%x + 0.5_DR*dps%pos(n)%x
        rk%pos(n)%y = ps%pos(n)%y + 0.5_DR*dps%pos(n)%y
        rk%vel(n)%x = ps%vel(n)%x + 0.5_DR*dps%vel(n)%x
        rk%vel(n)%y = ps%vel(n)%y + 0.5_DR*dps%vel(n)%y
      end do
    else if ( step=='3rd' ) then
      do n = 1, N_PARTICLES
        rk%pos(n)%x = ps%pos(n)%x + dps%pos(n)%x
        rk%pos(n)%y = ps%pos(n)%y + dps%pos(n)%y
        rk%vel(n)%x = ps%vel(n)%x + dps%vel(n)%x
        rk%vel(n)%y = ps%vel(n)%y + dps%vel(n)%y
      end do
    else
      print *,'<particles/rk_step_01_02_03> *** Bad arg. ***'
      stop
    end if

  end subroutine rk_step_01_02_03


  subroutine rk_step_04( ps, d1, d2, d3, d4 )
    type(particles_t), intent(inout) :: ps
    type(particles_t), intent(in) :: d1, d2, d3, d4
!
!   * The final step of the 4-th order Runge-Kutta method.
!
    real(DR), parameter :: ONE6TH = 1.0_DR / 6.0_DR
    real(DR), parameter :: ONE3RD = 1.0_DR / 3.0_DR

    integer(SI) :: n

    do n = 1, N_PARTICLES
      ps%pos(n)%x = ps%pos(n)%x  &
                  + ONE6TH*( d1%pos(n)%x &
                           + d4%pos(n)%x )  &
                  + ONE3RD*( d2%pos(n)%x &
                           + d3%pos(n)%x )
      ps%pos(n)%y = ps%pos(n)%y  &
                  + ONE6TH*( d1%pos(n)%y &
                           + d4%pos(n)%y )  &
                  + ONE3RD*( d2%pos(n)%y &
                           + d3%pos(n)%y )
      ps%vel(n)%x = ps%vel(n)%x  &
                  + ONE6TH*( d1%vel(n)%x &
                           + d4%vel(n)%x )  &
                  + ONE3RD*( d2%vel(n)%x &
                           + d3%vel(n)%x )
      ps%vel(n)%y = ps%vel(n)%y  &
                  + ONE6TH*( d1%vel(n)%y &
                           + d4%vel(n)%y )  &
                  + ONE3RD*( d2%vel(n)%y &
                           + d3%vel(n)%y )
    end do

  end subroutine rk_step_04


  subroutine zero_set( work, d1, d2, d3, d4 )
    type(particles_t), intent(out) :: work, d1, d2, d3, d4

    call iZero(work)
    call iZero(d1)
    call iZero(d2)
    call iZero(d3)
    call iZero(d4)

  contains

    subroutine iZero(ps)
      type(particles_t), intent(out) :: ps
      ps%pos(:)%x = 0.0_DR
      ps%pos(:)%y = 0.0_DR
      ps%vel(:)%x = 0.0_DR
      ps%vel(:)%y = 0.0_DR
    end subroutine iZero

  end subroutine zero_set


!
! Private
!==========
! Public
!

  subroutine particles__move( t0, dt )
    real(DR), intent(in) :: t0, dt

    integer(SI) :: n
    type(particles_t), save :: work, d1, d2, d3, d4
    logical, save :: first_time = .true.

    if ( first_time ) then
      call zero_set( work, d1, d2, d3, d4 )
      first_time = .false.
    end if

    !--<step 1>--!
    call equation_of_motion( particles, t0, dt, d1 )
    call runge_kutta( '1st', particles, d1, work )

    !--<step 2>--!
    call equation_of_motion( work, t0+dt/2, dt, d2 )
    call runge_kutta( '2nd', particles, d2, work )

    !--<step 3>--!
    call equation_of_motion( work, t0+dt/2, dt, d3 )
    call runge_kutta( '3rd', particles, d3, work )

    !--<step 4>--!
    call equation_of_motion( work, t0+dt, dt, d4 )
    call runge_kutta( particles, d1, d2, d3, d4 )

    call bounce(particles)
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
    real(DR), parameter :: VEL_MAX = 2.0_DR

    call random_number(rand)
    particles%pos(:)%x = rand(:)
    call random_number(rand)
    particles%vel(:)%x = ( 2*rand(:) - 1.0_DR ) * VEL_MAX

    call random_number(rand)
    particles%pos(:)%y = rand(:)
    call random_number(rand)
    particles%vel(:)%y = ( 2*rand(:) - 1.0_DR ) * VEL_MAX
  end subroutine particles__set

end module particles_m
