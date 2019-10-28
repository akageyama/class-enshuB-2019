!=========================================================
!
!  rk.f90
!    * The 4-th order 4-step Runge-Kutta integrator
!
!=========================================================

module rk_m
  use const_m
  use particles_m
  implicit none
  private
  public :: rk__integrate

  interface rk__integrate
     module procedure step_01_02_03, &
                      step_04
  end interface


contains

  subroutine step_01_02_03( step, ps, dps, rk )
    character(len=3), intent(in) :: step
    type(particles_t), intent(in) :: ps, dps
    type(particles_t), intent(out) :: rk
!
!   * 1st to 3rd steps in 4-stage-Runge-Kutta method.
!
    integer(SI) :: n

    if ( step=='1st' .or. step=='2nd' ) then
      do n = 1, p%get_n_particles()
        rk%pos%x(n) = ps%pos%x(n) + 0.5_DR*dps%pos%x(n)
        rk%pos%y(n) = ps%pos%y(n) + 0.5_DR*dps%pos%y(n)
        rk%vel%x(n) = ps%vel%x(n) + 0.5_DR*dps%vel%x(n)
        rk%vel%y(n) = ps%vel%y(n) + 0.5_DR*dps%vel%y(n)
      end do
    else if ( step=='3rd' ) then
      do n = 1, p%get_n_particles()
        rk%pos%x(n) = ps%pos%x(n) + dps%pos%x(n)
        rk%pos%y(n) = ps%pos%y(n) + dps%pos%y(n)
        rk%vel%x(n) = ps%vel%x(n) + dps%vel%x(n)
        rk%vel%y(n) = ps%vel%y(n) + dps%vel%y(n)
      end do
    else
      print *,'<rk/step_01_02_03> *** ERROR *** bad arg.'
      stop
    end if

  end subroutine step_01_02_03


  subroutine step_04( ps, dps1, dps2, dps3, dps4 )
    type(particles_t), intent(inout) :: ps
    type(particles_t), intent(in) :: dps1, dps2, dps3, dps4
!
!   * The final step of the 4-th order Runge-Kutta method.
!
    real(DR), parameter :: ONE6TH = 1.0_DR / 6.0_DR
    real(DR), parameter :: ONE3RD = 1.0_DR / 3.0_DR

    integer(SI) :: n

    do n = 1, parts%get_n_particles()
      ps%pos(n)%x = ps%pos(n)%x  &
                  + ONE6TH*( dps1%pos(n)%x &
                           + dps4%pos(n)%x )  &
                  + ONE3RD*( dps2%pos(n)%x &
                           + dps3%pos(n)%x )
      ps%pos(n)%y = ps%pos(n)%y  &
                  + ONE6TH*( dps1%pos(n)%y &
                           + dps4%pos(n)%y )  &
                  + ONE3RD*( dps2%pos(n)%y &
                           + dps3%pos(n)%y )
      ps%vel(n)%x = ps%vel(n)%x  &
                  + ONE6TH*( dps1%vel(n)%x &
                           + dps4%vel(n)%x )  &
                  + ONE3RD*( dps2%vel(n)%x &
                           + dps3%vel(n)%x )
      ps%vel(n)%y = ps%vel(n)%y  &
                  + ONE6TH*( dps1%vel(n)%y &
                           + dps4%vel(n)%y )  &
                  + ONE3RD*( dps2%vel(n)%y &
                           + dps3%vel(n)%y )
    end do

  end subroutine step_04

end module rk_m
