!=========================================================
!
!  class-enshu/template-image/src/main.f90
!
!=========================================================

program main
  use const_m
  use pgm_m
  use sml_m ! SmoothLife
  implicit none

  integer(SI) :: n

  call const__print

  ! call sml__set_by_image('kobe.pgm')
  call sml__set_by_program
  ! call gol__set_by_image('kobe_bitmap.pbm')
  ! call gol__set_by_image('pentomino_bitmap.pbm')
  ! call gol__set_by_image('start.pbm')
  ! call gol__set_by_program

  call sml__print_summary
  call sml__save


  do n = 1 , 200
    call sml__advance
    call sml__print_summary
    call sml__save
  end do

end program
