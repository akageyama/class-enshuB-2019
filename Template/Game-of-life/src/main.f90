!=========================================================
!
!  class-enshu/template-image/src/main.f90
!
!=========================================================

program main
  use const_m
  use pbm_m
  use gol_m ! game of life
  implicit none

  integer(SI) :: n

  call const__print

  call gol__set_by_image('kobe_bitmap.pbm')
  ! call gol__set_by_image('pentomino_bitmap.pbm')
  ! call gol__set_by_image('start.pbm')
  ! call gol__set_by_program

  call gol__print_summary
  call gol__save

  do n = 1 , 200
    call gol__advance
    call gol__print_summary
    call gol__save
  end do

end program
