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
  type(sml_t) :: sml

  ! call const__print

  !call sml__set_by_image( sml, 'sample_face.pgm' )
  call sml__set_by_program( sml )

  call sml__print_summary( sml )
  call sml__save( sml )

  do n = 1 , 2000
    call sml__advance( sml )
    call sml__print_summary( sml )
    call sml__save( sml )
  end do

end program
