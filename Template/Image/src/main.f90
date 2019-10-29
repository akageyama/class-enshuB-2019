!=========================================================
!
!  class-enshu/template-image/src/main.f90
!
!=========================================================

program main
  use const_m
  use ppm_m
  use pbm_m
  implicit none

  type(ppm_t) :: face
  type(pbm_t) :: kobe

  call const__print
  call ppm__read(face,'sample_face.ppm')
  call ppm__convert_grayscale(face)
  call ppm__save(face,'temp_grayscale.ppm')

  call pbm__read(kobe,'kobe_bitmap.pbm')
  call pbm__revert(kobe)
  call pbm__save(kobe,'temp_kobe_rev.pbm')
end program
