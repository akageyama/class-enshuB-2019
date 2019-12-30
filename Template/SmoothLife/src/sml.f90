!=========================================================
!
!  sml.f90
!    * SmoothLife
!    * Reference:
!        "Generalization of Conway's "Game of Life" to a 
!         continuous domain - SmoothLife" by Stephan Rafler
!         https://arxiv.org/abs/1111.1567
!
!=========================================================

module sml_m
  use const_m
  use pgm_m
  !$use omp_lib
  implicit none
  private
  public :: sml__advance,  &
            sml__print_summary,  &
            sml__reset,  &
            sml__save,  &
            sml__set_by_image,  &
            sml__set_by_program

  type, public :: sml_t
    integer(SI) :: nstep
    integer(SI) :: width, height
    real(DR), dimension(:,:), allocatable :: f
  end type sml_t

  integer(SI), parameter :: FILE_NUM = 10

contains

  subroutine assert( must_be_true, last_will )
    logical, intent(in) :: must_be_true
    character(len=*), intent(in) :: last_will

    if ( .not. must_be_true ) then
      print *,'*** Error *** <sml_m>: ' // trim(last_will)
      stop
    end if
  end subroutine assert


  subroutine boundary_condition(sml)
    type(sml_t), intent(inout) :: sml

    ! Let n = number of boundary overlap
    ! and w = pgm%width
    !
    ! When n=3
    !           1       2       3       4
    !           o-------o-------o-------o--...
    !           |       |       |
    !           |       |       |
    !           |       |       |
    !...-o------o-------o-------o-------o-------o-------o
    !   w-6    w-5     w-4     w-3     w-2     w-1      w
    integer(SI) :: i, width, height
    integer(SI) :: nbo ! number of grid points for the bounary overlap
    integer(SI) :: j

    width  = sml%width
    height = sml%height

    sml%f(    1,     :) = sml%f(width-1,        :)
    sml%f(width,     :) = sml%f(      2,        :)
    sml%f(    :,     1) = sml%f(      :, height-1)
    sml%f(    :,height) = sml%f(      :,        2)
  end subroutine boundary_condition


  function int_to_str6(i) result(str)
    integer(SI), intent(in) :: i
    character(len=6) :: str
    !===============
    !  Convert an integer into 6 characters.
    !          i =   10 --> str="000010"
    !          i =  123 --> str="000123"
    !          i = -123 --> str="XXXXXX"
    !===============
    if ( i < 0 .or. i > 999999 ) then
      str = 'XXXXXX'
    else
      write(str, '(i6.6)') i
    end if
  end function int_to_str6


  function double_sigmoid( x, a, b )
    real(DR), intent(in) :: x, a, b
    real(DR) :: double_sigmoid

    real(DR), parameter :: DELTA = 0.41_DR
    real(DR), parameter :: DELTA_INV = 1.0_DR / DELTA

    !                oooooooo
    !               o|      |o
    !              o |      | o
    !    oooooooooo--+------+--oooooooooo--> x
    !             |  |      |  |
    !            a1  a2    b1  b2

    real(DR) :: a1, a2, b1, b2

    a1 = a - DELTA/2
    a2 = a + DELTA/2
    b1 = b - DELTA/2
    b2 = b + DELTA/2

    if ( x < a1 ) then
      double_sigmoid = 0.0_DR
    else if ( x < a2 ) then
      double_sigmoid = ( x - a1 ) * DELTA_INV
    else if ( x < b1 ) then
      double_sigmoid = 1.0_DR
    else if ( x < b2 ) then
      double_sigmoid = 1.0_DR - ( x - b1 ) * DELTA_INV
    else
      double_sigmoid = 0.0_DR
    end if
  end function double_sigmoid



!  private
!=================
!  public


  subroutine sml__advance(sml)
    type(sml_t), intent(inout) :: sml

    integer(SI) :: i, j
    real(DR) :: s, df
    type(sml_t), save :: sml_copy

    sml_copy = sml

    ! Let n = number of boundary overlap
    ! and w = pgm%width
    !
    ! When n=3
    !           1       2       3       4
    !           o-------o-------o-------o--...
    !           |       |       |
    !           |       |       |
    !           |       |       |
    !...-o------o-------o-------o-------o-------o-------o
    !   w-6    w-5     w-4     w-3     w-2     w-1      w

    !$omp parallel do
    do j = 2, sml%height-1
      do i = 2, sml%width-1
        df = ( sml_copy%f(i-1,j-1) + sml_copy%f(i  ,j-1)  &
             + sml_copy%f(i+1,j-1) + sml_copy%f(i+1,j  )  &
             + sml_copy%f(i+1,j+1) + sml_copy%f(i  ,j+1)  &
             + sml_copy%f(i-1,j+1) + sml_copy%f(i-1,j  ) )  &
             - ( 9 * sml_copy%f(i,j) )
        if ( df >=0.0_DR ) then
          s = double_sigmoid( df,  2.40_DR ,  3.60_DR )
        else
          s = double_sigmoid( df, -7.60_DR , -5.40_DR )
        end if
        sml%f(i,j) = 0.05_DR * sml%f(i,j) + 0.95_DR*s ! Vertical, horizontal,
                                                      ! and diagonal gliders
                                                      ! odd/even oscillation.
                                                      ! Covers the whole space.
        !=sml%f(i,j) = 0.10_DR * sml%f(i,j) + 0.90_DR*s  ! Vertical glider in some
        !=                                               ! cases, and vortex-like
        !=                                               ! object.
        !!!!!sml%f(i,j) = 0.15_DR * sml%f(i,j) + 0.85_DR*s ! Vertical glider
        !!!!sml%f(i,j) = 0.17_DR * sml%f(i,j) + 0.83_DR*s ! Vertical glider
        !!! sml%f(i,j) = 0.20_DR * sml%f(i,j) + 0.80_DR*s ! Vertical glider

        call assert ( sml%f(i,j) >= 0.0_DR .and. sml%f(i,j) <= 1.0_DR,  &
                     "<sml__advance> sml%f(i,j) out of range.")
      end do
    end do
    !$omp end parallel do

    call boundary_condition( sml )

    sml%nstep = sml%nstep + 1
  end subroutine sml__advance


  subroutine sml__print_summary( sml )
    type(sml_t), intent(in) :: sml

    integer(SI) :: width, height, ngrids
    width  = sml%width
    height = sml%height
    ngrids = width * height
    print *, '(width,height,nstep) : ', width, height, sml%nstep
    print *, ' Average vital level: ', sum(sml%f(:,:)/ngrids)
  end subroutine sml__print_summary


  subroutine sml__set_by_program( sml )
    type(sml_t), intent(out) :: sml
    integer(SI) :: width  = 140
    integer(SI) :: height = 140

    integer(SI) :: i, j, i2, j2, skip
    integer(SI) :: some_non_negative_int
    real(DR) :: random

    sml%nstep  = 0
    sml%width  = width
    sml%height = height
    allocate ( sml%f( width, height ) )

    sml%f(:,:) = 0.0_DR  ! default zero

    skip = 1
    do j = 1 , height, skip
      do i = 1 , width, skip
        call random_number(random)  ! 0.0 to 1.0
        do j2 = 0 , skip-1
          do i2 = 0 , skip-1
            if ( .not. ( i+i2 >= 0.4*sml%width .and.  &
                         i+i2 <= 0.6*sml%width .and.  &
                         j+j2 >= 0.4*sml%height .and.  &
                         j+j2 <= 0.6*sml%height ) ) cycle
            if ( i+i2 <= sml%width .and. &
                 j+j2 <= sml%height ) then
              sml%f( i+i2, j+j2 ) = random
            end if
          end do
        end do
      end do
    end do

    call boundary_condition( sml )
  end subroutine sml__set_by_program


  subroutine sml__set_by_image( sml, filename )
    character(len=*), intent(in) :: filename
    type(sml_t), intent(out) :: sml

    type(pgm_t) :: pgm
    integer(SI) :: i, j, f_int

    call pgm__read( pgm, filename )

    sml%nstep  = 0
    sml%width  = pgm%width
    sml%height = pgm%height

    allocate( sml%f( sml%width, sml%height ) )

    do j = 1 , sml%height
      do i = 1 , sml%width
        ! f_int = pgm%max - pgm%whitelevel(i,j) ! conv white <-> black
        f_int = pgm%whitelevel(i,j)
        call assert( f_int >= 0 .and. f_int <= pgm%max,  &
                    "<sml__set_by_image> img%f(i,j) out of range." )
        sml%f(i,j) = real(f_int, DR) / pgm%max
      end do
    end do

    call boundary_condition( sml )
  end subroutine sml__set_by_image


  subroutine sml__reset( sml )
    type(sml_t), intent(out) :: sml
    sml%nstep = 0
    sml%f(:,:) = 0.0_DR
  end subroutine sml__reset


  subroutine sml__save( sml )
    type(sml_t), intent(in) :: sml

    type(pgm_t) :: pgm
    integer(SI) :: i, j, width, height
    integer(SI), parameter :: PGM_MAX = 100

    width  = sml%width
    height = sml%height

    pgm%header = 'P2'
    pgm%width  = width
    pgm%height = height
    pgm%max    = PGM_MAX

    allocate( pgm%whitelevel(width, height) )

    do j = 1, height
      do i = 1 , width
        pgm%whitelevel(i,j) = int( sml%f(i,j)*PGM_MAX )
      end do
    end do

    call pgm__save( pgm,  &
                    '../data/' // int_to_str6(sml%nstep) // '.pgm' )
  end subroutine sml__save

end module sml_m
