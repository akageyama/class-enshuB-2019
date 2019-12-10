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
  implicit none
  private
  public :: sml
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
    real(DR), dimension(:,:), allocatable :: f_copy
  end type sml_t

  type(sml_t) :: sml

  integer(SI), parameter :: FILE_NUM = 10

  integer(SI), parameter :: PAPERS_PARAM_RA_INT = 21
  integer(SI), parameter :: PAPERS_PARAM_RI_INT = 7

  real(DR), parameter :: PAPERS_PARAM_RA = real(PAPERS_PARAM_RA_INT, DR)
  real(DR), parameter :: PAPERS_PARAM_RI = real(PAPERS_PARAM_RI_INT, DR)

  real(DR), parameter :: PAPERS_PARAM_B = 1.0_DR

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
    !           1     2     3     4   ..  n-1    n    n+1   n+2
    !           o-----o-----o-----o-- .. --o-----o-----o-----o--...
    !           |     |     |     |        |     |
    !           |   To implement the periodic boundary condition,
    !           |   we assume that n grids are overlapped.
    !           |     |     |     |        |     |
    !...--o-----o-----o-----o-----o-- .. --o-----o
    !    w-n  w-n+1 w-n+2 w-n+3 w-n+4     w-1    w
    integer(SI) :: i, width, height
    integer(SI) :: nbo ! number of grid points for the bounary overlap

    nbo = PAPERS_PARAM_RA_INT
    width  = sml%width
    height = sml%height

    call assert( nbo <= width .and. nbo <= height, &
                "<sml/boundary_condition> nbo is too large" )

    do i = 1 , nbo
      sml%f(i,           :) = sml%f(width-nbo+i, :)
      sml%f(width-nbo+i, :) = sml%f(          i, :)
      sml%f(:,           i) = sml%f(:,height-nbo+i)
      sml%f(:,height-nbo+i) = sml%f(:,           i)
    end do
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


  function distance(i1, j1, i2, j2)
    integer(SI), intent(in) :: i1, j1, i2, j2

    real(DR) :: dx, dy

    dx = real(i1-i2, DR)
    dy = real(j1-j2, DR)

    distance = sqrt( dx*dx + dy*dy )
  end function distance


  function integrate_circle(i, j, rad)
    integer(SI), intent(in) :: i, j
    real(DR), intent(in) :: rad
    real(DR) :: integrate_circle

    integer(SI) :: ii, jj, delta
    real(DR) :: papers_variable_ell, darea
    real(DR) :: sum_area, sum_graylevel
    real(DR) :: rad_plus_b_half, rad_minus_b_half

    sum_area      = 0.0_DR    ! reset
    sum_graylevel = 0.0_DR

    !      *     |     |     |     |     |
    !      |     *     |     |     |     |
    !      o-----o---*-o-----o-----o-----o
    !      |     |     |     |     |     |
    !      |     |     |  *  |     |     |
    !      o-----o-----o-----o-----o-----o
    !      |     |     |     |*    |     |
    !      |     |     |     |     |     |
    !      o-----o-----o-----o-----o-----o
    !      |     |     |     |    *|     |
    !      |     |     |     |     |     |
    !      o-----o-----o-----o-----o-----o
    !      |     |     |     |     | *   |
    !      |     |     |     |     |     |
    !      X-----o-----o-----o-----o--*--o
    !       \                             \
    !        grid:(i,j)                   grid:(i+delta,j)
    !                     b/2     b/2
    !                    __|__   __|__
    !                   /     \ /     \
    !      +-----------x-------O-------x
    !      |                   |       |
    !      |<----------------->|       |
    !      |        rad                |
    !      |                           |
    !      +---------------------------*

    delta = int(rad+1) ! rad=1.01 -> delta=2
                       ! rad=3.14 -> delta=4
                       ! rad=4.99 -> delta=5

    call assert( rad > 0.0_DR,  &
                "<sml/integrate_circle> rad <= 0?!" )
    call assert ( i-delta >= 1,          "i-delta out of range" )
    call assert ( i+delta <= sml%width,  "i+delta out of range" )
    call assert ( j-delta >= 1,          "j-delta out of range" )
    call assert ( j+delta <= sml%height, "j+delta out of range" )

    rad_plus_b_half  = rad + PAPERS_PARAM_B / 2
    rad_minus_b_half = rad - PAPERS_PARAM_B / 2

    do jj = j-delta, j+delta
      do ii = i-delta, i+delta
        papers_variable_ell = distance(i, j, ii, jj)
        if ( papers_variable_ell <= rad_minus_b_half ) then
          darea = 1.0_DR
        else if ( rad_minus_b_half < papers_variable_ell .and. &
                  rad_plus_b_half >= papers_variable_ell ) then
          darea = (  rad_plus_b_half  &
                   - papers_variable_ell ) / PAPERS_PARAM_B
        else
          darea = 0.0_DR
        end if

        sum_graylevel = sum_graylevel + sml%graylevel_copy(i,j)*darea
        sum_area      = sum_area   + darea
      end do
    end do

    call assert( sum_area > 0.0_DR, &
                "<sml/calc_papers_variable_m> sum_surfae <= 0?!" )

    integrate_circle = sum_graylevel / sum_area

  end function integrate_circle



!  private
!=================
!  public


  subroutine sml__advance
    integer(SI) :: i, j, nbo

    sml%f_copy(:,:) = sml%f(:,:)

    ! Let n = number of boundary overlap
    ! and w = pgm%width
    !           1     2     3     4   ..  n-1    n    n+1   n+2
    !           o-----o-----o-----o-- .. --o-----o-----o-----o--...
    !           |     |     |     |        |     |
    !           |   To implement the periodic boundary condition,
    !           |   we assume that n grids are overlapped.
    !           |     |     |     |        |     |
    !...--o-----o-----o-----o-----o-- .. --o-----o
    !    w-n  w-n+1 w-n+2 w-n+3 w-n+4     w-1    w

    nbo = PAPERS_PARAM_RA_INT
    do j = nbo+1, sml%height-nbo
      do i = nbo+1, sml%width-nbo
        papers_variable_m = integrate_circle( i, j, PAPERS_PARAM_RI )
        papers_bariable_n = integrate_circle( i, j, PAPERS_PARAM_RA )  &
                          - papers_variable_m
        sml%f(i,j) = papers_function_s( papers_variable_n, papers_variable_m )
      end do
    end do

    call boundary_condition( sml )

    sml%nstep = sml%nstep + 1
  end subroutine sml__advance


  subroutine sml__print_summary
    print *, '(width,height) : ', sml%width, sml%height
    print *, ' nstep: ', sml%nstep
    print *, ' Num of alive cell: ', sum(sml%pgm%graylevel(:,:))
  end subroutine sml__print_summary


  subroutine sml__set_by_program
    integer(SI) :: width  = 101
    integer(SI) :: height = 71

    integer(SI) :: i, j, i2, j2
    integer(SI) :: some_non_negative_int
    real(DR) :: random

    sml%nstep  = 0
    sml%width  = width
    sml%height = height
    allocate ( sml%f( width, height ) )
    allocate ( sml%f_copy( width, height ) )

    sml%f(:,:) = 0.0_DR  ! default zero

    do j = 1 , height
      do i = 1 , width
        call random_number(random)
          sml%f( i, j ) = random
        end if
      end do
    end do

    sml%width  = width
    sml%height = height

    call boundary_condition( sml )
  end subroutine sml__set_by_program


  subroutine sml__set_by_image( filename )
    character(len=*), intent(in) :: filename

    type(pgm_t) :: pgm
    integer(SI) :: gray

    call pgm__read( pgm, filename )

    sml%nstep  = 0
    sml%width  = pgm%width
    sml%height = pgm%height

    allocate( sml%f( sml%width, sml%height ) )
    allocate( sml%f_copy( sml%width, sml%height ) )

    do j = 1 , sml%height
      do i = 1 , sml%width
        gray = pgm%graylevel(i,j)
        call assert( gray >= 0 .and. gray <= pgm%max,  &
        sml%f(i,j) = real(gray, DR) / pgm%max
      end do
    end do

    call boundary_condition( sml )
  end subroutine sml__set_by_image


  subroutine sml__reset
    sml%nstep = 0
    sml%f(:,:) = 0.0_DR
  end subroutine sml__reset

  subroutine sml__save
    type(pgm_t) :: pgm
    integer(SI) :: widht, height
    integer(SI), parameter :: PGM_MAX = 255

    width  = sml%width
    height = sml%height

    pgm%width  = width
    pgm%height = height
    pgm%max    = PGM_MAX

    allocate( pgm%graylevel(width, height) )

    do j = 1, height
      do i = 1 , width
        pgm%grayscale(i,j) = int(sml%f(i,j)*PGM_MAX)
      end do
    end do

    call pgm__save( pgm,  &
                    '../data/' // int_to_str6(sml%nstep) // '.pgm' )
  end subroutine sml__save

end module sml_m
