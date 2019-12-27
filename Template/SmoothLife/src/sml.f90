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

  integer(SI), parameter :: PAPERS_RA = 21
  integer(SI), parameter :: PAPERS_RI = 7

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

    nbo = PAPERS_RA
    width  = sml%width
    height = sml%height

    call assert( width-2*nbo > 0 .and. height-2*nbo > 0,  &
                "<sml/boundary_condition> nbo is too large" )

    do i = 1, nbo
      sml%f(i,           :) = sml%f(width-2*nbo+i, :)
      sml%f(width-nbo+i, :) = sml%f(        nbo+i, :)
      sml%f(:,           i) = sml%f(:,height-2*nbo+i)
      sml%f(:,height-nbo+i) = sml%f(:,         nbo+i)
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


  function integral( i, j, irad, sheet )
    integer(SI), intent(in) :: i, j, irad
    real(DR), dimension(-irad:irad,-irad:irad), intent(in) :: sheet
    real(DR) :: integral

    integer(SI) :: ii, jj
    real(DR) :: sum_f

    sum_f    = 0.0_DR ! reset

    !      *     |     |     |     |     |
    !      |     *     |     |     |     |
    !      o-----o---*-o-----o-----o-----o
    !      |     |     |     |     |     |
    !      |     |     |  *  |     |     |
    !      o-----o-----o-----o-----o-----o
    !      |     |     |     | *   |     |
    !      |     |     |     |     |     |
    !      o-----o-----o-----o-----o-----o
    !      |     |     |     |     *     |
    !      |     |     |     |     |     |
    !      o-----o-----o-----o-----o-----o
    !      |     |     |     |     |  *  |
    !      |     |     |     |     |     |
    !      X-----o-----o-----o-----o-----*
    !       \                             \
    !        grid=(i,j)                   grid=(i+irad,j)
    !                     b/2     b/2
    !                    __|__   __|__
    !                   /     \ /     \
    !      +-----------x-------O-------x
    !      |                   |       |
    !      |<----------------->|       |
    !      |       irad                |
    !      |                           |
    !      +---------------------------*

    call assert( irad > 0,  &
                 "<sml/integral> irad <= 0 ?!" )
    call assert( irad <= PAPERS_RA,  &
                 "<sml/integral> irad too large." )

    call assert ( i-irad >= 1,          "i-irad out of range" )
    call assert ( i+irad <= sml%width,  "i+irad out of range" )
    call assert ( j-irad >= 1,          "j-irad out of range" )
    call assert ( j+irad <= sml%height, "j+irad out of range" )

    do jj = -irad, +irad
      do ii = -irad, +irad
        sum_f = sum_f + sml%f_copy(i+ii,j+jj) * sheet(ii,jj)
      end do
    end do

    integral = sum_f
  end function integral


  subroutine make_sheet( irad, sheet )
    integer(SI), intent(in) :: irad
    real(DR), dimension(-irad:irad,-irad:irad), intent(out) :: sheet
    
    real(DR), parameter :: PAPERS_B = 1.0_DR
    real(DR) :: papers_ell, rad_plus_b_half, rad_minus_b_half
    integer(SI) :: i, j

    !      *-----o-----o-----o-----o-----o
    !      |   * |     |     |     |     |
    !      |     |     |     |     |     |
    !      o-----o---*-o-----o-----o-----o
    !      |     |     |     |     |     |
    !      |     |     |  *  |     |     |
    !      o-----o-----o-----o-----o-----o
    !      |     |     |     | *   |     |
    !      |     |     |     |     |     |
    !      o-----o-----o-----o-----o-----o
    !      |     |     |     |     |*    |
    !      |     |     |     |     |     |
    !      o-----o-----o-----o-----o-----o
    !      |     |     |     |     |   * |
    !      |     |     |     |     |     |
    !      X-----o-----o-----o-----o-----*
    !       \                             \
    !        center of interest           radius 
    !                     b/2     b/2
    !                    __|__   __|__
    !                   /     \ /     \
    !      +-----------x-------O-------x
    !      |                   |       |
    !      |<----------------->|       |
    !      |       irad                |
    !      |                           |
    !      +---------------------------*


    rad_plus_b_half  = real(irad,DR) + PAPERS_B / 2.0_DR
    rad_minus_b_half = real(irad,DR) - PAPERS_B / 2.0_DR

    do j = -irad, irad
      do i = -irad, irad
        papers_ell = iDistance( i, j )
        if ( papers_ell <= rad_minus_b_half ) then
          sheet(i,j) = 1.0_DR
        else if ( rad_minus_b_half < papers_ell .and. &
                  rad_plus_b_half >= papers_ell ) then
          sheet(i,j) = (  rad_plus_b_half  &
                        - papers_ell ) / PAPERS_B
        else
          sheet(i,j) = 0.0_DR
        end if
      end do
    end do

  contains

    function iDistance( i, j )
      integer(SI), intent(in) :: i, j
      real(DR) :: iDistance

      real(DR) :: x, y

      x = real(i, DR)
      y = real(j, DR)

      iDistance = sqrt( x*x + y*y )
    end function iDistance
    
  end subroutine make_sheet


  function papers_function_sigma1( x, a, alpha )
    real(DR), intent(in) :: x, a, alpha
    real(DR) :: papers_function_sigma1

    real(DR) :: denominator

    denominator = 1.0_DR + exp( -4.0_DR*(x-a)/alpha )
    papers_function_sigma1 = 1.0_DR / denominator
  end function papers_function_sigma1


  function papers_function_sigma2( x, a, b, alpha )
    real(DR), intent(in) :: x, a, b, alpha
    real(DR) :: papers_function_sigma2

    real(DR) :: sigma1a, sigma1b

    sigma1a = papers_function_sigma1( x, a, alpha )
    sigma1b = papers_function_sigma1( x, b, alpha )

    papers_function_sigma2 = sigma1a * ( 1.0_DR - sigma1b )
  end function papers_function_sigma2


  function papers_function_sigma_m( x, y, m )
    real(DR), intent(in) :: x, y, m
    real(DR) :: papers_function_sigma_m

    real(DR) :: sigma1
    real(DR), parameter :: PAPERS_ALPHA_M = 0.147_DR

    sigma1 = papers_function_sigma1( m, 0.5_DR, &
                                     PAPERS_ALPHA_M )

    papers_function_sigma_m = x * ( 1 - sigma1 ) + y * sigma1
  end function papers_function_sigma_m


  function papers_function_s( n, m )
    real(DR), intent(in) :: n, m
    real(DR) :: papers_function_s

    real(DR), parameter :: PAPERS_B1 = 0.278_DR
    real(DR), parameter :: PAPERS_B2 = 0.365_DR
    real(DR), parameter :: PAPERS_D1 = 0.267_DR
    real(DR), parameter :: PAPERS_D2 = 0.445_DR
    real(DR), parameter :: PAPERS_ALPHA_N = 0.028_DR
    real(DR) :: sigma_m_1, sigma_m_2

    sigma_m_1 = papers_function_sigma_m( PAPERS_B1,  &
                                         PAPERS_D1,  &
                                         m )
    sigma_m_2 = papers_function_sigma_m( PAPERS_B2,  &
                                         PAPERS_D2,  &
                                         m )
    papers_function_s = papers_function_sigma2( n,  &
                                                sigma_m_1,  &
                                                sigma_m_2,  &
                                                PAPERS_ALPHA_N )
  end function papers_function_s



!  private
!=================
!  public


  subroutine sml__advance
    integer(SI) :: i, j
    real(DR) :: ra_sq, ri_sq, factor_n, factor_m
    real(DR) :: s, integral_ri, integral_ra
    real(DR) :: papers_m, papers_n
    logical, save :: first_time = .true.
    real(DR), dimension( -PAPERS_RA:PAPERS_RA,  &
                         -PAPERS_RA:PAPERS_RA ), save :: sheet_ra
    real(DR), dimension( -PAPERS_RI:PAPERS_RI,  &
                         -PAPERS_RI:PAPERS_RI ), save :: sheet_ri

    if ( first_time ) then
      call make_sheet( PAPERS_RI, sheet_ri )
      call make_sheet( PAPERS_RA, sheet_ra )
      first_time = .false.
    end if

    sml%f_copy(:,:) = sml%f(:,:)

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

    ri_sq = real(PAPERS_RI, DR)**2
    ra_sq = real(PAPERS_RA, DR)**2

    factor_m = 1.0_DR / ( PI*ri_sq )
    factor_n = 1.0_DR / ( PI*(ra_sq-ri_sq) )

    do j = PAPERS_RA+1, sml%height-PAPERS_RA
      do i = PAPERS_RA+1, sml%width-PAPERS_RA
        integral_ri = integral( i, j, PAPERS_RI, sheet_ri )
        integral_ra = integral( i, j, PAPERS_RA, sheet_ra )
        papers_m = integral_ri * factor_m
        papers_n = ( integral_ra - integral_ri ) * factor_n
        s = papers_function_s( papers_n, papers_m )
        sml%f(i,j) = s
        call assert ( s >= 0.0_DR .and. s <= 1.0_DR,  &
                     "<sml__advance> s out of range.")
      end do
    end do

    call boundary_condition( sml )

    sml%nstep = sml%nstep + 1
  end subroutine sml__advance


  subroutine sml__print_summary
    integer(SI) :: width, height, ngrids
    width  = sml%width
    height = sml%height
    ngrids = width * height
    print *, '(width,height,nstep) : ', width, height, sml%nstep
    print *, ' Average vital level: ', sum(sml%f(:,:)/ngrids)
  end subroutine sml__print_summary


  subroutine sml__set_by_program
    integer(SI) :: width  = 400
    integer(SI) :: height = 400

    integer(SI) :: i, j, i2, j2, skip
    integer(SI) :: some_non_negative_int
    real(DR) :: random

    sml%nstep  = 0
    sml%width  = width
    sml%height = height
    allocate ( sml%f( width, height ) )
    allocate ( sml%f_copy( width, height ) )

    sml%f(:,:) = 0.0_DR  ! default zero

    skip = PAPERS_RA*0.8
    do j = 1 , height, skip
      do i = 1 , width, skip
        call random_number(random)  ! 0.0 to 1.0
        do j2 = 0 , skip-1
          do i2 = 0 , skip-1
            sml%f( i+i2, j+j2 ) = random
          end do
        end do
      end do
    end do

    call boundary_condition( sml )
  end subroutine sml__set_by_program


  subroutine sml__set_by_image( filename )
    character(len=*), intent(in) :: filename

    type(pgm_t) :: pgm
    integer(SI) :: i, j, f_int

    call pgm__read( pgm, filename )

    sml%nstep  = 0
    sml%width  = pgm%width
    sml%height = pgm%height

    allocate( sml%f( sml%width, sml%height ) )
    allocate( sml%f_copy( sml%width, sml%height ) )

    do j = 1 , sml%height
      do i = 1 , sml%width
        f_int = pgm%max - pgm%whitelevel(i,j)
        ! convert white <--> black
        call assert( f_int >= 0 .and. f_int <= pgm%max,  &
                    "<sml__set_by_image> img%f(i,j) out of range." )
        sml%f(i,j) = real(f_int, DR) / pgm%max
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
    integer(SI) :: i, j, width, height
    integer(SI), parameter :: PGM_MAX = 3

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
