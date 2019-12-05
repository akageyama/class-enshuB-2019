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
            sml__revert,  &
            sml__save,  &
            sml__set_by_image,  &
            sml__set_by_program

  type, public :: sml_t
    type(pgm_t) :: pgm
    type(pgm_t) :: pgm_copy
    integer(SI) :: nstep
    integer(SI) :: max
    integer(SI) :: width, height
    integer, dimension(:,:), allocatable :: graylevel_copy
  end type sml_t

  type(sml_t) :: sml

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


  subroutine boundary_condition(pgm)
    type(pgm_t), intent(inout) :: pgm
    !  0 0 0 0 0 0 0 0
    !  0 * * * * * * 0
    !  0 * * * * * * 0
    !  0 * * * * * * 0
    !  0 * * * * * * 0
    !  0 0 0 0 0 0 0 0
    pgm%graylevel(        1, :) = 0
    pgm%graylevel(pgm%width, :) = 0
    pgm%graylevel(:,         1) = 0
    pgm%graylevel(:,pgm%height) = 0
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


!  private
!=================
!  public


  subroutine sml__advance

    integer(SI) :: i, j, ii, jj
    integer(SI) :: neighbours

    sml%graylevel_copy(:,:) = sml%pgm%graylevel(:,:)

    call assert( maxval( sml%pgm%graylevel(:,:) ) <= 1 .and.  &
                 minval( sml%pgm%graylevel(:,:) ) >= 0,  &
                 'sml%pgm%graylevel value out of range.' )

    do j = 2, sml%height-1
      do i = 2, sml%width-1
        ! 1st step: Count alive cells near the target cell
        neighbours = 0  ! Number of alive ones.
        do jj = j-1, j+1
          do ii = i-1, i+1
            if (  (ii/=i) .or. (jj/=j) ) then
              if ( sml%graylevel_copy(ii,jj)==1 ) then
                neighbours = neighbours + 1
              end if
            end if
          end do
        end do

        ! 2nd step: Change the target cell state
        if ( sml%graylevel_copy(i,j)==1 ) then ! The cell is alive.
          if ( neighbours < 2 .or. neighbours > 3 ) then
            sml%pgm%graylevel(i,j) = 0  ! Die unless it has 2 or 3 neighbours
          end if
        else  ! The cell is dead.
          if ( neighbours == 3 ) then
            sml%pgm%graylevel(i,j) = 1  ! Reborn if it has 3 neighbours
          end if
        end if
      end do
    end do

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
    integer(SI) :: max = 255

    integer(SI) :: i, j, i2, j2
    integer(SI) :: some_non_negative_int
    real(DR) :: random

    sml%pgm%header = 'P2'
    sml%pgm%max    = max
    sml%pgm%width  = width
    sml%pgm%height = height
    allocate ( sml%pgm%graylevel( width, height ) )

    sml%pgm%graylevel(:,:) = 0  ! default zero

    do j = 1 , height
      do i = 1 , width
        call random_number(random)
        if ( random > 0.9_DR ) then
          sml%pgm%graylevel( i, j ) = 1
        end if
      end do
    end do

    allocate( sml%graylevel_copy( width, height ) )

    sml%width  = width
    sml%height = height
    sml%max    = max

    call boundary_condition( sml%pgm )
  end subroutine sml__set_by_program


  subroutine sml__set_by_image(filename)
    character(len=*), intent(in) :: filename

    call pgm__read( sml%pgm, filename )

    sml%width  = sml%pgm%width
    sml%height = sml%pgm%height
    sml%max    = sml%pgm%max

    allocate( sml%graylevel_copy( sml%width, sml%height ) )

    call boundary_condition(sml%pgm)
  end subroutine sml__set_by_image


  subroutine sml__reset
    sml%nstep = 0
    sml%pgm%graylevel(:,:) = 0
  end subroutine sml__reset


  subroutine sml__revert
    call pgm__revert( sml%pgm )
    call boundary_condition(sml%pgm)
  end subroutine sml__revert


  subroutine sml__save
    call pgm__save( sml%pgm,  &
                    '../data/' // int_to_str6(sml%nstep) // '.pgm' )
  end subroutine sml__save

end module sml_m
