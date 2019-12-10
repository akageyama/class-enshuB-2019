!=========================================================
!
!  gol.f90
!    * Game of Life
!
!=========================================================

module gol_m
  use const_m
  use pbm_m
  implicit none
  private
  public :: gol
  public :: gol__advance,  &
            gol__print_summary,  &
            gol__reset,  &
            gol__revert,  &
            gol__save,  &
            gol__set_by_image,  &
            gol__set_by_program

  type, public :: gol_t
    type(pbm_t) :: pbm
    type(pbm_t) :: pbm_copy
    integer(SI) :: nstep
    integer(SI) :: width, height
    integer, dimension(:,:), allocatable :: bitmap_copy
  end type gol_t

  type(gol_t) :: gol

  integer(SI), parameter :: FILE_NUM = 10


contains

  subroutine assert( must_be_true, last_will )
    logical, intent(in) :: must_be_true
    character(len=*), intent(in) :: last_will

    if ( .not. must_be_true ) then
      print *,'*** Error *** <pbm_m>: ' // trim(last_will)
      stop
    end if
  end subroutine assert


  subroutine boundary_condition(pbm)
    type(pbm_t), intent(inout) :: pbm
    !  0 0 0 0 0 0 0 0
    !  0 * * * * * * 0
    !  0 * * * * * * 0
    !  0 * * * * * * 0
    !  0 * * * * * * 0
    !  0 0 0 0 0 0 0 0
    pbm%bitmap(        1, :) = 0
    pbm%bitmap(pbm%width, :) = 0
    pbm%bitmap(:,         1) = 0
    pbm%bitmap(:,pbm%height) = 0
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


  subroutine gol__advance

    integer(SI) :: i, j, ii, jj
    integer(SI) :: neighbours

    gol%bitmap_copy(:,:) = gol%pbm%bitmap(:,:)

    call assert( maxval( gol%pbm%bitmap(:,:) ) <= 1 .and.  &
                 minval( gol%pbm%bitmap(:,:) ) >= 0,  &
                 'gol%pbm%bitmap value out of range.' )

    do j = 2, gol%height-1
      do i = 2, gol%width-1
        ! 1st step: Count alive cells near the target cell
        neighbours = 0  ! Number of alive ones.
        do jj = j-1, j+1
          do ii = i-1, i+1
            if (  (ii/=i) .or. (jj/=j) ) then
              if ( gol%bitmap_copy(ii,jj)==1 ) then
                neighbours = neighbours + 1
              end if
            end if
          end do
        end do

        ! 2nd step: Change the target cell state
        if ( gol%bitmap_copy(i,j)==1 ) then ! The cell is alive.
          if ( neighbours < 2 .or. neighbours > 3 ) then
            gol%pbm%bitmap(i,j) = 0  ! Die unless it has 2 or 3 neighbours
          end if
        else  ! The cell is dead.
          if ( neighbours == 3 ) then
            gol%pbm%bitmap(i,j) = 1  ! Reborn if it has 3 neighbours
          end if
        end if
      end do
    end do

    gol%nstep = gol%nstep + 1
  end subroutine gol__advance


  subroutine gol__print_summary
    print *, '(width,height) : ', gol%width, gol%height
    print *, ' nstep: ', gol%nstep
    print *, ' Num of alive cell: ', sum(gol%pbm%bitmap(:,:))
  end subroutine gol__print_summary


  subroutine gol__set_by_program
    integer(SI) :: width  = 101
    integer(SI) :: height = 71

    integer(SI) :: i, j, i2, j2
    integer(SI) :: some_non_negative_int
    real(DR) :: random

    gol%pbm%header = 'P1'
    gol%pbm%width  = width
    gol%pbm%height = height
    allocate ( gol%pbm%bitmap( width, height ) )

    gol%pbm%bitmap(:,:) = 0  ! default zero

    do j = 1 , height
      do i = 1 , width
        call random_number(random)
        if ( random > 0.9_DR ) then
          gol%pbm%bitmap( i, j ) = 1
        end if
      end do
    end do

    allocate( gol%bitmap_copy( width, height ) )

    gol%width  = width
    gol%height = height

    call boundary_condition( gol%pbm )
  end subroutine gol__set_by_program


  subroutine gol__set_by_image(filename)
    character(len=*), intent(in) :: filename

    call pbm__read( gol%pbm, filename )

    gol%width  = gol%pbm%width
    gol%height = gol%pbm%height

    allocate( gol%bitmap_copy( gol%width, gol%height ) )

    call boundary_condition(gol%pbm)
  end subroutine gol__set_by_image


  subroutine gol__reset
    gol%nstep = 0
    gol%pbm%bitmap(:,:) = 0
  end subroutine gol__reset


  subroutine gol__revert
    call pbm__revert( gol%pbm )
    call boundary_condition(gol%pbm)
  end subroutine gol__revert


  subroutine gol__save
    call pbm__save( gol%pbm,  &
                    '../data/' // int_to_str6(gol%nstep) // '.pbm' )
  end subroutine gol__save

end module gol_m
