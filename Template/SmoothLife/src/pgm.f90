!=========================================================
!
!  pgm.f90
!    * Grayscale image file handler in the pgm format.
!
!=========================================================

module pgm_m
  use const_m
  implicit none
  private
  public :: pgm__read,  &
            pgm__revert,  &
            pgm__save

  type, public :: pgm_t
    character(len=2) :: header ! must be 'P2'
    integer(SI) :: width, height
    integer(SI) :: max
    integer(SI), allocatable :: whitelevel(:,:) ! [0...max]
  end type pgm_t

  integer(SI), parameter :: FILE_NUM = 10


contains

  subroutine assert( must_be_true, last_will )
    logical, intent(in) :: must_be_true
    character(len=*), intent(in) :: last_will

    if ( .not. must_be_true ) then
      print *,'*** Error *** <pgm_m>: ' // trim(last_will)
      stop
    end if
  end subroutine assert


!  private
!=================
!  public


  subroutine pgm__read( image, filename )
    type(pgm_t), intent(out) :: image
    character(len=*), intent(in) :: filename
    integer(SI) :: i, j

    open(FILE_NUM, file=filename, form='formatted')
      read(FILE_NUM,*) image%header
      call assert( image%header == 'P2', &
                   '<pgm__read> '//trim(filename)//': not pgm?' )
      read(FILE_NUM,*) image%width, image%height
      read(FILE_NUM,*) image%max
      call assert( image%width > 0 .and. image%height > 0, &
                   '<pgm__read> width/height must be positive.' )
      call assert( image%max > 0, &
                   '<pgm__read> max must be positive.' )
      if ( allocated( image%whitelevel ) ) then
        deallocate( image%whitelevel )
      end if
      allocate( image%whitelevel( image%width, image%height ) )
      do j = 1 , image%height
        read(FILE_NUM,*) ( image%whitelevel(i,j), i=1, image%width )
      end do
    close(FILE_NUM)

    print *,' reading a pgm file:',  filename
    print *,' header=',  image%header
    print *,' width=',  image%width, ' height=',image%height
    print *,'    max=',  image%max
    ! do j = 1 , image%height
    !   write(*, '(i1,1x)') ( image%whitelevel(i,j), i=1, image%width )
    ! end do
  end subroutine pgm__read


  subroutine pgm__revert(image)
    type(pgm_t), intent(inout) :: image

    integer(SI) :: i, j
    do j = 1 , image%height
      do i = 1 , image%width
        image%whitelevel(i,j) = image%max - image%whitelevel(i,j)
      end do
    end do
  end subroutine pgm__revert


  subroutine pgm__save( image, filename )
    type(pgm_t), intent(in) :: image
    character(len=*), intent(in) :: filename
    integer(SI) :: i, j

    open(FILE_NUM, file=filename, form='formatted', status='replace')
      write(FILE_NUM,'(a)')        image%header
      write(FILE_NUM,'(i5,1x,i5)') image%width, image%height
      write(FILE_NUM,'(i5)')       image%max
      do j = 1 , image%height
        !write(FILE_NUM,'(i3,1x)') ( image%whitelevel(i,j), i=1, image%width )
        do i = 1 , image%width
          write(FILE_NUM,'(i3,1x)') image%whitelevel(i,j)
        end do
      end do
    close(FILE_NUM)
  end subroutine pgm__save

end module pgm_m
