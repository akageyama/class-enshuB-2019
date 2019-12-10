!=========================================================
!
!  pbm.f90
!    * Bitmap image file handler in the pbm format.
!
!=========================================================

module pbm_m
  use const_m
  implicit none
  private
  public :: pbm__read,  &
            pbm__revert,  &
            pbm__save

  type, public :: pbm_t
    character(len=2) :: header ! must be 'P1'
    integer(SI) :: width, height
    integer(SI), allocatable :: bitmap(:,:) ! 1 or 0
  end type pbm_t

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


!  private
!=================
!  public


  subroutine pbm__read(image, filename)
    type(pbm_t), intent(out) :: image
    character(len=*), intent(in) :: filename
    integer(SI) :: i, j

    open(FILE_NUM, file=filename, form='formatted')
      read(FILE_NUM,*) image%header
      call assert( image%header == 'P1', &
                   '<pbm__read> '//trim(filename)//': not PBM?' )
      read(FILE_NUM,*) image%width, image%height
      call assert( image%width > 0 .and. image%height > 0, &
                   '<pbm__read> width/height must be positive.' )
      call assert( .not. allocated(image%bitmap), &
                   '<pbm__read> bitmap already allocated.')
      allocate( image%bitmap(image%width,image%height) )
print *,' in the pbm__read, width, height = ', image%width, image%height
      do j = 1 , image%height
        read(FILE_NUM,*) ( image%bitmap(i,j), i=1, image%width )
      end do
    close(FILE_NUM)

    print *,' reading file:',  filename
    print *,' header=',  image%header
    print *,' width=',  image%width, ' height=',image%height
    ! do j = 1 , image%height
    !   write(*, '(i1,1x)') ( image%bitmap(i,j), i=1, image%width )
    ! end do
  end subroutine pbm__read


  subroutine pbm__revert(image)
    type(pbm_t), intent(inout) :: image

    integer(SI) :: i, j
    do j = 1 , image%height
      do i = 1 , image%width
        image%bitmap(i,j) = 1 - image%bitmap(i,j)
      end do
    end do
  end subroutine pbm__revert


  subroutine pbm__save(image, filename)
    type(pbm_t), intent(in) :: image
    character(len=*), intent(in) :: filename
    integer(SI) :: i, j

    open(FILE_NUM, file=filename, form='formatted', status='replace')
      write(FILE_NUM,'(a)')     image%header
      write(FILE_NUM,'(i3,1x,i3)') image%width, image%height
      do j = 1 , image%height
        write(FILE_NUM,'(i1,1x)') ( image%bitmap(i,j), i=1, image%width )
      end do
    close(FILE_NUM)
  end subroutine pbm__save

end module pbm_m
