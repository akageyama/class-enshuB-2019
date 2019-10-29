!=========================================================
!
!  ppm.f90
!    * Image file handler in the ppm format.
!
!=========================================================

module ppm_m
  use const_m
  implicit none
  private
  public :: ppm__convert_grayscale,  &
            ppm__read,  &
            ppm__save

  type, public :: ppm_t
    character(len=2) :: header ! must be 'P3'
    integer(SI) :: width, height
    integer(SI) :: max ! must be 255
    integer(SI), allocatable :: red(:,:), green(:,:), blue(:,:)
  end type ppm_t

  integer(SI), parameter :: FILE_NUM = 10


contains

  subroutine assert( must_be_true, last_will )
    logical, intent(in) :: must_be_true
    character(len=*), intent(in) :: last_will

    if ( .not. must_be_true ) then
      print *,'*** Error *** <ppm_m>: ' // trim(last_will)
      stop
    end if
  end subroutine assert


!  private
!=================
!  public


  subroutine ppm__read(image, filename)
    class(ppm_t), intent(out) :: image
    character(len=*), intent(in) :: filename
    integer(SI) :: i, j

    open(FILE_NUM, file=filename, form='formatted')
      read(FILE_NUM,*) image%header
      call assert( image%header == 'P3', &
                   trim(filename)//'is not a ppm?' )
      read(FILE_NUM,*) image%width, image%height
      call assert( image%width > 0 .and. image%height > 0, &
                   ': width/height must be positive.' )
      read(FILE_NUM,*) image%max
      call assert( image%max<=999, &
                   'Here we assume 0 <= (r,g,b) <= 999.' )
      allocate( image%red  (image%width,image%height), &
                image%green(image%width,image%height), &
                image%blue (image%width,image%height) )
      do j = 1 , image%height
        read(FILE_NUM,*) ( image%red  (i,j), &
                           image%green(i,j), &
                           image%blue (i,j), i=1, image%width )
      end do
    close(FILE_NUM)

    !> print *,' header=',  image%header
    !> print *,'  width=',  image%width, ' height=',image%height
    !> print *,'   max= ',  image%max
    !> do j = 1, image%height
    !>   do i = 1, image%width
    !>     print *, image%red(i,j), image%green(i,j), image%blue(i,j)
    !>   end do
    !> end do
  end subroutine ppm__read

  subroutine ppm__convert_grayscale(image)
    class(ppm_t), intent(inout) :: image

    integer(SI) :: i, j
    integer(SI) :: gray
    do j = 1 , image%height
      do i = 1 , image%width
        gray = 0.299_SR*image%red  (i,j) &
             + 0.587_SR*image%green(i,j) &
             + 0.114_SR*image%blue (i,j)
        gray = min(gray, image%max) ! just in case...
        image%red  (i,j) = gray
        image%green(i,j) = gray
        image%blue (i,j) = gray
      end do
    end do
  end subroutine ppm__convert_grayscale


  subroutine ppm__save(image, filename)
    class(ppm_t), intent(in) :: image
    character(len=*), intent(in) :: filename
    integer(SI) :: i, j

    open(FILE_NUM, file=filename, form='formatted')
      write(FILE_NUM,'(a)')     image%header
      write(FILE_NUM,'(i3xi3)') image%width, image%height
      write(FILE_NUM,'(i3)')    image%max
      do j = 1 , image%height
        write(FILE_NUM,'(*(i3xi3xi3x))') ( image%red  (i,j), &
                                           image%green(i,j), &
                                           image%blue (i,j), i=1, image%width )
      end do
    close(FILE_NUM)
  end subroutine ppm__save

end module ppm_m
