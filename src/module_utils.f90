module utils_mod

contains
 !!Interfaz a "date"
 !function date(date_str, fmt_str) result(output)
 !  implicit none
 !  character(*), intent(in) :: date_str, fmt_str
 !  character(256)           :: command
 !  character(20)            :: output
 !  command="date -d "//trim(date_str)//" '+"//trim(fmt_str)//"'  > tmp_date.txt"
 !  call system( trim(command) )
 !  !print*,trim(command)
 !  open(9, file='tmp_date.txt', status='old',action='read'); read(9, '(A)') output;  close(9)
 !  call system('rm tmp_date.txt')
 !end function

 function atoi(str)     !string -> int
   implicit none
   character(len=*), intent(in) :: str
   integer :: atoi
   read(str,*) atoi
 end function
 function itoa(i)       !int -> string
    implicit none
    integer, intent(in) :: i
    character(len=20) :: itoa
    write(itoa, '(i0)') i
    itoa = adjustl(itoa)
 end function
 function rtoa(r)       !real -> string
    implicit none
    real, intent(in) :: r
    character(len=16) :: rtoa
    write(rtoa, '(F16.3)') r
    rtoa = adjustl(rtoa)
 end function


function Mode(arr) 
  implicit none
  
  integer, intent(in) :: arr(:,:)  ! Input array
  integer :: mode            ! Most frequent value
  integer :: modeCount       ! Count of the most frequent value
  integer :: i, j, N, count
  integer,allocatable :: arr1d(:)

  !Reshape arr to 1-D
  N = product(shape(arr))
  allocate(arr1d(N))
  arr1d = reshape(arr, [N])
        
  ! Initialize the mode and its count to zero
  mode = 0
  modeCount = 0
  
  ! Find the mode
  do i = 1, N
    count = 0
    do j = 1, N
      if (arr1d(j) == arr1d(i)) then
        count = count + 1
      endif
    end do
    
    if (count > modeCount) then
      mode = arr1d(i)
      modeCount = count
    endif
  end do
  
  !! Print the mode and its count
  !write(*,*) "Mode:", mode
  !write(*,*) "Count:", modeCount
  
end function








end module
