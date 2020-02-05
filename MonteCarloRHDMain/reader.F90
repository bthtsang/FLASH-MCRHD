module reader
interface read_array_From_file
  subroutine read_array_from_file(myPE, path, arr, rows, cols)
    integer,intent(IN) :: myPE
    character(len=*),intent(IN) :: path
    real,pointer,intent(INOUT) :: arr(:,:)
    integer,intent(OUT) :: rows
    integer,intent(OUT) :: cols
  end subroutine read_array_from_file
end interface read_array_from_file
end module reader

subroutine read_array_from_file(myPE, path, arr, rows, cols)

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"

  integer, intent(IN) :: myPE
  character(len=*), intent(IN) :: path
  real,pointer,intent(INOUT) :: arr(:,:)
  integer, intent(OUT) :: rows
  integer, intent(OUT) :: cols

  integer :: i, status, num_cols, line_len, size, init_size
  real,allocatable,dimension(:,:) :: tmp_data, tmp_data2
  real,allocatable,dimension(:) :: tmp_array
  real,pointer,dimension(:,:) :: ptr_data
  character(len=1024) :: line
  logical :: blank
  
  if (myPE==MASTER_PE) print *, 'Reading from ', trim(path)
  
  !init_size = 250
  init_size = 5000

  open(9, file=path, form='formatted', action='read')
  rewind(9)

  ! read the first line, count columns
  read(9,'(a)',IOSTAT=status), line
  
  line_len = len(line)
  num_cols = 0
  blank = .true.

  do i=1,line_len
  
    if (line(i:i) == ' ') then
      blank = .true.
    else if (blank) then
      num_cols = num_cols + 1
      blank = .false.
    end if
  
  end do
  
  cols = num_cols
  
  rewind(9)
  
  allocate(tmp_data(init_size, num_cols))
  allocate(tmp_array(num_cols))

  size = init_size
  i = 0

  do
    ! check if tmp_data needs to be resized
    ! copy data to tmp_data2, then resize tmp_data and copy back
    if (i >= size) then
      allocate(tmp_data2(size+init_size, num_cols))
      
      tmp_data2 = 0.0
      tmp_data2(1:size,:) = tmp_data(1:size,:)
      
      deallocate(tmp_data)
      allocate(tmp_data(size+init_size, num_cols))
      
      tmp_data = 0.0
      tmp_data(1:size,:) = tmp_data2(1:size,:)
      
      deallocate(tmp_data2)
      
      size = size + init_size
    endif

    ! read array values
    read(9,*,IOSTAT=status), tmp_array(:)
    
    ! end of file, read error
    if (status /= 0) exit
    
    ! we have a good read, increment the count and add it
    i = i + 1
    tmp_data(i,:) = tmp_array(:)

  enddo

  close(9)
  
  rows = i
  
  if (myPE==MASTER_PE) print *, '    : rows = ', rows, ' cols = ', cols
  
  allocate(ptr_data(rows,cols))
  
  ptr_data(:,:) = tmp_data(1:rows,:)
  
  deallocate(tmp_data)
  deallocate(tmp_array)
  
  arr => ptr_data

  return

end subroutine read_array_from_file
