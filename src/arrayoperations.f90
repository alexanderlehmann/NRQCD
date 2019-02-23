!>@brief Contains various array operations, like sorting and removing duplicates
!!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!!@date 17.02.2019
!!@version 1.0
module arrayoperations
  use, intrinsic :: iso_fortran_env
  implicit none

  PRIVATE

  public :: &
       RemoveDuplicates,&
       Sort

  !>@brief Removes duplicates from 1-dimensional array containing integers, real or complex numbers
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
  interface RemoveDuplicates
     module procedure RemoveDuplicates_int8
     module procedure RemoveDuplicates_int16
     module procedure RemoveDuplicates_int32
     module procedure RemoveDuplicates_int64

     module procedure RemoveDuplicates_real32
     module procedure RemoveDuplicates_real64
     module procedure RemoveDuplicates_real128

     module procedure RemoveDuplicates_complex32
     module procedure RemoveDuplicates_complex64
     module procedure RemoveDuplicates_complex128
  end interface RemoveDuplicates

  !>@brief Sorts 1-dimensional integer or real array
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.02.2019
  !!@version 1.0
  interface Sort
     module procedure sort_int8
     module procedure sort_int16
     module procedure sort_int32
     module procedure sort_int64

     module procedure sort_real32
     module procedure sort_real64
     module procedure sort_real128
  end interface Sort
contains

  pure subroutine RemoveDuplicates_int8(array,nduplicates)
    use, intrinsic :: iso_fortran_env
    implicit none
    integer(int8), parameter :: kind=int8
    integer(kind), allocatable, intent(inout) :: array(:)
    integer(int64),allocatable, intent(out), optional :: nduplicates(:)

    integer(kind), allocatable :: reduced_array(:)

    integer(kind), parameter :: empty=huge(1_kind)
    
    integer(int64) :: i, original_size, reduced_size, j
    integer(int64), allocatable :: countduplicates(:)
    
    original_size = size(array)
    allocate(reduced_array(original_size))
    reduced_array = empty

    if(present(nduplicates)) then
       allocate(countduplicates(original_size))
       countduplicates = 0
    end if
    
    reduced_size=0
    do i=1,original_size
       ! if the entry already exists, check next entry
       if(any(reduced_array == array(i))) then
          if(present(nduplicates)) then
             j = FindLoc(reduced_array,DIM=1,VALUE=array(i))
             countduplicates(j) = countduplicates(j) + 1_int64
          end if
          cycle
       else
          ! no match found, therefore add it to the output
          reduced_size = reduced_size + 1_int64

          reduced_array(reduced_size) = array(i)
          if(present(nduplicates)) then
             countduplicates(reduced_size) = 1_int64
          end if
       end if
    end do

    ! Now reallocate input-array and assign reduced content
    deallocate(array)
    allocate(array(reduced_size))
    array = reduced_array(1:reduced_size)

    if(present(nduplicates)) then
       allocate(nduplicates(reduced_size))
       nduplicates = countduplicates(1:reduced_size)
    end if
  end subroutine RemoveDuplicates_Int8
  
  pure subroutine RemoveDuplicates_int16(array,nduplicates)
    use, intrinsic :: iso_fortran_env
    implicit none
    integer(int8), parameter :: kind=int16
    integer(kind), allocatable, intent(inout) :: array(:)
    integer(int64),allocatable, intent(out), optional :: nduplicates(:)

    integer(kind), allocatable :: reduced_array(:)

    integer(kind), parameter :: empty=huge(1_kind)
    
    integer(int64) :: i, original_size, reduced_size, j
    integer(int64), allocatable :: countduplicates(:)
    
    original_size = size(array)
    allocate(reduced_array(original_size))
    reduced_array = empty

    if(present(nduplicates)) then
       allocate(countduplicates(original_size))
       countduplicates = 0
    end if
    
    reduced_size=0
    do i=1,original_size
       ! if the entry already exists, check next entry
       if(any(reduced_array == array(i))) then
          if(present(nduplicates)) then
             j = FindLoc(reduced_array,DIM=1,VALUE=array(i))
             countduplicates(j) = countduplicates(j) + 1_int64
          end if
          cycle
       else
          ! no match found, therefore add it to the output
          reduced_size = reduced_size + 1_int64

          reduced_array(reduced_size) = array(i)
          if(present(nduplicates)) then
             countduplicates(reduced_size) = 1_int64
          end if
       end if
    end do

    ! Now reallocate input-array and assign reduced content
    deallocate(array)
    allocate(array(reduced_size))
    array = reduced_array(1:reduced_size)

    if(present(nduplicates)) then
       allocate(nduplicates(reduced_size))
       nduplicates = countduplicates(1:reduced_size)
    end if
  end subroutine RemoveDuplicates_Int16
  
  pure subroutine RemoveDuplicates_int32(array,nduplicates)
    use, intrinsic :: iso_fortran_env
    implicit none
    integer(int8), parameter :: kind=int32
    integer(kind), allocatable, intent(inout) :: array(:)
    integer(int64),allocatable, intent(out), optional :: nduplicates(:)

    integer(kind), allocatable :: reduced_array(:)

    integer(kind), parameter :: empty=huge(1_kind)
    
    integer(int64) :: i, original_size, reduced_size, j
    integer(int64), allocatable :: countduplicates(:)
    
    original_size = size(array)
    allocate(reduced_array(original_size))
    reduced_array = empty

    if(present(nduplicates)) then
       allocate(countduplicates(original_size))
       countduplicates = 0
    end if
    
    reduced_size=0
    do i=1,original_size
       ! if the entry already exists, check next entry
       if(any(reduced_array == array(i))) then
          if(present(nduplicates)) then
             j = FindLoc(reduced_array,DIM=1,VALUE=array(i))
             countduplicates(j) = countduplicates(j) + 1_int64
          end if
          cycle
       else
          ! no match found, therefore add it to the output
          reduced_size = reduced_size + 1_int64

          reduced_array(reduced_size) = array(i)
          if(present(nduplicates)) then
             countduplicates(reduced_size) = 1_int64
          end if
       end if
    end do

    ! Now reallocate input-array and assign reduced content
    deallocate(array)
    allocate(array(reduced_size))
    array = reduced_array(1:reduced_size)

    if(present(nduplicates)) then
       allocate(nduplicates(reduced_size))
       nduplicates = countduplicates(1:reduced_size)
    end if
  end subroutine RemoveDuplicates_Int32

  pure subroutine RemoveDuplicates_int64(array,nduplicates)
    use, intrinsic :: iso_fortran_env
    implicit none
    integer(int8), parameter :: kind=int64
    integer(kind), allocatable, intent(inout) :: array(:)
    integer(int64),allocatable, intent(out), optional :: nduplicates(:)

    integer(kind), allocatable :: reduced_array(:)

    integer(kind), parameter :: empty=huge(1_kind)
    
    integer(int64) :: i, original_size, reduced_size, j
    integer(int64), allocatable :: countduplicates(:)
    
    original_size = size(array)
    allocate(reduced_array(original_size))
    reduced_array = empty

    if(present(nduplicates)) then
       allocate(countduplicates(original_size))
       countduplicates = 0
    end if
    
    reduced_size=0
    do i=1,original_size
       ! if the entry already exists, check next entry
       if(any(reduced_array == array(i))) then
          if(present(nduplicates)) then
             j = FindLoc(reduced_array,DIM=1,VALUE=array(i))
             countduplicates(j) = countduplicates(j) + 1_int64
          end if
          cycle
       else
          ! no match found, therefore add it to the output
          reduced_size = reduced_size + 1_int64

          reduced_array(reduced_size) = array(i)
          if(present(nduplicates)) then
             countduplicates(reduced_size) = 1_int64
          end if
       end if
    end do

    ! Now reallocate input-array and assign reduced content
    deallocate(array)
    allocate(array(reduced_size))
    array = reduced_array(1:reduced_size)

    if(present(nduplicates)) then
       allocate(nduplicates(reduced_size))
       nduplicates = countduplicates(1:reduced_size)
    end if
  end subroutine RemoveDuplicates_Int64

  pure subroutine RemoveDuplicates_Real32(array,nduplicates)
    use, intrinsic :: iso_fortran_env
    implicit none
    integer(int8), parameter :: kind=real32
    real(kind), allocatable, intent(inout) :: array(:)
    integer(int64),allocatable, intent(out), optional :: nduplicates(:)

    real(kind), allocatable :: reduced_array(:)

    real(kind), parameter :: empty=huge(1._kind)
    
    integer(int64) :: i, original_size, reduced_size, j
    integer(int64), allocatable :: countduplicates(:)
    
    original_size = size(array)
    allocate(reduced_array(original_size))
    reduced_array = empty

    if(present(nduplicates)) then
       allocate(countduplicates(original_size))
       countduplicates = 0
    end if
    
    reduced_size=0
    do i=1,original_size
       ! if the entry already exists, check next entry
       if(any(reduced_array == array(i))) then
          if(present(nduplicates)) then
             j = FindLoc(reduced_array,DIM=1,VALUE=array(i))
             countduplicates(j) = countduplicates(j) + 1_int64
          end if
          cycle
       else
          ! no match found, therefore add it to the output
          reduced_size = reduced_size + 1_int64

          reduced_array(reduced_size) = array(i)
          if(present(nduplicates)) then
             countduplicates(reduced_size) = 1_int64
          end if
       end if
    end do

    ! Now reallocate input-array and assign reduced content
    deallocate(array)
    allocate(array(reduced_size))
    array = reduced_array(1:reduced_size)

    if(present(nduplicates)) then
       allocate(nduplicates(reduced_size))
       nduplicates = countduplicates(1:reduced_size)
    end if
  end subroutine RemoveDuplicates_Real32
  
  pure subroutine RemoveDuplicates_Real64(array,nduplicates)
    use, intrinsic :: iso_fortran_env
    implicit none
    integer(int8), parameter :: kind=real64
    real(kind), allocatable, intent(inout) :: array(:)
    integer(int64),allocatable, intent(out), optional :: nduplicates(:)

    real(kind), allocatable :: reduced_array(:)

    real(kind), parameter :: empty=huge(1._kind)
    
    integer(int64) :: i, original_size, reduced_size, j
    integer(int64), allocatable :: countduplicates(:)
    
    original_size = size(array)
    allocate(reduced_array(original_size))
    reduced_array = empty

    if(present(nduplicates)) then
       allocate(countduplicates(original_size))
       countduplicates = 0
    end if
    
    reduced_size=0
    do i=1,original_size
       ! if the entry already exists, check next entry
       if(any(reduced_array == array(i))) then
          if(present(nduplicates)) then
             j = FindLoc(reduced_array,DIM=1,VALUE=array(i))
             countduplicates(j) = countduplicates(j) + 1_int64
          end if
          cycle
       else
          ! no match found, therefore add it to the output
          reduced_size = reduced_size + 1_int64

          reduced_array(reduced_size) = array(i)
          if(present(nduplicates)) then
             countduplicates(reduced_size) = 1_int64
          end if
       end if
    end do

    ! Now reallocate input-array and assign reduced content
    deallocate(array)
    allocate(array(reduced_size))
    array = reduced_array(1:reduced_size)

    if(present(nduplicates)) then
       allocate(nduplicates(reduced_size))
       nduplicates = countduplicates(1:reduced_size)
    end if
  end subroutine RemoveDuplicates_Real64
  
  pure subroutine RemoveDuplicates_Real128(array,nduplicates)
    use, intrinsic :: iso_fortran_env
    implicit none
    integer(int8), parameter :: kind=real128
    real(kind), allocatable, intent(inout) :: array(:)
    integer(int64),allocatable, intent(out), optional :: nduplicates(:)

    real(kind), allocatable :: reduced_array(:)

    real(kind), parameter :: empty=huge(1._kind)
    
    integer(int64) :: i, original_size, reduced_size, j
    integer(int64), allocatable :: countduplicates(:)
    
    original_size = size(array)
    allocate(reduced_array(original_size))
    reduced_array = empty

    if(present(nduplicates)) then
       allocate(countduplicates(original_size))
       countduplicates = 0
    end if
    
    reduced_size=0
    do i=1,original_size
       ! if the entry already exists, check next entry
       if(any(reduced_array == array(i))) then
          if(present(nduplicates)) then
             j = FindLoc(reduced_array,DIM=1,VALUE=array(i))
             countduplicates(j) = countduplicates(j) + 1_int64
          end if
          cycle
       else
          ! no match found, therefore add it to the output
          reduced_size = reduced_size + 1_int64

          reduced_array(reduced_size) = array(i)
          if(present(nduplicates)) then
             countduplicates(reduced_size) = 1_int64
          end if
       end if
    end do

    ! Now reallocate input-array and assign reduced content
    deallocate(array)
    allocate(array(reduced_size))
    array = reduced_array(1:reduced_size)

    if(present(nduplicates)) then
       allocate(nduplicates(reduced_size))
       nduplicates = countduplicates(1:reduced_size)
    end if
  end subroutine RemoveDuplicates_Real128


  pure subroutine RemoveDuplicates_Complex32(array,nduplicates)
    use, intrinsic :: iso_fortran_env
    implicit none
    integer(int8), parameter :: kind=real32
    complex(kind), allocatable, intent(inout) :: array(:)
    integer(int64),allocatable, intent(out), optional :: nduplicates(:)

    complex(kind), allocatable :: reduced_array(:)

    complex(kind), parameter :: empty=huge(1._kind)
    
    integer(int64) :: i, original_size, reduced_size, j
    integer(int64), allocatable :: countduplicates(:)
    
    original_size = size(array)
    allocate(reduced_array(original_size))
    reduced_array = empty

    if(present(nduplicates)) then
       allocate(countduplicates(original_size))
       countduplicates = 0
    end if
    
    reduced_size=0
    do i=1,original_size
       ! if the entry already exists, check next entry
       if(any(reduced_array == array(i))) then
          if(present(nduplicates)) then
             j = FindLoc(reduced_array,DIM=1,VALUE=array(i))
             countduplicates(j) = countduplicates(j) + 1_int64
          end if
          cycle
       else
          ! no match found, therefore add it to the output
          reduced_size = reduced_size + 1_int64

          reduced_array(reduced_size) = array(i)
          if(present(nduplicates)) then
             countduplicates(reduced_size) = 1_int64
          end if
       end if
    end do

    ! Now reallocate input-array and assign reduced content
    deallocate(array)
    allocate(array(reduced_size))
    array = reduced_array(1:reduced_size)

    if(present(nduplicates)) then
       allocate(nduplicates(reduced_size))
       nduplicates = countduplicates(1:reduced_size)
    end if
  end subroutine RemoveDuplicates_Complex32
  
  pure subroutine RemoveDuplicates_Complex64(array,nduplicates)
    use, intrinsic :: iso_fortran_env
    implicit none
    integer(int8), parameter :: kind=real64
    complex(kind), allocatable, intent(inout) :: array(:)
    integer(int64),allocatable, intent(out), optional :: nduplicates(:)

    complex(kind), allocatable :: reduced_array(:)

    complex(kind), parameter :: empty=huge(1._kind)
    
    integer(int64) :: i, original_size, reduced_size, j
    integer(int64), allocatable :: countduplicates(:)
    
    original_size = size(array)
    allocate(reduced_array(original_size))
    reduced_array = empty

    if(present(nduplicates)) then
       allocate(countduplicates(original_size))
       countduplicates = 0
    end if
    
    reduced_size=0
    do i=1,original_size
       ! if the entry already exists, check next entry
       if(any(reduced_array == array(i))) then
          if(present(nduplicates)) then
             j = FindLoc(reduced_array,DIM=1,VALUE=array(i))
             countduplicates(j) = countduplicates(j) + 1_int64
          end if
          cycle
       else
          ! no match found, therefore add it to the output
          reduced_size = reduced_size + 1_int64

          reduced_array(reduced_size) = array(i)
          if(present(nduplicates)) then
             countduplicates(reduced_size) = 1_int64
          end if
       end if
    end do

    ! Now reallocate input-array and assign reduced content
    deallocate(array)
    allocate(array(reduced_size))
    array = reduced_array(1:reduced_size)

    if(present(nduplicates)) then
       allocate(nduplicates(reduced_size))
       nduplicates = countduplicates(1:reduced_size)
    end if
  end subroutine RemoveDuplicates_Complex64
  
  pure subroutine RemoveDuplicates_Complex128(array,nduplicates)
    use, intrinsic :: iso_fortran_env
    implicit none
    integer(int8), parameter :: kind=real128
    complex(kind), allocatable, intent(inout) :: array(:)
    integer(int64),allocatable, intent(out), optional :: nduplicates(:)

    complex(kind), allocatable :: reduced_array(:)

    complex(kind), parameter :: empty=huge(1._kind)
    
    integer(int64) :: i, original_size, reduced_size, j
    integer(int64), allocatable :: countduplicates(:)
    
    original_size = size(array)
    allocate(reduced_array(original_size))
    reduced_array = empty

    if(present(nduplicates)) then
       allocate(countduplicates(original_size))
       countduplicates = 0
    end if
    
    reduced_size=0
    do i=1,original_size
       ! if the entry already exists, check next entry
       if(any(reduced_array == array(i))) then
          if(present(nduplicates)) then
             j = FindLoc(reduced_array,DIM=1,VALUE=array(i))
             countduplicates(j) = countduplicates(j) + 1_int64
          end if
          cycle
       else
          ! no match found, therefore add it to the output
          reduced_size = reduced_size + 1_int64

          reduced_array(reduced_size) = array(i)
          if(present(nduplicates)) then
             countduplicates(reduced_size) = 1_int64
          end if
       end if
    end do

    ! Now reallocate input-array and assign reduced content
    deallocate(array)
    allocate(array(reduced_size))
    array = reduced_array(1:reduced_size)

    if(present(nduplicates)) then
       allocate(nduplicates(reduced_size))
       nduplicates = countduplicates(1:reduced_size)
    end if
  end subroutine RemoveDuplicates_Complex128

  pure subroutine sort_int8(array)
    use, intrinsic :: iso_fortran_env
    implicit none
    integer(int8), parameter :: kind=int8
    integer(kind), intent(inout) :: array(:)

    integer(kind) :: a,b
    integer(int64) :: n, i,j

    n = size(array)

    do j=2_int64, N
       a=array(j)
       do i=j-1_int64,1_int64,-1_int64
          if (array(i)<=a) then
             exit
          else
             array(i+1_int64) = array(i)
          end if
       end do
       array(i+1_int64) =a
    end do
  end subroutine sort_int8

  pure subroutine sort_int16(array)
    use, intrinsic :: iso_fortran_env
    implicit none
    integer(int8), parameter :: kind=int16
    integer(kind), intent(inout) :: array(:)

    integer(kind) :: a,b
    integer(int64) :: n, i,j

    n = size(array)

    do j=2_int64, N
       a=array(j)
       do i=j-1_int64,1_int64,-1_int64
          if (array(i)<=a) then
             exit
          else
             array(i+1_int64) = array(i)
          end if
       end do
       array(i+1_int64) =a
    end do
  end subroutine sort_int16

  pure subroutine sort_int32(array)
    use, intrinsic :: iso_fortran_env
    implicit none
    integer(int8), parameter :: kind=int32
    integer(kind), intent(inout) :: array(:)

    integer(kind) :: a,b
    integer(int64) :: n, i,j

    n = size(array)

    do j=2_int64, N
       a=array(j)
       do i=j-1_int64,1_int64,-1_int64
          if (array(i)<=a) then
             exit
          else
             array(i+1_int64) = array(i)
          end if
       end do
       array(i+1_int64) =a
    end do
  end subroutine sort_int32

  pure subroutine sort_int64(array)
    use, intrinsic :: iso_fortran_env
    implicit none
    integer(int8), parameter :: kind=int64
    integer(kind), intent(inout) :: array(:)

    integer(kind) :: a,b
    integer(int64) :: n, i,j

    n = size(array)

    do j=2_int64, N
       a=array(j)
       do i=j-1_int64,1_int64,-1_int64
          if (array(i)<=a) then
             exit
          else
             array(i+1_int64) = array(i)
          end if
       end do
       array(i+1_int64) =a
    end do
  end subroutine sort_int64

  pure subroutine sort_real32(array)
    use, intrinsic :: iso_fortran_env
    implicit none
    integer(int8), parameter :: kind=real32
    real(kind),    intent(inout) :: array(:)

    real(kind) :: a,b
    integer(int64) :: n, i,j

    n = size(array)

    do j=2_int64, N
       a=array(j)
       do i=j-1_int64,1_int64,-1_int64
          if (array(i)<=a) then
             exit
          else
             array(i+1_int64) = array(i)
          end if
       end do
       array(i+1_int64) =a
    end do
  end subroutine sort_real32

  pure subroutine sort_real64(array)
    use, intrinsic :: iso_fortran_env
    implicit none
    integer(int8), parameter :: kind=real64
    real(kind),    intent(inout) :: array(:)

    real(kind) :: a,b
    integer(int64) :: n, i,j

    n = size(array)

    do j=2_int64, N
       a=array(j)
       do i=j-1_int64,1_int64,-1_int64
          if (array(i)<=a) then
             exit
          else
             array(i+1_int64) = array(i)
          end if
       end do
       array(i+1_int64) =a
    end do
  end subroutine sort_real64

  pure subroutine sort_real128(array)
    use, intrinsic :: iso_fortran_env
    implicit none
    integer(int8), parameter :: kind=real128
    real(kind),    intent(inout) :: array(:)

    real(kind) :: a,b
    integer(int64) :: n, i,j

    n = size(array)

    do j=2_int64, N
       a=array(j)
       do i=j-1_int64,1_int64,-1_int64
          if (array(i)<=a) then
             exit
          else
             array(i+1_int64) = array(i)
          end if
       end do
       array(i+1_int64) =a
    end do
  end subroutine sort_real128
end module arrayoperations

  
