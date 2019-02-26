!------------------------------------------------------------------------------
! IO, Input-Output-Managing Module
!------------------------------------------------------------------------------
!
! MODULE: io
!> @brief
!! Manages input and output using files
!! @author
!! Alexander Lehmann,
!! UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!! @date
!! 03.09.2018
!! @version
!! 1.0
! REVISION HISTORY:
! 03 09 2018 - Initial Version
module io
  USE,INTRINSIC :: ISO_FORTRAN_ENV ! defines kinds

  implicit none
  private
  public OpenFile,CloseFile
  
  !> Mininum for the file ID numbers
  integer(int8), private, parameter :: minID=10
  !> Maximum for the file ID numbers
  integer(int8), private, parameter :: maxID=99
  !> Number of total available file ID numbers
  integer(int8), private, parameter :: nID = maxID-minID+1
  !> Current occupation state of file ID numbers
  logical, private, dimension(minID:maxID) :: ID_used = .false.

contains
  !> @brief
  !! Opens a file
  !! @returns
  !! File ID number
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 03.09.2018
  !! @version 1.0
  impure function OpenFile(filename,st,fm,act,ios) result(FileID)
    implicit none
    !> File name
    character(len=*), intent(in)            :: filename
    !> Status (write, read, etc.)
    character(len=*), intent(in),  optional :: st
    !> Format
    character(len=*), intent(in),  optional :: fm
    !> Access
    character(len=*), intent(in),  optional :: act
    !> IOS status
    integer(int8),   intent(out), optional :: ios
    !> File ID number of the opened file
    integer(int8) :: FileID
    
    character(len=260) :: open_st,open_fm,open_acc,open_act
    integer(int8)      :: open_ios

    FileID = GetFileID()
    ! DEFAULT VALUES
    open_st     = 'REPLACE'
    open_fm     = 'FORMATTED'
    open_act    = 'WRITE'
    !open_pos    = 'ASIS'
    !open_pad    = 'NO'
    !open_del    = 'NONE'

    ! OPTIONAL INPUT VALUES
    if(present(st)) open_st = st
    if(present(fm)) open_fm = fm
    if(present(act)) open_act = act
    !if(present(pos)) open_pos = pos
    !if(present(pad)) open_pad = pad
    !if(present(del)) open_del = del
    
    Open(UNIT=FileID,FILE=filename,STATUS=st,FORM=fm,ACTION=act,IOSTAT=open_ios)
    if(present(ios)) ios = open_ios
  end function OpenFile
  
  !> @brief
  !! Closes a file
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 03.09.2018
  !! @version 1.0
  impure subroutine CloseFile(fileID)
    implicit none
    !> File ID number of the to be closed file
    integer(int8), intent(in) :: fileID

    close(fileID)
    call FreeFileID(fileID)
  end subroutine CloseFile

  !> @brief
  !! Takes an unused file-ID
  !! @returns
  !! File-ID
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 03.09.2018
  !! @version 1.0
  impure function GetFileID() result(fileID)
    use mpiinterface, only: ThisProc,MPIstop
    implicit none
    !> Free file-ID
    integer(int8) :: fileID
    
    integer(int8) :: i

    do i=minID,maxID
       if(.not.ID_used(i)) then
          ID_used(i) = .true.
          fileID = i
          return
       end if
    end do

    call MPIstop("GetFileID failed to return a free ID: All file-ID's occupied.")
  end function GetFileID

  !> @brief
  !! Frees an used file-ID
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 03.09.2018
  !! @version 1.0  
  impure subroutine FreeFileID(fileID)
    implicit none
    !> To-be freed file-ID
    integer(int8), intent(in) :: fileID

    ID_used(fileID) = .false.
  end subroutine FreeFileID
end module io
