program main
  use, intrinsic :: iso_fortran_env
  IMPLICIT NONE

  character(len=100) :: FileName_input, FileName_output
  integer(int8) :: FileID_input, FileID_output

  integer :: ix,iy,iz
  integer(int8), parameter :: nDIM=3_int8
  real :: EigenVector(nDIM,nDIM), EigenValue(nDIM)
  
  ! Transformation into cylindrical coordinates
  !real :: rho, phi, z

  integer :: i
  real, parameter :: pi = acos(-1.)

  call get_command_argument(1,FileName_input)
  print*,'INPUT-FILE:'
  print*,FileName_input
  
  call get_command_argument(2,FileName_output)
  print*,'OUTPUT-FILE:'
  print*,FileName_output
  
  FileID_input  = 10
  FileID_output = 11
  OPEN(UNIT=FileID_input, FILE=FileName_input, STATUS='OLD',    FORM='FORMATTED',ACTION='READ')
  OPEN(UNIT=FileID_output,FILE=FileName_output,STATUS='REPLACE',FORM='FORMATTED',ACTION='WRITE')
  
  ! Getting rid of header first
  READ(FileID_input,*,END=1)
  do ! Process until end of file
     READ(FileID_input,*,END=1) ix,iy,iz,&
          EigenValue(1),EigenVector(1:nDIM,1),EigenValue(2),EigenVector(1:nDIM,2),EigenValue(3),EigenVector(1:nDIM,nDIM)

     WRITE(FileID_output,'(3(I2.2,1X))',advance='no') ix, iy, iz
     ! Computing cylindrical coordinates, projecting and then computing cartesian ones
     do i=1,nDIM
        !rho = norm2(EigenVector(1:2,i))
        !phi = atan2(EigenVector(2,i),EigenVector(1,i))
        !z   = EigenVector(nDIM,i)

        !if(abs(phi)>pi/2) then
        !   EigenVector(1:nDIM,i) = -EigenVector(1:nDIM,i)
        !else
        !end if
        if(EigenVector(3,i)<0) EigenVector(1:nDIM,i) = -EigenVector(1:nDIM,i)

        write(FileID_output,'(3(SP,E13.6,1X))',advance='no') &
             EigenValue(i)*EigenVector(1:nDIM,i)
             ![rho*cos(phi_new),rho*sin(phi_new),z_new]
     end do
     write(FileID_output,*)
  end do
1 continue
  CLOSE(FileID_input)
  CLOSE(FileID_output)

end program main
