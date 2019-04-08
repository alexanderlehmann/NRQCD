!-------------------------------------------------------------------------
! Wilsonline
!-------------------------------------------------------------------------
!
! MODULE: wilsonline
!>@brief Contains several definition of the wilson line
!!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!!@date 02.04.2019
!!@version 1.0
! REVISION HISTORY:
! 02 04 2019 - Initial Version
!-------------------------------------------------------------------------
module WilsonLine
  use, intrinsic :: iso_fortran_env
  use precision, only: fp

  use gaugeconfiguration_su3, only: &
       SU3GaugeConfiguration => GaugeConfiguration
  use su3, only: &
       SU3Generators => Generators, &
       nColours      => nSUN, &
       nGluons       => nGen
  use su2, only: &
       SU2Generators => Generators, &
       nSpins        => nSUN
  use nrqcd, only:&
       S2CS,C2CS,GetLinkCS_M,GetLinkCS_G,&
       NRQCDField,&
       nDoF

  use mpiinterface, only: mpistop
  
  implicit none
  
  PRIVATE

  public GetFermionicWilsonLoop
  
contains ! Module procedures

  impure complex(fp) function GetFermionicWilsonLoop(GaugeField_t1,GaugeField_t2,HeavyField_x0,HeavyField_xr,x0,r,messdir)
    use lattice, only: nDim, GetLatticeExtension, GetLatticeSize, GetLatticeIndex,GetProc_G,GetNeib_g
    use mpi
    use mpiinterface, only: intmpi, GetComplexSendType, ThisProc
    use matrixoperations, only: GetTrace
    implicit none
    !> SU(3)-Gauge configuration at \f$t=t_1\f$
    type(SU3GaugeConfiguration), intent(in) :: GaugeField_t1
    !> SU(3)-Gauge configuration at \f$t=t_2\f$
    type(SU3GaugeConfiguration), intent(in) :: GaugeField_t2
    !> NRQCD heavy field with quarkpropagator initialised at x0
    type(NRQCDField),            intent(in) :: HeavyField_x0
    !> NRQCD heavy field with antiquarkpropagator initialised at xr
    type(NRQCDField),            intent(in) :: HeavyField_xr
    !> Lattice index x0
    integer(int64), intent(in) :: x0
    !> Distance in lattice units in i'th direction
    integer(int64), intent(in) :: r
    !> Direction
    integer(int8), intent(in) :: messdir

    complex(fp), dimension(ndof,ndof) :: WilsonLoop, WilsonLine_t1, WilsonLine_t2
    complex(fp), dimension(ndof,ndof) :: quark_x0, antiq_xr

    integer(int64) :: xr, k

    ! MPI stuff
    integer(intmpi) :: proc_x0, proc_xr, mpierr
    
    WilsonLine_t1 = GetWilsonLine(GaugeField_t1,x0,r,messdir)
    WilsonLine_t2 = GetWilsonLine(GaugeField_t2,x0,r,messdir)

    proc_x0 = GetProc_G(x0)
    if(ThisProc()==proc_x0) then
       quark_x0 = HeavyField_x0%GetQuarkProp_G(x0)
    end if
    call mpi_bcast(quark_x0,size(quark_x0),GetComplexSendType(),proc_x0,mpi_comm_world,mpierr)
    
    xr = x0
    do k=1,r
       xr = GetNeib_G(messdir,xr)
    end do
    proc_xr = GetProc_G(xr)
    if(ThisProc()==proc_xr) then
       antiq_xr = HeavyField_xr%GetAntiQProp_G(xr)
    end if
    call mpi_bcast(antiq_xr,size(antiq_xr),GetComplexSendType(),proc_xr,mpi_comm_world,mpierr)

    WilsonLoop = matmul(matmul(matmul(&
         conjg(transpose(WilsonLine_t1(:,:))),&
         quark_x0),&
         WilsonLine_t2(:,:)),&
         conjg(antiq_xr))

    GetFermionicWilsonLoop = GetTrace(WilsonLoop)
  end function GetFermionicWilsonLoop

  impure function GetWilsonLine(GaugeField,x0,r,messdir)
    use mpiinterface, only: ThisProc, intmpi, GetComplexSendType
    use mpi
    use lattice, only: nDim, GetLatticeIndex, GetProc_G, GetNeib_G
    use matrixoperations, only: GetUnitMatrix
    implicit none
    !> SU(3)-Gauge configuration
    type(SU3GaugeConfiguration), intent(in) :: GaugeField
    !> Lattice index
    integer(int64),              intent(in) :: x0
    !> Distance in lattice units in i'th direction
    integer(int64),              intent(in) :: r
    !> Direction
    integer(int8),               intent(in) :: messdir
    !> Wilson line
    complex(fp), dimension(nDoF,nDoF) :: GetWilsonLine
    
    complex(fp), dimension(nDoF,nDoF) :: Link

    ! Indices
    integer(int64) :: k, xk, latticeindex

    ! MPI stuff
    integer(intmpi) :: proc, mpierr

    GetWilsonLine = GetUnitMatrix(nDoF)
    xk = x0
    do k=0,r-1,+1
       proc = GetProc_G(xk)

       if(ThisProc()==proc) then
          Link = GetLinkCS_G(GaugeField,messdir,xk)
          GetWilsonLine = matmul(GetWilsonLine,Link)
       end if
       call mpi_bcast(GetWilsonLine,size(GetWilsonLine),GetComplexSendType(),&
            proc,mpi_comm_world,mpierr)
       
       xk = GetNeib_G(messdir,xk)
    end do
  end function GetWilsonLine


  
  !>@brief Computes the gluonic on-axis-wilson-loop
  !!@returns Gluonic on-axis-wilson loop
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 02.04.2019
  !!@version 1.0
  impure function GetGluonicWilsonLoops(GaugeField_t1, GaugeField_t2, rmax) result(WilsonLoops)
    use lattice, only: nDim, GetLatticeExtension, GetLatticeSize
    use mpi
    use mpiinterface, only: intmpi, GetComplexSendType
    use matrixoperations, only: GetTrace
    implicit none
    !> SU(3)-Gauge configuration at \f$t=t_1\f$
    type(SU3GaugeConfiguration), intent(in) :: GaugeField_t1
    !> SU(3)-Gauge configuration at \f$t=t_2\f$
    type(SU3GaugeConfiguration), intent(in) :: GaugeField_t2
    !> Maximal distance in lattice units
    integer(int64),              intent(in) :: rmax

    complex(fp), dimension(0:rmax) :: WilsonLoops
    
    integer(int64) :: x(nDim), x1,x2,x3, n1,n2,n3, LatticeIndex, r
    integer(int8) :: i
    complex(fp), dimension(ndof,ndof,0:rmax) :: Lines_t1,Lines_t2
    complex(fp), dimension(ndof,ndof) :: WilsonLoop
    
    WilsonLoops = 0
    
    n1 = GetLatticeExtension(1)
    n2 = GetLatticeExtension(2)
    n3 = GetLatticeExtension(3)

    do i=1,nDim
       do x3=1,n3
          do x2=1,n2
             do x1=1,n1
                x  = [x1,x2,x3]

                Lines_t1 = GetWilsonLines(GaugeField_t1,x,rmax,i)
                Lines_t2 = GetWilsonLines(GaugeField_t2,x,rmax,i)
                
                do concurrent(r=0:rmax)
                   WilsonLoop = matmul(conjg(transpose(lines_t1(:,:,r))), lines_t2(:,:,r))

                   WilsonLoops(r) = WilsonLoops(r)&
                        + GetTrace(WilsonLoop)/GetLatticeSize()/nDoF/nDim
                end do
             end do
          end do
       end do
    end do

  end function GetGluonicWilsonLoops

  !>@brief Computes the fermionic on-axis-wilson-loop
  !!@returns  Fermionic on-axis-wilson loop
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 03.04.2019
  !!@version 1.0
  impure function GetFermionicWilsonLoops_Approximated(GaugeField_t1,GaugeField_t2,HeavyField,rmax)&
       result(WilsonLoops)
    use lattice, only: nDim, GetLatticeExtension, GetLatticeSize, GetLatticeIndex, GetProc_G
    use mpi
    use mpiinterface, only: intmpi, GetComplexSendType, ThisProc
    use matrixoperations, only: GetTrace
    implicit none
    !> SU(3)-Gauge configuration at \f$t=t_1\f$
    type(SU3GaugeConfiguration), intent(in) :: GaugeField_t1
    !> SU(3)-Gauge configuration at \f$t=t_2\f$
    type(SU3GaugeConfiguration), intent(in) :: GaugeField_t2
    !> NRQCD heavy field
    type(NRQCDField),            intent(in) :: HeavyField
    !> Maximal distance in lattice units
    integer(int64),              intent(in) :: rmax
    
    complex(fp), dimension(0:rmax) :: WilsonLoops
    
    integer(int64) :: x(nDim), x1,x2,x3, n1,n2,n3, LatticeIndex, r
    integer(int8) :: i
    complex(fp), dimension(ndof,ndof,0:rmax) :: Lines_t1,Lines_t2
    complex(fp), dimension(ndof,ndof) :: quark_x, antiq_x
    complex(fp), dimension(ndof,ndof) :: WilsonLoop
    
    integer(intmpi) :: proc, src, tag, mpierr, recvstatus(MPI_STATUS_SIZE)

    WilsonLoops = 0
    
    n1 = GetLatticeExtension(1)
    n2 = GetLatticeExtension(2)
    n3 = GetLatticeExtension(3)

    do i=1,nDim
       do x3=1,n3
          do x2=1,n2
             do x1=1,n1
                x = [x1,x2,x3]

                ! Getting the wilson lines at t1 and t2 from 0 to x
                Lines_t1 = GetWilsonLines(GaugeField_t1,x,rmax,i)
                Lines_t2 = GetWilsonLines(GaugeField_t2,x,rmax,i)

                LatticeIndex = GetLatticeIndex(x)
                proc = GetProc_G(LatticeIndex)
                if(ThisProc()==proc) then
                   quark_x = HeavyField%GetQuarkProp_G(LatticeIndex)
                   antiq_x = HeavyField%GetAntiQProp_G(LatticeIndex)
                   do r=0,rmax
                      WilsonLoop = matmul(matmul(matmul(&
                           conjg(transpose(lines_t1(:,:,r))), &
                           quark_x), &
                           lines_t2(:,:,r)),&
                           conjg(antiq_x))

                      WilsonLoops(r) = WilsonLoops(r)&
                           + GetTrace(WilsonLoop)/nDim
                   end do
                end if
                call mpi_bcast(WilsonLoops,size(wilsonloops),GetComplexSendType(),&
                     proc,mpi_comm_world,mpierr)
             end do
          end do
       end do
    end do

  end function GetFermionicWilsonLoops_Approximated

  !>@brief Computes the fermionic on-axis-wilson-loop
  !!@returns Fermionic on-axis-wilson loop
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 02.04.2019
  !!@version 1.0
  impure function GetFermionicWilsonLoops(GaugeField_t1,GaugeField_t2,HeavyField,rmax) result(WilsonLoops)
    use lattice, only: nDim, GetLatticeExtension, GetLatticeSize, GetLatticeIndex, GetProc_G
    use mpi
    use mpiinterface, only: intmpi, GetComplexSendType, ThisProc
    use matrixoperations, only: GetTrace
    implicit none
    !> SU(3)-Gauge configuration at \f$t=t_1\f$
    type(SU3GaugeConfiguration), intent(in) :: GaugeField_t1
    !> SU(3)-Gauge configuration at \f$t=t_2\f$
    type(SU3GaugeConfiguration), intent(in) :: GaugeField_t2
    !> NRQCD heavy field
    type(NRQCDField),            intent(in) :: HeavyField
    !> Maximal distance in lattice units
    integer(int64),              intent(in) :: rmax
    
    complex(fp), dimension(0:rmax) :: WilsonLoops
    
    integer(int64) :: x(nDim), x1,x2,x3, n1,n2,n3, LatticeIndex, r, xshifted(ndim)
    integer(int8) :: i
    complex(fp), dimension(ndof,ndof,0:rmax) :: Lines_t1,Lines_t2
    complex(fp), dimension(ndof,ndof) :: quark_x, quark_r, antiq_x, antiq_r
    complex(fp), dimension(ndof,ndof) :: WilsonLoop

    integer(intmpi), parameter :: buffersize=nDoF**2
    integer(intmpi), parameter :: root=0
    integer(intmpi) :: proc, src, tag, mpierr, recvstatus(MPI_STATUS_SIZE)
    
    WilsonLoops = 0
    
    n1 = GetLatticeExtension(1)
    n2 = GetLatticeExtension(2)
    n3 = GetLatticeExtension(3)

    do i=1,nDim
       do x3=1,n3
          do x2=1,n2
             do x1=1,n1
                x  = [x1,x2,x3]
                xshifted = x

                Lines_t1 = GetWilsonLines(GaugeField_t1,x,rmax,i)
                Lines_t2 = GetWilsonLines(GaugeField_t2,x,rmax,i)

                do r=0,rmax
                   latticeindex = GetLatticeIndex(x)
                   proc = GetProc_G(LatticeIndex)
                   if(ThisProc()==proc) then
                      quark_x = HeavyField%GetQuarkProp_G(LatticeIndex)
                      antiq_x = HeavyField%GetAntiQProp_G(LatticeIndex)
                   end if

                   if(proc.ne.root) then
                      tag = LatticeIndex
                      src = proc

                      if(ThisProc()==proc) then
                         call mpi_send(&
                              quark_x(              & ! What to send ...
                              1,1),                 & ! ... and it's first index
                              buffersize,           & ! How many points
                              GetComplexSendType(), & ! What type to send
                              root,                 & ! Destination
                              tag,                  & ! Tag
                              MPI_COMM_WORLD,       & ! Communicator
                              mpierr)                 ! Error code

                         tag = LatticeIndex + GetLatticeSize()
                         
                         call mpi_send(&
                              antiq_x(              & ! What to send ...
                              1,1),                 & ! ... and it's first index
                              buffersize,           & ! How many points
                              GetComplexSendType(), & ! What type to send
                              root,                 & ! Destination
                              tag,                  & ! Tag
                              MPI_COMM_WORLD,       & ! Communicator
                              mpierr)                 ! Error code
                      elseif(ThisProc()==root) then
                         call mpi_recv(&
                              quark_x(              & ! What to recieve ...
                              1,1),                 & ! ... and it's first index
                              buffersize,           & ! How many points
                              GetComplexSendType(), & ! What type to recive
                              src,                  & ! Source
                              tag,                  & ! Tag
                              MPI_COMM_WORLD,       & ! Communicator
                              recvstatus,           & ! Recieve-status
                              mpierr)                 ! Error code

                         tag = LatticeIndex + GetLatticeSize()
                         
                         call mpi_recv(&
                              antiq_x(              & ! What to recieve ...
                              1,1),                 & ! ... and it's first index
                              buffersize,           & ! How many points
                              GetComplexSendType(), & ! What type to recive
                              src,                  & ! Source
                              tag,                  & ! Tag
                              MPI_COMM_WORLD,       & ! Communicator
                              recvstatus,           & ! Recieve-status
                              mpierr)                 ! Error code
                      end if
                   end if
                   
                   xshifted(i) = x(i) + r
                   LatticeIndex = GetLatticeIndex(xshifted)
                   proc = GetProc_G(LatticeIndex)
                   if(ThisProc()==proc) then
                      quark_r = HeavyField%GetQuarkProp_G(LatticeIndex)
                      antiq_r = HeavyField%GetAntiQProp_G(LatticeIndex)
                   end if

                   if(proc.ne.root) then
                      tag = LatticeIndex + 2*GetLatticeSize()
                      src = proc

                      if(ThisProc()==proc) then
                         call mpi_send(&
                              quark_r(              & ! What to send ...
                              1,1),                 & ! ... and it's first index
                              buffersize,           & ! How many points
                              GetComplexSendType(), & ! What type to send
                              root,                 & ! Destination
                              tag,                  & ! Tag
                              MPI_COMM_WORLD,       & ! Communicator
                              mpierr)                 ! Error code

                         tag = LatticeIndex + 3*GetLatticeSize()
                         
                         call mpi_send(&
                              antiq_r(              & ! What to send ...
                              1,1),                 & ! ... and it's first index
                              buffersize,           & ! How many points
                              GetComplexSendType(), & ! What type to send
                              root,                 & ! Destination
                              tag,                  & ! Tag
                              MPI_COMM_WORLD,       & ! Communicator
                              mpierr)                 ! Error code
                      elseif(ThisProc()==root) then
                         call mpi_recv(&
                              quark_r(              & ! What to recieve ...
                              1,1),                 & ! ... and it's first index
                              buffersize,           & ! How many points
                              GetComplexSendType(), & ! What type to recive
                              src,                  & ! Source
                              tag,                  & ! Tag
                              MPI_COMM_WORLD,       & ! Communicator
                              recvstatus,           & ! Recieve-status
                              mpierr)                 ! Error code

                         tag = LatticeIndex + 3*GetLatticeSize()
                         
                         call mpi_recv(&
                              antiq_r(              & ! What to recieve ...
                              1,1),                 & ! ... and it's first index
                              buffersize,           & ! How many points
                              GetComplexSendType(), & ! What type to recive
                              src,                  & ! Source
                              tag,                  & ! Tag
                              MPI_COMM_WORLD,       & ! Communicator
                              recvstatus,           & ! Recieve-status
                              mpierr)                 ! Error code
                      end if
                   end if

                   if(ThisProc()==root) then
                      !WilsonLoop = &
                      !     ( &
                      !     + matmul(matmul(matmul(&
                      !     lines_t1(:,:,r), &
                      !     conjg(transpose(quark_x))),&
                      !     conjg(transpose(lines_t2(:,:,r)))),&
                      !     antiq_r)&
                      !     + matmul(matmul(matmul(&
                      !     conjg(transpose(lines_t1(:,:,r))), &
                      !     quark_r),&
                      !     lines_t2(:,:,r)),&
                      !     conjg(transpose(antiq_x)))&
                      !     ) / 2

                      WilsonLoop = matmul(matmul(matmul(&
                           lines_t1(:,:,r), &
                           conjg(transpose(quark_x))),&
                           conjg(transpose(lines_t2(:,:,r)))),&
                           antiq_r) !&
                           !+ matmul(matmul(&
                           !lines_t1(:,:,r), &
                           !conjg(transpose(quark_x))),&
                           !conjg(transpose(lines_t2(:,:,r))))
                      
                      WilsonLoops(r) = WilsonLoops(r)&
                           + GetTrace(WilsonLoop)/nDoF/nDim!/GetLatticeSize()/nDoF/nDim!*cmplx(0,-1,fp)
                   end if
                end do
             end do
          end do
       end do
    end do
    
    call mpi_bcast(WilsonLoops,size(wilsonloops),GetComplexSendType(),root,mpi_comm_world,mpierr)
  end function GetFermionicWilsonLoops

  !>@brief Computes Wilson line from x to x + r * î
  !!@returns Gluonic wilson loop
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 02.04.2019
  !!@version 1.0
  impure function GetWilsonLines(GaugeField,x,rmax,i) result(WilsonLines)
    use mpiinterface, only: ThisProc, intmpi, GetComplexSendType
    use mpi
    use lattice, only: nDim, GetLatticeIndex, GetProc_G
    use matrixoperations, only: GetUnitMatrix
    implicit none
    !> SU(3)-Gauge configuration
    type(SU3GaugeConfiguration), intent(in) :: GaugeField
    !> Lattice position
    integer(int64),              intent(in) :: x(nDim)
    !> Distance in lattice units in i'th direction
    integer(int64),              intent(in) :: rmax
    !> Direction
    integer(int8),               intent(in) :: i
    
    complex(fp), dimension(nDoF,nDoF,0:rmax) :: WilsonLines

    integer(int64) :: r, xshifted(ndim), latticeindex
    
    complex(fp),dimension(ndof,ndof) :: Link
   
    integer(intmpi), parameter :: buffersize=nDoF**2
    integer(intmpi), parameter :: root=0
    integer(intmpi) :: proc, src, tag, mpierr, recvstatus(MPI_STATUS_SIZE)
    
    xshifted    = x
    do r=0,rmax
       if(r==0) then
          WilsonLines(:,:,0) = GetUnitMatrix(nDoF)
       else
          WilsonLines(:,:,r) = WilsonLines(:,:,r-1)
          
          xshifted(i) = x(i) + r

          LatticeIndex = GetLatticeIndex(xshifted)
          proc = GetProc_G(LatticeIndex)

          if(ThisProc()==proc) then
             Link = GetLinkCS_G(GaugeField,i,LatticeIndex)
          end if

          if(proc.ne.root) then
             tag = LatticeIndex
             src = proc

             if(ThisProc()==proc) then
                call mpi_send(&
                     Link(                 & ! What to send ...
                     1,1),                 & ! ... and it's first index
                     buffersize,           & ! How many points
                     GetComplexSendType(), & ! What type to send
                     root,                 & ! Destination
                     tag,                  & ! Tag
                     MPI_COMM_WORLD,       & ! Communicator
                     mpierr)                 ! Error code
             elseif(ThisProc()==root) then
                call mpi_recv(&
                     Link(                 & ! What to recieve ...
                     1,1),                 & ! ... and it's first index
                     buffersize,           & ! How many points
                     GetComplexSendType(), & ! What type to recive
                     src,                  & ! Source
                     tag,                  & ! Tag
                     MPI_COMM_WORLD,       & ! Communicator
                     recvstatus,           & ! Recieve-status
                     mpierr)                 ! Error code
             end if
          end if

          if(proc==root) then
             WilsonLines(:,:,r) = matmul(WilsonLines(:,:,r),Link)
          end if

          call mpi_bcast(WilsonLines(1,1,r),size(WilsonLines(:,:,r)),GetComplexSendType(),root,mpi_comm_world,mpierr)
       end if

    end do
  end function GetWilsonLines
end module WilsonLine
