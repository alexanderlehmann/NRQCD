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
       S2CS,C2CS,GetLinkCS_M,GetLinkCS_G,nDoF

  implicit none

  PRIVATE

  public GetGluonicWilsonLoop,GetWilsonLoop
       

contains ! Module procedures
  impure complex(fp) function GetGluonicWilsonLoop(&
       GaugeField_t1, GaugeField_t2, origin, distance, proddir)
    use matrixoperations, only: GetTrace
    implicit none
    !> SU(3)-Gauge configuration at \f$t=t_1\f$
    type(SU3GaugeConfiguration), intent(in) :: GaugeField_t1
    !> SU(3)-Gauge configuration at \f$t=t_2\f$
    type(SU3GaugeConfiguration), intent(in) :: GaugeField_t2
    !> Lattice index at origin of link product
    integer(int64),              intent(in) :: origin
    !> Distance in lattice units in direction (proddir)
    integer(int64),              intent(in) :: distance
    !> Direction into which to take the product
    integer(int8),               intent(in) :: proddir

    complex(fp), dimension(nDoF,nDoF) :: LinkProduct_t1, LinkProduct_t2, WholeProduct
    
    ! Getting link product
    LinkProduct_t1 = GetLinkProduct(GaugeField_t1,origin,distance,proddir)
    LinkProduct_t2 = GetLinkProduct(GaugeField_t2,origin,distance,proddir)

    ! Final multiplication step
    WholeProduct = matmul(LinkProduct_t1,conjg(transpose(LinkProduct_t2)))

    GetGluonicWilsonLoop = GetTrace(WholeProduct)/nDoF ! (sqrt(ndof))² due to quark-normalization
  end function GetGluonicWilsonLoop
  
  impure function GetLinkProduct(GaugeField,origin,distance,proddir)
    use mpiinterface, only: ThisProc, intmpi, GetComplexSendType
    use mpi
    use lattice, only: nDim, GetLatticeIndex, GetProc_G, GetNeib_G
    use matrixoperations, only: GetUnitMatrix
    implicit none
    !> SU(3)-Gauge configuration
    type(SU3GaugeConfiguration), intent(in) :: GaugeField
    !> Lattice index at origin of link product
    integer(int64),              intent(in) :: origin
    !> Distance in lattice units in direction (proddir)
    integer(int64),              intent(in) :: distance
    !> Direction into which to take the product
    integer(int8),               intent(in) :: proddir
    !> Link product
    complex(fp), dimension(nDoF,nDoF) :: GetLinkProduct

    ! Link at a site
    complex(fp), dimension(nDoF,nDoF) :: Link
    
    ! Indices
    integer(int64) :: k, xk

    ! MPI stuff
    integer(intmpi) :: proc, mpierr

    GetLinkProduct = GetUnitMatrix(nDoF)
    xk = origin
    do k=1,distance,+1
       proc = GetProc_G(xk)

       if(ThisProc()==proc) then
          Link = GetLinkCS_G(GaugeField,proddir,xk)
          GetLinkProduct = matmul(GetLinkProduct,Link)
       end if
       
       call mpi_bcast(GetLinkProduct,size(GetLinkProduct),GetComplexSendType(),&
            proc,mpi_comm_world,mpierr)
       
       xk = GetNeib_G(proddir,xk)
    end do
  end function GetLinkProduct

  impure function GetWilsonLine(GaugeField,origin,distance,proddir)
    use mpiinterface, only: ThisProc, intmpi, GetComplexSendType
    use mpi
    use lattice, only: nDim, GetLatticeIndex, GetProc_G, GetNeib_G
    use matrixoperations, only: GetUnitMatrix
    implicit none
    !> SU(3)-Gauge configuration
    type(SU3GaugeConfiguration), intent(in) :: GaugeField
    !> Lattice index at origin of link product
    integer(int64),              intent(in) :: origin
    !> Distance in lattice units in direction (proddir)
    integer(int64),              intent(in) :: distance
    !> Direction into which to take the product
    integer(int8),               intent(in) :: proddir
    !> Link product
    complex(fp), dimension(nColours,nColours) :: GetWilsonLine

    ! Link at a site
    complex(fp), dimension(nColours,nColours) :: Link
    
    ! Indices
    integer(int64) :: k, xk

    ! MPI stuff
    integer(intmpi) :: proc, mpierr

    GetWilsonLine = GetUnitMatrix(nColours)
    xk = origin
    do k=1,distance,+1
       proc = GetProc_G(xk)

       if(ThisProc()==proc) then
          Link = GaugeField%GetLink_G(proddir,xk)
          GetWilsonLine = matmul(GetWilsonLine,Link)
       end if
       
       call mpi_bcast(GetWilsonLine,size(GetWilsonLine),GetComplexSendType(),&
            proc,mpi_comm_world,mpierr)
       
       xk = GetNeib_G(proddir,xk)
    end do
  end function GetWilsonLine

  impure function GetWilsonLoop(&
       GaugeField_t1, GaugeField_t2,distance, proddir, origin_in) result(res)
    use matrixoperations, only: GetTrace
    use lattice, only: GetLatticeSize
    use mpiinterface, only: ThisProc
    implicit none
    !> SU(3)-Gauge configuration at \f$t=t_1\f$
    type(SU3GaugeConfiguration), intent(in) :: GaugeField_t1
    !> SU(3)-Gauge configuration at \f$t=t_2\f$
    type(SU3GaugeConfiguration), intent(in) :: GaugeField_t2
    !> Distance in lattice units in direction (proddir)
    integer(int64),              intent(in) :: distance
    !> Direction into which to take the product
    integer(int8),               intent(in) :: proddir
    complex(fp) :: res
    
    complex(fp), dimension(nColours,nColours) :: LinkProduct_t1, LinkProduct_t2, WholeProduct
    
    integer(int64), optional :: origin_in
    integer(int64) :: origin

    if(present(origin_in)) then
       origin = origin_in
       
       ! Getting link product
       LinkProduct_t1 = GetWilsonLine(GaugeField_t1,origin,distance,proddir)
       LinkProduct_t2 = GetWilsonLine(GaugeField_t2,origin,distance,proddir)

       ! Final multiplication step
       WholeProduct = matmul(LinkProduct_t1,conjg(transpose(LinkProduct_t2)))

       res = GetTrace(WholeProduct)/nColours
    else
       res = 0

       do origin=1,GetLatticeSize()

          ! Getting link product
          LinkProduct_t1 = GetWilsonLine(GaugeField_t1,origin,distance,proddir)
          LinkProduct_t2 = GetWilsonLine(GaugeField_t2,origin,distance,proddir)

          ! Final multiplication step
          WholeProduct = matmul(LinkProduct_t1,conjg(transpose(LinkProduct_t2)))

          res = res + GetTrace(WholeProduct)/nColours
       end do
       res = res / GetLatticeSize()
    end if
  end function GetWilsonLoop
end module WilsonLine
