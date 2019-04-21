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

  implicit none

  PRIVATE

  public GetFermionicWilsonLoop, GetGluonicWilsonLoop, GetPointSplitCorrelator

contains ! Module procedures
  impure complex(fp) function GetPointSplitCorrelator(&
       GaugeField,Heavyfield,origin,distance,proddir)
    use mpiinterface, only: ThisProc, intmpi, GetComplexSendType
    use mpi
    use lattice, only: nDim, GetLatticeIndex, GetProc_G, GetNeib_G
    use matrixoperations, only: GetTrace
    use su2, only: SU2Generators => Generators
    implicit none
    !> SU(3)-Gauge configuration
    type(SU3GaugeConfiguration), intent(in) :: GaugeField
    !> NRQCD heavy field with quark- and antiquark-propagator initialised at origin
    type(NRQCDField),            intent(in) :: HeavyField
    !> Lattice index at origin of link product
    integer(int64),              intent(in) :: origin
    !> Distance in lattice units in direction (proddir)
    integer(int64),              intent(in) :: distance
    !> Direction into which to take the product
    integer(int8),               intent(in) :: proddir

    integer(int64) :: quark_index, antiq_index

    ! Indices
    integer(int64) :: LatticeIndex_distance, shiftstep
    integer(int8) :: k
    
    ! MPI stuff
    integer(intmpi) :: proc, mpierr
    
    complex(fp), dimension(nDoF,nDoF) :: LinkProduct, WholeProduct,&
         quarkprop, antiqprop,PauliMatrix
    
    
    
    quark_index = origin
    antiq_index = origin
    do shiftstep=1,distance/2
       quark_index = GetNeib_G(+proddir,quark_index)
       antiq_index = GetNeib_G(-proddir,antiq_index)
    end do
    
    
    ! Getting link product
    LinkProduct = GetLinkProduct(GaugeField,antiq_index,distance,proddir)
    LinkProduct = transpose(LinkProduct)
    
    ! Getting propagators
    proc = GetProc_G(quark_index)
    if(ThisProc()==proc) then
       QuarkProp = HeavyField%GetQuarkProp_G(quark_index)
    end if
    call mpi_bcast(QuarkProp,size(QuarkProp),GetComplexSendType(),&
         proc,mpi_comm_world,mpierr)

    proc = GetProc_G(antiq_index)
    if(ThisProc()==proc) then
       AntiQProp = HeavyField%GetAntiQProp_G(antiq_index)
    end if
    call mpi_bcast(AntiQProp,size(AntiQProp),GetComplexSendType(),&
         proc,mpi_comm_world,mpierr)

    GetPointSplitCorrelator = 0
    
    ! Final multiplication step
    do k=1,3
       PauliMatrix = S2CS(SU2Generators(:,:,k))

       WholeProduct = matmul(matmul(matmul(matmul(&
            quarkprop,&
            paulimatrix),&
            conjg(antiqprop)),&
            conjg(paulimatrix)),&
            conjg(LinkProduct))

       GetPointSplitCorrelator = GetPointSplitCorrelator + GetTrace(WholeProduct)
    end do

  end function GetPointSplitCorrelator
  

  impure complex(fp) function GetFermionicWilsonLoop(&
       GaugeField_t1, GaugeField_t2, HeavyField, origin, distance, proddir)
    use mpiinterface, only: ThisProc, intmpi, GetComplexSendType
    use mpi
    use lattice, only: nDim, GetLatticeIndex, GetProc_G, GetNeib_G
    use matrixoperations, only: GetTrace
    implicit none
    !> SU(3)-Gauge configuration at \f$t=t_1\f$
    type(SU3GaugeConfiguration), intent(in) :: GaugeField_t1
    !> SU(3)-Gauge configuration at \f$t=t_2\f$
    type(SU3GaugeConfiguration), intent(in) :: GaugeField_t2
    !> NRQCD heavy field with quarkpropagator initialised at origin and antiquark at origin+r
    type(NRQCDField),            intent(in) :: HeavyField
    !> Lattice index at origin of link product
    integer(int64),              intent(in) :: origin
    !> Distance in lattice units in direction (proddir)
    integer(int64),              intent(in) :: distance
    !> Direction into which to take the product
    integer(int8),               intent(in) :: proddir

    complex(fp), dimension(nDoF,nDoF) :: LinkProduct_t1, LinkProduct_t2, WholeProduct,&
         quarkprop, antiqprop
    
    ! Indices
    integer(int64) :: LatticeIndex_distance, shiftstep
    
    ! MPI stuff
    integer(intmpi) :: proc, mpierr

    ! Getting link product
    LinkProduct_t1 = GetLinkProduct(GaugeField_t1,origin,distance,proddir)
    LinkProduct_t2 = GetLinkProduct(GaugeField_t2,origin,distance,proddir)

    ! Getting propagators
    proc = GetProc_G(origin)
    if(ThisProc()==proc) then
       AntiQProp = HeavyField%GetAntiQProp_G(origin)
    end if
    call mpi_bcast(AntiQProp,size(AntiQProp),GetComplexSendType(),&
         proc,mpi_comm_world,mpierr)
    
    LatticeIndex_distance = origin
    do shiftstep=1,distance
       LatticeIndex_distance=GetNeib_G(proddir,LatticeIndex_distance)
    end do
    proc = GetProc_G(LatticeIndex_distance)
    if(ThisProc()==proc) then
       QuarkProp = HeavyField%GetQuarkProp_G(LatticeIndex_distance)
    end if
    call mpi_bcast(QuarkProp,size(QuarkProp),GetComplexSendType(),&
         proc,mpi_comm_world,mpierr)

    ! Final multiplication step
    WholeProduct = matmul(matmul(matmul(&
         LinkProduct_t1,&
         QuarkProp),&
         conjg(transpose(LinkProduct_t2))),&
         ! Transforming from t1->t2 with c* into t2->t1 propagator with c [c=wilsoncoefficients]
         conjg(transpose(AntiQProp)))

    GetFermionicWilsonLoop = GetTrace(WholeProduct)
  end function GetFermionicWilsonLoop

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

end module WilsonLine
