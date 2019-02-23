!----------------------------------------------------------------------
! RTLQCD, Real-Time-Lattice-QCD Simulation of Gauge Fields
!----------------------------------------------------------------------
!
!>@brief SU(3)-gauge-link-configuration in temporal gauge
!!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!!@date 22.01.2019
!!@version 1.0
!----------------------------------------------------------------------
module gaugeconfiguration_su3
  use, intrinsic :: iso_fortran_env
  use precision, only: fp

  use su3

  implicit none

  PRIVATe

  public gaugeconfiguration
  
  !> Module name
  character(len=22), parameter, public ::  modulename='gaugeconfiguration_su3'
  
  !>@brief Gauge configuration in temporal gauge
  !!@details Electric field is safe in lattice units as
  !!\f$\tilde E_{i,\vec{v}}=ga_ta_S(i)E_{i}(\vec{v})\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
  type GaugeConfiguration
     !> Gauge links \f$U_\mu(x)\f$
     complex(fp), allocatable, private :: links(:,:,:,:)
     !> Electric field
     real(fp),    allocatable, private :: efield(:,:,:)
     
   contains ! Member functions
     ! Allocation and deallocation
     procedure, public :: Allocate
     procedure, public :: Deallocate

     ! Boundary value communication
     procedure, public :: CommunicateBoundary
     procedure, public :: CommunicateBoundary_Links
     procedure, public :: CommunicateBoundary_Efield
     
     ! Return and setting of member variables
     procedure, public :: GetLink
     procedure, private:: SetLink
     procedure, private:: GetEfield
     procedure, private:: SetEfield

     ! Initialisation routines
     procedure, public :: ColdInit

     ! Gauss law deviation
     procedure, public :: GetDeviationFromGausslaw

     ! Energy
     procedure, public :: GetEnergy

     ! Field strength tensor and plaquettes
     procedure, public :: GetFieldStrengthTensor
     procedure, public :: GetPlaquette
     procedure, private :: GetSpatialPlaquette
     procedure, private :: GetTemporalPlaquette
  end type GaugeConfiguration


  
contains ! Module procedures
  
  !>@brief Allocation of gaugeconfiguration
  !!@details Allocates link variables and electric field in gaugeconfiguration
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
  pure subroutine Allocate(GaugeConf)
    use lattice, only: nDim, GetLocalLatticeSize_includingHalo
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(out) :: GaugeConf

    allocate(GaugeConf%links(nSUN,nSUN,nDim,&
         GetLocalLatticeSize_includingHalo()))
    allocate(GaugeConf%efield(nGen,nDim,&
         GetLocalLatticeSize_includingHalo()))
  end subroutine Allocate
  
  !>@brief Deallocation of gaugeconfiguration
  !!@details Deallocates link variables and electric field in gaugeconfiguration
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
  pure subroutine Deallocate(GaugeConf)
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(out) :: GaugeConf
    if(allocated(GaugeConf%links))  deallocate(GaugeConf%links)
    if(allocated(GaugeConf%efield)) deallocate(GaugeConf%efield)
  end subroutine Deallocate

  !>@brief Access routine to links
  !!@returns Link at given index
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
  pure function GetLink(GaugeConf,i,LatticeIndex)
    use matrixoperations, only: GetUnitmatrix
    use lattice, only: GetLocalIndex
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Direction
    integer(int8),  intent(in) :: i
    !> Lattice index
    integer(int64), intent(in) :: LatticeIndex
    !> Link variable
    complex(fp) :: GetLink(nSUN,nSUN)
    if(i/=0) then
       GetLink = GaugeConf%Links(:,:,i,GetLocalIndex(LatticeIndex))
    else
       GetLink = GetUnitmatrix(nSUN)
    end if
  end function GetLink

  !>@brief Setting routine for links
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
  pure subroutine SetLink(GaugeConf,i,LatticeIndex,Link)
    use lattice, only: GetLocalIndex
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(inout) :: GaugeConf
    !> Direction
    integer(int8),  intent(in) :: i
    !> Lattice index
    integer(int64), intent(in) :: LatticeIndex
    !> Value for link variable
    complex(fp),    intent(in) :: Link(nSUN,nSUN)
    GaugeConf%Links(:,:,i,GetLocalIndex(LatticeIndex)) = Link
  end subroutine SetLink

  !>@brief Access routine to the electric field
  !!@returns Efield at given index
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
  pure elemental function GetEfield(GaugeConf,a,i,LatticeIndex)
    use lattice, only: GetLocalIndex
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Generator index
    integer(int8),  intent(in) :: a
    !> Direction
    integer(int8),  intent(in) :: i
    !> Lattice index
    integer(int64), intent(in) :: LatticeIndex
    !> Efield
    real(fp) :: GetEfield
    if(i/=0) then
       GetEfield = GaugeConf%Efield(a,i,GetLocalIndex(LatticeIndex))
    else
       GetEfield = 0._fp
    end if
  end function GetEfield

  !>@brief Setting routine for links
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
  pure subroutine SetEfield(GaugeConf,a,i,LatticeIndex,Efield)
    use lattice, only: GetLocalIndex
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(inout) :: GaugeConf
    !> Generator index
    integer(int8),  intent(in) :: a
    !> Direction
    integer(int8),  intent(in) :: i
    !> Lattice index
    integer(int64), intent(in) :: LatticeIndex
    !> Value for link variable
    real(fp),       intent(in) :: Efield
    GaugeConf%Efield(a,i,GetLocalIndex(LatticeIndex)) = Efield
  end subroutine SetEfield

  !>@brief Communication routine for boundary values
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
  impure subroutine CommunicateBoundary(GaugeConf)
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(inout) :: GaugeConf
    call GaugeConf%CommunicateBoundary_Links
    call GaugeConf%CommunicateBoundary_Efield
  end subroutine CommunicateBoundary
  
  !>@brief Communication routine for boundary values
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
  impure subroutine CommunicateBoundary_Efield(GaugeConf)
    use HaloComm, only: HaloComm_CommunicateBoundary => CommunicateBoundary
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(inout) :: GaugeConf
    call HaloComm_CommunicateBoundary(GaugeConf%efield)
  end subroutine CommunicateBoundary_Efield
  
  !>@brief Communication routine for boundary values
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
  impure subroutine CommunicateBoundary_Links(GaugeConf)
    use HaloComm, only: HaloComm_CommunicateBoundary => CommunicateBoundary
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(inout) :: GaugeConf
    call HaloComm_CommunicateBoundary(GaugeConf%links)
  end subroutine CommunicateBoundary_Links

  !>@brief Initialises the links with unit matrices and electric field with zeroes
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
  pure subroutine ColdInit(GaugeConf)
    use matrixoperations, only: GetUnitMatrix
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(out) :: GaugeConf

    integer(int64) :: LocalIndex
    integer(int8)  :: i
    complex(fp) :: UnitMatrix(nSUN,nSUN)
    
    call GaugeConf%Allocate

    GaugeConf%Efield = 0
    UnitMatrix = GetUnitMatrix(nSUN)
    forall(i=1:size(GaugeConf%Links,3), LocalIndex=1:size(GaugeConf%Links,4))
       GaugeConf%Links(:,:,i,LocalIndex) = UnitMatrix
    end forall
  end subroutine ColdInit

  !>@brief Total deviation from Gauss-law of the total MPI-distributed configuration
  !!@returns Total deviation from Gauss-law of the total MPI-distributed configuration
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
  impure real(fp) function GetDeviationFromGaussLaw(GaugeConf)
    use mpiinterface, only: intmpi, GetRealSendType
    use mpi
    use lattice, only: nDim, GetLatticeSpacing,GetNeib,GetLocalIndex,&
         GetLocalLatticeIndices_allocatable, GetNeib
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf

    integer(intmpi) :: mpierr
    real(fp) :: local_contribution

    real(fp), dimension(ngen) :: efield, efield_neib
    complex(fp), dimension(nsun,nsun) :: Mefield, Mefield_neib, derivative, Link_neib
    
    integer(int8)  :: i, a
    integer(int64) :: is, neib, latticeindex
    integer(int64), allocatable :: LocalLatticeIndices(:)

    ! 1. Calculation of local contribution
    local_contribution = 0
    call GetLocalLatticeIndices_allocatable(LocalLatticeIndices)
    do concurrent (is=1:size(LocalLatticeIndices))
       LatticeIndex = LocalLatticeIndices(is)
       do concurrent(i=1_int8:ndim)
          efield = GaugeConf%GetEfield([1_int8:ngen],i,LatticeIndex)
          
          neib = GetNeib(-i,LatticeIndex)
          efield_neib = GaugeConf%GetEfield([1_int8:ngen],i,Neib)

          Mefield = GetAlgebraMatrix(efield)
          Mefield_neib = GetAlgebraMatrix(efield_neib)

          Link_neib = GaugeConf%GetLink(i,Neib)

          derivative = (Mefield &
               - matmul(matmul(&
               conjg(transpose(Link_Neib)),&
               Mefield_neib),&
               Link_neib))&
               /GetLatticeSpacing(i)**2/GetLatticeSpacing(0)
          do concurrent(a=1_int8:ngen)
             local_contribution = local_contribution &
                  + 2*Abs(Aimag(GetTraceWithGenerator(a,derivative)))
          end do
       end do
    end do

    ! 2. MPI-Sum over all partitions
    call MPI_ALLREDUCE(&
         local_contribution,&
         GetDeviationFromGaussLaw,&
         1_intmpi,&
         GetRealSendType(),&
         MPI_SUM,&
         MPI_COMM_WORLD,mpierr)
  end function GetDeviationFromGaussLaw

  !>@brief Returns field strength tensor
  !!@details: Field strength tensor in physical units \f$g F_{ij}(\vec{v})^aT^a\f$
  !!in order \f$O(a_i^2,a_j^2,a_i\cdot a_j)\f$
  !!@returns Spatial field strength tensor \f$F_{i,j,\vec{v}}\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 23.02.2019
  !!@version 1.0
  pure function GetFieldStrengthTensor(GaugeConf,i,j,LatticeIndex)
    use lattice, only: GetLatticeSpacing
    use matrixoperations, only: LogU
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Direction
    integer(int8),             intent(in) :: i,j
    !> Lattice index
    integer(int64),            intent(in) :: LatticeIndex
    !> Field strength tensor
    complex(fp) :: GetFieldStrengthTensor(nSUN,nSUN)

    complex(fp) :: Plaquette(nSUN,nSUN)

    Plaquette = GaugeConf%GetPlaquette(i,j,LatticeIndex)
    
    GetFieldStrengthTensor = &
         LogU(Plaquette)/GetLatticeSpacing(i)/GetLatticeSpacing(j)/cmplx(0,1,fp)
  end function GetFieldStrengthTensor

  !> @brief Returns spatial plaquette \f$U_{ij,\vec{v}} =U_{i,\vec{v}}\cdot U_{j,\vec{v}+\hat{i}}
  !! \cdot U_{i,\vec{v}+\hat{j}}^\dagger\cdot U_{j,\vec{v}}^\dagger\f$
  !! @returns Spatial plaquette
  !! \f$U_{i,j,v}=U_{i,v}\cdot U_{j,v+\hat{i}}\cdot U_{i,v+\hat{j}}^\dagger\cdot U_{j,v}^\dagger\f$
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 23.02.2019
  !! @version 1.0
  pure function GetSpatialPlaquette(GaugeConf,i,j,LatticeIndex)
    use lattice, only: GetNeib
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Spatial direction
    integer(int8),             intent(in) :: i,j
    !> Lattice index
    integer(int64),            intent(in) :: LatticeIndex
    !> Spatial plaquette
    complex(fp) :: GetSpatialPlaquette(nSUN,nSUN)

    complex(fp), dimension(nSUN,nSUN) :: link1,link2,link3,link4

    link1 = GaugeConf%GetLink(i,LatticeIndex)                             ! U_i(x)
    link2 = GaugeConf%GetLink(j,GetNeib(i,LatticeIndex))                  ! U_j(x+î)
    link3 = conjg(transpose(GaugeConf%GetLink(i,GetNeib(j,LatticeIndex))))! U_i(x+ĵ)†
    link4 = conjg(transpose(GaugeConf%GetLink(j,LatticeIndex)))           ! U_j(x)†

    GetSpatialPlaquette = matmul(matmul(matmul(&
         link1,&  ! U_i(x)
         link2),& ! U_j(x+î)
         link3),& ! U_i(x+ĵ)†
         link4)   ! U_j(x)†
  end function GetSpatialPlaquette

  !> @brief Returns temporal plaquette
  !! @details The temporal plaquette is computed via the chromo-electric field
  !! due to the knowledge of the evolution scheme. Using that
  !! \f$U_{i,\vec{v}}(t+a_t) =\exp\left(i \tilde{E}_{C,i,\vec{v}}(t)^aT^a\right)
  !! \cdot U_{i,\vec{v}}(t)\f$
  !! one can exactly compute the temporal plaquette with
  !! \f$U_{0j,\vec{v}}(t)=U_{j,\vec{v}}(t+a_t)
  !! \cdot U_{j,\vec{v}}(t)^\dagger
  !! =\exp\left(i \tilde{E}_{C,i,\vec{v}}(t)^aT^a\right)\f$
  !! @returns Temporal plaquette
  !! \f$U_{0,j,v}=U_{0,v}\cdot U_{j,v+\hat{0}}
  !! \cdot U_{0,v+\hat{j}}^\dagger\cdot U_{j,v}^\dagger
  !! = U_{j,v+\hat{0}}\cdot U_{j,v}^\dagger\f$
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 23.02.2019
  !! @version 1.0
  pure function GetTemporalPlaquette(GaugeConf,i,j,LatticeIndex)
    use matrixoperations, only: GetUnitMatrix
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Direction
    integer(int8),             intent(in) :: i,j
    !> Lattice index
    integer(int64),            intent(in) :: LatticeIndex
    !> Temporal plaquette
    complex(fp) :: GetTemporalPlaquette(nSUN,nSUN)

    real(fp) :: efield(ngen)
    
    if(     i==0 .and. j/=0 ) then
       efield = +GaugeConf%GetEfield([1_int8:ngen],j,LatticeIndex)
       GetTemporalPlaquette = GetGroupExp(efield)
    elseif( i/=0 .and. j==0 ) then
       efield = -GaugeConf%GetEfield([1_int8:ngen],i,LatticeIndex)
       GetTemporalPlaquette = GetGroupExp(efield)
    else
       GetTemporalPlaquette = GetUnitMatrix(nSUN)
    end if
  end function GetTemporalPlaquette
  
  !> @brief Returns plaquette \f$U_{jk,\vec{v}}
  !! =U_{j,\vec{v}}\cdot U_{k,\vec{v}+\hat{j}}
  !! \cdot U_{j,\vec{v}+\hat{k}}^\dagger\cdot U_{k,\vec{v}}^\dagger\f$
  !! @returns Plaquette
  !! \f$U_{j,k,v}=U_{j,v}\cdot U_{k,v+\hat{j}}
  !! \cdot U_{j,v+\hat{k}}^\dagger\cdot U_{k,v}^\dagger\f$
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 23.02.2019
  !! @version 1.0
  pure function GetPlaquette(GaugeConf,i,j,LatticeIndex)
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Direction
    integer(int8),             intent(in) :: i,j
    !> Lattice index
    integer(int64),            intent(in) :: LatticeIndex
    !> Temporal plaquette
    complex(fp) :: GetPlaquette(nSUN,nSUN)
    if(i==0 .or. j==0) then
       GetPlaquette = GaugeConf%GetTemporalPlaquette(i,j,LatticeIndex)
    else
       GetPlaquette = GaugeConf%GetSpatialPlaquette(i,j,LatticeIndex)
    end if
  end function GetPlaquette

  !>@brief Energy functional
  !!@returns Energy functional
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
  impure real(fp) function GetEnergy(GaugeConf)
    use precision, only: fp
    use matrixoperations, only: GetTrace
    use lattice, only: nDim, GetLatticeSpacing, GetLocalLatticeIndices_allocatable
    use mpiinterface, only: intmpi, GetRealSendType
    use mpi
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf

    integer(intmpi) :: mpierr
    real(fp) :: local_contribution

    real(fp), dimension(ngen) :: Efield
    complex(fp), dimension(nsun,nsun) :: Plaquette

    real(fp) :: PotentialTerm
    integer(int8)  :: i, j
    integer(int64) :: is, latticeindex
    integer(int64), allocatable :: LocalLatticeIndices(:)

    ! 1. Calculation of local contribution
    local_contribution = 0
    call GetLocalLatticeIndices_allocatable(LocalLatticeIndices)
    do concurrent (is=1:size(LocalLatticeIndices))
       LatticeIndex = LocalLatticeIndices(is)
       do concurrent(i=1_int8:ndim)
          efield = GaugeConf%GetEfield([1_int8:ngen],i,LatticeIndex)
          
          local_contribution = local_contribution &
               ! Electric energy
               + sum(efield**2)/GetLatticeSpacing(i)/GetLatticeSpacing(0_int8)/2
          ! Magnetic energy
          do concurrent(j=i+1_int8:ndim)
             Plaquette = GaugeConf%GetPlaquette(i,j,LatticeIndex)
             PotentialTerm = (1-real(GetTrace(Plaquette),fp)/nSUN)&
                  /GetLatticeSpacing(i)/GetLatticeSpacing(j)

             Local_contribution = Local_Contribution &
                  + 2*nSUN*PotentialTerm
          end do
       end do
    end do

    ! 2. MPI-Sum over all partitions
    call MPI_ALLREDUCE(&
         local_contribution,&
         GetEnergy,&
         1_intmpi,&
         GetRealSendType(),&
         MPI_SUM,&
         MPI_COMM_WORLD,mpierr)
  end function GetEnergy
end module gaugeconfiguration_su3
