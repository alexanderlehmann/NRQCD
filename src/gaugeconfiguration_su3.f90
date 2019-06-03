!----------------------------------------------------------------------
! RTLQCD, Real-Time-Lattice-QCD Simulation of Gauge Fields
!----------------------------------------------------------------------
!
!>@brief SU(3)-gauge-link-configuration in temporal gauge
!!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!!@date 22.02.2019
!!@version 1.0
!----------------------------------------------------------------------
module gaugeconfiguration_su3
  use, intrinsic :: iso_fortran_env
  use precision, only: fp
  
  use su3

  implicit none

  PRIVATE

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
     procedure, public :: GetLink_G
     procedure, public :: GetLink_M
     procedure, private:: SetLink_G
     procedure, private:: SetLink_M
     procedure, private:: GetEfield_G
     procedure, private:: GetEfield_M
     procedure, private:: SetEfield_G
     procedure, private:: SetEfield_M
     
     ! Initialisation routines
     procedure, public :: EquilibriumInit
     procedure, public :: SemicoldInit
     procedure, public :: HotInit
     procedure, public :: ColdInit
     procedure, public :: TransversePolarisedOccupiedInit_Box
     procedure, private :: TransversePolarisedOccupiedInit
     procedure, private :: TransversePolarisedOccupiedInit_specificSender

     ! Gauss law deviation
     procedure, public :: GetDeviationFromGausslaw

     ! Energy
     procedure, public :: GetEnergy

     ! Field strength tensor and plaquettes
     procedure, public :: GetFieldStrengthTensor_G
     procedure, public :: GetPlaquette_G
     procedure, private :: GetSpatialPlaquette_G
     procedure, private :: GetTemporalPlaquette_G
     
     procedure, public :: GetFieldStrengthTensor_M
     procedure, public :: GetPlaquette_M
     procedure, private :: GetSpatialPlaquette_M
     procedure, private :: GetTemporalPlaquette_M

     ! Return of semi-physical fields
     procedure, public :: GetGaugefield_AlgebraCoordinate
     procedure, public :: GetGaugefield_AlgebraCoordinates
     generic, public :: GetElectricField => &
          GetElectricField_AlgebraCoordinate, GetElectricField_AlgebraMatrix
     procedure, private :: GetElectricField_AlgebraMatrix
     procedure, private :: GetElectricField_AlgebraCoordinate
     
     ! AA- and EE-correlator
     procedure, public :: GetTransverseAACorrelator
     procedure, public :: GetTransverseEECorrelator
     !procedure, public :: GetAACorrelator

     ! Update routines
     generic, public :: Update => Update_Leapfrog
     procedure, private :: Update_Leapfrog
     procedure, private :: Update_Efield_Leapfrog
     procedure, private :: Update_Links_Leapfrog
     
     ! Gauge fixing
     procedure, public :: CoulombGaugefixing
     procedure, private :: CoulombGaugefixing_Links

     ! Utility routines for gaugixing
     procedure, private :: GetDivergenceOfGaugefield_M
     
  end type GaugeConfiguration


  
contains ! Module procedures
  
  !>@brief Allocation of gaugeconfiguration
  !!@details Allocates link variables and electric field in gaugeconfiguration
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
 pure subroutine Allocate(GaugeConf)
    use lattice, only: nDim, GetMemorySize
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(out) :: GaugeConf

    allocate(GaugeConf%links(nSUN,nSUN,nDim,&
         GetMemorySize()))
    allocate(GaugeConf%efield(nGen,nDim,&
         GetMemorySize()))
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
  !!@date 28.02.2019
  !!@version 1.1
 pure function GetLink_G(GaugeConf,i,LatticeIndex)
    use lattice, only: GetMemoryIndex
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Direction
    integer(int8),  intent(in) :: i
    !> Lattice index
    integer(int64), intent(in) :: LatticeIndex
    !> Link variable
    complex(fp) :: GetLink_G(nSUN,nSUN)
    GetLink_G = GaugeConf%GetLink_M(i,GetMemoryIndex(LatticeIndex))
  end function GetLink_G

  !>@brief Access routine to links
  !!@returns Link at given memory index
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 28.02.2019
  !!@version 1.0
 pure function GetLink_M(GaugeConf,i,MemoryIndex)
    use matrixoperations, only: GetUnitmatrix
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Direction
    integer(int8),  intent(in) :: i
    !> Memory index
    integer(int64), intent(in) :: MemoryIndex
    !> Link variable
    complex(fp) :: GetLink_M(nSUN,nSUN)
    if(i/=0) then
       GetLink_M = GaugeConf%Links(:,:,i,MemoryIndex)
    else
       GetLink_M = GetUnitmatrix(nSUN)
    end if
  end function GetLink_M

  !>@brief Setting routine for links
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 28.02.2019
  !!@version 1.1
 pure subroutine SetLink_G(GaugeConf,i,LatticeIndex,Link)
    use lattice, only: GetMemoryIndex
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(inout) :: GaugeConf
    !> Direction
    integer(int8),  intent(in) :: i
    !> Lattice index
    integer(int64), intent(in) :: LatticeIndex
    !> Value for link variable
    complex(fp),    intent(in) :: Link(nSUN,nSUN)
    GaugeConf%Links(:,:,i,GetMemoryIndex(LatticeIndex)) = Link
  end subroutine SetLink_G

  !>@brief Setting routine for links
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 28.02.2019
  !!@version 1.0
 pure subroutine SetLink_M(GaugeConf,i,MemoryIndex,Link)
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(inout) :: GaugeConf
    !> Direction
    integer(int8),  intent(in) :: i
    !> Memory index
    integer(int64), intent(in) :: MemoryIndex
    !> Value for link variable
    complex(fp),    intent(in) :: Link(nSUN,nSUN)
    GaugeConf%Links(:,:,i,MemoryIndex) = Link
  end subroutine SetLink_M

  !>@brief Access routine to the electric field
  !!@returns Efield at given index
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 28.02.2019
  !!@version 1.0
 pure elemental real(fp) function GetEfield_G(GaugeConf,a,i,LatticeIndex)
    use lattice, only: GetMemoryIndex
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Generator index
    integer(int8),  intent(in) :: a
    !> Direction
    integer(int8),  intent(in) :: i
    !> Lattice index
    integer(int64), intent(in) :: LatticeIndex
    GetEfield_G = GaugeConf%GetEfield_M(a,i,GetMemoryIndex(LatticeIndex))
  end function GetEfield_G

  !>@brief Access routine to the electric field
  !!@returns Efield at given memory index
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 28.02.2019
  !!@version 1.0
 pure elemental real(fp) function GetEfield_M(GaugeConf,a,i,MemoryIndex)
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Generator index
    integer(int8),  intent(in) :: a
    !> Direction
    integer(int8),  intent(in) :: i
    !> Memory index
    integer(int64), intent(in) :: MemoryIndex
    GetEfield_M = GaugeConf%Efield(a,i,MemoryIndex)
  end function GetEfield_M

  !>@brief Setting routine for links
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
 pure subroutine SetEfield_G(GaugeConf,a,i,LatticeIndex,Efield)
    use lattice, only: GetMemoryIndex
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
    GaugeConf%Efield(a,i,GetMemoryIndex(LatticeIndex)) = Efield
  end subroutine SetEfield_G
  
  !>@brief Setting routine for links
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
 pure subroutine SetEfield_M(GaugeConf,a,i,MemoryIndex,Efield)
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(inout) :: GaugeConf
    !> Generator index
    integer(int8),  intent(in) :: a
    !> Direction
    integer(int8),  intent(in) :: i
    !> Memory index
    integer(int64), intent(in) :: MemoryIndex
    !> Value for link variable
    real(fp),       intent(in) :: Efield
    GaugeConf%Efield(a,i,MemoryIndex) = Efield
  end subroutine SetEfield_M

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

    integer(int64) :: MemoryIndex
    integer(int8)  :: i
    
    call GaugeConf%Allocate

    GaugeConf%Efield = 0
    forall(i=1:size(GaugeConf%Links,3), MemoryIndex=1:size(GaugeConf%Links,4))
       GaugeConf%Links(:,:,i,MemoryIndex) = GetUnitMatrix(nSUN)
    end forall
  end subroutine ColdInit

  impure subroutine EquilibriumInit(GaugeConf,Beta,nefieldinit,nequilibrium,MeasureEnergy,filename)
    use lattice, only: GetLatticeSpacing, GetMemorySize, GetProc_M,&
         ndim, GetNeib_M, getvolume
    use random, only: GetRandomNormalCmplx
    use mpiinterface, only: ThisProc, mpistop
    use io
    use matrixoperations
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(out) :: GaugeConf
    !> \f$\beta_{\text{L}}\f$
    real(fp), intent(in) :: beta
    integer(int64), intent(in) :: nefieldinit
    integer(int64), intent(in) :: nequilibrium

    logical, optional, intent(in) :: MeasureEnergy
    character(len=*), optional, intent(in) :: filename

    real(fp) :: sigma
    
    integer(int8) :: k, a
    integer(int64) :: MemoryIndex
    complex(fp) :: r(ngen)

    integer :: iefieldinit
    integer :: iequibstep

    real(fp) :: blub
    real(fp) :: energy, deviation
    real(fp), parameter :: kappa=0.12_fp
    real(fp), dimension(ngen) :: ec
    complex(fp), dimension(nsun,nsun) :: U, G, Gneib, E, term, uconjg


    integer(int8) :: fileid

    integer(int64) :: i
    
    sigma = sqrt(nsun/beta)

    call GaugeConf%ColdInit

    ! Start thermalizing gauge links

    if(present(MeasureEnergy).and.MeasureEnergy) then
       if(MeasureEnergy.and..not.present(filename)) then
          call MPIStop('Filename for energy output missing')
       end if
       if(thisproc()==0) &
            fileID = OpenFile(filename=filename,&
            st='REPLACE',fm='FORMATTED',act='WRITE')
    end if
    do iefieldinit=1,nefieldinit
       ! Redrawing E-field
       gaugeconf%efield = 0
       do MemoryIndex=1,GetMemorySize()
          if(ThisProc()==GetProc_M(MemoryIndex)) then
             do k=1,nDim
                r = sigma*GetRandomNormalCmplx(int(ngen,int64))

                GaugeConf%Efield(:,k,MemoryIndex) = real(r,fp)
             end do
          end if
       end do
       
       call GaugeConf%CommunicateBoundary()
       
       deviation = GetGdeviation(GaugeConf)
       
       ! Projection of the electric field
       i = 0
       do while(deviation>1E-12)
          i = i + 1
          do MemoryIndex=1,GetMemorySize()
             if(ThisProc()==GetProc_M(MemoryIndex)) then
                G = GetG(GaugeConf,MemoryIndex)
                do k=1,nDim
                   ec = GaugeConf%GetEfield_M([1_int8:ngen],k,MemoryIndex)
                   U = GaugeConf%Links(:,:,k,MemoryIndex)
                   E = GetAlgebraMatrix(ec)

                   Gneib = GetG(GaugeConf,GetNeib_M(+k,MemoryIndex))
                   
                   term = &
                        kappa*(matmul(matmul(matmul(&
                        U,G),conjg(transpose(u))),Gneib)&
                        - G) &
                        + E

                   do a=1,ngen
                      GaugeConf%Efield(a,k,MemoryIndex) &
                           = GetAlgebraCoordinate(a,term)
                      
                   end do
                end do
             end if
          end do
          call GaugeConf%CommunicateBoundary()
          
          deviation = GetGdeviation(GaugeConf)/GetVolume()
          !if(ThisProc()==0) write(output_unit,*) deviation
       end do
       
       ! Evolution of the gauge links
       do iequibstep=1,nequilibrium
          call GaugeConf%Update
       end do

       blub = GaugeConf%GetDeviationFromGaussLaw()
       if(ThisProc()==0) print*,blub
       
       ! Print energy to terminal
       if(present(MeasureEnergy).and.MeasureEnergy) then
          energy = GaugeConf%GetEnergy()
          if(ThisProc()==0) then
             write(fileid,*) iefieldinit,energy
             write(output_unit,*) iefieldinit,energy
          end if
       end if
    end do
    if(present(MeasureEnergy).and.MeasureEnergy) then
       if(ThisProc()==0) then
          call Closefile(fileid)
       end if
    end if
  contains
    impure real(fp) function GetGdeviation(GaugeConf)
    use mpiinterface, only: intmpi, GetRealSendType, ThisProc
    use mpi
    use lattice, only: nDim, GetLatticeSpacing,GetNeib_M,GetLatticeIndex_M, GetMemorySize, GetProc_M
    implicit none
    !> Gauge configuration
    type(GaugeConfiguration), intent(in) :: GaugeConf

    integer(intmpi) :: mpierr
    real(fp) :: local_contribution

    real(fp), dimension(ngen) :: efield, efield_neib
    complex(fp), dimension(nsun,nsun) :: Mefield, Mefield_neib, derivative, Link_neib, G
    
    integer(int8)  :: i, a
    integer(int64) :: MemoryIndex, neib

    ! 1. Calculation of local contribution
    local_contribution = 0
    
    do concurrent(MemoryIndex=1:GetMemorySize(),i=1:ndim,ThisProc()==GetProc_M(MemoryIndex))
       efield = GaugeConf%GetEfield_M([1_int8:ngen],i,MemoryIndex)

       neib = GetNeib_M(-i,MemoryIndex)
       efield_neib = GaugeConf%GetEfield_M([1_int8:ngen],i,Neib)

       Mefield = GetAlgebraMatrix(efield)
       Mefield_neib = GetAlgebraMatrix(efield_neib)

       Link_neib = GaugeConf%GetLink_M(i,Neib)

       G = GetG(GaugeConf,MemoryIndex)
       do concurrent(a=1_int8:ngen)
          local_contribution = local_contribution &
               + 2*Abs(GetTraceWithGenerator(a,G))
       end do
    end do

    ! 2. MPI-Sum over all partitions
    call MPI_ALLREDUCE(&
         local_contribution,&
         GetGdeviation,&
         1_intmpi,&
         GetRealSendType(),&
         MPI_SUM,&
         MPI_COMM_WORLD,mpierr)
    end function GetGdeviation
    
    pure function GetG(GaugeConf,MemoryIndex)
      use lattice, only: nDim, GetNeib_M
      implicit none
      type(GaugeConfiguration), intent(in) :: GaugeConf
      integer(int64), intent(in) :: MemoryIndex

      complex(fp), dimension(nsun,nsun) :: GetG, Uneib, E, Eneib
      real(fp), dimension(ngen) :: ec,ecneib

      integer(int8) :: i
      integer(int64) :: neib

      GetG = 0
      do i=1,ndim
         neib = GetNeib_m(-i,MemoryIndex)

         Uneib = GaugeConf%Links(:,:,i,neib)
         
         Eneib = GetAlgebraMatrix(GaugeConf%Efield(:,i,neib))

         E = GetAlgebraMatrix(GaugeConf%Efield(:,i,MemoryIndex))
         
         GetG = GetG &
              + E - matmul(matmul(conjg(transpose(Uneib)),Eneib),Uneib)
      end do
    end function GetG
  end subroutine EquilibriumInit
  
  !>@brief Initialises the links with unit matrices on almost all lattice sites
  !! and electric field with zeroes
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 06.03.2019
  !!@version 1.0
  impure subroutine SemiColdInit(GaugeConf)
    use lattice, only: GetProc_G, GetMemoryIndex
    use matrixoperations, only: GetUnitMatrix
    use mpiinterface, only: ThisProc
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(out) :: GaugeConf

    integer(int64) :: MemoryIndex, LatticeIndex
    integer(int8)  :: i

    real(fp) :: afield(ngen)
    complex(fp) :: link(nsun,nsun)


    real(fp) :: retr_plaquette
    complex(fp) :: plaquette(nsun,nsun)
    
    call GaugeConf%ColdInit

    i = 1
    LatticeIndex = 1

    if(GetProc_G(LatticeIndex)==ThisProc()) then
       afield = 0
       afield(1) = 0.1_fp
       afield(2) = 0.2_fp
       afield(3) = 1.3_fp
       afield(4) = 2.4_fp
       afield(5) = 3.1_fp
       afield(6) = -1.1_fp
       afield(7) = -2.1_fp
       afield(8) = -0.1_fp
       link = GetGroupExp(afield)
       GaugeConf%Links(:,:,i,GetMemoryIndex(LatticeIndex)) = link
       !call GaugeConf%SetLink_G(i,LatticeIndex,link)
    end if

    call GaugeConf%CommunicateBoundary
  end subroutine SemiColdInit

  !>@brief Initialises the links and efield with random entries (still unitary links)
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 28.02.2019
  !!@version 1.0
  impure subroutine HotInit(GaugeConf)
    use mpiinterface, only: ThisProc
    use random, only: GetRandomUniformReal
    use lattice, only: ndim, GetProc_M, GetMemorySize, GetLatticeSpacing
    implicit none
    class(GaugeConfiguration), intent(out) :: GaugeConf

    integer(int64) :: MemoryIndex
    integer(int8)  :: i
    real(fp) :: r(ngen)

    call GaugeConf%Allocate
    
    do MemoryIndex=1,GetMemorySize()
       if(ThisProc()==GetProc_M(MemoryIndex)) then
          do i=1,ndim
             ! Link
             r = GetRandomUniformReal(int(ngen,int64))*GetLatticeSpacing(i)
             GaugeConf%Links(:,:,i,MemoryIndex) = GetGroupExp(r)

             ! E-field
             r = GetRandomUniformReal(int(ngen,int64))*GetLatticeSpacing(i)
             GaugeConf%Efield(:,i,MemoryIndex) = r
          end do
       end if
    end do

    call GaugeConf%CommunicateBoundary
  end subroutine HotInit

  !>@brief Initialises the configuration as in the shape of a box highly occupied,
  !! transverse polarised fields
  !!@details Initialises the configuration with gaussian initial conditions which fills
  !! the gluon occupation up to a saturation scale \f$q_s\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 24.02.2019
  !!@version 1.0
  impure subroutine TransversePolarisedOccupiedInit_Box(GaugeConf,SaturationScale,Amplitude,Coupling,&
       aa_correlator,ee_correlator)
    use, intrinsic :: iso_fortran_env
    use precision, only: fp
    use lattice, only: GetNorm2Momentum_M,GetMemorySize
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(out) :: GaugeConf
    !> Saturation scale \f$q_s\f$
    real(fp),                  intent(in)  :: SaturationScale
    !> Occupation of the box for \f$p<q_s\f$
    real(fp),                  intent(in)  :: Amplitude
    !> Coupling
    real(fp),                  intent(in)  :: Coupling
    !> AA-Correlator
    real(fp), optional, allocatable, intent(out) :: aa_correlator(:)
    !> EE-Correlator
    real(fp), optional, allocatable, intent(out) :: ee_correlator(:)
    
    real(fp), allocatable :: Occupation(:)
    integer(int64) :: MemoryIndex
    
    allocate(Occupation(GetMemorySize()))
    
    forall(MemoryIndex=1:GetMemorySize())
       Occupation(MemoryIndex)&
            = GetBoxOccupation(&
            GetNorm2Momentum_M(MemoryIndex),& !|p|
            SaturationScale,&
            Amplitude,&
            Coupling)
    end forall

    if(present(aa_correlator).and.present(ee_correlator)) then
       call GaugeConf%TransversePolarisedOccupiedInit(Occupation,aa_correlator,ee_correlator)
    elseif(present(aa_correlator).and..not.present(ee_correlator)) then
       call GaugeConf%TransversePolarisedOccupiedInit(Occupation,aa_correlator=aa_correlator)
    elseif(.not.present(aa_correlator).and.present(ee_correlator)) then
       call GaugeConf%TransversePolarisedOccupiedInit(Occupation,ee_correlator=ee_correlator)
    else
       call GaugeConf%TransversePolarisedOccupiedInit(Occupation)
    end if

    deallocate(Occupation)
    
  end subroutine TransversePolarisedOccupiedInit_Box

  !>@brief Initialises the configuration as highly occupied, transverse polarised fields
  !!@details Initialises the configuration according to a given gluon occupation
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 24.02.2019
  !!@version 1.0
  impure subroutine TransversePolarisedOccupiedInit_specificSender(GaugeConf,Occupation,aa_correlator,ee_correlator)
    use, intrinsic :: iso_fortran_env
    use precision, only: fp
    use tolerances, only: GetZeroTol
    use mpiinterface, only: intmpi, ThisProc, mpistop, SyncAll
    use lattice, only: nDim, GetMemorySize, GetMemoryIndex, GetLatticeSize,&
         GetProc_G, GetVolume, GetPolarisationVectors, GetNorm2Momentum_G, GetMomentum_G,&
         GetLatticeSpacing, GetLatticeIndex_M, GetLocalLatticeSize
    use random, only: IsModuleInitialised_random=>IsModuleInitialised,&
         GetRandomNormalCmplx_specificProcess,&
         modulename_random=>modulename
    use xpfft, only: p2x, x2p
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(out) :: GaugeConf
    !> Gauge particle occupation
    real(fp),                  intent(in)  :: Occupation(:)
    !> AA-Correlator
    real(fp), optional, allocatable, intent(out) :: aa_correlator(:)
    !> EE-Correlator
    real(fp), optional, allocatable, intent(out) :: ee_correlator(:)

    ! Gauge and electric field
    complex(fp), allocatable :: afield(:,:,:), efield(:,:,:)

    ! Prefactors
    complex(fp) :: prefactor_afield, prefactor_efield

    ! Number of transverse polarisations
    integer(int8), parameter :: nPol_transverse = ndim - 1_int8
    ! Transverse polarisation vectors
    complex(fp), dimension(nDim,nDim) :: PolarisationVectors
    ! Random numbers for fluctuations in mode decomposition
    complex(real64) :: r_afield(npol_transverse), r_efield(npol_transverse)
    ! Momentum
    complex(fp), dimension(nDim) :: Momentum(nDim)
    ! Norm of momentum
    real(fp) :: MomentumNorm

    ! Lattice
    integer(int64) :: LatticeIndex, MemoryIndex, is
    integer(int8) :: i,a
    ! MPI
    integer(intmpi) :: RecvProc
    integer(intmpi), parameter :: SendProc=0_intmpi

    ! Link
    complex(fp) :: Link(nSUN,nSUN)
    ! E-field
    real(fp) :: Efield_site(ngen)
    ! A-field
    real(fp) :: Afield_site(ngen)
    character(len=105) :: errormessage

    complex(fp) :: afield_p(ndim),efield_p(ndim)
    
    ! Check if random numbers are initialised
    if(.not.isModuleInitialised_random()) then
       errormessage = 'Error in TransversePolarisedOccupiedInit of '//modulename&
               //':'// modulename_random //' is not initialised.'
       call MPISTOP(errormessage)
    end if
    
    ! ..--** START: Allocations **--..
    allocate(afield(GetMemorySize(),nGen,nDim))
    allocate(efield(GetMemorySize(),nGen,nDim))
    ! ..--**  END : Allocations **--..

    ! ..--** START: Drawing random numbers and setting the field in p-space **--..
    !
    !        Note : All processes are doing this simultaneously together
    !               in order to ensure that the results do not depend on
    !               the number of used processes
    !
    GlobalLattice: do LatticeIndex=1,GetLatticeSize()
       RecvProc = GetProc_G(LatticeIndex)

       if(ThisProc()==RecvProc .or. ThisProc()==SendProc) then

          MomentumNorm = GetNorm2Momentum_G(LatticeIndex)
          if(MomentumNorm > GetZeroTol()) then
             if(ThisProc()==RecvProc) then
                Momentum = GetMomentum_G(LatticeIndex)
                call GetPolarisationVectors(Momentum,PolarisationVectors)

                MemoryIndex = GetMemoryIndex(LatticeIndex)

                prefactor_afield = &
                     sqrt(GetVolume()*Occupation(MemoryIndex)/MomentumNorm/2)

                prefactor_efield = &
                     cmplx(0,sqrt(GetVolume()*Occupation(MemoryIndex)*MomentumNorm/2),fp)
             end if ! is recieving process

             generator: do a=1_int8,ngen
                r_afield = GetRandomNormalCmplx_specificProcess(int(ngen,int64),SendProc,RecvProc)

                if(ThisProc()==RecvProc) then
                   r_efield = r_afield
                   
                   dims: do concurrent(i=1_int8:ndim)
                      afield(MemoryIndex,a,i) = &
                           prefactor_afield*(&
                                ! First transversal polarisation
                           + r_afield(1)*PolarisationVectors(i,1) &
                                ! Second transversal polarisation
                           + r_afield(2)*PolarisationVectors(i,2) )

                      efield(MemoryIndex,a,i) = &
                           prefactor_efield*(&
                                ! First transversal polarisation
                           + r_efield(1)*PolarisationVectors(i,1) &
                                ! Second transversal polarisation
                           + r_efield(2)*PolarisationVectors(i,2) )
                   end do dims
                end if ! is recieving process
             end do generator ! Generators
          
          else !p=0
             if(ThisProc()==RecvProc) then
                MemoryIndex = GetMemoryIndex(LatticeIndex)
                afield(MemoryIndex,:,:) = 0
                efield(MemoryIndex,:,:) = 0
             end if !recieving process
          end if !p>0
       end if !recieving or sending process

       call SyncAll
    end do GlobalLattice
    ! ..--**  END : Drawing random numbers and setting the field in p-space **--..

    ! Symmetrisation A(p) = A(-p)* in order to make A(x) real
    do i=1,ndim
       do a=1,ngen
          call SymmetriseInPspace(afield(:,a,i))
          call SymmetriseInPspace(efield(:,a,i))
       end do !a
    end do !i

    if(present(aa_correlator)) then
       allocate(aa_correlator(GetLocalLatticeSize()))
       is=0
       do MemoryIndex=1,GetMemorySize()
          LatticeIndex = GetLatticeIndex_M(MemoryIndex)
          if(GetProc_G(LatticeIndex)==ThisProc()) then
             is = is+1
             aa_correlator(is) = 0
             do concurrent(a=1:ngen)
                afield_p = afield(MemoryIndex,a,:)
                aa_correlator(is) = aa_correlator(is) &
                     + GetTransverseField_usingProjector(afield_p,LatticeIndex)&
                     /GetVolume()/nGen
             end do
          end if
       end do
    end if
    if(present(ee_correlator)) then
       allocate(ee_correlator(GetLocalLatticeSize()))
       is=0
       do MemoryIndex=1,GetMemorySize()
          LatticeIndex = GetLatticeIndex_M(MemoryIndex)
          if(GetProc_G(LatticeIndex)==ThisProc()) then
             is = is+1
             ee_correlator(is) = 0
             do concurrent(a=1:ngen)
                efield_p = efield(MemoryIndex,a,:)
                ee_correlator(is) = ee_correlator(is) &
                     + GetTransverseField_usingProjector(efield_p,LatticeIndex)&
                     /GetVolume()/nGen
             end do
          end if
       end do
    end if
    
    !..--** START: FFT p-->x **--..
    do i=1,ndim
       do a=1,ngen
          call p2x(afield(:,a,i))
          call p2x(efield(:,a,i))
       end do !a
    end do !i
    !..--**  END : FFT p-->x **--..

    
    !..--** START: Writing fields to configuration **--..
    call GaugeConf%Allocate
    !do concurrent(&
    !     is=1:size(LocalLatticeIndices),&
    !     i =1:ndim)
    do concurrent(MemoryIndex=1:GetMemorySize(), i=1:ndim)
       LatticeIndex = GetLatticeIndex_M(MemoryIndex)
       if(GetProc_G(LatticeIndex)==ThisProc()) then
          ! Link
          afield_site = real(afield(MemoryIndex,:,i),fp) &
                                ! Translation to lattice units
               *GetLatticeSpacing(i)
          GaugeConf%Links(:,:,i,MemoryIndex) = GetGroupExp(afield_site)

          ! E-field
          efield_site = real(efield(MemoryIndex,:,i),fp) &
                                ! Translation to lattice units
               *GetLatticeSpacing(i)
          GaugeConf%Efield(:,i,MemoryIndex) = efield_site
       end if
    end do
    !..--**  END : Writing fields to configuration **--..
    call GaugeConf%CommunicateBoundary()
  contains
    impure subroutine SymmetriseInPspace(field)
      use precision, only: fp
      use lattice, only: nDim, GetProc_G, GetNegativeLatticeIndex, GetLatticeSize
      use mpiinterface, only: ThisProc, SyncAll, intmpi, GetComplexSendType
      use mpi
      implicit none
      complex(fp), intent(inout) :: field(:)


      ! MPI
      complex(fp) :: a
      integer(intmpi) :: status(mpi_status_size)
      integer(intmpi), parameter :: buffersize=1_intmpi
      integer(intmpi) :: dest, source
      integer(intmpi) :: tag
      integer(intmpi) :: mpierr

      ! Indices
      integer(intmpi) :: PositiveProc, NegativeProc
      integer(int64)  :: PositiveLatticeIndex, NegativeLatticeIndex
      integer(int64)  :: PositiveMemoryIndex, NegativeMemoryIndex
      integer(int64)  :: LatticeSize

      LatticeSize = GetLatticeSize()

      do PositiveLatticeIndex=1,LatticeSize
         ! Getting negative lattice index, corresponding to -p
         NegativeLatticeIndex = GetNegativeLatticeIndex(PositiveLatticeIndex)

         ! Getting process numbers for positive and negative momentum
         PositiveProc = GetProc_g(PositiveLatticeIndex)
         NegativeProc = GetProc_G(NegativeLatticeIndex)

         ! Setting field such that f(p) = f(-p)*

         if(PositiveProc == NegativeProc) then
            ! positive and negative process are the same --> No MPI-communication needed
            if(ThisProc()==PositiveProc) then
               PositiveMemoryIndex = GetMemoryIndex(PositiveLatticeIndex)
               NegativeMemoryIndex = GetMemoryIndex(NegativeLatticeIndex)
               
               a = field(PositiveMemoryIndex)
               if(PositiveMemoryIndex==NegativeMemoryIndex) then
                  field(PositiveMemoryIndex) = real(a,fp)
               else
                  field(PositiveMemoryIndex) = a
                  field(NegativeMemoryIndex) = conjg(a)
               end if
            end if
         else 
            ! positive and negative process are not the same --> MPI-communication
            tag  = PositiveLatticeIndex
            source = PositiveProc
            dest = NegativeProc
            
            if(ThisProc()==PositiveProc) then
               PositiveMemoryIndex = GetMemoryIndex(PositiveLatticeIndex)
               a = field(PositiveMemoryIndex)

               call MPI_Send(&
                    a,                   & ! What to send
                    buffersize,          & ! How many points to send
                    GetComplexSendType(),& ! What type to send
                    dest,                & ! Destination
                    tag,                 & ! Tag
                    MPI_COMM_WORLD,      & ! Communicator
                    mpierr)                ! Error code
            elseif(ThisProc()==NegativeProc) then
               call MPI_Recv(&
                    a,                   & ! What to send
                    buffersize,          & ! How many points to send
                    GetComplexSendType(),& ! What type to send
                    source,              & ! Destination
                    tag,                 & ! Tag
                    MPI_COMM_WORLD,      & ! Communicator
                    status,              & ! Status
                    mpierr)                ! Error code
               
               NegativeMemoryIndex = GetMemoryIndex(NegativeLatticeIndex)
               field(NegativeMemoryIndex) = conjg(a)
            else
               ! do nothing
            end if
         end if

         call SyncAll
      end do
    end subroutine SymmetriseInPspace
  end subroutine TransversePolarisedOccupiedInit_SpecificSender

  !>@brief Initialises the configuration as highly occupied, transverse polarised fields
  !!@details Initialises the configuration according to a given gluon occupation
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 24.02.2019
  !!@version 1.0
  impure subroutine TransversePolarisedOccupiedInit(GaugeConf,Occupation,aa_correlator,ee_correlator)
    use, intrinsic :: iso_fortran_env
    use precision, only: fp
    use tolerances, only: GetZeroTol
    use mpiinterface, only: intmpi, ThisProc, mpistop
    use lattice, only: nDim, GetMemorySize,&
         GetProc_M, GetVolume, GetPolarisationVectors, GetNorm2Momentum_M, GetMomentum_M,&
         GetLatticeSpacing, GetLatticeSize, GetLocalLatticeSize, GetLatticeIndex_M, GetProc_G
    use random, only: IsModuleInitialised_random=>IsModuleInitialised,&
         GetRandomNormalCmplx,&
         modulename_random=>modulename
    use xpfft, only: p2x, x2p
    use mathconstants, only: pi
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(out) :: GaugeConf
    !> Gauge particle occupation
    real(fp),                  intent(in)  :: Occupation(:)
    !> AA-Correlator
    real(fp), optional, allocatable, intent(out) :: aa_correlator(:)
    !> EE-Correlator
    real(fp), optional, allocatable, intent(out) :: ee_correlator(:)

    ! Gauge and electric field
    complex(fp), allocatable :: afield(:,:,:), efield(:,:,:)

    ! Prefactors
    complex(fp) :: prefactor_afield, prefactor_efield

    ! Number of transverse polarisations
    integer(int8), parameter :: nPol_transverse = ndim - 1_int8
    ! Transverse polarisation vectors
    complex(fp), dimension(nDim,nDim) :: PolarisationVectors
    ! Random numbers for fluctuations in mode decomposition
    complex(real64) :: r_afield(npol_transverse), r_efield(npol_transverse)
    ! Momentum
    complex(fp), dimension(nDim) :: Momentum(nDim)
    ! Norm of momentum
    real(fp) :: MomentumNorm

    ! Lattice
    integer(int64) :: MemoryIndex, is, LatticeIndex
    integer(int8) :: i,a,pol

    ! Link
    complex(fp) :: Link(nSUN,nSUN)
    ! E-field
    real(fp) :: Efield_site(ngen)
    ! A-field
    real(fp) :: Afield_site(ngen)
    character(len=105) :: errormessage

    complex(fp) :: afield_p(ndim),efield_p(ndim)

    complex(fp), allocatable :: r(:,:,:)
    
    ! Check if random numbers are initialised
    if(.not.isModuleInitialised_random()) then
       errormessage = 'Error in TransversePolarisedOccupiedInit of '//modulename&
               //':'// modulename_random //' is not initialised.'
       call MPISTOP(errormessage)
    end if
    
    ! ..--** START: Allocations **--..
    allocate(afield(GetMemorySize(),nGen,nDim))
    allocate(efield(GetMemorySize(),nGen,nDim))
    ! ..--**  END : Allocations **--..

    ! Drawing random numbers in x-space ...
    allocate(r(GetMemorySize(),nGen,nPol_transverse))
    
    do MemoryIndex=1,GetMemorysize()
       if(GetProc_M(MemoryIndex)==ThisProc()) then
          do pol=1,nPol_transverse
             r(MemoryIndex,:,pol) = real(GetRandomNormalCmplx(int(ngen,int64)),fp)
          end do
       end if
    end do

    ! ... and FFT into p-space
    do pol=1,nPol_transverse
       do a=1,nGen
          call x2p(r(:,a,pol))
       end do
    end do
    
    ! ..--** START: Setting the field in p-space **--..
    afield=0
    efield=0
    
    do concurrent(&
         MemoryIndex=1:GetMemorySize(),a=1:nGen,i=1:nDim,&
         ThisProc()==GetProc_M(MemoryIndex) .and. GetNorm2Momentum_M(MemoryIndex)>GetZeroTol())

       MomentumNorm = GetNorm2Momentum_M(MemoryIndex)
       Momentum = GetMomentum_M(MemoryIndex)
       call GetPolarisationVectors(Momentum,PolarisationVectors)
       prefactor_afield = &
            sqrt(Occupation(MemoryIndex)/MomentumNorm*2)
       prefactor_efield = &
            cmplx(0,sqrt(Occupation(MemoryIndex)*MomentumNorm*2),fp)

       r_afield = r(MemoryIndex,a,:)
       r_efield = r_afield

       afield(MemoryIndex,a,i) = &
            prefactor_afield*(&
                                ! First transversal polarisation
            + r_afield(1)*PolarisationVectors(i,1) &
                                ! Second transversal polarisation
            + r_afield(2)*PolarisationVectors(i,2) )

       efield(MemoryIndex,a,i) = &
            prefactor_efield*(&
                                ! First transversal polarisation
            + r_efield(1)*PolarisationVectors(i,1) &
                                ! Second transversal polarisation
            + r_efield(2)*PolarisationVectors(i,2) )
    end do
    ! ..--**  END : Drawing random numbers and setting the field in p-space **--..

    if(present(aa_correlator)) then
       allocate(aa_correlator(GetLocalLatticeSize()))
       is=0
       do MemoryIndex=1,GetMemorySize()
          LatticeIndex = GetLatticeIndex_M(MemoryIndex)
          if(GetProc_G(LatticeIndex)==ThisProc()) then
             is = is+1
             aa_correlator(is) = 0
             do concurrent(a=1:ngen)
                afield_p = afield(MemoryIndex,a,:)
                aa_correlator(is) = aa_correlator(is) &
                     + GetTransverseField_usingProjector(afield_p,LatticeIndex)&
                     /GetVolume()/nGen
             end do
          end if
       end do
    end if
    if(present(ee_correlator)) then
       allocate(ee_correlator(GetLocalLatticeSize()))
       is=0
       do MemoryIndex=1,GetMemorySize()
          LatticeIndex = GetLatticeIndex_M(MemoryIndex)
          if(GetProc_G(LatticeIndex)==ThisProc()) then
             is = is+1
             ee_correlator(is) = 0
             do concurrent(a=1:ngen)
                efield_p = efield(MemoryIndex,a,:)
                ee_correlator(is) = ee_correlator(is) &
                     + GetTransverseField_usingProjector(efield_p,LatticeIndex)&
                     /GetVolume()/nGen
             end do
          end if
       end do
    end if
    
    !..--** START: FFT p-->x **--..
    do i=1,ndim
       do a=1,ngen
          call p2x(afield(:,a,i))
          call p2x(efield(:,a,i))
       end do !a
    end do !i
    !..--**  END : FFT p-->x **--..
    
    !..--** START: Writing fields to configuration **--..
    call GaugeConf%Allocate
    
    do concurrent(MemoryIndex=1:GetMemorySize(), i=1:ndim, ThisProc()== GetProc_M(MemoryIndex))
       ! Link
       afield_site = real(afield(MemoryIndex,:,i),fp) &
                                ! Translation to lattice units
            *GetLatticeSpacing(i)
       GaugeConf%Links(:,:,i,MemoryIndex) = GetGroupExp(afield_site)

       ! E-field
       efield_site = real(efield(MemoryIndex,:,i),fp) &
                                ! Translation to lattice units
            *GetLatticeSpacing(i)
       GaugeConf%Efield(:,i,MemoryIndex) = efield_site
    end do
    !..--**  END : Writing fields to configuration **--..
    call GaugeConf%CommunicateBoundary()
  end subroutine TransversePolarisedOccupiedInit

  !>@brief Box occupation
  !!@details Step function with amplitude until saturation scale. Afterwards coupling\f$^2/2\f$
  !!@returns Box occupation
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 24.02.2019
  !!@version 1.0
 pure elemental real(fp) function GetBoxOccupation(Momentum,SaturationScale,Amplitude,Coupling)
    use precision, only: fp
    implicit none
    !> Lattice momentum norm
    real(fp), intent(in) :: Momentum
    !> Saturation scale \f$q_s\f$
    real(fp), intent(in)  :: SaturationScale
    !> Occupation of the box for \f$p<q_s\f$
    real(fp), intent(in)  :: Amplitude
    !> Coupling
    real(fp), intent(in)  :: Coupling
    if( Momentum < SaturationScale ) then
       GetBoxOccupation = Amplitude
    else
       GetBoxOccupation = Coupling**2/2
    end if
  end function GetBoxOccupation
    

  !>@brief Total deviation from Gauss-law of the total MPI-distributed configuration
  !!@returns Total deviation from Gauss-law of the total MPI-distributed configuration
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
  impure real(fp) function GetDeviationFromGaussLaw(GaugeConf)
    use mpiinterface, only: intmpi, GetRealSendType, ThisProc
    use mpi
    use lattice, only: nDim, GetLatticeSpacing,GetNeib_M,GetLatticeIndex_M, GetMemorySize, GetProc_M
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf

    integer(intmpi) :: mpierr
    real(fp) :: local_contribution

    real(fp), dimension(ngen) :: efield, efield_neib
    complex(fp), dimension(nsun,nsun) :: Mefield, Mefield_neib, derivative, Link_neib
    
    integer(int8)  :: i, a
    integer(int64) :: MemoryIndex, neib

    ! 1. Calculation of local contribution
    local_contribution = 0
    
    do concurrent(MemoryIndex=1:GetMemorySize(),i=1:ndim,ThisProc()==GetProc_M(MemoryIndex))
       efield = GaugeConf%GetEfield_M([1_int8:ngen],i,MemoryIndex)

       neib = GetNeib_M(-i,MemoryIndex)
       efield_neib = GaugeConf%GetEfield_M([1_int8:ngen],i,Neib)

       Mefield = GetAlgebraMatrix(efield)
       Mefield_neib = GetAlgebraMatrix(efield_neib)

       Link_neib = GaugeConf%GetLink_M(i,Neib)

       derivative = (Mefield &
            - matmul(matmul(&
            conjg(transpose(Link_Neib)),&
            Mefield_neib),&
            Link_neib))&
            /GetLatticeSpacing(i)**2
       do concurrent(a=1_int8:ngen)
          local_contribution = local_contribution &
               + 2*Abs(Aimag(GetTraceWithGenerator(a,derivative)))
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
 pure function GetFieldStrengthTensor_G(GaugeConf,i,j,LatticeIndex)
    use lattice, only: GetMemoryIndex
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Direction
    integer(int8),             intent(in) :: i,j
    !> Lattice index
    integer(int64),            intent(in) :: LatticeIndex
    !> Field strength tensor
    complex(fp) :: GetFieldStrengthTensor_G(nSUN,nSUN)

    complex(fp) :: Plaquette(nSUN,nSUN)
    
    GetFieldStrengthTensor_G = GaugeConf%GetFieldStrengthTensor_M(i,j,GetMemoryIndex(LatticeIndex))
  end function GetFieldStrengthTensor_G

  !>@brief Returns field strength tensor
  !!@details: Field strength tensor in physical units \f$g F_{ij}(\vec{v})^aT^a\f$
  !!in order \f$O(a_i^2,a_j^2,a_i\cdot a_j)\f$
  !!@returns Spatial field strength tensor \f$F_{i,j,\vec{v}}\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 28.02.2019
  !!@version 1.0
 pure function GetFieldStrengthTensor_M(GaugeConf,i,j,MemoryIndex)
    use lattice, only: GetLatticeSpacing
    use matrixoperations, only: LogU
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Direction
    integer(int8),             intent(in) :: i,j
    !> Memory index
    integer(int64),            intent(in) :: MemoryIndex
    !> Field strength tensor
    complex(fp) :: GetFieldStrengthTensor_M(nSUN,nSUN)

    complex(fp) :: Plaquette(nSUN,nSUN)

    Plaquette = GaugeConf%GetPlaquette_M(i,j,MemoryIndex)
    
    GetFieldStrengthTensor_M = &
         LogU(Plaquette)/GetLatticeSpacing(i)/GetLatticeSpacing(j)/cmplx(0,1,fp)
  end function GetFieldStrengthTensor_M
  
  !> @brief Returns spatial plaquette \f$U_{ij,\vec{v}} =U_{i,\vec{v}}\cdot U_{j,\vec{v}+\hat{i}}
  !! \cdot U_{i,\vec{v}+\hat{j}}^\dagger\cdot U_{j,\vec{v}}^\dagger\f$
  !! @returns Spatial plaquette
  !! \f$U_{i,j,v}=U_{i,v}\cdot U_{j,v+\hat{i}}\cdot U_{i,v+\hat{j}}^\dagger\cdot U_{j,v}^\dagger\f$
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 28.02.2019
  !! @version 1.0
 pure function GetSpatialPlaquette_M(GaugeConf,i,j,MemoryIndex)
    use lattice, only: GetNeib_M
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Spatial direction
    integer(int8),             intent(in) :: i,j
    !> Memory index
    integer(int64),            intent(in) :: MemoryIndex
    !> Spatial plaquette
    complex(fp) :: GetSpatialPlaquette_M(nSUN,nSUN)

    complex(fp), dimension(nSUN,nSUN) :: link1,link2,link3,link4

    link1 = GaugeConf%GetLink_M(i,MemoryIndex)                               ! U_i(x)
    link2 = GaugeConf%GetLink_M(j,GetNeib_M(i,MemoryIndex))                  ! U_j(x+î)
    link3 = conjg(transpose(GaugeConf%GetLink_M(i,GetNeib_M(j,MemoryIndex))))! U_i(x+ĵ)†
    link4 = conjg(transpose(GaugeConf%GetLink_M(j,MemoryIndex)))             ! U_j(x)†

    GetSpatialPlaquette_M = matmul(matmul(matmul(&
         link1,&  ! U_i(x)
         link2),& ! U_j(x+î)
         link3),& ! U_i(x+ĵ)†
         link4)   ! U_j(x)†
  end function GetSpatialPlaquette_M

  !> @brief Returns spatial plaquette \f$U_{ij,\vec{v}} =U_{i,\vec{v}}\cdot U_{j,\vec{v}+\hat{i}}
  !! \cdot U_{i,\vec{v}+\hat{j}}^\dagger\cdot U_{j,\vec{v}}^\dagger\f$
  !! @returns Spatial plaquette
  !! \f$U_{i,j,v}=U_{i,v}\cdot U_{j,v+\hat{i}}\cdot U_{i,v+\hat{j}}^\dagger\cdot U_{j,v}^\dagger\f$
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 28.02.2019
  !! @version 1.0
  pure function GetSpatialPlaquette_G(GaugeConf,i,j,LatticeIndex)
    use lattice, only: GetMemoryIndex
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Spatial direction
    integer(int8),             intent(in) :: i,j
    !> Lattice index
    integer(int64),            intent(in) :: LatticeIndex
    !> Spatial plaquette
    complex(fp) :: GetSpatialPlaquette_G(nSUN,nSUN)
    GetSpatialPlaquette_G=GaugeConf%GetSpatialPlaquette_M(i,j,GetMemoryIndex(LatticeIndex))
  end function GetSpatialPlaquette_G

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
  !! @date 28.02.2019
  !! @version 1.0
 pure function GetTemporalPlaquette_M(GaugeConf,i,j,MemoryIndex)
    use matrixoperations, only: GetUnitMatrix
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Direction
    integer(int8),             intent(in) :: i,j
    !> Memory index
    integer(int64),            intent(in) :: MemoryIndex
    !> Temporal plaquette
    complex(fp) :: GetTemporalPlaquette_M(nSUN,nSUN)

    real(fp) :: efield(ngen)
    
    if(     i==0 .and. j/=0 ) then
       efield = +GaugeConf%GetEfield_M([1_int8:ngen],j,MemoryIndex)
       GetTemporalPlaquette_M = GetGroupExp(efield)
    elseif( i/=0 .and. j==0 ) then
       efield = -GaugeConf%GetEfield_M([1_int8:ngen],i,MemoryIndex)
       GetTemporalPlaquette_M = GetGroupExp(efield)
    else
       GetTemporalPlaquette_M = GetUnitMatrix(nSUN)
    end if
  end function GetTemporalPlaquette_M

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
  !! @date 28.02.2019
  !! @version 1.0
  pure function GetTemporalPlaquette_G(GaugeConf,i,j,LatticeIndex)
    use lattice, only: GetMemoryIndex
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Direction
    integer(int8),             intent(in) :: i,j
    !> Lattice index
    integer(int64),            intent(in) :: LatticeIndex
    !> Temporal plaquette
    complex(fp) :: GetTemporalPlaquette_G(nSUN,nSUN)
    GetTemporalPlaquette_G = GaugeConf%GetTemporalPlaquette_M(i,j,GetMemoryIndex(LatticeIndex))
  end function GetTemporalPlaquette_G
  
  !> @brief Returns plaquette \f$U_{jk,\vec{v}}
  !! =U_{j,\vec{v}}\cdot U_{k,\vec{v}+\hat{j}}
  !! \cdot U_{j,\vec{v}+\hat{k}}^\dagger\cdot U_{k,\vec{v}}^\dagger\f$
  !! @returns Plaquette
  !! \f$U_{j,k,v}=U_{j,v}\cdot U_{k,v+\hat{j}}
  !! \cdot U_{j,v+\hat{k}}^\dagger\cdot U_{k,v}^\dagger\f$
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 28.02.2019
  !! @version 1.0
 pure function GetPlaquette_M(GaugeConf,i,j,MemoryIndex)
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Direction
    integer(int8),             intent(in) :: i,j
    !> Memory index
    integer(int64),            intent(in) :: MemoryIndex
    !> Temporal plaquette
    complex(fp) :: GetPlaquette_M(nSUN,nSUN)
    if(i==0 .or. j==0) then
       GetPlaquette_M = GaugeConf%GetTemporalPlaquette_M(i,j,MemoryIndex)
    else
       GetPlaquette_M = GaugeConf%GetSpatialPlaquette_M(i,j,MemoryIndex)
    end if
  end function GetPlaquette_M

  !> @brief Returns plaquette \f$U_{jk,\vec{v}}
  !! =U_{j,\vec{v}}\cdot U_{k,\vec{v}+\hat{j}}
  !! \cdot U_{j,\vec{v}+\hat{k}}^\dagger\cdot U_{k,\vec{v}}^\dagger\f$
  !! @returns Plaquette
  !! \f$U_{j,k,v}=U_{j,v}\cdot U_{k,v+\hat{j}}
  !! \cdot U_{j,v+\hat{k}}^\dagger\cdot U_{k,v}^\dagger\f$
  !! @author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 28.02.2019
  !! @version 1.0
  pure function GetPlaquette_G(GaugeConf,i,j,LatticeIndex)
    use lattice, only: GetMemoryIndex
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Direction
    integer(int8),             intent(in) :: i,j
    !> Lattice index
    integer(int64),            intent(in) :: LatticeIndex
    !> Temporal plaquette
    complex(fp) :: GetPlaquette_G(nSUN,nSUN)
    GetPlaquette_G = GaugeConf%GetPlaquette_M(i,j,GetMemoryIndex(LatticeIndex))
  end function GetPlaquette_G

  !>@brief Energy functional
  !!@returns Energy functional
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.02.2019
  !!@version 1.0
  impure real(fp) function GetEnergy(GaugeConf)
    use precision, only: fp
    use matrixoperations, only: GetTrace
    use lattice, only: nDim, GetLatticeSpacing, GetMemorySize, GetProc_M, GetLatticeIndex_M
    use mpiinterface, only: intmpi, GetRealSendType, ThisProc
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
    integer(int64) :: MemoryIndex

    ! 1. Calculation of local contribution
    local_contribution = 0
    do concurrent(MemoryIndex=1:GetMemorySize(),i=1:ndim,ThisProc()==GetProc_M(MemoryIndex))
       
       efield = GaugeConf%GetElectricField_AlgebraCoordinate([1_int8:ngen],i,GetLatticeIndex_M(MemoryIndex))

       local_contribution = local_contribution &
                                ! Electric energy
            + sum(efield**2)/2
       ! Magnetic energy
       do concurrent(j=i+1_int8:ndim)
          Plaquette = GaugeConf%GetPlaquette_M(i,j,MemoryIndex)
          PotentialTerm = (1-real(GetTrace(Plaquette),fp)/nSUN)&
               /GetLatticeSpacing(i)/GetLatticeSpacing(j)

          Local_contribution = Local_Contribution &
               + 2*nSUN*PotentialTerm
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

  !>@brief Computes the transverse AA-correlator (MPI-distributed)
  !!@returns Transverse AA-correlator (MPI-distributed)
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 25.02.2019
  !!@version 1.0
  impure subroutine GetTransverseAACorrelator(GaugeConf,Correlator)
    use, intrinsic :: iso_fortran_env
    use precision, only: fp
    use lattice, only: nDim, GetLatticeIndex_M, GetMemoryIndex, GetProc_G, GetMemorySize,&
         GetLocalLatticeSize
    use mpiinterface, only: ThisProc
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in)  :: GaugeConf
    !> Transverse AA-correlator (MPI-distributed)
    real(fp), allocatable, intent(out) :: Correlator(:)

    integer(int8)  :: a, i
    integer(int64) :: MemoryIndex, LatticeIndex

    real(fp), allocatable :: field(:,:,:), correlator_a(:)
    
    ! Extracting the gauge field in position space
    allocate(field(GetMemorySize(),nDim,nGen))
    forall(MemoryIndex=1:GetMemorySize(),i=1:ndim)
       field(MemoryIndex,i,:)&
            = GaugeConf%GetGaugeField_AlgebraCoordinates(&
            i,GetLatticeIndex_M(MemoryIndex))
    end forall

    ! Computing average over all generators
    allocate(correlator(GetLocalLatticeSize()))
    correlator = 0
    
    do a=1,nGen
       call GetTransverseCorrelator(field(:,:,a),correlator_a)
       correlator = correlator + correlator_a/nGen
    end do
    deallocate(correlator_a,field)
  end subroutine GetTransverseAACorrelator

  !>@brief Computes the transverse EE-correlator (MPI-distributed)
  !!@returns Transverse EE-correlator (MPI-distributed)
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 25.02.2019
  !!@version 1.0
  impure subroutine GetTransverseEECorrelator(GaugeConf,Correlator)
    use, intrinsic :: iso_fortran_env
    use precision, only: fp
    use lattice, only: nDim, GetLatticeIndex_M, GetMemoryIndex, GetProc_G, GetMemorySize,&
         GetLocalLatticeSize
    use mpiinterface, only: ThisProc
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in)  :: GaugeConf
    !> Transverse AA-correlator (MPI-distributed)
    real(fp), allocatable, intent(out) :: Correlator(:)

    integer(int8)  :: a, i
    integer(int64) :: MemoryIndex, LatticeIndex

    real(fp), allocatable :: field(:,:,:), correlator_a(:)
    
    ! Extracting the gauge field in position space
    allocate(field(GetMemorySize(),nDim,nGen))
    forall(MemoryIndex=1:GetMemorySize(),i=1:ndim)
       field(MemoryIndex,i,:)&
            = GaugeConf%GetElectricField_AlgebraCoordinate(&
            [1_int8:ngen],i,GetLatticeIndex_M(MemoryIndex))
    end forall

    ! Computing average over all generators
    allocate(correlator(GetLocalLatticeSize()))
    correlator = 0
    
    do a=1,nGen
       call GetTransverseCorrelator(field(:,:,a),correlator_a)
       correlator = correlator + correlator_a/nGen
    end do
    deallocate(correlator_a,field)
  end subroutine GetTransverseEECorrelator

  !>@brief Computes the transverse correlator (MPI-distributed)
  !!@returns Transverse correlator
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 17.01.2019
  !!@version 1.0
  impure subroutine GetTransverseCorrelator(x_field,correlator)
    use lattice, only: nDim, GetVolume, GetLatticeIndex_M, GetMemorySize, GetLocalLatticeSize, GetProc_G
    use xpfft,     only: x2p
    use precision, only: fp
    use mpiinterface, only: ThisProc
    implicit none
    !> Field in position space (e.g. A(x), E(x)), indexed as: field(latticeindex,direction)
    real(fp),              intent(in)  :: x_field(:,:)
    !> Transverse correlator (MPI-distributed)
    real(fp), allocatable, intent(out) :: correlator(:)


    integer(int8) :: i
    integer(int64) :: LatticeIndex, MemoryIndex, is

    complex(fp), allocatable :: p_field(:,:)
    complex(fp) :: field_site(nDim)

    ! Allocating fourier transform of input
    allocate(p_field(GetMemorySize(),nDim))
    p_field = cmplx(x_field,0,fp)

    ! Performing fourier transform for all spatial components
    do i=1,ndim
       call x2p(p_field(:,i))
    end do
    
    ! Allocating output
    allocate(correlator(GetLocalLatticeSize()))
    correlator = 0
    
    is=0
    do MemoryIndex=1,GetMemorySize()
       LatticeIndex = GetLatticeIndex_M(MemoryIndex)
       if(GetProc_G(LatticeIndex)==ThisProc()) then
          is=is+1
          field_site = p_field(MemoryIndex,:)
          Correlator(is) = &
               GetTransverseField_usingProjector&
               (field_site,LatticeIndex)/GetVolume()
       end if
    end do
    deallocate(p_field)
  end subroutine GetTransverseCorrelator

  pure real(fp) function GetTransverseField_usingProjector(field,LatticeIndex)
    use precision, only: fp
    use lattice, only: nDim, GetMomentum_G, GetTransverseProjector
    implicit none
    !> Field in p-space
    complex(fp),    intent(in) :: field(nDim)
    !> Lattice index
    integer(int64), intent(in) :: LatticeIndex

    integer(int8) :: i,j
    complex(fp) :: res, momentum(nDim)

    res = 0
    do concurrent(i=1:ndim, j=1:ndim)
       momentum = GetMomentum_G(LatticeIndex)
       res = res + field(i)*GetTransverseProjector(i,j,momentum)*conjg(field(j))
    end do
    GetTransverseField_usingProjector = real(res,fp)
  end function GetTransverseField_usingProjector

  !>@brief Returns algebra-component of the gauge field
  !!@returns \f$A_{i,\vec{v}}^a\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 25.02.2019
  !!@version 1.0
 pure elemental function GetGaugefield_AlgebraCoordinate(GaugeConf,a,i,LatticeIndex)
    use precision, only: fp
    use lattice, only: GetLatticeSpacing
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Generator index
    integer(int8),  intent(in) :: a
    !> Direction
    integer(int8),  intent(in) :: i
    !> Latticeindex
    integer(int64), intent(in) :: LatticeIndex
    !> Algebra coordinates of electric field
    real(fp) :: GetGaugefield_AlgebraCoordinate

    real(fp) :: algebra_coordinates(ngen)

    algebra_coordinates = GaugeConf%GetGaugefield_AlgebraCoordinates(i,LatticeIndex)
    GetGaugefield_AlgebraCoordinate = algebra_coordinates(a)/GetLatticeSpacing(i)
  end function GetGaugefield_AlgebraCoordinate

  !>@brief Returns all algebra-components of the gauge field
  !!@returns  All \f$A_{i,\vec{v}}^a\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 25.02.2019
  !!@version 1.0
 pure function GetGaugefield_AlgebraCoordinates(GaugeConf,i,latticeindex)
    use lattice, only: GetLatticeSpacing
    implicit none
    !> Gauge link configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Direction
    integer(int8),  intent(in) :: i
    !> Latticeindex
    integer(int64), intent(in) :: latticeindex
    !> Algebra coordinates of electric field
    real(fp) :: GetGaugefield_AlgebraCoordinates(ngen)

    complex(fp) :: Link(nSUN,nSUN)

    Link = GaugeConf%GetLink_G(i,LatticeIndex)

    GetGaugefield_AlgebraCoordinates = GetGroupLog(Link)/GetLatticeSpacing(i)
  end function GetGaugefield_AlgebraCoordinates

  !>@brief Returns algebra-component of electric field
  !!@returns \f$E_{i,\vec{v}}^a\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 25.02.2019
  !!@version 1.0
 pure elemental function GetElectricField_AlgebraCoordinate(conf,a,i,latticeindex)
    use precision, only: fp
    use lattice, only: GetLatticeSpacing
    implicit none
    !> Gauge link configuration
    class(GaugeConfiguration), intent(in) :: conf
    !> Generator index
    integer(int8),  intent(in) :: a
    !> Direction
    integer(int8),  intent(in) :: i
    !> Latticeindex
    integer(int64), intent(in) :: latticeindex
    !> Algebra coordinates of electric field
    real(fp) :: GetElectricField_AlgebraCoordinate
    
    GetElectricField_AlgebraCoordinate &
         = conf%GetEfield_G(a,i,latticeindex)/GetLatticeSpacing(i)
  end function GetElectricField_AlgebraCoordinate

  !>@brief Returns electric field
  !!@returns \f$E_{i,\vec{v}}\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 25.02.2019
  !!@version 1.0
 pure function GetElectricField_AlgebraMatrix(conf,i,latticeindex)
    use precision, only: fp
    use lattice, only: GetLatticeSpacing
    implicit none
    !> Gauge link configuration
    class(GaugeConfiguration), intent(in) :: conf
    !> Direction
    integer(int8),  intent(in) :: i
    !> Latticeindex
    integer(int64), intent(in) :: latticeindex
    !> Electric field
    complex(fp) :: GetElectricField_AlgebraMatrix(nSUN,nSUN)
    
    real(fp) :: efield(ngen)
    efield = conf%GetEfield_G([1_int8:ngen],i,LatticeIndex)
    
    GetElectricField_AlgebraMatrix &
         = GetAlgebraMatrix(efield)/GetLatticeSpacing(i)
  end function GetElectricField_AlgebraMatrix

  !>@brief Update routine using the Leapfrog algorithm
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 27.02.2019
  !!@version 1.0
  impure subroutine Update_Leapfrog(GaugeConf,StepWidth_)
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(inout) :: GaugeConf
    !> Step width in units of \f$a_0\f$
    real(fp),                  intent(in), optional :: Stepwidth_

    real(fp) :: Stepwidth
    
    ! ..--** START: Optional parameters **--..
    if(present(Stepwidth_)) then
       Stepwidth = Stepwidth_
    else
       ! Default value
       Stepwidth = 1._fp ! in units of a0
    end if
    ! ..--**  END : Optional parameters **--..
    
    if(Stepwidth.gt.0) then
       call GaugeConf%Update_Links_Leapfrog(Stepwidth)
       call GaugeConf%CommunicateBoundary_Links
       call GaugeConf%Update_Efield_Leapfrog(Stepwidth)
       call GaugeConf%CommunicateBoundary_Efield
    else
       call GaugeConf%Update_Efield_Leapfrog(Stepwidth)
       call GaugeConf%Update_Links_Leapfrog(Stepwidth)
       call GaugeConf%CommunicateBoundary
       !call GaugeConf%CommunicateBoundary_Efield
    end if
  end subroutine Update_Leapfrog

  !>@brief Update routine for the links as part of the Leapfrog algorithm
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 27.02.2019
  !!@version 1.0
  pure subroutine Update_Links_Leapfrog(GaugeConf,StepWidth)
    use lattice, only: ndim, GetMemorySize, GetProc_M, getlatticespacing
    use mpiinterface, only: ThisProc
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(inout) :: GaugeConf
    !> Step width in units of \f$a_0\f$
    real(fp),                  intent(in)    :: StepWidth
    
    integer(int8) :: i
    integer(int64):: MemoryIndex

    complex(fp) :: TimeEvolutionOperator(nSUN,nSUN)
    real(fp)    :: efield_times_dt(nGen)

    do concurrent(MemoryIndex=1:GetMemorySize(),i=1:ndim,ThisProc()==GetProc_M(MemoryIndex))
       efield_times_dt       = GaugeConf%Efield(:,i,MemoryIndex)*StepWidth*GetLatticeSpacing(0_int8)
       TimeEvolutionOperator = GetGroupExp(efield_times_dt)

       GaugeConf%Links(:,:,i,MemoryIndex) =&
            matmul(TimeEvolutionOperator,GaugeConf%Links(:,:,i,MemoryIndex))
    end do
  end subroutine Update_Links_Leapfrog

  !>@brief Update routine for the electric field as part of the Leapfrog algorithm
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 27.02.2019
  !!@version 1.0
  pure subroutine Update_Efield_Leapfrog(GaugeConf,StepWidth)
    use mpiinterface, only: ThisProc
    use lattice, only: ndim, GetMemorySize, GetProc_M
    implicit none
    !> Gauge link configuration
    class(GaugeConfiguration), intent(inout) :: GaugeConf
    !> Step width in units of \f$a_0\f$
    real(fp),                  intent(in)    :: StepWidth
    
    integer(int8) :: i
    integer(int64):: MemoryIndex
    
    do concurrent(MemoryIndex=1:GetMemorySize(),i=1:ndim,ThisProc()==GetProc_M(MemoryIndex))
       call Update_Efield_Leapfrog_atSite_Direction(GaugeConf,StepWidth,i,MemoryIndex)
    end do

  contains
    pure subroutine Update_Efield_Leapfrog_atSite_Direction(GaugeConf,StepWidth,i,MemoryIndex)
      use lattice, only: ndim, GetLatticeSpacing, GetMemoryIndex
      implicit none
      !> Gauge link configuration
      type(GaugeConfiguration), intent(inout) :: GaugeConf
      !> Step width in units of \f$a_0\f$
      real(fp),                 intent(in)    :: StepWidth
      !> Direction
      integer(int8),            intent(in)    :: i
      !> Memory index
      integer(int64),           intent(in)    :: MemoryIndex

      integer(int8) :: k,a
      complex(fp) :: staplesum(nsun,nsun), link_times_staplesum(nsun,nsun)
      
      staplesum = 0
      do concurrent(k=1:ndim, k/=i)
         staplesum = staplesum + &
              (&
              GetUStaple(GaugeConf,i,k,MemoryIndex) + &
              GetDStaple(GaugeConf,i,k,MemoryIndex) &
              )/GetLatticeSpacing(k)**2
      end do !k

      link_times_staplesum = matmul(GaugeConf%Links(:,:,i,MemoryIndex),staplesum)

      forall(a=1:ngen)
         GaugeConf%efield(a,i,MemoryIndex)&
           = GaugeConf%efield(a,i,MemoryIndex)&
           - 2*StepWidth*GetLatticeSpacing(0)*Aimag(GetTraceWithGenerator(a,link_times_staplesum))
      end forall
      !GaugeConf%efield(:,i,MemoryIndex) = &
      !     GetWrappedAlgebraCoordinates(GaugeConf%efield(:,i,MemoryIndex))
    end subroutine Update_Efield_Leapfrog_AtSite_Direction

    !>@brief Returns staple
    !!@details The U-staple is defined as
    !! \f$ S^{\text{U}}_{ik,x}=U_{k,x+\hat{i}}\cdot U_{i,x+\hat{k}}^\dagger\cdot U_{k,x}^\dagger\f$
    !!@autho Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 27.02.2019
    !!@version 1.0
    pure function GetUStaple(GaugeConf,i,k,MemoryIndex)
      use lattice, only: GetNeib_M
      implicit none
      !> Gauge configuration
      type(GaugeConfiguration), intent(in) :: GaugeConf
      !> Direction
      integer(int8),            intent(in) :: i,k
      !> Memory index
      integer(int64),           intent(in) :: MemoryIndex
      !> Upwards staple
      complex(fp)                          :: GetUStaple(nsun,nsun)

      integer(int64) :: neib_i, neib_k

      neib_i = GetNeib_M(+i,MemoryIndex)
      neib_k = GetNeib_M(+k,MemoryIndex)

      GetUStaple = matmul(matmul(&
           GaugeConf%links(                :,:,k,neib_i),&
           conjg(transpose(GaugeConf%links(:,:,i,neib_k)))),&
           conjg(transpose(GaugeConf%links(:,:,k,MemoryIndex))))
    end function GetUStaple
    
    !>@brief Returns staple
    !!@details The U-staple is defined as \f$ S^{\text{D}}_{ik,x}=
    !! U_{k,x+\hat{i}-\hat{k}}^\dagger\cdot U_{i,x-\hat{k}}^\dagger\cdot U_{k,x-\hat{k}}\f$
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 27.02.2019
    !!@version 1.0
    pure function GetDStaple(GaugeConf,i,k,MemoryIndex)
      use lattice, only: GetNeib_M
      implicit none
      !> Gauge configuration
      type(GaugeConfiguration), intent(in) :: GaugeConf
      !> Direction
      integer(int8),            intent(in) :: i,k
      !> Memory index
      integer(int64),           intent(in) :: MemoryIndex
      !> Downwards staple
      complex(fp)                          :: GetDStaple(nsun,nsun)

      integer(int64) :: neib_k, neib_ik

      neib_k  = GetNeib_M(-k,MemoryIndex)
      neib_ik = GetNeib_M(+i,neib_k)

      GetDStaple = matmul(matmul(&
           conjg(transpose(GaugeConf%links(:,:,k,neib_ik))),&
           conjg(transpose(GaugeConf%links(:,:,i,neib_k)))),&
           GaugeConf%links(:,:,k,neib_k))
    end function GetDStaple
  end subroutine Update_Efield_Leapfrog


  !>@brief Performs gauge-fixing to Coulomb gauge
  !!@details Exactly Davies' algorithm. See <a href="https://doi.org/10.1103/PhysRevD.37.1581">
  !! Fourier acceleration in lattice gauge theories. I. Landau gauge fixing</a>
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 28.02.2019
  !!@version 1.0
  impure subroutine CoulombGaugefixing(GaugeConf,Tolerance_,alpha_,MaxIterations_)
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(inout) :: GaugeConf
    !> Tolerance for gauge fixing
    real(fp),                  intent(in), optional :: Tolerance_
    !> Alpha (slope) in Coulomb gauge fixing
    real(fp),                  intent(in), optional :: alpha_
    !> Maximum number of gauge fixing iteration steps
    integer(int64),            intent(in), optional :: MaxIterations_

    real(fp) :: Tolerance, alpha
    integer(int64) :: MaxIterations

    type(GaugeConfiguration) :: UpdatedConf

    ! ..--** START: Optional parameters **--..
    if(present(Tolerance_)) then
       Tolerance = Tolerance_
    else
       ! Default value
       Tolerance = 1.E-8
    end if

    if(present(alpha_)) then
       alpha = alpha_
    else
       ! Default value from
       ! Davies et al. in https://doi.org/10.1103/PhysRevD.37.1581
       alpha = 0.08_fp
    end if

    if(present(MaxIterations_)) then
       MaxIterations = MaxIterations_
    else
       ! Default value
       MaxIterations = 100
    end if
    ! ..--**  END : Optional parameters **--..

    call GaugeConf%CoulombGaugefixing_Links(Tolerance,Alpha,MaxIterations)
  end subroutine CoulombGaugefixing

  !>@brief Performs gauge-fixing to Coulomb gauge for the link variables
  !!@details Exactly Davies' algorithm. See <a href="https://doi.org/10.1103/PhysRevD.37.1581">
  !! Fourier acceleration in lattice gauge theories. I. Landau gauge fixing</a>
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 28.02.2019
  !!@version 1.0
  impure subroutine CoulombGaugefixing_Links(GaugeConf,Tolerance,Alpha,MaxIterations)
    use lattice, only: GetMaxNorm2Momentum, GetNorm2Momentum_M, nDim, GetNeib_M,&
         GetProc_M, GetMemorySize
    use xpfft, only: x2p, p2x
    use tolerances, only: GetZeroTol
    use mpiinterface, only: ThisProc
    use halocomm, only: MPI_CommunicateBoundary => CommunicateBoundary
    use matrixoperations, only: ExpAH
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(inout) :: GaugeConf
    !> Tolerance for gauge fixing
    real(fp),                  intent(in)    :: Tolerance
    !> Alpha (slope) in Coulomb gauge fixing
    real(fp),                  intent(in)    :: alpha
    !> Maximum number of gauge fixing iteration steps
    integer(int64),            intent(in)    :: MaxIterations

    ! Divergence term
    complex(fp), allocatable :: Divergence(:,:,:)
    ! Gauge transformation
    complex(fp), allocatable :: GaugeTransform(:,:,:)
    ! FFT-Inout
    complex(fp), allocatable :: fft_inout(:)

    ! Indices
    integer(int8) :: i,j
    integer(int64) :: MemoryIndex, iteration

    ! Acceleration step
    real(fp) :: p, pmax

    !..--** START: Allocation **--..
    allocate(Divergence(nsun,nsun,GetMemorySize()))
    allocate(GaugeTransform(nsun,nsun,GetMemorySize()))
    allocate(fft_inout(GetMemorySize()))
    !..--**  END : Allocation **--..

    GaugeFixingIteration: do iteration=1,MaxIterations

       ! Computing divergence
       Divergence=0
       forall(MemoryIndex=1:GetMemorySize(),ThisProc()==GetProc_M(MemoryIndex))
          Divergence(:,:,MemoryIndex) = GaugeConf%GetDivergenceOfGaugeField_M(MemoryIndex)
       end forall
       call MPI_CommunicateBoundary(Divergence)

       ! Checking if Coulomb gauge condition is fulfilled
       if(isCoulombGauged(Divergence,Tolerance)) exit GaugeFixingIteration

       ! ..--** START: Constructing the gauge transformation **--..
       
       ! ..--** START: Fourier Acceleration **--..
       do j=1,nsun
          do i=1,nsun
             fft_inout = Divergence(i,j,:)
             call x2p(fft_inout)
             Divergence(i,j,:) = fft_inout
          end do
       end do

       ! Acceleration step
       pmax = GetMaxNorm2Momentum()
       do concurrent(MemoryIndex=1:GetMemorySize(),      & !Indices
            GetNorm2Momentum_M(MemoryIndex)>GetZeroTol() & !|p|>0
            .and. ThisProc()==GetProc_M(MemoryIndex))      !only this process' lattice indices
          Divergence(:,:,MemoryIndex) = &
               Divergence(:,:,MemoryIndex)&
               * (pmax/GetNorm2Momentum_M(MemoryIndex))**2
       end do
       
       ! FFT p-->x
       do j=1,nsun
          do i=1,nsun
             fft_inout = Divergence(i,j,:)
             call p2x(fft_inout)
             Divergence(i,j,:) = fft_inout
          end do
       end do
       ! ..--**  END : Fourier Acceleration **--..

       Divergence = alpha*Divergence

       ! Algebra --> Group
       do concurrent(MemoryIndex=1:GetMemorySize(),ThisProc()==GetProc_M(MemoryIndex))
          gaugetransform(:,:,MemoryIndex) = ExpAH(Divergence(:,:,MemoryIndex))
       end do
       ! ..--**  END : Constructing the gauge transformation **--..

       ! Communicating boundary values of gauge transform
       call MPI_CommunicateBoundary(GaugeTransform)

       ! Performing gauge transformation on the lnks
       call PerformGaugeTransformation(GaugeConf,GaugeTransform)
    end do GaugeFixingIteration
    if(ThisProc()==0) write(output_unit,*)&
         iteration,'iterations in Coulomb gauge fixing performed';&
         call flush(output_unit)
  contains
    impure logical function IsCoulombGauged(Divergence,Tolerance)
      use mpiinterface, only: intmpi, ThisProc, NumProcs, SyncAll, GetRealSendType
      use lattice, only: GetProc_M
      use matrixoperations, only: FrobeniusNorm
      use mpi
      implicit none
      complex(fp), intent(in) :: Divergence(:,:,:)
      real(fp),    intent(in) :: Tolerance

      real(fp), allocatable :: Deviation(:)

      integer(int64) :: MemoryIndex
      logical :: IsCoulombGauged_local
      logical, allocatable :: IsCoulombGauged_allProcs(:)
      integer(intmpi) :: mpierr

      real(fp), allocatable :: maxdeviation(:)

      allocate(deviation(size(Divergence,rank(Divergence))))
      
      deviation = -1 ! Default value, important for mask in local check
      forall(MemoryIndex=1:size(deviation),ThisProc()==GetProc_M(MemoryIndex))
         Deviation(MemoryIndex) = FrobeniusNorm(Divergence(:,:,MemoryIndex))/nsun**2
      end forall

      ! Local check
      IsCoulombGauged_local = .true.
      do concurrent(MemoryIndex=1:size(deviation),&
           ThisProc()==GetProc_M(MemoryIndex).and.deviation(MemoryIndex)>Tolerance)
         IsCoulombGauged_local = .false.
      end do
      
      allocate(maxdeviation(NumProcs()))
      call MPI_ALLGATHER(maxval(deviation),1_intmpi,GetRealSendType(),maxdeviation,1_intmpi,GetRealSendType(),MPI_COMM_WORLD,mpierr)
      if(thisproc()==0) write(output_unit,*) maxval(maxdeviation)
      deallocate(deviation)

      ! Gathering information about local check from the other processes
      allocate(IsCoulombGauged_allProcs(NumProcs()))
      call MPI_ALLGATHER(&
           IsCoulombGauged_local,&
           1_intmpi, MPI_LOGICAL,&
           IsCoulombGauged_allProcs,&
           1_intmpi, MPI_LOGICAL,&
           MPI_COMM_WORLD,mpierr)
      
      ! Coulomb gauge condition fulfilled on all local lattices?
      if(all(IsCoulombGauged_allProcs)) then
         ! yes :)
         IsCoulombGauged = .true.
      else
         ! no  :(
         IsCoulombGauged = .false.
      end if
      deallocate(IsCoulombGauged_allProcs)
    end function IsCoulombGauged

    impure subroutine PerformGaugeTransformation(GaugeConf,GaugeTransformation)
      use lattice, only: GetMemorySize, nDim, GetProc_M, GetNeib_M
      use mpiinterface, only: ThisProc
      implicit none
      type(GaugeConfiguration), intent(inout) :: GaugeConf
      complex(fp),              intent(in)    :: GaugeTransformation(:,:,:)

      integer(int64) :: MemoryIndex
      integer(int8)  :: i

      do concurrent(MemoryIndex=1:GetMemorySize(), i=1:nDim, ThisProc()==GetProc_M(MemoryIndex))
         GaugeConf%Links(:,:,i,MemoryIndex) = matmul(matmul(&
              GaugeTransformation(:,:,MemoryIndex),&
              GaugeConf%Links(:,:,i,MemoryIndex)),&
              conjg(transpose(GaugeTransformation(:,:,GetNeib_M(+i,MemoryIndex)))))
      end do

      call GaugeConf%CommunicateBoundary_Links
    end subroutine PerformGaugeTransformation
  end subroutine CoulombGaugefixing_Links

  pure function GetDivergenceOfGaugeField_M(GaugeConf,MemoryIndex)
    use lattice, only: nDim, GetNeib_M
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Memory index
    integer(int64),            intent(in) :: MemoryIndex
    !> Laplace operator on links
    complex(fp) :: GetDivergenceOfGaugefield_M(nSUN,nSUN)

    integer(int8) :: i
    complex(fp) :: diff(nSUN,nSUN)

    GetDivergenceOfGaugeField_M = 0

    do concurrent(i=1:ndim)
       diff = &
            - GaugeConf%Links(:,:,i,MemoryIndex) &
            + GaugeConf%Links(:,:,i,GetNeib_M(-i,MemoryIndex))

       ! Skew-hermitisation
       diff = (diff - conjg(transpose(diff)))/2

       GetDivergenceOfGaugefield_M = GetDivergenceOfGaugefield_M + diff
    end do
  end function GetDivergenceOfGaugeField_M
end module gaugeconfiguration_su3
