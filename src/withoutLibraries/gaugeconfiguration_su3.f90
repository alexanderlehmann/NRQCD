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

     ! Gauss law deviation
     procedure, public :: GetDeviationFromGausslaw

     ! Energy
     procedure, public :: GetEnergy

     ! Field strength tensor and plaquettes
     procedure, public :: GetFieldStrengthTensorAlgebraCoordinate_G
     procedure, public :: GetFieldStrengthTensor_G
     procedure, public :: GetPlaquette_G
     procedure, private :: GetSpatialPlaquette_G
     procedure, private :: GetTemporalPlaquette_G

     procedure, public :: GetFieldStrengthTensor_M
     procedure, public :: GetPlaquette_M
     procedure, private :: GetSpatialPlaquette_M
     procedure, private :: GetTemporalPlaquette_M

     ! Clover leaf
     procedure, public :: GetFieldStrengthTensor_CloverLeaf_G
     procedure, public :: GetFieldStrengthTensor_CloverLeaf_M

     ! Return of semi-physical fields
     procedure, public :: GetGaugefield_AlgebraCoordinate
     procedure, public :: GetGaugefield_AlgebraCoordinates
     generic, public :: GetElectricField => &
          GetElectricField_AlgebraCoordinate, GetElectricField_AlgebraMatrix
     procedure, private :: GetElectricField_AlgebraMatrix
     procedure, private :: GetElectricField_AlgebraCoordinate

     ! Update routines
     generic, public :: Update => Update_Leapfrog
     procedure, private :: Update_Leapfrog
     procedure, private :: Update_Efield_Leapfrog
     procedure, private :: Update_Links_Leapfrog

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
    use lattice, only: GetNeib_M
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Direction
    integer(int8),  intent(in) :: i
    !> Memory index
    integer(int64), intent(in) :: MemoryIndex
    !> Link variable
    complex(fp) :: GetLink_M(nSUN,nSUN)

    if(i==0) then
       GetLink_M = GetUnitmatrix(nSUN)
    else if(i>0) then
       GetLink_M = GaugeConf%Links(:,:,i,MemoryIndex)
    else if(i<0) then
       GetLink_M = conjg(transpose(GaugeConf%Links(:,:,-i,GetNeib_M(i,MemoryIndex))))
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

  impure subroutine EquilibriumInit(GaugeConf,Beta,nefieldinit,nequilibrium,tolerance_Eprojection,ChargeDensity,MeasureEnergy,filename)
    use lattice!, only: GetLatticeSpacing, GetMemorySize, GetProc_M,&
    !ndim, GetNeib_M, GetLatticeSize
    use random, only: GetRandomNormalCmplx
    use mpiinterface, only: ThisProc, mpistop, intmpi, GetRealSendType
    use mpi
    use io
    use matrixoperations

    ! use gaugeobservables
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(out) :: GaugeConf
    !> \f$\beta_{\text{L}}\f$
    real(fp), intent(in) :: beta
    integer(int64), intent(in) :: nefieldinit
    integer(int64), intent(in) :: nequilibrium

    real(fp), optional, intent(in) :: tolerance_Eprojection
    real(fp), optional, intent(in) :: ChargeDensity(:,:)
    logical, optional, intent(in) :: MeasureEnergy
    character(len=*), optional, intent(in) :: filename

    real(fp) :: sigma

    integer(int8) :: k, a
    integer(int64) :: LatticeIndex, MemoryIndex, x,y,z
    complex(fp) :: r(ngen)

    integer :: iefieldinit
    integer :: iequibstep


    real(fp) :: energy, deviation
    real(fp), parameter :: kappa=0.01_fp!0.12_fp
    real(fp) :: kappa_times_dt
    real(fp), dimension(ngen) :: ec
    complex(fp), dimension(nsun,nsun) :: U, G, Gneib, E, term, uconjg


    integer(int8) :: fileid

    integer(int8) :: i,j
    integer(int64) :: iprojection

    character(len=80) :: energydensity_filename='energydensity.txt'
    integer(int8) :: energydensity_fileid

    integer(intmpi) :: proc, mpierr,mpistatus(mpi_status_size)

    real(fp) :: tolerance_Eprojection_

    type(GaugeConfiguration) :: GaugeConf_previous


    if(present(tolerance_Eprojection)) then
       tolerance_Eprojection_ = tolerance_Eprojection
    else
       ! Default tolerance for projection of E-fields
       tolerance_Eprojection_ = 1E-3
    end if

    sigma = sqrt(1/(GetLatticeSpacing(1_int8)*GetLatticeSpacing(2_int8)*GetLatticeSpacing(3_int8)*beta))
    !sigma = sqrt(1/beta)

    kappa_times_dt = kappa*GetLatticeSpacing(0_int8)

    call GaugeConf%ColdInit

    ! Start thermalizing gauge links
    if(present(MeasureEnergy)) then
       if(MeasureEnergy) then
          if(.not.present(filename)) then
             call MPIStop('Filename for energy output missing')
          end if
          if(thisproc()==0) &
               fileID = OpenFile(filename=filename,&
               st='REPLACE',fm='FORMATTED',act='WRITE')
       end if
    end if

    do iefieldinit=1,nefieldinit
       ! Redrawing E-field
       gaugeconf%efield = 0
       do MemoryIndex=1,GetMemorySize()
          if(ThisProc()==GetProc_M(MemoryIndex)) then
             do k=1,nDim
                r = sigma*GetRandomNormalCmplx(int(ngen,int64))

                GaugeConf%Efield(:,k,MemoryIndex) = real(r,fp)*GetLatticeSpacing(0_int8)
             end do
          end if
       end do

       call GaugeConf%CommunicateBoundary()
       if(present(ChargeDensity)) then
          deviation = GetGdeviation(GaugeConf,ChargeDensity)
       else
          deviation = GetGdeviation(GaugeConf)
       end if
       
       ! Projection of the electric field
       iprojection = 0
       do while(deviation>tolerance_Eprojection_)
          iprojection = iprojection + 1
          GaugeConf_previous = GaugeConf

          do MemoryIndex=1,GetMemorySize()
             if(ThisProc()==GetProc_M(MemoryIndex)) then
                if(present(ChargeDensity)) then
                   G = GetG(GaugeConf_previous,MemoryIndex,ChargeDensity)
                else
                   G = GetG(GaugeConf_previous,MemoryIndex)
                end if
                do k=1,nDim
                   ec = GaugeConf_previous%GetEfield_M([1_int8:ngen],k,MemoryIndex)
                   U = GaugeConf_previous%Links(:,:,k,MemoryIndex)
                   E = GetAlgebraMatrix(ec)

                   if(present(ChargeDensity)) then
                      Gneib = GetG(GaugeConf_previous,GetNeib_M(+k,MemoryIndex),ChargeDensity)
                   else
                      Gneib = GetG(GaugeConf_previous,GetNeib_M(+k,MemoryIndex))
                   end if

                   term = &
                        kappa_times_dt*(matmul(matmul(matmul(&
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

          if(present(ChargeDensity)) then
             deviation = GetGdeviation(GaugeConf,ChargeDensity)
          else
             deviation = GetGdeviation(GaugeConf)
          end if
       end do
       
       ! Evolution of the gauge links
       do iequibstep=1,nequilibrium
          call GaugeConf%Update
       end do

       ! Print energy to terminal
       if(present(MeasureEnergy)) then
          if(MeasureEnergy) then
             energy = GaugeConf%GetEnergy()
             if(ThisProc()==0) then
                write(fileid,*) iefieldinit,energy
                write(output_unit,*) iefieldinit,energy
             end if
          end if
       end if
       if(present(ChargeDensity)) then
          deviation = GetGdeviation(GaugeConf,ChargeDensity)
       else
          deviation = GetGdeviation(GaugeConf)
       end if

       if(ThisProc()==0) print*,'G[U,E]=',deviation
    end do

    if(present(MeasureEnergy)) then
       if(MeasureEnergy) then
          if(ThisProc()==0) then
             call Closefile(fileid)
          end if
       end if
    end if

  contains
    impure real(fp) function GetGdeviation(GaugeConf,ChargeDensity)
      use mpiinterface, only: intmpi, GetRealSendType, ThisProc
      use mpi
      use lattice, only: nDim, GetLatticeSpacing,GetNeib_M,GetLatticeIndex_M, GetMemorySize, GetProc_M
      implicit none
      !> Gauge configuration
      type(GaugeConfiguration), intent(in) :: GaugeConf

      real(fp), optional, intent(in) :: ChargeDensity(:,:)

      integer(intmpi) :: mpierr
      real(fp) :: local_contribution

      real(fp), dimension(ngen) :: efield, efield_neib
      complex(fp), dimension(nsun,nsun) :: Mefield, Mefield_neib, derivative, Link_neib, G

      integer(int8)  :: i, a
      integer(int64) :: MemoryIndex, neib

      ! 1. Calculation of local contribution
      local_contribution = 0

      do MemoryIndex=1,GetMemorySize()
         if(ThisProc()==GetProc_M(MemoryIndex)) then
            if(present(ChargeDensity)) then
               G = GetG(GaugeConf,MemoryIndex,ChargeDensity)
            else
               G = GetG(GaugeConf,MemoryIndex)
            end if
            do concurrent(a=1_int8:ngen)
               if(2*Abs(GetTraceWithGenerator(a,G))>local_contribution) then
                  local_contribution = 2*Abs(GetTraceWithGenerator(a,G))
               end if
            end do
         end if
      end do
      ! 2. MPI-Sum over all partitions
      call MPI_ALLREDUCE(&
           local_contribution,&
           GetGdeviation,&
           1_intmpi,&
           GetRealSendType(),&
           MPI_MAX,&
           MPI_COMM_WORLD,mpierr)
    end function GetGdeviation

    impure function GetG(GaugeConf,MemoryIndex,ChargeDensity)
      use lattice!, only: nDim, GetNeib_M
      implicit none
      type(GaugeConfiguration), intent(in) :: GaugeConf
      integer(int64), intent(in) :: MemoryIndex
      real(fp), intent(in), optional :: ChargeDensity(:,:)

      complex(fp), dimension(nsun,nsun) :: GetG, Uneib, E, Eneib
      real(fp), dimension(ngen) :: ec,ecneib

      integer(int8) :: i,a,colour
      integer(int64) :: neib

      complex(fp) :: ChargeMatrix(nsun,nsun)

      GetG = 0
      do i=1,ndim
         neib = GetNeib_m(-i,MemoryIndex)

         Uneib = GaugeConf%Links(:,:,i,neib)

         Eneib = GetAlgebraMatrix(GaugeConf%Efield(:,i,neib))

         E = GetAlgebraMatrix(GaugeConf%Efield(:,i,MemoryIndex))

         GetG = GetG &
              + (E - matmul(matmul(conjg(transpose(Uneib)),Eneib),Uneib))/GetLatticeSpacing(i)
      end do

      if(present(ChargeDensity)) then
         ChargeMatrix = 0
         do colour=1,nsun
            ChargeMatrix(colour,colour) = ChargeDensity(colour,MemoryIndex)&
                 /product(GetLatticeSpacing([1_int8,2_int8,3_int8]))
         end do

         do a=1,ngen
            GetG = GetG &
                 - GetGenerator(a)*GetTrace(matmul(GetGenerator(a),ChargeMatrix))
         end do
      end if

      GetG = GetG/GetLatticeSpacing(0_int8)
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
             !r = GetRandomUniformReal(int(ngen,int64))!*GetLatticeSpacing(i)*GetLatticeSpacing(0_int8)
             GaugeConf%Efield(:,i,MemoryIndex) = 0
          end do
       end if
    end do

    call GaugeConf%CommunicateBoundary
  end subroutine HotInit

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
            /GetLatticeSpacing(i)/GetLatticeSpacing(0)
       do concurrent(a=1_int8:ngen)
          local_contribution = local_contribution &
               + Abs(Aimag(GetTraceWithGenerator(a,derivative)))
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
    real(fp) :: GetFieldStrengthTensor_G(nSUN,nSUN)

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
    real(fp) :: GetFieldStrengthTensor_M(nSUN,nSUN)

    complex(fp) :: Plaquette(nSUN,nSUN)

    Plaquette = GaugeConf%GetPlaquette_M(i,j,MemoryIndex)

    GetFieldStrengthTensor_M = &
         real(LogU(Plaquette)/GetLatticeSpacing(i)/GetLatticeSpacing(j)/cmplx(0,1,fp),fp)
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
          !PotentialTerm = (1-real(GetTrace(Plaquette),fp)/nSUN)&
          !     /(GetLatticeSpacing(i)*GetLatticeSpacing(j))**2

          !Local_contribution = Local_Contribution &
          !     + 2*nSUN*PotentialTerm
          PotentialTerm = (NSUN-real(GetTrace(Plaquette)))&
               /(GetLatticeSpacing(i)*GetLatticeSpacing(j))**2
          Local_Contribution = Local_Contribution + PotentialTerm
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

    GetEnergy = GetEnergy&
         *GetLatticeSpacing(1_int8)&
         *GetLatticeSpacing(2_int8)&
         *GetLatticeSpacing(3_int8)
         
  end function GetEnergy

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
         = conf%GetEfield_G(a,i,latticeindex)!/GetLatticeSpacing(i)/GetLatticeSpacing(0_int8)
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
         = GetAlgebraMatrix(efield)!/GetLatticeSpacing(i)/GetLatticeSpacing(0_int8)
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
    use lattice, only: ndim, GetMemorySize, GetProc_M, GetLatticeSpacing
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
       efield_times_dt       = GaugeConf%Efield(:,i,MemoryIndex)*StepWidth*GetLatticeSpacing(0_int8)*GetLatticeSpacing(i)
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
              !StepWidth*(GetLatticeSpacing(0)/GetLatticeSpacing(k))**2 &
              StepWidth*GetLatticeSpacing(0)/GetLatticeSpacing(i)/GetLatticeSpacing(k)**2&
              *(&
              GetUStaple(GaugeConf,i,k,MemoryIndex) + &
              GetDStaple(GaugeConf,i,k,MemoryIndex) &
              )
      end do !k

      link_times_staplesum = matmul(GaugeConf%Links(:,:,i,MemoryIndex),staplesum)

      forall(a=1:ngen)
         GaugeConf%efield(a,i,MemoryIndex)&
              = GaugeConf%efield(a,i,MemoryIndex)&
              - 2*Aimag(GetTraceWithGenerator(a,link_times_staplesum))
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

  pure function GetDivergenceOfGaugeField_M(GaugeConf,MemoryIndex) result(res)
    use lattice, only: nDim, GetNeib_M
    implicit none
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Memory index
    integer(int64),            intent(in) :: MemoryIndex
    !> Laplace operator on links
    complex(fp) :: res(nSUN,nSUN)

    integer(int8) :: i
    complex(fp) :: diff(nSUN,nSUN)

    res = 0

    do concurrent(i=1:ndim)
       diff = &
            - GaugeConf%Links(:,:,i,MemoryIndex) &
            + GaugeConf%Links(:,:,i,GetNeib_M(-i,MemoryIndex))

       ! Skew-hermitisation
       diff = (diff - conjg(transpose(diff)))/2

       res = res + diff
    end do
  end function GetDivergenceOfGaugeField_M

  !>@brief Computes the eigenvalues and -vectors (principal axes) for an Energy-Momentum-Tensor
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 03.09.2019
  !!@version 1.0
  impure subroutine GetPrincipalAxes(EMTensor,EigenValues,EigenVectors)
    use lattice, only: nDIM
    use matrixoperations, only: EigenH
    implicit none
    real(fp), intent(in)  :: EMTensor(nDIM,nDIM)

    real(fp), intent(out) :: EigenValues(nDIM)
    real(fp), intent(out) :: EigenVectors(nDIM,nDIM)

    complex(fp) :: cEMTensor(nDIM,nDIM)
    complex(fp) :: cEigenVectors(nDIM,nDIM)

    cEMTensor = cmplx(EMTensor,0,fp)

    call EigenH(cEMTensor,EigenValues,cEigenVectors,sort=.TRUE.)

    EigenVectors = real(cEigenVectors,fp)
  end subroutine GetPrincipalAxes


  !>@brief Returns field strength tensor
  !!@details: Field strength tensor in clover leaf approximation
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 04.09.2019
  !!@version 1.0
  impure function GetFieldStrengthTensor_CloverLeaf_G(GaugeConf,i,j,LatticeIndex) result(res)
    use lattice, only: GetMemoryIndex
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Direction
    integer(int8),             intent(in) :: i,j
    !> Lattice index
    integer(int64),            intent(in) :: LatticeIndex
    !> Field strength tensor
    real(fp) :: res(nSUN,nSUN)

    res = GaugeConf%GetFieldStrengthTensor_CloverLeaf_M&
         (i,j,GetMemoryIndex(LatticeIndex))
  end function GetFieldStrengthTensor_CloverLeaf_G

  !>@brief Returns field strength tensor
  !!@details: Field strength tensor in physical units \f$g F_{ij}(\vec{v})^aT^a\f$ in clover-leaf
  !! approximation
  !!@returns Field strength tensor \f$F_{i,j,\vec{v}}\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 28.02.2019
  !!@version 1.0
  impure function GetFieldStrengthTensor_CloverLeaf_M(GaugeConf,i,j,MemoryIndex) result(res)
    use lattice, only: GetLatticeSpacing
    use lattice
    use matrixoperations, only: LogU
    use mpiinterface
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Direction
    integer(int8),             intent(in) :: i,j
    !> Memory index
    integer(int64),            intent(in) :: MemoryIndex
    !> Field strength tensor
    real(fp) :: res(nSUN,nSUN)

    res= &
         real(&
         (GetQ(GaugeConf,i,j,MemoryIndex) - GetQ(GaugeConf,j,i,MemoryIndex))/(cmplx(0,8,fp))) &
         /GetLatticeSpacing(i)/GetLatticeSpacing(j)
    
    print*,GetQ(GaugeConf,0_int8,1_int8,GetMemoryIndex(GetLatticeIndex(int([1,1,5],int64))))! - GetQ(GaugeConf,j,i,MemoryIndex))/(cmplx(0,8,fp)
    call mpistop
  contains
    pure function GetQ(GaugeConf,i,j,MemoryIndex) result(res)
      implicit none
      type(GaugeConfiguration), intent(in) :: GaugeConf
      integer(int8),  intent(in) :: i,j
      integer(int64), intent(in) :: MemoryIndex

      complex(fp) :: res(nSUN,nSUN)

      res = &
           + GaugeConf%GetPlaquette_M(i,j,MemoryIndex)   & !Pij,x
           + GaugeConf%GetPlaquette_M(j,-i,MemoryIndex)  & !Pj-i,x
           + GaugeConf%GetPlaquette_M(-i,-j,MemoryIndex) & !P-i-j,x
           + GaugeConf%GetPlaquette_M(-j,i,MemoryIndex)    !P-ji,x
    end function GetQ
  end function GetFieldStrengthTensor_CloverLeaf_M

  !>@brief Returns field strength tensor
  !!@details: Field strength tensor in physical units \f$g F_{ij}(\vec{v})^aT^a\f$
  !!in order \f$O(a_i^2,a_j^2,a_i\cdot a_j)\f$
  !!@returns Spatial field strength tensor \f$F_{i,j,\vec{v}}\f$
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !!and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 23.08.2019
  !!@version 1.0
  impure function GetFieldStrengthTensorAlgebraCoordinate_G(GaugeConf,a,i,j,LatticeIndex) result(res)
    use lattice, only: GetMemoryIndex
    !> Gauge configuration
    class(GaugeConfiguration), intent(in) :: GaugeConf
    !> Gluon index
    integer(int8),             intent(in) :: a
    !> Direction
    integer(int8),             intent(in) :: i,j
    !> Lattice index
    integer(int64),            intent(in) :: LatticeIndex
    !> Field strength tensor
    real(fp) :: res

    complex(fp), dimension(nsun,nsun) :: FieldStrengthTensor

    FieldStrengthTensor = cmplx(GaugeConf%GetFieldStrengthTensor_CloverLeaf_G(i,j,LatticeIndex),0,fp)
    res = GetAlgebraCoordinate(a,FieldStrengthTensor)
  end function GetFieldStrengthTensorAlgebraCoordinate_G
end module gaugeconfiguration_su3
