!----------------------------------------------------------------------
! PROGRAMS for Lattice-NRQCD
!----------------------------------------------------------------------
!
! MODULE: programs
!>@brief Program/mains
!!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!!@date 07.03.2019
!!@version 1.0
! REVISION HISTORY:
! 07 03 2019 - Initial Version
!----------------------------------------------------------------------
module programs
  PUBLIC

contains
  impure subroutine MeasureHeavyQuarkoniumCorrelators_Equilibrium
    use, intrinsic :: iso_fortran_env
    use precision
    use mpiinterface

    use lattice
    use gaugeconfiguration_su3
    use mpi
    use io

    use nrqcd

    use tolerances, only: GetZeroTol
    use, intrinsic :: ieee_arithmetic

    implicit none

    ! Simulation parameters
    integer(int64) :: LatticeExtensions(ndim)
    real(fp)       :: LatticeSpacings(0:ndim)
    integer(int64) :: RandomNumberSeed
    
    real(fp) :: HeavyQuarkmass
    complex(fp) :: WilsonCoefficients(nWilsonCoefficients)

    ! Time coordinates
    real(fp) :: smax
    integer(int64) :: it, tsteps
    
    ! Physical fields
    type(GaugeConfiguration) :: GaugeConf_t
    type(NRQCDField)         :: HeavyField_t

    real(fp) :: Beta
    integer(int64) :: nefieldinit,nequilibriumtimesteps
    logical :: MeasureEnergy
    character(len=80) :: EnergyFilename
    
    ! Counting
    integer :: i
    
    ! Output
    integer(int8) :: FileID_Norm, FileID_Correlator

    character(len=80) :: FileMesonCorrelator, FileNorm
    
    ! Measurement of progression
    integer(int64) :: iwork, nwork
    
    call InitSimulation

    ! Measure amount of work
    tsteps = nint(abs(smax)/LatticeSpacings(0),int64)
    nwork = tsteps
    
    ! Initialisation, defining the point t=0
    call Gaugeconf_t%EquilibriumInit(Beta,nefieldinit,nequilibriumtimesteps,MeasureEnergy=MeasureEnergy,Filename=EnergyFilename)
    call HeavyField_t%InitSinglePoint(&
         latticeindex_quark=1,&
         latticeindex_antiq=1)

    ! Evolution to t
    call OpenObservableFiles
    call PrintObservables(s=0._fp,HeavyField=HeavyField_t,GaugeConf=GaugeConf_t)
    
    iwork = 0
    do it=1,tsteps
       ! Updating quarks only every "UpdateQuarksEveryNsteps" steps
       call HeavyField_t%Update(GaugeConf_t,HeavyQuarkMass,WilsonCoefficients)
       
       call GaugeConf_t%Update

       call PrintObservables(s=LatticeSpacings(0)*it,&
            HeavyField=HeavyField_t,GaugeConf=GaugeConf_t)

       iwork = iwork + 1
       if(ThisProc()==0) write(output_unit,'(F7.3,A1)') real(iwork)/nwork*100,'%'
    end do

    if(ThisProc()==0) then
       ! Closing files
       call CloseFile(FileID_Correlator)
       call CloseFile(FileID_Norm)
    end if
    
    call EndSimulation
  contains
    impure subroutine OpenObservableFiles
      ! Opening files and printing header
      if(ThisProc()==0) then
         ! Opening files
         fileID_Correlator = OpenFile(filename=FileMesonCorrelator,&
              st='REPLACE',fm='FORMATTED',act='WRITE')
         fileID_Norm = OpenFile(filename=FileNorm,&
              st='REPLACE',fm='FORMATTED',act='WRITE')

         write(FileID_Correlator,'(A1)',advance='no') 's'
         write(FileID_Correlator,'(1X,A11,1X,A11)', advance='no') 'Re(O1(1S0))','Im(O1(1S0))'
         write(FileID_Correlator,'(1X,A11,1X,A11)', advance='no') 'Re(O1(3S1))','Im(O1(3S1))'
         write(FileID_Correlator,'(1X,A11,1X,A11)', advance='no') 'Re(O8(1S0))','Im(O8(1S0))'
         write(FileID_Correlator,'(1X,A11,1X,A11)', advance='yes') 'Re(O8(3S1))','Im(O8(3S1))'

         write(fileID_Norm,'(A1,1X,A5,1X,A5)') 's','Quark','Antiq'
      end if
    end subroutine OpenObservableFiles
    
    impure subroutine PrintObservables(s,HeavyField,GaugeConf)
      implicit none
      !> relative time difference s in Wigner coordinates
      real(fp), intent(in) :: s
      !> Heavy quark field
      type(NRQCDField), intent(in) :: HeavyField
      !> Gauge configuration
      type(GaugeConfiguration), intent(in) :: GaugeConf
      
      real(fp) :: norm_quark, norm_antiq
      complex(fp) :: O1_1S0, O1_3S1, O8_1S0, O8_3S1
      
      
      norm_quark = HeavyField%GetNorm_Quark()
      norm_antiq = HeavyField%GetNorm_AntiQ()

      O1_1S0 = HeavyField%GetMesonCorrelator_O1_1s0_ZeroMomentum()
      O1_3S1 = HeavyField%GetMesonCorrelator_O1_3s1_ZeroMomentum()
      O8_1S0 = HeavyField%GetMesonCorrelator_O8_1s0_ZeroMomentum()
      O8_3S1 = HeavyField%GetMesonCorrelator_O8_3s1_ZeroMomentum()

      if(ThisProc()==0) then
         ! Norm
         write(FileID_Norm,'(3(SP,E16.9,1X))') s,norm_quark,norm_antiq


         ! Correlators
         write(FileID_Correlator,'(9(SP,E16.9,1X))') &
              s,&
              real(O1_1S0,fp), aimag(O1_1S0),&
              real(O1_3S1,fp), aimag(O1_3S1),&
              real(O8_1S0,fp), aimag(O8_1S0),&
              real(O8_3S1,fp), aimag(O8_3S1)
      end if
    end subroutine PrintObservables
    
    !>@brief Initialisation of the simulation
    !!@details
    !! MPI\n
    !! Lattice-module\n
    !! Random number generator\n
    !! etc.
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 15.02.2019
    !!@version 1.0
    impure subroutine InitSimulation
      use precision, only: fp
      use, intrinsic :: iso_fortran_env

      use mpiinterface,       only: InitModule_MPIinterface       => InitModule, ThisProc, SyncAll
      use lattice,            only: InitModule_Lattice            => InitModule, nDim
      use halocomm,           only: InitModule_HaloComm           => InitModule
      use random,             only: InitModule_Random             => InitModule
      use xpfft,              only: InitModule_xpFFT              => InitModule
      use tolerances,         only: InitModule_tolerances         => InitModule
      use NRQCD,              only: InitModule_NRQCD              => InitModule
      implicit none

      integer(int64) :: arg_count
      character(len=80) :: arg
      integer(int8) :: i, j

      real(fp) :: c_re, c_im
      real(fp) :: kspTol
      
      !..--** Reading simulation parameters **--..
      arg_count = 1

      ! Spatial lattice parameters (extensions, spacings)
      do i=1,ndim
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(I4)') LatticeExtensions(i)
      end do
      
      do i=0,ndim
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(F10.13)') LatticeSpacings(i)
      end do

      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') smax

      ! Seed for random number generator
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I4)') RandomNumberSeed
      
      ! Initial gluon distribution (box): Saturation scale
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') Beta

      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I6)') nefieldinit
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I6)') nequilibriumtimesteps
      
      ! Tolerance for iterative PETSc solver
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(E15.7)') kspTol
      
      ! Heavy quark mass
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') HeavyQuarkmass
      
      ! Wilson coefficents
      do i=1,nWilsonCoefficients
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(F10.13)') c_re
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(F10.13)') c_im

         WilsonCoefficients(i) = cmplx(c_re,c_im,fp)
      end do
      
      arg_count = arg_count +1; call get_command_argument(arg_count,FileMesonCorrelator);
      arg_count = arg_count +1; call get_command_argument(arg_count,FileNorm);
      
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I1)') j
      if(j==0) then
         MeasureEnergy=.false.
      else
         MeasureEnergy=.true.
         arg_count = arg_count +1; call get_command_argument(arg_count,EnergyFilename);
      end if

      !..--** Module initialisations **--..
      call InitModule_MPIinterface
      call InitModule_Lattice(LatticeExtensions(1:ndim),LatticeSpacings(0:ndim))
      call InitModule_HaloComm
      call InitModule_xpFFT
      call InitModule_Random(RandomNumberSeed + ThisProc())
      call InitModule_tolerances
      call InitModule_NRQCD(HeavyQuarkMass,WilsonCoefficients,1._fp,kspTol)

      call SyncAll
    end subroutine InitSimulation

    !>@brief Ending of the simulation
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 15.02.2019
    !!@version 1.0
    subroutine EndSimulation
      use mpiinterface,   only: ThisProc, FinalizeModule_MPIinterface   => FinalizeModule
      use xpfft,          only: FinalizeModule_xpFFT          => FinalizeModule
      use NRQCD,          only: FinalizeModule_NRQCD          => FinalizeModule
      implicit none

      call FinalizeModule_NRQCD
      call FinalizeModule_xpFFT
      call FinalizeModule_MPIinterface
      
      if(ThisProc()==0) write(output_unit,*) "Simulation completed"
      STOP
    end subroutine EndSimulation
  end subroutine MeasureHeavyQuarkoniumCorrelators_Equilibrium
  
  impure subroutine MeasureHeavyQuarkoniumCorrelators_free
    use, intrinsic :: iso_fortran_env
    use precision
    use mpiinterface

    use lattice
    use gaugeconfiguration_su3
    use mpi
    use io

    use nrqcd

    use tolerances, only: GetZeroTol
    use, intrinsic :: ieee_arithmetic

    implicit none

    ! Simulation parameters
    integer(int64) :: LatticeExtensions(ndim)
    real(fp)       :: LatticeSpacings(0:ndim)
    integer(int64) :: RandomNumberSeed
    
    real(fp) :: HeavyQuarkmass
    complex(fp) :: WilsonCoefficients(nWilsonCoefficients)

    ! Time coordinates
    real(fp) :: smax
    integer(int64) :: it, tsteps
    
    ! Physical fields
    type(GaugeConfiguration) :: GaugeConf_t
    type(NRQCDField)         :: HeavyField_t

    ! Counting
    integer :: i
    
    ! Output
    integer(int8) :: FileID_Norm, FileID_Correlator

    character(len=80) :: FileMesonCorrelator, FileNorm
    
    ! Measurement of progression
    integer(int64) :: iwork, nwork
    
    call InitSimulation


    ! Measure amount of work
    tsteps = nint(abs(smax)/LatticeSpacings(0),int64)
    nwork = tsteps
    
    ! Initialisation, defining the point t=0
    call Gaugeconf_t%ColdInit
    call HeavyField_t%InitSinglePoint(&
         latticeindex_quark=1,&
         latticeindex_antiq=1)

    ! Evolution to t
    call OpenObservableFiles
    call PrintObservables(s=0._fp,HeavyField=HeavyField_t,GaugeConf=GaugeConf_t)
    
    iwork = 0
    do it=1,tsteps
       ! Updating quarks only every "UpdateQuarksEveryNsteps" steps
       call HeavyField_t%Update(GaugeConf_t,HeavyQuarkMass,WilsonCoefficients)
       
       call GaugeConf_t%Update

       call PrintObservables(s=LatticeSpacings(0)*it,&
            HeavyField=HeavyField_t,GaugeConf=GaugeConf_t)

       iwork = iwork + 1
       if(ThisProc()==0) write(output_unit,'(F7.3,A1)') real(iwork)/nwork*100,'%'
    end do

    if(ThisProc()==0) then
       ! Closing files
       call CloseFile(FileID_Correlator)
       call CloseFile(FileID_Norm)
    end if
    
    call EndSimulation
  contains
    impure subroutine OpenObservableFiles
      ! Opening files and printing header
      if(ThisProc()==0) then
         ! Opening files
         fileID_Correlator = OpenFile(filename=FileMesonCorrelator,&
              st='REPLACE',fm='FORMATTED',act='WRITE')
         fileID_Norm = OpenFile(filename=FileNorm,&
              st='REPLACE',fm='FORMATTED',act='WRITE')

         write(FileID_Correlator,'(A1)',advance='no') 's'
         write(FileID_Correlator,'(1X,A11,1X,A11)', advance='no') 'Re(O1(1S0))','Im(O1(1S0))'
         write(FileID_Correlator,'(1X,A11,1X,A11)', advance='no') 'Re(O1(3S1))','Im(O1(3S1))'
         write(FileID_Correlator,'(1X,A11,1X,A11)', advance='no') 'Re(O8(1S0))','Im(O8(1S0))'
         write(FileID_Correlator,'(1X,A11,1X,A11)', advance='yes') 'Re(O8(3S1))','Im(O8(3S1))'

         write(fileID_Norm,'(A1,1X,A5,1X,A5)') 's','Quark','Antiq'
      end if
    end subroutine OpenObservableFiles
    
    impure subroutine PrintObservables(s,HeavyField,GaugeConf)
      implicit none
      !> relative time difference s in Wigner coordinates
      real(fp), intent(in) :: s
      !> Heavy quark field
      type(NRQCDField), intent(in) :: HeavyField
      !> Gauge configuration
      type(GaugeConfiguration), intent(in) :: GaugeConf
      
      real(fp) :: norm_quark, norm_antiq
      complex(fp) :: O1_1S0, O1_3S1, O8_1S0, O8_3S1
      
      
      norm_quark = HeavyField%GetNorm_Quark()
      norm_antiq = HeavyField%GetNorm_AntiQ()

      O1_1S0 = HeavyField%GetMesonCorrelator_O1_1s0_ZeroMomentum()
      O1_3S1 = HeavyField%GetMesonCorrelator_O1_3s1_ZeroMomentum()
      O8_1S0 = HeavyField%GetMesonCorrelator_O8_1s0_ZeroMomentum()
      O8_3S1 = HeavyField%GetMesonCorrelator_O8_3s1_ZeroMomentum()

      if(ThisProc()==0) then
         ! Norm
         write(FileID_Norm,'(3(SP,E16.9,1X))') s,norm_quark,norm_antiq


         ! Correlators
         write(FileID_Correlator,'(9(SP,E16.9,1X))') &
              s,&
              real(O1_1S0,fp), aimag(O1_1S0),&
              real(O1_3S1,fp), aimag(O1_3S1),&
              real(O8_1S0,fp), aimag(O8_1S0),&
              real(O8_3S1,fp), aimag(O8_3S1)
      end if
    end subroutine PrintObservables
    
    !>@brief Initialisation of the simulation
    !!@details
    !! MPI\n
    !! Lattice-module\n
    !! Random number generator\n
    !! etc.
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 15.02.2019
    !!@version 1.0
    impure subroutine InitSimulation
      use precision, only: fp
      use, intrinsic :: iso_fortran_env

      use mpiinterface,       only: InitModule_MPIinterface       => InitModule, ThisProc, SyncAll
      use lattice,            only: InitModule_Lattice            => InitModule, nDim
      use halocomm,           only: InitModule_HaloComm           => InitModule
      use random,             only: InitModule_Random             => InitModule
      use xpfft,              only: InitModule_xpFFT              => InitModule
      use tolerances,         only: InitModule_tolerances         => InitModule
      use NRQCD,              only: InitModule_NRQCD              => InitModule
      implicit none

      integer(int64) :: arg_count
      character(len=80) :: arg
      integer(int8) :: i

      real(fp) :: c_re, c_im
      real(fp) :: kspTol
      
      !..--** Reading simulation parameters **--..
      arg_count = 1

      ! Spatial lattice parameters (extensions, spacings)
      do i=1,ndim
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(I4)') LatticeExtensions(i)
      end do
      
      do i=0,ndim
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(F10.13)') LatticeSpacings(i)
      end do

      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') smax

      ! Seed for random number generator
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I4)') RandomNumberSeed

      ! Tolerance for iterative PETSc solver
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(E15.7)') kspTol
      
      ! Heavy quark mass
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') HeavyQuarkmass
      
      ! Wilson coefficents
      do i=1,nWilsonCoefficients
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(F10.13)') c_re
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(F10.13)') c_im

         WilsonCoefficients(i) = cmplx(c_re,c_im,fp)
      end do
      
      arg_count = arg_count +1; call get_command_argument(arg_count,FileMesonCorrelator);
      arg_count = arg_count +1; call get_command_argument(arg_count,FileNorm);

      !..--** Module initialisations **--..
      call InitModule_MPIinterface
      call InitModule_Lattice(LatticeExtensions(1:ndim),LatticeSpacings(0:ndim))
      call InitModule_HaloComm
      call InitModule_xpFFT
      call InitModule_Random(RandomNumberSeed + ThisProc())
      call InitModule_tolerances
      call InitModule_NRQCD(HeavyQuarkMass,WilsonCoefficients,1._fp,kspTol)

      call SyncAll
    end subroutine InitSimulation

    !>@brief Ending of the simulation
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 15.02.2019
    !!@version 1.0
    subroutine EndSimulation
      use mpiinterface,   only: ThisProc, FinalizeModule_MPIinterface   => FinalizeModule
      use xpfft,          only: FinalizeModule_xpFFT          => FinalizeModule
      use NRQCD,          only: FinalizeModule_NRQCD          => FinalizeModule
      implicit none

      call FinalizeModule_NRQCD
      call FinalizeModule_xpFFT
      call FinalizeModule_MPIinterface
      
      if(ThisProc()==0) write(output_unit,*) "Simulation completed"
      STOP
    end subroutine EndSimulation
  end subroutine MeasureHeavyQuarkoniumCorrelators_free
  
  impure subroutine MeasureHeavyQuarkoniumCorrelators_oneT

    use, intrinsic :: iso_fortran_env
    use precision
    use mpiinterface

    use lattice
    use gaugeconfiguration_su3
    use mpi
    use io

    use nrqcd

    use tolerances, only: GetZeroTol
    use, intrinsic :: ieee_arithmetic

    implicit none

    ! Simulation parameters
    integer(int64) :: LatticeExtensions(ndim)
    real(fp)       :: LatticeSpacings(0:ndim)
    integer(int64) :: RandomNumberSeed

    real(fp) :: GluonSaturationScale !qs
    real(fp) :: GluonOccupationAmplitude ! Amplitude of box in units of 1/g^2
    real(fp) :: GluonCoupling
    real(fp) :: HeavyQuarkmass
    complex(fp) :: WilsonCoefficients(nWilsonCoefficients)

    ! Time coordinates
    real(fp) :: t, smax, s, t1, t2, trange
    integer(int64) :: it, is, it1, it2, t1steps, t2steps
    
    ! Physical fields
    type(GaugeConfiguration) :: GaugeConf_t1, GaugeConf_t2
    type(NRQCDField)         :: HeavyField_t1, HeavyField_t2

    ! Counting
    integer :: i
    
    ! Output
    integer(int8) :: FileID_Norm, FileID_Correlator, FileID_Correlator_gconjg

    character(len=80) :: FileMesonCorrelator, FileNorm, FileMesonCorrelator_gconjg
    
    ! Measurement of progression
    integer(int64) :: iwork, nwork
    
    call InitSimulation

    ! Measure amount of work
    trange = smax/2
    t1steps = nint(abs(trange)/LatticeSpacings(0),int64)
    nwork = 0
    do it1=1,t1steps
       t2steps = 2*it1
       do it2=1,t2steps
          nwork = nwork + 1
       end do
    end do
    
    ! Initialisation, defining the point t=0
    call GaugeConf_t1%TransversePolarisedOccupiedInit_Box(&
         GluonSaturationScale,GluonOccupationAmplitude,GluonCoupling)
    !call Gaugeconf_t1%ColdInit
    call HeavyField_t1%InitSinglePoint(&
         latticeindex_quark=1,&
         latticeindex_antiq=1)

    ! Evolution to t
    do it=1,nint(abs(t)/LatticeSpacings(0))
       call GaugeConf_t1%Update(sign(1._fp,t))
    end do

    call OpenObservableFiles
    call PrintObservables(s=0._fp,HeavyField=HeavyField_t1,GaugeConf=GaugeConf_t1)

    iwork = 0
    do it1=1,t1steps
       t2steps = 2*it1
       HeavyField_t2 = HeavyField_t1
       GaugeConf_t2 = GaugeConf_t1

       do it2=1,t2steps
          iwork = iwork + 1
          call HeavyField_t2%Update(GaugeConf_t2,HeavyQuarkMass,WilsonCoefficients)
          call GaugeConf_t2%Update
          if(ThisProc()==0) write(output_unit,'(F7.3,A1)') real(iwork)/nwork*100,'%'
       end do

       call PrintObservables(s=LatticeSpacings(0)*t2steps,&
            HeavyField=HeavyField_t2,GaugeConf=GaugeConf_t2)

       ! Bring gauge configuration to next t1
       call GaugeConf_t1%Update(-1._fp)
    end do

    if(ThisProc()==0) then
       ! Closing files
       call CloseFile(FileID_Correlator)
       call CloseFile(FileID_Norm)
       call CloseFile(FileID_Correlator_gconjg)
    end if
    
    call EndSimulation
  contains
    impure subroutine OpenObservableFiles
      ! Opening files and printing header
      if(ThisProc()==0) then
         ! Opening files
         fileID_Correlator = OpenFile(filename=FileMesonCorrelator,&
              st='REPLACE',fm='FORMATTED',act='WRITE')
         fileID_Norm = OpenFile(filename=FileNorm,&
              st='REPLACE',fm='FORMATTED',act='WRITE')

         write(FileID_Correlator,'(A1)',advance='no') 's'
         write(FileID_Correlator,'(1X,A11,1X,A11)', advance='no') 'Re(O1(1S0))','Im(O1(1S0))'
         write(FileID_Correlator,'(1X,A11,1X,A11)', advance='no') 'Re(O1(3S1))','Im(O1(3S1))'
         write(FileID_Correlator,'(1X,A11,1X,A11)', advance='no') 'Re(O8(1S0))','Im(O8(1S0))'
         write(FileID_Correlator,'(1X,A11,1X,A11)', advance='yes') 'Re(O8(3S1))','Im(O8(3S1))'

         write(fileID_Norm,'(A1,1X,A5,1X,A5)') 's','Quark','Antiq'

         
         fileID_Correlator_gconjg = OpenFile(filename=FileMesonCorrelator_gconjg,&
              st='REPLACE',fm='FORMATTED',act='WRITE')
         write(FileID_Correlator_gconjg,'(A1)',advance='no') 's'
         write(FileID_Correlator_gconjg,'(1X,A11,1X,A11)', advance='no') 'Re(O1(1S0))','Im(O1(1S0))'
         write(FileID_Correlator_gconjg,'(1X,A11,1X,A11)', advance='no') 'Re(O1(3S1))','Im(O1(3S1))'
         write(FileID_Correlator_gconjg,'(1X,A11,1X,A11)', advance='no') 'Re(O8(1S0))','Im(O8(1S0))'
         write(FileID_Correlator_gconjg,'(1X,A11,1X,A11)', advance='yes') 'Re(O8(3S1))','Im(O8(3S1))'
      end if
    end subroutine OpenObservableFiles
    
    impure subroutine PrintObservables(s,HeavyField,GaugeConf)
      implicit none
      !> relative time difference s in Wigner coordinates
      real(fp), intent(in) :: s
      !> Heavy quark field
      type(NRQCDField), intent(in) :: HeavyField
      !> Gauge configuration
      type(GaugeConfiguration), intent(in) :: GaugeConf
      
      real(fp) :: norm_quark, norm_antiq
      complex(fp) :: O1_1S0, O1_3S1, O8_1S0, O8_3S1
      complex(fp) :: O1_1S0_gconjg, O1_3S1_gconjg, O8_1S0_gconjg, O8_3S1_gconjg
      
      
      norm_quark = HeavyField%GetNorm_Quark()
      norm_antiq = HeavyField%GetNorm_AntiQ()

      O1_1S0 = HeavyField%GetMesonCorrelator_O1_1s0_ZeroMomentum()
      O1_3S1 = HeavyField%GetMesonCorrelator_O1_3s1_ZeroMomentum()
      O8_1S0 = HeavyField%GetMesonCorrelator_O8_1s0_ZeroMomentum()
      O8_3S1 = HeavyField%GetMesonCorrelator_O8_3s1_ZeroMomentum()


      O1_1S0_gconjg = HeavyField%GetMesonCorrelator_O1_1s0_ZeroMomentum_gconjg()
      O1_3S1_gconjg = HeavyField%GetMesonCorrelator_O1_3s1_ZeroMomentum_gconjg()
      O8_1S0_gconjg = HeavyField%GetMesonCorrelator_O8_1s0_ZeroMomentum_gconjg()
      O8_3S1_gconjg = HeavyField%GetMesonCorrelator_O8_3s1_ZeroMomentum_gconjg()
      
      if(ThisProc()==0) then
         ! Norm
         write(FileID_Norm,'(3(SP,E16.9,1X))') s,norm_quark,norm_antiq

         ! Correlators
         write(FileID_Correlator,'(9(SP,E16.9,1X))') &
              s,&
              real(O1_1S0,fp), aimag(O1_1S0),&
              real(O1_3S1,fp), aimag(O1_3S1),&
              real(O8_1S0,fp), aimag(O8_1S0),&
              real(O8_3S1,fp), aimag(O8_3S1)

         ! Correlators
         write(FileID_Correlator_gconjg,'(9(SP,E16.9,1X))') &
              s,&
              real(O1_1S0_gconjg,fp), aimag(O1_1S0_gconjg),&
              real(O1_3S1_gconjg,fp), aimag(O1_3S1_gconjg),&
              real(O8_1S0_gconjg,fp), aimag(O8_1S0_gconjg),&
              real(O8_3S1_gconjg,fp), aimag(O8_3S1_gconjg)
      end if
    end subroutine PrintObservables
    
    !>@brief Initialisation of the simulation
    !!@details
    !! MPI\n
    !! Lattice-module\n
    !! Random number generator\n
    !! etc.
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 15.02.2019
    !!@version 1.0
    impure subroutine InitSimulation
      use precision, only: fp
      use, intrinsic :: iso_fortran_env

      use mpiinterface,       only: InitModule_MPIinterface       => InitModule, ThisProc, SyncAll
      use lattice,            only: InitModule_Lattice            => InitModule, nDim
      use halocomm,           only: InitModule_HaloComm           => InitModule
      use random,             only: InitModule_Random             => InitModule
      use xpfft,              only: InitModule_xpFFT              => InitModule
      use tolerances,         only: InitModule_tolerances         => InitModule
      use NRQCD,              only: InitModule_NRQCD              => InitModule
      implicit none

      integer(int64) :: arg_count
      character(len=80) :: arg
      integer(int8) :: i

      real(fp) :: c_re, c_im
      real(fp) :: kspTol
      
      !..--** Reading simulation parameters **--..
      arg_count = 1

      ! Spatial lattice parameters (extensions, spacings)
      do i=1,ndim
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(I4)') LatticeExtensions(i)
      end do
      
      do i=0,ndim
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(F10.13)') LatticeSpacings(i)
      end do


      ! Start time for quarkonium correlator (t1)
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') t
      ! t2-t1
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') smax
      
      ! Seed for random number generator
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I4)') RandomNumberSeed

      ! Initial gluon distribution (box): Saturation scale
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') GluonSaturationScale

      ! Initial gluon distribution (box): Amplitude
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') GluonOccupationAmplitude

      ! Coupling (only relevant in initialisation)
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') GluonCoupling

      ! Tolerance for iterative PETSc solver
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(E15.7)') kspTol
      
      ! Heavy quark mass
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') HeavyQuarkmass
      
      ! Wilson coefficents
      do i=1,nWilsonCoefficients
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(F10.13)') c_re
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(F10.13)') c_im

         WilsonCoefficients(i) = cmplx(c_re,c_im,fp)
      end do
      
      arg_count = arg_count +1; call get_command_argument(arg_count,FileMesonCorrelator);
      arg_count = arg_count +1; call get_command_argument(arg_count,FileNorm);
      arg_count = arg_count +1; call get_command_argument(arg_count,FileMesonCorrelator_gconjg);

      !..--** Module initialisations **--..
      call InitModule_MPIinterface
      call InitModule_Lattice(LatticeExtensions(1:ndim),LatticeSpacings(0:ndim))
      call InitModule_HaloComm
      call InitModule_xpFFT
      call InitModule_Random(RandomNumberSeed + ThisProc())
      call InitModule_tolerances
      call InitModule_NRQCD(HeavyQuarkMass,WilsonCoefficients,1._fp,kspTol)

      call SyncAll
    end subroutine InitSimulation

    !>@brief Ending of the simulation
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 15.02.2019
    !!@version 1.0
    subroutine EndSimulation
      use mpiinterface,   only: ThisProc, FinalizeModule_MPIinterface   => FinalizeModule
      use xpfft,          only: FinalizeModule_xpFFT          => FinalizeModule
      use NRQCD,          only: FinalizeModule_NRQCD          => FinalizeModule
      implicit none

      call FinalizeModule_NRQCD
      call FinalizeModule_xpFFT
      call FinalizeModule_MPIinterface
      
      if(ThisProc()==0) write(output_unit,*) "Simulation completed"
      STOP
    end subroutine EndSimulation

  end subroutine MeasureHeavyQuarkoniumCorrelators_oneT

  impure subroutine MeasureWilsonAndHybridLines_Equilibrium
    use, intrinsic :: iso_fortran_env
    use precision
    use mpiinterface
    use lattice
    use gaugeconfiguration_su3
    use mpi
    use io
    use halocomm
    use wilsonline
    use random
    use nrqcd
    use statistics

    implicit none

    ! Simulation parameters
    integer(int64) :: LatticeExtensions(ndim)
    real(fp)       :: LatticeSpacings(0:ndim)
    integer(int64) :: TimeSteps
    integer(int64) :: RandomNumberSeed

    real(fp) :: Beta
    real(fp) :: tstart
    real(fp) :: TimeRange
    
    real(fp) :: HeavyQuarkmass
    complex(fp) :: WilsonCoefficients(nWilsonCoefficients)

    real(fp) :: kspTol
    logical :: MeasureHybridLoops

    logical :: MeasureEnergy
    character(len=80) :: EnergyFilename
    integer(int64) :: nefieldinit,nequilibriumtimesteps
    
    ! Physical fields
    type(GaugeConfiguration) :: GaugeConf_t1, GaugeConf_t2, GaugeConf_initial, ColdGaugeConf, GaugeConf_old
    type(NRQCDField)         :: HeavyField, HeavyField_freeStep, HeavyField_free
    
    ! Counting
    integer :: i

    integer(int64) :: it

    ! Wilson loop parameters
    integer(int8), parameter :: messdir=nDim
    integer(int64) :: meanMaxLatticeIndex=20
    integer(int64) :: x0, xr
    integer(int64) :: rmax, r
    ! Indices
    integer(int64) :: LatticeIndex_distance, shiftstep

    ! Observables
    complex(fp), allocatable :: WilsonLoops(:,:), HybridLoops(:,:), FreeStepHybridLoops(:,:),&
         FreeHybridLoops(:,:)
    real(fp), allocatable :: qnorm(:,:), anorm(:,:)

    real(fp) :: time
    integer(int64) :: TimePoints
    ! Output
    integer(int8) :: FileID_WilsonLoops, FileID_HybridLoops, FileID_FreeStepHybridLoops,&
         FileID_FreeHybridLoops, FileID_Anorm, FileID_Qnorm

    character(len=80) :: &
         FileName_WilsonLoops, FileName_HybridLoops, FileName_FreeStepHybridLoops,&
         FileName_FreeHybridLoops, FileName_Anorm, FileName_Qnorm

    integer(intmpi) :: proc
    
    ! Measurement of progression
    integer(int64) :: iwork, nwork

    call InitSimulation

    TimePoints = nint(TimeRange/LatticeSpacings(0))

    allocate(WilsonLoops(0:TimePoints,0:rmax))
    if(MeasureHybridLoops) then
       allocate(HybridLoops(0:TimePoints,0:rmax))
       allocate(FreeStepHybridLoops(0:TimePoints,0:rmax))
       allocate(FreeHybridLoops(0:TimePoints,0:rmax))
       allocate(qnorm(0:TimePoints,0:rmax))
       allocate(anorm(0:TimePoints,0:rmax))
    end if
    
    ! initialisation of config ....
    if(MeasureHybridLoops) then
       call ColdGaugeConf%ColdInit
    end if
    
    call GaugeConf_initial%EquilibriumInit(Beta,nefieldinit,nequilibriumtimesteps,MeasureEnergy=MeasureEnergy,FileName=EnergyFilename)
    if(thisproc()==0) write(output_unit,*) 'Done: Equilibration'
    GaugeConf_t1 = GaugeConf_initial

    do it=1,abs(NINT(tstart/LatticeSpacings(0)))
       call GaugeConf_t1%Update(sign(+1._real64,tstart))
    end do

    ! Measure amount of work
    if(MeasureHybridLoops) then
       meanMaxLatticeIndex=1
       nwork = meanMaxLatticeIndex*(TimePoints+1)*(rmax+1)

       iwork = 0
       do r=0,rmax
          if(ThisProc()==0) then
             write(output_unit,*) 'r=',r,'of',rmax
          end if

          ! Initialising quark-antiquark-pair
          ! with quark at variable remote point xr
          ! and antiquark at fixed (origin) point x0
          do x0=1,meanMaxLatticeIndex
             ! Setting time of gauge conf to t1
             GaugeConf_t2 = GaugeConf_t1

             ! Setting time of quarks to t1, using distance r between quark and antiquark
             ! Getting xr
             xr = x0
             do shiftstep=1,r
                xr=GetNeib_G(messdir,xr)
             end do
             call HeavyField%InitSinglePoint(&
                  latticeindex_quark=xr,&
                  latticeindex_antiq=x0)
             HeavyField_free = HeavyField
             HeavyField_freeStep = HeavyField
             GaugeConf_old = GaugeConf_t2
             do it=0,TimeSteps,+1
                iwork = iwork + 1
                if(ThisProc()==0) write(output_unit,'(F7.3,A1)') real(iwork)/nwork*100,'%'

                if(x0==1) then
                   ! No additional average over lattice points needed
                   ! because it is already included in the routine
                   WilsonLoops(it,r) = GetWilsonLoop(GaugeConf_t1,GaugeConf_t2,r,messdir)
                end if

                ! Hybrid loop with update under dynamic gauge conf
                HybridLoops(it,r) &
                     = GetHybridLoop(GaugeConf_t1,GaugeConf_t2,HeavyField,x0,r,messdir)
                ! Hybrid loop with one update step under cold gauge conf (free step)
                FreeStepHybridLoops(it,r) &
                     = GetHybridLoop(GaugeConf_t1,GaugeConf_old,HeavyField_freeStep,x0,r,messdir)
                ! Hybrid loop with full updates under cold gauge conf only (free)
                FreeHybridLoops(it,r) &
                     = GetHybridLoop(ColdGaugeConf,ColdGaugeConf,HeavyField_free,x0,r,messdir)
                ! Norm of heavy field, evolved under interacting gauge field
                qnorm(it,r) = HeavyField%GetNorm_Quark()
                anorm(it,r) = HeavyField%GetNorm_AntiQ()

                call HeavyField%Update(GaugeConf_t2,HeavyQuarkMass,WilsonCoefficients)
                call HeavyField_freeStep%Update(ColdGaugeConf,HeavyQuarkMass,WilsonCoefficients)
                call HeavyField_free%Update(ColdGaugeConf,HeavyQuarkMass,WilsonCoefficients)
                GaugeConf_old = GaugeConf_t2
                call GaugeConf_t2%Update
             end do
          end do
       end do
    else
       nwork = TimePoints+1

       iwork = 0
       ! Setting time of gauge conf to t1
       GaugeConf_t2 = GaugeConf_t1

       do it=0,TimeSteps,+1
          do r=0,rmax

             WilsonLoops(it,r) = GetWilsonLoop(GaugeConf_t1,GaugeConf_t2,r,messdir)
          end do
          call GaugeConf_t2%Update
          iwork = iwork + 1
          if(ThisProc()==0) write(output_unit,'(F7.3,A1)') real(iwork)/nwork*100,'%'
       end do
    end if

    call PrintObservables

    call EndSimulation

  contains

    !>@brief Output of observables
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 29.05.2019
    !!@version 1.0
    impure subroutine PrintObservables
      implicit none
      
      real(fp) :: t

      character(len=3) :: adv
      
      ! Opening files and printing header
      if(ThisProc()==0) then
         ! Opening files
         FileID_WilsonLoops = OpenFile(filename=FileName_WilsonLoops,&
              st='REPLACE',fm='FORMATTED',act='WRITE')
         if(MeasureHybridLoops) then
            FileID_HybridLoops = OpenFile(filename=FileName_HybridLoops,&
                 st='REPLACE',fm='FORMATTED',act='WRITE')
            FileID_FreeStepHybridLoops = OpenFile(filename=FileName_FreeStepHybridLoops,&
                 st='REPLACE',fm='FORMATTED',act='WRITE')
            FileID_FreeHybridLoops = OpenFile(filename=FileName_FreeHybridLoops,&
                 st='REPLACE',fm='FORMATTED',act='WRITE')
            FileID_qnorm = OpenFile(filename=FileName_qnorm,&
                 st='REPLACE',fm='FORMATTED',act='WRITE')
            FileID_anorm = OpenFile(filename=FileName_anorm,&
                 st='REPLACE',fm='FORMATTED',act='WRITE')
         end if

         write(FileID_WilsonLoops,'(A1,1X)',advance='no') 't'

         if(MeasureHybridLoops) then
            write(FileID_HybridLoops,'(A1,1X)',advance='no') 't'
            write(FileID_FreeStepHybridLoops,'(A1,1X)',advance='no') 't'
            write(FileID_FreeHybridLoops,'(A1,1X)',advance='no') 't'
            write(FileID_qnorm,'(A1,1X)',advance='no') 't'
            write(FileID_anorm,'(A1,1X)',advance='no') 't'
         end if
         
         do r=0,rmax
            if(r<rmax) then
               adv='no'
            else
               adv='yes'
            end if
            
            write(FileID_WilsonLoops,'(A5,I0.2,A1,1X)',advance='no') 'Re(r=',r,')'
            write(FileID_WilsonLoops,'(A5,I0.2,A1,1X)',advance=trim(adv)) 'Im(r=',r,')'

            if(MeasureHybridLoops) then
               write(FileID_HybridLoops,'(A5,I0.2,A1,1X)',advance='no') 'Re(r=',r,')'
               write(FileID_HybridLoops,'(A5,I0.2,A1,1X)',advance=trim(adv)) 'Im(r=',r,')'
               write(FileID_FreeStepHybridLoops,'(A5,I2.2,A1,1X)',advance='no') 'Re(r=',r,')'
               write(FileID_FreeStepHybridLoops,'(A5,I2.2,A1,1X)',advance=trim(adv)) 'Im(r=',r,')'
               write(FileID_FreeHybridLoops,'(A5,I2.2,A1,1X)',advance='no') 'Re(r=',r,')'
               write(FileID_FreeHybridLoops,'(A5,I2.2,A1,1X)',advance=trim(adv)) 'Im(r=',r,')'
               write(FileID_qnorm,'(A2,I2.2,1X)',advance=trim(adv)) 'r=',r
               write(FileID_anorm,'(A2,I2.2,1X)',advance=trim(adv)) 'r=',r
            end if
         end do

         do it=lbound(WilsonLoops,1),ubound(WilsonLoops,1)
            time = tstart + it*LatticeSpacings(0)
            
            write(FileID_WilsonLoops,'(SP,E16.9,1X)',advance='no') time

            if(MeasureHybridLoops) then
               write(FileID_HybridLoops,'(SP,E16.9,1X)',advance='no') time
               write(FileID_FreeStepHybridLoops,'(SP,E16.9,1X)',advance='no') time
               write(FileID_FreeHybridLoops,'(SP,E16.9,1X)',advance='no') time
               write(FileID_qnorm,'(SP,E16.9,1X)',advance='no') time
               write(FileID_anorm,'(SP,E16.9,1X)',advance='no') time
            end if
            do r=lbound(WilsonLoops,2),ubound(WilsonLoops,2)
               if(r<ubound(WilsonLoops,2)) then
                  adv='no'
               else
                  adv='yes'
               end if
               
               write(FileID_WilsonLoops,'(2(SP,E16.9,1X))',advance=trim(adv)) &
                    real(WilsonLoops(it,r)),aimag(WilsonLoops(it,r))

               if(MeasureHybridLoops) then
                  write(FileID_HybridLoops,'(2(SP,E16.9,1X))',advance=trim(adv)) &
                       real(HybridLoops(it,r)),aimag(HybridLoops(it,r))

                  write(FileID_FreeStepHybridLoops,'(2(SP,E16.9,1X))',advance=trim(adv)) &
                       real(FreeStepHybridLoops(it,r)),aimag(FreeStepHybridLoops(it,r))

                  write(FileID_FreeHybridLoops,'(2(SP,E16.9,1X))',advance=trim(adv)) &
                       real(FreeHybridLoops(it,r)),aimag(FreeHybridLoops(it,r))

                  write(FileID_qnorm,'(SP,E16.9,1X)',advance=trim(adv)) &
                       qnorm(it,r)
                  write(FileID_anorm,'(SP,E16.9,1X)',advance=trim(adv)) &
                       anorm(it,r)
               end if
            end do
         end do

         call CloseFile(FileID_WilsonLoops)
         if(MeasureHybridLoops) then
            call CloseFile(FileID_HybridLoops)
            call CloseFile(FileID_FreeStepHybridLoops)
            call CloseFile(FileID_FreeHybridLoops)
            call CloseFile(FileID_qnorm)
            call CloseFile(FileID_anorm)
         end if
      end if
    end subroutine PrintObservables

    !>@brief Initialisation of the simulation
    !!@details
    !! MPI\n
    !! Lattice-module\n
    !! Random number generator\n
    !! etc.
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 15.02.2019
    !!@version 1.0
    impure subroutine InitSimulation
      use precision, only: fp
      use, intrinsic :: iso_fortran_env

      use mpiinterface,       only: InitModule_MPIinterface       => InitModule, ThisProc, SyncAll
      use lattice,            only: InitModule_Lattice            => InitModule, nDim
      use halocomm,           only: InitModule_HaloComm           => InitModule
      use random,             only: InitModule_Random             => InitModule
      use xpfft,              only: InitModule_xpFFT              => InitModule
      use tolerances,         only: InitModule_tolerances         => InitModule
      use NRQCD,              only: InitModule_NRQCD              => InitModule
      implicit none

      integer(int64) :: arg_count
      character(len=80) :: arg
      integer(int8) :: i, j
      
      real(fp) :: c_re, c_im
      
      !..--** Reading simulation parameters **--..
      arg_count = 1

      ! Spatial lattice parameters (extensions, spacings)
      do i=1,ndim
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(I4)') LatticeExtensions(i)
      end do

      do i=0,ndim
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(F10.13)') LatticeSpacings(i)
      end do

      ! start time
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') tstart
      ! Center of mass time-range smax
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') TimeRange

      TimeSteps=ceiling(TimeRange/LatticeSpacings(0))

      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I4)') rmax
      
      ! Seed for random number generator
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I4)') RandomNumberSeed
      ! Initial gluon distribution (box): Saturation scale
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') Beta

      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I6)') nefieldinit
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I6)') nequilibriumtimesteps
      
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I1)') j
      if(j==0) then
         MeasureHybridLoops=.false.
      else
         MeasureHybridLoops=.true.
      end if
      
      ! Tolerance for iterative PETSc solver
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(E15.7)') kspTol
      
      ! Heavy quark mass
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') HeavyQuarkmass
      
      ! Wilson coefficents
      do i=1,nWilsonCoefficients
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(F10.13)') c_re
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(F10.13)') c_im

         WilsonCoefficients(i) = cmplx(c_re,c_im,fp)
      end do
      
      ! Output filenames
      arg_count = arg_count +1; call get_command_argument(arg_count,FileName_WilsonLoops);
      arg_count = arg_count +1; call get_command_argument(arg_count,FileName_HybridLoops);
      arg_count = arg_count +1; call get_command_argument(arg_count,FileName_FreeStepHybridLoops);
      arg_count = arg_count +1; call get_command_argument(arg_count,FileName_FreeHybridLoops);
      arg_count = arg_count +1; call get_command_argument(arg_count,FileName_qnorm);
      arg_count = arg_count +1; call get_command_argument(arg_count,FileName_anorm);
      
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I1)') j
      if(j==0) then
         MeasureEnergy=.false.
      else
         MeasureEnergy=.true.
         arg_count = arg_count +1; call get_command_argument(arg_count,EnergyFilename);
      end if
      
      !..--** Module initialisations **--..
      call InitModule_MPIinterface
      call InitModule_Lattice(LatticeExtensions(1:ndim),LatticeSpacings(0:ndim))
      call InitModule_HaloComm
      call InitModule_xpFFT
      call InitModule_Random(RandomNumberSeed + ThisProc())
      call InitModule_tolerances
      call InitModule_NRQCD(HeavyQuarkMass,WilsonCoefficients,1._fp,kspTol)

      call SyncAll
    end subroutine InitSimulation

    !>@brief Ending of the simulation
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 15.02.2019
    !!@version 1.0
    subroutine EndSimulation
      use mpiinterface,   only: ThisProc, FinalizeModule_MPIinterface   => FinalizeModule
      use xpfft,          only: FinalizeModule_xpFFT          => FinalizeModule
      use NRQCD,          only: FinalizeModule_NRQCD          => FinalizeModule
      implicit none

      call FinalizeModule_NRQCD
      call FinalizeModule_xpFFT
      call FinalizeModule_MPIinterface

      if(ThisProc()==0) write(output_unit,*) "Simulation completed"
      STOP
    end subroutine EndSimulation
  end subroutine MeasureWilsonAndHybridLines_Equilibrium

  impure subroutine MeasureWilsonLines_Equilibrium_StaticMeson
    use, intrinsic :: iso_fortran_env
    use precision
    use mpiinterface
    use lattice
    use gaugeconfiguration_su3
    use mpi
    use io
    use halocomm
    use wilsonline
    use random
    use nrqcd
    use statistics
    use su3, only: nsun


    use matrixoperations

    implicit none

    ! Simulation parameters
    integer(int64) :: LatticeExtensions(ndim)
    real(fp)       :: LatticeSpacings(0:ndim)
    integer(int64) :: TimeSteps
    integer(int64) :: RandomNumberSeed

    real(fp) :: Beta
    real(fp) :: tstart
    real(fp) :: TimeRange

    logical :: MeasureEnergy
    character(len=80) :: EnergyFilename
    integer(int64) :: nefieldinit,nequilibriumtimesteps
    
    ! Physical fields
    type(GaugeConfiguration) :: GaugeConf
    
    ! Counting
    integer :: i

    integer(int64) :: it

    integer(intmpi) :: mpierr
    
    ! Wilson loop parameters
    integer(int8), parameter :: messdir=nDim
    integer(int64) :: meanMaxLatticeIndex=20
    integer(int64) :: position_quark, position_antiquark
    integer(int8) :: rmax, r, k
    ! Indices
    integer(int64) :: LatticeIndex_distance, shiftstep
    integer(int8) :: colour

    real(fp) :: time
    integer(int64) :: TimePoints

    integer(intmpi) :: proc
    
    ! Measurement of progression
    integer(int64) :: iwork, nwork

    ! Charge density of static meson
    real(fp), allocatable :: ChargeDensity(:,:)
    real(fp) :: tolerance_Eprojection
    
    real(fp) :: localcharge,totalcharge
    
    call InitSimulation

    TimePoints = nint(TimeRange/LatticeSpacings(0))

    ! Initialising the gauge configuration for the given charge density
    call GaugeConf%EquilibriumInit(&
         Beta,nefieldinit,nequilibriumtimesteps,tolerance_Eprojection,&
         ChargeDensity,MeasureEnergy,EnergyFilename)

    ! Performing update steps while computing energy tensor

    call EndSimulation

  contains

    !>@brief Output of observables
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 29.05.2019
    !!@version 1.0
    impure subroutine PrintObservables
      implicit none
      
      real(fp) :: t

      character(len=3) :: adv

    end subroutine PrintObservables

    !>@brief Initialisation of the simulation
    !!@details
    !! MPI\n
    !! Lattice-module\n
    !! Random number generator\n
    !! etc.
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 15.02.2019
    !!@version 1.0
    impure subroutine InitSimulation
      use precision, only: fp
      use, intrinsic :: iso_fortran_env

      use mpiinterface,       only: InitModule_MPIinterface       => InitModule, ThisProc, SyncAll
      use lattice,            only: InitModule_Lattice            => InitModule, nDim
      use halocomm,           only: InitModule_HaloComm           => InitModule
      use random,             only: InitModule_Random             => InitModule
      use xpfft,              only: InitModule_xpFFT              => InitModule
      use tolerances,         only: InitModule_tolerances         => InitModule
      implicit none

      integer(int64) :: arg_count
      character(len=80) :: arg
      integer(int8) :: i, j
      
      real(fp) :: c_re, c_im
      character(len=100) :: FileName_ChargeDensity
      
      !..--** Reading simulation parameters **--..
      arg_count = 1

      ! Spatial lattice parameters (extensions, spacings)
      do i=1,ndim
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(I4)') LatticeExtensions(i)
      end do

      do i=0,ndim
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(F10.13)') LatticeSpacings(i)
      end do

      ! start time
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') tstart
      ! Center of mass time-range smax
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') TimeRange

      TimeSteps=ceiling(TimeRange/LatticeSpacings(0))

      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I4)') rmax
      
      ! Seed for random number generator
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I4)') RandomNumberSeed
      ! Initial gluon distribution (box): Saturation scale
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') Beta

      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(E15.7)') tolerance_Eprojection
      
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I6)') nefieldinit
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I6)') nequilibriumtimesteps
      
      ! Output filenames
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I1)') j
      if(j==0) then
         MeasureEnergy=.false.
      else
         MeasureEnergy=.true.
         arg_count = arg_count +1; call get_command_argument(arg_count,EnergyFilename);
      end if

       arg_count = arg_count +1; call get_command_argument(arg_count,FileName_ChargeDensity);
      
      !..--** Module initialisations **--..
      call InitModule_MPIinterface
      call InitModule_Lattice(LatticeExtensions(1:ndim),LatticeSpacings(0:ndim))
      call InitModule_HaloComm
      call InitModule_xpFFT
      call InitModule_Random(RandomNumberSeed + ThisProc())
      call InitModule_tolerances

      ! Read charge density
      call ReadChargeDensity(ChargeDensity,FileName_ChargeDensity)

      call SyncAll
    end subroutine InitSimulation
    
    impure subroutine ReadChargeDensity(ChargeDensity,FileName_ChargeDensity)
      use mpiinterface
      use io
      implicit none
      real(fp), allocatable, intent(out) :: ChargeDensity(:,:)
      character(len=*), intent(in) :: FileName_ChargeDensity

      integer(int64) :: position(nDim), latticeindex
      integer(int8) :: colour
      real(fp), dimension(nSUN) :: charges

      integer(int8) :: FileID_ChargeDensity
      integer(int64) :: nLines, iLine
    
      allocate(ChargeDensity(nsun,GetMemorySize()))
      ChargeDensity = 0

      FileID_ChargeDensity = OpenFile(FileName_ChargeDensity,st='OLD',fm='FORMATTED',act='READ')
      ! Determine number of entries in chargedensity-file
      nLines = 0
      do
         READ(FileID_ChargeDensity,*,END=10)
         nLines = nLines + 1
      end do
10    REWIND(FileID_ChargeDensity)
      
      ! Read in charge density
      do iLine=1,nLines
         READ(FileID_ChargeDensity,fmt=*) position, charges

         LatticeIndex = GetLatticeIndex(position)

         if(GetProc_G(LatticeIndex)==ThisProc()) then
            ChargeDensity(:,GetMemoryIndex(LatticeIndex)) = &
                 ChargeDensity(:,GetMemoryIndex(LatticeIndex)) + Charges
         end if
      end do
    end subroutine ReadChargeDensity

    !>@brief Ending of the simulation
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 15.02.2019
    !!@version 1.0
    subroutine EndSimulation
      use mpiinterface,   only: ThisProc, FinalizeModule_MPIinterface   => FinalizeModule
      use xpfft,          only: FinalizeModule_xpFFT          => FinalizeModule
      implicit none

      call FinalizeModule_xpFFT
      call FinalizeModule_MPIinterface

      if(ThisProc()==0) write(output_unit,*) "Simulation completed"
      STOP
    end subroutine EndSimulation
  end subroutine MeasureWilsonLines_Equilibrium_StaticMeson
  
  
  !>@brief Program for measuring the various wilson lines
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 09.06.2019
  !!@version 1.0
  impure subroutine MeasureWilsonLines_nonEquilibrium
    use, intrinsic :: iso_fortran_env
    use precision
    use mpiinterface
    use lattice
    use gaugeconfiguration_su3
    use mpi
    use io
    use halocomm
    use wilsonline
    use random
    use nrqcd
    
    implicit none

    ! Simulation parameters
    integer(int64) :: LatticeExtensions(ndim)
    real(fp)       :: LatticeSpacings(0:ndim)
    
    integer(int64) :: RandomNumberSeed
    
    real(fp) :: GluonSaturationScale !qs
    real(fp) :: GluonOccupationAmplitude ! Amplitude of box in units of 1/g^2
    real(fp) :: GluonCoupling

    complex(fp) :: WilsonCoefficients(nWilsonCoefficients)

    ! Time coordinates
    real(fp) :: t, smax, s, t1, t2, trange
    integer(int64) :: it, is, it1, it2, t1steps, t2steps

    ! Physical fields
    type(GaugeConfiguration) :: GaugeConf_t1, GaugeConf_t2

    ! Counting
    integer :: i

    integer(int64) :: rmax
    
    ! Measurement of progression
    integer(int64) :: iwork, nwork

    character(len=80) :: FileName_WilsonLoops
    integer(int8) :: FileID_WilsonLoops

    call InitSimulation
    
    ! Measure amount of work
    trange = smax/2
    t1steps = nint(abs(trange)/LatticeSpacings(0),int64)
    nwork = 0
    do it1=1,t1steps
       t2steps = 2*it1
       do it2=1,t2steps
          nwork = nwork + 1
       end do
    end do
    
    ! Initialisation, defining the point t=0
    call GaugeConf_t1%TransversePolarisedOccupiedInit_Box(&
         GluonSaturationScale,GluonOccupationAmplitude,GluonCoupling)
    
    ! Evolution to t
    do it=1,nint(abs(t)/LatticeSpacings(0))
       call GaugeConf_t1%Update(sign(1._fp,t))
    end do

    call OpenObservableFiles
    call PrintObservables(s=0._fp, GaugeConf_t1=GaugeConf_t1, GaugeConf_t2=GaugeConf_t1)

    iwork=0

    do it1=1,t1steps
       t2steps = 2*it1
       GaugeConf_t2 = GaugeConf_t1

       do it2=1,t2steps
          iwork = iwork + 1

          call GaugeConf_t2%Update
          if(ThisProc()==0) write(output_unit,'(F7.3,A1)') real(iwork)/nwork*100,'%'
       end do
       call PrintObservables(s=LatticeSpacings(0)*t2steps,&
            GaugeConf_t1=GaugeConf_t1, GaugeConf_t2=GaugeConf_t2)
       
       ! Bring gauge configuration to next t1
       call GaugeConf_t1%Update(-1._fp)
    end do

    call CloseObservableFiles
    call EndSimulation
  contains
    impure subroutine OpenObservableFiles
      implicit none
      integer(int8) :: r
      character(len=3) :: adv
      
      ! Opening files and printing header
      if(ThisProc()==0) then
         ! Opening files
         fileID_WilsonLoops = OpenFile(filename=FileName_WilsonLoops,&
              st='REPLACE',fm='FORMATTED',act='WRITE')

         write(FileID_WilsonLoops,'(A1,1X)',advance='no') 's'
         do r=0,rmax
            if(r<rmax) then
               adv='no'
            else
               adv='yes'
            end if
            
            write(FileID_WilsonLoops,'(A5,I0.2,A1,1X)',advance='no') 'Re(r=',r,')'
            write(FileID_WilsonLoops,'(A5,I0.2,A1,1X)',advance=trim(adv)) 'Im(r=',r,')'
         end do
      end if
    end subroutine OpenObservableFiles

    impure subroutine CloseObservableFiles
      implicit none
      if(ThisProc()==0) then
         ! Closing files
         call CloseFile(FileID_WilsonLoops)
      end if
    end subroutine CloseObservableFiles

    impure subroutine PrintObservables(s,GaugeConf_t1,GaugeConf_t2)
      implicit none
      real(fp), intent(in) :: s
      type(GaugeConfiguration), intent(in) :: GaugeConf_t1, GaugeConf_t2
      integer(int64) :: r
      integer(int8) :: messdir
      character(len=3) :: adv

      complex(fp) :: WilsonLoop
      
      
      if(ThisProc()==0) then
         write(FileID_WilsonLoops,'(SP,E16.9,1X)',advance='no') s
      end if
      
      do r=0,rmax
         WilsonLoop = 0
         do messdir=1,ndim
            WilsonLoop = WilsonLoop + GetWilsonLoop(GaugeConf_t1,GaugeConf_t2,r,messdir)/ndim
         end do

         if(ThisProc()==0) then
            if(r<rmax) then
               adv='no'
            else
               adv='yes'
            end if
            write(FileID_WilsonLoops,'(2(SP,E16.9,1X))',advance=trim(adv)) &
                 real(WilsonLoop),aimag(WilsonLoop)
         end if
      end do
    end subroutine PrintObservables

    !>@brief Initialisation of the simulation
    !!@details
    !! MPI\n
    !! Lattice-module\n
    !! Random number generator\n
    !! etc.
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 15.02.2019
    !!@version 1.0
    impure subroutine InitSimulation
      use precision, only: fp
      use, intrinsic :: iso_fortran_env

      use mpiinterface,       only: InitModule_MPIinterface       => InitModule, ThisProc, SyncAll
      use lattice,            only: InitModule_Lattice            => InitModule, nDim
      use halocomm,           only: InitModule_HaloComm           => InitModule
      use random,             only: InitModule_Random             => InitModule
      use xpfft,              only: InitModule_xpFFT              => InitModule
      use tolerances,         only: InitModule_tolerances         => InitModule
      implicit none

      integer(int64) :: arg_count
      character(len=80) :: arg
      integer(int8) :: i
      
      !..--** Reading simulation parameters **--..
      arg_count = 1

      ! Spatial lattice parameters (extensions, spacings)
      do i=1,ndim
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(I4)') LatticeExtensions(i)
      end do
      
      do i=0,ndim
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(F10.13)') LatticeSpacings(i)
      end do

      ! Mean time t in (s,t) diagram
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') t
      ! t2-t1
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') smax

      ! rmax
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I4)') rmax

      ! Seed for random number generator
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I4)') RandomNumberSeed

      ! Initial gluon distribution (box): Saturation scale
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') GluonSaturationScale

      ! Initial gluon distribution (box): Amplitude
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') GluonOccupationAmplitude

      ! Coupling (only relevant in initialisation)
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') GluonCoupling
      
      arg_count = arg_count +1; call get_command_argument(arg_count,FileName_WilsonLoops);

      !..--** Module initialisations **--..
      call InitModule_MPIinterface
      call InitModule_Lattice(LatticeExtensions(1:ndim),LatticeSpacings(0:ndim))
      call InitModule_HaloComm
      call InitModule_xpFFT
      call InitModule_Random(RandomNumberSeed + ThisProc())
      call InitModule_tolerances

      call SyncAll
    end subroutine InitSimulation

    !>@brief Ending of the simulation
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 15.02.2019
    !!@version 1.0
    subroutine EndSimulation
      use mpiinterface,   only: ThisProc, FinalizeModule_MPIinterface   => FinalizeModule
      use xpfft,          only: FinalizeModule_xpFFT          => FinalizeModule
      implicit none

      call FinalizeModule_xpFFT
      call FinalizeModule_MPIinterface

      if(ThisProc()==0) write(output_unit,*) "Simulation completed"
      STOP
    end subroutine EndSimulation
  end subroutine MeasureWilsonLines_NonEquilibrium

  !>@brief Program for computing the spectral function out of a correlator
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.03.2019
  !!@version 1.0
  impure subroutine ComputeSpectrum
    use, intrinsic :: iso_fortran_env
    use precision
    use mpiinterface
    use io
    use windowing
    use twfft

    use mathconstants, only: pi
    implicit none
    character(len=200) :: FileName_input, FileName_output
    integer :: nHeaderLines
    complex(fp), allocatable :: Correlator(:)
    real(fp), allocatable :: Spectrum(:)

    real(fp) :: dt, MaxTime

    integer :: it
    integer(int8) :: fileID
    
    call InitSimulation

    if(ThisProc()==0) then
       write(output_unit,*)'**********************'
       write(output_unit,*)'Extraction of spectrum'
       write(output_unit,*)'**********************'
       write(output_unit,*)
       write(output_unit,*)'Input file:  ',trim(FileName_input)
       write(output_unit,*)'Output file: ',trim(FileName_output)
       call flush(output_unit)

       call ReadSignal(FileName_input,nHeaderLines,Correlator,dt)
       if(ThisProc()==0) then
          fileID = OpenFile(filename='input.txt',fm='formatted',act='write',st='replace')

          do it=1,size(correlator)
             write(FileID,*) it,real(correlator(it)),aimag(correlator(it))
          end do

          call CloseFile(fileID)
       end if

       
       call Symmetrisation(correlator)
       if(ThisProc()==0) then
          fileID = OpenFile(filename='symm.txt',fm='formatted',act='write',st='replace')

          do it=1,size(correlator)
             write(FileID,*) it,real(correlator(it)),aimag(correlator(it))
          end do

          call CloseFile(fileID)
       end if

       
       call ApplyHannWindow(Correlator)
       if(ThisProc()==0) then
          fileID = OpenFile(filename='hanned.txt',fm='formatted',act='write',st='replace')

          do it=1,size(correlator)
             write(FileID,*) it,real(correlator(it)),aimag(correlator(it))
          end do

          call CloseFile(fileID)
       end if

       call Symm2FFTformat(Correlator)
       if(ThisProc()==0) then
          fileID = OpenFile(filename='fftformat.txt',fm='formatted',act='write',st='replace')
          do it=1,size(correlator)
             write(FileID,*) it,real(correlator(it)),aimag(correlator(it))
          end do
          call CloseFile(fileID)
       end if

       call t2w(Correlator,dt)

       ! Extract spectrum
       allocate(spectrum(size(correlator)))
       
       spectrum = 2*aimag(correlator)

       call WriteSpectrum2File(Spectrum,FileName_output,dt)
    end if
    call EndSimulation

  contains
    
    subroutine Symmetrisation(signal)
      implicit none
      complex(real64), allocatable, intent(inout) :: signal(:)
      complex(real64), allocatable :: signal_input(:)

      integer :: i,k

      if(.not.allocated(signal)) then
         STOP 'Signal not allocated'
      end if
      allocate(signal_input(size(signal)))
      signal_input = signal
      deallocate(signal)
      allocate(signal(size(signal_input)*2-2))

      
      do i=2,size(signal_input)
         k = size(signal_input) -(i-1)
         signal(k) = -conjg(signal_input(i))
      end do
      do i=1,size(signal_input)-1
         k = i + size(signal_input) - 1
         signal(k) = signal_input(i)
      end do
    end subroutine Symmetrisation
    
    !>@brief Transforms symmetrized signal into FFT format
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 22.03.2019
    !!@version 1.0
    impure subroutine Symm2FFTformat(signal)
      implicit none
      complex(fp), allocatable, intent(inout) :: signal(:)
      complex(fp), allocatable :: signal_input(:)

      integer :: i,k

      allocate(signal_input(size(signal)))
      signal_input = signal

      signal(1:size(signal)/2) = signal_input(size(signal)/2+1:size(signal))

      do i=1,size(signal)/2
         k = size(signal)/2 + i 

         signal(k) = signal_input(i)
      end do
    end subroutine Symm2FFTformat
    
    !>@brief Reads signal/correlator from file
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 22.03.2019
    !!@version 1.0
    impure subroutine ReadSignal(FileName,nHeaderLines,Signal,dt)
      use io
      implicit none
      character(len=*),         intent(in)  :: FileName
      integer,                  intent(in)  :: nHeaderLines
      complex(fp), allocatable, intent(out) :: Signal(:)
      real(fp),                 intent(out) :: dt

      integer :: i, SignalSize
      integer(int8) :: FileID
      integer :: stat
      real(fp) :: t,signal_re,signal_im
      real(fp) :: t2
      
      ! 1. Open file
      FileID = OpenFile(filename=FileName,fm='formatted',act='read',st='old')
      
      ! 2. Determine length of signal
      do i=1,nHeaderLines
         read(unit=FileID,fmt=*)
      end do
      SignalSize = 0
      do
         read(unit=FileID,fmt=*,iostat=stat) t, signal_re, signal_im
         if(t>MaxTime) exit
         if(stat/=0) exit
         SignalSize = SignalSize + 1
      end do
      call CloseFile(FileID)
      
      ! 3. Read signal
      allocate(Signal(SignalSize))
      FileID = OpenFile(filename=FileName,fm='formatted',act='read',st='old')
      
      do i=1,nHeaderLines
         read(unit=FileID,fmt=*)
      end do

      do i=1,SignalSize
         read(unit=FileID,fmt=*,iostat=stat) t, signal_re, signal_im
         if(i==1) then
            t2 = t
         elseif(i==2) then
            dt = t - t2
         end if

         Signal(i) = cmplx(signal_re,signal_im,fp)
      end do
      
      call CloseFile(FileID) 
    end subroutine ReadSignal

    !>@brief Writes spectrum to file
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 22.03.2019
    !!@version 1.0
    impure subroutine WriteSpectrum2File(Spectrum,FileName,dt)
      use io
      use mathconstants, only: pi
      implicit none
      real(fp),        intent(in) :: Spectrum(:)
      character(len=*),intent(in) :: FileName
      real(fp),        intent(in) :: dt

      integer(int8) :: FileID

      integer(int64) :: i,k
      
      real(fp) :: wmin, wmax, dw
      real(fp), parameter :: twopi = 2*pi
      
      fileID = OpenFile(filename=filename,fm='formatted',act='write',st='replace')

      wmin = -pi/dt
      wmax = +pi/dt
      dw   = twopi/(dt*size(spectrum))

      ! Starting with negative frequencies...
      do i=size(spectrum)/2+1,size(spectrum)
         k = i-size(spectrum)/2 - 1
         write(fileID,fmt='(2(SP,E13.6,1X))') wmin + k*dw,spectrum(i)
         !print*,'-',i,k,wmin+k*dw
      end do

      ! ...continueing with positive frequencies
      do i=1,size(spectrum)/2
         k = i -1
         write(fileID,fmt='(2(SP,E13.6,1X))') k*dw,spectrum(i)
         !print*,'+',i,k,k*dw
      end do
      call CloseFile(fileID)
    end subroutine WriteSpectrum2File

    !>@brief Initialisation of the simulation
    !!@details
    !! MPI\n
    !! Lattice-module\n
    !! Random number generator\n
    !! etc.
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 22.03.2019
    !!@version 1.0
    impure subroutine InitSimulation
      use precision, only: fp
      use, intrinsic :: iso_fortran_env

      use mpiinterface, only: InitModule_MPIinterface       => InitModule, ThisProc, SyncAll
      implicit none

      integer(int64) :: arg_count
      character(len=200) :: arg
      
      !..--** Reading simulation parameters **--..
      arg_count = 1

      ! Number of lines in header
      arg_count = arg_count + 1; call get_command_argument(arg_count,arg);
      read(arg,'(I4)') nheaderlines

      ! Input filename
      arg_count = arg_count + 1; call get_command_argument(arg_count,arg);
      FileName_input = arg

      ! Output filename
      arg_count = arg_count + 1; call get_command_argument(arg_count,arg);
      FileName_output = arg

      ! maximum time
      arg_count = arg_count + 1; call get_command_argument(arg_count,arg);
      read(arg,'(F6.1)') MaxTime
      

      !..--** Module initialisations **--..
      call InitModule_MPIinterface

      call SyncAll
    end subroutine InitSimulation

    !>@brief Ending of the simulation
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 22.03.2019
    !!@version 1.0
    subroutine EndSimulation
      use mpiinterface, only: ThisProc, FinalizeModule_MPIinterface   => FinalizeModule
      implicit none

      call FinalizeModule_MPIinterface

      if(ThisProc()==0) write(output_unit,*) "Simulation completed"
      STOP
    end subroutine EndSimulation
  end subroutine ComputeSpectrum
  
  !>@brief Program for measuring the heavy quarkonium correlator
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 07.03.2019
  !!@version 1.0
  impure subroutine MeasureHeavyQuarkoniumCorrelators
    use, intrinsic :: iso_fortran_env
    use precision
    use mpiinterface

    use lattice
    use gaugeconfiguration_su3
    use mpi
    use io

    use nrqcd

    use tolerances, only: GetZeroTol
    use, intrinsic :: ieee_arithmetic

    implicit none

    ! Simulation parameters
    integer(int64) :: LatticeExtensions(ndim)
    real(fp)       :: LatticeSpacings(0:ndim)
    integer(int64) :: RandomNumberSeed

    real(fp) :: GluonSaturationScale !qs
    real(fp) :: GluonOccupationAmplitude ! Amplitude of box in units of 1/g^2
    real(fp) :: GluonCoupling
    real(fp) :: tmin
    real(fp) :: tmax
    real(fp) :: smax
    real(fp) :: HeavyQuarkmass
    complex(fp) :: WilsonCoefficients(nWilsonCoefficients)

    ! Time coordinates
    real(fp) :: t, s, t1, t2, ds
    integer(int64) :: is, it, it1, it2
    integer(int64) :: itmin, itmax, idt
    integer(int64) :: ismin, ismax, ids

    real(fp) :: t2min, t2max, t2range
    integer(int64) :: it2min, it2max, it2range

    real(fp) :: t1min, t1max, t1range
    integer(int64) :: it1min, it1max, it1range
    
    ! Physical fields
    type(GaugeConfiguration) :: GaugeConf_t1, GaugeConf_t2
    type(NRQCDField)         :: HeavyField_t1, HeavyField_t2

    ! Counting
    integer :: i

    ! Output
    real(fp) :: norm_quark, norm_antiq
    complex(fp) :: mesoncorrelator
    integer(int8) :: FileID_Norm, FileID_Correlator

    character(len=80) :: FileMesonCorrelator, FileNorm

    complex(fp), allocatable :: correlator(:,:)
    real(fp),allocatable :: qnorm(:,:),anorm(:,:)

    real(fp) :: nan

    real(fp) :: dttest, dstest
    integer(int64) :: ittest, istest
    logical :: tTest, sTest

    ! Measurement of progression
    integer(int64) :: iprogress, nprogress

    ! Format builder
    character(len=100) :: form
    character(len=1) :: auxd, auxw
    integer :: d, p, w

    nan = ieee_value(nan, ieee_quiet_nan)
    
    call InitSimulation

    ! Allocating result arrays
    ismin = 0
    ismax = nint(smax/ds)
    itmin = nint(tmin/ds)
    itmax = nint(tmax/ds)
    
    allocate(correlator(ismin:ismax,itmin:itmax))
    correlator = cmplx(nan,nan,fp)
    correlator(0,:) = cmplx(0,1,fp)

    allocate(qnorm(ismin:ismax,itmin:itmax))
    qnorm = nan
    qnorm(0,:) = 1
    allocate(anorm(ismin:ismax,itmin:itmax))
    anorm = qnorm

    !goto 1
    ! Initialisation, defining the point t=0
    call GaugeConf_t1%TransversePolarisedOccupiedInit_Box(&
         GluonSaturationScale,GluonOccupationAmplitude,GluonCoupling)
    !call Gaugeconf_t1%ColdInit
    call HeavyField_t1%InitSinglePoint(&
         latticeindex_quark=1,&
         latticeindex_antiq=1)

    ! Evolve gauge configuration from t1=0 (initialisation time) to minimum t1-value
    t1min = tmin - smax/2
    t1max = tmax ! ... srange=0 at highest t-value
    it1min = nint(t1min/LatticeSpacings(0))
    it1max = nint(t1max/LatticeSpacings(0))

    t1 = 0
    do it1=1,nint(abs(t1min)/LatticeSpacings(0))
       call GaugeConf_t1%Update(sign(1._fp,t1min))
       t1 = t1 + sign(LatticeSpacings(0),t1min)

       !if(ThisProc()==0) print*,t1
    end do

    ! Counting number of (t,s) points --> Needed for output of progression
    nprogress = 0
    do it1=nint(t1min/LatticeSpacings(0)),nint(t1max/LatticeSpacings(0))
       t1 = it1*LatticeSpacings(0)

       it2range = nint(smax/LatticeSpacings(0))
       do it2=it1+1,it1+it2range
          nprogress = nprogress + 1
          
          t2 = it2*LatticeSpacings(0)
          t = (t2+t1)/2
          s =  t2-t1
          ! Check if t and s on the (t,s)-grid
          ittest = nint(t/ds)
          istest = nint(s/ds)

          dttest = ittest*ds
          dstest = istest*ds

          tTest = abs(dttest-t).le.GetZeroTol()
          sTest = abs(dstest-s).le.GetZeroTol()
          
          if(t > tmax) then
             exit
          end if
       end do
    end do

    ! Evoling along t1 ...
    iprogress = 0
    do it1=nint(t1min/LatticeSpacings(0)),nint(t1max/LatticeSpacings(0))
       t1 = it1*LatticeSpacings(0)
       
       ! ... and t2 --> getting (t,s)-coordinates
       GaugeConf_t2 = GaugeConf_t1
       HeavyField_t2 = HeavyField_t1
       it2range = nint(smax/LatticeSpacings(0))
       do it2=it1+1,it1+it2range
          iprogress = iprogress + 1
          if(ThisProc()==0) then
             write(output_unit,'(E16.9,A1)') real(iprogress)/nprogress*100,'%'
          end if
          
          call HeavyField_t2%Update(GaugeConf_t2,HeavyQuarkMass,WilsonCoefficients)
          call GaugeConf_t2%Update

          t2 = it2*LatticeSpacings(0)
          t = (t2+t1)/2
          s =  t2-t1

          ! Check if t and s on the (t,s)-grid
          ittest = nint(t/ds)
          istest = nint(s/ds)

          dttest = ittest*ds
          dstest = istest*ds

          tTest = abs(dttest-t).le.GetZeroTol()
          sTest = abs(dstest-s).le.GetZeroTol()
          
          if(tmin.le.t .and. t.le.tmax .and. tTest .and. sTest) then
             it = nint(t/ds)
             is = nint(s/ds)
             
             correlator(is,it) = HeavyField_t2%GetMesonCorrelator_O1_3s1_ZeroMomentum()
             qnorm(is,it) = HeavyField_t2%GetNorm_Quark()
             anorm(is,it) = HeavyField_t2%GetNorm_AntiQ()
             
             !if(ThisProc()==0) then
             !print*,real([t1,t2,t,s],real32),int([it,is],int16)
             !print*,int([it,is],int16)
             !end if
          else if(t > tmax) then
             exit
          end if
       end do

       call GaugeConf_t1%Update
    end do

1   continue
    if(ThisProc()==0) then
       ! Opening files
       fileID_Correlator = OpenFile(filename=FileMesonCorrelator,&
            st='REPLACE',fm='FORMATTED',act='WRITE')
       fileID_Norm = OpenFile(filename=FileNorm,&
            st='REPLACE',fm='FORMATTED',act='WRITE')


       ! Write header
       write(FileID_Correlator,'(A1,1X)',advance='no') 's'
       write(FileID_Norm,'(A1,1X)',advance='no') 's'
       do it=lbound(correlator,2),ubound(correlator,2)
          t = it*ds

          d = nint(-log10(ds))
          if(abs(t)>1) then
             p = nint(log10(abs(t)))
          else
             p = 0
          end if
          w = d + p + 1 + 1
          
          write(auxw,'(I1)') w
          write(auxd,'(I1)') d
          form = '(F'// auxw // '.' // auxd // ',1X)'
          
          write(FileID_Correlator,'(A2)',advance='no') 'Re'
          write(FileID_Correlator,form,advance='no')  t
          write(FileID_Correlator,'(A2)',advance='no') 'Im'
          write(FileID_Correlator,form,advance='no')  t
          
          write(FileID_Norm,'(A1)',advance='no') 'Q'
          write(FileID_Norm,form,advance='no')  t
          write(FileID_Norm,'(A1)',advance='no') 'A'
          write(FileID_Norm,form,advance='no')  t

          if(it==ubound(correlator,2)) then
             write(FileID_Correlator,*)
             write(FileID_Norm,*)            
          end if
          
       end do
       
       do is=lbound(correlator,1),ubound(correlator,1)
          s = is*ds
          write(FileID_Correlator,'(SP,E16.9,1X)',advance='no') s
          write(FileID_Norm,'(SP,E16.9,1X)',advance='no') s
          do it=lbound(correlator,2),ubound(correlator,2)
             if(it<ubound(correlator,2)) then
                write(FileID_Correlator,'(2(SP,E16.9,1X))',advance='no') &
                     real(correlator(is,it)),aimag(correlator(is,it))

                write(FileID_Norm,'(2(SP,E16.9,1X))',advance='no') &
                     qnorm(is,it),anorm(is,it)
             else
                write(FileID_Correlator,'((SP,E16.9,1X),(SP,E16.9))',advance='yes') &
                     real(correlator(is,it)),aimag(correlator(is,it))
                
                write(FileID_Norm,'(2(SP,E16.9,1X))',advance='yes') &
                     qnorm(is,it),anorm(is,it)
             end if
          end do
       end do
    
       ! Closing files
       call CloseFile(FileID_Correlator)
       call CloseFile(FileID_Norm)
    end if
    
    call EndSimulation
  contains
    !>@brief Initialisation of the simulation
    !!@details
    !! MPI\n
    !! Lattice-module\n
    !! Random number generator\n
    !! etc.
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 15.02.2019
    !!@version 1.0
    impure subroutine InitSimulation
      use precision, only: fp
      use, intrinsic :: iso_fortran_env

      use mpiinterface,       only: InitModule_MPIinterface       => InitModule, ThisProc, SyncAll
      use lattice,            only: InitModule_Lattice            => InitModule, nDim
      use halocomm,           only: InitModule_HaloComm           => InitModule
      use random,             only: InitModule_Random             => InitModule
      use xpfft,              only: InitModule_xpFFT              => InitModule
      use tolerances,         only: InitModule_tolerances         => InitModule
      use NRQCD,              only: InitModule_NRQCD              => InitModule
      implicit none

      integer(int64) :: arg_count
      character(len=80) :: arg
      integer(int8) :: i

      real(fp) :: c_re, c_im
      real(fp) :: kspTol
      
      !..--** Reading simulation parameters **--..
      arg_count = 1

      ! Spatial lattice parameters (extensions, spacings)
      do i=1,ndim
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(I4)') LatticeExtensions(i)
      end do

      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') ds
      LatticeSpacings(0) = ds/2
      
      do i=1,ndim
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(F10.13)') LatticeSpacings(i)
      end do


      ! Minimal CoM time
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') tmin
      ! Maximum CoM time
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') tmax
      ! Upper bound for s
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') smax
      !if(abs(smax).lt.LatticeSpacings(0)) smax=huge(smax)

      ! Seed for random number generator
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I4)') RandomNumberSeed

      ! Initial gluon distribution (box): Saturation scale
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') GluonSaturationScale

      ! Initial gluon distribution (box): Amplitude
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') GluonOccupationAmplitude

      ! Coupling (only relevant in initialisation)
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') GluonCoupling

      ! Tolerance for iterative PETSc solver
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(E15.7)') kspTol
      
      ! Heavy quark mass
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') HeavyQuarkmass
      
      ! Wilson coefficents
      do i=1,nWilsonCoefficients
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(F10.13)') c_re
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(F10.13)') c_im

         WilsonCoefficients(i) = cmplx(c_re,c_im,fp)
      end do
      
      arg_count = arg_count +1; call get_command_argument(arg_count,FileMesonCorrelator);
      arg_count = arg_count +1; call get_command_argument(arg_count,FileNorm);

      !..--** Module initialisations **--..
      call InitModule_MPIinterface
      call InitModule_Lattice(LatticeExtensions(1:ndim),LatticeSpacings(0:ndim))
      call InitModule_HaloComm
      call InitModule_xpFFT
      call InitModule_Random(RandomNumberSeed + ThisProc())
      call InitModule_tolerances
      call InitModule_NRQCD(HeavyQuarkMass,WilsonCoefficients,1._fp,kspTol)

      call SyncAll
    end subroutine InitSimulation

    !>@brief Ending of the simulation
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 15.02.2019
    !!@version 1.0
    subroutine EndSimulation
      use mpiinterface,   only: ThisProc, FinalizeModule_MPIinterface   => FinalizeModule
      use xpfft,          only: FinalizeModule_xpFFT          => FinalizeModule
      use NRQCD,          only: FinalizeModule_NRQCD          => FinalizeModule
      implicit none

      call FinalizeModule_NRQCD
      call FinalizeModule_xpFFT
      call FinalizeModule_MPIinterface
      
      if(ThisProc()==0) write(output_unit,*) "Simulation completed"
      STOP
    end subroutine EndSimulation
  end subroutine MeasureHeavyQuarkoniumCorrelators
  
  !>@brief Program for measuring the time-evolution of the energy and gauss law deviation
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 07.03.2019
  !!@version 1.0
  impure subroutine MeasureEnergyAndGaussLawDeviation
    use, intrinsic :: iso_fortran_env
    use precision
    use mpiinterface

    use lattice
    use gaugeconfiguration_su3
    use mpi
    use io

    implicit none

    ! Simulation parameters
    integer(int64) :: LatticeExtensions(ndim)
    real(fp)       :: LatticeSpacings(0:ndim)
    integer(int64) :: TimeSteps
    integer(int64) :: RandomNumberSeed

    real(fp)   :: GluonSaturationScale !qs
    real(fp)   :: GluonOccupationAmplitude ! Amplitude of box in units of 1/g^2
    real(fp)   :: GluonCoupling
    real(fp)   :: CoMTime
    real(fp)   :: TimeRange

    ! Physical fields
    type(GaugeConfiguration) :: GaugeConf

    ! Monitoring variables
    integer(int64) :: it
    real(fp) :: time
    
    ! Observable output (Gauss, Energy)
    integer(int8) :: fileID_eg
    real(fp) :: gauss, energy

    call InitSimulation

    call GaugeConf%TransversePolarisedOccupiedInit_Box(&
         GluonSaturationScale,GluonOccupationAmplitude,GluonCoupling)
    if(ThisProc()==0) &
         fileID_eg = OpenFile(filename="energy_gauss.txt",st='REPLACE',fm='FORMATTED',act='WRITE')
    do it=0,TimeSteps
       time = it*GetLatticeSpacing(0_int8)
       gauss = GaugeConf%GetDeviationFromGaussLaw()
       energy= GaugeConf%GetEnergy()
       if(ThisProc()==0) write(fileID_eg,'(3(SP,E13.6,1X))') time,energy,gauss
       if(ThisProc()==0) write(output_unit,*) time
       call GaugeConf%Update
    end do
    if(ThisProc()==0) call CloseFile(fileID_eg)

    call EndSimulation
  contains

    !>@brief Initialisation of the simulation
    !!@details
    !! MPI\n
    !! Lattice-module\n
    !! Random number generator\n
    !! etc.
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 15.02.2019
    !!@version 1.0
    impure subroutine InitSimulation
      use precision, only: fp
      use, intrinsic :: iso_fortran_env

      use mpiinterface,       only: InitModule_MPIinterface       => InitModule, ThisProc, SyncAll
      use lattice,            only: InitModule_Lattice            => InitModule, nDim
      use halocomm,           only: InitModule_HaloComm           => InitModule
      use random,             only: InitModule_Random             => InitModule
      use xpfft,              only: InitModule_xpFFT              => InitModule
      use tolerances,         only: InitModule_tolerances         => InitModule
      implicit none

      integer(int64) :: arg_count
      character(len=80) :: arg
      integer(int8) :: i

      !..--** Reading simulation parameters **--..
      arg_count = 1

      ! Spatial lattice parameters (extensions, spacings)
      do i=1,ndim
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(I4)') LatticeExtensions(i)
      end do

      do i=0,ndim
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(F10.13)') LatticeSpacings(i)
      end do

      ! Center of mass time T
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') CoMTime
      ! Center of mass time-range smax
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') TimeRange

      TimeSteps=ceiling(TimeRange/LatticeSpacings(0))

      ! Seed for random number generator
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I4)') RandomNumberSeed

      ! Initial gluon distribution (box): Saturation scale
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') GluonSaturationScale

      ! Initial gluon distribution (box): Amplitude
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') GluonOccupationAmplitude

      ! Coupling (only relevant in initialisation)
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') GluonCoupling

      !..--** Module initialisations **--..
      call InitModule_MPIinterface
      call InitModule_Lattice(LatticeExtensions(1:ndim),LatticeSpacings(0:ndim))
      call InitModule_HaloComm
      call InitModule_xpFFT
      call InitModule_Random(RandomNumberSeed + ThisProc())
      call InitModule_tolerances

      call SyncAll
    end subroutine InitSimulation

    !>@brief Ending of the simulation
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 15.02.2019
    !!@version 1.0
    subroutine EndSimulation
      use mpiinterface, only: ThisProc, FinalizeModule_MPIinterface => FinalizeModule
      use xpfft,        only: FinalizeModule_xpFFT        => FinalizeModule
      implicit none

      call FinalizeModule_xpFFT
      call FinalizeModule_MPIinterface

      if(ThisProc()==0) write(output_unit,*) "Simulation completed"
      STOP
    end subroutine EndSimulation
  end subroutine MeasureEnergyAndGaussLawDeviation

  !>@brief Program for measuring the gluon distribution function
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 04.06.2019
  !!@version 1.0
  impure subroutine MeasureGluondistribution_Equilibrium
    use, intrinsic :: iso_fortran_env
    use precision
    use mpiinterface

    use lattice
    use gaugeconfiguration_su3
    use mpi
    use io
    implicit none

    ! Simulation parameters
    integer(int64) :: LatticeExtensions(ndim)
    real(fp)       :: LatticeSpacings(0:ndim)
    integer(int64) :: RandomNumberSeed

    ! Physical fields
    type(GaugeConfiguration) :: GaugeConf

    real(fp) :: TimeRange

    real(fp) :: Beta
    integer(int64) :: nefieldinit,nequilibriumtimesteps
    logical :: MeasureThermilisationEnergy
    character(len=80) :: Filename_ThermilisationEnergy

    character(len=80) :: Filename_Gluondist
    
    integer(int64) :: GaugeFixingMaxIterations
    real(fp) :: GaugefixingCoefficient, GaugeFixingTolerance
    integer(int64) :: Number_of_Measurements_of_Gluondistribution
    integer(int64) :: TimePoints_between_Measurement_of_Gluondistribution
    real(fp) :: TimeBetweenGluonMeasurements

    integer(int64) :: TimeSteps

    call InitSimulation
    
    call Gaugeconf%EquilibriumInit(Beta,nefieldinit,nequilibriumtimesteps,MeasureEnergy=MeasureThermilisationEnergy,FileName=Filename_ThermilisationEnergy)

    call PrintObservables(GaugeConf=GaugeConf)

    call EndSimulation

  contains
    impure subroutine PrintObservables(GaugeConf)
      use lattice
      implicit none
      
      type(GaugeConfiguration), intent(in) :: GaugeConf

      type(GaugeConfiguration) :: GaugeConf_gaugefixed

      real(fp), allocatable :: aa_correlator(:),ee_correlator(:), gluondist(:)

      real(fp), allocatable :: momenta(:)
      
      integer(int64) :: MemoryIndex, LatticeIndex, i, is

      ! MPI
      integer(intmpi) :: mpierr, mpistatus(mpi_status_size), tag, src, buffersize
      
      ! Output
      integer(int8) :: FileID
      character(len=128) :: filename,time_tag
      integer(intmpi), allocatable :: mpisendrequest(:), mpisendstatus(:,:)

      integer(intmpi) :: proc

      integer(intmpi), parameter :: dest=0_intmpi
      
      GaugeConf_gaugefixed = GaugeConf
      
      call gaugeconf_gaugefixed%CoulombGaugefixing(&
           Tolerance_    =GaugeFixingTolerance,&
           alpha_        =GaugefixingCoefficient,&
           MaxIterations_=GaugefixingMaxIterations)

      
      call gaugeconf_gaugefixed%GetTransverseAACorrelator(aa_correlator)
      call gaugeconf_gaugefixed%GetTransverseEECorrelator(ee_correlator)

      allocate(gluondist(size(aa_correlator)))
      allocate(Momenta(GetLocalLatticeSize()))
      allocate(mpisendrequest(2))
      
      if(ThisProc()==dest) then
          fileID = OpenFile(filename=Filename_Gluondist,st='REPLACE',fm='FORMATTED',act='WRITE')
       end if
       do src=0,NumProcs()-1
          if(ThisProc()==dest .or. ThisProc()==src) then
             
             ! Momenta
             buffersize = size(Momenta)
             tag = (0+NumProcs())*src
             if(ThisProc()==src) then
                is = 0
                do MemoryIndex=1,GetMemorySize()
                   LatticeIndex = GetLatticeIndex_M(MemoryIndex)
                   if(ThisProc()==GetProc_G(LatticeIndex)) then
                      is = is + 1
                      Momenta(is) = GetNorm2Momentum_G(LatticeIndex)
                   end if
                end do
                call mpi_isend(&
                     Momenta,                      & ! What to send
                     buffersize,                   & ! How many points
                     MPI_DOUBLE,                   & ! What type
                     dest,                         & ! Recieving process
                     tag,                          & ! Tag
                     MPI_COMM_WORLD,               & ! Communicator
                     mpisendrequest(1),            & ! Request handle
                     mpierr)                         ! Error code
             end if
             if(ThisProc()==dest) then
                call mpi_recv(&
                     Momenta,                      & ! What to recieve
                     buffersize,                   & ! How many points
                     MPI_DOUBLE,                   & ! What type
                     src,                          & ! Sending process
                     tag,                          & ! Tag
                     MPI_COMM_WORLD,               & ! Communicator
                     mpistatus,                    & ! Status
                     mpierr)                         ! Error code
             end if

             ! Gluon distribution
             buffersize = size(gluondist)
             tag = (1+NumProcs())*src
             if(src .ne. dest) then
                if(ThisProc()==src) then
                   gluondist = sqrt(aa_correlator*ee_correlator)
                   call mpi_isend(&
                        GluonDist,& ! What to send
                        buffersize,                   & ! How many points
                        MPI_DOUBLE,                   & ! What type
                        dest,                         & ! Recieving process
                        tag,                          & ! Tag
                        MPI_COMM_WORLD,               & ! Communicator
                        mpisendrequest(2),            & ! Request handle
                        mpierr)                         ! Error code
                end if
                if(ThisProc()==dest) then
                   call mpi_recv(&
                        GluonDist,& ! What to recieve
                        buffersize,                   & ! How many points
                        MPI_DOUBLE,                   & ! What type
                        src,                          & ! Sending process
                        tag,                          & ! Tag
                        MPI_COMM_WORLD,               & ! Communicator
                        mpistatus,                    & ! Status
                        mpierr)                         ! Error code
                end if
             else
                gluondist = sqrt(aa_correlator*ee_correlator)
             end if

             if(ThisProc()==dest) then
                do is=1,size(Momenta)
                   write(fileID,'(2(SP,E13.6,1X))') &
                        Momenta(is),&
                        GluonDist(is)
                end do
             end if
          end if

          call SyncAll
       end do
       if(ThisProc()==dest) call CloseFile(FileID)

      
    end subroutine PrintObservables
    
    !>@brief Initialisation of the simulation
    !!@details
    !! MPI\n
    !! Lattice-module\n
    !! Random number generator\n
    !! etc.
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 15.02.2019
    !!@version 1.0
    impure subroutine InitSimulation
      use precision, only: fp
      use, intrinsic :: iso_fortran_env

      use mpiinterface,       only: InitModule_MPIinterface       => InitModule, ThisProc, SyncAll
      use lattice,            only: InitModule_Lattice            => InitModule, nDim
      use halocomm,           only: InitModule_HaloComm           => InitModule
      use random,             only: InitModule_Random             => InitModule
      use xpfft,              only: InitModule_xpFFT              => InitModule
      use tolerances,         only: InitModule_tolerances         => InitModule
      implicit none

      integer(int64) :: arg_count
      character(len=80) :: arg
      integer(int8) :: i, j

      !..--** Reading simulation parameters **--..
      arg_count = 1

      ! Spatial lattice parameters (extensions, spacings)
      do i=1,ndim
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(I4)') LatticeExtensions(i)
      end do

      do i=0,ndim
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(F10.13)') LatticeSpacings(i)
      end do

      ! Center of mass time-range smax
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') TimeRange

      TimeSteps=ceiling(TimeRange/LatticeSpacings(0))

      ! Seed for random number generator
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I4)') RandomNumberSeed

      ! Initial gluon distribution (box): Saturation scale
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') beta

      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I6)') nefieldinit
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I6)') nequilibriumtimesteps      

      ! Time between gluon measurements
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') TimeBetweenGluonMeasurements
      Number_of_Measurements_of_Gluondistribution = aint(TimeRange/TimeBetweenGluonMeasurements,int64)
      TimePoints_between_Measurement_of_Gluondistribution&
           = int(TimeBetweenGluonMeasurements/LatticeSpacings(0))

      ! Maximum number of iterations for coulomb gauge fixing
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I7)') GaugefixingMaxIterations
      ! Tolerance for coulomb gauge fixing
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(E15.7)') GaugefixingTolerance
      ! Tolerance for coulomb gauge fixing
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') GaugefixingCoefficient

      arg_count = arg_count +1; call get_command_argument(arg_count,Filename_Gluondist);
      
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I1)') j
      if(j==0) then
         MeasureThermilisationEnergy=.false.
      else
         MeasureThermilisationEnergy=.true.
         arg_count = arg_count +1; call get_command_argument(arg_count,Filename_ThermilisationEnergy);
      end if
      
      !..--** Module initialisations **--..
      call InitModule_MPIinterface
      call InitModule_Lattice(LatticeExtensions(1:ndim),LatticeSpacings(0:ndim))
      call InitModule_HaloComm
      call InitModule_xpFFT
      call InitModule_Random(RandomNumberSeed + ThisProc())
      call InitModule_tolerances

      call SyncAll
    end subroutine InitSimulation

    !>@brief Ending of the simulation
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 15.02.2019
    !!@version 1.0
    subroutine EndSimulation
      use mpiinterface, only: ThisProc, FinalizeModule_MPIinterface => FinalizeModule
      use xpfft,        only: FinalizeModule_xpFFT        => FinalizeModule
      implicit none
      
      call FinalizeModule_xpFFT
      call FinalizeModule_MPIinterface

      if(ThisProc()==0) write(output_unit,*) "Simulation completed"
      STOP
    end subroutine EndSimulation
  end subroutine MeasureGluondistribution_Equilibrium

  
  !>@brief Program for measuring the gluon distribution function
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 07.03.2019
  !!@version 1.0
  impure subroutine MeasureGluondistribution
    use, intrinsic :: iso_fortran_env
    use precision
    use mpiinterface

    use lattice
    use gaugeconfiguration_su3
    use mpi
    use io

    implicit none

    ! Simulation parameters
    integer(int64) :: LatticeExtensions(ndim)
    real(fp)       :: LatticeSpacings(0:ndim)
    integer(int64) :: RandomNumberSeed

    real(fp)   :: GluonSaturationScale !qs
    real(fp)   :: GluonOccupationAmplitude ! Amplitude of box in units of 1/g^2
    real(fp)   :: GluonCoupling
    real(fp)   :: CoMTime

    character(len=80) :: Filename_Gluondist
    
    ! Physical fields
    type(GaugeConfiguration) :: GaugeConf

    integer(int64) :: GaugeFixingMaxIterations
    real(fp) :: GaugefixingCoefficient, GaugeFixingTolerance

    integer(int64) :: TimeSteps, it

    call InitSimulation

    call GaugeConf%TransversePolarisedOccupiedInit_Box(&
         GluonSaturationScale,GluonOccupationAmplitude,GluonCoupling)

    do it=1,TimeSteps
       call GaugeConf%Update(sign(+1._fp,CoMTime))
    end do

    ! Measuring
    call PrintObservables(GaugeConf=GaugeConf)
    
    call EndSimulation
  contains

    impure subroutine PrintObservables(GaugeConf)
      use lattice
      implicit none
      
      type(GaugeConfiguration), intent(in) :: GaugeConf

      type(GaugeConfiguration) :: GaugeConf_gaugefixed

      real(fp), allocatable :: aa_correlator(:),ee_correlator(:), gluondist(:)

      real(fp), allocatable :: momenta(:)
      
      integer(int64) :: MemoryIndex, LatticeIndex, i, is

      ! MPI
      integer(intmpi) :: mpierr, mpistatus(mpi_status_size), tag, src, buffersize
      
      ! Output
      integer(int8) :: FileID
      integer(intmpi), allocatable :: mpisendrequest(:), mpisendstatus(:,:)

      integer(intmpi) :: proc

      integer(intmpi), parameter :: dest=0_intmpi


      GaugeConf_gaugefixed = GaugeConf
      
      call gaugeconf_gaugefixed%CoulombGaugefixing(&
           Tolerance_    =GaugeFixingTolerance,&
           alpha_        =GaugefixingCoefficient,&
           MaxIterations_=GaugefixingMaxIterations)

      
      call gaugeconf_gaugefixed%GetTransverseAACorrelator(aa_correlator)
      call gaugeconf_gaugefixed%GetTransverseEECorrelator(ee_correlator)

      allocate(gluondist(size(aa_correlator)))
      allocate(Momenta(GetLocalLatticeSize()))
      allocate(mpisendrequest(2))
      
      if(ThisProc()==dest) then
          fileID = OpenFile(filename=Filename_Gluondist,st='REPLACE',fm='FORMATTED',act='WRITE')
       end if
       do src=0,NumProcs()-1
          if(ThisProc()==dest .or. ThisProc()==src) then
             
             ! Momenta
             buffersize = size(Momenta)
             tag = (0+NumProcs())*src
             if(ThisProc()==src) then
                is = 0
                do MemoryIndex=1,GetMemorySize()
                   LatticeIndex = GetLatticeIndex_M(MemoryIndex)
                   if(ThisProc()==GetProc_G(LatticeIndex)) then
                      is = is + 1
                      Momenta(is) = GetNorm2Momentum_G(LatticeIndex)
                   end if
                end do
                call mpi_isend(&
                     Momenta,                      & ! What to send
                     buffersize,                   & ! How many points
                     MPI_DOUBLE,                   & ! What type
                     dest,                         & ! Recieving process
                     tag,                          & ! Tag
                     MPI_COMM_WORLD,               & ! Communicator
                     mpisendrequest(1),            & ! Request handle
                     mpierr)                         ! Error code
             end if
             if(ThisProc()==dest) then
                call mpi_recv(&
                     Momenta,                      & ! What to recieve
                     buffersize,                   & ! How many points
                     MPI_DOUBLE,                   & ! What type
                     src,                          & ! Sending process
                     tag,                          & ! Tag
                     MPI_COMM_WORLD,               & ! Communicator
                     mpistatus,                    & ! Status
                     mpierr)                         ! Error code
             end if

             ! Gluon distribution
             buffersize = size(gluondist)
             tag = (1+NumProcs())*src
             if(src .ne. dest) then
                if(ThisProc()==src) then
                   gluondist = sqrt(aa_correlator*ee_correlator)
                   call mpi_isend(&
                        GluonDist,& ! What to send
                        buffersize,                   & ! How many points
                        MPI_DOUBLE,                   & ! What type
                        dest,                         & ! Recieving process
                        tag,                          & ! Tag
                        MPI_COMM_WORLD,               & ! Communicator
                        mpisendrequest(2),            & ! Request handle
                        mpierr)                         ! Error code
                end if
                if(ThisProc()==dest) then
                   call mpi_recv(&
                        GluonDist,& ! What to recieve
                        buffersize,                   & ! How many points
                        MPI_DOUBLE,                   & ! What type
                        src,                          & ! Sending process
                        tag,                          & ! Tag
                        MPI_COMM_WORLD,               & ! Communicator
                        mpistatus,                    & ! Status
                        mpierr)                         ! Error code
                end if
             else
                gluondist = sqrt(aa_correlator*ee_correlator)
             end if

             if(ThisProc()==dest) then
                do is=1,size(Momenta)
                   write(fileID,'(2(SP,E13.6,1X))') &
                        Momenta(is),&
                        GluonDist(is)
                end do
             end if
          end if

          call SyncAll
       end do
       if(ThisProc()==dest) call CloseFile(FileID)
      
    end subroutine PrintObservables

    !>@brief Initialisation of the simulation
    !!@details
    !! MPI\n
    !! Lattice-module\n
    !! Random number generator\n
    !! etc.
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 15.02.2019
    !!@version 1.0
    impure subroutine InitSimulation
      use precision, only: fp
      use, intrinsic :: iso_fortran_env

      use mpiinterface,       only: InitModule_MPIinterface       => InitModule, ThisProc, SyncAll
      use lattice,            only: InitModule_Lattice            => InitModule, nDim
      use halocomm,           only: InitModule_HaloComm           => InitModule
      use random,             only: InitModule_Random             => InitModule
      use xpfft,              only: InitModule_xpFFT              => InitModule
      use tolerances,         only: InitModule_tolerances         => InitModule
      implicit none

      integer(int64) :: arg_count
      character(len=80) :: arg
      integer(int8) :: i

      !..--** Reading simulation parameters **--..
      arg_count = 1

      ! Spatial lattice parameters (extensions, spacings)
      do i=1,ndim
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(I4)') LatticeExtensions(i)
      end do

      do i=0,ndim
         arg_count = arg_count +1; call get_command_argument(arg_count,arg);
         read(arg,'(F10.13)') LatticeSpacings(i)
      end do

      ! Center of mass time T
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') CoMTime

      TimeSteps=ceiling(abs(CoMTime)/LatticeSpacings(0))

      ! Seed for random number generator
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I4)') RandomNumberSeed

      ! Initial gluon distribution (box): Saturation scale
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') GluonSaturationScale

      ! Initial gluon distribution (box): Amplitude
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') GluonOccupationAmplitude

      ! Coupling (only relevant in initialisation)
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') GluonCoupling

      ! Maximum number of iterations for coulomb gauge fixing
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I7)') GaugefixingMaxIterations
      ! Tolerance for coulomb gauge fixing
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(E15.7)') GaugefixingTolerance
      ! Tolerance for coulomb gauge fixing
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') GaugefixingCoefficient

      ! Name of output file
      arg_count = arg_count +1; call get_command_argument(arg_count,Filename_Gluondist);
      
      !..--** Module initialisations **--..
      call InitModule_MPIinterface
      call InitModule_Lattice(LatticeExtensions(1:ndim),LatticeSpacings(0:ndim))
      call InitModule_HaloComm
      call InitModule_xpFFT
      call InitModule_Random(RandomNumberSeed + ThisProc())
      call InitModule_tolerances

      call SyncAll
    end subroutine InitSimulation

    !>@brief Ending of the simulation
    !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
    !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
    !!@date 15.02.2019
    !!@version 1.0
    subroutine EndSimulation
      use mpiinterface, only: ThisProc, FinalizeModule_MPIinterface => FinalizeModule
      use xpfft,        only: FinalizeModule_xpFFT        => FinalizeModule
      implicit none
      
      call FinalizeModule_xpFFT
      call FinalizeModule_MPIinterface

      if(ThisProc()==0) write(output_unit,*) "Simulation completed"
      STOP
    end subroutine EndSimulation
  end subroutine MeasureGluondistribution
end module programs



program simulation
  use, intrinsic :: iso_fortran_env
  use mpiinterface, only: MPIstop
  use programs
  implicit none

  integer(int8) :: mode

  character(len=3) :: arg

  call get_command_argument(1,arg);
  read(arg,'(I3)') mode

  select case(mode)
  case(1)
     call MeasureGluondistribution
  case(11)
     call MeasureGluondistribution_Equilibrium
  case(2)
     call MeasureEnergyAndGaussLawDeviation
  case(3)
     call MeasureHeavyQuarkoniumCorrelators
  case (31)
     call MeasureHeavyQuarkoniumCorrelators_oneT
  case (32)
     call MeasureHeavyQuarkoniumCorrelators_free
  case (33)
     call MeasureHeavyQuarkoniumCorrelators_Equilibrium
  case(4)
     call ComputeSpectrum
  case (5)
     call MeasureWilsonLines_nonEquilibrium
  case (51)
     call MeasureWilsonAndHybridLines_Equilibrium
  case (52)
     call MeasureWilsonLines_Equilibrium_StaticMeson
  case default
     call MPIStop('Invalid simulation mode selected')
  end select
end program simulation
