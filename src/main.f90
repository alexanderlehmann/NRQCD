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
         write(FileID_Norm,'(3(SP,E16.9))') s,norm_quark,norm_antiq


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
    integer :: i, quarkskips, UpdateQuarksEveryNsteps
    
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
       
       quarkskips = 0
       do it2=1,t2steps
          iwork = iwork + 1
          if(modulo(quarkskips,UpdateQuarksEveryNsteps)==0) then
             if(t2steps-it2+1.lt.UpdateQuarksEveryNsteps) then
                call HeavyField_t2%Update(GaugeConf_t2,HeavyQuarkMass,WilsonCoefficients,&
                   real(t2steps-it2+1,fp))
             else
                call HeavyField_t2%Update(GaugeConf_t2,HeavyQuarkMass,WilsonCoefficients,&
                   real(UpdateQuarksEveryNsteps,fp))
             end if
             quarkskips = 0
             if(ThisProc()==0) write(output_unit,'(F7.3,A1)') real(iwork)/nwork*100,'%'
          end if
          call GaugeConf_t2%Update
          quarkskips = quarkskips + 1
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
         write(FileID_Norm,'(3(SP,E16.9))') s,norm_quark,norm_antiq

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

      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I3)') UpdateQuarksEveryNsteps
      
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
    
    ! Physical fields
    type(GaugeConfiguration) :: GaugeConf_t1, GaugeConf_t2, GaugeConf_initial
    type(NRQCDField)         :: HeavyField
    
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
    complex(fp) :: WilsonLoop, HybridLoop

    complex(fp), allocatable :: WilsonLoops(:,:,:), HybridLoops(:,:,:)
    real(fp), allocatable :: rObservable(:), iObservable(:)
    real(fp) :: rMean, rStderr,iMean,iStderr
    real(fp) :: time
    integer(int64) :: TimePoints
    ! Output
    integer(int8) :: FileID_WilsonLoops, FileID_HybridLoops

    character(len=80) :: &
         FileName_WilsonLoops, FileName_HybridLoops

    integer(intmpi) :: proc

    integer(int64) :: measurement
    integer(int64) :: nMeasurement

    call InitSimulation

    rmax = LatticeExtensions(messdir)/2

    TimePoints = nint(TimeRange/LatticeSpacings(0))
    
    allocate(WilsonLoops(nMeasurement,0:rmax,0:TimePoints))
    WilsonLoops = 0
    allocate(HybridLoops(nMeasurement,0:rmax,0:TimePoints))
    HybridLoops = 0

    do measurement=1,nMeasurement

       if(ThisProc()==0) write(output_unit,*)&
            int(measurement,int16),'of',&
            int(nMeasurement,int16),'configurations';&
            call flush(output_unit)

       ! initialisation of config ....
       call GaugeConf_initial%EquilibriumInit(Beta)
       if(thisproc()==0) write(output_unit,*) 'Done: Equilibration'
       GaugeConf_t1 = GaugeConf_initial
       
       do it=1,abs(NINT(tstart/LatticeSpacings(0)))
          call GaugeConf_t1%Update(sign(+1._real64,tstart))
       end do

       
       do r=0,rmax
          if(ThisProc()==0) then
             write(output_unit,*) 'r=',r,'of',rmax
          end if
          
          ! Initialising quark-antiquark-pair
          ! with quark at variable remote point xr
          ! and antiquark at fixed (origin) point x0
          meanMaxLatticeIndex=1
          do x0=1,meanMaxLatticeIndex
             ! Getting xr
             xr = x0
             do shiftstep=1,r
                xr=GetNeib_G(messdir,xr)
             end do
             
             GaugeConf_t2 = GaugeConf_t1
             call HeavyField%InitSinglePoint(&
                  latticeindex_quark=xr,&
                  latticeindex_antiq=x0)
             do it=0,nint(TimeRange/LatticeSpacings(0)),+1
                if(ThisProc()==0) write(output_unit,'(I5,1X,A2,1X,I5)')it,'of',nint(TimeRange/LatticeSpacings(0))
                if(x0==1) then
                   WilsonLoop = GetWilsonLoop(GaugeConf_t1,GaugeConf_t2,r,messdir)
                   WilsonLoops(measurement,r,it) = WilsonLoop
                end if
                
                HybridLoop = GetHybridLoop(GaugeConf_t1,GaugeConf_t2,HeavyField,x0,r,messdir)
                HybridLoops(measurement,r,it) = HybridLoops(measurement,r,it) &
                     + HybridLoop/meanMaxLatticeIndex

                call HeavyField%Update(GaugeConf_t2,HeavyQuarkMass,WilsonCoefficients)
                call GaugeConf_t2%Update
             end do
          end do
       end do
    end do
1   continue
    if(ThisProc()==0) then
       fileID_WilsonLoops = OpenFile(filename=FileName_WilsonLoops,&
            st='REPLACE',fm='FORMATTED',act='WRITE')
       fileID_HybridLoops = OpenFile(filename=FileName_HybridLoops,&
            st='REPLACE',fm='FORMATTED',act='WRITE')

       allocate(rObservable(size(WilsonLoops,1)))
       allocate(iObservable(size(WilsonLoops,1)))

       do it=0,nint(TimeRange/LatticeSpacings(0)),+1

          time = tstart + it*LatticeSpacings(0)
          do r=0,rmax
             if(r==0) then
                write(FileID_WilsonLoops,'(1(SP,E16.9,1X))',advance='no') time
             end if
             rObservable = real(WilsonLoops(:,r,it))
             iObservable = aimag(WilsonLoops(:,r,it))

             rMean = GetMean(rObservable)
             iMean = GetMean(iObservable)
             rStdErr = GetStdError(rObservable)
             iStdErr = GetStdError(iObservable)

             !WilsonLoop = cmplx(rMean,iMean)
             
             if(r<rmax) then
                write(FileID_WilsonLoops,'(4(SP,E16.9,1X))',advance='no') &
                     rMean,rStdErr,iMean,iStderr
             else
                write(FileID_WilsonLoops,'(4(SP,E16.9,1X))',advance='yes') &
                     rMean,rStdErr,iMean,iStderr
             end if


             
             if(r==0) then
                write(FileID_HybridLoops,'(1(SP,E16.9,1X))',advance='no') time
             end if
             rObservable = real(HybridLoops(:,r,it))
             iObservable = aimag(HybridLoops(:,r,it))

             rMean = GetMean(rObservable)
             iMean = GetMean(iObservable)
             rStdErr = GetStdError(rObservable)
             iStdErr = GetStdError(iObservable)

             !HybridLoop = cmplx(rMean,iMean)
             
             if(r<rmax) then
                write(FileID_HybridLoops,'(4(SP,E16.9,1X))',advance='no') &
                     rMean,rStdErr,iMean,iStderr
             else
                write(FileID_HybridLoops,'(4(SP,E16.9,1X))',advance='yes') &
                     rMean,rStdErr,iMean,iStderr
             end if
          end do
       end do

       call CloseFile(FileID_WilsonLoops)
       call CloseFile(FileID_HybridLoops)
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

      ! Seed for random number generator
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I4)') RandomNumberSeed
      ! Initial gluon distribution (box): Saturation scale
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') Beta

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
      
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I4)') nMeasurement
      
      ! Output filenames
      arg_count = arg_count +1; call get_command_argument(arg_count,FileName_WilsonLoops);
      arg_count = arg_count +1; call get_command_argument(arg_count,FileName_HybridLoops);

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
  
  impure subroutine MeasureWilsonLines_Equilibrium

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

    ! Physical fields
    type(GaugeConfiguration) :: GaugeConf, GaugeConf_initial

    ! Counting
    integer :: i

    integer(int64) :: it

    ! Wilson loop parameters
    integer(int8), parameter :: messdir=nDim
    integer(int64) :: rmax, r

    ! Observables
    complex(fp) :: WilsonLoop, TimeDerivativeWilsonLoop, Potential

    complex(fp), allocatable :: WilsonLoops(:,:,:), TimeDerivativeWilsonLoops(:,:,:)
    real(fp), allocatable :: rObservable(:), iObservable(:)
    real(fp) :: rMean, rStderr,iMean,iStderr
    real(fp) :: time
    integer(int64) :: TimePoints
    ! Output
    integer(int8) :: FileID_WilsonLoops!, FileID_TimeDerivativeWilsonLoops, FileID_Potential

    character(len=80) :: &
         FileName_WilsonLoops!, FileName_TimeDerivativeWilsonLoops, FileName_Potential

    integer(intmpi) :: proc

    integer(int64) :: measurement
    integer(int64) :: nMeasurement

    call InitSimulation

    rmax = LatticeExtensions(messdir)/2

    TimePoints = nint(TimeRange/LatticeSpacings(0))
    
    allocate(WilsonLoops(nMeasurement,rmax,0:TimePoints))
    WilsonLoops = 0
    !allocate(TimeDerivativeWilsonLoops(nMeasurement,rmax,0:TimePoints))
    !TimeDerivativeWilsonLoops = 0

    do measurement=1,nMeasurement

       if(ThisProc()==0) write(output_unit,*)&
            int(measurement,int16),'of',&
            int(nMeasurement,int16),'configurations';&
            call flush(output_unit)

       ! initialisation of config ....
       call GaugeConf_initial%EquilibriumInit(Beta)
       
       if(thisproc()==0) write(output_unit,*) 'Done: Equilibration'
       GaugeConf = GaugeConf_initial

       do it=1,abs(NINT(tstart/LatticeSpacings(0)))
          call GaugeConf%Update(sign(+1._real64,tstart))
       end do

       do it=0,nint(TimeRange/LatticeSpacings(0)),+1
          if(thisproc()==0) write(output_unit,*) it
          do r=1,rmax
             WilsonLoop = GetWilsonLoop(GaugeConf_initial,GaugeConf,r,messdir)
             !TimeDerivativeWilsonLoop = &
             !     GetTimeDerivativeWilsonLoop(GaugeConf_initial,GaugeConf,r,messdir)
             !GetPotentialWilsonLoop(GaugeConf_initial,GaugeConf,r,messdir)

             WilsonLoops(measurement,r,it) = WilsonLoop
             !TimeDerivativeWilsonLoops(measurement,r,it) = TimeDerivativeWilsonLoop
          end do

          call GaugeConf%Update
       end do
    end do
    if(ThisProc()==0) then
       fileID_WilsonLoops = OpenFile(filename=FileName_WilsonLoops,&
            st='REPLACE',fm='FORMATTED',act='WRITE')
       !fileID_TimeDerivativeWilsonLoops = OpenFile(filename=FileName_TimeDerivativeWilsonLoops,&
       !     st='REPLACE',fm='FORMATTED',act='WRITE')
       !fileID_Potential = OpenFile(filename=FileName_Potential,&
       !     st='REPLACE',fm='FORMATTED',act='WRITE')
       do it=0,nint(TimeRange/LatticeSpacings(0)),+1

          time = tstart + it*LatticeSpacings(0)
          if(ThisProc()==0) then
             write(output_unit,*) real(time,real32)
             call flush(output_unit)
          end if
          do r=1,rmax
             if(r==1) then
                write(FileID_WilsonLoops,'(1(SP,E16.9,1X))',advance='no') time
             end if
             rObservable = real(WilsonLoops(:,r,it))
             iObservable = aimag(WilsonLoops(:,r,it))

             rMean = GetMean(rObservable)
             iMean = GetMean(iObservable)
             rStdErr = GetStdError(rObservable)
             iStdErr = GetStdError(iObservable)

             WilsonLoop = cmplx(rMean,iMean)
             
             if(r<rmax) then
                write(FileID_WilsonLoops,'(4(SP,E16.9,1X))',advance='no') &
                     rMean,rStdErr,iMean,iStderr
             else
                write(FileID_WilsonLoops,'(4(SP,E16.9,1X))',advance='yes') &
                     rMean,rStdErr,iMean,iStderr
             end if

             
             !if(r==1) then
             !   write(FileID_TimeDerivativeWilsonLoops,'(1(SP,E16.9,1X))',advance='no') time
             !end if
             !rObservable = real(TimeDerivativeWilsonLoops(:,r,it))
             !iObservable = aimag(TimeDerivativeWilsonLoops(:,r,it))

             !rMean = GetMean(rObservable)
             !iMean = GetMean(iObservable)
             !rStdErr = GetStdError(rObservable)
             !iStdErr = GetStdError(iObservable)

             !TimeDerivativeWilsonLoop = cmplx(rMean,iMean)
             
             !if(r<rmax) then
             !   write(FileID_TimeDerivativeWilsonLoops,'(4(SP,E16.9,1X))',advance='no') &
             !        rMean,rStdErr,iMean,iStderr
             !else
             !   write(FileID_TimeDerivativeWilsonLoops,'(4(SP,E16.9,1X))',advance='yes') &
             !        rMean,rStdErr,iMean,iStderr
             !end if

             ! Assume uncorrelated data (incorrect, but anyhow)
             !potential = cmplx(0,1)*TimeDerivativeWilsonLoop/WilsonLoop
             !if(r==1) then
             !   write(FileID_potential,'(1(SP,E16.9,1X))',advance='no') time
             !end if
             !if(r<rmax) then
             !   write(FileID_potential,'(2(SP,E16.9,1X))',advance='no') &
             !        real(potential),aimag(potential)
             !else
             !   write(FileID_potential,'(2(SP,E16.9,1X))',advance='yes') &
             !        real(potential),aimag(potential)
             !end if

          end do
       end do

       call CloseFile(FileID_WilsonLoops)
       !call CloseFile(FileID_TimeDerivativeWilsonLoops)
       !call CloseFile(FileID_Potential)
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

      ! start time
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') tstart
      ! Center of mass time-range smax
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') TimeRange

      TimeSteps=ceiling(TimeRange/LatticeSpacings(0))

      ! Seed for random number generator
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I4)') RandomNumberSeed
      if(thisproc()==0) print*,'seed'
      ! Initial gluon distribution (box): Saturation scale
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') Beta
      if(thisproc()==0) print*,'beta'
      
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I4)') nMeasurement
      if(thisproc()==0) print*,'nmeasurement'
      
      ! Output filenames
      arg_count = arg_count +1; call get_command_argument(arg_count,FileName_WilsonLoops);
     !arg_count = arg_count +1; call get_command_argument(arg_count,FileName_TimeDerivativeWilsonLoops);
      !arg_count = arg_count +1; call get_command_argument(arg_count,FileName_Potential)

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
  end subroutine MeasureWilsonLines_Equilibrium
  
  impure subroutine DeterminePotential_oneTimePoint
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
    use statistics, only: GetMean, GetStdError
    implicit none

    ! Simulation parameters
    integer(int64) :: LatticeExtensions(ndim)
    real(fp)       :: LatticeSpacings(0:ndim)
    integer(int64) :: RandomNumberSeed

    real(fp) :: GluonSaturationScale !qs
    real(fp) :: GluonOccupationAmplitude ! Amplitude of box in units of 1/g^2
    real(fp) :: GluonCoupling
    real(fp) :: tstart
    
    ! Physical fields
    type(GaugeConfiguration) :: GaugeConf, GaugeConf_initial
    
    ! Counting
    integer :: i

    real(fp) :: tEnd
    integer(int64) :: it

    ! Wilson loop parameters
    integer(int8), parameter :: messdir=nDim
    integer(int64) :: rmax, r

    ! Observables
    complex(fp), allocatable :: WilsonLoops(:,:), PotentialWilsonLoops(:,:)
    real(fp), allocatable :: rObservable(:), iObservable(:)
    real(fp) :: rMean, rStderr,iMean,iStderr
    real(fp) :: time
    
    ! Output
    integer(int8) :: FileID

    character(len=80) :: &
         FileName_WilsonLoops, FileName_PotentialWilsonLoops

    integer(intmpi) :: proc

    integer(int64) :: measurement
    integer(int64), parameter :: nMeasurement=100

    call InitSimulation

    rmax = LatticeExtensions(messdir)/2

    allocate(WilsonLoops(nMeasurement,rmax))
    WilsonLoops = 0
    allocate(PotentialWilsonLoops(nMeasurement,rmax))
    PotentialWilsonLoops = 0

    do measurement=1,nMeasurement
       if(ThisProc()==0) write(output_unit,*)&
            int(measurement,int16),'of',&
            int(nMeasurement,int16),'configurations';&
            call flush(output_unit)
       
       ! initialisation of config ....
       call GaugeConf_initial%TransversePolarisedOccupiedInit_Box(&
            GluonSaturationScale,GluonOccupationAmplitude,GluonCoupling)
       GaugeConf = GaugeConf_initial

       do it=1,abs(NINT(tstart/LatticeSpacings(0)))
          call GaugeConf%Update(sign(+1._real64,tstart))
       end do

       do r=1,rmax
          WilsonLoops(measurement,r) = &
               GetWilsonLoop(GaugeConf_initial,GaugeConf,r,messdir)
          
          PotentialWilsonLoops(measurement,r) = &
               GetPotentialWilsonLoop(GaugeConf_initial,GaugeConf,r,messdir)
       end do
    end do

    if(ThisProc()==0) then
       allocate(rObservable(rmax))
       allocate(iObservable(rmax))
       
       fileID = OpenFile(filename=FileName_WilsonLoops,&
            st='REPLACE',fm='FORMATTED',act='WRITE')
       do r=1,rmax
          rObservable = real(WilsonLoops(:,r))
          iObservable = aimag(WilsonLoops(:,r))
          
          rMean = GetMean(rObservable)
          iMean = GetMean(iObservable)
          rStdErr = GetStdError(rObservable)
          iStdErr = GetStdError(iObservable)
          
          write(FileID,'(1(I2,1X),4(SP,E16.9,1X))') r,rMean,rStdErr,iMean,iStderr
       end do
       call CloseFile(FileID)

       
       fileID = OpenFile(filename=FileName_PotentialWilsonLoops,&
            st='REPLACE',fm='FORMATTED',act='WRITE')
       do r=1,rmax
          rObservable = real(PotentialWilsonLoops(:,r))
          iObservable = aimag(PotentialWilsonLoops(:,r))
          
          rMean = GetMean(rObservable)
          iMean = GetMean(iObservable)
          rStdErr = GetStdError(rObservable)
          iStdErr = GetStdError(iObservable)
          
          write(FileID,'(1(I2,1X),4(SP,E16.9,1X))') r,rMean,rStdErr,iMean,iStderr
       end do
       call CloseFile(FileID)
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

      ! start time
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') tstart

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

      ! Output filenames
      arg_count = arg_count +1; call get_command_argument(arg_count,FileName_WilsonLoops);
      arg_count = arg_count +1; call get_command_argument(arg_count,FileName_PotentialWilsonLoops);

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
  end subroutine DeterminePotential_oneTimePoint
  
  impure subroutine DeterminePotential
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
    use statistics

    implicit none

    ! Simulation parameters
    integer(int64) :: LatticeExtensions(ndim)
    real(fp)       :: LatticeSpacings(0:ndim)
    integer(int64) :: TimeSteps
    integer(int64) :: RandomNumberSeed

    real(fp) :: GluonSaturationScale !qs
    real(fp) :: GluonOccupationAmplitude ! Amplitude of box in units of 1/g^2
    real(fp) :: GluonCoupling
    real(fp) :: tstart
    real(fp) :: TimeRange

    ! Physical fields
    type(GaugeConfiguration) :: GaugeConf, GaugeConf_initial

    ! Counting
    integer :: i

    integer(int64) :: it

    ! Wilson loop parameters
    integer(int8), parameter :: messdir=nDim
    integer(int64) :: rmax, r

    ! Observables
    complex(fp) :: WilsonLoop, TimeDerivativeWilsonLoop, Potential

    complex(fp), allocatable :: WilsonLoops(:,:,:), TimeDerivativeWilsonLoops(:,:,:)
    real(fp), allocatable :: rObservable(:), iObservable(:)
    real(fp) :: rMean, rStderr,iMean,iStderr
    real(fp) :: time
    integer(int64) :: TimePoints
    ! Output
    integer(int8) :: FileID_WilsonLoops, FileID_TimeDerivativeWilsonLoops, FileID_Potential

    character(len=80) :: &
         FileName_WilsonLoops, FileName_TimeDerivativeWilsonLoops, FileName_Potential

    integer(intmpi) :: proc

    integer(int64) :: measurement
    integer(int64), parameter :: nMeasurement=10

    call InitSimulation

    rmax = LatticeExtensions(messdir)/2

    TimePoints = nint(TimeRange/LatticeSpacings(0))
    
    allocate(WilsonLoops(nMeasurement,rmax,0:TimePoints))
    WilsonLoops = 0
    allocate(TimeDerivativeWilsonLoops(nMeasurement,rmax,0:TimePoints))
    TimeDerivativeWilsonLoops = 0

    do measurement=1,nMeasurement

       if(ThisProc()==0) write(output_unit,*)&
            int(measurement,int16),'of',&
            int(nMeasurement,int16),'configurations';&
            call flush(output_unit)

       ! initialisation of config ....
       call GaugeConf_initial%TransversePolarisedOccupiedInit_Box(&
            GluonSaturationScale,GluonOccupationAmplitude,GluonCoupling)
       GaugeConf = GaugeConf_initial

       do it=1,abs(NINT(tstart/LatticeSpacings(0)))
          call GaugeConf%Update(sign(+1._real64,tstart))
       end do

       do it=0,nint(TimeRange/LatticeSpacings(0)),+1

          do r=1,rmax
             WilsonLoop = GetWilsonLoop(GaugeConf_initial,GaugeConf,r,messdir)
             TimeDerivativeWilsonLoop = &
                  GetTimeDerivativeWilsonLoop(GaugeConf_initial,GaugeConf,r,messdir)
             !GetPotentialWilsonLoop(GaugeConf_initial,GaugeConf,r,messdir)

             WilsonLoops(measurement,r,it) = WilsonLoop
             TimeDerivativeWilsonLoops(measurement,r,it) = TimeDerivativeWilsonLoop
          end do

          call GaugeConf%Update
       end do
    end do
    if(ThisProc()==0) then
       fileID_WilsonLoops = OpenFile(filename=FileName_WilsonLoops,&
            st='REPLACE',fm='FORMATTED',act='WRITE')
       fileID_TimeDerivativeWilsonLoops = OpenFile(filename=FileName_TimeDerivativeWilsonLoops,&
            st='REPLACE',fm='FORMATTED',act='WRITE')
       fileID_Potential = OpenFile(filename=FileName_Potential,&
            st='REPLACE',fm='FORMATTED',act='WRITE')
       do it=0,nint(TimeRange/LatticeSpacings(0)),+1

          time = tstart + it*LatticeSpacings(0)
          if(ThisProc()==0) then
             write(output_unit,*) real(time,real32)
             call flush(output_unit)
          end if
          do r=1,rmax
             if(r==1) then
                write(FileID_WilsonLoops,'(1(SP,E16.9,1X))',advance='no') time
             end if
             rObservable = real(WilsonLoops(:,r,it))
             iObservable = aimag(WilsonLoops(:,r,it))

             rMean = GetMean(rObservable)
             iMean = GetMean(iObservable)
             rStdErr = GetStdError(rObservable)
             iStdErr = GetStdError(iObservable)

             WilsonLoop = cmplx(rMean,iMean)
             
             if(r<rmax) then
                write(FileID_WilsonLoops,'(4(SP,E16.9,1X))',advance='no') &
                     rMean,rStdErr,iMean,iStderr
             else
                write(FileID_WilsonLoops,'(4(SP,E16.9,1X))',advance='yes') &
                     rMean,rStdErr,iMean,iStderr
             end if

             
             if(r==1) then
                write(FileID_TimeDerivativeWilsonLoops,'(1(SP,E16.9,1X))',advance='no') time
             end if
             rObservable = real(TimeDerivativeWilsonLoops(:,r,it))
             iObservable = aimag(TimeDerivativeWilsonLoops(:,r,it))

             rMean = GetMean(rObservable)
             iMean = GetMean(iObservable)
             rStdErr = GetStdError(rObservable)
             iStdErr = GetStdError(iObservable)

             TimeDerivativeWilsonLoop = cmplx(rMean,iMean)
             
             if(r<rmax) then
                write(FileID_TimeDerivativeWilsonLoops,'(4(SP,E16.9,1X))',advance='no') &
                     rMean,rStdErr,iMean,iStderr
             else
                write(FileID_TimeDerivativeWilsonLoops,'(4(SP,E16.9,1X))',advance='yes') &
                     rMean,rStdErr,iMean,iStderr
             end if

             ! Assume uncorrelated data (incorrect, but anyhow)
             potential = cmplx(0,1)*TimeDerivativeWilsonLoop/WilsonLoop
             if(r==1) then
                write(FileID_potential,'(1(SP,E16.9,1X))',advance='no') time
             end if
             if(r<rmax) then
                write(FileID_potential,'(2(SP,E16.9,1X))',advance='no') &
                     real(potential),aimag(potential)
             else
                write(FileID_potential,'(2(SP,E16.9,1X))',advance='yes') &
                     real(potential),aimag(potential)
             end if

          end do
       end do

       call CloseFile(FileID_WilsonLoops)
       call CloseFile(FileID_TimeDerivativeWilsonLoops)
       call CloseFile(FileID_Potential)
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

      ! start time
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(F10.13)') tstart
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

      ! Output filenames
      arg_count = arg_count +1; call get_command_argument(arg_count,FileName_WilsonLoops);
      arg_count = arg_count +1; call get_command_argument(arg_count,FileName_TimeDerivativeWilsonLoops);
      arg_count = arg_count +1; call get_command_argument(arg_count,FileName_Potential)

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


    
  end subroutine DeterminePotential

  
  impure subroutine DetermineCustomPotentials
    use, intrinsic :: iso_fortran_env
    use precision
    use mpiinterface
    use lattice
    use gaugeconfiguration_su3
    use mpi
    use io
    use halocomm
    use nrqcd
    use wilsonline
    use random

    implicit none

    ! Simulation parameters
    integer(int64) :: LatticeExtensions(ndim)
    real(fp)       :: LatticeSpacings(0:ndim)
    integer(int64) :: TimeSteps
    integer(int64) :: RandomNumberSeed

    real(fp) :: GluonSaturationScale !qs
    real(fp) :: GluonOccupationAmplitude ! Amplitude of box in units of 1/g^2
    real(fp) :: GluonCoupling
    real(fp) :: tstart
    real(fp) :: TimeRange
    real(fp) :: HeavyQuarkmass
    complex(fp) :: WilsonCoefficients(nWilsonCoefficients)

    ! Physical fields
    type(GaugeConfiguration) :: GaugeConf, GaugeConf_atStart, ColdGaugeConf
    type(NRQCDField)         :: HeavyField,HeavyField_atStart,HeavyField_freeStep

    ! Counting
    integer :: i

    real(fp) :: tEnd
    integer(int64) :: it

    ! Wilson line parameters
    integer(int8), parameter :: messdir=nDim
    integer(int64), parameter :: x0=1
    integer(int64) :: rmax, r

    ! Observable-arrays
    real(fp), allocatable :: time(:)
    complex(fp),allocatable :: &
         GluonicWilsonLoops(:,:), PointSplitLoops(:,:), FreePointSplitLoops(:,:),&
         GluonicPotentials(:,:), PointSplitPotentials(:,:)
    
    ! Output
    integer(int8) :: &
         FileID_GluonicWilsonLoops,&
         FileID_PointSplitLoops,&
         FileID_FreePointSplitLoops,&
         FileID_GluonicPotentials,&
         FileID_PointSplitPotentials

    character(len=80) ::&
         FileGluonicWilsonLoops,&
         FilePointSplitLoops,&
         FileFreePointSplitLoops,&
         FileGluonicPotentials,&
         FilePointSplitPotentials
    
    integer(intmpi) :: proc

    call InitSimulation
    
    rmax = LatticeExtensions(messdir)/2/2
    
    ! Allocating observable arrays
    allocate(Time(0:nint(TimeRange/LatticeSpacings(0))))
    allocate(GluonicWilsonLoops(0:nint(TimeRange/LatticeSpacings(0)),0:rmax))
    allocate(PointSplitLoops(0:nint(TimeRange/LatticeSpacings(0)),0:rmax))
    allocate(FreePointSplitLoops(0:nint(TimeRange/LatticeSpacings(0)),0:rmax))
    allocate(GluonicPotentials(0:nint(TimeRange/LatticeSpacings(0)),0:rmax))
    allocate(PointSplitPotentials(0:nint(TimeRange/LatticeSpacings(0)),0:rmax))

    call GaugeConf_atStart%TransversePolarisedOccupiedInit_Box(&
         GluonSaturationScale,GluonOccupationAmplitude,GluonCoupling)
    call ColdGaugeConf%ColdInit
    
    ! Evolve gauge configuration to t0
    do it=1,abs(NINT(tstart/LatticeSpacings(0)))
       call GaugeConf_atStart%Update(sign(+1._real64,tstart))
    end do
    
    
    ! Initialising quark-antiquark-pair
    call HeavyField_atStart%InitSinglePoint(&
         latticeindex_quark=x0,&
         latticeindex_antiq=x0)

    ! Positive time evolution
    GaugeConf =GaugeConf_atStart
    HeavyField=HeavyField_atStart

    ! Opening files
    if(ThisProc()==0) then
       ! Opening files
       fileID_GluonicWilsonLoops  = OpenFile(filename=FileGluonicWilsonLoops,&
            st='REPLACE',fm='FORMATTED',act='WRITE')
       fileID_PointSplitLoops= OpenFile(filename=FilePointSplitLoops,&
            st='REPLACE',fm='FORMATTED',act='WRITE')
       fileID_FreePointSplitLoops= OpenFile(filename=FileFreePointSplitLoops,&
            st='REPLACE',fm='FORMATTED',act='WRITE')
       FileID_GluonicPotentials= OpenFile(filename=FileGluonicPotentials,&
            st='REPLACE',fm='FORMATTED',act='WRITE')
       FileID_PointSplitPotentials= OpenFile(filename=FilePointSplitPotentials,&
            st='REPLACE',fm='FORMATTED',act='WRITE')
    end if

    
    do it=0,nint(TimeRange/LatticeSpacings(0)),+1
       time(it) = tstart + it*LatticeSpacings(0)
       if(thisproc()==0) then
          write(output_unit,*) real(time(it),real32)
          call flush(6)
       end if

       HeavyField_freeStep = HeavyField
       
       do r=0,rmax
          GluonicWilsonLoops(it,r) &
               = GetGluonicWilsonLoop(GaugeConf_atStart,GaugeConf,x0, 2*r, messdir)
          PointSplitLoops(it,r) &
               = GetPointSplitCorrelator(GaugeConf,HeavyField, x0, 2*r, messdir)
          FreePointSplitLoops(it,r) &
               = GetPointSplitCorrelator(ColdGaugeConf,HeavyField_freeStep,x0,2*r,messdir)

          if(it>0) then
             ! Gluonic potential in analogy to Wilson's static potential definition
             ! combined with heuristic potential definition
             GluonicPotentials(it-1,r) &
                  =cmplx(0,1,fp)*(&
                  GluonicWilsonLoops(it,r) - GluonicWilsonLoops(it-1,r))&
                  /GluonicWilsonLoops(it-1,r)&
                  /LatticeSpacings(0)

             ! Fermionic wilsonloop from Rothkopf's definition
             ! used as input for heuristic definition of a potential
             PointSplitPotentials(it-1,r) &
                  =cmplx(0,1,fp)*(&
                  PointSplitLoops(it,r) - FreePointSplitLoops(it,r))&
                  /PointSplitLoops(it-1,r)&
                  /LatticeSpacings(0)
          end if
          
          if(ThisProc()==0.and.it>0) then
             if(r==0) then
                write(fileID_GluonicWilsonLoops,'(1(SP,E16.9,1X))',advance='no') time(it-1)
                write(fileID_PointSplitLoops,'(1(SP,E16.9,1X))',advance='no') time(it-1)
                write(fileID_FreePointSplitLoops,'(1(SP,E16.9,1X))',advance='no') time(it-1)
             end if
             if(r<rmax) then
                write(fileID_GluonicWilsonLoops,'(2(SP,E16.9,1X))',advance='no') &
                     real(GluonicWilsonLoops(it-1,r)), aimag(GluonicWilsonLoops(it-1,r))

                write(fileID_PointSplitLoops,'(2(SP,E16.9,1X))',advance='no') &
                     real(PointSplitLoops(it-1,r)), aimag(PointSplitLoops(it-1,r))

                write(fileID_FreePointSplitLoops,'(2(SP,E16.9,1X))',advance='no') &
                     real(FreePointSplitLoops(it-1,r)), aimag(FreePointSplitLoops(it-1,r))

                if(it<ubound(time,1)) then
                   if(r==0) then
                      write(fileID_GluonicPotentials,'(1(SP,E16.9,1X))',advance='no') time(it-1)
                      write(fileID_PointSplitPotentials,'(1(SP,E16.9,1X))',advance='no') time(it-1)
                   end if
                   write(fileID_GluonicPotentials,'(2(SP,E16.9,1X))',advance='no') &
                        real(GluonicPotentials(it-1,r)), aimag(GluonicPotentials(it-1,r))

                   write(fileID_PointSplitPotentials,'(2(SP,E16.9,1X))',advance='no') &
                        real(PointSplitPotentials(it-1,r)), aimag(PointSplitPotentials(it-1,r))
                end if
             else
                if(r==0) then
                   write(fileID_GluonicWilsonLoops,'(1(SP,E16.9,1X))',advance='no') time(it-1)
                   write(fileID_PointSplitLoops,'(1(SP,E16.9,1X))',advance='no') time(it-1)
                   write(fileID_FreePointSplitLoops,'(1(SP,E16.9,1X))',advance='no') time(it-1)
                end if
                write(fileID_GluonicWilsonLoops,'(2(SP,E16.9,1X))',advance='yes') &
                     real(GluonicWilsonLoops(it-1,r)), aimag(GluonicWilsonLoops(it-1,r))

                write(fileID_PointSplitLoops,'(2(SP,E16.9,1X))',advance='yes') &
                     real(PointSplitLoops(it-1,r)), aimag(PointSplitLoops(it-1,r))

                write(fileID_FreePointSplitLoops,'(2(SP,E16.9,1X))',advance='yes') &
                     real(FreePointSplitLoops(it-1,r)), aimag(FreePointSplitLoops(it-1,r))

                if(it<ubound(time,1)) then
                   
                   if(r==0) then
                      write(fileID_GluonicPotentials,'(1(SP,E16.9,1X))',advance='no') time(it-1)
                      write(fileID_PointSplitPotentials,'(1(SP,E16.9,1X))',advance='no') time(it-1)
                   end if
                   write(fileID_GluonicPotentials,'(2(SP,E16.9,1X))',advance='yes') &
                        real(GluonicPotentials(it-1,r)), aimag(GluonicPotentials(it-1,r))

                   write(fileID_PointSplitPotentials,'(2(SP,E16.9,1X))',advance='yes') &
                        real(PointSplitPotentials(it-1,r)), aimag(PointSplitPotentials(it-1,r))
                end if
             end if
          end if
       end do
       
       call HeavyField_freeStep%Update(ColdGaugeConf,HeavyQuarkMass,WilsonCoefficients,+1._fp)
       call HeavyField%Update(GaugeConf,HeavyQuarkMass,WilsonCoefficients,+1._fp)
       call GaugeConf%Update(+1._fp)
    end do

    if(ThisProc()==0) then

       ! Closing files
       call CloseFile(FileID_GluonicWilsonloops)
       call CloseFile(FileID_PointSplitLoops)
       call CloseFile(FileID_FreePointSplitLoops)
       call CloseFile(FileID_GluonicPotentials)
       call CloseFile(FileID_PointSplitPotentials)
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

      ! Output filenames
      arg_count = arg_count +1; call get_command_argument(arg_count,FileGluonicWilsonLoops);
      arg_count = arg_count +1; call get_command_argument(arg_count,FilePointSplitLoops);
      arg_count = arg_count +1; call get_command_argument(arg_count,FileFreePointSplitLoops);
      arg_count = arg_count +1; call get_command_argument(arg_count,FileGluonicPotentials);
      arg_count = arg_count +1; call get_command_argument(arg_count,FilePointSplitPotentials);

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
  end subroutine DetermineCustomPotentials
    
  
  impure subroutine DeterminePotentials
    use, intrinsic :: iso_fortran_env
    use precision
    use mpiinterface
    use lattice
    use gaugeconfiguration_su3
    use mpi
    use io
    use halocomm
    use nrqcd
    use wilsonline
    use random

    implicit none

    ! Simulation parameters
    integer(int64) :: LatticeExtensions(ndim)
    real(fp)       :: LatticeSpacings(0:ndim)
    integer(int64) :: TimeSteps
    integer(int64) :: RandomNumberSeed

    real(fp) :: GluonSaturationScale !qs
    real(fp) :: GluonOccupationAmplitude ! Amplitude of box in units of 1/g^2
    real(fp) :: GluonCoupling
    real(fp) :: tstart
    real(fp) :: TimeRange
    real(fp) :: HeavyQuarkmass
    complex(fp) :: WilsonCoefficients(nWilsonCoefficients)

    ! Physical fields
    type(GaugeConfiguration) :: GaugeConf, GaugeConf_atStart, ColdGaugeConf
    type(NRQCDField)         :: HeavyField,HeavyField_atStart,HeavyField_freeStep

    ! Counting
    integer :: i

    real(fp) :: tEnd
    integer(int64) :: it

    ! Wilson line parameters
    integer(int8), parameter :: messdir=nDim
    integer(int64), parameter :: x0=1
    integer(int64) :: rmax, r, xr

    ! Observable-arrays
    real(fp), allocatable :: time(:)
    complex(fp),allocatable :: &
         GluonicWilsonLoops(:,:), FermionicWilsonLoops(:,:), FreeFermionicWilsonLoops(:,:),&
         GluonicPotentials(:,:), FermionicPotentials(:,:)
    
    ! Output
    integer(int8) :: &
         FileID_GluonicWilsonLoops,&
         FileID_FermionicWilsonLoops,&
         FileID_FreeFermionicWilsonLoops,&
         FileID_GluonicPotentials,&
         FileID_FermionicPotentials

    character(len=80) ::&
         FileGluonicWilsonLoops,&
         FileFermionicWilsonLoops,&
         FileFreeFermionicWilsonLoops,&
         FileGluonicPotentials,&
         FileFermionicPotentials
    
    integer(intmpi) :: proc

    call InitSimulation
    
    rmax = LatticeExtensions(messdir)/2
    
    ! Allocating observable arrays
    allocate(Time(0:nint(TimeRange/LatticeSpacings(0))))
    allocate(GluonicWilsonLoops(0:nint(TimeRange/LatticeSpacings(0)),0:rmax))
    allocate(FermionicWilsonLoops(0:nint(TimeRange/LatticeSpacings(0)),0:rmax))
    allocate(FreeFermionicWilsonLoops(0:nint(TimeRange/LatticeSpacings(0)),0:rmax))
    allocate(GluonicPotentials(0:nint(TimeRange/LatticeSpacings(0)),0:rmax))
    allocate(FermionicPotentials(0:nint(TimeRange/LatticeSpacings(0)),0:rmax))

    call GaugeConf_atStart%TransversePolarisedOccupiedInit_Box(&
         GluonSaturationScale,GluonOccupationAmplitude,GluonCoupling)
    call ColdGaugeConf%ColdInit
    
    ! Evolve gauge configuration to t0
    do it=1,abs(NINT(tstart/LatticeSpacings(0)))
       call GaugeConf_atStart%Update(sign(+1._real64,tstart))
    end do
    
    xr = x0
    do r=0,rmax
       if(ThisProc()==0) then
          write(output_unit,*) 'r=',r,'of',rmax
       end if
       
       ! Initialising quark-antiquark-pair
       ! with quark at variable remote point xr
       ! and antiquark at fixed (origin) point x0
       call HeavyField_atStart%InitSinglePoint(&
            latticeindex_quark=xr,&
            latticeindex_antiq=x0)
       
       ! Positive time evolution
       GaugeConf =GaugeConf_atStart
       HeavyField=HeavyField_atStart

       do it=0,nint(TimeRange/LatticeSpacings(0)),+1
          if(r==0) then
             time(it) = tstart + it*LatticeSpacings(0)
          end if

          
          HeavyField_freeStep = HeavyField

          
          GluonicWilsonLoops(it,r) &
               = GetGluonicWilsonLoop(GaugeConf_atStart,GaugeConf,x0, r, messdir)
          FermionicWilsonLoops(it,r) &
               = GetFermionicWilsonLoop(GaugeConf_atStart,GaugeConf,HeavyField, x0, r, messdir)
          FreeFermionicWilsonLoops(it,r) &
               = GetFermionicWilsonLoop(ColdGaugeConf,ColdGaugeConf,HeavyField_freeStep,x0,r,messdir)
          
          call HeavyField_freeStep%Update(ColdGaugeConf,HeavyQuarkMass,WilsonCoefficients,+1._fp)
          call HeavyField%Update(GaugeConf,HeavyQuarkMass,WilsonCoefficients,+1._fp)
          call GaugeConf%Update(+1._fp)

          if(it>0) then
             ! Gluonic potential in analogy to Wilson's static potential definition
             ! combined with heuristic potential definition
             GluonicPotentials(it-1,r) &
                  =cmplx(0,1,fp)*(&
                  GluonicWilsonLoops(it,r) - GluonicWilsonLoops(it-1,r))&
                  /GluonicWilsonLoops(it-1,r)&
                  /LatticeSpacings(0)

             ! Fermionic wilsonloop from Rothkopf's definition
             ! used as input for heuristic definition of a potential
             FermionicPotentials(it-1,r) &
                  =cmplx(0,1,fp)*(&
                  FermionicWilsonLoops(it,r) - FreeFermionicWilsonLoops(it,r))&
                  /FermionicWilsonLoops(it-1,r)&
                  /LatticeSpacings(0)
          end if
       end do

       ! Preparing next step
       xr = GetNeib_G(messdir,xr)
    end do

    if(ThisProc()==0) then
       ! Opening files
       fileID_GluonicWilsonLoops  = OpenFile(filename=FileGluonicWilsonLoops,&
            st='REPLACE',fm='FORMATTED',act='WRITE')
       fileID_FermionicWilsonLoops= OpenFile(filename=FileFermionicWilsonLoops,&
            st='REPLACE',fm='FORMATTED',act='WRITE')
       fileID_FreeFermionicWilsonLoops= OpenFile(filename=FileFreeFermionicWilsonLoops,&
            st='REPLACE',fm='FORMATTED',act='WRITE')
       FileID_GluonicPotentials= OpenFile(filename=FileGluonicPotentials,&
            st='REPLACE',fm='FORMATTED',act='WRITE')
       FileID_FermionicPotentials= OpenFile(filename=FileFermionicPotentials,&
            st='REPLACE',fm='FORMATTED',act='WRITE')
       
       do it=lbound(time,1),ubound(time,1)
          write(fileID_GluonicWilsonLoops,'(1(SP,E16.9,1X))',advance='no') time(it)
          write(fileID_FermionicWilsonLoops,'(1(SP,E16.9,1X))',advance='no') time(it)
          write(fileID_FreeFermionicWilsonLoops,'(1(SP,E16.9,1X))',advance='no') time(it)
          write(fileID_GluonicPotentials,'(1(SP,E16.9,1X))',advance='no') time(it)
          write(fileID_FermionicPotentials,'(1(SP,E16.9,1X))',advance='no') time(it)
          
          do r=0,rmax
             if(r<rmax) then
                write(fileID_GluonicWilsonLoops,'(2(SP,E16.9,1X))',advance='no') &
                     real(GluonicWilsonLoops(it,r)), aimag(GluonicWilsonLoops(it,r))
                
                write(fileID_FermionicWilsonLoops,'(2(SP,E16.9,1X))',advance='no') &
                     real(FermionicWilsonLoops(it,r)), aimag(FermionicWilsonLoops(it,r))
                
                write(fileID_FreeFermionicWilsonLoops,'(2(SP,E16.9,1X))',advance='no') &
                     real(FreeFermionicWilsonLoops(it,r)), aimag(FreeFermionicWilsonLoops(it,r))

                if(it<ubound(time,1)) then
                   write(fileID_GluonicPotentials,'(2(SP,E16.9,1X))',advance='no') &
                        real(GluonicPotentials(it,r)), aimag(GluonicPotentials(it,r))

                   write(fileID_FermionicPotentials,'(2(SP,E16.9,1X))',advance='no') &
                        real(FermionicPotentials(it,r)), aimag(FermionicPotentials(it,r))
                end if
             else
                write(fileID_GluonicWilsonLoops,'(2(SP,E16.9,1X))',advance='yes') &
                     real(GluonicWilsonLoops(it,r)), aimag(GluonicWilsonLoops(it,r))
                
                write(fileID_FermionicWilsonLoops,'(2(SP,E16.9,1X))',advance='yes') &
                     real(FermionicWilsonLoops(it,r)), aimag(FermionicWilsonLoops(it,r))
                
                write(fileID_FreeFermionicWilsonLoops,'(2(SP,E16.9,1X))',advance='yes') &
                     real(FreeFermionicWilsonLoops(it,r)), aimag(FreeFermionicWilsonLoops(it,r))

                if(it<ubound(time,1)) then
                   write(fileID_GluonicPotentials,'(2(SP,E16.9,1X))',advance='yes') &
                        real(GluonicPotentials(it,r)), aimag(GluonicPotentials(it,r))

                   write(fileID_FermionicPotentials,'(2(SP,E16.9,1X))',advance='yes') &
                        real(FermionicPotentials(it,r)), aimag(FermionicPotentials(it,r))
                end if
             end if
          end do
       end do

       ! Closing files
       call CloseFile(FileID_GluonicWilsonloops)
       call CloseFile(FileID_FermionicWilsonLoops)
       call CloseFile(FileID_FreeFermionicWilsonLoops)
       call CloseFile(FileID_GluonicPotentials)
       call CloseFile(FileID_FermionicPotentials)
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

      ! Output filenames
      arg_count = arg_count +1; call get_command_argument(arg_count,FileGluonicWilsonLoops);
      arg_count = arg_count +1; call get_command_argument(arg_count,FileFermionicWilsonLoops);
      arg_count = arg_count +1; call get_command_argument(arg_count,FileFreeFermionicWilsonLoops);
      arg_count = arg_count +1; call get_command_argument(arg_count,FileGluonicPotentials);
      arg_count = arg_count +1; call get_command_argument(arg_count,FileFermionicPotentials);

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
    
  end subroutine DeterminePotentials

  
  !>@brief Program for measuring the various wilson lines
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 02.04.2019
  !!@version 1.0
  impure subroutine MeasureWilsonLines
    use, intrinsic :: iso_fortran_env
    use precision
    use mpiinterface
    use lattice
    use gaugeconfiguration_su3
    use mpi
    use io
    use halocomm
    use nrqcd
    use wilsonline
    use random
    
    implicit none

    ! Simulation parameters
    integer(int64) :: LatticeExtensions(ndim)
    real(fp)       :: LatticeSpacings(0:ndim)
    integer(int64) :: TimeSteps
    integer(int64) :: RandomNumberSeed

    real(fp) :: GluonSaturationScale !qs
    real(fp) :: GluonOccupationAmplitude ! Amplitude of box in units of 1/g^2
    real(fp) :: GluonCoupling
    real(fp) :: tstart
    real(fp) :: TimeRange
    real(fp) :: HeavyQuarkmass
    complex(fp) :: WilsonCoefficients(nWilsonCoefficients)

    ! Physical fields
    type(GaugeConfiguration) :: GaugeConf, GaugeConf_at0
    type(NRQCDField)         :: HeavyField,HeavyField_at0

    ! Counting
    integer :: i

    real(fp) :: tEnd
    integer(int64) :: it

    ! Wilson line parameters
    integer(int8), parameter :: messdir=nDim
    integer(int64), parameter :: x0=1
    integer(int64) :: rmax, r, xr

    ! Observable
    complex(fp) :: GluonicWilsonLoop, FermionicWilsonLoop, MesonCorrelator
    real(fp) :: quarknorm, antiqnorm

    ! Observable-arrays
    complex(fp),allocatable :: mesoncorrelators(:), GluonicWilsonLoops(:,:), FermionicWilsonLoops(:,:)
    real(fp),   allocatable :: time(:), quarknorms(:,:), antiqnorms(:,:)

    ! Output
    integer(int8) :: FileID_GluonicWilsonLoops,FileID_FermionicWilsonLoops,FileID_MesonCorrelator,&
         FileID_Norm

    character(len=80) :: FileGluonicWilsonLoops, FileFermionicWilsonLoops, FileMesonCorrelator,&
         FileNorm
    
    integer(intmpi) :: proc

    call InitSimulation
    
    rmax = LatticeExtensions(messdir)/2

    ! Allocating observable arrays
    allocate(Time(0:nint(TimeRange/LatticeSpacings(0))))
    allocate(MesonCorrelators(0:nint(TimeRange/LatticeSpacings(0))))
    allocate(AntiQNorms(0:nint(TimeRange/LatticeSpacings(0)),0:rmax))
    allocate(QuarkNorms(0:nint(TimeRange/LatticeSpacings(0)),0:rmax))
    allocate(GluonicWilsonLoops(0:nint(TimeRange/LatticeSpacings(0)),0:rmax))
    allocate(FermionicWilsonLoops(0:nint(TimeRange/LatticeSpacings(0)),0:rmax))
    
    call GaugeConf_at0%TransversePolarisedOccupiedInit_Box(&
         GluonSaturationScale,GluonOccupationAmplitude,GluonCoupling)

    ! Evolve gauge configuration to t0
    do it=1,abs(NINT(tstart/LatticeSpacings(0)))
       call GaugeConf_at0%Update(sign(+1._real64,tstart))
    end do
    
    xr = x0
    do r=0,rmax
       if(ThisProc()==0) then
          write(output_unit,*) 'r=',r,'of',rmax
       end if
       
       ! Initialising quark-antiquark-pair
       ! with quark at variable remote point xr
       ! and antiquark at fixed (origin) point x0
       call HeavyField_at0%InitSinglePoint(&
            latticeindex_quark=xr,&
            latticeindex_antiq=x0)
       
       ! Positive time evolution
       GaugeConf =GaugeConf_at0
       HeavyField=HeavyField_at0
       do it=0,nint(TimeRange/LatticeSpacings(0)),+1

          if(r==0) then
             time(it) = it*LatticeSpacings(0)
             
             mesoncorrelators(it) &
                  = HeavyField%GetMesonCorrelator_O1_3s1_ZeroMomentum()
          end if

          quarknorms(it,r) = HeavyField%GetNorm_Quark()
          antiqnorms(it,r) = HeavyField%GetNorm_AntiQ()
          
          GluonicWilsonLoops(it,r) &
               = GetGluonicWilsonLoop(GaugeConf_at0, GaugeConf, x0, r, messdir)
          FermionicWilsonLoops(it,r) &
               = GetFermionicWilsonLoop(GaugeConf_at0, GaugeConf, HeavyField, x0, r, messdir)
          
          call HeavyField%Update(GaugeConf,HeavyQuarkMass,WilsonCoefficients,+1._fp)
          call GaugeConf%Update(+1._fp)
       end do

       ! Preparing next step
       xr = GetNeib_G(messdir,xr)
    end do

    if(ThisProc()==0) then
       ! Opening files
       fileID_MesonCorrelator = OpenFile(filename=FileMesonCorrelator,&
            st='REPLACE',fm='FORMATTED',act='WRITE')
       fileID_Norm = OpenFile(filename=FileNorm,&
            st='REPLACE',fm='FORMATTED',act='WRITE')
       fileID_GluonicWilsonLoops  = OpenFile(filename=FileGluonicWilsonLoops,&
            st='REPLACE',fm='FORMATTED',act='WRITE')
       fileID_FermionicWilsonLoops= OpenFile(filename=FileFermionicWilsonLoops,&
            st='REPLACE',fm='FORMATTED',act='WRITE')
       
       do it=lbound(time,1),ubound(time,1)
          write(FileID_MesonCorrelator,'(3(SP,E16.9,1X))') &
               time(it),real(mesoncorrelators(it),real64),aimag(mesoncorrelators(it))

          write(fileID_Norm,'(1(SP,E16.9,1X))',advance='no') time(it)
          write(fileID_GluonicWilsonLoops,'(1(SP,E16.9,1X))',advance='no') time(it)
          write(fileID_FermionicWilsonLoops,'(1(SP,E16.9,1X))',advance='no') time(it)
          
          do r=0,rmax
             if(r<rmax) then
                write(fileID_Norm,'(2(SP,E16.9,1X))',advance='no') &
                     quarknorms(it,r), antiqnorms(it,r)
                
                write(fileID_GluonicWilsonLoops,'(2(SP,E16.9,1X))',advance='no') &
                     real(GluonicWilsonLoops(it,r)), aimag(GluonicWilsonLoops(it,r))
                
                write(fileID_FermionicWilsonLoops,'(2(SP,E16.9,1X))',advance='no') &
                     real(FermionicWilsonLoops(it,r)), aimag(FermionicWilsonLoops(it,r))
             else
                write(fileID_Norm,'(2(SP,E16.9,1X))',advance='yes') &
                     quarknorms(it,r), antiqnorms(it,r)
                
                write(fileID_GluonicWilsonLoops,'(2(SP,E16.9,1X))',advance='yes') &
                     real(GluonicWilsonLoops(it,r)), aimag(GluonicWilsonLoops(it,r))
                
                write(fileID_FermionicWilsonLoops,'(2(SP,E16.9,1X))',advance='yes') &
                     real(FermionicWilsonLoops(it,r)), aimag(FermionicWilsonLoops(it,r))
             end if
          end do
       end do

       ! Closing files
       call CloseFile(fileID_MesonCorrelator)
       call CloseFile(FileID_Norm)
       call CloseFile(FileID_GluonicWilsonloops)
       call CloseFile(FileID_FermionicWilsonLoops)
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

      ! Output filenames
      arg_count = arg_count +1; call get_command_argument(arg_count,FileGluonicWilsonLoops);
      arg_count = arg_count +1; call get_command_argument(arg_count,FileFermionicWilsonLoops);
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
  end subroutine MeasureWilsonLines
  
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
    integer(int64) :: TimeSteps
    integer(int64) :: RandomNumberSeed

    real(fp)   :: GluonSaturationScale !qs
    real(fp)   :: GluonOccupationAmplitude ! Amplitude of box in units of 1/g^2
    real(fp)   :: GluonCoupling
    integer(int64) :: EnsembleSize
    real(fp)   :: CoMTime
    real(fp)   :: TimeRange

    ! Physical fields
    type(GaugeConfiguration) :: GaugeConf

    ! Monitoring variables
    integer(int64) :: it, idistmeas
    real(fp) :: time

    integer(int64) :: MemoryIndex, LatticeIndex, i, is

    ! MPI
    integer(intmpi) :: mpierr, mpistatus(mpi_status_size), tag, src, dest, buffersize

    real(fp), allocatable :: aa_correlator_opt(:),ee_correlator_opt(:)
    ! Ensemble variables
    real(fp), allocatable :: AA_correlator_ensemble(:,:,:), EE_correlator_ensemble(:,:,:)
    real(fp), allocatable :: AA_correlator(:,:), AA_correlator_stderror(:,:)
    real(fp), allocatable :: EE_correlator(:,:), EE_correlator_stderror(:,:)
    real(fp), allocatable :: GluonDistribution(:,:), GluonDistribution_stderror(:,:)
    real(fp), allocatable :: &
         GluonDistribution_AtTimePoint(:), GluonDistribution_AtTimePoint_stderror(:)
    integer :: iensemble
    real(fp), allocatable :: momenta(:)

    integer(intmpi) :: proc

    ! Output
    integer(int8) :: FileID
    character(len=128) :: filename,time_tag
    integer(intmpi), allocatable :: mpisendrequest(:), mpisendstatus(:,:)
    
    integer(int64) :: GaugeFixingMaxIterations
    real(fp) :: GaugefixingCoefficient, GaugeFixingTolerance
    integer(int64) :: Number_of_Measurements_of_Gluondistribution
    integer(int64) :: TimePoints_between_Measurement_of_Gluondistribution
    real(fp) :: TimeBetweenGluonMeasurements

    type(GaugeConfiguration) :: gaugeconf_gaugefixed

    call InitSimulation

    if(ThisProc()==0) print*,GaugefixingCoefficient
    
    allocate(AA_correlator_ensemble(GetLocalLatticeSize(),&
         0:Number_of_Measurements_of_Gluondistribution,EnsembleSize))
    allocate(EE_correlator_ensemble(GetLocalLatticeSize(),&
         0:Number_of_Measurements_of_Gluondistribution,EnsembleSize))

    ensemble: do iensemble=1,EnsembleSize
       if(ThisProc()==0) write(output_unit,*)&
            int(iensemble,int16),'of',&
            int(EnsembleSize,int16),'configurations';&
            call flush(output_unit)

       !call GaugeConf%SemicoldInit
       call GaugeConf%TransversePolarisedOccupiedInit_Box(&
            GluonSaturationScale,GluonOccupationAmplitude,GluonCoupling&
       !,aa_correlator_opt,ee_correlator_opt&
            )
       TimeEvolution: do it=0,TimeSteps

          if(modulo(it,TimePoints_between_Measurement_of_Gluondistribution)==0) then
             if(ThisProc()==0) write(output_unit,*)&
                  it,'of',&
                  TimeSteps,'time steps';&
                  call flush(output_unit)

             idistmeas = it/TimePoints_between_Measurement_of_Gluondistribution

             gaugeconf_gaugefixed = GaugeConf
             call gaugeconf_gaugefixed%CoulombGaugefixing(&
                  Tolerance_    =GaugeFixingTolerance,&
                  alpha_        =GaugefixingCoefficient,&
                  MaxIterations_=GaugefixingMaxIterations)
             call gaugeconf_gaugefixed%GetTransverseAACorrelator(aa_correlator_opt)
             call gaugeconf_gaugefixed%GetTransverseEECorrelator(ee_correlator_opt)

             AA_correlator_ensemble(:,idistmeas,iensemble) = aa_correlator_opt
             EE_correlator_ensemble(:,idistmeas,iensemble) = ee_correlator_opt

          end if

          call GaugeConf%Update
       end do TimeEvolution
    end do ensemble

    ! ..--** Start: Statistics and output **--..
    it = 0
    idistmeas = 0
    call StatisticsCorrelator(AA_correlator_ensemble,AA_correlator,AA_correlator_stderror,0)
    call StatisticsCorrelator(EE_correlator_ensemble,EE_correlator,EE_correlator_stderror,0)

    allocate(GluonDistribution(&
         lbound(AA_correlator,1):ubound(AA_correlator,1),&
         lbound(AA_correlator,2):ubound(AA_correlator,2)))
    allocate(GluonDistribution_stderror(&
         lbound(AA_correlator,1):ubound(AA_correlator,1),&
         lbound(AA_correlator,2):ubound(AA_correlator,2)))

    GluonDistribution = &
                                !AA_correlator
         sqrt(AA_correlator*EE_correlator)
    GluonDistribution_stderror = &
                                !AA_correlator_stderror
         sqrt(&
         + (AA_correlator_stderror*EE_correlator)**2 &
         + (EE_correlator_stderror*AA_correlator)**2&
         )

    ! Printing to file
    dest=0

    allocate(GluonDistribution_AtTimePoint(GetLocalLatticeSize()))
    allocate(GluonDistribution_AtTimePoint_stderror(GetLocalLatticeSize()))
    allocate(Momenta(GetLocalLatticeSize()))
    is=0
    allocate(mpisendrequest(3))
    do idistmeas=lbound(GluonDistribution,2),ubound(GluonDistribution,2)
       if(ThisProc()==dest) then
          write(filename,"(A18,I0.3,A1,I0.3,A1,I0.3)") 'gluondistribution_', &
               LatticeExtensions(1),'x',LatticeExtensions(2),'x',&
               LatticeExtensions(3)

          time = (idistmeas - lbound(GluonDistribution,2)) &
               * TimeBetweenGluonMeasurements
          write(time_tag,"(F12.3)") time
          time_tag = '_t' // trim(ADJUSTL(time_tag))
          filename = trim(filename) // trim(time_tag) // '.txt'

          fileID = OpenFile(filename=filename,st='REPLACE',fm='FORMATTED',act='WRITE')
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
             buffersize = size(GluonDistribution_atTimePoint)
             tag = (1+NumProcs())*src
             if(ThisProc()==src) then
                GluonDistribution_atTimePoint = GluonDistribution(:,idistmeas)
                call mpi_isend(&
                     GluonDistribution_atTimePoint,& ! What to send
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
                     GluonDistribution_atTimePoint,& ! What to recieve
                     buffersize,                   & ! How many points
                     MPI_DOUBLE,                   & ! What type
                     src,                          & ! Sending process
                     tag,                          & ! Tag
                     MPI_COMM_WORLD,               & ! Communicator
                     mpistatus,                    & ! Status
                     mpierr)                         ! Error code
             end if

             ! Standard error of gluon distribution
             buffersize = size(GluonDistribution_atTimePoint)
             tag = (2+NumProcs())*src
             if(ThisProc()==src) then
                GluonDistribution_atTimePoint_stderror = GluonDistribution_stderror(:,idistmeas)
                call mpi_isend(&
                     GluonDistribution_atTimePoint_stderror,& ! What to send
                     buffersize,                   & ! How many points
                     MPI_DOUBLE,                   & ! What type
                     dest,                         & ! Recieving process
                     tag,                          & ! Tag
                     MPI_COMM_WORLD,               & ! Communicator
                     mpisendrequest(3),            & ! Request handle
                     mpierr)                         ! Error code
             end if
             if(ThisProc()==dest) then
                call mpi_recv(&
                     GluonDistribution_atTimePoint_stderror,& ! What to recieve
                     buffersize,                   & ! How many points
                     MPI_DOUBLE,                   & ! What type
                     src,                          & ! Sending process
                     tag,                          & ! Tag
                     MPI_COMM_WORLD,               & ! Communicator
                     mpistatus,                    & ! Status
                     mpierr)                         ! Error code
             end if


             if(ThisProc()==dest) then
                do is=1,size(Momenta)
                   write(fileID,'(3(SP,E13.6,1X))') &
                        Momenta(is),&
                        GluonDistribution_AtTimePoint(is),&
                        GluonDistribution_AtTimePoint_stderror(is)
                end do
             end if
          end if

          call SyncAll
       end do
       if(ThisProc()==0) call CloseFile(FileID)
    end do
    call endsimulation
  contains
    pure subroutine StatisticsCorrelator(&
         correlator_ensemble,correlator,correlator_stderror,lbound_time)
      use precision, only: fp
      use statistics, only: GetMean, GetStdError

      implicit none
      real(fp),intent(in) :: correlator_ensemble(:,lbound_time:,:)
      real(fp),allocatable,intent(out):: correlator(:,:)
      real(fp),allocatable,intent(out):: correlator_stderror(:,:)
      integer, intent(in) :: lbound_time

      integer :: latticeindex, it

      allocate(correlator(&
           lbound(correlator_ensemble,1):ubound(correlator_ensemble,1),&
           lbound(correlator_ensemble,2):ubound(correlator_ensemble,2)))

      allocate(correlator_stderror(&
           lbound(correlator_ensemble,1):ubound(correlator_ensemble,1),&
           lbound(correlator_ensemble,2):ubound(correlator_ensemble,2)))

      forall(&
           it          =lbound(correlator_ensemble,2):ubound(correlator_ensemble,2),&
           latticeindex=lbound(correlator_ensemble,1):ubound(correlator_ensemble,1))&
           correlator(latticeindex,it) = GetMean(correlator_ensemble(latticeindex,it,:))

      if(size(correlator_ensemble,3)>1) then
         forall(&
              it          =lbound(correlator_ensemble,2):ubound(correlator_ensemble,2),&
              latticeindex=lbound(correlator_ensemble,1):ubound(correlator_ensemble,1))&
              correlator_stderror(latticeindex,it) = GetStdError(correlator_ensemble(latticeindex,it,:))
      else
         correlator_stderror = 0
      end if
    end subroutine StatisticsCorrelator

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

      ! Ensemble size (statistical average of classical statistical simulation)
      arg_count = arg_count +1; call get_command_argument(arg_count,arg);
      read(arg,'(I6)') EnsembleSize

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
  case(2)
     call MeasureEnergyAndGaussLawDeviation
  case(3)
     call MeasureHeavyQuarkoniumCorrelators
  case (31)
     call MeasureHeavyQuarkoniumCorrelators_oneT
  case (32)
     call MeasureHeavyQuarkoniumCorrelators_free
  case(4)
     call ComputeSpectrum
  case (5)
     call MeasureWilsonLines
  case (51)
     call MeasureWilsonLines_Equilibrium
  case (52)
     call MeasureWilsonAndHybridLines_Equilibrium
  case (6)
     call DeterminePotentials
  case (7)
     call DetermineCustomPotentials
  case (8)
     call DeterminePotential
  case (9)
     call DeterminePotential_oneTimePoint
  case default
     call MPIStop('Invalid simulation mode selected')
  end select
end program simulation
