!------------------------------------------------------------------------------
! Fast Fourier Transforms between time and frequency domain
!------------------------------------------------------------------------------
!
! MODULE: dft
!>@brief Providing interfaces for fast fourier transforms between real and momentum space
!!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!!@date 22.03.2019
!!@version 1.0
! REVISION HISTORY:
! 22 03 2019 - Initial Version
!------------------------------------------------------------------------------
module twfft
  use, intrinsic :: iso_fortran_env ! defines kinds
  use, intrinsic :: iso_c_binding
  use mkl_dfti
  use mathconstants, only: pi
  use precision, only: fp
  implicit none

  PRIVATE

  public &
       t2w, w2t

contains
  !>@brief DFT for 1-d signal from frequency to time domain using MKL-CDFT
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.03.2019
  !!@version 1.0
  subroutine w2t(data_,dt)
    implicit none
    !> Fourier-transform input and output
    complex(fp), intent(inout) :: data_(:)
    !> Temporal spacing
    real(fp), intent(in) :: dt
    

    complex(real64), allocatable :: data(:)

    integer(int64) status

    integer(int64), parameter :: tw_nrank=1
    integer(int64) :: tw_length, tw_size
    TYPE(DFTI_DESCRIPTOR), POINTER :: tw_desc

    tw_length = size(data_)

    status = DftiCreateDescriptor(tw_desc,DFTI_DOUBLE, &
         DFTI_COMPLEX,tw_nrank,tw_length)
    status = DftiSetValue(tw_desc,&
         DFTI_FORWARD_SCALE,&
         real(1/(dt*tw_length),real64))
    status = DftiCommitDescriptor(tw_desc)

    allocate(data(size(data_)))
    data = data_
    status = DftiComputeForward(tw_desc,data)
    data_ = data
    deallocate(data)

    status = DftiFreeDescriptor(tw_desc)

  end subroutine w2t

  !>@brief DFT for 1-d signal from time to frequency domain using MKL-CDFT
  !!@author Alexander Lehmann, UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !!@date 22.03.2019
  !!@version 1.0
  subroutine t2w(data_,dt)
    implicit none
    !> Fourier-transform input and output
    complex(fp), intent(inout) :: data_(:)
    !> Temporal spacing
    real(fp), intent(in) :: dt

    complex(real64), allocatable :: data(:)
    
    integer(int64) status

    integer(int64), parameter :: tw_nrank=1
    integer(int64) :: tw_length, tw_size
    TYPE(DFTI_DESCRIPTOR), POINTER :: tw_desc

    tw_length = size(data_)

    status = DftiCreateDescriptor(tw_desc,DFTI_DOUBLE, &
         DFTI_COMPLEX,tw_nrank,tw_length)
    status = DftiSetValue(tw_desc,&
         DFTI_BACKWARD_SCALE,&
         real(dt,real64))
    status = DftiCommitDescriptor(tw_desc)

    allocate(data(size(data_)))
    data = data_
    status = DftiComputeBackward(tw_desc,data)
    data_ = data
    deallocate(data)
    
    status = DftiFreeDescriptor(tw_desc)

  end subroutine t2w
end module twfft
