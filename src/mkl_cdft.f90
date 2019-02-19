!===============================================================================
! Copyright 2002-2018 Intel Corporation.
!
! This software and the related documents are Intel copyrighted  materials,  and
! your use of  them is  governed by the  express license  under which  they were
! provided to you (License).  Unless the License provides otherwise, you may not
! use, modify, copy, publish, distribute,  disclose or transmit this software or
! the related documents without Intel's prior written permission.
!
! This software and the related documents  are provided as  is,  with no express
! or implied  warranties,  other  than those  that are  expressly stated  in the
! License.
!===============================================================================

!  Content:
!      Intel(R) Math Kernel Library (Intel(R) MKL) interface for Cluster DFT routines
!*******************************************************************************

! Include to build module MKL_DFTI
!INCLUDE 'mkl_dfti.f90'

! Definition of module MKL_CDFT_DM_TYPE. It is used just to define type DFTI_DESCRIPTOR_DM
MODULE MKL_CDFT_DM_TYPE
  use, intrinsic :: iso_fortran_env

  ! Definition of descriptor.
  ! Structure of this type is not used in Fortran code. The pointer to this type is used only
  TYPE DFTI_DESCRIPTOR_DM
     PRIVATE
     INTEGER(INT64) DESCRIPTOR
  END TYPE DFTI_DESCRIPTOR_DM

END MODULE MKL_CDFT_DM_TYPE

! Definition of module MKL_CDFT. It is used to define constants and interfaces of routines
MODULE MKL_CDFT
  use, intrinsic :: iso_fortran_env

  ! Module MKL_CDFT includes definitions from module MKL_DFTI and MKL_CDFT_DM_TYPE
  USE MKL_DFTI
  USE MKL_CDFT_DM_TYPE

  IMPLICIT NONE

  ! Codes of parameters for DftiGetValueDM / DftiSetValueDM
  INTEGER(INT64), PARAMETER :: CDFT_LOCAL_SIZE        =1000
  INTEGER(INT64), PARAMETER :: CDFT_LOCAL_X_START     =1001
  INTEGER(INT64), PARAMETER :: CDFT_LOCAL_NX          =1002
  INTEGER(INT64), PARAMETER :: CDFT_MPI_COMM          =1003
  INTEGER(INT64), PARAMETER :: CDFT_WORKSPACE         =1004
  INTEGER(INT64), PARAMETER :: CDFT_LOCAL_OUT_X_START =1005
  INTEGER(INT64), PARAMETER :: CDFT_LOCAL_OUT_NX      =1006

  ! Codes of errors
  INTEGER(INT64), PARAMETER :: CDFT_MPI_ERROR     =1000
  INTEGER(INT64), PARAMETER :: CDFT_SPREAD_ERROR  =1001

  ! Interfaces of routines
  INTERFACE DftiCreateDescriptorDM
     MODULE PROCEDURE DftiCreateDescriptorDM1
     MODULE PROCEDURE DftiCreateDescriptorDMn
     MODULE PROCEDURE DftiCreateDescriptorDM1_s
     MODULE PROCEDURE DftiCreateDescriptorDM1_d
     MODULE PROCEDURE DftiCreateDescriptorDMn_s
     MODULE PROCEDURE DftiCreateDescriptorDMn_d
  END INTERFACE DftiCreateDescriptorDM
  PRIVATE DftiCreateDescriptorDM1
  PRIVATE DftiCreateDescriptorDMn
  PRIVATE DftiCreateDescriptorDM1_s
  PRIVATE DftiCreateDescriptorDM1_d
  PRIVATE DftiCreateDescriptorDMn_s
  PRIVATE DftiCreateDescriptorDMn_d

  INTERFACE DftiGetValueDM
     MODULE PROCEDURE DftiGetValueDMs
     MODULE PROCEDURE DftiGetValueDMd
     MODULE PROCEDURE DftiGetValueDMi
  END INTERFACE DftiGetValueDM
  PRIVATE DftiGetValueDMs
  PRIVATE DftiGetValueDMd
  PRIVATE DftiGetValueDMi

  INTERFACE DftiSetValueDM
     MODULE PROCEDURE DftiSetValueDMs
     MODULE PROCEDURE DftiSetValueDMd
     MODULE PROCEDURE DftiSetValueDMi
     MODULE PROCEDURE DftiSetValueDMpc
     MODULE PROCEDURE DftiSetValueDMpz
     MODULE PROCEDURE DftiSetValueDMps
     MODULE PROCEDURE DftiSetValueDMpd
     MODULE PROCEDURE DftiSetValueDMpi
  END INTERFACE DftiSetValueDM
  PRIVATE DftiSetValueDMs
  PRIVATE DftiSetValueDMd
  PRIVATE DftiSetValueDMi
  PRIVATE DftiSetValueDMpc
  PRIVATE DftiSetValueDMpz
  PRIVATE DftiSetValueDMps
  PRIVATE DftiSetValueDMpd
  PRIVATE DftiSetValueDMpi

  INTERFACE DftiCommitDescriptorDM
     INTEGER(INT64) FUNCTION DftiCommitDescriptorDM_internal(H)
       USE MKL_CDFT_DM_TYPE
       TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
     END FUNCTION DftiCommitDescriptorDM_internal
  END INTERFACE DftiCommitDescriptorDM
  PRIVATE DftiCommitDescriptorDM_internal

  INTERFACE DftiComputeForwardDM
     MODULE PROCEDURE DftiComputeForwardDMos
     MODULE PROCEDURE DftiComputeForwardDMod
     MODULE PROCEDURE DftiComputeForwardDMoc
     MODULE PROCEDURE DftiComputeForwardDMoz
     MODULE PROCEDURE DftiComputeForwardDMis
     MODULE PROCEDURE DftiComputeForwardDMid
     MODULE PROCEDURE DftiComputeForwardDMic
     MODULE PROCEDURE DftiComputeForwardDMiz
  END INTERFACE DftiComputeForwardDM
  PRIVATE DftiComputeForwardDMos
  PRIVATE DftiComputeForwardDMod
  PRIVATE DftiComputeForwardDMoc
  PRIVATE DftiComputeForwardDMoz
  PRIVATE DftiComputeForwardDMis
  PRIVATE DftiComputeForwardDMid
  PRIVATE DftiComputeForwardDMic
  PRIVATE DftiComputeForwardDMiz

  INTERFACE DftiComputeBackwardDM
     MODULE PROCEDURE DftiComputeBackwardDMos
     MODULE PROCEDURE DftiComputeBackwardDMod
     MODULE PROCEDURE DftiComputeBackwardDMoc
     MODULE PROCEDURE DftiComputeBackwardDMoz
     MODULE PROCEDURE DftiComputeBackwardDMis
     MODULE PROCEDURE DftiComputeBackwardDMid
     MODULE PROCEDURE DftiComputeBackwardDMic
     MODULE PROCEDURE DftiComputeBackwardDMiz
  END INTERFACE DftiComputeBackwardDM
  PRIVATE DftiComputeBackwardDMos
  PRIVATE DftiComputeBackwardDMod
  PRIVATE DftiComputeBackwardDMoc
  PRIVATE DftiComputeBackwardDMoz
  PRIVATE DftiComputeBackwardDMis
  PRIVATE DftiComputeBackwardDMid
  PRIVATE DftiComputeBackwardDMic
  PRIVATE DftiComputeBackwardDMiz

  INTERFACE DftiFreeDescriptorDM
     INTEGER(INT64) FUNCTION DftiFreeDescriptorDM_internal(H)
       USE MKL_CDFT_DM_TYPE
       TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
     END FUNCTION DftiFreeDescriptorDM_internal
  END INTERFACE DftiFreeDescriptorDM
  PRIVATE DftiFreeDescriptorDM_internal

CONTAINS

  !INTERFACE DftiCreateDescriptorDM

  ! overloading of DftiCreateDescriptorDM for nD DFT
  INTEGER(INT64) FUNCTION DftiCreateDescriptorDMn(C,H,P1,P2,D,L)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    INTEGER(INT64) C
    INTEGER(INT64) P1,P2,D,L(*)
    INTENT(IN) :: C,P1,P2,D,L
    INTERFACE
       INTEGER(INT64) FUNCTION DftiCreateDescriptorDMn_internal(C,H,P1,P2,D,L)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         INTEGER(INT64) C
         INTEGER(INT64) P1,P2,D,L(*)
         INTENT(IN) :: C,P1,P2,D,L
       END FUNCTION DftiCreateDescriptorDMn_internal
    END INTERFACE
    DftiCreateDescriptorDMn = DftiCreateDescriptorDMn_internal(C,H,P1,P2,D,L)
  END FUNCTION DftiCreateDescriptorDMn

  ! overloading of DftiCreateDescriptorDM for 1D DFT
  INTEGER(INT64) FUNCTION DftiCreateDescriptorDM1(C,H,P1,P2,D,L)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    INTEGER(INT64) C
    INTEGER(INT64) P1,P2,D,L
    INTENT(IN) :: C,P1,P2,D,L
    INTERFACE
       INTEGER(INT64) FUNCTION DftiCreateDescriptorDM1_internal(C,H,P1,P2,D,L)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         INTEGER(INT64) C
         INTEGER(INT64) P1,P2,D,L
         INTENT(IN) :: C,P1,P2,D,L
       END FUNCTION DftiCreateDescriptorDM1_internal
    END INTERFACE
    DftiCreateDescriptorDM1 = DftiCreateDescriptorDM1_internal(C,H,P1,P2,D,L)
  END FUNCTION DftiCreateDescriptorDM1
  INTEGER(INT64) FUNCTION DftiCreateDescriptorDM1_s(C,H,P1R,P2,D,L)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    REAL(4), INTENT(IN) :: P1R
    INTEGER(INT64) C
    INTEGER(INT64) P1,P2,D,L
    INTENT(IN) :: C,P2,D,L
    P1 = INT(P1R)
    DftiCreateDescriptorDM1_s = DftiCreateDescriptorDM1(C,H,P1,P2,D,L)
  END FUNCTION DftiCreateDescriptorDM1_s
  INTEGER(INT64) FUNCTION DftiCreateDescriptorDM1_d(C,H,P1R,P2,D,L)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    REAL(8), INTENT(IN) :: P1R
    INTEGER(INT64) C
    INTEGER(INT64) P1,P2,D,L
    INTENT(IN) :: C,P2,D,L
    P1 = INT(P1R)
    DftiCreateDescriptorDM1_d = DftiCreateDescriptorDM1(C,H,P1,P2,D,L)
  END FUNCTION DftiCreateDescriptorDM1_d
  INTEGER(INT64) FUNCTION DftiCreateDescriptorDMn_s(C,H,P1R,P2,D,L)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    REAL(4), INTENT(IN) :: P1R
    INTEGER(INT64) C
    INTEGER(INT64) P1,P2,D,L(*)
    INTENT(IN) :: C,P2,D,L
    P1 = INT(P1R)
    DftiCreateDescriptorDMn_s = DftiCreateDescriptorDMn(C,H,P1,P2,D,L)
  END FUNCTION DftiCreateDescriptorDMn_s
  INTEGER(INT64) FUNCTION DftiCreateDescriptorDMn_d(C,H,P1R,P2,D,L)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    REAL(8), INTENT(IN) :: P1R
    INTEGER(INT64) C
    INTEGER(INT64) P1,P2,D,L(*)
    INTENT(IN) :: C,P2,D,L
    P1 = INT(P1R)
    DftiCreateDescriptorDMn_d = DftiCreateDescriptorDMn(C,H,P1,P2,D,L)
  END FUNCTION DftiCreateDescriptorDMn_d
  !END INTERFACE DftiCreateDescriptorDM
  !INTERFACE DftiSetValueDM

  ! overloading of DftiSetValueDM for SP value
  INTEGER(INT64) FUNCTION DftiSetValueDMs(H,P,V)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    INTEGER(INT64) P
    REAL(4) V
    INTENT(IN) :: P,V
    INTERFACE
       INTEGER(INT64) FUNCTION DftiSetValueDMf_internal(H,P,V)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         INTEGER(INT64) P
         REAL(4) V
         INTENT(IN) :: P,V
       END FUNCTION DftiSetValueDMf_internal
    END INTERFACE
    DftiSetValueDMs = DftiSetValueDMf_internal(H,P,V)
  END FUNCTION DftiSetValueDMs

  ! overloading of DftiSetValueDM for DP value
  INTEGER(INT64) FUNCTION DftiSetValueDMd(H,P,V)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    INTEGER(INT64) P
    REAL(8) V
    INTENT(IN) :: P,V
    INTERFACE
       INTEGER(INT64) FUNCTION DftiSetValueDMd_internal(H,P,V)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         INTEGER(INT64) P
         REAL(8) V
         INTENT(IN) :: P,V
       END FUNCTION DftiSetValueDMd_internal
    END INTERFACE
    DftiSetValueDMd = DftiSetValueDMd_internal(H,P,V)
  END FUNCTION DftiSetValueDMd

  ! overloading of DftiSetValueDM for integer(int64) value
  INTEGER(INT64) FUNCTION DftiSetValueDMi(H,P,V)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    INTEGER(INT64) P,V
    INTENT(IN) :: P,V
    INTERFACE
       INTEGER(INT64) FUNCTION DftiSetValueDMi_internal(H,P,V)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         INTEGER(INT64) P,V
         INTENT(IN) :: P,V
       END FUNCTION DftiSetValueDMi_internal
    END INTERFACE
    DftiSetValueDMi = DftiSetValueDMi_internal(H,P,V)
  END FUNCTION DftiSetValueDMi

  ! overloading of DftiSetValueDM for SP complex array
  INTEGER(INT64) FUNCTION DftiSetValueDMpc(H,P,V)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    INTEGER(INT64) P
    COMPLEX(4) V(*)
    INTENT(IN) :: P,V
    INTERFACE
       INTEGER(INT64) FUNCTION DftiSetValueDMp_internal(H,P,V)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         INTEGER(INT64) P
         COMPLEX(4) V(*)
         INTENT(IN) :: P,V
       END FUNCTION DftiSetValueDMp_internal
    END INTERFACE
    DftiSetValueDMpc = DftiSetValueDMp_internal(H,P,V)
  END FUNCTION DftiSetValueDMpc

  ! overloading of DftiSetValueDM for DP complex array
  INTEGER(INT64) FUNCTION DftiSetValueDMpz(H,P,V)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    INTEGER(INT64) P
    COMPLEX(8) V(*)
    INTENT(IN) :: P,V
    INTERFACE
       INTEGER(INT64) FUNCTION DftiSetValueDMp_internal(H,P,V)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         INTEGER(INT64) P
         COMPLEX(8) V(*)
         INTENT(IN) :: P,V
       END FUNCTION DftiSetValueDMp_internal
    END INTERFACE
    DftiSetValueDMpz = DftiSetValueDMp_internal(H,P,V)
  END FUNCTION DftiSetValueDMpz

  ! overloading of DftiSetValueDM for SP array
  INTEGER(INT64) FUNCTION DftiSetValueDMps(H,P,V)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    INTEGER(INT64) P
    REAL(4) V(*)
    INTENT(IN) :: P,V
    INTERFACE
       INTEGER(INT64) FUNCTION DftiSetValueDMp_internal(H,P,V)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         INTEGER(INT64) P
         REAL(4) V(*)
         INTENT(IN) :: P,V
       END FUNCTION DftiSetValueDMp_internal
    END INTERFACE
    DftiSetValueDMps = DftiSetValueDMp_internal(H,P,V)
  END FUNCTION DftiSetValueDMps

  ! overloading of DftiSetValueDM for DP array
  INTEGER(INT64) FUNCTION DftiSetValueDMpd(H,P,V)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    INTEGER(INT64) P
    REAL(8) V(*)
    INTENT(IN) :: P,V
    INTERFACE
       INTEGER(INT64) FUNCTION DftiSetValueDMp_internal(H,P,V)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         INTEGER(INT64) P
         REAL(8) V(*)
         INTENT(IN) :: P,V
       END FUNCTION DftiSetValueDMp_internal
    END INTERFACE
    DftiSetValueDMpd = DftiSetValueDMp_internal(H,P,V)
  END FUNCTION DftiSetValueDMpd

  ! overloading of DftiSetValueDM for integer(int64) array
  INTEGER(INT64) FUNCTION DftiSetValueDMpi(H,P,V)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    INTEGER(INT64) P,V(*)
    INTENT(IN) :: P,V
    INTERFACE
       INTEGER(INT64) FUNCTION DftiSetValueDMp_internal(H,P,V)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         INTEGER(INT64) P,V(*)
         INTENT(IN) :: P,V
       END FUNCTION DftiSetValueDMp_internal
    END INTERFACE
    DftiSetValueDMpi = DftiSetValueDMp_internal(H,P,V)
  END FUNCTION DftiSetValueDMpi
  !END INTERFACE DftiSetValueDM
  !INTERFACE DftiGetValueDM

  ! overloading of DftiGetValueDM for SP value
  INTEGER(INT64) FUNCTION DftiGetValueDMs(H,P,V)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    INTEGER(INT64) P
    REAL(4) V
    INTENT(IN)  :: P
    INTENT(OUT) :: V
    INTERFACE
       INTEGER(INT64) FUNCTION DftiGetValueDM_internal(H,P,V)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         INTEGER(INT64) P
         REAL(4) V
         INTENT(IN)  :: P
         INTENT(OUT) :: V
       END FUNCTION DftiGetValueDM_internal
    END INTERFACE
    DftiGetValueDMs = DftiGetValueDM_internal(H,P,V)
  END FUNCTION DftiGetValueDMs

  ! overloading of DftiGetValueDM for DP value
  INTEGER(INT64) FUNCTION DftiGetValueDMd(H,P,V)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    INTEGER(INT64) P
    REAL(8) V
    INTENT(IN)  :: P
    INTENT(OUT) :: V
    INTERFACE
       INTEGER(INT64) FUNCTION DftiGetValueDM_internal(H,P,V)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         INTEGER(INT64) P
         REAL(8) V
         INTENT(IN)  :: P
         INTENT(OUT) :: V
       END FUNCTION DftiGetValueDM_internal
    END INTERFACE
    DftiGetValueDMd = DftiGetValueDM_internal(H,P,V)
  END FUNCTION DftiGetValueDMd

  ! overloading of DftiGetValueDM for integer(int64) value
  INTEGER(INT64) FUNCTION DftiGetValueDMi(H,P,V)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    INTEGER(INT64) P,V
    INTENT(IN)  :: P
    INTENT(OUT) :: V
    INTERFACE
       INTEGER(INT64) FUNCTION DftiGetValueDM_internal(H,P,V)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         INTEGER(INT64) P,V
         INTENT(IN)  :: P
         INTENT(OUT) :: V
       END FUNCTION DftiGetValueDM_internal
    END INTERFACE
    DftiGetValueDMi = DftiGetValueDM_internal(H,P,V)
  END FUNCTION DftiGetValueDMi
  !END INTERFACE DftiGetValueDM
  !INTERFACE DftiComputeForwardDM

  ! overloading of DftiComputeForwardDM for SP R2C DFT (out-of-place)
  INTEGER(INT64) FUNCTION DftiComputeForwardDMos(H,IN,OUT)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    REAL(4) IN(*)
    COMPLEX(4) OUT(*)
    INTENT(IN)  :: IN
    INTENT(OUT) :: OUT
    INTERFACE
       INTEGER(INT64) FUNCTION DftiComputeForwardDMo_internal(H,IN,OUT)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         REAL(4) IN(*)
         COMPLEX(4) OUT(*)
         INTENT(IN)  :: IN
         INTENT(OUT) :: OUT
       END FUNCTION DftiComputeForwardDMo_internal
    END INTERFACE
    DftiComputeForwardDMos = DftiComputeForwardDMo_internal(H,IN,OUT)
  END FUNCTION DftiComputeForwardDMos

  ! overloading of DftiComputeForwardDM for DP R2C DFT (out-of-place)
  INTEGER(INT64) FUNCTION DftiComputeForwardDMod(H,IN,OUT)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    REAL(8) IN(*)
    COMPLEX(8) OUT(*)
    INTENT(IN)  :: IN
    INTENT(OUT) :: OUT
    INTERFACE
       INTEGER(INT64) FUNCTION DftiComputeForwardDMo_internal(H,IN,OUT)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         REAL(8) IN(*)
         COMPLEX(8) OUT(*)
         INTENT(IN)  :: IN
         INTENT(OUT) :: OUT
       END FUNCTION DftiComputeForwardDMo_internal
    END INTERFACE
    DftiComputeForwardDMod = DftiComputeForwardDMo_internal(H,IN,OUT)
  END FUNCTION DftiComputeForwardDMod

  ! overloading of DftiComputeForwardDM for SP C2C DFT (out-of-place)
  INTEGER(INT64) FUNCTION DftiComputeForwardDMoc(H,IN,OUT)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    COMPLEX(4) IN(*),OUT(*)
    INTENT(IN)  :: IN
    INTENT(OUT) :: OUT
    INTERFACE
       INTEGER(INT64) FUNCTION DftiComputeForwardDMo_internal(H,IN,OUT)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         COMPLEX(4) IN(*),OUT(*)
         INTENT(IN)  :: IN
         INTENT(OUT) :: OUT
       END FUNCTION DftiComputeForwardDMo_internal
    END INTERFACE
    DftiComputeForwardDMoc = DftiComputeForwardDMo_internal(H,IN,OUT)
  END FUNCTION DftiComputeForwardDMoc

  ! overloading of DftiComputeForwardDM for DP C2C DFT (out-of-place)
  INTEGER(INT64) FUNCTION DftiComputeForwardDMoz(H,IN,OUT)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    COMPLEX(8) IN(*),OUT(*)
    INTENT(IN)  :: IN
    INTENT(OUT) :: OUT
    INTERFACE
       INTEGER(INT64) FUNCTION DftiComputeForwardDMo_internal(H,IN,OUT)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         COMPLEX(8) IN(*),OUT(*)
         INTENT(IN)  :: IN
         INTENT(OUT) :: OUT
       END FUNCTION DftiComputeForwardDMo_internal
    END INTERFACE
    DftiComputeForwardDMoz = DftiComputeForwardDMo_internal(H,IN,OUT)
  END FUNCTION DftiComputeForwardDMoz

  ! overloading of DftiComputeForwardDM for SP R2C DFT (inplace)
  INTEGER(INT64) FUNCTION DftiComputeForwardDMis(H,IN)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    REAL(4) IN(*)
    INTERFACE
       INTEGER(INT64) FUNCTION DftiComputeForwardDMi_internal(H,IN)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         REAL(4) IN(*)
       END FUNCTION DftiComputeForwardDMi_internal
    END INTERFACE
    DftiComputeForwardDMis = DftiComputeForwardDMi_internal(H,IN)
  END FUNCTION DftiComputeForwardDMis

  ! overloading of DftiComputeForwardDM for DP R2C DFT (inplace)
  INTEGER(INT64) FUNCTION DftiComputeForwardDMid(H,IN)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    REAL(8) IN(*)
    INTERFACE
       INTEGER(INT64) FUNCTION DftiComputeForwardDMi_internal(H,IN)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         REAL(8) IN(*)
       END FUNCTION DftiComputeForwardDMi_internal
    END INTERFACE
    DftiComputeForwardDMid = DftiComputeForwardDMi_internal(H,IN)
  END FUNCTION DftiComputeForwardDMid

  ! overloading of DftiComputeForwardDM for SP C2C DFT (inplace)
  INTEGER(INT64) FUNCTION DftiComputeForwardDMic(H,IN)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    COMPLEX(4) IN(*)
    INTERFACE
       INTEGER(INT64) FUNCTION DftiComputeForwardDMi_internal(H,IN)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         COMPLEX(4) IN(*)
       END FUNCTION DftiComputeForwardDMi_internal
    END INTERFACE
    DftiComputeForwardDMic = DftiComputeForwardDMi_internal(H,IN)
  END FUNCTION DftiComputeForwardDMic

  ! overloading of DftiComputeForwardDM for DP C2C DFT (inplace)
  INTEGER(INT64) FUNCTION DftiComputeForwardDMiz(H,IN)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    COMPLEX(8) IN(*)
    INTERFACE
       INTEGER(INT64) FUNCTION DftiComputeForwardDMi_internal(H,IN)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         COMPLEX(8) IN(*)
       END FUNCTION DftiComputeForwardDMi_internal
    END INTERFACE
    DftiComputeForwardDMiz = DftiComputeForwardDMi_internal(H,IN)
  END FUNCTION DftiComputeForwardDMiz
  !END INTERFACE DftiComputeForwardDM
  !INTERFACE DftiComputeBackwardDM

  ! overloading of DftiComputeBackwardDM for SP R2C DFT (out-of-place)
  INTEGER(INT64) FUNCTION DftiComputeBackwardDMos(H,IN,OUT)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    COMPLEX(4) IN(*)
    REAL(4) OUT(*)
    INTENT(IN)  :: IN
    INTENT(OUT) :: OUT
    INTERFACE
       INTEGER(INT64) FUNCTION DftiComputeBackwardDMo_internal(H,IN,OUT)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         COMPLEX(4) IN(*)
         REAL(4) OUT(*)
         INTENT(IN)  :: IN
         INTENT(OUT) :: OUT
       END FUNCTION DftiComputeBackwardDMo_internal
    END INTERFACE
    DftiComputeBackwardDMos = DftiComputeBackwardDMo_internal(H,IN,OUT)
  END FUNCTION DftiComputeBackwardDMos

  ! overloading of DftiComputeBackwardDM for DP R2C DFT (out-of-place)
  INTEGER(INT64) FUNCTION DftiComputeBackwardDMod(H,IN,OUT)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    COMPLEX(8) IN(*)
    REAL(8) OUT(*)
    INTENT(IN)  :: IN
    INTENT(OUT) :: OUT
    INTERFACE
       INTEGER(INT64) FUNCTION DftiComputeBackwardDMo_internal(H,IN,OUT)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         COMPLEX(8) IN(*)
         REAL(8) OUT(*)
         INTENT(IN)  :: IN
         INTENT(OUT) :: OUT
       END FUNCTION DftiComputeBackwardDMo_internal
    END INTERFACE
    DftiComputeBackwardDMod = DftiComputeBackwardDMo_internal(H,IN,OUT)
  END FUNCTION DftiComputeBackwardDMod

  ! overloading of DftiComputeBackwardDM for SP C2C DFT (out-of-place)
  INTEGER(INT64) FUNCTION DftiComputeBackwardDMoc(H,IN,OUT)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    COMPLEX(4) IN(*),OUT(*)
    INTENT(IN)  :: IN
    INTENT(OUT) :: OUT
    INTERFACE
       INTEGER(INT64) FUNCTION DftiComputeBackwardDMo_internal(H,IN,OUT)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         COMPLEX(4) IN(*),OUT(*)
         INTENT(IN)  :: IN
         INTENT(OUT) :: OUT
       END FUNCTION DftiComputeBackwardDMo_internal
    END INTERFACE
    DftiComputeBackwardDMoc = DftiComputeBackwardDMo_internal(H,IN,OUT)
  END FUNCTION DftiComputeBackwardDMoc

  ! overloading of DftiComputeBackwardDM for DP C2C DFT (out-of-place)
  INTEGER(INT64) FUNCTION DftiComputeBackwardDMoz(H,IN,OUT)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    COMPLEX(8) IN(*),OUT(*)
    INTENT(IN)  :: IN
    INTENT(OUT) :: OUT
    INTERFACE
       INTEGER(INT64) FUNCTION DftiComputeBackwardDMo_internal(H,IN,OUT)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         COMPLEX(8) IN(*),OUT(*)
         INTENT(IN)  :: IN
         INTENT(OUT) :: OUT
       END FUNCTION DftiComputeBackwardDMo_internal
    END INTERFACE
    DftiComputeBackwardDMoz = DftiComputeBackwardDMo_internal(H,IN,OUT)
  END FUNCTION DftiComputeBackwardDMoz

  ! overloading of DftiComputeBackwardDM for SP R2C DFT (inplace)
  INTEGER(INT64) FUNCTION DftiComputeBackwardDMis(H,IN)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    REAL(4) IN(*)
    INTERFACE
       INTEGER(INT64) FUNCTION DftiComputeBackwardDMi_internal(H,IN)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         REAL(4) IN(*)
       END FUNCTION DftiComputeBackwardDMi_internal
    END INTERFACE
    DftiComputeBackwardDMis = DftiComputeBackwardDMi_internal(H,IN)
  END FUNCTION DftiComputeBackwardDMis

  ! overloading of DftiComputeBackwardDM for DP R2C DFT (inplace)
  INTEGER(INT64) FUNCTION DftiComputeBackwardDMid(H,IN)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    REAL(8) IN(*)
    INTERFACE
       INTEGER(INT64) FUNCTION DftiComputeBackwardDMi_internal(H,IN)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         REAL(8) IN(*)
       END FUNCTION DftiComputeBackwardDMi_internal
    END INTERFACE
    DftiComputeBackwardDMid = DftiComputeBackwardDMi_internal(H,IN)
  END FUNCTION DftiComputeBackwardDMid

  ! overloading of DftiComputeBackwardDM for SP C2C DFT (inplace)
  INTEGER(INT64) FUNCTION DftiComputeBackwardDMic(H,IN)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    COMPLEX(4) IN(*)
    INTERFACE
       INTEGER(INT64) FUNCTION DftiComputeBackwardDMi_internal(H,IN)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         COMPLEX(4) IN(*)
       END FUNCTION DftiComputeBackwardDMi_internal
    END INTERFACE
    DftiComputeBackwardDMic = DftiComputeBackwardDMi_internal(H,IN)
  END FUNCTION DftiComputeBackwardDMic

  ! overloading of DftiComputeBackwardDM for DP C2C DFT (inplace)
  INTEGER(INT64) FUNCTION DftiComputeBackwardDMiz(H,IN)
    USE MKL_CDFT_DM_TYPE
    TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    COMPLEX(8) IN(*)
    INTERFACE
       INTEGER(INT64) FUNCTION DftiComputeBackwardDMi_internal(H,IN)
         USE MKL_CDFT_DM_TYPE
         TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
         COMPLEX(8) IN(*)
       END FUNCTION DftiComputeBackwardDMi_internal
    END INTERFACE
    DftiComputeBackwardDMiz = DftiComputeBackwardDMi_internal(H,IN)
  END FUNCTION DftiComputeBackwardDMiz
  !END INTERFACE DftiComputeBackwardDM
END MODULE MKL_CDFT
