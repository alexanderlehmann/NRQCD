!----------------------------------------------------------------------
! MATRIXOPERATIONS, Matrix functions and routines
!----------------------------------------------------------------------
!
! MODULE: matrixoperations
!> @brief
!! Matrix functions and routines
!! @author
!! Alexander Lehmann,
!! UiS (<alexander.lehmann@uis.no>)
!! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
!! @date
!! 03.09.2018
!! @version
!! 1.0
! REVISION HISTORY:
! 03 09 2018
module matrixoperations
  USE,INTRINSIC :: ISO_FORTRAN_ENV ! defines kinds

  implicit none

  public

  interface ExpAH
     module procedure ExpAH_diaglibrary
     !module procedure ExpAH_mkl
  end interface ExpAH

  interface FrobeniusNorm
     module procedure FrobeniusNorm_cmplx64
     module procedure FrobeniusNorm_real64
  end interface FrobeniusNorm
contains
  !> @brief
  !! Frobenius-norm of a complex matrix
  !! @returns
  !! \f$\sqrt{\sum\limits_{i,j}|A_{ij}|^2}\f$
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 17.11.2018
  !! @version 1.0
  pure real(real64) function FrobeniusNorm_cmplx64(A)
    implicit none
    complex(real64), intent(in) :: A(:,:)

    FrobeniusNorm_cmplx64 = sqrt(sum(real(A*conjg(A),real64)))
  end function FrobeniusNorm_Cmplx64

  !> @brief
  !! Frobenius-norm of a real matrix
  !! @returns
  !! \f$\sqrt{\sum\limits_{i,j}|A_{ij}|^2}\f$
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 17.11.2018
  !! @version 1.0
  pure real(real64) function FrobeniusNorm_real64(A)
    implicit none
    real(real64), intent(in) :: A(:,:)

     FrobeniusNorm_Real64 = sqrt(sum(A**2))
  end function FrobeniusNorm_Real64
  
  !> @brief
  !! Computes Kronecker product
  !! @returns
  !! Kronecker product \f$A\otimes B\f$
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 03.09.2018
  !! @version 1.0
  pure function GetKronProd(A,B) result(AB)
    IMPLICIT NONE
    !> Matrix A
    complex(real64), dimension (:,:), intent(in) :: A
    !> Matrix B
    complex(real64), dimension (:,:), intent(in) :: B
    !> Kronecker product \f$A\otimes B\f$
    complex(real64), dimension (size(A,1)*size(B,1),size(A,2)*size(B,2)) :: AB
    
    integer R,RA,RB,C,CA,CB,I,J !Assistants.
    RA = UBOUND(A,DIM = 1)      !Ascertain the upper bounds of the incoming arrays.
    CA = UBOUND(A,DIM = 2)      !Their lower bounds will be deemed one,
    RB = UBOUND(B,DIM = 1)      !And the upper bound as reported will correspond.
    CB = UBOUND(B,DIM = 2)      !UBOUND(A) would give an array of two values, RA and CA, more for higher dimensionality.
    R = 0 !Syncopation: start the row offset.
    do I = 1,RA !Step down the rows of A.
       C = 0    !For each row, start the column offset.
       do J = 1,CA !Step along the columns of A.
          AB(R + 1:R + RB,C + 1:C + CB) = A(I,J)*B !Place a block of B values.
          C = C + CB !Advance a block of columns.
       end do !On to the next column of A.
       R = R + RB !Advance a block of rows.
    end do !On to the next row of A.
  end function GetKronProd

  !> @brief
  !! Computes trace
  !! @returns
  !! Trace \f$\Tr(A)\f$
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 03.09.2018
  !! @version 1.0
  pure function GetTrace(matrix) result(trace)
    implicit none
    !> Matrix of which the trace is to be computed
    complex(real64), intent(in) :: matrix(:,:)
    !> Trace
    complex(real64)             :: trace

    integer :: i

    trace = 0._real64
    do i=1,size(matrix,1)
       trace = trace + matrix(i,i)
    end do
  end function GetTrace
  
  !> @brief
  !! Returns a unit matrix of size \f$n\f$
  !! @returns
  !! Unit matrix of size \f$n\f$
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 03.09.2018
  !! @version 1.0
  pure function GetUnitmatrix(n) result(res)
    implicit none
    !> Size of the unit matrix
    integer, intent(in) :: n
    !> Unit matrix
    complex(real64)           :: res(n,n)

    integer :: i

    res = 0._real64
    forall(i=1:n) res(i,i) = 1._real64
  end function GetUnitmatrix

  !> @brief
  !! Computes hermitian conjugate
  !! @returns
  !! Hermitian conjugate
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 03.09.2018
  !! @version 1.0
  pure function GetHerm(matrix) result(res)
    implicit none
    !> Matrix
    complex(real64), dimension(:,:),               intent(in) :: matrix
    !> Hermitian conjugate
    complex(real64), dimension(size(matrix,1),size(matrix,2)) :: res

    res = conjg(transpose(matrix))
  end function GetHerm

  !> @brief
  !! Performs eigen-system decomposition of hermitian matrix using the Jacobi algorithm
  !! @details
  !! Unitary transformation matrix \f$U\f$ fulfills
  !! \f$d=U\cdot A\cdot U^\dagger\Leftrightarrow A=U^\dagger\cdot d\cdot U\Leftrightarrow U\cdot A=d\cdot U
  !! @author
  !! Taken from diag library (then modernized),
  !! adapted from Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202\n
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 03.09.2018
  !! @version 1.0
  pure subroutine EigenH(A_, d, U, sort)
    implicit none
    !> Hermitian input \f$A\f$
    complex(real64), intent(in)  :: A_(:,:)
    !> Eigen values
    real(real64),    intent(out) :: d(size(A_,1))
    !> Unitary transform which diagonalizes \f$A\f$
    complex(real64), intent(out), optional :: U(size(A_,1),size(A_,1))
    !> If given, sorting of eigen-values performed according to \f$d_i>d_{i+1}\forall i\f$
    integer,  intent(in),  optional :: sort

    complex(real64)  :: A(size(A_,1),size(A_,1))
    integer    :: n,p, q, j
    Real(real64)     :: red, off, thresh
    Real(real64)     :: t, delta, invc, s
    Complex(real64)  :: x, y, Apq
    Real(real64)     :: ev(size(A_,1),2)

    integer :: sweep
    REAL(REAL64), parameter :: SYM_EPS=2._real64**(-103)

    n = size(A_,1)

    A = A_
    ev = 0
    forall(p=1:n) ev(p,2) = real(A(p,p),real64)
    d = ev(:,2)

    U = GetUnitmatrix(n)

    red = .04_real64/n**4

    do sweep = 1, 50
       off = 0._real64
       do q = 2, n
          do p = 1, q - 1
             off = off + real(A(p,q)*Conjg(A(p,q)),real64)
          end do
       end do
       if(  off .le. SYM_EPS ) exit
       
       thresh = 0._real64
       if( sweep .lt. 4 ) thresh = off*red

       do q = 2, n
          do p = 1, q - 1
             Apq = A(p,q)
             off = real(Apq*Conjg(Apq),real64)
             if( (sweep .gt. 4) .and. (off .lt. SYM_EPS*(ev(p,2)**2 + ev(q,2)**2)) ) then
                A(p,q) = 0._real64
             else if( off .gt. thresh ) then
                t = .5_real64*(ev(p,2) - ev(q,2))
                t = 1._real64/(t + sign(sqrt(t**2 + off), t))

                delta = t*off
                ev(p,1) = ev(p,1) + delta
                ev(p,2) = d(p) + ev(p,1)
                ev(q,1) = ev(q,1) - delta
                ev(q,2) = d(q) + ev(q,1)

                invc = sqrt(delta*t + 1._real64)
                s = t/invc
                t = delta/(invc + 1._real64)

                do j = 1, p - 1
                   x = A(j,p)
                   y = A(j,q)
                   A(j,p) = x + s*(Conjg(Apq)*y - t*x)
                   A(j,q) = y - s*(Apq*x + t*y)
                end do

                do j = p + 1, q - 1
                   x = A(p,j)
                   y = A(j,q)
                   A(p,j) = x + s*(Apq*Conjg(y) - t*x)
                   A(j,q) = y - s*(Apq*Conjg(x) + t*y)
                end do

                do j = q + 1, n
                   x = A(p,j)
                   y = A(q,j)
                   A(p,j) = x + s*(Apq*y - t*x)
                   A(q,j) = y - s*(Conjg(Apq)*x + t*y)
                end do

                A(p,q) = 0._real64

                do j = 1, n
                   x = U(p,j)
                   y = U(q,j)
                   U(p,j) = x + s*(Apq*y - t*x)
                   U(q,j) = y - s*(Conjg(Apq)*x + t*y)
                end do
             end if
          end do
       end do

       ev(:,1) = 0
       d = ev(:,2)
    end do

    if( present(sort) ) then
       do p = 1, n - 1
          j = p
          t = d(p)
          do q = p + 1, n
             if( (t - d(q)) .gt. 0 ) then
                j = q
                t = d(q)
             end if
          end do

          if( j .ne. p ) then
             d(j) = d(p)
             d(p) = t
             do q = 1, n
                x = U(p,q)
                U(p,q) = U(j,q)
                U(j,q) = x
             end do
          end if
       end do
    end if
  end subroutine EigenH

  !> @brief
  !! Performs eigen-system decomposition of general complex matrix using the Jacobi algorithm
  !! @details
  !! Unitary transformation matrix \f$U\f$ fulfills
  !! \f$d=U\cdot A\cdot U^\dagger\Leftrightarrow A=U^\dagger\cdot d\cdot U\Leftrightarrow U\cdot A=d\cdot U
  !! @author
  !! Taken from diag library (then modernized),
  !! adapted from Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202\n
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 03.09.2018
  !! @version 1.0
  pure subroutine EigenC(A_, d, U, sort)
    implicit none
    !> Input \f$A\f$
    complex(real64), intent(in)  :: A_(:,:)
    !> Eigen values
    complex(real64), intent(out) :: d(size(A_,1))
    !> Unitary transform which diagonalizes \f$A\f$
    complex(real64), intent(out), optional :: U(size(A_,1),size(A_,1))
    !> If given, sorting of eigen-values performed according to \f$d_i>d_{i+1}\forall i\f$
    integer,  intent(in),  optional :: sort

    complex(real64)  :: A(size(A_,1),size(A_,1))

    integer  :: p, q, j, n
    Real(real64)    :: red, off, thresh, norm
    Complex(real64) :: delta, t, s, invc, sx, sy, tx, ty
    Complex(real64) :: x, y
    Complex(real64) :: ev(size(A_,1),2)
    REAL(REAL64), parameter :: EPS=2._real64**(-103)
    integer :: sweep

    n = size(A_,1)
    A = A_

    ev=0
    forall(p=1:n) ev(p,2)=A(p,p)
    d= ev(:,2)

    U = GetUnitmatrix(n)

    red = .01_real64/n**4

    do sweep = 1, 50
       off = 0
       do q = 2, n
          do p = 1, q - 1
             off = off + abs(A(p,q))**2 + abs(A(q,p))**2
          enddo
       enddo
       if( off .le. EPS ) exit

       thresh = 0
       if( sweep .lt. 4 ) thresh = off*red

       do q = 2, n
          do p = 1, q - 1
             off = abs(A(p,q))**2 + abs(A(q,p))**2
             if( sweep .gt. 4 .and. off .lt. &
                  EPS*(abs(ev(p,2))**2 + abs(ev(q,2))**2) ) then
                A(p,q) = 0
                A(q,p) = 0
             else if( off .gt. thresh ) then
                delta = A(p,q)*A(q,p)
                x = .5_real64*(ev(p,2) - ev(q,2))
                y = sqrt(x**2 + delta)
                t = x - y
                s = x + y
                if( abs(t) .lt. abs(s) ) t = s

                t = 1/t
                delta = delta*t
                ev(p,1) = ev(p,1) + delta
                ev(p,2) = d(p) + ev(p,1)
                ev(q,1) = ev(q,1) - delta
                ev(q,2) = d(q) + ev(q,1)

                invc = sqrt(delta*t + 1)
                s = t/invc
                t = t/(invc + 1)
                sx = s*A(p,q)
                ty = t*A(p,q)
                sy = s*A(q,p)
                tx = t*A(q,p)

                do j = 1, n
                   x = A(j,p)
                   y = A(j,q)
                   A(j,p) = x + sy*(y - ty*x)
                   A(j,q) = y - sx*(x + tx*y)
                   x = A(p,j)
                   y = A(q,j)
                   A(p,j) = x + sx*(y - tx*x)
                   A(q,j) = y - sy*(x + ty*y)
                end do
                A(p,q) = 0
                A(q,p) = 0
                do j = 1, n
                   x = U(p,j)
                   y = U(q,j)
                   U(p,j) = x + sx*(y - tx*x)
                   U(q,j) = y - sy*(x + ty*y)
                end do
             end if
          end do
       end do

       ev(:,1) = 0
       d = ev(:,2)
    end do

    ! normalize the eigenvectors
    do p = 1, n
       norm = 0
       do q = 1, n
          norm = norm + abs(U(p,q))**2
       end do
       norm = 1/sqrt(norm)
       U(p,:) = U(p,:)*norm
    end do
    if( present(sort) ) then
       ! sort the eigenvalues by their real part
       do p = 1, n - 1
          j = p
          t = d(p)
          do q = p + 1, n
             if( (Real(t) - Real(d(q))) .gt. 0 ) then
                j = q
                t = d(q)
             end if
          end do

          if( j .ne. p ) then
             d(j) = d(p)
             d(p) = t
             do q = 1, n
                x = U(p,q)
                U(p,q) = U(j,q)
                U(j,q) = x
             end do
          end if
       end do
    end if
  end subroutine EigenC

  !> @brief
  !! Computes the logarithm of a U(n) or SU(n) group element, returning the algebra element
  !! @details
  !! Via diagonalisation of the input, the complex logarithm into the first branch
  !! is taken and then the back-wards-similarity transform performed (diagonalisation reverted)
  !! @returns
  !! \f$H=\log(U)/\text{i}\f$
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 03.09.2018
  !! @version 1.0
  pure function ArgLogU(U) result(res)
    implicit none
    !> Group element
    complex(real64), dimension(:,:),     intent(in) :: U
    !> Algebra element
    complex(real64), dimension(size(U,1),size(U,2)) :: res

    complex(real64), dimension(size(U,1)) :: ev
    complex(real64), dimension(size(U,1),size(U,2)) :: S

    integer :: i

    ! Eigen-decomposition with U = S^+ ev S
    call EigenC(U,ev,S)

    ! Complex logarithm via atan2
    res = 0._real64
    forall(i=1:size(ev)) res(i,i) = ATAN2(aimag(ev(i)),real(ev(i),real64))

    ! Multiply similarity transform matrices to the result
    res = matmul(matmul(GetHerm(S),res),S)
  end function ArgLogU
  
  !> @brief
  !! Computes the logarithm of a U(n) or SU(n) group element
  !! @details
  !! Via diagonalisation of the input, the complex logarithm into the first branch
  !! is taken and then the back-wards-similarity transform performed (diagonalisation reverted)
  !! @returns
  !! \f$H=\log(U)\f$
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 03.09.2018
  !! @version 1.0
  pure function LogU(U) result(res)
    implicit none
    !> Group element
    complex(real64), dimension(:,:),     intent(in) :: U
    !> Algebra element
    complex(real64), dimension(size(U,1),size(U,2)) :: res

    complex(real64), dimension(size(U,1)) :: ev
    complex(real64), dimension(size(U,1),size(U,2)) :: S

    integer :: i

    ! Eigen-decomposition with U = S^+ ev S
    call EigenC(U,ev,S)

    ! Complex logarithm via atan2
    res = 0._real64
    forall(i=1:size(ev)) res(i,i) = cmplx(0,ATAN2(aimag(ev(i)),real(ev(i),real64)),real64)

    ! Multiply similarity transform matrices to the result
    res = matmul(matmul(GetHerm(S),res),S)
  end function LogU
  
  !> @brief
  !! Matrix exponential of skew-hermitian input
  !! @details
  !! \f$ \exp(A) = \exp( \text{i}S\cdot (-\text{i})A\cdot S^\dagger) = S\cdot \exp(\text{i} D )\cdot S^\dagger\f$
  !! First computes eigenvalues and eigenvectors, thus obtaining \f$S\f$, then computes exponential
  !! of the eigen values and then in a final step reverts the similarity transform
  !! @returns
  !! \f$H=\exp(A)\f$
  !! @author
  !! Alexander Lehmann,
  !! UiS (<alexander.lehmann@uis.no>)
  !! and ITP Heidelberg (<lehmann@thpys.uni-heidelberg.de>)
  !! @date 03.09.2018
  !! @version 1.0
  pure function ExpAH_diaglibrary(A_) result(res)
    implicit none
    !> Skew-hermitian input A
    complex(real64), dimension(:,:), intent(in) :: A_
    !> Exponential of A
    complex(real64), dimension(size(A_,1),size(A_,2)) :: res

    real(real64),    dimension(size(A_,1)) :: ev
    complex(real64), dimension(size(A_,1),size(A_,2)) :: s
    integer :: i

    ! Eigen-decomposition with U = S^+ ev S
    call EigenH(A_*cmplx(0._real64,-1._real64,real64),ev,S)

    ! Exponential
    res = 0._real64
    forall(i=1:size(A_,1)) res(i,i) = exp(cmplx(0._real64,ev(i),real64))

    ! Multiply similarity transform matrices to the result
    res = matmul(matmul(GetHerm(S),res),S)
  end function ExpAH_Diaglibrary

  !function ExpAH_mkl(A_) result(res)
  !  implicit none
  !  !> Skew-hermitian input A
  !  complex(real64), dimension(:,:), intent(in) :: A_
  !  !> Exponential of A
  !  complex(real64), dimension(size(A_,1),size(A_,2)) :: res
  !
  !  !.. External Subroutines ..
  !  EXTERNAL ZHEEVD, ZCOPY, ZSCAL, ZGEMM
  !
  !  !.. Parameters ..
  !  INTEGER N
  !  INTEGER LDA
  !  INTEGER, parameter :: LWMAX = 1000
  !  
  !  !.. Local Scalars ..
  !  INTEGER          INFO, LWORK, LIWORK, LRWORK
  !  
  !  !.. Local Arrays ..
  !  INTEGER          IWORK( LWMAX )
  !  REAL(REAL64),    allocatable :: ev( : )
  !  REAL(REAL64) :: RWORK( LWMAX )
  !  COMPLEX(REAL64), allocatable :: A( :,: )
  !  complex(real64) :: WORK( LWMAX )
  !  
  !  complex(real64), allocatable :: S(:,:)
  !
  !  integer :: i
  !  
  !  !--**START**-- SIZE PARAMETERS
  !  N = size(A_,1)
  !  LDA=N
  !  
  !  LWORK  = 2*N + N**2
  !  LIWORK = 3 + 5*N
  !  LRWORK = 1 + 5*N + 2*N**2
  !  !--** END **-- SIZE PARAMETERS
  !
  !  !--**START**-- Allocation
  !  allocate(ev(N))
  !  allocate(A( LDA, N ))
  !  !--** END **-- Allocation
  !
  !  !--**START**-- Transforming skew-hermitian input to hermitian matrix:A -> -i.A
  !  A = A_*cmplx(0._real64,-1._real64,real64)
  !  !--** END **-- Transforming skew-hermitian input to hermitian matrix:A -> -i.A
  !
  !  ! Diagonalisation for hermitian matrix -i.A
  !  call zheevd('V','L',N,A,LDA,ev, WORK, LWORK, RWORK,LRWORK, IWORK, LIWORK, INFO )
  !
  !  ! Exponentiation of eigenvalues of A = i.(-i.A)
  !  res = 0._real64
  !  forall(i=1:N) res(i,i) = exp(cmplx(0._real64,ev(i),real64))
  !
  !  S = A
  !  !--**START**-- Inverse similar transform: exp(A) = exp(S.A.S†) = S.exp(A).S†
  !  res = matmul(matmul(S,res),GetHerm(S))
  !  !--** END **-- Inverse similar transform: exp(A) = exp(S.A.S†) = S.exp(A).S†
  !end function ExpAH_mkl
end module matrixoperations
