   ! This routine provides an interface between an old LINPACK routines
   ! and the LAPACK routines.  
   
   ! Currently, the following Linpack routines are mapped
   ! LINPACK -> LAPACK
   !  DGBFA -> DGBTRF
   !  DGBSL -> DGBTRS
   !  DGEFA -> DGETRF
   !  DGESL -> DGETRS
   !  DPOFA -> dpotrf
   !  dtrsl -> dtrtrs

   ! This uses LAPACK DGETRS to replace LINPACK DGESL
   SUBROUTINE DGESL (A, LDA, N, IPVT, B, JOB)
      implicit none
      INTEGER LDA,N,IPVT(N),JOB
      real(8) A(LDA,N),B(N)
   !***BEGIN PROLOGUE  DGESL
   !***PURPOSE  Solve the real system A*X=B or TRANS(A)*X=B using the
   !            factors computed by DGECO or DGEFA.
   !***LIBRARY   SLATEC (LINPACK)
   !***CATEGORY  D2A1
   !***TYPE      DOUBLE PRECISION (SGESL-S, DGESL-D, CGESL-C)
   !***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
   !***AUTHOR  Moler, C. B., (U. of New Mexico)
   !***DESCRIPTION
   !
   !     DGESL solves the double precision system
   !     A * X = B  or  TRANS(A) * X = B
   !     using the factors computed by DGECO or DGEFA.
   !
   !     On Entry
   !
   !        A       DOUBLE PRECISION(LDA, N)
   !                the output from DGECO or DGEFA.
   !
   !        LDA     INTEGER
   !                the leading dimension of the array  A .
   !
   !        N       INTEGER
   !                the order of the matrix  A .
   !
   !        IPVT    INTEGER(N)
   !                the pivot vector from DGECO or DGEFA.
   !
   !        B       DOUBLE PRECISION(N)
   !                the right hand side vector.
   !
   !        JOB     INTEGER
   !                = 0         to solve  A*X = B ,
   !                = nonzero   to solve  TRANS(A)*X = B  where
   !                            TRANS(A)  is the transpose.
   !
   !     On Return
   !
   !        B       the solution vector  X .
   !
   !     Error Condition
   !
   !        A division by zero will occur if the input factor contains a
   !        zero on the diagonal.  Technically this indicates singularity
   !        but it is often caused by improper arguments or improper
   !        setting of LDA .  It will not occur if the subroutines are
   !        called correctly and if DGECO has set RCOND .GT. 0.0
   !        or DGEFA has set INFO .EQ. 0 .
   !
   !     To compute  INVERSE(A) * C  where  C  is a matrix
   !     with  P  columns
   !           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
   !           IF (RCOND is too small) GO TO ...
   !           DO 10 J = 1, P
   !              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
   !        10 CONTINUE
   !
   !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
   !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
   !***ROUTINES CALLED  DAXPY, DDOT
   !***REVISION HISTORY  (YYMMDD)
   !   780814  DATE WRITTEN
   !   890831  Modified array declarations.  (WRB)
   !   890831  REVISION DATE from Version 3.2
   !   891214  Prologue converted to Version 4.0 format.  (BAB)
   !   900326  Removed duplicate information from DESCRIPTION section.
   !           (WRB)
   !   920501  Reformatted the REFERENCES section.  (WRB)
   !***END PROLOGUE  DGESL
      
   !      SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
   !
   !  -- LAPACK routine (version 3.1) --
   !     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
   !     November 2006
   !
   !  Purpose
   !  =======
   !
   !  DGETRS solves a system of linear equations
   !     A * X = B  or  A' * X = B
   !  with a general N-by-N matrix A using the LU factorization computed
   !  by DGETRF.
   !
   !  Arguments
   !  =========
   !
   !  TRANS   (input) CHARACTER*1
   !          Specifies the form of the system of equations:
   !          = 'N':  A * X = B  (No transpose)
   !          = 'T':  A'* X = B  (Transpose)
   !          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
   !
   !  N       (input) INTEGER
   !          The order of the matrix A.  N >= 0.
   !
   !  NRHS    (input) INTEGER
   !          The number of right hand sides, i.e., the number of columns
   !          of the matrix B.  NRHS >= 0.
   !
   !  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
   !          The factors L and U from the factorization A = P*L*U
   !          as computed by DGETRF.
   !
   !  LDA     (input) INTEGER
   !          The leading dimension of the array A.  LDA >= max(1,N).
   !
   !  IPIV    (input) INTEGER array, dimension (N)
   !          The pivot indices from DGETRF; for 1<=i<=N, row i of the
   !          matrix was interchanged with row IPIV(i).
   !
   !  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
   !          On entry, the right hand side matrix B.
   !          On exit, the solution matrix X.
   !
   !  LDB     (input) INTEGER
   !          The leading dimension of the array B.  LDB >= max(1,N).
   !
   !  INFO    (output) INTEGER
   !          = 0:  successful exit
   !          < 0:  if INFO = -i, the i-th argument had an illegal value
   !
   !  =====================================================================

      CHARACTER          TRANS
      INTEGER            LDB, NRHS, INFO

      trans = 'N'
      if (job /= 0) trans = 'T'
      LDB = N
      NRHS = 1
      call DGETRS( TRANS, N, NRHS, A, LDA, IPVT, B, LDB, INFO )
      if (info < 0) then
         write(6,*) 'Error in DGETRS call via DGESL, INFO = ', info
         stop
      end if
      return
   end subroutine dgesl

   ! This uses LAPACK DGBTRF to replace LINPACK DGBFA
   SUBROUTINE DGBFA (ABD, LDA, N, ML, MU, IPVT, INFO)
      implicit none
      INTEGER LDA,N,ML,MU,IPVT(N),INFO
      real(8) ABD(LDA,N)
   !***BEGIN PROLOGUE  DGBFA
   !***PURPOSE  Factor a band matrix using Gaussian elimination.
   !***LIBRARY   SLATEC (LINPACK)
   !***CATEGORY  D2A2
   !***TYPE      DOUBLE PRECISION (SGBFA-S, DGBFA-D, CGBFA-C)
   !***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION
   !***AUTHOR  Moler, C. B., (U. of New Mexico)
   !***DESCRIPTION
   !
   !     DGBFA factors a double precision band matrix by elimination.
   !
   !     DGBFA is usually called by DGBCO, but it can be called
   !     directly with a saving in time if  RCOND  is not needed.
   !
   !     On Entry
   !
   !        ABD     DOUBLE PRECISION(LDA, N)
   !                contains the matrix in band storage.  The columns
   !                of the matrix are stored in the columns of  ABD  and
   !                the diagonals of the matrix are stored in rows
   !                ML+1 through 2*ML+MU+1 of  ABD .
   !                See the comments below for details.
   !
   !        LDA     INTEGER
   !                the leading dimension of the array  ABD .
   !                LDA must be .GE. 2*ML + MU + 1 .
   !
   !        N       INTEGER
   !                the order of the original matrix.
   !
   !        ML      INTEGER
   !                number of diagonals below the main diagonal.
   !                0 .LE. ML .LT.  N .
   !
   !        MU      INTEGER
   !                number of diagonals above the main diagonal.
   !                0 .LE. MU .LT.  N .
   !                More efficient if  ML .LE. MU .
   !     On Return
   !
   !        ABD     an upper triangular matrix in band storage and
   !                the multipliers which were used to obtain it.
   !                The factorization can be written  A = L*U  where
   !                L  is a product of permutation and unit lower
   !                triangular matrices and  U  is upper triangular.
   !
   !        IPVT    INTEGER(N)
   !                an integer vector of pivot indices.
   !
   !        INFO    INTEGER
   !                = 0  normal value.
   !                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
   !                     condition for this subroutine, but it does
   !                     indicate that DGBSL will divide by zero if
   !                     called.  Use  RCOND  in DGBCO for a reliable
   !                     indication of singularity.
   !
   !     Band Storage
   !
   !           If  A  is a band matrix, the following program segment
   !           will set up the input.
   !
   !                   ML = (band width below the diagonal)
   !                   MU = (band width above the diagonal)
   !                   M = ML + MU + 1
   !                   DO 20 J = 1, N
   !                      I1 = MAX(1, J-MU)
   !                      I2 = MIN(N, J+ML)
   !                      DO 10 I = I1, I2
   !                         K = I - J + M
   !                         ABD(K,J) = A(I,J)
   !                10    CONTINUE
   !                20 CONTINUE
   !
   !           This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
   !           In addition, the first  ML  rows in  ABD  are used for
   !           elements generated during the triangularization.
   !           The total number of rows needed in  ABD  is  2*ML+MU+1 .
   !           The  ML+MU by ML+MU  upper left triangle and the
   !           ML by ML  lower right triangle are not referenced.
   !
   !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
   !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
   !***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX
   !***REVISION HISTORY  (YYMMDD)
   !   780814  DATE WRITTEN
   !   890531  Changed all specific intrinsics to generic.  (WRB)
   !   890831  Modified array declarations.  (WRB)
   !   890831  REVISION DATE from Version 3.2
   !   891214  Prologue converted to Version 4.0 format.  (BAB)
   !   900326  Removed duplicate information from DESCRIPTION section.
   !           (WRB)
   !   920501  Reformatted the REFERENCES section.  (WRB)
   !***END PROLOGUE  DGBFA
   !
   !====================================================================      
   !      SUBROUTINE DGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
   !
   !  -- LAPACK routine (version 3.1) --
   !     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
   !     November 2006
   !
   !     .. Scalar Arguments ..
   !      INTEGER            INFO, KL, KU, LDAB, M, N
   !     ..
   !     .. Array Arguments ..
   !      INTEGER            IPIV( * )
   !      DOUBLE PRECISION   AB( LDAB, * )
   !     ..
   !
   !  Purpose
   !  =======
   !
   !  DGBTRF computes an LU factorization of a real m-by-n band matrix A
   !  using partial pivoting with row interchanges.
   !
   !  This is the blocked version of the algorithm, calling Level 3 BLAS.
   !
   !  Arguments
   !  =========
   !
   !  M       (input) INTEGER
   !          The number of rows of the matrix A.  M >= 0.
   !
   !  N       (input) INTEGER
   !          The number of columns of the matrix A.  N >= 0.
   !
   !  KL      (input) INTEGER
   !          The number of subdiagonals within the band of A.  KL >= 0.
   !
   !  KU      (input) INTEGER
   !          The number of superdiagonals within the band of A.  KU >= 0.
   !
   !  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
   !          On entry, the matrix A in band storage, in rows KL+1 to
   !          2*KL+KU+1; rows 1 to KL of the array need not be set.
   !          The j-th column of A is stored in the j-th column of the
   !          array AB as follows:
   !          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
   !
   !          On exit, details of the factorization: U is stored as an
   !          upper triangular band matrix with KL+KU superdiagonals in
   !          rows 1 to KL+KU+1, and the multipliers used during the
   !          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
   !          See below for further details.
   !
   !  LDAB    (input) INTEGER
   !          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
   !
   !  IPIV    (output) INTEGER array, dimension (min(M,N))
   !          The pivot indices; for 1 <= i <= min(M,N), row i of the
   !          matrix was interchanged with row IPIV(i).
   !
   !  INFO    (output) INTEGER
   !          = 0: successful exit
   !          < 0: if INFO = -i, the i-th argument had an illegal value
   !          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
   !               has been completed, but the factor U is exactly
   !               singular, and division by zero will occur if it is used
   !               to solve a system of equations.
   !
   !  Further Details
   !  ===============
   !
   !  The band storage scheme is illustrated by the following example, when
   !  M = N = 6, KL = 2, KU = 1:
   !
   !  On entry:                       On exit:
   !
   !      *    *    *    +    +    +       *    *    *   u14  u25  u36
   !      *    *    +    +    +    +       *    *   u13  u24  u35  u46
   !      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
   !     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
   !     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
   !     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
   !
   !  Array elements marked * are not used by the routine; elements marked
   !  + need not be set on entry, but are required by the routine to store
   !  elements of U because of fill-in resulting from the row interchanges.
   !
   !  =====================================================================
      ! Since call to DGBFA does not specify M, use N as the row dim of A,
      ! as I belive A is square in linpack routine.
      integer m
      m = n
      call DGBTRF( M, N, ML, MU, ABD, LDA, IPVT, INFO )  
      if (info < 0) then
         write(6,*) 'Error in DGBTRF call via DGBFA, INFO = ', info
         stop
      end if
      return
   end subroutine DGBFA

   ! This uses LAPACK DGBTRS to replace LINPACK DGBSL
   SUBROUTINE DGBSL (ABD, LDA, N, ML, MU, IPVT, B, JOB)
      implicit none
      INTEGER LDA,N,ML,MU,IPVT(N),JOB
      real(8) ABD(LDA,N),B(N)
   !***BEGIN PROLOGUE  DGBSL
   !***PURPOSE  Solve the real band system A*X=B or TRANS(A)*X=B using
   !            the factors computed by DGBCO or DGBFA.
   !***LIBRARY   SLATEC (LINPACK)
   !***CATEGORY  D2A2
   !***TYPE      DOUBLE PRECISION (SGBSL-S, DGBSL-D, CGBSL-C)
   !***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
   !***AUTHOR  Moler, C. B., (U. of New Mexico)
   !***DESCRIPTION
   !
   !     DGBSL solves the double precision band system
   !     A * X = B  or  TRANS(A) * X = B
   !     using the factors computed by DGBCO or DGBFA.
   !
   !     On Entry
   !
   !        ABD     DOUBLE PRECISION(LDA, N)
   !                the output from DGBCO or DGBFA.
   !
   !        LDA     INTEGER
   !                the leading dimension of the array  ABD .
   !
   !        N       INTEGER
   !                the order of the original matrix.
   !
   !        ML      INTEGER
   !                number of diagonals below the main diagonal.
   !
   !        MU      INTEGER
   !                number of diagonals above the main diagonal.
   !
   !        IPVT    INTEGER(N)
   !                the pivot vector from DGBCO or DGBFA.
   !
   !        B       DOUBLE PRECISION(N)
   !                the right hand side vector.
   !
   !        JOB     INTEGER
   !                = 0         to solve  A*X = B ,
   !                = nonzero   to solve  TRANS(A)*X = B , where
   !                            TRANS(A)  is the transpose.
   !
   !     On Return
   !
   !        B       the solution vector  X .
   !
   !     Error Condition
   !
   !        A division by zero will occur if the input factor contains a
   !        zero on the diagonal.  Technically this indicates singularity
   !        but it is often caused by improper arguments or improper
   !        setting of LDA .  It will not occur if the subroutines are
   !        called correctly and if DGBCO has set RCOND .GT. 0.0
   !        or DGBFA has set INFO .EQ. 0 .
   !
   !     To compute  INVERSE(A) * C  where  C  is a matrix
   !     with  P  columns
   !           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
   !           IF (RCOND is too small) GO TO ...
   !           DO 10 J = 1, P
   !              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
   !        10 CONTINUE
   !
   !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
   !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
   !***ROUTINES CALLED  DAXPY, DDOT
   !***REVISION HISTORY  (YYMMDD)
   !   780814  DATE WRITTEN
   !   890531  Changed all specific intrinsics to generic.  (WRB)
   !   890831  Modified array declarations.  (WRB)
   !   890831  REVISION DATE from Version 3.2
   !   891214  Prologue converted to Version 4.0 format.  (BAB)
   !   900326  Removed duplicate information from DESCRIPTION section.
   !           (WRB)
   !   920501  Reformatted the REFERENCES section.  (WRB)
   !***END PROLOGUE  DGBSL
   !
   !======================================================================
   !
   !      SUBROUTINE DGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB,
   !     $                   INFO )
   !
   !  -- LAPACK routine (version 3.1) --
   !     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
   !     November 2006
   !
   !     .. Scalar Arguments ..
   !      CHARACTER          TRANS
   !      INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
   !     ..
   !     .. Array Arguments ..
   !      INTEGER            IPIV( * )
   !      DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )
   !     ..
   !
   !  Purpose
   !  =======
   !
   !  DGBTRS solves a system of linear equations
   !     A * X = B  or  A' * X = B
   !  with a general band matrix A using the LU factorization computed
   !  by DGBTRF.
   !
   !  Arguments
   !  =========
   !
   !  TRANS   (input) CHARACTER*1
   !          Specifies the form of the system of equations.
   !          = 'N':  A * X = B  (No transpose)
   !          = 'T':  A'* X = B  (Transpose)
   !          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
   !
   !  N       (input) INTEGER
   !          The order of the matrix A.  N >= 0.
   !
   !  KL      (input) INTEGER
   !          The number of subdiagonals within the band of A.  KL >= 0.
   !
   !  KU      (input) INTEGER
   !          The number of superdiagonals within the band of A.  KU >= 0.
   !
   !  NRHS    (input) INTEGER
   !          The number of right hand sides, i.e., the number of columns
   !          of the matrix B.  NRHS >= 0.
   !
   !  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
   !          Details of the LU factorization of the band matrix A, as
   !          computed by DGBTRF.  U is stored as an upper triangular band
   !          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
   !          the multipliers used during the factorization are stored in
   !          rows KL+KU+2 to 2*KL+KU+1.
   !
   !  LDAB    (input) INTEGER
   !          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
   !
   !  IPIV    (input) INTEGER array, dimension (N)
   !          The pivot indices; for 1 <= i <= N, row i of the matrix was
   !          interchanged with row IPIV(i).
   !
   !  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
   !          On entry, the right hand side matrix B.
   !          On exit, the solution matrix X.
   !
   !  LDB     (input) INTEGER
   !          The leading dimension of the array B.  LDB >= max(1,N).
   !
   !  INFO    (output) INTEGER
   !          = 0:  successful exit
   !          < 0: if INFO = -i, the i-th argument had an illegal value
   !
   !  =====================================================================
      character trans
      integer nrhs, ldb, info
   
      trans = 'N'
      if (job /= 0) trans = 'T'
      LDB = N
      NRHS = 1   
      call DGBTRS( TRANS, N, ML, MU, NRHS, ABD, LDA, IPVT, B, LDB, INFO )
      if (info < 0) then
         write(6,*) 'Error in DGBTRS call via DGBSL, INFO = ', info
         stop
      end if
      return
   end subroutine dgbsl   

   ! This uses LAPACK DGETRF to replace LINPACK DGEFA
   SUBROUTINE DGEFA (A, LDA, N, IPVT, INFO)
      implicit none
      INTEGER LDA,N,IPVT(N),INFO
      real(8) A(LDA,N)
      
   !***BEGIN PROLOGUE  DGEFA
   !***PURPOSE  Factor a matrix using Gaussian elimination.
   !***LIBRARY   SLATEC (LINPACK)
   !***CATEGORY  D2A1
   !***TYPE      DOUBLE PRECISION (SGEFA-S, DGEFA-D, CGEFA-C)
   !***KEYWORDS  GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
   !             MATRIX FACTORIZATION
   !***AUTHOR  Moler, C. B., (U. of New Mexico)
   !***DESCRIPTION
   !
   !     DGEFA factors a double precision matrix by Gaussian elimination.
   !
   !     DGEFA is usually called by DGECO, but it can be called
   !     directly with a saving in time if  RCOND  is not needed.
   !     (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) .
   !
   !     On Entry
   !
   !        A       DOUBLE PRECISION(LDA, N)
   !                the matrix to be factored.
   !
   !        LDA     INTEGER
   !                the leading dimension of the array  A .
   !
   !        N       INTEGER
   !                the order of the matrix  A .
   !
   !     On Return
   !
   !        A       an upper triangular matrix and the multipliers
   !                which were used to obtain it.
   !                The factorization can be written  A = L*U  where
   !                L  is a product of permutation and unit lower
   !                triangular matrices and  U  is upper triangular.
   !
   !        IPVT    INTEGER(N)
   !                an integer vector of pivot indices.
   !
   !        INFO    INTEGER
   !                = 0  normal value.
   !                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
   !                     condition for this subroutine, but it does
   !                     indicate that DGESL or DGEDI will divide by zero
   !                     if called.  Use  RCOND  in DGECO for a reliable
   !                     indication of singularity. 
   !
   !  Replace the above call with
   !      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
   !
   !  -- LAPACK routine (version 3.1) --
   !     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
   !     November 2006
   !
   !     .. Scalar Arguments ..
   !      INTEGER            INFO, LDA, M, N
   !     ..
   !     .. Array Arguments ..
   !      INTEGER            IPIV( * )
   !      DOUBLE PRECISION   A( LDA, * )
   !     ..
   !
   !  Purpose
   !  =======
   !
   !  DGETRF computes an LU factorization of a general M-by-N matrix A
   !  using partial pivoting with row interchanges.
   !
   !  The factorization has the form
   !     A = P * L * U
   !  where P is a permutation matrix, L is lower triangular with unit
   !  diagonal elements (lower trapezoidal if m > n), and U is upper
   !  triangular (upper trapezoidal if m < n).
   !
   !  This is the right-looking Level 3 BLAS version of the algorithm.
   !
   !  Arguments
   !  =========
   !
   !  M       (input) INTEGER
   !          The number of rows of the matrix A.  M >= 0.
   !
   !  N       (input) INTEGER
   !          The number of columns of the matrix A.  N >= 0.
   !
   !  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   !          On entry, the M-by-N matrix to be factored.
   !          On exit, the factors L and U from the factorization
   !          A = P*L*U; the unit diagonal elements of L are not stored.
   !
   !  LDA     (input) INTEGER
   !          The leading dimension of the array A.  LDA >= max(1,M).
   !
   !  IPIV    (output) INTEGER array, dimension (min(M,N))
   !          The pivot indices; for 1 <= i <= min(M,N), row i of the
   !          matrix was interchanged with row IPIV(i).
   !
   !  INFO    (output) INTEGER
   !          = 0:  successful exit
   !          < 0:  if INFO = -i, the i-th argument had an illegal value
   !          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
   !                has been completed, but the factor U is exactly
   !                singular, and division by zero will occur if it is used
   !                to solve a system of equations.
   !
   !  =====================================================================


      ! The problem faced here is that the LAPACK routine requires the
      ! column dimension of A; i.e., the M in A(M,N).  Since the LINPACK
      ! routine does not require this, just set M equal to N, since I 
      ! believe the linpack routines assumes it's square
      integer M
      M = N
      call DGETRF( M, N, A, LDA, IPVT, INFO )
      if (info < 0) then
         write(6,*) 'Error in DGETRF call via DGEFA, INFO = ', info
         stop
      end if
      return
   end subroutine DGEFA

   ! This uses LAPACK DPOTRF to replace LINPACK DPOFA
   subroutine dpofa(a,lda,n,info)
      implicit none
      integer lda,n,info
      double precision a(lda,*)
      !
      !     dpofa factors a double precision symmetric positive definite
      !     matrix.
      !
      !     dpofa is usually called by dpoco, but it can be called
      !     directly with a saving in time if  rcond  is not needed.
      !     (time for dpoco) = (1 + 18/n)*(time for dpofa) .
      !
      !     on entry
      !
      !        a       double precision(lda, n)
      !                the symmetric matrix to be factored.  only the
      !                diagonal and upper triangle are used.
      !
      !        lda     integer
      !                the leading dimension of the array  a .
      !
      !        n       integer
      !                the order of the matrix  a .
      !
      !     on return
      !
      !        a       an upper triangular matrix  r  so that  a = trans(r)*r
      !                where  trans(r)  is the transpose.
      !                the strict lower triangle is unaltered.
      !                if  info .ne. 0 , the factorization is not complete.
      !
      !        info    integer
      !                = 0  for normal return.
      !                = k  signals an error condition.  the leading minor
      !                     of order  k  is not positive definite.
      !
      !     linpack.  this version dated 08/14/78 .
      !     cleve moler, university of new mexico, argonne national lab.
      !
      !     subroutines and functions
      !
      !     blas ddot
      !     fortran sqrt
      !
      !     internal variables
      !     none
      call dpotrf ('U',n,a,lda,info) 
      !           
      return
   end subroutine dpofa     
   !====================== The end of dpofa ===============================

   ! This uses LAPACK DTRTRS to replace LINPACK DTRSL
   subroutine dtrsl(t,ldt,n,b,job,info)
      implicit none
      integer ldt,n,job,info
      double precision t(ldt,*),b(*)
      !
      !
      !     dtrsl solves systems of the form
      !
      !                   t * x = b
      !     or
      !                   trans(t) * x = b
      !
      !     where t is a triangular matrix of order n. here trans(t)
      !     denotes the transpose of the matrix t.
      !
      !     on entry
      !
      !         t         double precision(ldt,n)
      !                   t contains the matrix of the system. the zero
      !                   elements of the matrix are not referenced, and
      !                   the corresponding elements of the array can be
      !                   used to store other information.
      !
      !         ldt       integer
      !                   ldt is the leading dimension of the array t.
      !
      !         n         integer
      !                   n is the order of the system.
      !
      !         b         double precision(n).
      !                   b contains the right hand side of the system.
      !
      !         job       integer
      !                   job specifies what kind of system is to be solved.
      !                   if job is
      !
      !                        00   solve t*x=b, t lower triangular,
      !                        01   solve t*x=b, t upper triangular,
      !                        10   solve trans(t)*x=b, t lower triangular,
      !                        11   solve trans(t)*x=b, t upper triangular.
      !
      !     on return
      !
      !         b         b contains the solution, if info .eq. 0.
      !                   otherwise b is unaltered.
      !
      !         info      integer
      !                   info contains zero if the system is nonsingular.
      !                   otherwise info contains the index of
      !                   the first zero diagonal element of t.
      !
      !     linpack. this version dated 08/14/78 .
      !     g. w. stewart, university of maryland, argonne national lab.
      !
      !     subroutines and functions
      !
      !     blas daxpy,ddot
      !     fortran mod
      !
      !     internal variables
      !
      integer nrhs, ldb
      character uplo*1, trans*1, diag*1
      
      !     translate job for dtrtrs   
      select case (job)
      case (0)
         uplo = 'L'
         trans = 'N'
      case (1)
         uplo = 'U'
         trans = 'N'
      case (10)
         uplo = 'L'
         trans = 'T'
      case (11)
         uplo = 'U'
         trans = 'T'
      end select 
      diag = 'N' 
      nrhs = 1
      ldb = n
      call dtrtrs( uplo, trans, diag, n, nrhs, t, ldt, b, ldb, info )
      
      return
   end subroutine dtrsl
      
   !====================== The end of dtrsl ===============================


       