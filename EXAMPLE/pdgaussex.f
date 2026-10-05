      PROGRAM PDGAUSSEX
      IMPLICIT NONE
*
*  -- ScaLAPACK example code --
*
*     This program solves the dense linear system  A * x = b  with
*     Gaussian elimination with partial pivoting (PDGETRF + PDGETRS,
*     i.e. the two halves of PDGESV), where
*
*        A(i,j) = exp( -|| p_i - p_j ||**2 / ( 2 * SIGMA**2 ) )
*                 + LAMBDA * delta_ij,
*
*     p_1, ..., p_N are N points drawn uniformly at random in the unit
*     cube [0,1]^3 and b is the vector of all ones.
*
*     Default problem: N = 80000, distributed block-cyclically on a
*     4 x 4 process grid with 20000 x 20000 blocks, i.e. every process
*     owns exactly one 20000 x 20000 block of A (3.2 GB of local
*     storage per process, 51.2 GB in total).
*
*     Matrix assembly is fully distributed: every process generates
*     the same set of points from a fixed seed (3*N doubles, cheap)
*     and evaluates the kernel only for the entries it owns. No
*     global matrix is ever formed on a single process.
*
*     After the solve the matrix is re-assembled (instead of keeping
*     a 3.2 GB copy per process) to compute the scaled residual
*        || A*x - b ||_inf / ( ||A||_inf * ||x||_inf * N * eps ).
*     An estimate of the reciprocal condition number is also printed.
*
*     NOTE: Gaussian kernel matrices become extremely ill-conditioned
*     when SIGMA is large compared to the typical point spacing
*     (about N**(-1/3) = 0.023 for N = 80000). If RCOND is close to
*     machine precision, decrease SIGMA or add a regularization
*     LAMBDA > 0 to the diagonal.
*
*     N and the (square) block size NB may be given on the command
*     line:  xdgaussex [N [NB]]   (defaults: N = 80000, NB = 20000).
*
*     Run with at least NPROW*NPCOL (= 16) MPI processes, e.g.
*        mpirun -np 16 ./xdgaussex
*
*     .. Parameters ..
      INTEGER            NPROW, NPCOL
      PARAMETER          ( NPROW = 4, NPCOL = 4 )
      DOUBLE PRECISION   SIGMA, LAMBDA
      PARAMETER          ( SIGMA = 0.01D+0, LAMBDA = 0.0D+0 )
      INTEGER            DLEN_
      PARAMETER          ( DLEN_ = 9 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IAM, ICTXT, IL, INFO, JL, LIWORK, LLDA,
     $                   LLDB, LWORK, MYCOL, MYROW, NP, NPROCS, NQ,
     $                   NPROW0, NPCOL0, N, NB, ISTAT
      CHARACTER*32       ARG
      DOUBLE PRECISION   ANORM, EPS, RCOND, RESID, RNORM, T0, T1, T2,
     $                   T3, T4, XNORM
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA( DLEN_ ), DESCB( DLEN_ ), ISEED( 4 ),
     $                   IWQRY( 1 )
      DOUBLE PRECISION   WQRY( 1 )
      INTEGER, ALLOCATABLE :: IPIV( : ), IWORK( : ), IROW( : ),
     $                        ICOL( : )
      DOUBLE PRECISION, ALLOCATABLE :: A( :, : ), B( :, : ),
     $                                 R( :, : ), PTS( :, : ),
     $                                 WORK( : )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_ABORT, BLACS_BARRIER, BLACS_EXIT,
     $                   BLACS_GET, BLACS_GRIDEXIT, BLACS_GRIDINFO,
     $                   BLACS_GRIDINIT, BLACS_PINFO, DESCINIT, DLARNV,
     $                   PDGECON, PDGEMV, PDGETRF, PDGETRS
*     ..
*     .. External Functions ..
      INTEGER            INDXL2G, NUMROC
      DOUBLE PRECISION   DWALLTIME00, PDLAMCH, PDLANGE
      EXTERNAL           DWALLTIME00, INDXL2G, NUMROC, PDLAMCH, PDLANGE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          COMMAND_ARGUMENT_COUNT, DBLE,
     $                   GET_COMMAND_ARGUMENT, INT, MAX
*     ..
*     .. Executable Statements ..
*
      N = 80000
      NB = 20000
      ISTAT = 0
      IF( COMMAND_ARGUMENT_COUNT().GE.1 ) THEN
         CALL GET_COMMAND_ARGUMENT( 1, ARG )
         READ( ARG, *, IOSTAT = ISTAT ) N
      END IF
      IF( ISTAT.EQ.0 .AND. COMMAND_ARGUMENT_COUNT().GE.2 ) THEN
         CALL GET_COMMAND_ARGUMENT( 2, ARG )
         READ( ARG, *, IOSTAT = ISTAT ) NB
      END IF
      IF( ISTAT.NE.0 .OR. N.LT.1 .OR. NB.LT.1 ) THEN
         WRITE( *, FMT = 9988 )
         STOP
      END IF
*
      CALL BLACS_PINFO( IAM, NPROCS )
      IF( NPROCS.LT.NPROW*NPCOL ) THEN
         IF( IAM.EQ.0 )
     $      WRITE( *, FMT = 9999 ) NPROW*NPCOL, NPROCS
         CALL BLACS_EXIT( 0 )
         STOP
      END IF
*
*     Define the NPROW x NPCOL process grid
*
      CALL BLACS_GET( -1, 0, ICTXT )
      CALL BLACS_GRIDINIT( ICTXT, 'Row-major', NPROW, NPCOL )
      CALL BLACS_GRIDINFO( ICTXT, NPROW0, NPCOL0, MYROW, MYCOL )
*
*     Processes not in the grid skip to the end
*
      IF( MYROW.LT.0 .OR. MYROW.GE.NPROW .OR. MYCOL.LT.0 .OR.
     $    MYCOL.GE.NPCOL ) GO TO 10
*
      NP = NUMROC( N, NB, MYROW, 0, NPROW )
      NQ = NUMROC( N, NB, MYCOL, 0, NPCOL )
      LLDA = MAX( 1, NP )
      LLDB = MAX( 1, NP )
*
      CALL DESCINIT( DESCA, N, N, NB, NB, 0, 0, ICTXT, LLDA, INFO )
      CALL DESCINIT( DESCB, N, 1, NB, NB, 0, 0, ICTXT, LLDB, INFO )
*
      IF( IAM.EQ.0 ) THEN
         WRITE( *, FMT = 9998 ) N, NB, NB, NPROW, NPCOL, SIGMA, LAMBDA
         WRITE( *, FMT = 9997 ) NP, NQ,
     $      8.0D+0*DBLE( NP )*DBLE( NQ ) / 1.0D+9
      END IF
*
*     Allocate local storage
*
      ALLOCATE( A( LLDA, MAX( 1, NQ ) ), B( LLDB, 1 ), R( LLDB, 1 ),
     $          IPIV( NP+NB ), IROW( MAX( 1, NP ) ),
     $          ICOL( MAX( 1, NQ ) ), PTS( 3, N ), STAT = INFO )
      IF( INFO.NE.0 ) THEN
         WRITE( *, FMT = 9996 ) MYROW, MYCOL
         CALL BLACS_ABORT( ICTXT, 1 )
      END IF
*
*     Every process generates the same N random points in [0,1]^3
*     (DLARNV with identical seed, uniform (0,1) distribution)
*
      ISEED( 1 ) = 2026
      ISEED( 2 ) = 10
      ISEED( 3 ) = 1
      ISEED( 4 ) = 1
      CALL DLARNV( 1, ISEED, 3*N, PTS )
*
*     Global indices of the locally owned rows and columns
*
      DO 20 IL = 1, NP
         IROW( IL ) = INDXL2G( IL, NB, MYROW, 0, NPROW )
   20 CONTINUE
      DO 30 JL = 1, NQ
         ICOL( JL ) = INDXL2G( JL, NB, MYCOL, 0, NPCOL )
   30 CONTINUE
*
*     Distributed assembly of A and of the right-hand side b = ones
*
      CALL BLACS_BARRIER( ICTXT, 'All' )
      T0 = DWALLTIME00()
      CALL ASSEMBLE( NP, NQ, IROW, ICOL, PTS, SIGMA, LAMBDA, A,
     $               LLDA )
      DO 40 IL = 1, NP
         B( IL, 1 ) = ONE
   40 CONTINUE
      CALL BLACS_BARRIER( ICTXT, 'All' )
      T1 = DWALLTIME00()
*
*     ||A||_1 is needed for the condition number estimate
*     (A is symmetric so ||A||_1 = ||A||_inf)
*
      LWORK = MAX( NP, NQ ) + 2*NB
      ALLOCATE( WORK( LWORK ) )
      ANORM = PDLANGE( '1', N, N, A, 1, 1, DESCA, WORK )
      DEALLOCATE( WORK )
*
*     Gaussian elimination with partial pivoting: A = P * L * U
*
      CALL PDGETRF( N, N, A, 1, 1, DESCA, IPIV, INFO )
      CALL BLACS_BARRIER( ICTXT, 'All' )
      T2 = DWALLTIME00()
      IF( INFO.NE.0 ) THEN
         IF( IAM.EQ.0 )
     $      WRITE( *, FMT = 9995 ) INFO, INFO
         GO TO 50
      END IF
*
*     Forward and back substitution, b is overwritten by x
*
      CALL PDGETRS( 'No transpose', N, 1, A, 1, 1, DESCA, IPIV, B, 1,
     $              1, DESCB, INFO )
      CALL BLACS_BARRIER( ICTXT, 'All' )
      T3 = DWALLTIME00()
      IF( INFO.NE.0 ) THEN
         IF( IAM.EQ.0 )
     $      WRITE( *, FMT = 9989 ) 'PDGETRS', INFO
         GO TO 50
      END IF
*
*     Estimate the reciprocal condition number of A
*
      CALL PDGECON( '1', N, A, 1, 1, DESCA, ANORM, RCOND, WQRY, -1,
     $              IWQRY, -1, INFO )
      LWORK = MAX( 1, INT( WQRY( 1 ) ) )
      LIWORK = MAX( 1, IWQRY( 1 ) )
      ALLOCATE( WORK( LWORK ), IWORK( LIWORK ) )
      CALL PDGECON( '1', N, A, 1, 1, DESCA, ANORM, RCOND, WORK, LWORK,
     $              IWORK, LIWORK, INFO )
      DEALLOCATE( WORK, IWORK )
      IF( INFO.NE.0 ) THEN
         IF( IAM.EQ.0 )
     $      WRITE( *, FMT = 9989 ) 'PDGECON', INFO
         GO TO 50
      END IF
*
*     Re-assemble A and compute r = b - A*x with b = ones
*
      CALL ASSEMBLE( NP, NQ, IROW, ICOL, PTS, SIGMA, LAMBDA, A,
     $               LLDA )
      DO 60 IL = 1, NP
         R( IL, 1 ) = ONE
   60 CONTINUE
      CALL PDGEMV( 'No transpose', N, N, -ONE, A, 1, 1, DESCA, B, 1, 1,
     $             DESCB, 1, ONE, R, 1, 1, DESCB, 1 )
      CALL BLACS_BARRIER( ICTXT, 'All' )
      T4 = DWALLTIME00()
*
      LWORK = MAX( NP, NQ ) + 2*NB
      ALLOCATE( WORK( LWORK ) )
      EPS = PDLAMCH( ICTXT, 'Epsilon' )
      RNORM = PDLANGE( 'I', N, 1, R, 1, 1, DESCB, WORK )
      XNORM = PDLANGE( 'I', N, 1, B, 1, 1, DESCB, WORK )
      DEALLOCATE( WORK )
      RESID = RNORM / ( ANORM*XNORM*DBLE( N )*EPS )
*
      IF( IAM.EQ.0 ) THEN
         WRITE( *, FMT = 9994 ) T1 - T0, T2 - T1, T3 - T2, T4 - T3
         WRITE( *, FMT = 9993 ) ANORM, RCOND, XNORM, RNORM, RESID
         IF( RCOND.LT.DBLE( N )*EPS ) THEN
            WRITE( *, FMT = 9992 )
         ELSE IF( RESID.LT.10.0D+0 ) THEN
            WRITE( *, FMT = 9991 )
         ELSE
            WRITE( *, FMT = 9990 )
         END IF
      END IF
*
   50 CONTINUE
      DEALLOCATE( A, B, R, IPIV, IROW, ICOL, PTS )
      CALL BLACS_GRIDEXIT( ICTXT )
*
   10 CONTINUE
      CALL BLACS_EXIT( 0 )
*
 9999 FORMAT( 'This example needs at least ', I6, ' processes, ',
     $        'only ', I6, ' available.' )
 9998 FORMAT( /'Gaussian kernel system: N = ', I8, ', block = ', I6,
     $        ' x ', I6, ', grid = ', I3, ' x ', I3, /
     $        'SIGMA = ', ES10.3, ', LAMBDA = ', ES10.3 )
 9997 FORMAT( 'Local block of A on process (0,0): ', I8, ' x ', I8,
     $        ' (', F8.2, ' GB)' )
 9996 FORMAT( 'Process (', I3, ',', I3, '): memory allocation failed' )
 9995 FORMAT( 'PDGETRF failed: U(', I8, ',', I8, ') is exactly zero' )
 9994 FORMAT( /'Time assembly          (s) = ', F12.3, /
     $        'Time LU (PDGETRF)      (s) = ', F12.3, /
     $        'Time solve (PDGETRS)   (s) = ', F12.3, /
     $        'Time cond + residual   (s) = ', F12.3 )
 9993 FORMAT( /'||A||_1                          = ', ES12.5, /
     $        'RCOND (estimate)                 = ', ES12.5, /
     $        '||x||_inf                        = ', ES12.5, /
     $        '||b - A*x||_inf                  = ', ES12.5, /
     $        '||b-Ax||/(||A|| ||x|| N eps)     = ', ES12.5 )
 9992 FORMAT( /'WARNING: A is numerically singular, the solution is ',
     $        'not meaningful.', /'Decrease SIGMA or increase LAMBDA.' )
 9991 FORMAT( /'The answer is correct (backward stable solve).' )
 9990 FORMAT( /'The answer is suspicious.' )
 9988 FORMAT( 'Usage: xdgaussex [N [NB]]  (N, NB positive integers)' )
 9989 FORMAT( A, ' failed with INFO = ', I8 )
*
*     End of PDGAUSSEX
*
      END
*
      SUBROUTINE ASSEMBLE( NP, NQ, IROW, ICOL, PTS, SIGMA, LAMBDA,
     $                     A, LDA )
*
*     Fill the local part A(1:NP,1:NQ) of the distributed Gaussian
*     kernel matrix. IROW/ICOL hold the global row/column indices of
*     the local rows/columns, PTS(1:3,*) are the points.
*
*     .. Scalar Arguments ..
      INTEGER            LDA, NP, NQ
      DOUBLE PRECISION   LAMBDA, SIGMA
*     ..
*     .. Array Arguments ..
      INTEGER            ICOL( * ), IROW( * )
      DOUBLE PRECISION   A( LDA, * ), PTS( 3, * )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IL, J, JL
      DOUBLE PRECISION   DX, DY, DZ, S, XJ, YJ, ZJ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          EXP
*     ..
*     .. Executable Statements ..
*
      S = -1.0D+0 / ( 2.0D+0*SIGMA*SIGMA )
      DO 20 JL = 1, NQ
         J = ICOL( JL )
         XJ = PTS( 1, J )
         YJ = PTS( 2, J )
         ZJ = PTS( 3, J )
         DO 10 IL = 1, NP
            I = IROW( IL )
            DX = PTS( 1, I ) - XJ
            DY = PTS( 2, I ) - YJ
            DZ = PTS( 3, I ) - ZJ
            A( IL, JL ) = EXP( S*( DX*DX + DY*DY + DZ*DZ ) )
            IF( I.EQ.J )
     $         A( IL, JL ) = A( IL, JL ) + LAMBDA
   10    CONTINUE
   20 CONTINUE
      RETURN
*
*     End of ASSEMBLE
*
      END
