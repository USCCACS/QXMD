



SUBROUTINE INVERS( W, WI, N, WORK, IP, MX, IERR )
!---------- inverse matrix of W ----------------------------------------
!   input   W   : matrix
!           N   : order of the matrix
!   output  WI  : inverse matrix
!
!     WORK, IP  : work area
!-----------------------------------------------------------------------
IMPLICIT REAL*8 (A-H,O-Z)
DIMENSION  W(MX,MX)
DIMENSION  WI(MX,MX), WORK(MX,MX), IP(MX)
DATA DZERO / 1.0D-30 /

IERR = 0
elmax = 0.d0
DO I = 1, N
DO J = 1, N
   elmax = max( elmax, abs(W(I,J)) )
end do
end do
!--- DET : determinant of W --------------------------------------------
DO 10 I = 1, N
DO 10 J = 1, N
   WORK(I,J) = W(I,J)/elmax
10 CONTINUE
CALL DETMAT( DET, WORK, N, IP, MX, IERR )
IF( IERR.NE.0 ) RETURN
!      WRITE(*,*) DET
IF( ABS(DET).LT.DZERO ) THEN
    IERR = 1
    RETURN
ENDIF

do I = 1, N
do J = 1, N
   ISUM = 0
   do II = 1, N
      IF( II.EQ.J ) cycle
      ISUM = ISUM + 1
      JSUM = 0
      do JJ = 1, N
         IF( JJ.EQ.I ) cycle
         JSUM = JSUM + 1
         WORK(ISUM,JSUM) = W(II,JJ)/elmax
      end do
   end do
   CALL DETMAT( DETIJ, WORK, N-1, IP, MX, IERR )
   IF( IERR.NE.0 ) RETURN
   WI(I,J) = ( (-1.D0)**(I+J) )*DETIJ/DET
end do
end do
DO I = 1, N
DO J = 1, N
   WI(I,J) = WI(I,J)/elmax
end do
end do

RETURN
END


SUBROUTINE DETMAT( DET, W, N, IP, MX, IERR )
!--------- DET : determinant of W --------------------------------------
IMPLICIT REAL*8 (A-H,O-Z)
DIMENSION  W(MX,MX), IP(MX)
DATA DZERO / 1.0D-30 /

IERR = 0

DO 10 I = 1, N
   IP(I) = I
10 CONTINUE
FUGO = 1.D0
DET  = 1.D0

DO 20 I = 1, N
   PIVOT = 0.D0
   DO 21 J = 1, N
      PIVOT = MAX( PIVOT, ABS(W(I,J)) )
21    CONTINUE
   IF( PIVOT.LT.DZERO ) THEN
!             WRITE(*,*) ' *** error-1 in subroutine DETMAT'
!             IERR = 2
       DET = 0.d0
       RETURN
   ENDIF
   DO 22 J = 1, N
      W(I,J) = W(I,J)/PIVOT
22    CONTINUE
   DET = DET * PIVOT
20 CONTINUE
!--- reducing the matrix to triangular form ---
DO 1000 I = 1, N - 1
   PIVOT = 0.0
   DO 1100 J = I, N
      JJ = IP(J)
      IF( ABS(W(JJ,I)).GT.PIVOT ) THEN
          PIVOT = ABS(W(JJ,I))
          JPV   = J
      ENDIF
1100    CONTINUE
   IF( PIVOT.LT.DZERO ) THEN
!             WRITE(*,*) ' *** error-2 in subroutine DETMAT'
!             IERR = 3
       DET = 0.d0
       RETURN
   ENDIF
   IF( JPV.NE.I ) THEN
       IPOL    = IP(JPV)
       IP(JPV) = IP(I)
       IP(I)   = IPOL
       FUGO    = FUGO * ( -1.D0 )
   ENDIF
   IPV = IP(I)
   DO 1200 J = I + 1, N
      JJ  = IP(J)
      DMJ = W(JJ,I)/W(IPV,I)
      W(JJ,I) = 0.0
      DO 1200 J2 = I + 1, N
         W(JJ,J2) = W(JJ,J2) - DMJ*W(IPV,J2)
1200    CONTINUE
1000 CONTINUE
DO 1500 I = 1, N
   IV = IP(I)
   DET = DET * W(IV,I)
1500 CONTINUE
DET = DET * FUGO
!-----------------------------------------------------------------------
RETURN
END




SUBROUTINE SPLINS( L, N, X, Y, Q, B, A, W, IP, KQM, MM, MM1, NDV )
!-----------------------------------------------------------------------
!     determination of coefficients of spline function
!-----------------------------------------------------------------------
!     L    ..... order of B-spline function
!     N    ..... No. of data point ( X(I), Y(I) ) ( I = 1...N )
!     Q    ..... setten
!     B    ..... matrix whose components are values of B-spline function
!     A    ..... coefficients of spline function
!     W    ..... work area
!     IP   ..... work area for Gauss-Jordan elimination method
!-----------------------------------------------------------------------
IMPLICIT REAL*8 ( A-H, O-Z )
DIMENSION   X(MM), Y(MM)
DIMENSION   Q(KQM), B(MM,MM1), A(MM,0:NDV), W(MM), IP(MM)

K = L + 1

!-----------------------------------------------------------------------
!     determination of setten : Q
DO 10 I = 1, K
   Q(I) = X(1)
10 CONTINUE
DO 20 I = 1, N - K
   Q(I+K) = ( X(I) + X(I+K) )*0.5D0
20 CONTINUE
DO 30 I = 1, K
   Q(N+I) = X(N)
30 CONTINUE
!-----------------------------------------------------------------------
!     B-spline function : B     by   de Boor-Cox method
DO 90 I = 1, N
DO 90 J = 1, N
   B(J,I) = 0.0
90 CONTINUE
DO 100 J = 1, N
   XX = X(J)
   CALL DEBOOR( N, L, XX, Q, W, IB, KQM, MM )
   DO 140 I = IB - L, IB
      B(J,I) = W(I)
140    CONTINUE
100 CONTINUE
DO 150 I = 1, N
   B(I,N+1) = Y(I)
150 CONTINUE


!-----------------------------------------------------------------------
!     determination of coefficients : A
!       A(I,NDV) is the coefficients for NDV-th derivatives.
CALL GSSJOR( B, IP, N, MM, MM1, INDER )
DO 200 I = 1, N
   A(I,0) = B(I,N+1)
200 CONTINUE
DO 210 ND = 1, NDV
   DO 220 I = 1, ND
      A(I,ND) = 0.0
220    CONTINUE
   DO 230 I = ND + 1, N
      A(I,ND) = DBLE(K-ND)*( A(I,ND-1) - A(I-1,ND-1) )  &
&                         /( Q(I+K-ND) - Q(I) )
230    CONTINUE
210 CONTINUE


RETURN
END


SUBROUTINE SPLINV( L, N, XX, YY, Q, A, W, KQM, MM, NDV, LD )
!-----------------------------------------------------------------------
!     interpolation by spline function
!-----------------------------------------------------------------------
!     L    ..... order of B-spline function
!     N    ..... No. of data point ( X(I), Y(I) ) ( I = 1...N )
!     Q    ..... setten
!     A    ..... coefficients of spline function
!     W    ..... work area
!-----------------------------------------------------------------------
IMPLICIT REAL*8 ( A-H, O-Z )
DIMENSION   Q(KQM), A(MM,0:NDV), W(MM)

CALL DEBOOR( N, L-LD, XX, Q, W, IB, KQM, MM )

YY = 0.0
DO 140 I = IB - (L-LD), IB
   YY = YY + A(I,LD) * W(I)
140 CONTINUE


RETURN
END


SUBROUTINE DEBOOR( N, L, XX, Q, W, IB, KQM, MM )
!-----------------------------------------------------------------------
!     B-spline function : B     by   de Boor-Cox method
!-----------------------------------------------------------------------
IMPLICIT REAL*8 ( A-H, O-Z )
DIMENSION   Q(KQM), W(MM)
DATA  DZERO / 1.0D-10 /

K = L + 1

do I = 1, N
   W(I) = 0.0
end do
IB   = N
do I = K, N
   IF( Q(I).LE.XX .AND. Q(I+1).GT.XX ) THEN
       IB   = I
       exit
   END IF
end do
W(IB) = 1.D0
do NO = 1, L
   KK = NO + 1
   do I = IB - NO, IB
      A1 = 0.0
      A2 = 0.0
      IF( ABS( Q(I+KK-1) - Q(I) ).GT.DZERO )  &
&         A1 = ( XX - Q(I) )/( Q(I+KK-1) - Q(I) ) * W(I)
      IF( ABS( Q(I+KK) - Q(I+1) ).GT.DZERO )  &
&         A2 = ( Q(I+KK) -XX )/( Q(I+KK) - Q(I+1) ) * W(I+1)
      W(I) = A1 + A2
   end do

end do


RETURN
END




            SUBROUTINE GSSJOR( W, IP, N, MM, MM1, INDER )
!-----------------------------------------------------------------------
!        <<<   Gauss-Jurdan elemination method >>>   since  '90/09/07
!-----------------------------------------------------------------------
!  ( input )
!     N        ......  number of variables
!     W(i,j)   ......  coefficients
!     IP       ......  work area
!  ( output )
!     W(i,N1)  ......  solution x_i
!-----------------------------------------------------------------------
IMPLICIT REAL*8 ( A-H, O-Z )
DIMENSION  W(MM,MM1),  IP(MM)
DATA DZERO / 1.0D-30 /

N1 = N + 1

DO 10 I = 1, N
   IP(I) = I
10 CONTINUE

DO 20 I = 1, N
   PIVOT = 0.0
   DO 21 J = 1, N1
      PIVOT = MAX( PIVOT, ABS(W(I,J)) )
21    CONTINUE
   IF( PIVOT.LT.DZERO ) THEN
       INDER = 5
       RETURN
   ENDIF
   DO 22 J = 1, N1
      W(I,J) = W(I,J)/PIVOT
22    CONTINUE
20 CONTINUE


!--- reducing the matrix to triangular form ----------------------------
DO 1000 I = 1, N - 1

   PIVOT = 0.0
   DO 1100 J = I, N
      JJ = IP(J)
      IF( ABS(W(JJ,I)).GT.PIVOT ) THEN
          PIVOT = ABS(W(JJ,I))
          JPV   = J
      ENDIF
1100    CONTINUE
   IF( PIVOT.LT.DZERO ) THEN
       INDER = 6
       RETURN
   ENDIF
   IF( JPV.NE.I ) THEN
       IPOL    = IP(JPV)
       IP(JPV) = IP(I)
       IP(I)   = IPOL
   ENDIF

   IPV = IP(I)
   DO 1200 J = I + 1, N
      JJ  = IP(J)
      DMJ = W(JJ,I)/W(IPV,I)
!            W(JJ,I) = DMJ
      W(JJ,I) = 0.0
      DO 1200 J2 = I + 1, N1
         W(JJ,J2) = W(JJ,J2) - DMJ*W(IPV,J2)
1200    CONTINUE

1000 CONTINUE



!--- solutions ---------------------------------------------------------
DO 1500 I = 1, N
   IV = IP(I)
   DIEL = W(IV,I)
   W(IV,I) = 1.0D0
   DO 1500 J = I + 1, N1
      W(IV,J) = W(IV,J)/DIEL
1500 CONTINUE

IPN = IP(N)
DO 1550 I = 1, N
   IV = IP(I)
   W(IPN,I) = W(IV,N1)
1550 CONTINUE
DO 1600 I = N - 1, 1, -1
   IV = IP(I)
   DO 1600 J = I + 1, N
      W(IPN,I) = W(IPN,I) - W(IV,J)*W(IPN,J)
1600 CONTINUE
DO 1650 I = 1, N
   W(I,N1) = W(IPN,I)
1650 CONTINUE

!--- value for the normal return
 INDER = 0


 RETURN
 END




!#if ! PCHOLESKY
SUBROUTINE RsEIGQR( AR, EVR, N, MTRX, EPSI, IFVEC, IFSRT,  &
&                   CHECK, W1R, W1I, work )
!-----------------------------------------------------------------------
!    Eigenvalues and eigenvectors of real symmetric  matrix   1992/10/16
!
!     Householders method for reducing the matrix to tridiagonal form
!     QR method for obtaining the eigenvalues
!
!      up-date  1992/11/12
!-----------------------------------------------------------------------
!  ( input )
!     AR    ...... Real parts of the Hermitian matrix elements
!       ( AR_ij ) for i>=j must be set correctly.
!     N     ...... Order of the Hermitian matrix
!     IFVEC ......  = 1 : obtain eigenvectors
!                   = 0 : do not obtain eigenvectors
!     IFSRT ......  = 2 : sorting eigenvalues and eigenvectors
!                   = 1 : sorting eigenvalues only ( W1I: sorting list )
!                   = 0 : do not sorting
!     EPSI  ...... Accuracy
!  ( output )
!     AR    ...... Real parts of the eigenvectors
!      ( AR_ij ) for i=1...N is eigenvector of eigenvalue, EVR_j.
!     EVR   ...... The eigenvalues
!  *****  W1R, W1I ...... work area  *****
!-----------------------------------------------------------------------
implicit none
integer :: N, MTRX
real*8  :: AR(MTRX,N), EVR(N)
real*8  :: EPSI
integer :: IFVEC, IFSRT
character*5 :: CHECK
real*8  :: W1R(N), W1I(N), work(N)

!------declare local variables
integer :: k, k1, i, j, KN, KN1, KAI, KNTSF, KNM, IOR, IO, l, IT
real*8  :: T1, TS, A1, U1R, DMX, WMX, ESHIFT, evol, EVRSFT
real*8  :: COLD, AA, DCOSR, DSINR, EK1R, ARK, W1NR, BK1, DSOLR, ARK1, EVRI, ARI
LOGICAL   VEC
real*8,  parameter :: ZERO = 1.0D-60
integer, parameter :: ITMIN = 2, KMAX = 50
logical :: lzero, lconvg


!      WRITE(*,*) 'ITMIN, KMAX'
!      READ(*,*) ITMIN,KMAX
VEC = IFVEC .EQ. 1


!-----------------------------------------------------------------------
!----    Householders method for reducing the matrix to tridiagonal form


DO K = 1, N-2
   K1 = K + 1

   T1 = 0.0
   DO I = K1, N
      T1 = T1 + AR(I,K)*AR(I,K)
   end do
   lzero = .true.
   DO I = K1, N
      lzero = lzero .and. abs(AR(I,K)) < ZERO
   end do
!         IF( ABS(T1) > ZERO ) THEN
   if( .not.lzero ) then

   TS = SQRT(T1) * SIGN( 1.0D0, AR(K1,K) )

   A1 = 1.0D0/( T1 + AR(K1,K)*TS )
   A1 = SQRT(A1)

   W1R(K) = 0.0
   W1R(K1) = (AR(K1,K) + TS)*A1
   W1R(K+2:N) = AR(K+2:N,K)*A1
   AR(K,K1:N) = W1R(K1:N)

!         DO 1300 I = K, N
!            W1I(I) = 0.0
!            DO 1310 J = K1, I
!               W1I(I) = W1I(I) + AR(I,J)*W1R(J)
! 1310       CONTINUE
!            DO 1320 J = I+1, N
!               W1I(I) = W1I(I) + AR(J,I)*W1R(J)
! 1320       CONTINUE
! 1300    CONTINUE
   W1I(K:N) = 0.d0
   DO J = K1, N
      W1I(J:N) = W1I(J:N) + AR(J:N,J)*W1R(J)
   end do
   DO I = K, N
      DO J = I+1, N
         W1I(I) = W1I(I) + AR(J,I)*W1R(J)
      end do
   end do

   U1R = 0.0
   DO I = K1, N
      U1R = U1R + W1I(I)*W1R(I)
   end do
   U1R = -0.5D0*U1R

   W1I(K:N) = W1I(K:N) + U1R*W1R(K:N)

   AR(K:K+1,K) = AR(K:K+1,K) - W1R(K:K+1)*W1I(K) - W1I(K:K+1)*W1R(K)
   AR(K+2:N,K) = 0.0
   DO J = K1, N
      AR(J:N,J) = AR(J:N,J) - W1R(J:N)*W1I(J) - W1I(J:N)*W1R(J)
   end do

   else

       AR(K,K1:N) = 0.0

   end if

!check
!         WRITE(2,*) 'K=',K
!         DO 9700 L1 = 1, N
!            WRITE(2,9000) (AR(L1,L2),AI(L1,L2),L2=1,N)
! 9700    CONTINUE
! 9000    FORMAT(100('(',2E12.4') '))

end do




!-----------------------------------------------------------------------
!----   Store matrix elements
!----      EVR        :  diagonal elements ( real )
!----      W1R        :  non-diagonal (i,i+1)-elements


do I = 1, N
   EVR(I) = AR(I,I)
end do
do I = 1, N-1
   W1R(I) = AR(I+1,I)
end do


IF( VEC ) then
!-----------------------------------------------------------------------
!----  Set the Unitary matrix elements


do I = 1, N
   AR(I,I) = 1.0D0
end do
DO J = 1, N - 1
   AR(J+1:N,J) = 0.d0
end do


DO K = 1, N-2
   KN  = N - K
   KN1 = KN - 1

   AR(KN,KN+1:N) = 0.0
   w1i(KN:N) = AR(KN1,KN:N)
   DO I = KN, N
      work(i) = 0.d0
      DO J = KN, N
         work(i) = work(i) + w1i(j)*AR(J,I)
      end do
   end do
   DO I = KN, N
      AR(KN:N,I) = AR(KN:N,I) - w1i(KN:N)*work(i)
   end do
end do

AR(1,2:N) = 0.0


!check
!         WRITE(2,*) 'U'
!         DO 710 L1 = 1, N
!            WRITE(2,9000) (AR(L1,L2),L2=1,N)
!  710    CONTINUE
end if




!-----------------------------------------------------------------------
!----               Givens' QR method


DMX = 0.0
WMX = 0.0
DO I = 1, N
   DMX = MAX( DMX, ABS(EVR(I)) )
end do
DO I = 1, N-1
   WMX = MAX( WMX, ABS(W1R(I)) )
end do

ESHIFT = MAX( WMX, DMX ) * 1.1D0
EVR(1:N) = EVR(1:N) + ESHIFT

!check
!         WRITE(*,*) 'A_initial'
!         WRITE(*,9002) (EVR(L1),L1=1,N)
!         WRITE(*,9002) (W1R(L1),L1=1,N-1)
!         WRITE(*,9002) (W1I(L1),L1=1,N-1)
! 9002    FORMAT(6E12.4)


KAI   = 0
KNTSF = 0
KNM   = N
evol  = 1.d+10

EVRSFT = 0.0
iterdo: do
   !----   eigenvalue shift
   EVR(1:KNM) = EVR(1:KNM) - EVRSFT


KAI   = KAI   + 1
KNTSF = KNTSF + 1

COLD = 1.0D0
DO K = 2, KNM
   K1 = K - 1
   AA   = EVR(K1)*EVR(K1) + W1R(K1)*W1R(K1)
   AA   = SQRT(AA)
!           IF( ABS(AA).LT.ZERO ) THEN
!               WRITE(*,*) AA
!               CHECK = 'ZERO '
!               RETURN
!           ENDIF
   DCOSR = EVR(K1)/AA
   DSINR = W1R(K1)/AA

   EK1R = W1R(K1)*COLD

   ARK  = EVR(K)
   W1NR = EK1R*DCOSR + ARK*DSINR
   EVR(K) = -EK1R*DSINR + ARK*DCOSR

   BK1 = COLD*DCOSR

   EVR(K1) = AA*BK1 + W1NR*DSINR
   IF( K.GE.3 ) THEN
       W1R(K-2) = AA*DSOLR
   ENDIF

   COLD  = DCOSR
   DSOLR = DSINR
   IF( VEC ) THEN
       DO I = 1, N
          ARK1 = AR(I,K-1)
          ARK  = AR(I,K)
          AR(I,K-1) =  ARK1*DCOSR + ARK*DSINR
          AR(I,K)   = -ARK1*DSINR + ARK*DCOSR
       end do
   ENDIF
end do
ARK1 = EVR(KNM)
EVR(KNM)   = ARK1*COLD
W1R(KNM-1) = ARK1*DSOLR



!----   diagonal part

EVR(1:KNM) = EVR(1:KNM) + EVRSFT

!      WRITE(*,*) 'kaisuu',KAI
!      WRITE(*,9002) (EVR(L1),L1=1,N)
!      WRITE(*,9002) (EVR(L1)-ESHIFT,L1=1,N)
!C         WRITE(*,9002) (W1R(L1),L1=1,N-1)
!C         WRITE(*,9002) (W1I(L1),L1=1,N-1)


!----   hantei

IF( KNTSF.GT.KMAX ) THEN
    CHECK = 'KAISU'
    RETURN
ENDIF
DMX = 0.0
DO I = 1, KNM
   DMX = MAX( DMX, ABS(EVR(I)) )
end do
DMX = DMX*EPSI

IF( ABS(W1R(KNM-1)).LT.DMX .or.  &
&   kntsf > 2 .and. abs(evol-EVR(KNM)) < zero ) then
    KNM = KNM - 1
    KNTSF = 0
ENDIF
evol = EVR(KNM)

lconvg = .false.
do I = 1, KNM - 1
   lconvg = lconvg .or. ABS(W1R(I)) > DMX
end do
if( .not.lconvg ) exit
   IF( KAI > ITMIN ) then
!       if( abs(EVR(KNM)-EVR(KNM-1)) > zero ) then
           EVRSFT = EVR(KNM)
!       else
!           !---2012/09/04
!           EVRSFT = EVR(KNM) * (1.d0 + 1.d-10)
!       end if
       !---2014/03/09
       do I = 1, KNM - 1 
          if( abs(EVR(KNM)-EVR(i)) < zero ) then
              EVRSFT = EVR(KNM) * (1.d0 + 1.d-10)
              exit
          end if
       end do
   end if
end do iterdo

!      WRITE(*,*) '*** iteration :',KAI
!-----------------------------------------------------------------------

EVR(1:N) = EVR(1:N) - ESHIFT


IF( IFSRT.LE.0 ) RETURN
!-----------------------------------------------------------------------
!----    sorting eigenvalues


W1I(1)= 1

ido: do I = 2, N
   jdo: do J = 1, I-1
   IF( EVR(I).GE.EVR(J) ) cycle jdo

      EVRI = EVR(I)
      do L = I, J+1, -1
         EVR(L) = EVR(L-1)
         W1I(L) = W1I(L-1)
      end do
      EVR(J) = EVRI
      W1I(J) = I
      cycle ido

   end do jdo
   W1I(I) = I
end do ido



IF( .NOT.VEC .OR. IFSRT.LE.1 ) RETURN

!----    sorting eigenvectors


W1R(1:N) = 0
ido2: do I = 1, N
   IF( NINT(W1R(I)).NE.0 ) cycle ido2
   W1R(I) = 10.0
   IT = NINT(W1I(I))
   IF( I.EQ.IT ) cycle ido2
       IOR = I
       IO  = I

       do
       do L=1,N
          ARI = AR(L,IO)
          AR(L,IO)  = AR(L,IT)
          AR(L,IT) = ARI
       end do
       IO = IT
       W1R(IO) = 10.0
       IT = NINT(W1I(IO))
       IF( IT.EQ.IOR ) exit
       end do

end do ido2


RETURN
END




SUBROUTINE MsEIGQR( AR, AI, EVR, N, MTRX, EPSI, IFVEC, IFSRT,  &
& CHECK, W1R, W1I, W2R, work1, work2, bwork1, bwork2 )
!-----------------------------------------------------------------------
!     Eigenvalues and eigenvectors of Hermitian matrix        1992/10/16
!
!     Householders method for reducing the matrix to tridiagonal form
!     QR method for obtaining the eigenvalues
!
!      up-date  1992/10/20
!-----------------------------------------------------------------------
!  ( input )
!     AR    ...... Real parts of the Hermitian matrix elements
!     AI    ...... Imaginary parts of the Hermitian matrix elements
!       ( AR_ij, AI_ij ) for i>=j must be set correctly.
!     N     ...... Order of the Hermitian matrix
!     IFVEC ......  = 1 : obtain eigenvectors
!                   = 0 : do not obtain eigenvectors
!     IFSRT ......  = 2 : sorting eigenvalues and eigenvectors
!                   = 1 : sorting eigenvalues only ( W1I: sorting list )
!                   = 0 : do not sorting
!     EPSI  ...... Accuracy
!  ( output )
!     AR    ...... Real parts of the eigenvectors
!     AI    ...... Imaginary parts of the eigenvectors
!      ( AR_ij, AI_ij ) for i=1...N is eigenvector of eigenvalue, EVR_j.
!     EVR   ...... The eigenvalues
!  *****  W1R, W1I, W2R ...... work area  *****
!-----------------------------------------------------------------------
implicit none
integer :: N, MTRX
real*8  :: AR(MTRX,*), AI(MTRX,*), EVR(*)
real*8  :: EPSI
integer :: IFVEC, IFSRT
character*5 :: CHECK
real*8  :: W1R(*),W1I(*),W2R(*)
real*8  :: work1(*), work2(*), bwork1(*), bwork2(*)

!------declare local variables
integer :: k, k1, i, j, kn, kn1, KAI, KNTSF, KNM, l, it, IOR, IO
real*8  :: t1, ts, AAK1, AAKR, TSR, TSI, a1, U1R, U1I, DMX, WMX, ESHIFT
real*8  :: EVRSFT, COLD, aa, DCOSR, DSINR, DSINI, EK1R, EK1I, ARK, W1NR, W1NI, BK1
real*8  :: DSOLR, DSOLI, ARK1, AIK1, AIK, ARI, AII, EVRI, evol
LOGICAL :: VEC
real*8,  parameter :: ZERO = 1.0D-60
integer, parameter :: ITMIN = 2, KMAX = 50
logical :: lzero, lconvg


!      WRITE(*,*) 'ITMIN, KMAX'
!      READ(*,*) ITMIN,KMAX
VEC = IFVEC .EQ. 1


!-----------------------------------------------------------------------
!----    Householders method for reducing the matrix to tridiagonal form


kdo: do K = 1, N-2
   K1 = K + 1

   T1 = 0.0
   DO I = K1, N
      T1 = T1 + AR(I,K)*AR(I,K) + AI(I,K)*AI(I,K)
   end do
   lzero = .true.
   DO I = K1, N
      lzero = lzero .and. AR(I,K)*AR(I,K) + AI(I,K)*AI(I,K) < ZERO
   end do
!         IF( ABS(T1).LT.ZERO ) THEN
   if( lzero ) then
       AR(K,K1:N) = 0.0
       AI(K,K1:N) = 0.0
       cycle kdo
   ENDIF

   TS = SQRT(T1)

   AAK1 = AR(K1,K)*AR(K1,K) + AI(K1,K)*AI(K1,K)
   IF( ABS(AAK1).GT.ZERO ) THEN
       AAK1 = SQRT(AAK1)
       AAKR = TS/AAK1
       TSR  = AR(K1,K)*AAKR
       TSI  = AI(K1,K)*AAKR
       A1   = T1 + AAK1*TS
     ELSE
       TSR  = TS
       TSI  = 0.0
       A1   = T1
   ENDIF

   A1 = 1.0D0/A1
   A1 = SQRT(A1)

   W1R(K) = 0.0
   W1I(K) = 0.0
   W1R(K1) = (AR(K1,K) + TSR)*A1
   W1I(K1) = (AI(K1,K) + TSI)*A1
   W1R(K+2:N) = AR(K+2:N,K)*A1
   W1I(K+2:N) = AI(K+2:N,K)*A1
   AR(K,K1:N) = W1R(K1:N)
   AI(K,K1:N) = W1I(K1:N)

   W2R(K:N) = 0.0
   EVR(K:N) = 0.0
   DO J = K1, N
      W2R(J:N) = W2R(J:N) + AR(J:N,J)*W1R(J) - AI(J:N,J)*W1I(J)
      EVR(J:N) = EVR(J:N) + AR(J:N,J)*W1I(J) + AI(J:N,J)*W1R(J)
   end do
   DO I = K, N
      DO J = I+1, N
         W2R(I) = W2R(I) + AR(J,I)*W1R(J) + AI(J,I)*W1I(J)
         EVR(I) = EVR(I) + AR(J,I)*W1I(J) - AI(J,I)*W1R(J)
      end do
   end do

   U1R = 0.0
   U1I = 0.0
   DO I = K1, N
      U1R = U1R + W2R(I)*W1R(I) + EVR(I)*W1I(I)
      U1I = U1I + W2R(I)*W1I(I) - EVR(I)*W1R(I)
   end do
   U1R = -0.5D0*U1R
   U1I = -0.5D0*U1I

   W2R(K:N) = W2R(K:N) + U1R*W1R(K:N) - U1I*W1I(K:N)
   EVR(K:N) = EVR(K:N) + U1R*W1I(K:N) + U1I*W1R(K:N)

   AR(K:K+1,K) = AR(K:K+1,K)  &
&              - ( W1R(K:K+1)*W2R(K) + W1I(K:K+1)*EVR(K) )  &
&              - ( W2R(K:K+1)*W1R(K) + EVR(K:K+1)*W1I(K) )
   AI(K:K+1,K) = AI(K:K+1,K)  &
&              - ( W1I(K:K+1)*W2R(K) - W1R(K:K+1)*EVR(K) )  &
&              - ( EVR(K:K+1)*W1R(K) - W2R(K:K+1)*W1I(K) )
   AR(K+2:N,K) = 0.0
   AI(K+2:N,K) = 0.0
   DO J = K1, N
      AR(J:N,J) = AR(J:N,J)  &
&               - ( W1R(J:N)*W2R(J) + W1I(J:N)*EVR(J) )  &
&               - ( W2R(J:N)*W1R(J) + EVR(J:N)*W1I(J) )
      AI(J:N,J) = AI(J:N,J)  &
&               - ( W1I(J:N)*W2R(J) - W1R(J:N)*EVR(J) )  &
&               - ( EVR(J:N)*W1R(J) - W2R(J:N)*W1I(J) )
   end do


!check
!         WRITE(2,*) 'K=',K
!         DO 9700 L1 = 1, N
!            WRITE(2,9000) (AR(L1,L2),AI(L1,L2),L2=1,N)
! 9700    CONTINUE
! 9000    FORMAT(100('(',2E12.4') '))

end do kdo




!-----------------------------------------------------------------------
!----   Store matrix elements
!----      EVR        :  diagonal elements ( real )
!----    ( W1R, W1I ) :  non-diagonal (i,i+1)-elements ( complex )


DO I = 1, N
   EVR(I) = AR(I,I)
end do
DO I = 1, N-1
   W1R(I) = AR(I+1,I)
   W1I(I) = - AI(I+1,I)
   W2R(I) = W1R(I)*W1R(I) + W1I(I)*W1I(I)
end do


IF( VEC ) then
!-----------------------------------------------------------------------
!----  Set the Unitary matrix elements


DO I = 1, N
   AR(I,I) = 1.0D0
   AI(I,I) = 0.0
end do
DO I = 2, N
   AR(I,1:I-1) = 0.0
   AI(I,1:I-1) = 0.0
end do


DO K = 1, N-2
   KN  = N - K
   KN1 = KN - 1

   AR(KN,KN+1:N) = 0.0
   AI(KN,KN+1:N) = 0.0
   bwork1(KN:N) = AR(KN1,KN:N)
   bwork2(KN:N) = AI(KN1,KN:N)
   work1(KN:N) = 0.0
   work2(KN:N) = 0.0
   DO I = KN, N
      DO J = KN, N
         work1(i) = work1(i)  &
&                 + bwork1(J)*AR(J,I) + bwork2(J)*AI(J,I)
         work2(i) = work2(i)  &
&                 + bwork1(J)*AI(J,I) - bwork2(J)*AR(J,I)
      end do
   end do
   DO I = KN, N
      AR(KN:N,I) = AR(KN:N,I) - bwork1(KN:N)*work1(i) + bwork2(KN:N)*work2(i)
      AI(KN:N,I) = AI(KN:N,I) - bwork1(KN:N)*work2(i) - bwork2(KN:N)*work1(i)
   end do
end do

AR(1,2:N) = 0.0
AI(1,2:N) = 0.0


!check
!         WRITE(2,*) 'U'
!         DO 710 L1 = 1, N
!            WRITE(2,9000) (AR(L1,L2),AI(L1,L2),L2=1,N)
!  710    CONTINUE

end if




!-----------------------------------------------------------------------
!----               Givens' QR method


DMX = 0.0
WMX = 0.0
DO I = 1, N
   DMX = MAX( DMX, ABS(EVR(I)) )
end do
DO I = 1, N-1
   WMX = MAX( WMX, W2R(I) )
end do
WMX = SQRT(WMX)

ESHIFT = MAX( WMX, DMX ) * 1.1D0
EVR(1:N) = EVR(1:N) + ESHIFT

!check
!         WRITE(*,*) 'A_initial'
!         WRITE(*,9002) (EVR(L1),L1=1,N)
!         WRITE(*,9002) (W1R(L1),L1=1,N-1)
!         WRITE(*,9002) (W1I(L1),L1=1,N-1)
! 9002    FORMAT(6E12.4)


KAI   = 0
KNTSF = 0
KNM   = N
evol  = 1.d+10

EVRSFT = 0.0
iterdo: do
   !----   eigenvalue shift
   EVR(1:KNM) = EVR(1:KNM) - EVRSFT


KAI   = KAI   + 1
KNTSF = KNTSF + 1

COLD = 1.0D0
DO K = 2, KNM
   K1 = K - 1
   AA   = EVR(K1)*EVR(K1) + W2R(K1)
   AA   = SQRT(AA)
!           IF( ABS(AA).LT.ZERO ) THEN
!               WRITE(*,*) AA
!               CHECK = 'ZERO '
!               RETURN
!           ENDIF
   DCOSR = EVR(K1)/AA
   DSINR = W1R(K1)/AA
   DSINI = -W1I(K1)/AA

   EK1R = W1R(K1)*COLD
   EK1I = W1I(K1)*COLD

   ARK  = EVR(K)
   W1NR = EK1R*DCOSR + ARK*DSINR
   W1NI = EK1I*DCOSR - ARK*DSINI
   EVR(K) = -(EK1R*DSINR - EK1I*DSINI) + ARK*DCOSR

   BK1 = COLD*DCOSR

   EVR(K1) = AA*BK1 + W1NR*DSINR - W1NI*DSINI
   IF( K.GE.3 ) THEN
       W1R(K-2) = AA*DSOLR
       W1I(K-2) = -AA*DSOLI
   ENDIF

   COLD  = DCOSR
   DSOLR = DSINR
   DSOLI = DSINI
   IF( VEC ) THEN
       DO I = 1, N
          ARK1 = AR(I,K-1)
          AIK1 = AI(I,K-1)
          ARK  = AR(I,K)
          AIK  = AI(I,K)
          AR(I,K-1) = ARK1*DCOSR + ARK*DSINR - AIK*DSINI
          AI(I,K-1) = AIK1*DCOSR + ARK*DSINI + AIK*DSINR
          AR(I,K)   = -(ARK1*DSINR + AIK1*DSINI) + ARK*DCOSR
          AI(I,K)   = -(AIK1*DSINR - ARK1*DSINI) + AIK*DCOSR
       end do
   ENDIF
end do
ARK1 = EVR(KNM)
EVR(KNM)   = ARK1*COLD
W1R(KNM-1) = ARK1*DSOLR
W1I(KNM-1) = -ARK1*DSOLI


!----   diagonal part

EVR(1:KNM) = EVR(1:KNM) + EVRSFT
W2R(1:KNM-1) = W1R(1:KNM-1)*W1R(1:KNM-1) + W1I(1:KNM-1)*W1I(1:KNM-1)

!      WRITE(*,*) 'kaisuu',KAI
!      WRITE(*,9002) (EVR(L1),L1=1,N)
!      WRITE(*,9002) (EVR(L1)-ESHIFT,L1=1,N)
!      WRITE(*,9002) (W2R(L1),L1=1,N-1)
!C         WRITE(*,9002) (W1R(L1),L1=1,N-1)
!C         WRITE(*,9002) (W1I(L1),L1=1,N-1)


!----   hantei

IF( KNTSF.GT.KMAX ) THEN
    CHECK = 'KAISU'
    RETURN
ENDIF
DMX = 0.0
DO I = 1, KNM
   DMX = MAX( DMX, ABS(EVR(I)) )
end do
DMX = DMX*EPSI
DMX = DMX*DMX

IF( ABS(W2R(KNM-1)).LT.DMX .or.  & !KNTSF.GT.KMAX ) THEN
&   kntsf > 2 .and. abs(evol-EVR(KNM)) < zero ) then
    KNM = KNM - 1
    KNTSF = 0
ENDIF
evol = EVR(KNM)

lconvg = .false.
do I = 1, KNM - 1
   lconvg = lconvg .or. ABS(W2R(I)) > DMX
end do
if( .not.lconvg ) exit
   IF( KAI > ITMIN ) then
!       if( abs(EVR(KNM)-EVR(KNM-1)) > zero ) then
           EVRSFT = EVR(KNM)
!       else
!           !---2012/09/04
!           EVRSFT = EVR(KNM) * (1.d0 + 1.d-10)
!       end if
       !---2014/03/09
       do I = 1, KNM - 1 
          if( abs(EVR(KNM)-EVR(i)) < zero ) then
              EVRSFT = EVR(KNM) * (1.d0 + 1.d-10)
              exit
          end if
       end do
   end if
end do iterdo


!      WRITE(*,*) '*** iteration :',KAI
!-----------------------------------------------------------------------

EVR(1:N) = EVR(1:N) - ESHIFT


IF( IFSRT.LE.0 ) RETURN
!-----------------------------------------------------------------------
!----    sorting eigenvalues


W1I(1)= 1

ido: do I = 2, N
   jdo: do J = 1, I-1
   IF( EVR(I).GE.EVR(J) ) cycle jdo

      EVRI = EVR(I)
      do L = I, J+1, -1
         EVR(L) = EVR(L-1)
         W1I(L) = W1I(L-1)
      end do
      EVR(J) = EVRI
      W1I(J) = I
      cycle ido

   end do jdo
   W1I(I) = I
end do ido



IF( .NOT.VEC .OR. IFSRT.LE.1 ) RETURN

!----    sorting eigenvectors


do I = 1, N
   W1R(I) = 0
end do
ido2: do I = 1, N
   IF( NINT(W1R(I)).NE.0 ) cycle ido2
   W1R(I) = 10.0
   IT = NINT(W1I(I))
   IF( I.EQ.IT ) cycle ido2
       IOR = I
       IO  = I

       do
       do L=1,N
          ARI = AR(L,IO)
          AII = AI(L,IO)
          AR(L,IO)  = AR(L,IT)
          AI(L,IO)  = AI(L,IT)
          AR(L,IT) = ARI
          AI(L,IT) = AII
       end do
       IO = IT
       W1R(IO) = 10.0
       IT = NINT(W1I(IO))
       IF( IT.EQ.IOR ) cycle ido2
       end do

end do ido2


RETURN
END




subroutine vsgar( n, key, ix, maxn, dstkey, x1, x2, buff )
integer n, key(maxn), ix(maxn)
integer dstkey(maxn), x1(maxn), x2(maxn), buff(maxn)
integer minkey, maxkey, dstmax, keylen, mask, p
logical x1flg
!      equivalence (x2,buff)


if( n.le.1 ) return

dstmax = n - 1
minkey = key(1)
maxkey = key(1)
do i = 2, n
   minkey = min( minkey, key(i) )
   maxkey = max( maxkey, key(i) )
end do

!---error trap
if( maxkey.eq.minkey ) then
    maxkey = maxkey + 1
    key(n) = key(n) + 1
end if

do i = 1, n
   dstkey(i) = int( float(key(i) - minkey)  &
&               * ( float(n-1)/float(maxkey-minkey)*0.5 ) )
end do
keylen = 0
do
   if( dstmax.eq.0 ) exit
   keylen = keylen + 1
   dstmax = dstmax/2
end do
do i = 1, n
   x1(i) = i
end do
x1flg = .true.

do j = 1, keylen
   mask = 2**(j-1)
   p = 0
   if( x1flg ) then
       do i = 1, n
          if( iand(dstkey(x1(i)),mask).eq.0 ) then
              p = p + 1
              x2(p) = x1(i)
          end if
       end do
       do i = 1, n
          if( iand(dstkey(x1(i)),mask).eq.mask ) then
              p = p + 1
              x2(p) = x1(i)
          end if
       end do
     else
       do i = 1, n
          if( iand(dstkey(x2(i)),mask).eq.0 ) then
              p = p + 1
              x1(p) = x2(i)
          end if
       end do
       do i = 1, n
          if( iand(dstkey(x2(i)),mask).eq.mask ) then
              p = p + 1
              x1(p) = x2(i)
          end if
       end do
   end if
#ifdef CRAY_T3E
!---     dummy if statement to prevent optimization
   if( j.eq.0 ) write(*,*) j,p,(x1(i),x2(i),i=1,n)
#endif
   x1flg = .not.x1flg
end do

if( .not.x1flg ) then
    do i = 1, n
       x1(i) = x2(i)
    end do
end if
do i = 1, n
   buff(i) = key(i)
end do
do i = 1, n
   key(i) = buff(x1(i))
end do
do i = 1, n
   buff(i) = ix(i)
end do
do i = 1, n
   ix(i) = buff(x1(i))
end do

return
end




!----------------------------------------------------------------------
!     Odd-even transposition sort (vector)
!             from  'Iwanami Koza Software Science 9'  ISBN4-00-010349-0
!
!     n      : No. of records
!     key    : key of records
!     ix     : data pointer of record
!-----------------------------------------------------------------------
subroutine vsoe( n, key, ix, maxn )
integer n, key(maxn), ix(maxn)
integer mode, notrs, i, irev, buff


if( n.le.1 ) return

mode = 1
notrs = 0

do
   irev = 0
   do i = mode, n-1, 2
      if( key(i).gt.key(i+1) ) irev = irev + 1
   end do
   if( irev.eq.0 ) then
       notrs = notrs + 1
     else
       notrs = 0
       do i = mode, n-1, 2
          if( key(i).gt.key(i+1) )  then
              buff     = key(i)
              key(i)   = key(i+1)
              key(i+1) = buff
              buff     = ix(i)
              ix(i)    = ix(i+1)
              ix(i+1)  = buff
          end if
       end do
   end if
   if( mode.eq.1 ) then
       mode = 2
     else
       mode = 1
   end if
   if( notrs >= 2 ) exit
end do


return
end




!      DOUBLE PRECISION FUNCTION DERFNC( CD, ERERR )
!C-----------------------------------------------------------------------
!C     complementary error function   -- Simpson's quadrature --
!c                                                       date   '90/11/09
!C-----------------------------------------------------------------------
!C    CD    ......  maximum length of integration
!C    ERERR ......  accuracy of integration
!C-----------------------------------------------------------------------
!      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!      parameter( CNS1    = 0.3333333333333D0,
!     &           CNS2    = 1.333333333333D0,
!     &           CNS3    = 0.6666666666667D0,
!     &           PAIRUB  = 1.128379167096d0  )
!c      F(X)= DEXP(-X*X)
!
!         NST  = 0
!         WAKE = 1.0D0
!         HDEL = CD/WAKE
!c         ERFS = F(0.0D0) + F(CD)
!         ERFS = 1.0D0 + exp(-CD*CD)
!         ERF1 = HDEL*ERFS/2.0D0
!         ER2  = 0.0
!         do
!         NST  = NST + 1
!         WAKE = WAKE*2.0D0
!         HDEL = CD/WAKE
!         NSUM = WAKE
!         NSUM = NSUM/2
!         ER4  = 0.0
!         DO 420 I=1,NSUM
!            J   = 2*I - 1
!            DDI = DBLE(J)*HDEL
!c            ER4 = ER4 + F(DDI)
!            ER4 = ER4 + exp(-DDI*DDI)
!  420    CONTINUE
!         ERF2 = HDEL*(CNS1*ERFS + CNS2*ER4 + CNS3*ER2)
!         ER2  = ER2 + ER4
!         DSA   = ABS( ERF2 - ERF1 )
!         ERF1  = ERF2
!         IF( DSA < ERERR ) exit
!         end do
!         DERFNC = 1.0D0 - PAIRUB*ERF1
! 
!      RETURN
!      END




DOUBLE PRECISION FUNCTION derfnc( x, ererr )
!-----------------------------------------------------------------------
!     complementary error function
!                                           NUMERICAL RECIPES  Chapter 6
!-----------------------------------------------------------------------
implicit real*8 (a-h,o-z)

if( x.lt.0.d0 ) then
    derfnc = 1.d0 + gammp( 0.5d0, x*x, ererr, ierr )
  else
    derfnc = gammq( 0.5d0, x*x, ererr, ierr )
end if

return
end




DOUBLE PRECISION FUNCTION gammp( a, x, ererr, ierr )
!-----------------------------------------------------------------------
!     incomplete gamma function P(a,x)
!                                           NUMERICAL RECIPES  Chapter 6
!-----------------------------------------------------------------------
implicit real*8 (a-h,o-z)

if( x.lt.a+1.d0 ) then
!--- use the series representation
    call gser( gamser, a, x, gln, ererr, ierr )
    gammp = gamser
  else
!--- use the continued fraction representation
    call gcf( gammcf, a, x, gln, ererr, ierr )
    gammp = 1.d0 - gammcf
end if

return
end




DOUBLE PRECISION FUNCTION gammq( a, x, ererr, ierr )
!-----------------------------------------------------------------------
!     incomplete gamma function Q(a,x) = 1 - P(a,x)
!                                           NUMERICAL RECIPES  Chapter 6
!-----------------------------------------------------------------------
implicit real*8 (a-h,o-z)

if( x.lt.a+1.d0 ) then
!--- use the series representation
    call gser( gamser, a, x, gln, ererr, ierr )
    gammq = 1.d0 - gamser
  else
!--- use the continued fraction representation
    call gcf( gammcf, a, x, gln, ererr, ierr )
    gammq = gammcf
end if

return
end




subroutine gser( gamser, a, x, gln, eps, ierr )
!-----------------------------------------------------------------------
!     incomplete gamma function P(a,x)
!     evaluated by its series representation
!                                           NUMERICAL RECIPES  Chapter 6
!-----------------------------------------------------------------------
implicit real*8 (a-h,o-z)
parameter( itmax = 100 )

ierr = 0
gln = gammln( a )
if( x.le.0.d0 ) then
    gamser = 0.d0
    return
end if
ap = a
sum = 1.d0/a
del = sum
ierr = 1
do n = 1, itmax
   ap = ap + 1.d0
   del = del * x/ap
   sum = sum + del
   if( abs(del).lt.abs(sum)*eps ) then
       ierr = 0
       exit
   end if
end do
gamser = sum*exp(-x+a*log(x)-gln)

return
end




subroutine gcf( gammcf, a, x, gln, eps, ierr )
!-----------------------------------------------------------------------
!     incomplete gamma function Q(a,x)
!     evaluated by its continued fraction representation
!                                           NUMERICAL RECIPES  Chapter 6
!-----------------------------------------------------------------------
implicit real*8 (a-h,o-z)
parameter( itmax = 100 )

ierr = 0
gln = gammln( a )
gold = 0.d0
a0 = 1.d0
a1 = x
b0 = 0.d0
b1 = 1.d0
fac = 1.d0
ierr = 1
do n = 1, itmax
   an = dble(n)
   ana = an - a
   a0 = ( a1 + a0*ana )*fac
   b0 = ( b1 + b0*ana )*fac
   anf = an*fac
   a1 = x*a0 + anf*a1
   b1 = x*b0 + anf*b1
   if( a1.ne.0.d0 ) then
       fac = 1.d0/a1
       g = b1*fac
       if( abs((g-gold)/g).lt.eps ) then
           ierr = 0
           exit
       end if
       gold = g
   end if
end do
gammcf = exp(-x+a*log(x)-gln)*g

return
end




DOUBLE PRECISION FUNCTION gammln( xx )
!-----------------------------------------------------------------------
!     log(gamma function(xx)) for xx > 0.
!     Full accuracy is obtained for xx > 1.
!                                           NUMERICAL RECIPES  Chapter 6
!-----------------------------------------------------------------------
implicit real*8 (a-h,o-z)
dimension cof(6)
data cof, stp / 76.18009173d0, -86.50532033d0, 24.01409822d0,  &
& -1.231739516d0, .120858003d-02, -.536382d-5,  2.50662827465d0 /
data half, one, fpf / 0.5d0, 1.d0, 5.5d0 /

x = xx - one
tmp = x + fpf
tmp = ( x + half )*log(tmp) - tmp
ser = one
do j = 1, 6
   x = x + one
   ser = ser + cof(j)/x
end do
gammln = tmp + log(stp*ser)

return
end




subroutine rndgaus( rn )
!----------------------------------------------------------------------
!-- random number ( Gaussian distribution exp(-rn*rn/2)/sqrt(2*pi) ) --
!----------------------------------------------------------------------
implicit none
real*8  :: rn

real*8  :: gset, rr, v1, v2, rsq, fac
integer :: iset = 0
real*8  :: seed = 25.d0
save seed, iset, gset

if( iset == 0 ) then
    do
       call rnd00( rr, seed )
       v1 = 2.d0*rr - 1.d0
       call rnd00( rr, seed )
       v2 = 2.d0*rr - 1.d0
       rsq = v1*v1 + v2*v2
       if( rsq < 1.d0 ) exit
    end do

    fac = sqrt(-2.d0*log(rsq)/rsq)
    gset = v1*fac
    rn   = v2*fac
    iset = 1
else
    rn = gset
    iset = 0
end if


return
end




subroutine nrmdst( am, ds, rn )
!--- random number ( am : average,  ds : dispersion )  -----------------
implicit none
real*8  :: am, ds, rn

real*8  :: a, rr
integer :: i
real*8  :: seed = 25.d0
save seed

a = 0.d0
do i = 1, 12
   call rnd00( rr, seed )
!   rr = (rr + 1.d0) * .5d0
   a = a + rr
end do
rn = (a - 6.d0)*ds + am

return
end




subroutine rndexp( rn, rambda )
!----------------------------------------------------------------------
!-- random number ( Exponential distribution ) --
!
!   f(r) = rambda*exp(-rambda*r)
!
!   Average = 1/rambda
!----------------------------------------------------------------------
implicit none
real*8  :: rn
real*8  :: rambda

real*8  :: rr
real*8  :: seed = 25.d0
save seed

call rnd00( rr, seed )
rn = -log(1.d0 - rr)/rambda

return
end




subroutine rnd00( rn, rndx )
!--- uniform random number ( 0. < RN < 1. )  --------------------------
implicit none
real*8  :: rn, rndx

! real*8  :: rndx = 25d0
real*8  :: rnda = 16807d0
real*8  :: rndm = 0.2147483647D+10
real*8  :: rnd
save rnda, rndm

  rnd  = rnda*rndx
  rndx = dmod(rnd,rndm)
  rn    = rndx/rndm

return
end




subroutine dexchange( a, b, n )
!-----exchange a and b
implicit none
integer :: n
real*8,  dimension(n) :: a, b
integer :: i
real*8  :: c

do i = 1, n
   c    = a(i)
   a(i) = b(i)
   b(i) = c
end do

return
end subroutine




subroutine dsum_a_by_b( s, a, b, n )
!-----sum of a*b
implicit none
integer :: n
real*8  :: s
real*8,  dimension(n) :: a, b
integer :: i

s = 0.d0
do i = 1, n
   s = s + a(i)*b(i)
end do

return
end subroutine




subroutine dadd_b_to_a( a, b, n )
!-----add a and b
implicit none
integer :: n
real*8,  dimension(n) :: a, b

a(1:n) = a(1:n) + b(1:n)

return
end subroutine




subroutine dcopy_a_to_b( a, b, n )
!-----copy a to b
implicit none
integer :: n
real*8,  dimension(n) :: a, b

b(1:n) = a(1:n)

return
end subroutine




subroutine dmulti_b_to_a( a, b, n )
!-----copy a to b
implicit none
integer :: n
real*8,  dimension(n) :: a, b

a(1:n) = a(1:n)*b(1:n)

return
end subroutine




DOUBLE PRECISION FUNCTION vecratio( u, v )
implicit none
real*8,  dimension(3) :: u, v

vecratio = sqrt(u(1)*u(1) + u(2)*u(2) + u(3)*u(3))  &
&        / sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))

return
end function




subroutine matinv33( a, b )
!-----------------------------------------------------------------------
!  ( input )
!     a     ...... 3x3 matrix
!  ( output )
!     b     ...... inverse of a
!-----------------------------------------------------------------------
!   i.e.   sum_k a(i,k) * b(k,j) = delta_ij 
!-----------------------------------------------------------------------
implicit none
real*8,  dimension(3,3) :: a, b

!------declare local variables
real*8  :: v, vr


b(1,1)  = a(2,2)*a(3,3) - a(3,2)*a(2,3)
b(1,2)  = a(3,2)*a(1,3) - a(1,2)*a(3,3)
b(1,3)  = a(1,2)*a(2,3) - a(2,2)*a(1,3)
v  = a(1,1)*b(1,1) + a(2,1)*b(1,2) + a(3,1)*b(1,3)
vr = 1.0d0/v


b(1,1) = b(1,1)*vr
b(1,2) = b(1,2)*vr
b(1,3) = b(1,3)*vr

b(2,1) = (a(2,3)*a(3,1) - a(3,3)*a(2,1))*vr
b(2,2) = (a(3,3)*a(1,1) - a(1,3)*a(3,1))*vr
b(2,3) = (a(1,3)*a(2,1) - a(2,3)*a(1,1))*vr

b(3,1) = (a(2,1)*a(3,2) - a(3,1)*a(2,2))*vr
b(3,2) = (a(3,1)*a(1,2) - a(1,1)*a(3,2))*vr
b(3,3) = (a(1,1)*a(2,2) - a(2,1)*a(1,2))*vr


return
end subroutine
