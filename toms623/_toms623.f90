!     ALGORITHM 623 COLLECTED ALGORITHMS FROM ACM.
!     ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.10, NO. 4,
!     DEC., 1984, P. 437.
      SUBROUTINE ADNODE (KK,X,Y,Z, IADJ,IEND, IER)
      INTEGER KK, IADJ(1), IEND(KK), IER, INDX
      DOUBLE PRECISION    X(KK), Y(KK), Z(KK)
      LOGICAL SWPTST
      EXTERNAL INDX
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE ADDS NODE KK TO A TRIANGULATION OF THE
! CONVEX HULL OF NODES 1,...,KK-1, PRODUCING A TRIANGULATION
! OF THE CONVEX HULL OF NODES 1,...,KK.  A SEQUENCE OF EDGE
! SWAPS IS THEN APPLIED TO THE MESH, RESULTING IN AN OPTIMAL
! TRIANGULATION.  ADNODE IS PART OF AN INTERPOLATION PACKAGE
! WHICH ALSO PROVIDES ROUTINES TO INITIALIZE THE DATA STRUC-
! TURE, PLOT THE MESH, AND DELETE ARCS.
!
! INPUT PARAMETERS -   KK - INDEX OF THE NODE TO BE ADDED
!                           TO THE MESH.  KK .GE. 4.
!
!                   X,Y,Z - VECTORS OF LENGTH .GE. KK CON-
!                           TAINING CARTESIAN COORDINATES
!                           OF THE NODES.  (X(I),Y(I),Z(I))
!                           DEFINES NODE I FOR I = 1,...,KK.
!
!                    IADJ - SET OF ADJACENCY LISTS OF NODES
!                           1,...,KK-1.
!
!                    IEND - POINTERS TO THE ENDS OF
!                           ADJACENCY LISTS IN IADJ FOR
!                           EACH NODE IN THE MESH.
!
! IADJ AND IEND MAY BE CREATED BY TRMESH.
!
! KK, X, Y, AND Z ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS - IADJ,IEND - UPDATED WITH THE ADDITION
!                                 OF NODE KK AS THE LAST
!                                 ENTRY.
!
!                           IER - ERROR INDICATOR
!                                 IER = 0 IF NO ERRORS
!                                         WERE ENCOUNTERED.
!                                 IER = 1 IF ALL NODES
!                                         (INCLUDING KK) ARE
!                                         COLLINEAR.
!
! MODULES REFERENCED BY ADNODE - TRFIND, INTADD, BDYADD,
!                                COVSPH, SHIFTD, INDX,
!                                SWPTST, SWAP
!
!***********************************************************
!
      INTEGER K, KM1, I1, I2, I3, INDKF, INDKL, NABOR1, &
              IO1, IO2, IN1, INDK1, IND2F, IND21
      DOUBLE PRECISION    P(3), DUM
!
! LOCAL PARAMETERS -
!
! K =        LOCAL COPY OF KK
! KM1 =      K - 1
! I1,I2,I3 = VERTICES OF A TRIANGLE CONTAINING K
! INDKF =    IADJ INDEX OF THE FIRST NEIGHBOR OF K
! INDKL =    IADJ INDEX OF THE LAST NEIGHBOR OF K
! NABOR1 =   FIRST NEIGHBOR OF K BEFORE ANY SWAPS OCCUR
! IO1,IO2 =  ADJACENT NEIGHBORS OF K DEFINING AN ARC TO
!              BE TESTED FOR A SWAP
! IN1 =      VERTEX OPPOSITE K -- FIRST NEIGHBOR OF IO2
!              WHICH PRECEDES IO1.  IN1,IO1,IO2 ARE IN
!              COUNTERCLOCKWISE ORDER.
! INDK1 =    INDEX OF IO1 IN THE ADJACENCY LIST FOR K
! IND2F =    INDEX OF THE FIRST NEIGHBOR OF IO2
! IND21 =    INDEX OF IO1 IN THE ADJACENCY LIST FOR IO2
! P =        CARTESIAN COORDINATES OF NODE KK
! DUM =      DUMMY PARAMETER FOR CALL TO TRFIND
!
      IER = 0
      K = KK
!
! INITIALIZATION
!
      KM1 = K - 1
      P(1) = X(K)
      P(2) = Y(K)
      P(3) = Z(K)
!
! ADD NODE K TO THE MESH
!
      CALL TRFIND(KM1,P,X,Y,Z,IADJ,IEND, DUM,DUM,DUM, &
                  I1,I2,I3)
      IF (I1 .EQ. 0) GO TO 5
      IF (I3 .NE. 0) CALL INTADD(K,I1,I2,I3, IADJ,IEND )
      IF (I3 .EQ. 0  .AND.  I1 .NE. I2) &
        CALL BDYADD (K,I1,I2, IADJ,IEND )
      IF (I1 .EQ. I2) CALL COVSPH(K,I1, IADJ,IEND )
!
! INITIALIZE VARIABLES FOR OPTIMIZATION OF THE MESH
!
      INDKF = IEND(KM1) + 1
      INDKL = IEND(K)
      NABOR1 = IADJ(INDKF)
      IO2 = NABOR1
      INDK1 = INDKF + 1
      IO1 = IADJ(INDK1)
!
! BEGIN LOOP -- FIND THE VERTEX OPPOSITE K
!
    1 IND2F = 1
      IF (IO2 .NE. 1) IND2F = IEND(IO2-1) + 1
      IND21 = INDX(IO2,IO1,IADJ,IEND)
      IF (IND2F .EQ. IND21) GO TO 2
      IN1 = IADJ(IND21-1)
      GO TO 3
!
! IN1 IS THE LAST NEIGHBOR OF IO2
!
    2 IND21 = IEND(IO2)
      IN1 = IADJ(IND21)
      IF (IN1 .EQ. 0) GO TO 4
!
! SWAP TEST -- IF A SWAP OCCURS, TWO NEW ARCS ARE OPPOSITE K
!              AND MUST BE TESTED.  INDK1 AND INDKF MUST BE
!              DECREMENTED.
!
    3 IF ( .NOT. SWPTST(IO1,IO2,IN1,K,X,Y,Z) ) GO TO 4
      CALL SWAP(IN1,K,IO1,IO2, IADJ,IEND )
      IO1 = IN1
      INDK1 = INDK1 - 1
      INDKF = INDKF - 1
      GO TO 1
!
! NO SWAP OCCURRED.  RESET IO2 AND IO1, AND TEST FOR
!   TERMINATION.
!
    4 IF (IO1 .EQ. NABOR1) RETURN
      IO2 = IO1
      INDK1 = INDK1 + 1
      IF (INDK1 .GT. INDKL) INDK1 = INDKF
      IO1 = IADJ(INDK1)
      IF (IO1 .NE. 0) GO TO 1
      RETURN
!
! ALL NODES ARE COLLINEAR
!
    5 IER = 1
      RETURN
      END
      SUBROUTINE APLYR (X,Y,Z,CX,SX,CY,SY, XP,YP,ZP)
      DOUBLE PRECISION X, Y, Z, CX, SX, CY, SY, XP, YP, ZP
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE APPLIES THE ROTATION R DEFINED BY SUB-
! ROUTINE CONSTR TO THE UNIT VECTOR (X Y Z)**T, I.E. (X,Y,Z)
! IS ROTATED TO (XP,YP,ZP).  IF (XP,YP,ZP) LIES IN THE
! SOUTHERN HEMISPHERE (ZP .LT. 0), (XP,YP) ARE SET TO THE
! COORDINATES OF THE NEAREST POINT OF THE EQUATOR, ZP RE-
! MAINING UNCHANGED.
!
! INPUT PARAMETERS - X,Y,Z - COORDINATES OF A POINT ON THE
!                            UNIT SPHERE.
!
!              CX,SX,CY,SY - ELEMENTS OF THE ROTATION DE-
!                            FINED BY CONSTR.
!
! INPUT PARAMETERS ARE NOT ALTERED EXCEPT AS NOTED BELOW.
!
! OUTPUT PARAMETERS - XP,YP,ZP - COORDINATES OF THE ROTATED
!                                POINT ON THE SPHERE UNLESS
!                                ZP .LT. 0, IN WHICH CASE
!                                (XP,YP,0) IS THE CLOSEST
!                                POINT OF THE EQUATOR TO THE
!                                ROTATED POINT.  STORAGE FOR
!                                XP, YP, AND ZP MAY COINCIDE
!                                WITH STORAGE FOR X, Y, AND
!                                Z, RESPECTIVELY, IF THE
!                                LATTER NEED NOT BE SAVED.
!
! MODULES REFERENCED BY APLYR - NONE
!
! INTRINSIC FUNCTION CALLED BY APLYR - SQRT
!
!***********************************************************
!
      DOUBLE PRECISION T
!
! LOCAL PARAMETER -
!
! T = TEMPORARY VARIABLE
!
      T = SX*Y + CX*Z
      YP = CX*Y - SX*Z
      ZP = SY*X + CY*T
      XP = CY*X - SY*T
      IF (ZP .GE. 0.) RETURN
!
! MOVE (XP,YP,ZP) TO THE EQUATOR
!
      T = SQRT(XP*XP + YP*YP)
      IF (T .EQ. 0.) GO TO 1
      XP = XP/T
      YP = YP/T
      RETURN
!
! MOVE THE SOUTH POLE TO AN ARBITRARY POINT OF THE EQUATOR
!
    1 XP = 1.
      YP = 0.
      RETURN
      END
      SUBROUTINE APLYRT (G1P,G2P,CX,SX,CY,SY, G)
      DOUBLE PRECISION G1P, G2P, CX, SX, CY, SY, G(3)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE APPLIES THE INVERSE (TRANSPOSE) OF THE
! ROTATION DEFINED BY SUBROUTINE CONSTR TO THE VECTOR
! (G1P G2P 0)**T, I.E. THE GRADIENT (G1P,G2P,0) IN THE ROT-
! ATED COORDINATE SYSTEM IS MAPPED TO (G1,G2,G3) IN THE
! ORIGINAL COORDINATE SYSTEM.
!
! INPUT PARAMETERS - G1P,G2P - X- AND Y-COMPONENTS, RESPECT-
!                              IVELY, OF THE GRADIENT IN THE
!                              ROTATED COORDINATE SYSTEM.
!
!                CX,SX,CY,SY - ELEMENTS OF THE ROTATION R
!                              CONSTRUCTED BY SUBROUTINE
!                              CONSTR.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS - G - X-, Y-, AND Z-COMPONENTS (IN THAT
!                         ORDER) OF THE INVERSE ROTATION
!                         APPLIED TO (G1P,G2P,0) -- GRADIENT
!                         IN THE ORIGINAL COORDINATE SYSTEM.
!
! MODULES REFERENCED BY APLYRT - NONE
!
!***********************************************************
!
      DOUBLE PRECISION T
!
! LOCAL PARAMETERS -
!
! T = TEMPORARY VARIABLE
!
      T = SY*G1P
      G(1) = CY*G1P
      G(2) = CX*G2P - SX*T
      G(3) = -SX*G2P - CX*T
      RETURN
      END
      SUBROUTINE ARCINT (P,P1,P2,W1,W2,G1,G2, W,G,GN)
      DOUBLE PRECISION    P(3), P1(3), P2(3), W1, W2, G1(3), G2(3), &
              W, G(3), GN
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN 3 POINTS P, P1, AND P2 LYING ON A COMMON GEODESIC
! OF THE UNIT SPHERE WITH P BETWEEN P1 AND P2, ALONG WITH
! DATA VALUES AND GRADIENTS AT P1 AND P2, THIS SUBROUTINE
! COMPUTES AN INTERPOLATED VALUE W AND A GRADIENT VECTOR G
! AT P.  W IS COMPUTED BY HERMITE CUBIC INTERPOLATION REL-
! ATIVE TO ARC-LENGTH ALONG THE GEODESIC.  THE TANGENTIAL
! COMPONENT OF G IS THE DERIVATIVE (WITH RESPECT TO ARC-
! LENGTH) OF THE CUBIC INTERPOLANT AT P, WHILE THE NORMAL
! COMPONENT OF G IS OBTAINED BY LINEAR INTERPOLATION OF THE
! NORMAL COMPONENTS OF THE GRADIENTS AT P1 AND P2.  THIS
! ALGORITHM IS DUE TO C. L. LAWSON.
!
! INPUT PARAMETERS - P - CARTESIAN COORDINATES OF A POINT
!                        LYING ON THE ARC DEFINED BY P1 AND
!                        P2.
!
!                P1,P2 - COORDINATES OF DISTINCT POINTS ON
!                        THE UNIT SPHERE DEFINING AN ARC
!                        WITH LENGTH LESS THAN 180 DEGREES.
!
!                W1,W2 - DATA VALUES ASSOCIATED WITH P1 AND
!                        P2, RESPECTIVELY.
!
!                G1,G2 - GRADIENT VECTORS ASSOCIATED WITH P1
!                        AND P2.  G1 AND G2 ARE ORTHOGONAL
!                        TO P1 AND P2, RESPECTIVELY.
!
! THE ABOVE PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
!                    G - ARRAY OF LENGTH 3.
!
! OUTPUT PARAMETERS - W - INTERPOLATED VALUE AT P.
!
!                     G - INTERPOLATED GRADIENT AT P.
!
!                    GN - NORMAL COMPONENT OF G WITH THE
!                         DIRECTION P1 X P2 TAKEN TO BE
!                         POSITIVE.  THE EXTRAPOLATION
!                         PROCEDURE REQUIRES THIS COMPONENT.
!
! FOR EACH VECTOR V, V(1), V(2), AND V(3) CONTAIN X-, Y-,
!   AND Z-COMPONENTS, RESPECTIVELY.
!
! MODULES REFERENCED BY ARCINT - ARCLEN
!
! INTRINSIC FUNCTION CALLED BY ARCINT - SQRT
!
!***********************************************************
!
      INTEGER I, LUN
      DOUBLE PRECISION    UN(3), UNORM, TAU1, TAU2, A, AL, S, T, GT,&
              ARCLEN
      DATA    LUN/6/
!
! LOCAL PARAMETERS -
!
! I =         DO-LOOP INDEX
! LUN =       LOGICAL UNIT FOR ERROR MESSAGES
! UN =        UNIT NORMAL TO THE PLANE OF P, P1, AND P2
! UNORM =     EUCLIDEAN NORM OF P1 X P2 -- USED TO NORMALIZE
!               UN
! TAU1,TAU2 = TANGENTIAL DERIVATIVES (COMPONENTS OF G1,G2)
!               AT P1 AND P2
! A =         ANGLE IN RADIANS (ARC-LENGTH) BETWEEN P1 AND
!               P2
! AL =        ARC-LENGTH BETWEEN P1 AND P
! S =         NORMALIZED VALUE OF AL -- AS P VARIES FROM P1
!               TO P2, S VARIES FROM 0 TO 1
! T =         1-S -- S AND T ARE BARYCENTRIC COORDINATES OF
!               P WITH RESPECT TO THE ARC FROM P1 TO P2
! GT =        TANGENTIAL COMPONENT OF G -- COMPONENT IN THE
!               DIRECTION UN X P
!
!
! COMPUTE UNIT NORMAL UN
!
      UN(1) = P1(2)*P2(3) - P1(3)*P2(2)
      UN(2) = P1(3)*P2(1) - P1(1)*P2(3)
      UN(3) = P1(1)*P2(2) - P1(2)*P2(1)
      UNORM = SQRT(UN(1)*UN(1) + UN(2)*UN(2) + UN(3)*UN(3))
      IF (UNORM .EQ. 0.) GO TO 2
!
! NORMALIZE UN
!
      DO 1 I = 1,3
    1   UN(I) = UN(I)/UNORM
!
! COMPUTE TANGENTIAL DERIVATIVES AT THE ENDPOINTS --
!   TAU1 = (G1,UN X P1) = (G1,P2)/UNORM AND
!   TAU2 = (G2,UN X P2) = -(G2,P1)/UNORM.
!
      TAU1 = (G1(1)*P2(1) + G1(2)*P2(2) + G1(3)*P2(3))/UNORM
      TAU2 =-(G2(1)*P1(1) + G2(2)*P1(2) + G2(3)*P1(3))/UNORM
!
! COMPUTE ARC-LENGTHS A, AL
!
      A = ARCLEN(P1,P2)
      IF (A .EQ. 0.) GO TO 2
      AL = ARCLEN(P1,P)
!
! COMPUTE W BY HERMITE CUBIC INTERPOLATION
!
      S = AL/A
      T = 1. - S
      W = W1*(2.*S+1.)*T*T + W2*(3.-2.*S)*S*S + &
          A*S*T*(TAU1*T - TAU2*S)
!
! COMPUTE TANGENTIAL AND NORMAL DERIVATIVES AT P
!
      GT = 6.*S*T*(W2-W1)/A + &
           TAU1*T*(1.-3.*S) + TAU2*S*(3.*S-2.)
      GN = T*(UN(1)*G1(1) + UN(2)*G1(2) + UN(3)*G1(3)) + &
           S*(UN(1)*G2(1) + UN(2)*G2(2) + UN(3)*G2(3))
!
! COMPUTE G = GT*(UN X P) + GN*UN
!
      G(1) = GT*(UN(2)*P(3) - UN(3)*P(2)) + GN*UN(1)
      G(2) = GT*(UN(3)*P(1) - UN(1)*P(3)) + GN*UN(2)
      G(3) = GT*(UN(1)*P(2) - UN(2)*P(1)) + GN*UN(3)
      RETURN
!
! P1 X P2 = 0.  PRINT AN ERROR MESSAGE AND TERMINATE
!   PROCESSING.
!
    2 WRITE (LUN,100) (P1(I),I=1,3), (P2(I),I=1,3)
  100 FORMAT (1H1,24HERROR IN ARCINT -- P1 = ,2(F9.6,3H,  ), &
              F9.6/1H ,19X,5HP2 = ,2(F9.6,3H,  ),F9.6)
      STOP
      END
      DOUBLE PRECISION FUNCTION ARCLEN (P,Q)
      DOUBLE PRECISION P(3), Q(3)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS FUNCTION COMPUTES THE ARC-LENGTH (ANGLE IN RADIANS)
! BETWEEN A PAIR OF POINTS ON THE UNIT SPHERE.
!
! INPUT PARAMETERS - P,Q - VECTORS OF LENGTH 3 CONTAINING
!                          THE X-, Y-, AND Z-COORDINATES (IN
!                          THAT ORDER) OF POINTS ON THE UNIT
!                          SPHERE.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS FUNCTION.
!
! OUTPUT PARAMETER - ARCLEN - ANGLE IN RADIANS BETWEEN THE
!                             UNIT VECTORS P AND Q.  0 .LE.
!                             ARCLEN .LE. PI.
!
! MODULES REFERENCED BY ARCLEN - NONE
!
! INTRINSIC FUNCTIONS CALLED BY ARCLEN - ATAN, SQRT
!
!***********************************************************
!
      INTEGER I
      DOUBLE PRECISION    D
!
! LOCAL PARAMETERS -
!
! I = DO-LOOP INDEX
! D = EUCLIDEAN NORM SQUARED OF P+Q
!
      D = 0.
      DO 1 I = 1,3
    1   D = D + (P(I) + Q(I))**2
      IF (D .EQ. 0.) GO TO 2
      IF (D .GE. 4.) GO TO 3
      ARCLEN = 2.*ATAN(SQRT((4.-D)/D))
      RETURN
!
! P AND Q ARE SEPARATED BY 180 DEGREES
!
    2 ARCLEN = 4.*ATAN(1.)
      RETURN
!
! P AND Q COINCIDE
!
    3 ARCLEN = 0.
      RETURN
      END
      SUBROUTINE BDYADD (KK,I1,I2, IADJ,IEND )
      INTEGER KK, I1, I2, IADJ(1), IEND(KK)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE ADDS A BOUNDARY NODE TO A TRIANGULATION OF
! A SET OF KK-1 POINTS ON THE UNIT SPHERE.  IADJ AND IEND
! ARE UPDATED WITH THE INSERTION OF NODE KK.
!
! INPUT PARAMETERS -   KK - INDEX OF AN EXTERIOR NODE TO BE
!                           ADDED.  KK .GE. 4.
!
!                      I1 - FIRST (RIGHTMOST AS VIEWED FROM
!                           KK) BOUNDARY NODE IN THE MESH
!                           WHICH IS VISIBLE FROM KK - THE
!                           LINE SEGMENT KK-I1 INTERSECTS
!                           NO ARCS.
!
!                      I2 - LAST (LEFTMOST) BOUNDARY NODE
!                           WHICH IS VISIBLE FROM KK.
!
!                    IADJ - SET OF ADJACENCY LISTS OF NODES
!                           IN THE MESH.
!
!                    IEND - POINTERS TO THE ENDS OF
!                           ADJACENCY LISTS IN IADJ FOR
!                           EACH NODE IN THE MESH.
!
!   IADJ AND IEND MAY BE CREATED BY TRMESH AND MUST CONTAIN
! THE VERTICES I1 AND I2.  I1 AND I2 MAY BE DETERMINED BY
! TRFIND.
!
! KK, I1, AND I2 ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS - IADJ,IEND - UPDATED WITH THE ADDITION
!                                 OF NODE KK AS THE LAST
!                                 ENTRY.  NODE KK WILL BE
!                                 CONNECTED TO I1, I2, AND
!                                 ALL BOUNDARY NODES BETWEEN
!                                 THEM.  NO OPTIMIZATION OF
!                                 THE MESH IS PERFORMED.
!
! MODULE REFERENCED BY BDYADD - SHIFTD
!
! INTRINSIC FUNCTIONS CALLED BY BDYADD - MIN0, MAX0
!
!***********************************************************
!
      INTEGER K, KM1, NRIGHT, NLEFT, NF, NL, N1, N2, I, &
              IMIN, IMAX, KEND, NEXT, INDX
!
! LOCAL PARAMETERS -
!
! K =            LOCAL COPY OF KK
! KM1 =          K - 1
! NRIGHT,NLEFT = LOCAL COPIES OF I1, I2
! NF,NL =        INDICES OF IADJ BOUNDING THE PORTION OF THE
!                  ARRAY TO BE SHIFTED
! N1 =           IADJ INDEX OF THE FIRST NEIGHBOR OF NLEFT
! N2 =           IADJ INDEX OF THE LAST NEIGHBOR OF NRIGHT
! I =            DO-LOOP INDEX
! IMIN,IMAX =    BOUNDS ON DO-LOOP INDEX -- FIRST AND LAST
!                  ELEMENTS OF IEND TO BE INCREMENTED
! KEND =         POINTER TO THE LAST NEIGHBOR OF K IN IADJ
! NEXT =         NEXT BOUNDARY NODE TO BE CONNECTED TO KK
! INDX =         INDEX FOR IADJ
!
      K = KK
      KM1 = K - 1
      NRIGHT = I1
      NLEFT = I2
!
! INITIALIZE VARIABLES
!
      NL = IEND(KM1)
      N1 = 1
      IF (NLEFT .NE. 1) N1 = IEND(NLEFT-1) + 1
      N2 = IEND(NRIGHT)
      NF = MAX0(N1,N2)
!
! INSERT K AS A NEIGHBOR OF MAX(NRIGHT,NLEFT)
!
      CALL SHIFTD(NF,NL,2, IADJ)
      IADJ(NF+1) = K
      IMIN = MAX0(NRIGHT,NLEFT)
      DO 1 I = IMIN,KM1
    1   IEND(I) = IEND(I) + 2
!
! INITIALIZE KEND AND INSERT K AS A NEIGHBOR OF
!   MIN(NRIGHT,NLEFT)
!
      KEND = NL + 3
      NL = NF - 1
      NF = MIN0(N1,N2)
      CALL SHIFTD(NF,NL,1, IADJ)
      IADJ(NF) = K
      IMAX = IMIN - 1
      IMIN = MIN0(NRIGHT,NLEFT)
      DO 2 I = IMIN,IMAX
    2   IEND(I) = IEND(I) + 1
!
! INSERT NRIGHT AS THE FIRST NEIGHBOR OF K
!
      IADJ(KEND) = NRIGHT
!
! INITIALIZE INDX FOR LOOP ON BOUNDARY NODES BETWEEN NRIGHT
!   AND NLEFT
!
      INDX = IEND(NRIGHT) - 2
    3 NEXT = IADJ(INDX)
      IF (NEXT .EQ. NLEFT) GO TO 4
!
! CONNECT NEXT AND K
!
      KEND = KEND + 1
      IADJ(KEND) = NEXT
      INDX = IEND(NEXT)
      IADJ(INDX) = K
      INDX = INDX - 1
      GO TO 3
!
! INSERT NLEFT AND 0 AS THE LAST NEIGHBORS OF K
!
    4 IADJ(KEND+1) = NLEFT
      KEND = KEND + 2
      IADJ(KEND) = 0
      IEND(K) = KEND
      RETURN
      END
      SUBROUTINE BNODES (N,IADJ,IEND, NB,NA,NT,NODES)
      INTEGER N, IADJ(6*(N-1)), IEND(N), NB, NA, NT, NODES(1)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A TRIANGULATION OF N POINTS ON THE UNIT SPHERE,
! THIS SUBROUTINE RETURNS A VECTOR CONTAINING THE INDICES
! (IF ANY) OF THE COUNTERCLOCKWISE-ORDERED SEQUENCE OF NODES
! ON THE BOUNDARY OF THE CONVEX HULL OF THE SET OF POINTS.
! THE BOUNDARY IS EMPTY IF THE POINTS DO NOT LIE IN A SINGLE
! HEMISPHERE.  THE NUMBERS OF BOUNDARY NODES, ARCS, AND
! TRIANGLES ARE ALSO RETURNED.
!
! INPUT PARAMETERS -     N - NUMBER OF NODES IN THE MESH.
!
!                     IADJ - SET OF ADJACENCY LISTS OF
!                            NODES IN THE MESH.
!
!                     IEND - POINTERS TO THE ENDS OF
!                            ADJACENCY LISTS IN IADJ FOR
!                            EACH NODE IN THE MESH.
!
!                    NODES - VECTOR OF LENGTH .GE. NB.
!                            (NB .LE. N).
!
!   IADJ AND IEND MAY BE CREATED BY TRMESH AND ARE NOT
! ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS -   NB - NUMBER OF BOUNDARY NODES.
!
!                    NA,NT - NUMBER OF ARCS AND TRIANGLES,
!                            RESPECTIVELY, IN THE MESH.
!
!                    NODES - VECTOR OF NB BOUNDARY NODE
!                            INDICES RANGING FROM 1 TO N.
!
! MODULES REFERENCED BY BNODES - NONE
!
!***********************************************************
!
      INTEGER NN, NST, INDL, K, N0, INDF
!
! LOCAL PARAMETERS -
!
! NN =   LOCAL COPY OF N
! NST =  FIRST ELEMENT OF NODES -- ARBITRARILY CHOSEN
! INDL = IADJ INDEX OF THE LAST NEIGHBOR OF NST
! K =    NODES INDEX
! N0 =   BOUNDARY NODE TO BE ADDED TO NODES
! INDF = IADJ INDEX OF THE FIRST NEIGHBOR OF N0
!
      NN = N
!
! SEARCH FOR A BOUNDARY NODE
!
      DO 1 NST = 1,NN
        INDL = IEND(NST)
        IF (IADJ(INDL) .EQ. 0) GO TO 2
    1   CONTINUE
!
! NO BOUNDARY NODE EXISTS
!
      NB = 0
      NT = 2*(NN-2)
      NA = 3*(NN-2)
      RETURN
!
! NST IS THE FIRST BOUNDARY NODE ENCOUNTERED.  INITIALIZE
!   FOR BOUNDARY TRAVERSAL.
!
    2 NODES(1) = NST
      K = 1
      N0 = NST
!
! TRAVERSE THE BOUNDARY IN COUNTERCLOCKWISE ORDER
!
    3 INDF = 1
      IF (N0 .GT. 1) INDF = IEND(N0-1) + 1
      N0 = IADJ(INDF)
      IF (N0 .EQ. NST) GO TO 4
      K = K + 1
      NODES(K) = N0
      GO TO 3
!
! TERMINATION
!
    4 NB = K
      NT = 2*(NN-1) - NB
      NA = 3*(NN-1) - NB
      RETURN
      END
      SUBROUTINE CIRCLE (N, X,Y)
      INTEGER N
      DOUBLE PRECISION    X(N), Y(N)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE COMPUTES THE COORDINATES OF A SEQUENCE
! OF EQUISPACED POINTS ON A UNIT CIRCLE.  AN (N-1)-SIDED
! POLYGONAL APPROXIMATION TO THE CIRCLE MAY BE PLOTTED BY
! CONNECTING (X(I),Y(I)) TO (X(I+1),Y(I+1)) FOR I = 1,...,
! N-1.
!
! INPUT PARAMETERS -   N - NUMBER OF POINTS.  N .GE. 2.
!
!                    X,Y - VECTORS OF LENGTH .GE. N.
!
! N IS NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS - X,Y - COORDINATES OF POINTS ON THE
!                           CIRCLE WHERE (X(N),Y(N)) =
!                           (X(1),Y(1)) = (1,0).
!
! MODULES REFERENCED BY CIRCLE - NONE
!
! INTRINSIC FUNCTIONS CALLED BY CIRCLE - ATAN, FLOAT, COS,
!                                        SIN
!
!***********************************************************
!
      INTEGER NM1, I
      DOUBLE PRECISION    DTH, TH
!
! LOCAL PARAMETERS -
!
! NM1 = N - 1
! I =   DO-LOOP INDEX
! DTH = ANGLE BETWEEN ADJACENT POINTS
! TH =  POLAR COORDINATE ANGLE
!
      NM1 = N - 1
      IF (NM1 .LT. 1) RETURN
      DTH = 8.*ATAN(1.)/FLOAT(NM1)
      DO 1 I = 1,NM1
        TH = FLOAT(I-1)*DTH
        X(I) = COS(TH)
    1   Y(I) = SIN(TH)
      X(N) = X(1)
      Y(N) = Y(1)
      RETURN
      END
      SUBROUTINE CONSTR (XK,YK,ZK, CX,SX,CY,SY)
      DOUBLE PRECISION XK, YK, ZK, CX, SX, CY, SY
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE CONSTRUCTS THE ELEMENTS OF A 3 BY 3
! ORTHOGONAL MATRIX R WHICH ROTATES A POINT (XK,YK,ZK) ON
! THE UNIT SPHERE TO THE NORTH POLE, I.E.
!
!      (XK)     (CY  0 -SY)   (1   0   0)   (XK)     (0)
!  R * (YK)  =  ( 0  1   0) * (0  CX -SX) * (YK)  =  (0)
!      (ZK)     (SY  0  CY)   (0  SX  CX)   (ZK)     (1)
!
! INPUT PARAMETERS - XK,YK,ZK - COMPONENTS OF A UNIT VECTOR
!                               TO BE ROTATED TO (0,0,1).
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS - CX,SX,CY,SY - ELEMENTS OF R -- CX,SX
!                                   DEFINE A ROTATION ABOUT
!                                   THE X-AXIS AND CY,SY DE-
!                                   FINE A ROTATION ABOUT
!                                   THE Y-AXIS.
!
! MODULES REFERENCED BY CONSTR - NONE
!
! INTRINSIC FUNCTION CALLED BY CONSTR - SQRT
!
!***********************************************************
!
      CY = SQRT(YK*YK + ZK*ZK)
      SY = XK
      IF (CY .EQ. 0.) GO TO 1
      CX = ZK/CY
      SX = YK/CY
      RETURN
!
! (XK,YK,ZK) LIES ON THE X-AXIS
!
    1 CX = 1.
      SX = 0.
      RETURN
      END
      SUBROUTINE COVSPH (KK,NODE, IADJ,IEND )
      INTEGER KK, NODE, IADJ(1), IEND(KK)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE CONNECTS AN EXTERIOR NODE KK TO ALL
! BOUNDARY NODES OF A TRIANGULATION OF KK-1 POINTS ON THE
! UNIT SPHERE, PRODUCING A TRIANGULATION WHICH COVERS THE
! SPHERE.  IADJ AND IEND ARE UPDATED WITH THE ADDITION OF
! NODE KK, BUT NO OPTIMIZATION OF THE MESH IS PERFORMED.
! ALL BOUNDARY NODES MUST BE VISIBLE FROM KK.
!
! INPUT PARAMETERS -   KK - INDEX OF THE EXTERIOR NODE TO
!                           BE ADDED.  KK .GE. 4.
!
!                    NODE - BOUNDARY NODE INDEX IN THE
!                           RANGE 1,...,KK-1.
!
!                    IADJ - SET OF ADJACENCY LISTS FOR
!                           NODES 1,...,KK-1.
!
!                    IEND - POINTERS TO THE ENDS OF ADJA-
!                           CENCY LISTS IN IADJ FOR NODES
!                           1,...,KK-1.
!
! KK AND NODE ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS - IADJ,IEND - UPDATED WITH THE ADDITION
!                                 OF NODE KK AS THE LAST
!                                 ENTRY.  ALL NODES ARE
!                                 INTERIOR.
!
! MODULES REFERENCED BY COVSPH - NONE
!
!***********************************************************
!
      INTEGER K, ND, NEXT, KEND, INDX
!
! LOCAL PARAMETERS -
!
! K,ND = LOCAL COPIES OF KK AND NODE
! NEXT = BOUNDARY NODE TO BE CONNECTED TO K
! KEND = IADJ INDEX OF THE LAST NEIGHBOR OF K
! INDX = IADJ INDEX
!
      K = KK
      ND = NODE
!
! INITIALIZATION
!
      NEXT = ND
      KEND = IEND(K-1)
!
! WALK ALONG THE BOUNDARY CONNECTING NODE K AND NEXT.  K
!   K REPLACES 0 AS THE LAST NEIGHBOR OF NEXT.
!
    1 KEND = KEND + 1
      IADJ(KEND) = NEXT
      INDX = IEND(NEXT)
      IADJ(INDX) = K
      NEXT = IADJ(INDX-1)
      IF (NEXT .NE. ND) GO TO 1
      IEND(K) = KEND
      RETURN
      END
      SUBROUTINE DELETE (N,NOUT1,NOUT2, IADJ,IEND, IER)
      INTEGER N, NOUT1, NOUT2, IADJ(6*(N-1)), IEND(N), IER
      EXTERNAL INDX
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE DELETES A BOUNDARY EDGE FROM A TRIANGU-
! LATION OF A SET OF POINTS ON THE UNIT SPHERE.  IT MAY BE
! NECESSARY TO FORCE CERTAIN EDGES TO BE PRESENT BEFORE
! CALLING DELETE (SEE SUBROUTINE EDGE).  NOTE THAT SUBROU-
! TINES EDGE, TRFIND, AND THE ROUTINES WHICH CALL TRFIND
! (ADNODE, UNIF, INTRC1, AND INTRC0) SHOULD NOT BE CALLED
! FOLLOWING A DELETION.
!
! INPUT PARAMETERS -     N - NUMBER OF NODES IN THE TRIAN-
!                            GULATION.
!
!              NOUT1,NOUT2 - PAIR OF ADJACENT NODES ON THE
!                            BOUNDARY DEFINING THE ARC TO
!                            BE REMOVED.  NOUT2 MUST BE THE
!                            LAST NONZERO NEIGHBOR OF NOUT1.
!
! THE ABOVE PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
!                IADJ,IEND - DATA STRUCTURE DEFINING THE
!                            TRIANGULATION (SEE SUBROUTINE
!                            TRMESH).
!
! OUTPUT PARAMETERS - IADJ,IEND - UPDATED WITH THE REMOVAL
!                                 OF THE ARC NOUT1-NOUT2
!                                 IF IER .EQ. 0.
!
!                           IER - ERROR INDICATOR
!                                 IER = 0 IF NO ERRORS WERE
!                                         ENCOUNTERED.
!                                 IER = 1 IF NOUT1 OR NOUT2
!                                         IS NOT ON THE
!                                         BOUNDARY.
!                                 IER = 2 IF NOUT1 OR NOUT2
!                                         HAS ONLY 2 NONZERO
!                                         NEIGHBORS.
!                                 IER = 3 IF NOUT2 IS NOT
!                                         THE LAST NEIGHBOR
!                                         OF NOUT1.
!                                 IER = 4 IF A DELETION
!                                         WOULD DIVIDE THE
!                                         MESH INTO TWO
!                                         REGIONS.
!
! MODULES REFERENCED BY DELETE - SHIFTD, INDX
!
!***********************************************************
!
      INTEGER NN, IOUT1, IOUT2, IO1, IO2, IND12, IND21, & 
              ITEMP, IND1F, IND1L, IND2F, IND2L, NEWBD, &
              INDNF, INDNL, INDN0, INDFP2, INDLM3, NF, NL, &
              I, IMAX, INDX
!
! LOCAL PARAMETERS -
!
! NN =          LOCAL COPY OF N
! IOUT1,IOUT2 = LOCAL COPIES OF NOUT1 AND NOUT2
! IO1,IO2 =     NOUT1,NOUT2 IN ORDER OF INCREASING MAGNITUDE
! IND12 =       INDEX OF IO2 IN THE ADJACENCY LIST FOR IO1
! IND21 =       INDEX OF IO1 IN THE ADJACENCY LIST FOR IO2
! ITEMP =       TEMPORARY STORAGE LOCATION FOR PERMUTATIONS
! IND1F =       IADJ INDEX OF THE FIRST NEIGHBOR OF IO1
! IND1L =       IADJ INDEX OF THE LAST NEIGHBOR OF IO1
! IND2F =       IADJ INDEX OF THE FIRST NEIGHBOR OF IO2
! IND2L =       IADJ INDEX OF THE LAST NEIGHBOR OF IO2
! NEWBD =       THE NEIGHBOR COMMON TO NOUT1 AND NOUT2
! INDNF =       IADJ INDEX OF THE FIRST NEIGHBOR OF NEWBD
! INDNL =       IADJ INDEX OF THE LAST NEIGHBOR OF NEWBD
! INDN0 =       INDEX OF 0 IN THE ADJACENCY LIST FOR NEWBD
!                 BEFORE PERMUTING THE NEIGHBORS
! INDFP2 =      INDNF + 2
! INDLM3 =      INDNL - 3
! NF,NL =       BOUNDS ON THE PORTION OF IADJ TO BE SHIFTED
! I =           DO-LOOP INDEX
! IMAX =        UPPER BOUND ON DO-LOOP FOR SHIFTING IEND
!
      NN = N
      IOUT1 = NOUT1
      IOUT2 = NOUT2
!
! INITIALIZE INDICES
!
      IND1F = 1
      IF (IOUT1 .GT. 1) IND1F = IEND(IOUT1-1) + 1
      IND1L = IEND(IOUT1)
      IND2F = 1
      IF (IOUT2 .GT. 1) IND2F = IEND(IOUT2-1) + 1
      IND2L = IEND(IOUT2)
      NEWBD = IADJ(IND1L-2)
      INDN0 = INDX(NEWBD,IOUT2,IADJ,IEND)
      INDNL = IEND(NEWBD)
!
! ORDER VERTICES SUCH THAT THE NEIGHBORS OF IO1 PRECEDE
!   THOSE OF IO2
!
      IF (IOUT1 .GT. IOUT2) GO TO 1
      IO1 = IOUT1
      IO2 = IOUT2
      IND12 = IND1L - 1
      IND21 = IND2F
      GO TO 2
    1 IO1 = IOUT2
      IO2 = IOUT1
      IND12 = IND2F
      IND21 = IND1L - 1
!
! CHECK FOR ERRORS
!
    2 IF ( (IADJ(IND1L) .NE. 0) .OR. (IADJ(IND2L) .NE. 0) ) &
         GO TO 21
      IF ( (IND1L-IND1F .LE. 2) .OR. (IND2L-IND2F .LE. 2) ) &
         GO TO 22
      IF (IADJ(IND1L-1) .NE. IOUT2) GO TO 23
      IF (IADJ(INDNL) .EQ. 0) GO TO 24
!
! DELETE THE EDGE IO1-IO2 AND MAKE NEWBD A BOUNDARY NODE
!
      IF (NEWBD .LT. IO1) GO TO 8
      IF (NEWBD .LT. IO2) GO TO 6
!
! THE VERTICES ARE ORDERED IO1, IO2, NEWBD.
! DELETE IO2 AS A NEIGHBOR OF IO1.
!
      NF = IND12 + 1
      NL = IND21 - 1
      CALL SHIFTD(NF,NL,-1, IADJ)
      IMAX = IO2 - 1
      DO 3 I = IO1,IMAX
    3   IEND(I) = IEND(I) - 1
!
! DELETE IO1 AS A NEIGHBOR OF IO2
!
      NF = NL + 2
      NL = INDN0
      CALL SHIFTD(NF,NL,-2, IADJ)
      IMAX = NEWBD - 1
      DO 4 I = IO2,IMAX
    4   IEND(I) = IEND(I) - 2
!
! SHIFT THE BOTTOM OF IADJ UP 1 LEAVING ROOM FOR 0 AS A
!   NEIGHBOR OF NEWBD
!
      INDN0 = INDN0 - 1
      NF = NL + 1
      NL = IEND(NN)
      IF (NF .LE. NL) CALL SHIFTD(NF,NL,-1, IADJ)
      DO 5 I = NEWBD,NN
    5   IEND(I) = IEND(I) - 1
      GO TO 12
!
! THE VERTICES ARE ORDERED IO1, NEWBD, IO2.
! DELETE IO2 AS A NEIGHBOR OF IO1 LEAVING ROOM FOR 0 AS A
!   NEIGHBOR OF NEWBD.
!
    6 NF = IND12 + 1
      NL = INDN0
      CALL SHIFTD(NF,NL,-1, IADJ)
      IMAX = NEWBD - 1
      DO 7 I = IO1,IMAX
    7   IEND(I) = IEND(I) - 1
      GO TO 10
!
! THE VERTICES ARE ORDERED NEWBD, IO1, IO2.
! DELETE IO2 AS A NEIGHBOR OF IO1 LEAVING ROOM FOR 0 AS A
!   NEIGHBOR OF NEWBD.
!
    8 INDN0 = INDN0 + 1
      NF = INDN0
      NL = IND12 - 1
      IF (NF .LE. NL) CALL SHIFTD(NF,NL,1, IADJ)
      IMAX = IO1 - 1
      DO 9 I = NEWBD,IMAX
    9   IEND(I) = IEND(I) + 1
!
! DELETE IO1 AS A NEIGHBOR OF IO2
!
   10 NF = IND21 + 1
      NL = IEND(NN)
      CALL SHIFTD(NF,NL,-1, IADJ)
      DO 11 I = IO2,NN
   11   IEND(I) = IEND(I) - 1
!
! PERMUTE THE NEIGHBORS OF NEWBD WITH END-AROUND SHIFTS SO
!   THAT 0 IS THE LAST NEIGHBOR
!
   12 INDNF = 1
      IF (NEWBD .GT. 1) INDNF = IEND(NEWBD-1) + 1
      INDNL = IEND(NEWBD)
      IF (INDN0-INDNF .GE. INDNL-INDN0) GO TO 16
!
! SHIFT UPWARD
!
      IF (INDN0 .GT. INDNF) GO TO 13
      CALL SHIFTD(INDNF+1,INDNL,-1, IADJ)
      GO TO 20
   13 INDFP2 = INDNF + 2
      IF (INDN0 .LT. INDFP2) GO TO 15
      DO 14 I = INDFP2,INDN0
        ITEMP = IADJ(INDNF)
        CALL SHIFTD(INDNF+1,INDNL,-1, IADJ)
   14   IADJ(INDNL) = ITEMP
!
! THE LAST SHIFT IS BY 2
!
   15 ITEMP = IADJ(INDNF)
      CALL SHIFTD(INDFP2,INDNL,-2, IADJ)
      IADJ(INDNL-1) = ITEMP
      GO TO 20
!
! SHIFT DOWNWARD
!
   16 IF (INDN0 .EQ. INDNL) GO TO 20
      IF (INDN0 .LT. INDNL-1) GO TO 17
      CALL SHIFTD(INDNF,INDNL-2,1, IADJ)
      IADJ(INDNF) = IADJ(INDNL)
      GO TO 20
   17 INDLM3 = INDNL - 3
      IF (INDN0 .GT. INDLM3) GO TO 19
      DO 18 I = INDN0,INDLM3
        ITEMP = IADJ(INDNL)
        CALL SHIFTD(INDNF,INDNL-1,1, IADJ)
   18   IADJ(INDNF) = ITEMP
!
! THE LAST SHIFT IS BY 2
!
   19 ITEMP = IADJ(INDNL-1)
      CALL SHIFTD(INDNF,INDLM3,2, IADJ)
      IADJ(INDNF+1) = IADJ(INDNL)
      IADJ(INDNF) = ITEMP
!
! INSERT 0 AS THE LAST NEIGHBOR OF NEWBD
!
   20 IADJ(INDNL) = 0
      IER = 0
      RETURN
!
! ONE OF THE VERTICES IS NOT ON THE BOUNDARY
!
   21 IER = 1
      RETURN
!
! ONE OF THE VERTICES HAS ONLY TWO NONZERO NEIGHBORS -- THE
!   TRIANGULATION WOULD BE DESTROYED BY A DELETION
!
   22 IER = 2
      RETURN
!
! NOUT2 IS NOT THE LAST NONZERO NEIGHBOR OF NOUT1
!
   23 IER = 3
      RETURN
!
! A DELETION WOULD DIVIDE THE MESH INTO TWO REGIONS
!   CONNECTED AT A SINGLE NODE
!
   24 IER = 4
      RETURN
      END
      SUBROUTINE EDGE (IN1,IN2,X,Y,Z, LWK,IWK,IADJ, &
                       IEND, IER)
      LOGICAL SWPTST
      INTEGER IN1, IN2, LWK, IWK(2,1), IADJ(1), IEND(1), IER
      DOUBLE PRECISION    X(1), Y(1), Z(1)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A TRIANGULATION OF N NODES AND A PAIR OF NODAL
! INDICES IN1 AND IN2, THIS ROUTINE SWAPS ARCS AS NECESSARY
! TO FORCE IN1 AND IN2 TO BE ADJACENT.  ONLY ARCS WHICH
! INTERSECT IN1-IN2 ARE SWAPPED OUT.  IF A THIESSEN TRIANGU-
! LATION IS INPUT, THE RESULTING TRIANGULATION IS AS CLOSE
! AS POSSIBLE TO A THIESSEN TRIANGULATION IN THE SENSE THAT
! ALL ARCS OTHER THAN IN1-IN2 ARE LOCALLY OPTIMAL.
!   A SEQUENCE OF CALLS TO EDGE MAY BE USED TO FORCE THE
! PRESENCE OF A SET OF EDGES DEFINING THE BOUNDARY OF A NON-
! CONVEX REGION.  SUBSEQUENT DELETION OF EDGES OUTSIDE THIS
! REGION (BY SUBROUTINE DELETE) RESULTS IN A NONCONVEX TRI-
! ANGULATION WHICH MAY SERVE AS A FINITE ELEMENT GRID.
! (EDGE SHOULD NOT BE CALLED AFTER A CALL TO DELETE.)  IF,
! ON THE OTHER HAND, INTERPOLATION IS TO BE PERFORMED IN THE
! NONCONVEX REGION, EDGES MUST NOT BE DELETED, BUT IT IS
! STILL ADVANTAGEOUS TO HAVE THE NONCONVEX BOUNDARY PRESENT
! IF IT IS DESIRABLE THAT INTERPOLATED VALUES BE INFLUENCED
! BY THE GEOMETRY.  NOTE THAT SUBROUTINE GETNP WHICH IS USED
! TO SELECT THE NODES ENTERING INTO LOCAL DERIVATIVE ESTI-
! MATES WILL NOT NECESSARILY RETURN CLOSEST NODES IF THE
! TRIANGULATION HAS BEEN RENDERED NONOPTIMAL BY A CALL TO
! EDGE.  HOWEVER, THE EFFECT WILL BE MERELY TO FURTHER EN-
! HANCE THE INFLUENCE OF THE NONCONVEX GEOMETRY ON INTERPO-
! LATED VALUES.
!
! INPUT PARAMETERS - IN1,IN2 - INDICES (OF X,Y,Z) IN THE
!                              RANGE 1,...,N DEFINING A PAIR
!                              OF NODES TO BE CONNECTED BY
!                              AN ARC.
!
!                      X,Y,Z - N-VECTORS CONTAINING CARTE-
!                              SIAN COORDINATES OF THE
!                              NODES.
!
! THE ABOVE PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
!                        LWK - NUMBER OF COLUMNS RESERVED
!                              FOR IWK.  THIS MUST BE AT
!                              LEAST NI -- THE NUMBER OF
!                              ARCS WHICH INTERSECT IN1-IN2.
!                              (NI IS BOUNDED BY N-3).
!
!                        IWK - INTEGER WORK ARRAY DIMENSION-
!                              ED 2 BY LWK (OR VECTOR OF
!                              LENGTH .GE. 2*LWK).
!
!                  IADJ,IEND - DATA STRUCTURE DEFINING THE
!                              TRIANGULATION.  SEE SUBROU-
!                              TINE TRMESH.
!
! OUTPUT PARAMETERS - LWK - NUMBER OF IWK COLUMNS REQUIRED
!                           IF IER = 0 OR IER = 2.  LWK = 0
!                           IFF IN1 AND IN2 WERE ADJACENT
!                           ON INPUT.
!
!                     IWK - CONTAINS THE INDICES OF THE END-
!                           POINTS OF THE NEW ARCS OTHER
!                           THAN IN1-IN2 UNLESS IER .GT. 0
!                           OR LWK = 0.  NEW ARCS TO THE
!                           LEFT OF IN1->IN2 ARE STORED IN
!                           THE FIRST K-1 COLUMNS (LEFT POR-
!                           TION OF IWK), COLUMN K CONTAINS
!                           ZEROS, AND NEW ARCS TO THE RIGHT
!                           OF IN1->IN2 OCCUPY COLUMNS K+1,
!                           ...,LWK.  (K CAN BE DETERMINED
!                           BY SEARCHING IWK FOR THE ZEROS.)
!
!               IADJ,IEND - UPDATED IF NECESSARY TO REFLECT
!                           THE PRESENCE OF AN ARC CONNECT-
!                           ING IN1 AND IN2, UNALTERED IF
!                           IER .NE. 0.
!
!                     IER - ERROR INDICATOR
!                           IER = 0 IF NO ERRORS WERE EN-
!                                   COUNTERED.
!                           IER = 1 IF IN1 .LT. 1, IN2 .LT.
!                                   1, IN1 = IN2, OR LWK
!                                   .LT. 0 ON INPUT.
!                           IER = 2 IF MORE SPACE IS REQUIR-
!                                   ED IN IWK.  SEE LWK.
!                           IER = 3 IF IN1 AND IN2 COULD NOT
!                                   BE CONNECTED DUE TO AN
!                                   INVALID DATA STRUCTURE.
!
! MODULES REFERENCED BY EDGE - SWAP, INDEX, SHIFTD, SWPTST
!
!***********************************************************
!
      INTEGER N1, N2, IWEND, IWL, INDF, INDX, N1LST, NL, NR, &
              NEXT, IWF, LFT, N0, IWC, IWCP1, IWCM1, I, IO1, &
              IO2, INDL
      DOUBLE PRECISION    X1, Y1, Z1, X2, Y2, Z2, X0, Y0, Z0
      DOUBLE PRECISION    XA,YA,ZA,XB,YB,ZB,XP,YP,ZP
      LOGICAL SWP, LEFT
!
! LOCAL PARAMETERS -
!
! N1,N2 =    LOCAL COPIES OF IN1 AND IN2 OR NODES OPPOSITE
!              AN ARC IO1-IO2 TO BE TESTED FOR A SWAP IN
!              THE OPTIMIZATION LOOPS
! IWEND =    INPUT OR OUTPUT VALUE OF LWK
! IWL =      IWK (COLUMN) INDEX OF THE LAST (RIGHTMOST) ARC
!              WHICH INTERSECTS IN1->IN2
! INDF =     IADJ INDEX OF THE FIRST NEIGHBOR OF IN1 OR IO1
! INDX =     IADJ INDEX OF A NEIGHBOR OF IN1, NL, OR IO1
! N1LST =    LAST NEIGHBOR OF IN1
! NL,NR =    ENDPOINTS OF AN ARC WHICH INTERSECTS IN1-IN2
!              WITH NL LEFT IN1->IN2
! NEXT =     NODE OPPOSITE NL->NR
! IWF =      IWK (COLUMN) INDEX OF THE FIRST (LEFTMOST) ARC
!              WHICH INTERSECTS IN1->IN2
! LFT =      FLAG USED TO DETERMINE IF A SWAP RESULTS IN THE
!              NEW ARC INTERSECTING IN1-IN2 -- LFT = 0 IFF
!              N0 = IN1, LFT = -1 IMPLIES N0 LEFT IN1->IN2,
!              AND LFT = 1 IMPLIES N0 LEFT IN2->IN1
! N0 =       NODE OPPOSITE NR->NL
! IWC =      IWK INDEX BETWEEN IWF AND IWL -- NL->NR IS
!              STORED IN IWK(1,IWC)->IWK(2,IWC)
! IWCP1 =    IWC + 1
! IWCM1 =    IWC - 1
! I =        DO-LOOP INDEX AND COLUMN INDEX FOR IWK
! IO1,IO2 =  ENDPOINTS OF AN ARC TO BE TESTED FOR A SWAP IN
!              THE OPTIMIZATION LOOPS
! INDL =     IADJ INDEX OF THE LAST NEIGHBOR OF IO1
! X1,Y1,Z1 = COORDINATES OF IN1
! X2,Y2,Z2 = COORDINATES OF IN2
! X0,Y0,Z0 = COORDINATES OF N0
! SWP =      FLAG SET TO .TRUE. IFF A SWAP OCCURS IN AN OPT-
!              IMIZATION LOOP
! LEFT =     STATEMENT FUNCTION WHICH RETURNS THE VALUE
!              .TRUE. IFF (XP,YP,ZP) IS ON OR TO THE LEFT OF
!              THE VECTOR (XA,YA,ZA)->(XB,YB,ZB)
!
      LEFT(XA,YA,ZA,XB,YB,ZB,XP,YP,ZP) = XP*(YA*ZB-YB*ZA) &
            - YP*(XA*ZB-XB*ZA) + ZP*(XA*YB-XB*YA) .GE. 0.
!
! STORE IN1, IN2, AND LWK IN LOCAL VARIABLES AND CHECK FOR
!   ERRORS.
!
      N1 = IN1
      N2 = IN2
      IWEND = LWK
      IF (N1 .LT. 1  .OR.  N2 .LT. 1  .OR.  N1 .EQ. N2  .OR. &
          IWEND .LT. 0) GO TO 35
!
! STORE THE COORDINATES OF N1 AND N2 AND INITIALIZE IWL.
!
      X1 = X(N1)
      Y1 = Y(N1)
      Z1 = Z(N1)
      X2 = X(N2)
      Y2 = Y(N2)
      Z2 = Z(N2)
      IWL = 0
!
! SET NR AND NL TO ADJACENT NEIGHBORS OF N1 SUCH THAT
!   NR LEFT N2->N1 AND NL LEFT N1->N2.
!
!   SET INDF AND INDX TO THE INDICES OF THE FIRST AND LAST
!     NEIGHBORS OF N1 AND SET N1LST TO THE LAST NEIGHBOR.
!
      INDF = 1
      IF (N1 .GT. 1) INDF = IEND(N1-1) + 1
      INDX = IEND(N1)
      N1LST = IADJ(INDX)
      IF (N1LST .EQ. 0) INDX = INDX - 1
      IF (N1LST .EQ. 0) GO TO 2
!
!   N1 IS AN INTERIOR NODE.  LOOP THROUGH THE NEIGHBORS NL
!     IN REVERSE ORDER UNTIL NL LEFT N1->N2.
!
      NL = N1LST
    1 IF ( LEFT(X1,Y1,Z1,X2,Y2,Z2,X(NL),Y(NL),Z(NL)) ) &
           GO TO 2
      INDX = INDX - 1
      NL = IADJ(INDX)
      IF (INDX .GT. INDF) GO TO 1
!
!   NL IS THE FIRST NEIGHBOR OF N1.  SET NR TO THE LAST
!     NEIGHBOR AND TEST FOR AN ARC N1-N2.
!
      NR = N1LST
      IF (NL .EQ. N2) GO TO 34
      GO TO 4
!
!   NL = IADJ(INDX) LEFT N1->N2 AND INDX .GT. INDF.  SET
!     NR TO THE PRECEDING NEIGHBOR OF N1.
!
    2 INDX = INDX - 1
      NR = IADJ(INDX)
      IF ( LEFT(X2,Y2,Z2,X1,Y1,Z1,X(NR),Y(NR),Z(NR)) ) &
           GO TO 3
      IF (INDX .GT. INDF) GO TO 2
!
!   SET NL AND NR TO THE FIRST AND LAST NEIGHBORS OF N1 AND
!     TEST FOR AN INVALID DATA STRUCTURE (N1 CANNOT BE A
!     BOUNDARY NODE AND CANNOT BE ADJACENT TO N2).
!
      NL = NR
      NR = N1LST
      IF (NR .EQ. 0  .OR.  NR .EQ. N2) GO TO 37
      GO TO 4
!
!   SET NL TO THE NEIGHBOR FOLLOWING NR AND TEST FOR AN ARC
!     N1-N2.
!
    3 NL = IADJ(INDX+1)
      IF (NL .EQ. N2  .OR.  NR .EQ. N2) GO TO 34
!
! STORE THE ORDERED SEQUENCE OF INTERSECTING EDGES NL->NR IN
!   IWK(1,IWL)->IWK(2,IWL).
!
    4 IWL = IWL + 1
      IF (IWL .LE. IWEND) IWK(1,IWL) = NL
      IF (IWL .LE. IWEND) IWK(2,IWL) = NR
!
!   SET NEXT TO THE NEIGHBOR OF NL WHICH FOLLOWS NR.
!
      INDX = IEND(NL)
      IF (IADJ(INDX) .NE. NR) GO TO 5
!
!   NR IS THE LAST NEIGHBOR OF NL.  SET NEXT TO THE FIRST
!     NEIGHBOR.
!
      INDX = 0
      IF (NL .NE. 1) INDX = IEND(NL-1)
      GO TO 6
!
!   NR IS NOT THE LAST NEIGHBOR OF NL.  LOOP THROUGH THE
!     NEIGHBORS IN REVERSE ORDER.
!
    5 INDX = INDX - 1
      IF (IADJ(INDX) .NE. NR) GO TO 5
!
!   STORE NEXT, TEST FOR AN INVALID TRIANGULATION (NL->NR
!     CANNOT BE A BOUNDARY EDGE), AND TEST FOR TERMINATION
!     OF THE LOOP.
!
    6 NEXT = IADJ(INDX+1)
      IF (NEXT .EQ. 0) GO TO 37
      IF (NEXT .EQ. N2) GO TO 8
!
!   SET NL OR NR TO NEXT.
!
      IF ( LEFT(X1,Y1,Z1,X2,Y2,Z2,X(NEXT),Y(NEXT),Z(NEXT)) ) &
           GO TO 7
      NR = NEXT
      GO TO 4
    7 NL = NEXT
      GO TO 4
!
! IWL IS THE NUMBER OF ARCS WHICH INTERSECT N1-N2.  STORE
!   LWK AND TEST FOR SUFFICIENT SPACE.
!
    8 LWK = IWL
      IF (IWL .GT. IWEND) GO TO 36
      IWEND = IWL
!
! INITIALIZE FOR EDGE SWAPPING LOOP -- ALL POSSIBLE SWAPS
!   ARE APPLIED (EVEN IF THE NEW ARC AGAIN INTERSECTS
!   N1-N2), ARCS TO THE LEFT OF N1->N2 ARE STORED IN THE
!   LEFT PORTION OF IWK, AND ARCS TO THE RIGHT ARE STORED IN
!   THE RIGHT PORTION.  IWF AND IWL INDEX THE FIRST AND LAST
!   INTERSECTING ARCS.
!
      IER = 0
      IWF = 1
!
! TOP OF LOOP -- SET N0 TO N1 AND NL->NR TO THE FIRST EDGE.
!   IWC POINTS TO THE ARC CURRENTLY BEING PROCESSED.  LFT
!   .LE. 0 IFF N0 LEFT N1->N2.
!
    9 LFT = 0
      N0 = N1
      X0 = X1
      Y0 = Y1
      Z0 = Z1
      NL = IWK(1,IWF)
      NR = IWK(2,IWF)
      IWC = IWF
!
!   SET NEXT TO THE NODE OPPOSITE NL->NR UNLESS IWC IS THE
!     LAST ARC.
!
   10 IF (IWC .EQ. IWL) GO TO 21
      IWCP1 = IWC + 1
      NEXT = IWK(1,IWCP1)
      IF (NEXT .NE. NL) GO TO 15
      NEXT = IWK(2,IWCP1)
!
!   NEXT RIGHT N1->N2 AND IWC .LT. IWL.  TEST FOR A POSSIBLE
!     SWAP.
!
      IF ( .NOT. LEFT(X0,Y0,Z0,X(NR),Y(NR),Z(NR),X(NEXT), &
           Y(NEXT),Z(NEXT)) ) GO TO 13
      IF (LFT .GE. 0) GO TO 11 
      IF ( .NOT. LEFT(X(NL),Y(NL),Z(NL),X0,Y0,Z0,X(NEXT), &
           Y(NEXT),Z(NEXT)) ) GO TO 13
!
!   REPLACE NL->NR WITH N0->NEXT.
!
      CALL SWAP(NEXT,N0,NL,NR, IADJ,IEND )
      IWK(1,IWC) = N0
      IWK(2,IWC) = NEXT
      GO TO 14
!
!   SWAP NL-NR FOR N0-NEXT, SHIFT COLUMNS IWC+1,...,IWL TO
!     THE LEFT, AND STORE N0-NEXT IN THE RIGHT PORTION OF
!     IWK.
!
   11 CALL SWAP(NEXT,N0,NL,NR, IADJ,IEND )
      DO 12 I = IWCP1,IWL
        IWK(1,I-1) = IWK(1,I)
   12   IWK(2,I-1) = IWK(2,I)
      IWK(1,IWL) = N0
      IWK(2,IWL) = NEXT
      IWL = IWL - 1
      NR = NEXT
      GO TO 10
!
!   A SWAP IS NOT POSSIBLE.  SET N0 TO NR.
!
   13 N0 = NR
      X0 = X(N0)
      Y0 = Y(N0)
      Z0 = Z(N0)
      LFT = 1
!
!   ADVANCE TO THE NEXT ARC.
!
   14 NR = NEXT
      IWC = IWC + 1
      GO TO 10
!
!   NEXT LEFT N1->N2, NEXT .NE. N2, AND IWC .LT. IWL.
!     TEST FOR A POSSIBLE SWAP.
!
   15 IF ( .NOT. LEFT(X(NL),Y(NL),Z(NL),X0,Y0,Z0,X(NEXT), &
           Y(NEXT),Z(NEXT)) ) GO TO 19
      IF (LFT .LE. 0) GO TO 16
      IF ( .NOT. LEFT(X0,Y0,Z0,X(NR),Y(NR),Z(NR),X(NEXT), &
           Y(NEXT),Z(NEXT)) ) GO TO 19
!
!   REPLACE NL->NR WITH NEXT->N0.
!
      CALL SWAP(NEXT,N0,NL,NR, IADJ,IEND )
      IWK(1,IWC) = NEXT
      IWK(2,IWC) = N0
      GO TO 20
!
!   SWAP NL-NR FOR N0-NEXT, SHIFT COLUMNS IWF,...,IWC-1 TO
!     THE RIGHT, AND STORE N0-NEXT IN THE LEFT PORTION OF
!     IWK.
!
   16 CALL SWAP(NEXT,N0,NL,NR, IADJ,IEND )
      I = IWC
   17 IF (I .EQ. IWF) GO TO 18
      IWK(1,I) = IWK(1,I-1)
      IWK(2,I) = IWK(2,I-1)
      I = I - 1
      GO TO 17
   18 IWK(1,IWF) = N0
      IWK(2,IWF) = NEXT
      IWF = IWF + 1
      GO TO 20
!
!   A SWAP IS NOT POSSIBLE.  SET N0 TO NL.
!
   19 N0 = NL
      X0 = X(N0)
      Y0 = Y(N0)
      Z0 = Z(N0)
      LFT = -1
!
!   ADVANCE TO THE NEXT ARC.
!
   20 NL = NEXT
      IWC = IWC + 1
      GO TO 10
!
!   N2 IS OPPOSITE NL->NR (IWC = IWL).
!
   21 IF (N0 .EQ. N1) GO TO 24
      IF (LFT .LT. 0) GO TO 22
!
!   N0 RIGHT N1->N2.  TEST FOR A POSSIBLE SWAP.
!
      IF ( .NOT. LEFT(X0,Y0,Z0,X(NR),Y(NR),Z(NR),X2,Y2,Z2) ) &
           GO TO 9
!
!   SWAP NL-NR FOR N0-N2 AND STORE N0-N2 IN THE RIGHT
!     PORTION OF IWK.
!
      CALL SWAP(N2,N0,NL,NR, IADJ,IEND )
      IWK(1,IWL) = N0
      IWK(2,IWL) = N2
      IWL = IWL - 1
      GO TO 9
!
!   N0 LEFT N1->N2.  TEST FOR A POSSIBLE SWAP.
!
   22 IF ( .NOT. LEFT(X(NL),Y(NL),Z(NL),X0,Y0,Z0,X2,Y2,Z2) ) &
           GO TO 9
!
!   SWAP NL-NR FOR N0-N2, SHIFT COLUMNS IWF,...,IWL-1 TO THE
!     RIGHT, AND STORE N0-N2 IN THE LEFT PORTION OF IWK.
!
      CALL SWAP(N2,N0,NL,NR, IADJ,IEND )
      I = IWL
   23 IWK(1,I) = IWK(1,I-1)
      IWK(2,I) = IWK(2,I-1)
      I = I - 1
      IF (I .GT. IWF) GO TO 23
      IWK(1,IWF) = N0
      IWK(2,IWF) = N2
      IWF = IWF + 1
      GO TO 9
!
! IWF = IWC = IWL.  SWAP OUT THE LAST ARC FOR N1-N2 AND
!   STORE ZEROS IN IWK.
!
   24 CALL SWAP(N2,N1,NL,NR, IADJ,IEND )
      IWK(1,IWC) = 0
      IWK(2,IWC) = 0
      IF (IWC .EQ. 1) GO TO 29
!
! OPTIMIZATION LOOPS -- OPTIMIZE THE SET OF NEW ARCS TO THE
!   LEFT OF IN1->IN2.  THE LOOP IS REPEATED UNTIL NO SWAPS
!   ARE PERFORMED.
!
      IWCM1 = IWC - 1
   25 SWP = .FALSE.
      DO 28 I = 1,IWCM1
        IO1 = IWK(1,I)
        IO2 = IWK(2,I)
!
!   SET N1 TO THE NEIGHBOR OF IO1 WHICH FOLLOWS IO2 AND SET
!     N2 TO THE NEIGHBOR OF IO1 WHICH PRECEDES IO2.
!
        INDF = 1
        IF (IO1 .GT. 1) INDF = IEND(IO1-1) + 1
        INDL = IEND(IO1)
        INDX = INDL
        IF (IADJ(INDX) .NE. IO2) GO TO 26
!
!   IO2 IS THE LAST NEIGHBOR OF IO1.
!
        N1 = IADJ(INDF)
        N2 = IADJ(INDX-1)
        GO TO 27
!
!   IO2 IS NOT THE LAST NEIGHBOR OF IO1.  LOOP THROUGH THE
!     NEIGHBORS IN REVERSE ORDER.
!
   26   INDX = INDX - 1
        IF (IADJ(INDX) .NE. IO2) GO TO 26
        N1 = IADJ(INDX+1)
        IF (INDX .NE. INDF) N2 = IADJ(INDX-1)
        IF (INDX .EQ. INDF) N2 = IADJ(INDL)
!
!   TEST IO1-IO2 FOR A SWAP.
!
   27   IF ( .NOT. SWPTST(N1,N2,IO1,IO2,X,Y,Z) ) GO TO 28
        SWP = .TRUE.
        CALL SWAP(N1,N2,IO1,IO2, IADJ,IEND )
        IWK(1,I) = N1
        IWK(2,I) = N2
   28   CONTINUE
      IF (SWP) GO TO 25
!
! TEST FOR TERMINATION.
!
   29 IF (IWC .EQ. IWEND) RETURN
      IWCP1 = IWC + 1
!
! OPTIMIZE THE SET OF NEW ARCS TO THE RIGHT OF IN1->IN2.
!
   30 SWP = .FALSE.
      DO 33 I = IWCP1,IWEND
        IO1 = IWK(1,I)
        IO2 = IWK(2,I)
!
!   SET N1 AND N2 TO THE NODES OPPOSITE IO1->IO2 AND
!     IO2->IO1, RESPECTIVELY.
!
        INDF = 1
        IF (IO1 .GT. 1) INDF = IEND(IO1-1) + 1
        INDL = IEND(IO1)
        INDX = INDL
        IF (IADJ(INDX) .NE. IO2) GO TO 31
!
        N1 = IADJ(INDF)
        N2 = IADJ(INDX-1)
        GO TO 32
!
   31   INDX = INDX - 1
        IF (IADJ(INDX) .NE. IO2) GO TO 31
        N1 = IADJ(INDX+1)
        IF (INDX .NE. INDF) N2 = IADJ(INDX-1)
        IF (INDX .EQ. INDF) N2 = IADJ(INDL)
!
   32   IF ( .NOT. SWPTST(N1,N2,IO1,IO2,X,Y,Z) ) GO TO 33
        SWP = .TRUE.
        CALL SWAP(N1,N2,IO1,IO2, IADJ,IEND )
        IWK(1,I) = N1
        IWK(2,I) = N2
   33   CONTINUE
      IF (SWP) GO TO 30
      RETURN
!
! IN1 AND IN2 WERE ADJACENT ON INPUT.
!
   34 IER = 0
      LWK = 0
      RETURN
!
! PARAMETER OUT OF RANGE
!
   35 IER = 1
      RETURN
!
! INSUFFICIENT SPACE IN IWK
!
   36 IER = 2
      RETURN
!
! INVALID TRIANGULATION DATA STRUCTURE
!
   37 IER = 3
      RETURN
      END
      SUBROUTINE GETNP (X,Y,Z,IADJ,IEND,L, NPTS, DF,IER)
      INTEGER IADJ(1), IEND(1), L, NPTS(L), IER
      DOUBLE PRECISION    X(1), Y(1), Z(1), DF
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A THIESSEN TRIANGULATION OF N NODES ON THE UNIT
! SPHERE AND AN ARRAY NPTS CONTAINING THE INDICES OF L-1
! NODES ORDERED BY ANGULAR DISTANCE FROM NPTS(1), THIS SUB-
! ROUTINE SETS NPTS(L) TO THE INDEX OF THE NEXT NODE IN THE
! SEQUENCE -- THE NODE, OTHER THAN NPTS(1),...,NPTS(L-1),
! WHICH IS CLOSEST TO NPTS(1).  THUS, THE ORDERED SEQUENCE
! OF K CLOSEST NODES TO N1 (INCLUDING N1) MAY BE DETERMINED
! BY K-1 CALLS TO GETNP WITH NPTS(1) = N1 AND L = 2,3,...,K
! FOR K .GE. 2.
!   THE ALGORITHM USES THE FACT THAT, IN A THIESSEN TRIAN-
! GULATION, THE K-TH CLOSEST NODE TO A GIVEN NODE N1 IS A
! NEIGHBOR OF ONE OF THE K-1 CLOSEST NODES TO N1.
!
! INPUT PARAMETERS - X,Y,Z - VECTORS OF LENGTH N CONTAINING
!                            THE CARTESIAN COORDINATES OF
!                            THE NODES.
!
!                     IADJ - SET OF ADJACENCY LISTS OF NODES
!                            IN THE TRIANGULATION.
!
!                     IEND - POINTERS TO THE ENDS OF ADJA-
!                            CENCY LISTS FOR EACH NODE IN
!                            THE TRIANGULATION.
!
!                        L - NUMBER OF NODES IN THE SEQUENCE
!                            ON OUTPUT.  2 .LE. L .LE. N.
!
!                     NPTS - ARRAY OF LENGTH .GE. L CONTAIN-
!                            ING THE INDICES OF THE L-1
!                            CLOSEST NODES TO NPTS(1) IN THE
!                            FIRST L-1 LOCATIONS.
!
! IADJ AND IEND MAY BE CREATED BY SUBROUTINE TRMESH.
!
! INPUT PARAMETERS OTHER THAN NPTS ARE NOT ALTERED BY THIS
!   ROUTINE.
!
! OUTPUT PARAMETERS - NPTS - UPDATED WITH THE INDEX OF THE
!                            L-TH CLOSEST NODE TO NPTS(1) IN
!                            POSITION L UNLESS IER = 1.
!
!                       DF - INCREASING FUNCTION (NEGATIVE
!                            COSINE) OF THE ANGULAR DISTANCE
!                            BETWEEN NPTS(1) AND NPTS(L)
!                            UNLESS IER = 1.
!
!                      IER - ERROR INDICATOR
!                            IER = 0 IF NO ERRORS WERE EN-
!                                    COUNTERED.
!                            IER = 1 IF L IS OUT OF RANGE.
!
! MODULES REFERENCED BY GETNP - NONE
!
! INTRINSIC FUNCTION CALLED BY GETNP - IABS
!
!***********************************************************
!
      INTEGER LM1, N1, I, NI, NP, INDF, INDL, INDX, NB
      DOUBLE PRECISION    X1, Y1, Z1, DNP, DNB
!
! LOCAL PARAMETERS -
!
! LM1 =      L - 1
! N1 =       NPTS(1)
! I =        NPTS INDEX AND DO-LOOP INDEX
! NI =       NPTS(I)
! NP =       CANDIDATE FOR NPTS(L)
! INDF =     IADJ INDEX OF THE FIRST NEIGHBOR OF NI
! INDL =     IADJ INDEX OF THE LAST NEIGHBOR OF NI
! INDX =     IADJ INDEX IN THE RANGE INDF,...,INDL
! NB =       NEIGHBOR OF NI AND CANDIDATE FOR NP
! X1,Y1,Z1 = COORDINATES OF N1
! DNP,DNB =  NEGATIVE COSINES OF THE ANGULAR DISTANCES FROM
!              N1 TO NP AND TO NB, RESPECTIVELY
!
      NP = 0
      LM1 = L - 1
      IF (LM1 .LT. 1) GO TO 4
      IER = 0
      N1 = NPTS(1)
      X1 = X(N1)
      Y1 = Y(N1)
      Z1 = Z(N1)
!
! MARK THE ELEMENTS OF NPTS
!
      DO 1 I = 1,LM1
        NI = NPTS(I)
    1   IEND(NI) = -IEND(NI)
!
! CANDIDATES FOR NP = NPTS(L) ARE THE UNMARKED NEIGHBORS
!   OF NODES IN NPTS.  DNP IS INITIALIZED TO -COS(PI) --
!   THE MAXIMUM DISTANCE.
!
      DNP = 1.
!
! LOOP ON NODES NI IN NPTS
!
      DO 2 I = 1,LM1
        NI = NPTS(I)
        INDF = 1
        IF (NI .GT. 1) INDF = IABS(IEND(NI-1)) + 1
        INDL = -IEND(NI)
!
! LOOP ON NEIGHBORS NB OF NI
!
        DO 2 INDX = INDF,INDL
          NB = IADJ(INDX)
          IF (NB .EQ. 0  .OR.  IEND(NB) .LT. 0) GO TO 2
!
! NB IS AN UNMARKED NEIGHBOR OF NI.  REPLACE NP IF NB IS
!   CLOSER TO N1.
!
          DNB = -(X(NB)*X1 + Y(NB)*Y1 + Z(NB)*Z1)
          IF (DNB .GE. DNP) GO TO 2
          NP = NB
          DNP = DNB
    2     CONTINUE
      NPTS(L) = NP
      DF = DNP
!
! UNMARK THE ELEMENTS OF NPTS
!
      DO 3 I = 1,LM1
        NI = NPTS(I)
    3   IEND(NI) = -IEND(NI)
      RETURN
!
! L IS OUT OF RANGE
!
    4 IER = 1
      RETURN
      END
      SUBROUTINE GIVENS ( A,B, C,S)
      DOUBLE PRECISION A, B, C, S
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE CONSTRUCTS THE GIVENS PLANE ROTATION --
!     ( C  S)
! G = (     ) WHERE C*C + S*S = 1 -- WHICH ZEROS THE SECOND
!     (-S  C)
! ENTRY OF THE 2-VECTOR (A B)-TRANSPOSE.  A CALL TO GIVENS
! IS NORMALLY FOLLOWED BY A CALL TO ROTATE WHICH APPLIES
! THE TRANSFORMATION TO A 2 BY N MATRIX.  THIS ROUTINE WAS
! TAKEN FROM LINPACK.
!
! INPUT PARAMETERS - A,B - COMPONENTS OF THE 2-VECTOR TO BE
!                          ROTATED.
!
! OUTPUT PARAMETERS -  A - OVERWRITTEN BY R = +/-SQRT(A*A
!                          + B*B)
!
!                      B - OVERWRITTEN BY A VALUE Z WHICH
!                          ALLOWS C AND S TO BE RECOVERED
!                          AS FOLLOWS -
!                          C = SQRT(1-Z*Z), S=Z IF ABS(Z)
!                              .LE. 1.
!                          C = 1/Z, S = SQRT(1-C*C) IF
!                              ABS(Z) .GT. 1.
!
!                      C - +/-(A/R)
!
!                      S - +/-(B/R)
!
! MODULES REFERENCED BY GIVENS - NONE
!
! INTRINSIC FUNCTIONS CALLED BY GIVENS - ABS, SQRT
!
!***********************************************************
!
      DOUBLE PRECISION AA, BB, R, U, V
!
! LOCAL PARAMETERS -
!
! AA,BB = LOCAL COPIES OF A AND B
! R =     C*A + S*B = +/-SQRT(A*A+B*B)
! U,V =   VARIABLES USED TO SCALE A AND B FOR COMPUTING R
!
      AA = A
      BB = B
      IF (ABS(AA) .LE. ABS(BB)) GO TO 1
!
! ABS(A) .GT. ABS(B)
!
      U = AA + AA
      V = BB/U
      R = SQRT(.25 + V*V) * U
      C = AA/R
      S = V * (C + C)
!
! NOTE THAT R HAS THE SIGN OF A, C .GT. 0, AND S HAS
!   SIGN(A)*SIGN(B)
!
      B = S
      A = R
      RETURN
!
! ABS(A) .LE. ABS(B)
!
    1 IF (BB .EQ. 0.) GO TO 2
      U = BB + BB
      V = AA/U
!
! STORE R IN A
!
      A = SQRT(.25 + V*V) * U
      S = BB/A
      C = V * (S + S)
!
! NOTE THAT R HAS THE SIGN OF B, S .GT. 0, AND C HAS
!   SIGN(A)*SIGN(B)
!
      B = 1.
      IF (C .NE. 0.) B = 1./C
      RETURN
!
! A = B = 0.
!
    2 C = 1.
      S = 0.
      RETURN
      END
      SUBROUTINE INTADD (KK,I1,I2,I3, IADJ,IEND )
      INTEGER KK, I1, I2, I3, IADJ(1), IEND(KK)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE ADDS AN INTERIOR NODE TO A TRIANGULATION
! OF A SET OF KK-1 POINTS ON THE UNIT SPHERE.  IADJ AND IEND
! ARE UPDATED WITH THE INSERTION OF NODE KK IN THE TRIANGLE
! WHOSE VERTICES ARE I1, I2, AND I3.
!
! INPUT PARAMETERS -        KK - INDEX OF NODE TO BE
!                                INSERTED.  KK .GE. 4.
!
!                     I1,I2,I3 - INDICES OF THE VERTICES OF
!                                A TRIANGLE CONTAINING NODE
!                                KK - IN COUNTERCLOCKWISE
!                                ORDER.
!
!                         IADJ - SET OF ADJACENCY LISTS
!                                OF NODES IN THE MESH.
!
!                         IEND - POINTERS TO THE ENDS OF
!                                ADJACENCY LISTS IN IADJ FOR
!                                EACH NODE IN THE MESH.
!
!   IADJ AND IEND MAY BE CREATED BY TRMESH AND MUST CONTAIN
! THE VERTICES I1, I2, AND I3.  I1,I2,I3 MAY BE DETERMINED
! BY TRFIND.
!
! KK, I1, I2, AND I3 ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS - IADJ,IEND - UPDATED WITH THE ADDITION
!                                 OF NODE KK AS THE LAST
!                                 ENTRY.  NODE KK WILL BE
!                                 CONNECTED TO NODES I1, I2,
!                                 AND I3.  NO OPTIMIZATION
!                                 OF THE MESH IS PERFORMED.
!
! MODULE REFERENCED BY INTADD - SHIFTD
!
! INTRINSIC FUNCTION CALLED BY INTADD - MOD
!
!***********************************************************
!
      INTEGER K, KM1, N(3), NFT(3), IP1, IP2, IP3, INDX, NF, &
              NL, N1, N2, IMIN, IMAX, I, ITEMP
!
! LOCAL PARAMETERS -
!
! K =           LOCAL COPY OF KK
! KM1 =         K - 1
! N =           VECTOR CONTAINING I1, I2, I3
! NFT =         POINTERS TO THE TOPS OF THE 3 SETS OF IADJ
!                 ELEMENTS TO BE SHIFTED DOWNWARD
! IP1,IP2,IP3 = PERMUTATION INDICES FOR N AND NFT
! INDX =        INDEX FOR IADJ AND N
! NF,NL =       INDICES OF FIRST AND LAST ENTRIES IN IADJ
!                 TO BE SHIFTED DOWN
! N1,N2 =       FIRST 2 VERTICES OF A NEW TRIANGLE --
!                 (N1,N2,KK)
! IMIN,IMAX =   BOUNDS ON DO-LOOP INDEX -- FIRST AND LAST
!                 ELEMENTS OF IEND TO BE INCREMENTED
! I =           DO-LOOP INDEX
! ITEMP =       TEMPORARY STORAGE LOCATION
!
      K = KK
!
! INITIALIZATION
!
      N(1) = I1
      N(2) = I2
      N(3) = I3
!
! SET UP NFT
!
      DO 2 I = 1,3
        N1 = N(I)
        INDX = MOD(I,3) + 1
        N2 = N(INDX)
        INDX = IEND(N1) + 1
!
! FIND THE INDEX OF N2 AS A NEIGHBOR OF N1
!
    1   INDX = INDX - 1
        IF (IADJ(INDX) .NE. N2) GO TO 1
    2   NFT(I) = INDX + 1
!
! ORDER THE VERTICES BY DECREASING MAGNITUDE -
!   N(IP(I+1)) PRECEDES N(IP(I)) IN IEND FOR
!   I = 1,2.
!
      IP1 = 1
      IP2 = 2
      IP3 = 3
      IF ( N(2) .LE. N(1) ) GO TO 3
      IP1 = 2
      IP2 = 1
    3 IF ( N(3) .LE. N(IP1) ) GO TO 4
      IP3 = IP1
      IP1 = 3
    4 IF ( N(IP3) .LE. N(IP2) )  GO TO 5
      ITEMP = IP2
      IP2 = IP3
      IP3 = ITEMP
!
! ADD NODE K TO THE ADJACENCY LISTS OF EACH VERTEX AND
!   UPDATE IEND.  FOR EACH VERTEX, A SET OF IADJ ELEMENTS
!   IS SHIFTED DOWNWARD AND K IS INSERTED.  SHIFTING STARTS
!   AT THE END OF THE ARRAY.
!
    5 KM1 = K - 1
      NL = IEND(KM1)
      NF = NFT(IP1)
      IF (NF .LE. NL) CALL SHIFTD(NF,NL,3, IADJ)
      IADJ(NF+2) = K
      IMIN = N(IP1)
      IMAX = KM1
      DO 6 I = IMIN,IMAX
    6   IEND(I) = IEND(I) + 3
!
      NL = NF - 1
      NF = NFT(IP2)
      CALL SHIFTD(NF,NL,2, IADJ)
      IADJ(NF+1) = K
      IMAX = IMIN - 1
      IMIN = N(IP2)
      DO 7 I = IMIN,IMAX
    7   IEND(I) = IEND(I) + 2
!
      NL = NF - 1
      NF = NFT(IP3)
      CALL SHIFTD(NF,NL,1, IADJ)
      IADJ(NF) = K
      IMAX = IMIN - 1
      IMIN = N(IP3)
      DO 8 I = IMIN,IMAX
    8   IEND(I) = IEND(I) + 1
!
! ADD NODE K TO IEND AND ITS NEIGHBORS TO IADJ
!
      INDX = IEND(KM1)
      IEND(K) = INDX + 3
      DO 9 I = 1,3
        INDX = INDX + 1
    9   IADJ(INDX) = N(I)
      RETURN
      END
      SUBROUTINE INTRNN (N,PLAT,PLON,X,Y,Z,W,IADJ,IEND, IST, &
                         PW,IER)
      INTEGER N, IADJ(6*(N-1)), IEND(N), IST, IER
      DOUBLE PRECISION    PLAT, PLON, X(N), Y(N), Z(N), W(N), PW
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
! (hacked version of INTRC0, jeff whitaker 09/2010)
!
!   GIVEN A TRIANGULATION OF A SET OF NODES ON THE UNIT
! SPHERE, ALONG WITH DATA VALUES AT THE NODES, THIS SUB-
! ROUTINE COMPUTES THE VALUE AT A POINT P BY NEAREST
! NEIGHBOR INTERPOLATION.
!
! INPUT PARAMETERS -     N - NUMBER OF NODES IN THE TRIANGU-
!                            LATION.  N .GE. 3.
!
!                PLAT,PLON - LATITUDE AND LONGITUDE OF P IN
!                            RADIANS.
!
!                    X,Y,Z - VECTORS CONTAINING CARTESIAN
!                            COORDINATES OF THE NODES.
!
!                        W - VECTOR CONTAINING DATA VALUES
!                            AT THE NODES.  W(I) IS ASSOCI-
!                            ATED WITH (X(I),Y(I),Z(I)) FOR
!                            I = 1,...,N.
!
!                IADJ,IEND - TRIANGULATION DATA STRUCTURE
!                            CREATED BY SUBROUTINE TRMESH.
!
!                      IST - INDEX OF THE STARTING NODE IN
!                            THE SEARCH FOR A TRIANGLE CON-
!                            TAINING P.  1 .LE. IST .LE. N.
!                            THE OUTPUT VALUE OF IST FROM A
!                            PREVIOUS CALL MAY BE A GOOD
!                            CHOICE.
!
! INPUT PARAMETERS OTHER THAN IST ARE NOT ALTERED BY THIS
!   ROUTINE.
!
! OUTPUT PARAMETERS - IST - INDEX OF ONE OF THE VERTICES OF
!                           THE TRIANGLE CONTAINING P (OR
!                           NEAREST P) UNLESS IER = -1 OR
!                           IER = -2.
!
!                      PW - VALUE OF THE INTERPOLATORY
!                           FUNCTION AT P IF IER .LE. 0.
!
!                     IER - ERROR INDICATOR
!                           IER = 0 IF INTERPOLATION WAS
!                                   PERFORMED SUCCESSFULLY.
!                           IER = -3 IF POINT IS EXTERIOR
!                                    TO TRIANGULATION.
!                           IER = -1 IF N .LT. 3 OR IST IS
!                                    OUT OF RANGE.
!                           IER = -2 IF THE NODES ARE COL-
!                                    LINEAR.
!
! MODULES REFERENCED BY INTRNN - TRFIND
!
! INTRINSIC FUNCTIONS CALLED BY INTRNN - COS, SIN
!
!***********************************************************
!
      INTEGER I1, I2, I3, I(1)
      DOUBLE PRECISION    P(3), P1(3), P2(3), P3(3), B1,B2,B3, &
              DIST(3), ARCLEN
!
      IF (N .LT. 3  .OR.  IST .LT. 1  .OR.  IST .GT. N) &
          GO TO 5
!
! TRANSFORM (PLAT,PLON) TO CARTESIAN COORDINATES
!
      P(1) = COS(PLAT)*COS(PLON)
      P(2) = COS(PLAT)*SIN(PLON)
      P(3) = SIN(PLAT)
!
! FIND THE VERTEX INDICES OF A TRIANGLE CONTAINING P
!
      CALL TRFIND(IST,P,X,Y,Z,IADJ,IEND, B1,B2,B3,I1,I2,I3)
      IF (I1 .EQ. 0) GO TO 6
      IST = I1
      IF (I3 .EQ. 0) GO TO 1
!
! P IS CONTAINED IN THE TRIANGLE (I1,I2,I3).  STORE THE
!  VERTEX COORDINATES IN LOCAL VARIABLES
!
      P1(1) = X(I1)
      P1(2) = Y(I1)
      P1(3) = Z(I1)
      P2(1) = X(I2)
      P2(2) = Y(I2)
      P2(3) = Z(I2)
      P3(1) = X(I3)
      P3(2) = Y(I3)
      P3(3) = Z(I3)
      dist(1) =  ARCLEN (P,P1)
      dist(2) =  ARCLEN (P,P2)
      dist(3) =  ARCLEN (P,P3)
      I = minloc(dist)
      if (i(1) .eq. 1) pw = w(i1)
      if (i(1) .eq. 2) pw = w(i2)
      if (i(1) .eq. 3) pw = w(i3)
      IER = 0
      RETURN
!
! N OR IST OUT OF RANGE
!
    5 IER = -1
      RETURN
!
! COLLINEAR NODES
!
    6 IER = -2
      RETURN
!
! P EXTERIOR TO TRIANGULATION
!
    1 IER = -3
      RETURN
      END
      SUBROUTINE INTRC0 (N,PLAT,PLON,X,Y,Z,W,IADJ,IEND, IST, &
                         PW,IER)
      INTEGER N, IADJ(6*(N-1)), IEND(N), IST, IER
      DOUBLE PRECISION    PLAT, PLON, X(N), Y(N), Z(N), W(N), PW
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A TRIANGULATION OF A SET OF NODES ON THE UNIT
! SPHERE, ALONG WITH DATA VALUES AT THE NODES, THIS SUB-
! ROUTINE COMPUTES THE VALUE AT A POINT P OF A CONTINUOUS
! FUNCTION WHICH INTERPOLATES THE DATA VALUES.  THE INTERP-
! OLATORY FUNCTION IS LINEAR ON EACH UNDERLYING TRIANGLE
! (PLANAR TRIANGLE WITH THE SAME VERTICES AS A SPHERICAL
! TRIANGLE).  IF P IS NOT CONTAINED IN A TRIANGLE, AN EX-
! TRAPOLATED VALUE IS TAKEN TO BE THE INTERPOLATED VALUE AT
! THE NEAREST POINT OF THE TRIANGULATION BOUNDARY.
!
! INPUT PARAMETERS -     N - NUMBER OF NODES IN THE TRIANGU-
!                            LATION.  N .GE. 3.
!
!                PLAT,PLON - LATITUDE AND LONGITUDE OF P IN
!                            RADIANS.
!
!                    X,Y,Z - VECTORS CONTAINING CARTESIAN
!                            COORDINATES OF THE NODES.
!
!                        W - VECTOR CONTAINING DATA VALUES
!                            AT THE NODES.  W(I) IS ASSOCI-
!                            ATED WITH (X(I),Y(I),Z(I)) FOR
!                            I = 1,...,N.
!
!                IADJ,IEND - TRIANGULATION DATA STRUCTURE
!                            CREATED BY SUBROUTINE TRMESH.
!
!                      IST - INDEX OF THE STARTING NODE IN
!                            THE SEARCH FOR A TRIANGLE CON-
!                            TAINING P.  1 .LE. IST .LE. N.
!                            THE OUTPUT VALUE OF IST FROM A
!                            PREVIOUS CALL MAY BE A GOOD
!                            CHOICE.
!
! INPUT PARAMETERS OTHER THAN IST ARE NOT ALTERED BY THIS
!   ROUTINE.
!
! OUTPUT PARAMETERS - IST - INDEX OF ONE OF THE VERTICES OF
!                           THE TRIANGLE CONTAINING P (OR
!                           NEAREST P) UNLESS IER = -1 OR
!                           IER = -2.
!
!                      PW - VALUE OF THE INTERPOLATORY
!                           FUNCTION AT P IF IER .LE. 0.
!
!                     IER - ERROR INDICATOR
!                           IER = 0 IF INTERPOLATION WAS
!                                   PERFORMED SUCCESSFULLY.
!                           IER = 1 IF EXTRAPOLATION WAS
!                                   PERFORMED SUCCESSFULLY.
!                           IER = -1 IF N .LT. 3 OR IST IS
!                                    OUT OF RANGE.
!                           IER = -2 IF THE NODES ARE COL-
!                                    LINEAR.
!                           IER = -3 IF P IS NOT IN A TRI-
!                                    ANGLE AND THE ANGLE BE-
!                                    TWEEN P AND THE NEAREST
!                                    BOUNDARY POINT IS AT
!                                    LEAST 90 DEGREES.
!
! MODULES REFERENCED BY INTRC0 - TRFIND
!
! INTRINSIC FUNCTIONS CALLED BY INTRC0 - COS, SIN
!
!***********************************************************
!
      INTEGER I1, I2, I3, N1, N2, INDX
      DOUBLE PRECISION    P(3), P1(3), P2(3), P3(3), B1,B2,B3,SUMM,&
              S12, PTN1, PTN2
!
! LOCAL PARAMETERS -
!
! I1,I2,I3 = VERTEX INDICES RETURNED BY TRFIND
! N1,N2 =    ENDPOINTS OF A BOUNDARY ARC WHICH IS VISIBLE
!              FROM P WHEN P IS NOT CONTAINED IN A TRIANGLE
! INDX =     IADJ INDEX OF N1 AS A NEIGHBOR OF N2
! P =        CARTESIAN COORDINATES OF P
! P1,P2,P3 = CARTESIAN COORDINATES OF I1, I2, AND I3
! B1,B2,B3 = BARYCENTRIC COORDINATES OF THE CENTRAL PROJEC-
!              TION OF P ONTO THE UNDERLYING PLANAR TRIANGLE
!              OR (B1 AND B2) PROJECTION OF Q ONTO THE
!              UNDERLYING LINE SEGMENT N1-N2 WHEN P IS
!              EXTERIOR -- UNNORMALIZED COORDINATES ARE
!              COMPUTED BY TRFIND WHEN P IS IN A TRIANGLE
! SUMM =     QUANTITY USED TO NORMALIZE THE BARYCENTRIC
!              COORDINATES
! S12 =      SCALAR PRODUCT (N1,N2)
! PTN1 =     SCALAR PRODUCT (P,N1)
! PTN2 =     SCALAR PRODUCT (P,N2)
!
      IF (N .LT. 3  .OR.  IST .LT. 1  .OR.  IST .GT. N) &
          GO TO 5
!
! TRANSFORM (PLAT,PLON) TO CARTESIAN COORDINATES
!
      P(1) = COS(PLAT)*COS(PLON)
      P(2) = COS(PLAT)*SIN(PLON)
      P(3) = SIN(PLAT)
!
! FIND THE VERTEX INDICES OF A TRIANGLE CONTAINING P
!
      CALL TRFIND(IST,P,X,Y,Z,IADJ,IEND, B1,B2,B3,I1,I2,I3)
      IF (I1 .EQ. 0) GO TO 6
      IST = I1
      IF (I3 .EQ. 0) GO TO 1
!
! P IS CONTAINED IN THE TRIANGLE (I1,I2,I3).  STORE THE
!   VERTEX COORDINATES IN LOCAL VARIABLES
!
      P1(1) = X(I1)
      P1(2) = Y(I1)
      P1(3) = Z(I1)
      P2(1) = X(I2)
      P2(2) = Y(I2)
      P2(3) = Z(I2)
      P3(1) = X(I3)
      P3(2) = Y(I3)
      P3(3) = Z(I3)
!
! NORMALIZE THE BARYCENTRIC COORDINATES
!
      SUMM = B1 + B2 + B3
      B1 = B1/SUMM
      B2 = B2/SUMM
      B3 = B3/SUMM
      PW = B1*W(I1) + B2*W(I2) + B3*W(I3)
      IER = 0
      RETURN
!
! P IS EXTERIOR TO THE TRIANGULATION, AND I1 AND I2 ARE
!   BOUNDARY NODES WHICH ARE VISIBLE FROM P.  SET PW TO THE
!   INTERPOLATED VALUE AT Q WHERE Q IS THE CLOSEST BOUNDARY
!   POINT TO P.
!
! TRAVERSE THE BOUNDARY STARTING FROM THE RIGHTMOST VISIBLE
!   NODE I1.
!
    1 N1 = I1
      PTN1 = P(1)*X(N1) + P(2)*Y(N1) + P(3)*Z(N1)
      IF (I1 .NE. I2) GO TO 3
!
! ALL BOUNDARY NODES ARE VISIBLE FROM P.  FIND A BOUNDARY
!   ARC N1->N2 SUCH THAT P LEFT (N2 X N1)->N1
!
! COUNTERCLOCKWISE BOUNDARY TRAVERSAL --
!   SET N2 TO THE FIRST NEIGHBOR OF N1
!
    2 INDX = 1
      IF (N1 .NE. 1) INDX = IEND(N1-1) + 1
      N2 = IADJ(INDX)
!
! COMPUTE INNER PRODUCTS (N1,N2) AND (P,N2), AND COMPUTE
!   B2 = DET(P,N1,N2 X N1)
!
      S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
      PTN2 = P(1)*X(N2) + P(2)*Y(N2) + P(3)*Z(N2)
      B2 = PTN2 - S12*PTN1
      IF (B2 .LE. 0.) GO TO 3
!
! P RIGHT (N2 X N1)->N1 -- ITERATE
!
      N1 = N2
      I1 = N1
      PTN1 = PTN2
      GO TO 2
!
! P LEFT (N2 X N1)->N1 WHERE N2 IS THE FIRST NEIGHBOR OF P1
!   CLOCKWISE BOUNDARY TRAVERSAL --
!
    3 N2 = N1
      PTN2 = PTN1
!
! SET N1 TO THE LAST NEIGHBOR OF N2 AND TEST FOR TERMINATION
!
      INDX = IEND(N2) - 1
      N1 = IADJ(INDX)
      IF (N1 .EQ. I1) GO TO 7
!
! COMPUTE INNER PRODUCTS (N1,N2) AND (P,N1)
!
      S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
      PTN1 = P(1)*X(N1) + P(2)*Y(N1) + P(3)*Z(N1)
!
! COMPUTE B2 = DET(P,N1,N2 X N1) = DET(Q,N1,N2 X N1)*(P,Q)
!
      B2 = PTN2 - S12*PTN1
      IF (B2 .LE. 0.) GO TO 3
!
! COMPUTE B1 = DET(P,N2 X N1,N2) = DET(Q,N2 X N1,N2)*(P,Q)
!
      B1 = PTN1 - S12*PTN2
      IF (B1 .GT. 0.) GO TO 4
!
! Q = N2
!
      PW = W(N2)
      IER = 1
      RETURN
!
! P STRICTLY LEFT (N2 X N1)->N2 AND P STRICTLY LEFT
!   N1->(N2 X N1).  THUS Q LIES ON THE INTERIOR OF N1->N2.
!   NORMALIZE THE COORDINATES AND COMPUTE PW.
!
    4 SUMM = B1 + B2
      PW = (B1*W(N1) + B2*W(N2))/SUMM
      IER = 1
      RETURN
!
! N OR IST OUT OF RANGE
!
    5 IER = -1
      RETURN
!
! COLLINEAR NODES
!
    6 IER = -2
      RETURN
!
! THE ANGULAR DISTANCE BETWEEN P AND THE CLOSEST BOUNDARY
!   POINT TO P IS AT LEAST 90 DEGREES
!
    7 IER = -3
      RETURN
      END
      SUBROUTINE INTRC1 (N,PLAT,PLON,X,Y,Z,W,IADJ,IEND, &
                         IFLAG,GRAD, IST, PW,IER)
      INTEGER N, IADJ(6*(N-1)), IEND(N), IFLAG, IST, IER
      DOUBLE PRECISION    PLAT, PLON, X(N), Y(N), Z(N), W(N), &
              GRAD(3,N), PW
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A TRIANGULATION OF A SET OF NODES ON THE UNIT
! SPHERE, ALONG WITH DATA VALUES AT THE NODES, THIS SUB-
! ROUTINE CONSTRUCTS THE VALUE AT A POINT P OF A ONCE CON-
! TINUOUSLY DIFFERENTIABLE FUNCTION WHICH INTERPOLATES THE
! DATA VALUES.  IF THE TRIANGULATION DOES NOT COVER THE
! ENTIRE SPHERE, THE SURFACE IS EXTENDED CONTINUOUSLY BEYOND
! THE BOUNDARY ALLOWING EXTRAPOLATION.
!
! INPUT PARAMETERS -     N - NUMBER OF NODES.  N .GE. 3 AND
!                            N .GE. 7 IF IFLAG = 0.
!
!                PLAT,PLON - LATITUDE AND LONGITUDE OF P IN
!                            RADIANS.
!
!                    X,Y,Z - VECTORS CONTAINING CARTESIAN
!                            COORDINATES OF THE NODES.
!
!                        W - VECTOR CONTAINING DATA VALUES
!                            AT THE NODES.  W(I) IS ASSOCI-
!                            ATED WITH (X(I),Y(I),Z(I)) FOR
!                            I = 1,...,N.
!
!                IADJ,IEND - DATA STRUCTURE REPRESENTING THE
!                            TRIANGULATION.  SEE SUBROUTINE
!                            TRMESH.
!
!                    IFLAG - OPTION INDICATOR
!                            IFLAG = 0 IF INTRC1 IS TO PRO-
!                                      VIDE ESTIMATED GRAD-
!                                      IENTS (FROM GRADL).
!                                      N .GE. 7 IN THIS
!                                      CASE.
!                            IFLAG = 1 IF GRADIENTS ARE PRO-
!                                      VIDED IN GRAD.  THIS
!                                      IS MORE EFFICIENT IF
!                                      INTRC1 IS TO BE
!                                      CALLED SEVERAL TIMES.
!
!                     GRAD - ARRAY DIMENSIONED 3 BY N WHOSE
!                            I-TH COLUMN CONTAINS AN ESTI-
!                            MATED GRADIENT AT NODE I IF
!                            IFLAG = 1 (SEE SUBROUTINE
!                            GRADL).  GRAD MAY BE A DUMMY
!                            VARIABLE (NOT USED) IF IFLAG
!                            = 0.
!
!                      IST - INDEX OF THE STARTING NODE IN
!                            THE SEARCH FOR A TRIANGLE CON-
!                            TAINING P.  1 .LE. IST .LE. N.
!                            THE OUTPUT VALUE OF IST FROM A
!                            PREVIOUS CALL MAY BE A GOOD
!                            CHOICE.
!
! INPUT PARAMETERS OTHER THAN IST ARE NOT ALTERED BY THIS
!   ROUTINE.
!
! OUTPUT PARAMETERS - IST - INDEX OF ONE OF THE VERTICES OF
!                           THE TRIANGLE CONTAINING P (OR
!                           NEAREST P) UNLESS IER = -1 OR
!                           IER = -2.
!
! OUTPUT PARAMETERS -  PW - INTERPOLATED VALUE AT P IF
!                           IER .GE. 0.
!
!                     IER - ERROR INDICATOR
!                           IER = 0 IF INTERPOLATION WAS
!                                   PERFORMED SUCCESSFULLY.
!                           IER = 1 IF EXTRAPOLATION WAS
!                                   PERFORMED SUCCESSFULLY.
!                           IER = -1 IF N, IFLAG, OR IST IS
!                                    OUT OF RANGE.
!                           IER = -2 IF THE NODES ARE COL-
!                                    LINEAR.
!                           IER = -3 IF THE ANGULAR DISTANCE
!                                    BETWEEN P AND THE NEAR-
!                                    EST POINT OF THE TRIAN-
!                                    GULATION IS AT LEAST 90
!                                    DEGREES.
!
! MODULES REFERENCED BY INTRC1 - TRFIND, WVAL, ARCINT,
!                                ARCLEN,
!             (AND OPTIONALLY)   GRADL, GETNP, CONSTR,
!                                APLYR, SETUP, GIVENS,
!                                ROTATE, APLYRT
!
! INTRINSIC FUNCTIONS CALLED BY INTRC1 - COS, SIN, SQRT
!
!***********************************************************
!
      INTEGER NN, I1, I2, I3, I, IERR, N1, N2, INDX
      DOUBLE PRECISION    P(3), P1(3), P2(3), P3(3), W1, W2, W3, G1(3), &
              G2(3), G3(3), B1, B2, B3, SUMM, DUM(3), S12, &
              PTN1, PTN2, Q(3), QNORM, WQ, GQ(3), A, PTGQ, &
              GQN, ARCLEN
!
! LOCAL PARAMETERS -
!
! NN =       LOCAL COPY OF N
! I1,I2,I3 = VERTEX INDICES RETURNED BY TRFIND
! I =        DO-LOOP INDEX
! IERR =     ERROR FLAG FOR CALLS TO GRADL
! N1,N2 =    INDICES OF THE ENDPOINTS OF A BOUNDARY ARC WHEN
!              P IS EXTERIOR (NOT CONTAINED IN A TRIANGLE)
! INDX =     IADJ INDEX OF N2 AS A NEIGHBOR OF N1 OR VICE-
!              VERSA
! P =        CARTESIAN COORDINATES OF P
! P1,P2,P3 = CARTESIAN COORDINATES OF THE VERTICES I1, I2,
!              AND I3, OR (P1 AND P2) COORDINATES OF N1 AND
!              N2 IF P IS EXTERIOR
! W1,W2,W3 = DATA VALUES ASSOCIATED WITH I1, I2, AND I3, OR
!              (W1 AND W2) VALUES ASSOCIATED WITH N1 AND
!              N2 IF P IS EXTERIOR
! G1,G2,G3 = GRADIENTS AT I1, I2, AND I3, OR (G1 AND G2) AT
!              N1 AND N2
! B1,B2,B3 = BARYCENTRIC COORDINATES OF THE CENTRAL PROJEC-
!              TION OF P ONTO THE UNDERLYING PLANAR TRIANGLE
!              OR (B1 AND B2) PROJECTION OF Q ONTO THE
!              UNDERLYING LINE SEGMENT N1-N2 WHEN P IS
!              EXTERIOR -- UNNORMALIZED COORDINATES ARE
!              COMPUTED BY TRFIND WHEN P IS IN A TRIANGLE
! SUMM =     QUANTITY USED TO NORMALIZE THE BARYCENTRIC
!             COORDINATES
! DUM =      DUMMY PARAMETER FOR WVAL AND ARCINT
! S12 =      SCALAR PRODUCT (N1,N2) -- FACTOR OF B1 AND B2
! PTN1 =     SCALAR PRODUCT (P,N1) -- FACTOR OF B1 AND B2
! PTN2 =     SCALAR PRODUCT (P,N2) -- FACTOR OF B1 AND B2
! Q =        CLOSEST BOUNDARY POINT TO P WHEN P IS EXTERIOR
! QNORM =    FACTOR USED TO NORMALIZE Q
! WQ,GQ =    INTERPOLATED VALUE AND GRADIENT AT Q
! A =        ANGULAR SEPARATION BETWEEN P AND Q
! PTGQ =     SCALAR PRODUCT (P,GQ) -- FACTOR OF THE COMPONENT
!              OF GQ IN THE DIRECTION Q->P
! GQN =      NEGATIVE OF THE COMPONENT OF GQ IN THE DIRECTION
!              Q->P
!
      NN = N
      IF (NN .LT. 3  .OR.  (IFLAG .NE. 1  .AND.  NN .LT. 7) &
         .OR.  IFLAG .LT. 0  .OR.  IFLAG .GT. 1  .OR. &
         IST .LT. 1  .OR.  IST .GT. NN) GO TO 16
!
! TRANSFORM (PLAT,PLON) TO CARTESIAN COORDINATES
!
      P(1) = COS(PLAT)*COS(PLON)
      P(2) = COS(PLAT)*SIN(PLON)
      P(3) = SIN(PLAT)
!
! LOCATE P WITH RESPECT TO THE TRIANGULATION
!
      CALL TRFIND(IST,P,X,Y,Z,IADJ,IEND, B1,B2,B3,I1,I2,I3)
      IF (I1 .EQ. 0) GO TO 17
      IST = I1
      IF (I3 .EQ. 0) GO TO 4
!
! P IS CONTAINED IN THE TRIANGLE (I1,I2,I3).  STORE THE DATA
!   VALUES AND VERTEX COORDINATES IN LOCAL VARIABLES.
!
      W1 = W(I1)
      W2 = W(I2)
      W3 = W(I3)
      P1(1) = X(I1)
      P1(2) = Y(I1)
      P1(3) = Z(I1)
      P2(1) = X(I2)
      P2(2) = Y(I2)
      P2(3) = Z(I2)
      P3(1) = X(I3)
      P3(2) = Y(I3)
      P3(3) = Z(I3)
      IF (IFLAG .NE. 1) GO TO 2
!
! GRADIENTS ARE USER-PROVIDED
!
      DO 1 I = 1,3
        G1(I) = GRAD(I,I1)
        G2(I) = GRAD(I,I2)
    1   G3(I) = GRAD(I,I3)
      GO TO 3
!
! COMPUTE GRADIENT ESTIMATES AT THE VERTICES
!
    2 CALL GRADL(NN,I1,X,Y,Z,W,IADJ,IEND, G1,IERR)
      IF (IERR .LT. 0) GO TO 17
      CALL GRADL(NN,I2,X,Y,Z,W,IADJ,IEND, G2,IERR)
      IF (IERR .LT. 0) GO TO 17
      CALL GRADL(NN,I3,X,Y,Z,W,IADJ,IEND, G3,IERR)
      IF (IERR .LT. 0) GO TO 17
!
! NORMALIZE THE COORDINATES
!
    3 SUMM = B1 + B2 + B3
      B1 = B1/SUMM
      B2 = B2/SUMM
      B3 = B3/SUMM
      CALL WVAL(B1,B2,B3,P1,P2,P3,W1,W2,W3,G1,G2,G3,0, PW, &
                DUM)
      IER = 0
      RETURN
!
! P IS EXTERIOR TO THE TRIANGULATION, AND I1 AND I2 ARE
!   BOUNDARY NODES WHICH ARE VISIBLE FROM P.  EXTRAPOLATE TO
!   P BY LINEAR (WITH RESPECT TO ARC-LENGTH) INTERPOLATION
!   OF THE VALUE AND DIRECTIONAL DERIVATIVE (GRADIENT COMP-
!   ONENT IN THE DIRECTION Q->P) OF THE INTERPOLATORY
!   SURFACE AT Q WHERE Q IS THE CLOSEST BOUNDARY POINT TO P.
!
! DETERMINE Q BY TRAVERSING THE BOUNDARY STARTING FROM I1
!
    4 N1 = I1
      PTN1 = P(1)*X(N1) + P(2)*Y(N1) + P(3)*Z(N1)
      IF (I1 .NE. I2) GO TO 6
!
! ALL BOUNDARY NODES ARE VISIBLE FROM P.  FIND A BOUNDARY
!   ARC N1->N2 SUCH THAT P LEFT (N2 X N1)->N1
!
! COUNTERCLOCKWISE BOUNDARY TRAVERSAL --
!   SET N2 TO THE FIRST NEIGHBOR OF N1
!
    5 INDX = 1
      IF (N1 .NE. 1) INDX = IEND(N1-1) + 1
      N2 = IADJ(INDX)
!
! COMPUTE INNER PRODUCTS (N1,N2) AND (P,N2), AND COMPUTE
!   B2 = DET(P,N1,N2 X N1)
!
      S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
      PTN2 = P(1)*X(N2) + P(2)*Y(N2) + P(3)*Z(N2)
      B2 = PTN2 - S12*PTN1
      IF (B2 .LE. 0.) GO TO 6
!
! P RIGHT (N2 X N1)->N1 -- ITERATE
!
      N1 = N2
      I1 = N1
      PTN1 = PTN2
      GO TO 5
!
! P LEFT (N2 X N1)->N1 WHERE N2 IS THE FIRST NEIGHBOR OF N1
!   CLOCKWISE BOUNDARY TRAVERSAL --
!
    6 N2 = N1
      PTN2 = PTN1
!
! SET N1 TO THE LAST NEIGHBOR OF N2 AND TEST FOR TERMINATION
!
      INDX = IEND(N2) - 1
      N1 = IADJ(INDX)
      IF (N1 .EQ. I1) GO TO 18
!
! COMPUTE INNER PRODUCTS (N1,N2) AND (P,N1)
!
      S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
      PTN1 = P(1)*X(N1) + P(2)*Y(N1) + P(3)*Z(N1)
!
! COMPUTE B2 = DET(P,N1,N2 X N1) = DET(Q,N1,N2 X N1)*(P,Q)
!
      B2 = PTN2 - S12*PTN1
      IF (B2 .LE. 0.) GO TO 6
!
! COMPUTE B1 = DET(P,N2 X N1,N2) = DET(Q,N2 X N1,N2)*(P,Q)
!
      B1 = PTN1 - S12*PTN2
      IF (B1 .GT. 0.) GO TO 10
!
! Q = N2.  STORE VALUE, COORDINATES, AND AND GRADIENT AT Q.
!
      WQ = W(N2)
      Q(1) = X(N2)
      Q(2) = Y(N2)
      Q(3) = Z(N2)
      IF (IFLAG .NE. 1) GO TO 8
      DO 7 I = 1,3
    7   GQ(I) = GRAD(I,N2)
      GO TO 9
    8 CALL GRADL(NN,N2,X,Y,Z,W,IADJ,IEND, GQ,IERR)
      IF (IERR .LT. 0) GO TO 17
!
! EXTRAPOLATE TO P -- PW = WQ + A*(GQ,Q X (PXQ)/SIN(A))
!   WHERE A IS THE ANGULAR SEPARATION BETWEEN Q AND P,
!   AND SIN(A) IS THE MAGNITUDE OF P X Q.
!
    9 A = ARCLEN(Q,P)
      PTGQ = P(1)*GQ(1) + P(2)*GQ(2) + P(3)*GQ(3)
      PW = WQ
      IF (A .NE. 0.) PW = PW + PTGQ*A/SIN(A)
      IER = 1
      RETURN
!
! P STRICTLY LEFT (N2 X N1)->N2 AND P STRICTLY LEFT
!   N1->(N2 X N1).  THUS Q LIES ON THE INTERIOR OF N1->N2.
!   STORE DATA VALUES AND COORDINATES OF N1 AND N2 IN
!   LOCAL VARIABLES
!
   10 W1 = W(N1)
      W2 = W(N2)
      P1(1) = X(N1)
      P1(2) = Y(N1)
      P1(3) = Z(N1)
      P2(1) = X(N2)
      P2(2) = Y(N2)
      P2(3) = Z(N2)
!
! COMPUTE THE CENTRAL PROJECTION OF Q ONTO P2-P1 AND
!   NORMALIZE TO OBTAIN Q
!
      QNORM = 0.
      DO 11 I = 1,3
        Q(I) = B1*P1(I) + B2*P2(I)
   11   QNORM = QNORM + Q(I)*Q(I)
      QNORM = SQRT(QNORM)
      DO 12 I = 1,3
   12   Q(I) = Q(I)/QNORM
!
! STORE OR COMPUTE GRADIENTS AT N1 AND N2
!
      IF (IFLAG .NE. 1) GO TO 14
      DO 13 I = 1,3
        G1(I) = GRAD(I,N1)
   13   G2(I) = GRAD(I,N2)
      GO TO 15
   14 CALL GRADL(NN,N1,X,Y,Z,W,IADJ,IEND, G1,IERR)
      IF (IERR .LT. 0) GO TO 17
      CALL GRADL(NN,N2,X,Y,Z,W,IADJ,IEND, G2,IERR)
      IF (IERR .LT. 0) GO TO 17
!
! COMPUTE AN INTERPOLATED VALUE AND NORMAL GRADIENT-
!   COMPONENT AT Q
!
   15 CALL ARCINT(Q,P1,P2,W1,W2,G1,G2, WQ,DUM,GQN)
!
! EXTRAPOLATE TO P -- THE NORMAL GRADIENT COMPONENT GQN IS
!   THE NEGATIVE OF THE COMPONENT IN THE DIRECTION Q->P.
!
      PW = WQ - GQN*ARCLEN(Q,P)
      IER = 1
      RETURN
!
! N, IFLAG, OR IST OUT OF RANGE
!
   16 IER = -1
      RETURN
!
! COLLINEAR NODES
!
   17 IER = -2
      RETURN
!
! THE DISTANCE BETWEEN P AND THE CLOSEST BOUNDARY POINT
!   IS AT LEAST 90 DEGREES
!
   18 IER = -3
      RETURN
      END
      SUBROUTINE PERMUT (N,IP, A )
      INTEGER N, IP(N)
      DOUBLE PRECISION    A(N)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE APPLIES A SET OF PERMUTATIONS TO A VECTOR.
!
! INPUT PARAMETERS -  N - LENGTH OF A AND IP.
!
!                    IP - VECTOR CONTAINING THE SEQUENCE OF
!                         INTEGERS 1,...,N PERMUTED IN THE
!                         SAME FASHION THAT A IS TO BE PER-
!                         MUTED.
!
!                     A - VECTOR TO BE PERMUTED.
!
! N AND IP ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETER -  A - REORDERED VECTOR REFLECTING THE
!                         PERMUTATIONS DEFINED BY IP.
!
! MODULES REFERENCED BY PERMUT - NONE
!
!***********************************************************
!
      INTEGER NN, K, J, IPJ
      DOUBLE PRECISION    TEMP
!
! LOCAL PARAMETERS -
!
! NN =   LOCAL COPY OF N
! K =    INDEX FOR IP AND FOR THE FIRST ELEMENT OF A IN A
!          PERMUTATION
! J =    INDEX FOR IP AND A, J .GE. K
! IPJ =  IP(J)
! TEMP = TEMPORARY STORAGE FOR A(K)
!
      NN = N
      IF (NN .LT. 2) RETURN
      K = 1
!
! LOOP ON PERMUTATIONS
!
    1 J = K
      TEMP = A(K)
!
! APPLY PERMUTATION TO A.  IP(J) IS MARKED (MADE NEGATIVE)
!   AS BEING INCLUDED IN THE PERMUTATION.
!
    2 IPJ = IP(J)
      IP(J) = -IPJ
      IF (IPJ .EQ. K) GO TO 3
      A(J) = A(IPJ)
      J = IPJ
      GO TO 2
    3 A(J) = TEMP
!
! SEARCH FOR AN UNMARKED ELEMENT OF IP
!
    4 K = K + 1
      IF (K .GT. NN) GO TO 5
      IF (IP(K) .GT. 0) GO TO 1
      GO TO 4
!
! ALL PERMUTATIONS HAVE BEEN APPLIED.  UNMARK IP.
!
    5 DO 6 K = 1,NN
    6   IP(K) = -IP(K)
      RETURN
      END
      SUBROUTINE QSORT (N,X, IND)
      INTEGER N, IND(N)
      DOUBLE PRECISION    X(N)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE USES AN ORDER N*LOG(N) QUICK SORT TO
! SORT THE DOUBLE PRECISION ARRAY X INTO INCREASING ORDER.  THE ALGOR-
! ITHM IS AS FOLLOWS.  IND IS INITIALIZED TO THE ORDERED
! SEQUENCE OF INDICES 1,...,N, AND ALL INTERCHANGES ARE
! APPLIED TO IND.  X IS DIVIDED INTO TWO PORTIONS BY PICKING
! A CENTRAL ELEMENT T.  THE FIRST AND LAST ELEMENTS ARE COM-
! PARED WITH T, AND INTERCHANGES ARE APPLIED AS NECESSARY SO
! THAT THE THREE VALUES ARE IN ASCENDING ORDER.  INTER-
! CHANGES ARE THEN APPLIED SO THAT ALL ELEMENTS GREATER THAN
! T ARE IN THE UPPER PORTION OF THE ARRAY AND ALL ELEMENTS
! LESS THAN T ARE IN THE LOWER PORTION.  THE UPPER AND LOWER
! INDICES OF ONE OF THE PORTIONS ARE SAVED IN LOCAL ARRAYS,
! AND THE PROCESS IS REPEATED ITERATIVELY ON THE OTHER
! PORTION.  WHEN A PORTION IS COMPLETELY SORTED, THE PROCESS
! BEGINS AGAIN BY RETRIEVING THE INDICES BOUNDING ANOTHER
! UNSORTED PORTION.
!
! INPUT PARAMETERS -   N - LENGTH OF THE ARRAY X.
!
!                      X - VECTOR OF LENGTH N TO BE SORTED.
!
!                    IND - VECTOR OF LENGTH .GE. N.
!
! N AND X ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETER - IND - SEQUENCE OF INDICES 1,...,N
!                          PERMUTED IN THE SAME FASHION AS X
!                          WOULD BE.  THUS, THE ORDERING ON
!                          X IS DEFINED BY Y(I) = X(IND(I)).
!
! MODULES REFERENCED BY QSORT - NONE
!
! INTRINSIC FUNCTIONS CALLED BY QSORT - INT, FLOAT
!
!***********************************************************
!
! NOTE -- IU AND IL MUST BE DIMENSIONED .GE. LOG(N) WHERE
!         LOG HAS BASE 2.
!
!***********************************************************
!
      INTEGER IU(21), IL(21)
      INTEGER M, I, J, K, L, IJ, IT, ITT, INDX
      DOUBLE PRECISION    R, T
!
! LOCAL PARAMETERS -
!
! IU,IL =  TEMPORARY STORAGE FOR THE UPPER AND LOWER
!            INDICES OF PORTIONS OF THE ARRAY X
! M =      INDEX FOR IU AND IL
! I,J =    LOWER AND UPPER INDICES OF A PORTION OF X
! K,L =    INDICES IN THE RANGE I,...,J
! IJ =     RANDOMLY CHOSEN INDEX BETWEEN I AND J
! IT,ITT = TEMPORARY STORAGE FOR INTERCHANGES IN IND
! INDX =   TEMPORARY INDEX FOR X
! R =      PSEUDO RANDOM NUMBER FOR GENERATING IJ
! T =      CENTRAL ELEMENT OF X
!
      IF (N .LE. 0) RETURN
!
! INITIALIZE IND, M, I, J, AND R
!
      DO 1 I = 1,N
    1   IND(I) = I
      M = 1
      I = 1
      J = N
      R = .375
!
! TOP OF LOOP
!
    2 IF (I .GE. J) GO TO 10
      IF (R .GT. .5898437) GO TO 3
      R = R + .0390625
      GO TO 4
    3 R = R - .21875
!
! INITIALIZE K
!
    4 K = I
!
! SELECT A CENTRAL ELEMENT OF X AND SAVE IT IN T
!
      IJ = I + INT(R*FLOAT(J-I))
      IT = IND(IJ)
      T = X(IT)
!
! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
!   INTERCHANGE IT WITH T
!
      INDX = IND(I)
      IF (X(INDX) .LE. T) GO TO 5
      IND(IJ) = INDX
      IND(I) = IT
      IT = INDX
      T = X(IT)
!
! INITIALIZE L
!
    5 L = J
!
! IF THE LAST ELEMENT OF THE ARRAY IS LESS THAN T,
!   INTERCHANGE IT WITH T
!
      INDX = IND(J)
      IF (X(INDX) .GE. T) GO TO 7
      IND(IJ) = INDX
      IND(J) = IT
      IT = INDX
      T = X(IT)
!
! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
!   INTERCHANGE IT WITH T
!
      INDX = IND(I)
      IF (X(INDX) .LE. T) GO TO 7
      IND(IJ) = INDX
      IND(I) = IT
      IT = INDX
      T = X(IT)
      GO TO 7
!
! INTERCHANGE ELEMENTS K AND L
!
    6 ITT = IND(L)
      IND(L) = IND(K)
      IND(K) = ITT
!
! FIND AN ELEMENT IN THE UPPER PART OF THE ARRAY WHICH IS
!   NOT LARGER THAN T
!
    7 L = L - 1
      INDX = IND(L)
      IF (X(INDX) .GT. T) GO TO 7
!
! FIND AN ELEMENT IN THE LOWER PART OF THE ARRAY WHCIH IS
!   NOT SMALLER THAN T
!
    8 K = K + 1
      INDX = IND(K)
      IF (X(INDX) .LT. T) GO TO 8
!
! IF K .LE. L, INTERCHANGE ELEMENTS K AND L
!
      IF (K .LE. L) GO TO 6
!
! SAVE THE UPPER AND LOWER SUBSCRIPTS OF THE PORTION OF THE
!   ARRAY YET TO BE SORTED
!
      IF (L-I .LE. J-K) GO TO 9
      IL(M) = I
      IU(M) = L
      I = K
      M = M + 1
      GO TO 11
!
    9 IL(M) = K
      IU(M) = J
      J = L
      M = M + 1
      GO TO 11
!
! BEGIN AGAIN ON ANOTHER UNSORTED PORTION OF THE ARRAY
!
   10 M = M - 1
      IF (M .EQ. 0) RETURN
      I = IL(M)
      J = IU(M)
!
   11 IF (J-I .GE. 11) GO TO 4
      IF (I .EQ. 1) GO TO 2
      I = I - 1
!
! SORT ELEMENTS I+1,...,J.  NOTE THAT 1 .LE. I .LT. J AND
!   J-I .LT. 11.
!
   12 I = I + 1
      IF (I .EQ. J) GO TO 10
      INDX = IND(I+1)
      T = X(INDX)
      IT = INDX
      INDX = IND(I)
      IF (X(INDX) .LE. T) GO TO 12
      K = I
!
   13 IND(K+1) = IND(K)
      K = K - 1
      INDX = IND(K)
      IF (T .LT. X(INDX)) GO TO 13
      IND(K+1) = IT
      GO TO 12
      END
      SUBROUTINE REORDR (N,IFLAG, A,B,C,D, IND)
      INTEGER N, IFLAG, IND(N)
      DOUBLE PRECISION    A(N), B(N), C(N), D(N)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE USES A QUICK SORT TO REORDER THE DOUBLE PRECISION
! ARRAY A INTO INCREASING ORDER.  A RECORD OF THE PERMUTA-
! TIONS APPLIED TO A IS STORED IN IND, AND THESE PERMUTA-
! TIONS MAY BE APPLIED TO ONE, TWO, OR THREE ADDITIONAL
! VECTORS BY THIS ROUTINE.  ANY OTHER VECTOR V MAY BE PER-
! MUTED IN THE SAME FASHION BY CALLING SUBROUTINE PERMUT
! WITH N, IND, AND V AS PARAMETERS.
!   A SET OF NODES (X(I),Y(I),Z(I)) AND DATA VALUES W(I)
! MAY BE PREPROCESSED BY REORDR FOR INCREASED EFFICIENCY IN
! THE TRIANGULATION ROUTINE TRMESH.  THIS PREPROCESSING IS
! ALSO USEFUL FOR DETECTING DUPLICATE NODES.  NOTE THAT THE
! FIRST THREE NODES MUST NOT BE COLLINEAR ON INPUT TO
! TRMESH.  EITHER X, Y, OR Z MAY BE USED AS THE SORT KEY
! (ASSOCIATED WITH A).  IF THE NODAL COORDINATES ARE IN
! TERMS OF LATITUDE AND LONGITUDE, IT IS USUALLY ADVAN-
! TAGEOUS TO SORT ON LONGITUDE WITH THE SAME INTERCHANGES
! APPLIED TO LATITUDE AND THE DATA VALUES.  D WOULD BE A
! DUMMY PARAMETER IN THIS CASE.
!
! INPUT PARAMETERS -     N - NUMBER OF NODES.
!
!                    IFLAG - NUMBER OF VECTORS TO BE PER-
!                            MUTED.
!                            IFLAG .LE. 0 IF A, B, C, AND
!                              D ARE TO REMAIN UNALTERED.
!                            IFLAG = 1 IF ONLY A IS TO BE
!                              PERMUTED.
!                            IFLAG = 2 IF A AND B ARE TO BE
!                              PERMUTED.
!                            IFLAG = 3 IF A, B, AND C ARE
!                              TO BE PERMUTED.
!                            IFLAG .GE. 4 IF A, B, C, AND D
!                              ARE TO BE PERMUTED.
!
!                  A,B,C,D - VECTORS OF LENGTH N TO BE
!                            SORTED (ON THE COMPONENTS OF A)
!                            OR DUMMY PARAMETERS, DEPENDING
!                            ON IFLAG.
!
!                      IND - VECTOR OF LENGTH .GE. N.
!
! N, IFLAG, AND ANY DUMMY PARAMETERS ARE NOT ALTERED BY THIS
!   ROUTINE.
!
! OUTPUT PARAMETERS - A,B,C,D - SORTED OR UNALTERED VECTORS,
!                               DEPENDING ON IFLAG.
!
!                         IND - SEQUENCE OF INDICES 1,...,N
!                               PERMUTED IN THE SAME FASHION
!                               AS THE DOUBLE PRECISION VECTORS.  THUS
!                               THE ORDERING MAY BE APPLIED
!                               TO A VECTOR V AND STORED IN
!                               W BY SETTING W(I) =
!                               V(IND(I)), OR V MAY BE OVER-
!                               WRITTEN WITH THE ORDERING BY
!                               A CALL TO PERMUT.
!
! MODULES REFERENCED BY REORDR - QSORT, PERMUT
!
!***********************************************************
!
      INTEGER NN, NV
!
! LOCAL PARAMETERS -
!
! NN = LOCAL COPY OF N
! NV = LOCAL COPY OF IFLAG
!
      NN = N
      NV = IFLAG
      CALL QSORT (NN,A, IND)
      IF (NV .LE. 0) RETURN
      CALL PERMUT (NN,IND, A )
      IF (NV .EQ. 1) RETURN
      CALL PERMUT (NN,IND, B )
      IF (NV .EQ. 2) RETURN
      CALL PERMUT (NN,IND, C )
      IF (NV .EQ. 3) RETURN
      CALL PERMUT (NN,IND, D )
      RETURN
      END
      SUBROUTINE ROTATE (N,C,S, X,Y )
      INTEGER N
      DOUBLE PRECISION    C, S, X(N), Y(N)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!                                            ( C  S)
!   THIS ROUTINE APPLIES THE GIVENS ROTATION (     ) TO THE
!                                            (-S  C)
!               (X(1) ... X(N))
! 2 BY N MATRIX (             ).  THIS ROUTINE WAS TAKEN
!               (Y(1) ... Y(N))
! FROM LINPACK.
!
! INPUT PARAMETERS -   N - NUMBER OF COLUMNS TO BE ROTATED.
!
!                    C,S - ELEMENTS OF THE GIVENS ROTATION.
!                          THESE MAY BE DETERMINED BY
!                          SUBROUTINE GIVENS.
!
!                    X,Y - VECTORS OF LENGTH .GE. N
!                          CONTAINING THE 2-VECTORS TO BE
!                          ROTATED.
!
!   THE PARAMETERS N, C, AND S ARE NOT ALTERED BY THIS
! ROUTINE.
!
! OUTPUT PARAMETERS - X,Y - ROTATED VECTORS
!
! MODULES REFERENCED BY ROTATE - NONE
!
!***********************************************************
!
      INTEGER I
      DOUBLE PRECISION    XI, YI
!
! LOCAL PARAMETERS -
!
! I =     DO-LOOP INDEX
! XI,YI = X(I), Y(I)
!
      IF (N .LE. 0 .OR. (C .EQ. 1. .AND. S .EQ. 0.)) RETURN
      DO 1 I = 1,N
        XI = X(I)
        YI = Y(I)
        X(I) = C*XI + S*YI
    1   Y(I) = -S*XI + C*YI
      RETURN
      END
      SUBROUTINE SETUP (XI,YI,WI,WK,S1,S2,WT, ROW)
      DOUBLE PRECISION XI, YI, WI, WK, S1, S2, WT, ROW(6)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE SETS UP THE I-TH ROW OF AN AUGMENTED
! REGRESSION MATRIX FOR A WEIGHTED LEAST SQUARES FIT OF A
! QUADRATIC FUNCTION Q(X,Y) TO A SET OF DATA VALUES WI
! WHERE Q(0,0) = WK.  THE FIRST 3 COLUMNS (QUADRATIC TERMS)
! ARE SCALED BY 1/S2 AND THE FOURTH AND FIFTH COLUMNS (LIN-
! EAR TERMS) ARE SCALED BY 1/S1.
!
! INPUT PARAMETERS - XI,YI - COORDINATES OF NODE I.
!
!                       WI - DATA VALUE AT NODE I.
!
!                       WK - DATA VALUE INTERPOLATED BY Q AT
!                            THE ORIGIN.
!
!                    S1,S2 - INVERSE SCALE FACTORS.
!
!                       WT - WEIGHT FACTOR CORRESPONDING TO
!                            THE I-TH EQUATION.
!
!                      ROW - VECTOR OF LENGTH 6.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETER - ROW - VECTOR CONTAINING A ROW OF THE
!                          AUGMENTED REGRESSION MATRIX.
!
! MODULES REFERENCED BY SETUP - NONE
!
!***********************************************************
!
      DOUBLE PRECISION W1, W2
!
! LOCAL PARAMETERS -
!
! W1 = WEIGHTED SCALE FACTOR FOR THE LINEAR TERMS
! W2 = WEIGHTED SCALE FACTOR FOR THE QUADRATIC TERMS
!
      W1 = WT/S1
      W2 = WT/S2
      ROW(1) = XI*XI*W2
      ROW(2) = XI*YI*W2
      ROW(3) = YI*YI*W2
      ROW(4) = XI*W1
      ROW(5) = YI*W1
      ROW(6) = (WI-WK)*WT
      RETURN
      END
      SUBROUTINE SHIFTD (NFRST,NLAST,KK, IARR )
      INTEGER NFRST, NLAST, KK, IARR(1)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE SHIFTS A SET OF CONTIGUOUS ELEMENTS OF AN
! INTEGER ARRAY KK POSITIONS DOWNWARD (UPWARD IF KK .LT. 0).
! THE LOOPS ARE UNROLLED IN ORDER TO INCREASE EFFICIENCY.
!
! INPUT PARAMETERS - NFRST,NLAST - BOUNDS ON THE PORTION OF
!                                  IARR TO BE SHIFTED.  ALL
!                                  ELEMENTS BETWEEN AND
!                                  INCLUDING THE BOUNDS ARE
!                                  SHIFTED UNLESS NFRST .GT.
!                                  NLAST, IN WHICH CASE NO
!                                  SHIFT OCCURS.
!
!                             KK - NUMBER OF POSITIONS EACH
!                                  ELEMENT IS TO BE SHIFTED.
!                                  IF KK .LT. 0 SHIFT UP.
!                                  IF KK .GT. 0 SHIFT DOWN.
!
!                           IARR - INTEGER ARRAY OF LENGTH
!                                  .GE. NLAST + MAX(KK,0).
!
! NFRST, NLAST, AND KK ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETER -        IARR - SHIFTED ARRAY.
!
! MODULES REFERENCED BY SHIFTD - NONE
!
!***********************************************************
!
      INTEGER INC, K, NF, NL, NLP1, NS, NSL, I, IBAK, INDX, &
              IMAX
      DATA    INC/5/
!
! LOCAL PARAMETERS -
!
! INC =  DO-LOOP INCREMENT (UNROLLING FACTOR) -- IF INC IS
!          CHANGED, STATEMENTS MUST BE ADDED TO OR DELETED
!          FROM THE DO-LOOPS
! K =    LOCAL COPY OF KK
! NF =   LOCAL COPY OF NFRST
! NL =   LOCAL COPY OF NLAST
! NLP1 = NL + 1
! NS =   NUMBER OF SHIFTS
! NSL =  NUMBER OF SHIFTS DONE IN UNROLLED DO-LOOP (MULTIPLE
!          OF INC)
! I =    DO-LOOP INDEX AND INDEX FOR IARR
! IBAK = INDEX FOR DOWNWARD SHIFT OF IARR
! INDX = INDEX FOR IARR
! IMAX = BOUND ON DO-LOOP INDEX
!
      K = KK
      NF = NFRST
      NL = NLAST
      IF (NF .GT. NL  .OR.  K .EQ. 0) RETURN
      NLP1 = NL + 1
      NS = NLP1 - NF
      NSL = INC*(NS/INC)
      IF ( K .LT. 0) GO TO 4
!
! SHIFT DOWNWARD STARTING FROM THE BOTTOM
!
      IF (NSL .LE. 0) GO TO 2
      DO 1 I = 1,NSL,INC
        IBAK = NLP1 - I
        INDX = IBAK + K
        IARR(INDX) = IARR(IBAK)
        IARR(INDX-1) = IARR(IBAK-1)
        IARR(INDX-2) = IARR(IBAK-2)
        IARR(INDX-3) = IARR(IBAK-3)
        IARR(INDX-4) = IARR(IBAK-4)
    1   CONTINUE
!
! PERFORM THE REMAINING NS-NSL SHIFTS ONE AT A TIME
!
    2 IBAK = NLP1 - NSL
    3 IF (IBAK .LE. NF) RETURN
      IBAK = IBAK - 1
      INDX = IBAK + K
      IARR(INDX) = IARR(IBAK)
      GO TO 3
!
! SHIFT UPWARD STARTING FROM THE TOP
!
    4 IF (NSL .LE. 0) GO TO 6
      IMAX = NLP1 - INC
      DO 5 I = NF,IMAX,INC
        INDX = I + K
        IARR(INDX) = IARR(I)
        IARR(INDX+1) = IARR(I+1)
        IARR(INDX+2) = IARR(I+2)
        IARR(INDX+3) = IARR(I+3)
        IARR(INDX+4) = IARR(I+4)
    5   CONTINUE
!
! PERFORM THE REMAINING NS-NSL SHIFTS ONE AT A TIME
!
    6 I = NSL + NF
    7 IF (I .GT. NL) RETURN
      INDX = I + K
      IARR(INDX) = IARR(I)
      I = I + 1
      GO TO 7
      END
      SUBROUTINE SWAP (NIN1,NIN2,NOUT1,NOUT2, IADJ,IEND )
      INTEGER NIN1, NIN2, NOUT1, NOUT2, IADJ(1), IEND(1)
      EXTERNAL INDX
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE SWAPS THE DIAGONALS IN A CONVEX QUADRI-
! LATERAL.
!
! INPUT PARAMETERS -  NIN1,NIN2,NOUT1,NOUT2 - NODAL INDICES
!                            OF A PAIR OF ADJACENT TRIANGLES
!                            WHICH FORM A CONVEX QUADRILAT-
!                            ERAL.  NOUT1 AND NOUT2 ARE CON-
!                            NECTED BY AN ARC WHICH IS TO BE
!                            REPLACED BY THE ARC NIN1-NIN2.
!                            (NIN1,NOUT1,NOUT2) MUST BE TRI-
!                            ANGLE VERTICES IN COUNTERCLOCK-
!                            WISE ORDER.
!
! THE ABOVE PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
!                IADJ,IEND - TRIANGULATION DATA STRUCTURE
!                            (SEE SUBROUTINE TRMESH).
!
! OUTPUT PARAMETERS - IADJ,IEND - UPDATED WITH THE ARC
!                                 REPLACEMENT.
!
! MODULES REFERENCED BY SWAP - INDX, SHIFTD
!
!***********************************************************
!
      INTEGER IN(2), IO(2), IP1, IP2, J, K, NF, NL, I, &
              IMIN, IMAX
!
! LOCAL PARAMETERS -
!
! IN =        NIN1 AND NIN2 ORDERED BY INCREASING MAGNITUDE
!               (THE NEIGHBORS OF IN(1) PRECEDE THOSE OF
!               IN(2) IN IADJ)
! IO =        NOUT1 AND NOUT2 IN INCREASING ORDER
! IP1,IP2 =   PERMUTATION OF (1,2) SUCH THAT IO(IP1)
!               PRECEDES IO(IP2) AS A NEIGHBOR OF IN(1)
! J,K =       PERMUTATION OF (1,2) USED AS INDICES OF IN
!               AND IO
! NF,NL =     IADJ INDICES BOUNDARY A PORTION OF THE ARRAY
!               TO BE SHIFTED
! I =         IEND INDEX
! IMIN,IMAX = BOUNDS ON THE PORTION OF IEND TO BE INCRE-
!               MENTED OR DECREMENTED
!
      IN(1) = NIN1
      IN(2) = NIN2
      IO(1) = NOUT1
      IO(2) = NOUT2
      IP1 = 1
!
! ORDER THE INDICES SO THAT IN(1) .LT. IN(2) AND IO(1) .LT.
!   IO(2), AND CHOOSE IP1 AND IP2 SUCH THAT (IN(1),IO(IP1),
!   IO(IP2)) FORMS A TRIANGLE.
!
      IF (IN(1) .LT. IN(2)) GO TO 1
      IN(1) = IN(2)
      IN(2) = NIN1
      IP1 = 2
    1 IF (IO(1) .LT. IO(2)) GO TO 2
      IO(1) = IO(2)
      IO(2) = NOUT1
      IP1 = 3 - IP1
    2 IP2 = 3 - IP1
      IF (IO(2) .LT. IN(1)) GO TO 8
      IF (IN(2) .LT. IO(1)) GO TO 12
!
! IN(1) AND IO(1) PRECEDE IN(2) AND IO(2).  FOR (J,K) =
!   (1,2) AND (2,1), DELETE IO(K) AS A NEIGHBOR OF IO(J)
!   BY SHIFTING A PORTION OF IADJ EITHER UP OR DOWN AND
!   AND INSERT IN(K) AS A NEIGHBOR OF IN(J).
!
      DO 7 J = 1,2
        K = 3 - J
        IF (IN(J) .GT. IO(J)) GO TO 4
!
!   THE NEIGHBORS OF IN(J) PRECEDE THOSE OF IO(J) -- SHIFT
!     DOWN BY 1
!
        NF = 1 + INDX(IN(J),IO(IP1),IADJ,IEND)
        NL = -1 + INDX(IO(J),IO(K),IADJ,IEND)
        IF (NF .LE. NL) CALL SHIFTD(NF,NL,1, IADJ )
        IADJ(NF) = IN(K)
        IMIN = IN(J)
        IMAX = IO(J)-1
        DO 3 I = IMIN,IMAX
    3     IEND(I) = IEND(I) + 1
        GO TO 6
!
!   THE NEIGHBORS OF IO(J) PRECEDE THOSE OF IN(J) -- SHIFT
!     UP BY 1
!
    4   NF = 1 + INDX(IO(J),IO(K),IADJ,IEND)
        NL = -1 + INDX(IN(J),IO(IP2),IADJ,IEND)
        IF (NF .LE. NL) CALL SHIFTD(NF,NL,-1, IADJ )
        IADJ(NL) = IN(K)
        IMIN = IO(J)
        IMAX = IN(J) - 1
        DO 5 I = IMIN,IMAX
    5     IEND(I) = IEND(I) - 1
!
!   REVERSE (IP1,IP2) FOR (J,K) = (2,1)
!
    6   IP1 = IP2
        IP2 = 3 - IP1
    7   CONTINUE
      RETURN
!
! THE VERTICES ARE ORDERED (IO(1),IO(2),IN(1),IN(2)).
!   DELETE IO(2) BY SHIFTING UP BY 1
!
    8 NF = 1 + INDX(IO(1),IO(2),IADJ,IEND)
      NL = -1 + INDX(IO(2),IO(1),IADJ,IEND)
      IF (NF .LE. NL) CALL SHIFTD(NF,NL,-1, IADJ )
      IMIN = IO(1)
      IMAX = IO(2)-1
      DO 9 I = IMIN,IMAX
    9   IEND(I) = IEND(I) - 1
!
!   DELETE IO(1) BY SHIFTING UP BY 2 AND INSERT IN(2)
!
      NF = NL + 2
      NL = -1 + INDX(IN(1),IO(IP2),IADJ,IEND)
      IF (NF .LE. NL) CALL SHIFTD(NF,NL,-2, IADJ )
      IADJ(NL-1) = IN(2)
      IMIN = IO(2)
      IMAX = IN(1)-1
      DO 10 I = IMIN,IMAX
   10   IEND(I) = IEND(I) - 2
!
!   SHIFT UP BY 1 AND INSERT IN(1)
!
      NF = NL + 1
      NL = -1 + INDX(IN(2),IO(IP1),IADJ,IEND)
      CALL SHIFTD(NF,NL,-1, IADJ )
      IADJ(NL) = IN(1)
      IMIN = IN(1)
      IMAX = IN(2)-1
      DO 11 I = IMIN,IMAX
   11   IEND(I) = IEND(I) - 1
      RETURN
!
! THE VERTICES ARE ORDERED (IN(1),IN(2),IO(1),IO(2)).
!   DELETE IO(1) BY SHIFTING DOWN BY 1
!
   12 NF = 1 + INDX(IO(1),IO(2),IADJ,IEND)
      NL = -1 + INDX(IO(2),IO(1),IADJ,IEND)
      IF (NF .LE. NL) CALL SHIFTD(NF,NL,1, IADJ )
      IMIN = IO(1)
      IMAX = IO(2) - 1
      DO 13 I = IMIN,IMAX
   13   IEND(I) = IEND(I) + 1
!
!   DELETE IO(2) BY SHIFTING DOWN BY 2 AND INSERT IN(1)
!
      NL = NF - 2
      NF = 1 + INDX(IN(2),IO(IP2),IADJ,IEND)
      IF (NF .LE. NL) CALL SHIFTD(NF,NL,2, IADJ )
      IADJ(NF+1) = IN(1)
      IMIN = IN(2)
      IMAX = IO(1) - 1
      DO 14 I = IMIN,IMAX
   14   IEND(I) = IEND(I) + 2
!
!   SHIFT DOWN BY 1 AND INSERT IN(2)
!
      NL = NF - 1
      NF = 1 + INDX(IN(1),IO(IP1),IADJ,IEND)
      CALL SHIFTD(NF,NL,1, IADJ )
      IADJ(NF) = IN(2)
      IMIN = IN(1)
      IMAX = IN(2) - 1
      DO 15 I = IMIN,IMAX
   15   IEND(I) = IEND(I) + 1
      RETURN
      END
      LOGICAL FUNCTION SWPTST (N1,N2,N3,N4,X,Y,Z)
      INTEGER N1, N2, N3, N4
      DOUBLE PRECISION    X(1), Y(1), Z(1)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS FUNCTION DECIDES WHETHER OR NOT TO REPLACE A
! DIAGONAL ARC IN A QUADRILATERAL WITH THE OTHER DIAGONAL.
! THE DECISION WILL BE POSITIVE (SWPTST = .TRUE.) IF AND
! ONLY IF N4 LIES ABOVE THE PLANE (IN THE HALF-SPACE NOT
! CONTAINING THE ORIGIN) DEFINED BY (N1,N2,N3), OR EQUIV-
! ALENTLY, IF THE PROJECTION OF N4 ONTO THIS PLANE IS
! INTERIOR TO THE CIRCUMCIRCLE OF (N1,N2,N3).  THE DECISION
! WILL BE NEGATIVE IF THE QUADRILATERAL IS NOT STRICTLY
! CONVEX.
!
! INPUT PARAMETERS - N1,N2,N3,N4 - INDICES OF THE FOUR NODES
!                      DEFINING THE QUADRILATERAL.  N1 AND
!                      N2 ARE CURRENTLY CONNECTED BY A DIAG-
!                      ONAL ARC WHICH SHOULD BE REPLACED BY
!                      AN ARC CONNECTING N3 TO N4 IF THE
!                      DECISION IS MADE TO SWAP.  N1, N2, N3
!                      MUST BE IN COUNTERCLOCKWISE ORDER.
!
!              X,Y,Z - VECTORS OF NODAL COORDINATES.  (X(I),
!                      Y(I),Z(I)) ARE THE COORDINATES OF
!                      NODE I FOR I = N1, N2, N3, AND N4.
!
! NONE OF THE INPUT PARAMETERS ARE ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETER -  SWPTST - .TRUE. IFF THE ARC CONNECTING
!                              N1 AND N2 IS TO BE REPLACED.
!
! MODULES REFERENCED BY SWPTST - NONE
!
!***********************************************************
!
      DOUBLE PRECISION X4, Y4, Z4, DX1, DX2, DX3, DY1, DY2, DY3, &
           DZ1, DZ2, DZ3
!
! LOCAL PARAMETERS -
!
! X4,Y4,Z4 =    COORDINATES OF N4
! DX1,DY1,DZ1 = COORDINATES OF N1 - N4
! DX2,DY2,DZ2 = COORDINATES OF N2 - N4
! DX3,DY3,DZ3 = COORDINATES OF N3 - N4
!
      X4 = X(N4)
      Y4 = Y(N4)
      Z4 = Z(N4)
      DX1 = X(N1) - X4
      DX2 = X(N2) - X4
      DX3 = X(N3) - X4
      DY1 = Y(N1) - Y4
      DY2 = Y(N2) - Y4
      DY3 = Y(N3) - Y4
      DZ1 = Z(N1) - Z4
      DZ2 = Z(N2) - Z4
      DZ3 = Z(N3) - Z4
!
! N4 LIES ABOVE THE PLANE OF (N1,N2,N3) IFF N3 LIES ABOVE
!   THE PLANE OF (N2,N1,N4) IFF DET(N3-N4,N2-N4,N1-N4) =
!   (N3-N4,N2-N4 X N1-N4) .GT. 0.
!
      SWPTST = DX3*(DY2*DZ1 - DY1*DZ2) &
              -DY3*(DX2*DZ1 - DX1*DZ2) &
              +DZ3*(DX2*DY1 - DX1*DY2) .GT. 0.
      RETURN
      END
      SUBROUTINE TRANS (N,RLAT,RLON, X,Y,Z)
      INTEGER N
      DOUBLE PRECISION    RLAT(N), RLON(N), X(N), Y(N), Z(N)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE TRANSFORMS SPHERICAL COORDINATES INTO
! CARTESIAN COORDINATES ON THE UNIT SPHERE FOR INPUT TO
! SUBROUTINE TRMESH.  STORAGE FOR X AND Y MAY COINCIDE WITH
! STORAGE FOR RLAT AND RLON IF THE LATTER NEED NOT BE SAVED.
!
! INPUT PARAMETERS -    N - NUMBER OF NODES (POINTS ON THE
!                           UNIT SPHERE) WHOSE COORDINATES
!                           ARE TO BE TRANSFORMED.
!
!                    RLAT - N-VECTOR CONTAINING LATITUDINAL
!                           COORDINATES OF THE NODES IN
!                           RADIANS.
!
!                    RLON - N-VECTOR CONTAINING LONGITUDINAL
!                           COORDINATES OF THE NODES IN
!                           RADIANS.
!
!                   X,Y,Z - VECTORS OF LENGTH .GE. N.
!
! N, RLAT, AND RLON ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS - X,Y,Z - CARTESIAN COORDINATES IN THE
!                             RANGE (-1,1).  X(I)**2 +
!                             Y(I)**2 + Z(I)**2 = 1 FOR
!                             I = 1,...,N.
!
! MODULES REFERENCED BY TRANS - NONE
!
! INTRINSIC FUNCTIONS CALLED BY TRANS - COS, SIN
!
!***********************************************************
!
      INTEGER NN, I
      DOUBLE PRECISION    PHI, THETA, COSPHI
!
! LOCAL PARAMETERS -
!
! NN =     LOCAL COPY OF N
! I =      DO-LOOP INDEX
! PHI =    LATITUDE
! THETA =  LONGITUDE
! COSPHI = COS(PHI)
!
      NN = N
      DO 1 I = 1,NN
        PHI = RLAT(I)
        THETA = RLON(I)
        COSPHI = COS(PHI)
        X(I) = COSPHI*COS(THETA)
        Y(I) = COSPHI*SIN(THETA)
    1   Z(I) = SIN(PHI)
      RETURN
      END
      SUBROUTINE TRFIND (NST,P,X,Y,Z,IADJ,IEND, B1,B2,B3, &
                         I1,I2,I3)
      INTEGER NST, IADJ(1), IEND(1), I1, I2, I3
      DOUBLE PRECISION    P(3), X(1), Y(1), Z(1), B1, B2, B3
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE LOCATES A POINT P RELATIVE TO A TRIANGU-
! LATION OF THE CONVEX HULL OF A SET OF NODES (POINTS ON THE
! UNIT SPHERE).
!
! INPUT PARAMETERS -    NST - INDEX OF NODE AT WHICH TRFIND
!                             BEGINS SEARCH.  SEARCH TIME
!                             DEPENDS ON THE PROXIMITY OF
!                             NST TO P.
!
!                         P - X-, Y-, AND Z-COORDINATES (IN
!                             THAT ORDER) OF THE POINT TO BE
!                             LOCATED.
!
!                     X,Y,Z - VECTORS CONTAINING CARTESIAN
!                             COORDINATES OF THE NODES IN
!                             THE MESH.  (X(I),Y(I),Z(I))
!                             DEFINES NODE I FOR I = 1,...,N
!                             WHERE N .GE. 3.
!
!                      IADJ - SET OF ADJACENCY LISTS OF
!                             NODES IN THE MESH.
!
!                      IEND - POINTERS TO THE ENDS OF
!                             ADJACENCY LISTS IN IADJ FOR
!                             EACH NODE IN THE MESH.
!
! IADJ AND IEND MAY BE CREATED BY TRMESH.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS - B1,B2,B3 - UNNORMALIZED BARYCENTRIC
!                                COORDINATES OF THE CENTRAL
!                                PROJECTION OF P ONTO THE
!                                UNDERLYING PLANAR TRIANGLE
!                                IF P IS IN THE CONVEX HULL
!                                OF THE NODES, UNCHANGED IF
!                                I1 = 0.
!
!                     I1,I2,I3 - COUNTERCLOCKWISE-ORDERED
!                                VERTEX INDICES OF A TRI-
!                                ANGLE CONTAINING P IF P IS
!                                AN INTERIOR POINT.  IF P IS
!                                NOT CONTAINED IN THE CONVEX
!                                HULL OF THE NODES, I1 AND
!                                I2 ARE THE RIGHTMOST AND
!                                LEFTMOST NODES WHICH ARE
!                                VISIBLE FROM P, AND I3 = 0.
!                                IF ALL BOUNDARY NODES ARE
!                                VISIBLE FROM P, I1 AND I2
!                                COINCIDE.  IF P AND ALL OF
!                                THE NODES LIE ON A SINGLE
!                                GREAT CIRCLE THEN I1 = I2
!                                = I3 = 0.
!
! MODULES REFERENCED BY TRFIND - NONE
!
! INTRINSIC FUNCTIONS CALLED BY TRFIND - MAX0, ABS
!
!***********************************************************
!
      INTEGER N0, N1, N2, N3, N4, INDX, IND, NF, NL, N1S, &
              N2S, NEXT, NRST, NRMAX, LUN
      DOUBLE PRECISION    XP, YP, ZP, B3P1, S12, PTN1, PTN2, Q(3), DET
      DOUBLE PRECISION    X1,Y1,Z1,X2,Y2,Z2,X0,Y0,Z0
      DATA    NRMAX/5/,  LUN/6/
!
! LOCAL PARAMETERS -
!
! N0,N1,N2 = NODES IN COUNTERCLOCKWISE ORDER DEFINING A
!              CONE (WITH VERTEX N0) CONTAINING P, OR END-
!              POINTS OF A BOUNDARY EDGE SUCH THAT P RIGHT
!              N1->N2
! N3,N4 =    NODES OPPOSITE N1->N2 AND N2->N1, RESPECTIVELY
! INDX,IND = INDICES FOR IADJ
! NF,NL =    FIRST AND LAST NEIGHBORS OF N0 IN IADJ, OR
!              FIRST (RIGHTMOST) AND LAST (LEFTMOST) NODES
!              VISIBLE FROM P WHEN P IS EXTERIOR TO THE
!              BOUNDARY
! N1S,N2S =  INITIALLY-DETERMINED VALUES OF N1 AND N2 WHEN
!              P IS EXTERIOR
! NEXT =     CANDIDATE FOR I1 OR I2 WHEN P IS EXTERIOR
! NRST =     NUMBER OF RESTARTS WITH NEW N0
! NRMAX =    MAXIMUM ALLOWABLE NUMBER OF RESTARTS BEFORE
!              TERMINATING WITH AN ERROR MESSAGE
! LUN =      LOGICAL UNIT FOR ERROR MESSAGES
! XP,YP,ZP = LOCAL VARIABLES CONTAINING P(1), P(2), AND P(3)
! B3P1 =     B3 + 1 -- B3P1 = 1 IFF B3 = 0 TO WITHIN THE
!              MACHINE PRECISION
! S12 =      SCALAR PRODUCT <N1,N2>
! PTN1 =     SCALAR PRODUCT <P,N1>
! PTN2 =     SCALAR PRODUCT <P,N2>
! Q =        (N2 X N1) X N2  OR  N1 X (N2 X N1) -- USED IN
!              THE BOUNDARY TRAVERSAL WHEN P IS EXTERIOR
! DET =      STATEMENT FUNCTION WHICH COMPUTES A DETERMI-
!              NANT -- DET(X1,...,Z0) .GE. 0 IFF (X0,Y0,Z0)
!              IS IN THE (CLOSED) LEFT HEMISPHERE DEFINED BY
!              THE PLANE CONTAINING (0,0,0), (X1,Y1,Z1),
!              AND (X2,Y2,Z2) WHERE LEFT IS DEFINED RELA-
!              TIVE TO AN OBSERVER AT (X1,Y1,Z1) FACING
!              (X2,Y2,Z2).
      DET (X1,Y1,Z1,X2,Y2,Z2,X0,Y0,Z0) = X0*(Y1*Z2-Y2*Z1) &
           - Y0*(X1*Z2-X2*Z1) + Z0*(X1*Y2-X2*Y1)
!
! INITIALIZE VARIABLES
!
      XP = P(1)
      YP = P(2)
      ZP = P(3)
      NRST = 0
      N0 = MAX0(NST,1)
!
! FIND A CONE WITH VERTEX N0 CONTAINING P
!
    1 INDX = IEND(N0)
      NL = IADJ(INDX)
      INDX = 1
      IF (N0 .NE. 1) INDX = IEND(N0-1) + 1
      NF = IADJ(INDX)
      N1 = NF
      IF (NL .NE. 0) GO TO 3
!
! N0 IS A BOUNDARY NODE.  SET NL TO THE LAST NONZERO
!   NEIGHBOR OF N0.
!
      IND = IEND(N0) - 1
      NL = IADJ(IND)
      IF ( DET(X(N0),Y(N0),Z(N0),X(NF),Y(NF),Z(NF), &
               XP,YP,ZP) .GE. 0. ) GO TO 2
!
! P IS TO THE RIGHT OF THE BOUNDARY EDGE N0->NF
!
      N1 = N0
      N2 = NF
      GO TO 16
    2 IF ( DET(X(NL),Y(NL),Z(NL),X(N0),Y(N0),Z(N0), &
               XP,YP,ZP) .GE. 0. ) GO TO 4
!
! P IS TO THE RIGHT OF THE BOUNDARY EDGE NL->N0
!
      N1 = NL
      N2 = N0
      GO TO 16
!
! N0 IS AN INTERIOR NODE.  FIND N1.
!
    3 IF ( DET(X(N0),Y(N0),Z(N0),X(N1),Y(N1),Z(N1), &
               XP,YP,ZP) .GE. 0. ) GO TO 4
      INDX = INDX + 1
      N1 = IADJ(INDX)
      IF (N1 .EQ. NL) GO TO 7
      GO TO 3
!
! P IS TO THE LEFT OF ARCS N0->N1 AND NL->N0.  SET N2 TO THE
!   NEXT NEIGHBOR OF N0 (FOLLOWING N1).
!
    4 INDX = INDX + 1
      N2 = IADJ(INDX)
      IF ( DET(X(N0),Y(N0),Z(N0),X(N2),Y(N2),Z(N2), &
               XP,YP,ZP) .LT. 0. ) GO TO 8
      N1 = N2
      IF (N1 .NE. NL) GO TO 4
      IF ( DET(X(N0),Y(N0),Z(N0),X(NF),Y(NF),Z(NF), &
               XP,YP,ZP) .LT. 0. ) GO TO 7
!
! P IS LEFT OF OR ON ARCS N0->NB FOR ALL NEIGHBORS NB
!   OF N0.  TEST FOR P = +/-N0.
!
      IF (ABS(X(N0)*XP + Y(N0)*YP + Z(N0)*ZP) .GE. 1.) &
          GO TO 6
!
! ALL POINTS ARE COLLINEAR IFF P IS LEFT OF NB->N0 FOR
!   ALL NEIGHBORS NB OF N0.  SEARCH THE NEIGHBORS OF N0
!   IN REVERSE ORDER.  NOTE -- N1 = NL AND INDX POINTS TO
!   NL.
!
    5 IF ( DET(X(N1),Y(N1),Z(N1),X(N0),Y(N0),Z(N0), &
               XP,YP,ZP) .LT. 0. ) GO TO 6
      IF (N1 .EQ. NF) GO TO 24
      INDX = INDX - 1
      N1 = IADJ(INDX)
      GO TO 5
!
! P IS TO THE RIGHT OF N1->N0, OR P = +/-N0.  SET N0 TO N1
!   AND START OVER.
!
    6 N0 = N1
      NRST = NRST + 1
      IF (NRST .EQ. NRMAX) GO TO 25
      GO TO 1
!
! P IS BETWEEN ARCS N0->N1 AND N0->NF
!
    7 N2 = NF
!
! P IS CONTAINED IN A CONE DEFINED BY LINE SEGMENTS N0-N1
!   AND N0-N2 WHERE N1 IS ADJACENT TO N2
!
    8 N3 = N0
    9 B3 = DET(X(N1),Y(N1),Z(N1),X(N2),Y(N2),Z(N2),XP,YP,ZP)
      IF (B3 .GE. 0.) GO TO 13
!
! SET N4 TO THE FIRST NEIGHBOR OF N2 FOLLOWING N1
!
      INDX = IEND(N2)
      IF (IADJ(INDX) .NE. N1) GO TO 10
!
! N1 IS THE LAST NEIGHBOR OF N2.
! SET N4 TO THE FIRST NEIGHBOR.
!
      INDX = 1
      IF (N2 .NE. 1) INDX = IEND(N2-1) + 1
      N4 = IADJ(INDX)
      GO TO 11
!
! N1 IS NOT THE LAST NEIGHBOR OF N2
!
   10 INDX = INDX-1
      IF (IADJ(INDX) .NE. N1) GO TO 10
      N4 = IADJ(INDX+1)
      IF (N4 .EQ. 0) GO TO 16
!
! DEFINE A NEW ARC N1->N2 WHICH INTERSECTS THE LINE
!   SEGMENT N0-P
!
   11 IF ( DET(X(N0),Y(N0),Z(N0),X(N4),Y(N4),Z(N4), &
               XP,YP,ZP) .GE. 0. ) GO TO 12
      N3 = N2
      N2 = N4
      GO TO 9
   12 N3 = N1
      N1 = N4
      GO TO 9
!
! P IS IN (N1,N2,N3) UNLESS N0, N1, N2, AND P ARE COLLINEAR
!   OR P IS CLOSE TO -N0.
!
   13 B3P1 = B3 + 1.
      IF (B3P1 .LE. 1.) GO TO 14
!
! B3 .NE. 0.
!
      B1 = DET(X(N2),Y(N2),Z(N2),X(N3),Y(N3),Z(N3),XP,YP,ZP)
      B2 = DET(X(N3),Y(N3),Z(N3),X(N1),Y(N1),Z(N1),XP,YP,ZP)
      IF (B1 .GE. 0.  .AND.  B2 .GE. 0.) GO TO 15
!
! RESTART WITH N0 = N3
!
      NRST = NRST + 1
      IF (NRST .EQ. NRMAX) GO TO 25
      N0 = N3
      GO TO 1
!
! B3 = 0 AND THUS P LIES ON N1->N2. COMPUTE
!   B1 = DET(P,N2 X N1,N2) AND B2 = DET(P,N1,N2 X N1).
!
   14 B3 = 0.
      S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
      PTN1 = XP*X(N1) + YP*Y(N1) + ZP*Z(N1)
      PTN2 = XP*X(N2) + YP*Y(N2) + ZP*Z(N2)
      B1 = PTN1 - S12*PTN2
      B2 = PTN2 - S12*PTN1
      IF (B1 .GE. 0.  .AND.  B2 .GE. 0.) GO TO 15
!
! RESTART WITH N0 = 1
!
      NRST = NRST + 1
      IF (NRST .EQ. NRMAX) GO TO 25
      N0 = 1
      GO TO 1
!
! P IS IN (N1,N2,N3)
!
   15 I1 = N1
      I2 = N2
      I3 = N3
      RETURN
!
! P RIGHT N1->N2 WHERE N1->N2 IS A BOUNDARY EDGE.
!   SAVE N1 AND N2, AND SET NL = 0 TO INDICATE THAT
!   NL HAS NOT YET BEEN FOUND.
!
   16 N1S = N1
      N2S = N2
      NL = 0
!
!           COUNTERCLOCKWISE BOUNDARY TRAVERSAL
!
   17 INDX = 1
      IF (N2 .NE. 1) INDX = IEND(N2-1) + 1
      NEXT = IADJ(INDX)
      IF ( DET(X(N2),Y(N2),Z(N2),X(NEXT),Y(NEXT),Z(NEXT), &
               XP,YP,ZP) .LT. 0. ) GO TO 18
!
! N2 IS THE RIGHTMOST VISIBLE NODE IF P FORWARD N2->N1
!   OR NEXT FORWARD N2->N1 -- SET Q TO (N2 X N1) X N2
!
      S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
      Q(1) = X(N1) - S12*X(N2)
      Q(2) = Y(N1) - S12*Y(N2)
      Q(3) = Z(N1) - S12*Z(N2)
      IF (XP*Q(1) + YP*Q(2) + ZP*Q(3) .GE. 0.) GO TO 19
      IF (X(NEXT)*Q(1) + Y(NEXT)*Q(2) + Z(NEXT)*Q(3) &
          .GE. 0.) GO TO 19
!
! N1, N2, NEXT, AND P ARE NEARLY COLLINEAR, AND N2 IS
!   THE LEFTMOST VISIBLE NODE
!
      NL = N2
   18 N1 = N2
      N2 = NEXT
      IF (N2 .NE. N1S) GO TO 17
!
! ALL BOUNDARY NODES ARE VISIBLE
!
      I1 = N1S
      I2 = N1S
      I3 = 0
      RETURN
!
! N2 IS THE RIGHTMOST VISIBLE NODE
!
   19 NF = N2
      IF (NL .NE. 0) GO TO 23
!
! RESTORE INITIAL VALUES OF N1 AND N2, AND BEGIN SEARCH
!   FOR THE LEFTMOST VISIBLE NODE
!
      N2 = N2S
      N1 = N1S
!
!           CLOCKWISE BOUNDARY TRAVERSAL
!
   20 INDX = IEND(N1) - 1
      NEXT = IADJ(INDX)
      IF ( DET(X(NEXT),Y(NEXT),Z(NEXT),X(N1),Y(N1),Z(N1), &
               XP,YP,ZP) .LT. 0. ) GO TO 21
!
! N1 IS THE LEFTMOST VISIBLE NODE IF P OR NEXT IS
!   FORWARD OF N1->N2 -- COMPUTE Q = N1 X (N2 X N1)
!
      S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
      Q(1) = X(N2) - S12*X(N1)
      Q(2) = Y(N2) - S12*Y(N1)
      Q(3) = Z(N2) - S12*Z(N1)
      IF (XP*Q(1) + YP*Q(2) + ZP*Q(3) .GE. 0.) GO TO 22
      IF (X(NEXT)*Q(1) + Y(NEXT)*Q(2) + Z(NEXT)*Q(3) &
          .GE. 0.) GO TO 22
!
! P, NEXT, N1, AND N2 ARE NEARLY COLLINEAR AND N1 IS THE
!   RIGHTMOST VISIBLE NODE
!
      NF = N1
   21 N2 = N1
      N1 = NEXT
      GO TO 20
!
! N1 IS THE LEFTMOST VISIBLE NODE
!
   22 NL = N1
!
! NF AND NL HAVE BEEN FOUND
!
   23 I1 = NF
      I2 = NL
      I3 = 0
      RETURN
!
! ALL POINTS ARE COLLINEAR
!
   24 I1 = 0
      I2 = 0
      I3 = 0
      RETURN
!
! MORE THAN NRMAX RESTARTS.  PRINT AN ERROR MESSAGE AND
!   TERMINATE PROCESSING.
!
   25 WRITE (LUN,100) NST, N0, N1, N2, N3, N4
  100 FORMAT (1H1,37HERROR IN TRFIND -- POSSIBLE INFINITE , &
              23HLOOP DUE TO F.P. ERROR.//1H , &
              21HNST,N0,N1,N2,N3,N4 = ,5(I5,2H, ),I5//1H , &
              40HTRY REORDERING THE NODES OR CHANGING NST)
      STOP
      END
      SUBROUTINE TRMESH (N,X,Y,Z, IADJ,IEND,IER)
      INTEGER N, IADJ(6*(N-1)), IEND(N), IER
      DOUBLE PRECISION    X(N), Y(N), Z(N)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE CREATES A THIESSEN TRIANGULATION OF THE
! CONVEX HULL OF N ARBITRARILY DISTRIBUTED POINTS (NODES) ON
! THE UNIT SPHERE.  IF THE NODES ARE NOT CONTAINED IN A
! SINGLE HEMISPHERE, THEIR CONVEX HULL IS THE ENTIRE SPHERE,
! AND THERE ARE NO BOUNDARY NODES.  THE TRIANGULATION IS
! OPTIMAL IN THE SENSE THAT IT IS AS NEARLY EQUIANGULAR AS
! POSSIBLE.  TRMESH IS PART OF AN INTERPOLATION PACKAGE
! WHICH ALSO PROVIDES SUBROUTINES TO REORDER THE NODES,
! TRANSFORM COORDINATES, ADD A NEW NODE, DELETE AN ARC, PLOT
! THE MESH, AND PRINT THE DATA STRUCTURE.
!   UNLESS THE NODES ARE ALREADY ORDERED IN SOME REASONABLE
! FASHION, THEY SHOULD BE REORDERED BY SUBROUTINE REORDR FOR
! INCREASED EFFICIENCY BEFORE CALLING TRMESH.  SPHERICAL
! COORDINATES (LATITUDE AND LONGITUDE) MAY BE CONVERTED TO
! CARTESIAN COORDINATES BY SUBROUTINE TRANS.
!
! INPUT PARAMETERS -     N - NUMBER OF NODES IN THE MESH.
!                            N .GE. 3.
!
!                    X,Y,Z - N-VECTORS CONTAINING CARTESIAN
!                            COORDINATES OF DISTINCT NODES.
!                            (X(I),Y(I),Z(I)) DEFINES NODE
!                            I.  X(I)**2 + Y(I)**2 + Z(I)**2
!                            = 1 FOR ALL I.  THE FIRST THREE
!                            NODES MUST NOT BE COLLINEAR
!                            (LIE ON A SINGLE GREAT CIRCLE).
!
!                     IADJ - VECTOR OF LENGTH .GE. 6*(N-1).
!
!                     IEND - VECTOR OF LENGTH .GE. N.
!
! N, X, Y, AND Z ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS - IADJ - ADJACENCY LISTS OF NEIGHBORS IN
!                            COUNTERCLOCKWISE ORDER.  THE
!                            LIST FOR NODE I+1 FOLLOWS THAT
!                            FOR NODE I. THE VALUE 0 DENOTES
!                            THE BOUNDARY (OR A PSEUDO-NODE
!                            FROM WHICH ALL BOUNDARY NODES
!                            ARE VISIBLE) AND IS ALWAYS THE
!                            LAST NEIGHBOR OF A BOUNDARY
!                            NODE.  IADJ IS UNCHANGED IF IER
!                            .NE. 0.
!
!                     IEND - POINTERS TO THE ENDS OF
!                            ADJACENCY LISTS (SETS OF
!                            NEIGHBORS) IN IADJ.  THE
!                            NEIGHBORS OF NODE 1 BEGIN IN
!                            IADJ(1).  FOR K .GT. 1, THE
!                            NEIGHBORS OF NODE K BEGIN IN
!                            IADJ(IEND(K-1)+1) AND K HAS
!                            IEND(K) - IEND(K-1) NEIGHBORS
!                            INCLUDING (POSSIBLY) THE
!                            BOUNDARY.  IADJ(IEND(K)) .EQ. 0
!                            IFF NODE K IS ON THE BOUNDARY.
!                            IEND IS UNCHANGED IF IER .NE.
!                            0.
!
!                      IER - ERROR INDICATOR
!                            IER = 0 IF NO ERRORS WERE
!                                    ENCOUNTERED.
!                            IER = 1 IF N .LT. 3.
!                            IER = 2 IF THE FIRST THREE
!                                    NODES ARE COLLINEAR.
!
! MODULES REFERENCED BY TRMESH - ADNODE, TRFIND, INTADD,
!                                BDYADD, COVSPH, SHIFTD,
!                                INDEX, SWPTST, SWAP
!
!***********************************************************
!
      INTEGER NN, K, IERR
      LOGICAL LEFT
      DOUBLE PRECISION X1,Y1,Z1,X2,Y2,Z2,X0,Y0,Z0
!
! LOCAL PARAMETERS -
!
! NN =   LOCAL COPY OF N
! K =    INDEX OF NODE TO BE ADDED, DO-LOOP INDEX
! IERR = ERROR FLAG FOR CALL TO ADNODE - NOT CHECKED
! LEFT = STATEMENT FUNCTION -- TRUE IFF (X0,Y0,Z0) IS IN
!          THE (CLOSED) LEFT HEMISPHERE DEFINED BY THE
!          PLANE CONTAINING (0,0,0), (X1,Y1,Z1), AND
!          (X2,Y2,Z2) WHERE LEFT IS DEFINED RELATIVE TO
!          AN OBSERVER AT (X1,Y1,Z1) FACING (X2,Y2,Z2)
!
      LEFT (X1,Y1,Z1,X2,Y2,Z2,X0,Y0,Z0) = X0*(Y1*Z2-Y2*Z1) &
             - Y0*(X1*Z2-X2*Z1) + Z0*(X1*Y2-X2*Y1) .GE. 0.
      NN = N
      IER = 1
      IF (NN .LT. 3) RETURN
      IER = 0
      IF ( .NOT. LEFT (X(1),Y(1),Z(1),X(2),Y(2),Z(2), &
                       X(3),Y(3),Z(3)) ) GO TO 2
      IF ( .NOT. LEFT (X(2),Y(2),Z(2),X(1),Y(1),Z(1), &
                       X(3),Y(3),Z(3)) ) GO TO 1
!
! FIRST 3 NODES ARE COLLINEAR
!
      IER = 2
      RETURN
!
! FIRST TRIANGLE IS (1,2,3), I.E. 3 STRICTLY LEFT 1-2,
!   I.E. 3 LIES IN THE LEFT HEMISPHERE DEFINED BY ARC 1-2
!
    1 IADJ(1) = 2
      IADJ(2) = 3
      IADJ(3) = 0
      IADJ(4) = 3
      IADJ(5) = 1
      IADJ(6) = 0
      IADJ(7) = 1
      IADJ(8) = 2
      IADJ(9) = 0
      GO TO 3
!
! FIRST TRIANGLE IS (3,2,1) = (2,1,3) = (1,3,2).  INITIALIZE
!   IADJ.
!
    2 IADJ(1) = 3
      IADJ(2) = 2
      IADJ(3) = 0
      IADJ(4) = 1
      IADJ(5) = 3
      IADJ(6) = 0
      IADJ(7) = 2
      IADJ(8) = 1
      IADJ(9) = 0
!
! INITIALIZE IEND
!
    3 IEND(1) = 3
      IEND(2) = 6
      IEND(3) = 9
      IF (NN .EQ. 3) RETURN
!
! ADD NODES 4,...,N TO THE TRIANGULATION
!
      DO 4 K = 4,NN
        !print *,'node ',k,' out of',nn
        CALL ADNODE (K,X,Y,Z, IADJ,IEND, IERR)
    4   CONTINUE
      RETURN
      END
      SUBROUTINE GRADG (N,X,Y,Z,W,IADJ,IEND,EPS, NIT, &
                        GRAD, IER)
      INTEGER N, IADJ(1), IEND(N), NIT, IER
      DOUBLE PRECISION    X(N), Y(N), Z(N), W(N), EPS, GRAD(3,N)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A TRIANGULATION OF N NODES ON THE UNIT SPHERE WITH
! ASSOCIATED DATA VALUES, THIS ROUTINE USES A GLOBAL METHOD
! TO COMPUTE ESTIMATED GRADIENTS AT THE NODES.  THE METHOD
! CONSISTS OF MINIMIZING A QUADRATIC FUNCTIONAL Q(G) OVER
! THE N-VECTOR G OF GRADIENTS WHERE Q APPROXIMATES THE LIN-
! EARIZED CURVATURE OF AN INTERPOLANT F OVER THE TRIANGULA-
! TION.  THE RESTRICTION OF F TO AN ARC OF THE TRIANGULATION
! IS TAKEN TO BE THE HERMITE CUBIC (WITH RESPECT TO ARC-
! LENGTH) INTERPOLANT OF THE DATA VALUES AND TANGENTIAL
! GRADIENT COMPONENTS AT THE ENDPOINTS OF THE ARC, AND Q IS
! THE SUM OF THE LINEARIZED CURVATURES OF F ALONG THE ARCS
! -- THE INTEGRALS OVER THE ARCS OF D2F(A)**2 WHERE D2F(A)
! IS THE SECOND DERIVATIVE OF F WITH RESPECT TO ARC-LENGTH
! A.
!   SINCE THE GRADIENT AT NODE K LIES IN THE PLANE TANGENT
! TO THE SPHERE SURFACE AT K, IT IS DEFINED BY TWO COMPO-
! NENTS -- ITS X AND Y COMPONENTS IN THE COORDINATE SYSTEM
! OBTAINED BY ROTATING K TO THE NORTH POLE.  THUS, THE MIN-
! IMIZATION PROBLEM CORRESPONDS TO AN ORDER 2N SYMMETRIC
! POSITIVE-DEFINITE SPARSE LINEAR SYSTEM WHICH IS SOLVED BY
! THE BLOCK GAUSS-SEIDEL METHOD WITH 2 BY 2 BLOCKS.
!   AN ALTERNATIVE METHOD, SUBROUTINE GRADL, COMPUTES A
! LOCAL APPROXIMATION TO THE GRADIENT AT A SINGLE NODE AND,
! WHILE SLIGHTLY LESS EFFICIENT, WAS FOUND TO BE GENERALLY
! MORE ACCURATE WHEN THE NODAL DISTRIBUTION IS VERY DENSE,
! VARIES GREATLY, OR DOES NOT COVER THE SPHERE.  GRADG, ON
! THE OTHER HAND, WAS FOUND TO BE SLIGHTLY MORE ACCURATE ON
! A SOMEWHAT UNIFORM DISTRIBUTION OF 514 NODES.
!
! INPUT PARAMETERS - N - NUMBER OF NODES.  N .GE. 3.
!
!                X,Y,Z - CARTESIAN COORDINATES OF THE NODES.
!                        X(I)**2 + Y(I)**2 + Z(I)**2 = 1 FOR
!                        I = 1,...,N.
!
!                    W - DATA VALUES AT THE NODES.  W(I) IS
!                        ASSOCIATED WITH (X(I),Y(I),Z(I)).
!
!            IADJ,IEND - DATA STRUCTURE DEFINING THE TRIAN-
!                        GULATION.  SEE SUBROUTINE TRMESH.
!
!                  EPS - NONNEGATIVE CONVERGENCE CRITERION.
!                        THE METHOD IS TERMINATED WHEN THE
!                        MAXIMUM CHANGE IN A GRADIENT COMPO-
!                        NENT BETWEEN ITERATIONS IS AT MOST
!                        EPS.  EPS = 1.E-3 IS SUFFICIENT FOR
!                        EFFECTIVE CONVERGENCE.
!
! THE ABOVE PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
!                  NIT - MAXIMUM NUMBER OF GAUSS-SEIDEL
!                        ITERATIONS TO BE APPLIED.  THIS
!                        MAXIMUM WILL LIKELY BE ACHIEVED IF
!                        EPS IS SMALLER THAN THE MACHINE
!                        PRECISION.  OPTIMAL EFFICIENCY WAS
!                        ACHIEVED IN TESTING WITH EPS = 0
!                        AND NIT = 5 OR 6.
!
!                 GRAD - INITIAL SOLUTION ESTIMATES (ZERO
!                        VECTORS ARE SUFFICIENT).  GRAD(I,J)
!                        CONTAINS COMPONENT I OF THE GRADI-
!                        ENT AT NODE J FOR I = 1,2,3 (X,Y,Z)
!                        AND J = 1,...,N.  GRAD( ,J) MUST BE
!                        ORTHOGONAL TO NODE J -- GRAD(1,J)*
!                        X(J) + GRAD(2,J)*Y(J) + GRAD(3,J)*
!                        Z(J) = 0.
!
! OUTPUT PARAMETERS - NIT - NUMBER OF GAUSS-SEIDEL ITERA-
!                           TIONS EMPLOYED.
!
!                    GRAD - ESTIMATED GRADIENTS.  SEE THE
!                           DESCRIPTION UNDER INPUT PARAME-
!                           TERS.  GRAD IS NOT CHANGED IF
!                           IER = 2.
!
!                     IER - ERROR INDICATOR
!                           IER = 0 IF THE CONVERGENCE CRI-
!                                   TERION WAS ACHIEVED.
!                           IER = 1 IF CONVERGENCE WAS NOT
!                                   ACHIEVED WITHIN NIT
!                                   ITERATIONS.
!                           IER = 2 IF N OR EPS IS OUT OF
!                                   RANGE OR NIT .LT. 0 ON
!                                   INPUT.
!
! MODULES REFERENCED BY GRADG - CONSTR, APLYRT
!
! INTRINSIC FUNCTIONS CALLED BY GRADG - ATAN, SQRT, AMAX1,
!                                       ABS
!
!***********************************************************
!
      INTEGER NN, MAXIT, ITER, K, INDF, INDL, INDX, NB
      DOUBLE PRECISION    TOL, DGMAX, XK, YK, ZK, WK, G1, G2, G3, CX, &
              SX, CY, SY, A11, A12, A22, R1, R2, XNB, YNB, &
              ZNB, ALFA, XS, YS, SINAL, D, T, DG1, DG2,&
              DGK(3)
!
! LOCAL PARAMETERS -
!
! NN =          LOCAL COPY OF N
! MAXIT =       INPUT VALUE OF NIT
! ITER =        NUMBER OF ITERATIONS USED
! K =           DO-LOOP AND NODE INDEX
! INDF,INDL =   IADJ INDICES OF THE FIRST AND LAST NEIGHBORS
!                 OF K
! INDX =        IADJ INDEX IN THE RANGE INDF,...,INDL
! NB =          NEIGHBOR OF K
! TOL =         LOCAL COPY OF EPS
! DGMAX =       MAXIMUM CHANGE IN A GRADIENT COMPONENT BE-
!                 TWEEN ITERATIONS
! XK,YK,ZK,WK = X(K), Y(K), Z(K), W(K)
! G1,G2,G3 =    COMPONENTS OF GRAD( ,K)
! CX,SX,CY,SY = COMPONENTS OF A ROTATION MAPPING K TO THE
!                 NORTH POLE (0,0,1)
! A11,A12,A22 = MATRIX COMPONENTS OF THE 2 BY 2 BLOCK A*DG
!                 = R WHERE A IS SYMMETRIC, (DG1,DG2,0) IS
!                 THE CHANGE IN THE GRADIENT AT K, AND R IS
!                 THE RESIDUAL
! R1,R2 =       COMPONENTS OF THE RESIDUAL -- DERIVATIVES OF
!                 Q WITH RESPECT TO THE COMPONENTS OF THE
!                 GRADIENT AT NODE K
! XNB,YNB,ZNB = COORDINATES OF NODE NB IN THE ROTATED COOR-
!                 DINATE SYSTEM
! ALFA =        ARC-LENGTH BETWEEN NODES K AND NB
! XS,YS =       XNB**2, YNB**2
! SINAL =       SIN(ALFA) -- MAGNITUDE OF THE VECTOR CROSS
!                 PRODUCT K X NB
! D =           ALFA*SINAL**2 -- FACTOR IN THE 2 BY 2 SYSTEM
! T =           TEMPORARY STORAGE AND FACTOR OF R1 AND R2
! DG1,DG2 =     SOLUTION OF THE 2 BY 2 SYSTEM -- FIRST 2
!                 COMPONENTS OF DGK IN THE ROTATED COORDI-
!                 NATE SYSTEM
! DGK =         CHANGE IN GRAD( ,K) FROM THE PREVIOUS ESTI-
!                 MATE
!
      NN = N
      TOL = EPS
      MAXIT = NIT
!
! ERROR CHECKS AND INITIALIZATION
!
      IF (NN .LT. 3  .OR.  TOL .LT. 0.  .OR.  MAXIT .LT. 0) &
         GO TO 5
      ITER = 0
!
! TOP OF ITERATION LOOP
!
    1 IF (ITER .EQ. MAXIT) GO TO 4
      DGMAX = 0.
      INDL = 0
      DO 3 K = 1,NN
        XK = X(K)
        YK = Y(K)
        ZK = Z(K)
        WK = W(K)
        G1 = GRAD(1,K)
        G2 = GRAD(2,K)
        G3 = GRAD(3,K)
!
!   CONSTRUCT THE ROTATION MAPPING NODE K TO THE NORTH POLE
!
        CALL CONSTR (X(K),Y(K),Z(K), CX,SX,CY,SY)
!
!   INITIALIZE COMPONENTS OF THE 2 BY 2 SYSTEM
!
        A11 = 0.
        A12 = 0.
        A22 = 0.
        R1 = 0.
        R2 = 0.
!
!   LOOP ON NEIGHBORS NB OF K
!
        INDF = INDL + 1
        INDL = IEND(K)
        DO 2 INDX = INDF,INDL
          NB = IADJ(INDX)
          IF (NB .EQ. 0) GO TO 2
!
!   COMPUTE THE COORDINATES OF NB IN THE ROTATED SYSTEM
!
          T = SX*Y(NB) + CX*Z(NB)
          YNB = CX*Y(NB) - SX*Z(NB)
          ZNB = SY*X(NB) + CY*T
          XNB = CY*X(NB) - SY*T
!
!   COMPUTE ARC-LENGTH ALFA BETWEN NB AND K, SINAL =
!     SIN(ALFA), AND D = ALFA*SIN(ALFA)**2
!
          ALFA = 2.*ATAN(SQRT((1.-ZNB)/(1.+ZNB)))
          XS = XNB*XNB
          YS = YNB*YNB
          SINAL = SQRT(XS+YS)
          D = ALFA*(XS+YS)
!
!   UPDATE THE SYSTEM COMPONENTS FOR NODE NB
!
          A11 = A11 + XS/D
          A12 = A12 + XNB*YNB/D
          A22 = A22 + YS/D
          T = 1.5*(W(NB)-WK)/(ALFA*ALFA*SINAL) + &
              ((GRAD(1,NB)*XK + GRAD(2,NB)*YK + &
                GRAD(3,NB)*ZK)/2. - (G1*X(NB) + &
                G2*Y(NB) + G3*Z(NB)))/D
          R1 = R1 + T*XNB
          R2 = R2 + T*YNB
    2     CONTINUE
!
!   SOLVE THE 2 BY 2 SYSTEM AND UPDATE DGMAX
!
        DG2 = (A11*R2 - A12*R1)/(A11*A22 - A12*A12)
        DG1 = (R1 - A12*DG2)/A11
        DGMAX = AMAX1(DGMAX,ABS(DG1),ABS(DG2))
!
!   ROTATE (DG1,DG2,0) BACK TO THE ORIGINAL COORDINATE
!     SYSTEM AND UPDATE GRAD( ,K)
!
        CALL APLYRT (DG1,DG2,CX,SX,CY,SY, DGK)
        GRAD(1,K) = G1 + DGK(1)
        GRAD(2,K) = G2 + DGK(2)
    3   GRAD(3,K) = G3 + DGK(3)
!
!   INCREMENT ITER AND TEST FOR CONVERGENCE
!
      ITER = ITER + 1
      IF (DGMAX .GT. TOL) GO TO 1
!
! METHOD CONVERGED
!
      NIT = ITER
      IER = 0
      RETURN
!
! METHOD FAILED TO CONVERGE WITHIN NIT ITERATIONS
!
    4 IER = 1
      RETURN
!
! PARAMETER OUT OF RANGE
!
    5 NIT = 0
      IER = 2
      RETURN
      END
      SUBROUTINE GRADL (N,K,X,Y,Z,W,IADJ,IEND, G,IER)
      INTEGER N, K, IADJ(1), IEND(N), IER
      DOUBLE PRECISION    X(N), Y(N), Z(N), W(N), G(3)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A TRIANGULATION OF A SET OF NODES ON THE UNIT
! SPHERE WITH THEIR ASSOCIATED DATA VALUES W, THIS ROUTINE
! ESTIMATES A GRADIENT VECTOR AT NODE K AS FOLLOWS -- THE
! COORDINATE SYSTEM IS ROTATED SO THAT K BECOMES THE NORTH
! POLE, NODE K AND A SET OF NEARBY NODES ARE PROJECTED
! ORTHOGONALLY ONTO THE X-Y PLANE (IN THE NEW COORDINATE
! SYSTEM), A QUADRATIC IS FITTED IN A WEIGHTED LEAST-SQUARES
! SENSE TO THE DATA VALUES AT THE PROJECTED NODES SUCH THAT
! THE VALUE (ASSOCIATED WITH K) AT (0,0) IS INTERPOLATED, X-
! AND Y-PARTIAL DERIVATIVE ESTIMATES DX AND DY ARE COMPUTED
! BY DIFFERENTIATING THE QUADRATIC AT (0,0), AND THE ESTI-
! MATED GRADIENT G IS OBTAINED BY ROTATING (DX,DY,0) BACK TO
! THE ORIGINAL COORDINATE SYSTEM.  NOTE THAT G LIES IN THE
! PLANE TANGENT TO THE SPHERE AT NODE K, I.E. G IS ORTHOGO-
! NAL TO THE UNIT VECTOR REPRESENTED BY NODE K.  A MARQUARDT
! STABILIZATION FACTOR IS USED IF NECESSARY TO ENSURE A
! WELL-CONDITIONED LEAST SQUARES SYSTEM, AND A UNIQUE SOLU-
! TION EXISTS UNLESS THE NODES ARE COLLINEAR.
!
! INPUT PARAMETERS -    N - NUMBER OF NODES IN THE TRIANGU-
!                           LATION.  N .GE. 7.
!
!                       K - NODE AT WHICH THE GRADIENT IS
!                           SOUGHT.  1 .LE. K .LE. N.
!
!                   X,Y,Z - CARTESIAN COORDINATES OF THE
!                           NODES.
!
!                       W - DATA VALUES AT THE NODES.  W(I)
!                           IS ASSOCIATED WITH (X(I),Y(I),
!                           Z(I)) FOR I = 1,...,N.
!
!                    IADJ - SET OF ADJACENCY LISTS OF NODES
!                           IN THE TRIANGULATION.
!
!                    IEND - POINTERS TO THE ENDS OF
!                           ADJACENCY LISTS FOR EACH NODE.
!
! IADJ AND IEND MAY BE CREATED BY SUBROUTINE TRMESH.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS -    G - X-, Y-, AND Z-COMPONENTS (IN
!                            THAT ORDER) OF THE ESTIMATED
!                            GRADIENT AT NODE K UNLESS
!                            IER .LT. 0.
!
!                      IER - ERROR INDICATOR
!                            IER .GE. 6 IF NO ERRORS WERE
!                                       ENCOUNTERED.  IER
!                                       CONTAINS THE NUMBER
!                                       OF NODES (INCLUDING
!                                       K) USED IN THE LEAST
!                                       SQUARES FIT.
!                            IER = -1 IF N OR K IS OUT OF
!                                     RANGE.
!                            IER = -2 IF THE LEAST SQUARES
!                                     SYSTEM HAS NO UNIQUE
!                                     SOLUTION DUE TO DUP-
!                                     LICATE OR COLLINEAR
!                                     NODES.
!
! MODULES REFERENCED BY GRADL - GETNP, CONSTR, APLYR,
!                               SETUP, GIVENS, ROTATE,
!                               APLYRT
!
! INTRINSIC FUNCTIONS CALLED BY GRADL - MIN0, FLOAT, SQRT,
!                                       AMIN1, ABS
!
!***********************************************************
!
      INTEGER NN, KK, LMN, LMX, LMIN, LMAX, LM1, LNP, &
              NPTS(30), IERR, NP, I, J, IM1, IP1, JP1, L
      DOUBLE PRECISION    WK, SUM, DF, RF, RTOL, AVSQ, AV, RIN, CX, SX, &
              CY, SY, XP, YP, ZP, WT, A(6,6), C, S, DMIN, &
              DTOL, SF, DX, DY
      DATA    LMN/10/
      DATA    LMX/30/, RTOL/1.E-6/, DTOL/.01/, SF/1./
!
! LOCAL PARAMETERS -
!
! NN,KK =     LOCAL COPIES OF N AND K
! LMN,LMX =   MINIMUM AND MAXIMUM VALUES OF LNP FOR N
!               SUFFICIENTLY LARGE.  IN MOST CASES LMN-1
!               NODES ARE USED IN THE FIT.  7 .LE. LMN .LE.
!               LMX.
! LMIN,LMAX = MIN(LMN,N), MIN(LMX,N)
! LM1 =       LMIN-1
! LNP =       LENGTH OF NPTS OR LMAX+1
! NPTS =      ARRAY CONTAINING THE INDICES OF A SEQUENCE OF
!               NODES ORDERED BY ANGULAR DISTANCE FROM K.
!               NPTS(1)=K AND THE FIRST LNP-1 ELEMENTS OF
!               NPTS ARE USED IN THE LEAST SQUARES FIT.
!               UNLESS LNP = LMAX+1, NPTS(LNP) DETERMINES R
!               (SEE RIN).
! IERR =      ERROR FLAG FOR CALLS TO GETNP (NOT CHECKED)
! NP =        ELEMENT OF NPTS TO BE ADDED TO THE SYSTEM
! I,J =       LOOP INDICES
! IM1,IP1 =   I-1, I+1
! JP1 =       J+1
! L =         NUMBER OF COLUMNS OF A**T TO WHICH A ROTATION
!               IS APPLIED
! WK =        W(K) -- DATA VALUE AT NODE K
! SUM =       SUM OF SQUARED EUCLIDEAN DISTANCES (IN THE
!               ROTATED COORDINATE SYSTEM) BETWEEN THE
!               ORIGIN AND THE NODES USED IN THE LEAST
!               SQUARES FIT
! DF =        NEGATIVE Z-COMPONENT (IN THE ROTATED COORDI-
!               NATE SYSTEM) OF AN ELEMENT NP OF NPTS --
!               INCREASING FUNCTION OF THE ANGULAR DISTANCE
!               BETWEEN K AND NP.  DF LIES IN THE INTERVAL
!               (-1,1).
! RF =        VALUE OF DF ASSOCIATED WITH NPTS(LNP) UNLESS
!               LNP = LMAX+1 (SEE RIN)
! RTOL =      TOLERANCE FOR DETERMINING LNP (AND HENCE R) --
!               IF THE INCREASE IN DF BETWEEN TWO SUCCESSIVE
!               ELEMENTS OF NPTS IS LESS THAN RTOL, THEY ARE
!               TREATED AS BEING THE SAME DISTANCE FROM NODE
!               K AND AN ADDITIONAL NODE IS ADDED
! AVSQ =      AV*AV -- ACCUMULATED IN SUM
! AV =        ROOT-MEAN-SQUARE DISTANCE (IN THE ROTATED
!               COORDINATE SYSTEM) BETWEEN THE ORIGIN AND
!               THE NODES (OTHER THAN K) IN THE LEAST
!               SQUARES FIT.  THE FIRST 3 COLUMNS OF A**T
!               ARE SCALED BY 1/AVSQ, THE NEXT 2 BY 1/AV.
! RIN =       INVERSE OF A RADIUS OF INFLUENCE R WHICH
!               ENTERS INTO WT -- R = 1+RF UNLESS ALL ELE-
!               MENTS OF NPTS ARE USED IN THE FIT (LNP =
!               LMAX+1), IN WHICH CASE R IS THE DISTANCE
!               FUNCTION ASSOCIATED WITH SOME POINT MORE
!               DISTANT FROM K THAN NPTS(LMAX)
! CX,SX =     COMPONENTS OF A PLANE ROTATION ABOUT THE X-
!               AXIS WHICH, TOGETHER WITH CY AND SY, DEFINE
!               A MAPPING FROM NODE K TO THE NORTH POLE
!               (0,0,1)
! CY,SY =     COMPONENTS OF A PLANE ROTATION ABOUT THE Y-
!               AXIS
! XP,YP,ZP =  COORDINATES OF NP IN THE ROTATED COORDINATE
!               SYSTEM UNLESS ZP .LT. 0, IN WHICH CASE
!               (XP,YP,0) LIES ON THE EQUATOR
! WT =        WEIGHT FOR THE EQUATION CORRESPONDING TO NP --
!               WT = (R-D)/(R*D) = 1/D - RIN WHERE D = 1-ZP
!               IS ASSOCIATED WITH NP
! A =         TRANSPOSE OF THE (UPPER TRIANGLE OF THE) AUG-
!               MENTED REGRESSION MATRIX
! C,S =       COMPONENTS OF THE PLANE ROTATION USED TO
!               TRIANGULARIZE THE REGRESSION MATRIX
! DMIN =      MINIMUM OF THE MAGNITUDES OF THE DIAGONAL
!               ELEMENTS OF THE TRIANGULARIZED REGRESSION
!               MATRIX
! DTOL =      TOLERANCE FOR DETECTING AN ILL-CONDITIONED
!               SYSTEM -- DMIN IS REQUIRED TO BE AT LEAST
!               DTOL
! SF =        MARQUARDT STABILIZATION FACTOR USED TO DAMP
!               OUT THE FIRST 3 SOLUTION COMPONENTS (SECOND
!               PARTIALS OF THE QUADRATIC) WHEN THE SYSTEM
!               IS ILL-CONDITIONED.  INCREASING SF RESULTS
!               IN MORE DAMPING (A MORE NEARLY LINEAR FIT).
! DX,DY =     X AND Y COMPONENTS OF THE ESTIMATED GRADIENT
!               IN THE ROTATED COORDINATE SYSTEM
!
      NN = N
      KK = K
      WK = W(KK)
!
! CHECK FOR ERRORS AND INITIALIZE LMIN, LMAX
!
      IF (NN .LT. 7  .OR.  KK .LT. 1  .OR.  KK .GT. NN) &
         GO TO 13
      LMIN = MIN0(LMN,NN)
      LMAX = MIN0(LMX,NN)
!
! COMPUTE NPTS, LNP, AVSQ, AV, AND R.
!   SET NPTS TO THE CLOSEST LMIN-1 NODES TO K.  DF CONTAINS
!   THE NEGATIVE Z-COMPONENT (IN THE ROTATED COORDINATE
!   SYSTEM) OF THE NEW NODE ON RETURN FROM GETNP.
!
      SUM = 0.
      NPTS(1) = KK
      LM1 = LMIN - 1
      DO 1 LNP = 2,LM1
        CALL GETNP (X,Y,Z,IADJ,IEND,LNP, NPTS, DF,IERR)
    1   SUM = SUM + 1. - DF*DF
!
!   ADD ADDITIONAL NODES TO NPTS UNTIL THE INCREASE IN
!     R = 1+RF IS AT LEAST RTOL.
!
      DO 2 LNP = LMIN,LMAX
        CALL GETNP (X,Y,Z,IADJ,IEND,LNP, NPTS, RF,IERR)
        IF (RF-DF .GE. RTOL) GO TO 3
    2   SUM = SUM + 1. - RF*RF
!
!   USE ALL LMAX NODES IN THE LEAST SQUARES FIT.  R IS
!     ARBITRARILY INCREASED BY 5 PER CENT.
!
      RF = 1.05*RF + .05
      LNP = LMAX + 1
!
!   THERE ARE LNP-2 EQUATIONS CORRESPONDING TO NODES
!     NPTS(2),...,NPTS(LNP-1).
!
    3 AVSQ = SUM/FLOAT(LNP-2)
      AV = SQRT(AVSQ)
      RIN = 1./(1.+RF)
!
! CONSTRUCT THE ROTATION
!
      CALL CONSTR (X(KK),Y(KK),Z(KK), CX,SX,CY,SY)
!
! SET UP THE FIRST 5 EQUATIONS OF THE AUGMENTED REGRESSION
!   MATRIX (TRANSPOSED) AS THE COLUMNS OF A, AND ZERO OUT
!   THE LOWER TRIANGLE (UPPER TRIANGLE OF A) WITH GIVENS
!   ROTATIONS
!
      DO 5 I = 1,5
        NP = NPTS(I+1)
        CALL APLYR (X(NP),Y(NP),Z(NP),CX,SX,CY,SY, XP,YP,ZP)
        WT = 1./(1.-ZP) - RIN
        CALL SETUP (XP,YP,W(NP),WK,AV,AVSQ,WT, A(1,I))
        IF (I .EQ. 1) GO TO 5
        IM1 = I - 1
        DO 4 J = 1,IM1
          JP1 = J + 1
          L = 6 - J
          CALL GIVENS ( A(J,J),A(J,I), C,S)
    4     CALL ROTATE (L,C,S, A(JP1,J),A(JP1,I) )
    5   CONTINUE
!
! ADD THE ADDITIONAL EQUATIONS TO THE SYSTEM USING
!   THE LAST COLUMN OF A -- I .LE. LNP.
!
      I = 7
    6   IF (I .EQ. LNP) GO TO 8
        NP = NPTS(I)
        CALL APLYR (X(NP),Y(NP),Z(NP),CX,SX,CY,SY, XP,YP,ZP)
        WT = 1./(1.-ZP) - RIN
        CALL SETUP (XP,YP,W(NP),WK,AV,AVSQ,WT, A(1,6))
        DO 7 J = 1,5
          JP1 = J + 1
          L = 6 - J
          CALL GIVENS ( A(J,J),A(J,6), C,S)
    7     CALL ROTATE (L,C,S, A(JP1,J),A(JP1,6) )
        I = I + 1
        GO TO 6
!
! TEST THE SYSTEM FOR ILL-CONDITIONING
!
    8 DMIN = AMIN1( ABS(A(1,1)),ABS(A(2,2)),ABS(A(3,3)), &
                    ABS(A(4,4)),ABS(A(5,5)) )
      IF (DMIN .GE. DTOL) GO TO 12
      IF (LNP .GT. LMAX) GO TO 9
!
! ADD ANOTHER NODE TO THE SYSTEM AND INCREASE R --
!   I .EQ. LNP
!
      LNP = LNP + 1
      IF (LNP .LE. LMAX) CALL GETNP (X,Y,Z,IADJ,IEND,LNP, &
                                     NPTS, RF,IERR)
      RIN = 1./(1.05*(1.+RF))
      GO TO 6
!
! STABILIZE THE SYSTEM BY DAMPING SECOND PARTIALS --ADD
!   MULTIPLES OF THE FIRST THREE UNIT VECTORS TO THE FIRST
!   THREE EQUATIONS.
!
    9 DO 11 I = 1,3
        A(I,6) = SF
        IP1 = I + 1
        DO 10 J = IP1,6
   10     A(J,6) = 0.
        DO 11 J = I,5
          JP1 = J + 1
          L = 6 - J
          CALL GIVENS ( A(J,J),A(J,6), C,S)
   11     CALL ROTATE (L,C,S, A(JP1,J),A(JP1,6) )
!
! TEST THE LINEAR PORTION OF THE STABILIZED SYSTEM FOR
!   ILL-CONDITIONING
!
      DMIN = AMIN1( ABS(A(4,4)),ABS(A(5,5)) )
      IF (DMIN .LT. DTOL) GO TO 14
!
! SOLVE THE 2 BY 2 TRIANGULAR SYSTEM FOR THE ESTIMATED
!   PARTIAL DERIVATIVES
!
   12 DY = A(6,5)/A(5,5)
      DX = (A(6,4) - A(5,4)*DY)/A(4,4)/AV
      DY = DY/AV
!
! ROTATE THE GRADIENT (DX,DY,0) BACK INTO THE ORIGINAL
!   COORDINATE SYSTEM
!
      CALL APLYRT (DX,DY,CX,SX,CY,SY, G)
      IER = LNP - 1
      RETURN
!
! N OR K IS OUT OF RANGE
!
   13 IER = -1
      RETURN
!
! NO UNIQUE SOLUTION DUE TO COLLINEAR NODES
!
   14 IER = -2
      RETURN
      END
      INTEGER FUNCTION INDX (NVERTX,NABOR,IADJ,IEND)
      INTEGER NVERTX, NABOR, IADJ(1), IEND(1)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS FUNCTION RETURNS THE INDEX OF NABOR IN THE
! ADJACENCY LIST FOR NVERTX.
!
! INPUT PARAMETERS - NVERTX - NODE WHOSE ADJACENCY LIST IS
!                             TO BE SEARCHED.
!
!                     NABOR - NODE WHOSE INDEX IS TO BE
!                             RETURNED.  NABOR MUST BE
!                             CONNECTED TO NVERTX.
!
!                      IADJ - SET OF ADJACENCY LISTS.
!
!                      IEND - POINTERS TO THE ENDS OF
!                             ADJACENCY LISTS IN IADJ.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS FUNCTION.
!
! OUTPUT PARAMETER -  INDEX - IADJ(INDEX) = NABOR.
!
! MODULES REFERENCED BY INDEX - NONE
!
!***********************************************************
!
      INTEGER NB, INDXX
!
! LOCAL PARAMETERS -
!
! NB =   LOCAL COPY OF NABOR
! INDXX = INDEX FOR IADJ
!
      NB = NABOR
!
! INITIALIZATION
!
      INDXX = IEND(NVERTX) + 1
!
! SEARCH THE LIST OF NVERTX NEIGHBORS FOR NB
!
    1 INDXX = INDXX - 1
      IF (IADJ(INDXX) .NE. NB) GO TO 1
!
      INDX = INDXX
      RETURN
      END
      SUBROUTINE WVAL (B1,B2,B3,V1,V2,V3,W1,W2,W3,G1,G2,G3, &
                       IFLAG, PW,PG)
      INTEGER IFLAG
      DOUBLE PRECISION    B1, B2, B3, V1(3), V2(3), V3(3), W1, W2, W3, &
              G1(3), G2(3), G3(3), PW, PG(3)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN DATA VALUES AND GRADIENTS AT THE THREE VERTICES OF
! A SPHERICAL TRIANGLE CONTAINING A POINT P, THIS ROUTINE
! COMPUTES THE VALUE AND, OPTIONALLY, THE GRADIENT OF F AT P
! WHERE F INTERPOLATES THE VERTEX DATA.  ALONG THE TRIANGLE
! EDGES, THE INTERPOLATORY FUNCTION F IS THE HERMITE CUBIC
! (WITH RESPECT TO ARC-LENGTH) INTERPOLANT OF THE VALUES AND
! TANGENTIAL GRADIENT COMPONENTS AT THE ENDPOINTS, AND THE
! DERIVATIVE NORMAL TO THE ARC VARIES LINEARLY WITH RESPECT
! TO ARC-LENGTH BETWEEN THE NORMAL GRADIENT COMPONENTS AT
! THE ENDPOINTS.  THUS THE METHOD YIELDS C-1 CONTINUITY WHEN
! USED TO INTERPOLATE OVER A TRIANGULATION.  THE INTERPOLANT
! USES A FIRST-ORDER BLENDING METHOD ON THE UNDERLYING
! PLANAR TRIANGLE.
!
! INPUT PARAMETERS - B1,B2,B3 - BARYCENTRIC COORDINATES OF
!                               PP WITH RESPECT TO THE
!                               (PLANAR) UNDERLYING TRIANGLE
!                               (V1,V2,V3) WHERE PP IS THE
!                               CENTRAL PROJECTION OF P ONTO
!                               THIS TRIANGLE.
!
!                    V1,V2,V3 - CARTESIAN COORDINATES OF THE
!                               VERTICES OF A SPHERICAL TRI-
!                               ANGLE CONTAINING P.  V3 LEFT
!                               V1->V2.
!
!                    W1,W2,W3 - DATA VALUES ASSOCIATED WITH
!                               THE VERTICES.
!
!                    G1,G2,G3 - GRADIENTS ASSOCIATED WITH
!                               THE VERTICES.  GI IS ORTHOG-
!                               ONAL TO VI FOR I = 1,2,3.
!
!                       IFLAG - OPTION INDICATOR
!                               IFLAG = 0 IF ONLY PW IS TO
!                                         BE COMPUTED.
!                               IFLAG = 1 IF BOTH PW AND PG
!                                         ARE DESIRED.
!
!                          PG - VECTOR OF LENGTH 3 IF IFLAG
!                               = 1, NOT USED OTHERWISE.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS -      PW - INTERPOLATED VALUE AT P.
!
!                          PG - INTERPOLATED GRADIENT AT P
!                               (ORTHOGONAL TO P) IF IFLAG
!                               = 1.
!
! EACH VECTOR V ABOVE CONTAINS X-, Y-, AND Z-COMPONENTS IN
!   V(1), V(2), AND V(3), RESPECTIVELY.
!
! MODULES REFERENCED BY WVAL - ARCINT, ARCLEN
!
! INTRINSIC FUNCTION CALLED BY WVAL - SQRT
!
!***********************************************************
!
      INTEGER I
      DOUBLE PRECISION    C1,C2,C3, SUMM, U1(3), U2(3), U3(3), U1N, &
              U2N, U3N, Q1(3), Q2(3), Q3(3), VAL, W, G(3),&
              DUM
!
! LOCAL PARAMETERS -
!
! I =           DO-LOOP INDEX
! C1,C2,C3 =    COEFFICIENTS (WEIGHT FUNCTIONS) OF PARTIAL
!                 INTERPOLANTS.  C1 = 1 ON THE EDGE OPPOSITE
!                 V1 AND C1 = 0 ON THE OTHER EDGES.  SIMI-
!                 LARLY FOR C2 AND C3.  C1+C2+C3 = 1.
! SUMM =         QUANTITY USED TO NORMALIZE C1, C2, AND C3
! U1,U2,U3 =    POINTS ON THE BOUNDARY OF THE PLANAR TRIAN-
!                 GLE AND LYING ON THE LINES CONTAINING PP
!                 AND THE VERTICES.  U1 IS OPPOSITE V1, ETC.
! U1N,U2N,U3N = QUANTITIES USED TO COMPUTE Q1, Q2, AND Q3
!                 (MAGNITUDES OF U1, U2, AND U3)
! Q1,Q2,Q3 =    CENTRAL PROJECTIONS OF U1, U2, AND U3 ONTO
!                 THE SPHERE AND THUS LYING ON AN ARC OF THE
!                 SPHERICAL TRIANGLE
! VAL =         LOCAL VARIABLE USED TO ACCUMULATE THE CON-
!                 TRIBUTIONS TO PW
! W,G =         VALUE AND GRADIENT AT Q1, Q2, OR Q3 OBTAINED
!                 BY INTERPOLATION ALONG ONE OF THE ARCS OF
!                 THE SPHERICAL TRIANGLE
! DUM =         DUMMY VARIABLE FOR THE CALL TO ARCINT
!
!
! COMPUTE WEIGHT FUNCTIONS C1, C2, AND C3
!
      C1 = B2*B3
      C2 = B3*B1
      C3 = B1*B2
      SUMM = C1 + C2 + C3
      IF (SUMM .GT. 0.) GO TO 2
!
! P COINCIDES WITH A VERTEX
!
      PW = B1*W1 + B2*W2 + B3*W3
      IF (IFLAG .NE. 1) RETURN
      DO 1 I = 1,3
    1   PG(I) = B1*G1(I) + B2*G2(I) + B3*G3(I)
      RETURN
!
! NORMALIZE C1, C2, AND C3
!
    2 C1 = C1/SUMM
      C2 = C2/SUMM
      C3 = C3/SUMM
!
! COMPUTE (U1,U2,U3) AND (U1N,U2N,U3N)
!
      U1N = 0.
      U2N = 0.
      U3N = 0.
      DO 3 I = 1,3
        U1(I) = (B2*V2(I) + B3*V3(I))/(B2+B3)
        U2(I) = (B3*V3(I) + B1*V1(I))/(B3+B1)
        U3(I) = (B1*V1(I) + B2*V2(I))/(B1+B2)
        U1N = U1N + U1(I)*U1(I)
        U2N = U2N + U2(I)*U2(I)
    3   U3N = U3N + U3(I)*U3(I)
!
! COMPUTE Q1, Q2, AND Q3
!
      U1N = SQRT(U1N)
      U2N = SQRT(U2N)
      U3N = SQRT(U3N)
      DO 4 I = 1,3
        Q1(I) = U1(I)/U1N
        Q2(I) = U2(I)/U2N
    4   Q3(I) = U3(I)/U3N
!
! COMPUTE INTERPOLATED VALUE (VAL) AT P BY LOOPING ON
!   TRIANGLE SIDES
!
      VAL = 0.
!
! CONTRIBUTION FROM SIDE OPPOSITE V1 --
!
!   COMPUTE VALUE AND GRADIENT AT Q1 BY INTERPOLATING
!     BETWEEN V2 AND V3
!
      CALL ARCINT (Q1,V2,V3,W2,W3,G2,G3, W,G,DUM)
!
!   ADD IN THE CONTRIBUTION
!
      VAL = VAL + C1*( W + B1*B1*(3.-2.*B1)*(W1-W) + &
                  B1*(1.-B1)* &
                  (B1*(G1(1)*U1(1)+G1(2)*U1(2)+G1(3)*U1(3)) &
                   + (1.-B1)* &
                   (G(1)*V1(1)+G(2)*V1(2)+G(3)*V1(3))/U1N) )
!
! CONTRIBUTION FROM SIDE OPPOSITE V2 --
!
!   COMPUTE VALUE AND GRADIENT AT Q2 BY INTERPOLATING
!     BETWEEN V3 AND V1
!
      CALL ARCINT (Q2,V3,V1,W3,W1,G3,G1, W,G,DUM)
!
!   ADD IN THE CONTRIBUTION
!
      VAL = VAL + C2*( W + B2*B2*(3.-2.*B2)*(W2-W) + &
                  B2*(1.-B2)* &
                  (B2*(G2(1)*U2(1)+G2(2)*U2(2)+G2(3)*U2(3)) &
                   + (1.-B2)* &
                   (G(1)*V2(1)+G(2)*V2(2)+G(3)*V2(3))/U2N) )
!
! CONTRIBUTION FROM SIDE OPPOSITE V3 --
!
!   COMPUTE INTERPOLATED VALUE AND GRADIENT AT Q3
!     BY INTERPOLATING BETWEEN V1 AND V2
!
      CALL ARCINT (Q3,V1,V2,W1,W2,G1,G2, W,G,DUM)
!
!   ADD IN THE FINAL CONTRIBUTION
!
      VAL = VAL + C3*( W + B3*B3*(3.-2.*B3)*(W3-W) + &
                  B3*(1.-B3)* &
                  (B3*(G3(1)*U3(1)+G3(2)*U3(2)+G3(3)*U3(3)) &
                   + (1.-B3)* &
                   (G(1)*V3(1)+G(2)*V3(2)+G(3)*V3(3))/U3N) )
      PW = VAL
      IF (IFLAG .NE. 1) RETURN
      RETURN
      END
      subroutine intrpc1_n(npts,nptso,olats,olons,x,y,z,datain,iadj,&
                           iend,odata,ierr)
      integer, intent(in) :: npts, nptso
      integer, intent(out) :: ierr
      double precision, intent(in), dimension(nptso) :: olats,olons
      double precision, intent(in), dimension(npts) :: datain,x,y,z
      double precision, intent(out), dimension(nptso) :: odata
      double precision, dimension(3,nptso) :: grad
      integer, intent(in), dimension(npts) :: iend
      integer, intent(in), dimension(6*(npts-1)) :: iadj
      integer n,ierr1,ist
      ist = 1
      ierr = 0
      do n=1,nptso
         call intrc1(npts,olats(n),olons(n),x,y,z,datain,iadj,iend,&
                     0,grad,ist,odata(n),ierr1)
         if (ierr1 .ne. 0) then
           !print *,n,'warning: ierr = ',ierr1,' in intrc1_n'
           !print *,olats(n), olons(n), npts
           !stop
           ierr = ierr + ierr1
         endif
      enddo
      end subroutine intrpc1_n
      subroutine intrpc0_n(npts,nptso,olats,olons,x,y,z,datain,iadj,&
                           iend,odata,ierr)
      integer, intent(in) :: npts, nptso
      integer, intent(out) :: ierr
      double precision, intent(in), dimension(nptso) :: olats,olons
      double precision, intent(in), dimension(npts) :: datain,x,y,z
      double precision, intent(out), dimension(nptso) :: odata
      integer, intent(in), dimension(npts) :: iend
      integer, intent(in), dimension(6*(npts-1)) :: iadj
      integer n,ierr1,ist
      ist = 1
      ierr = 0
      do n=1,nptso
         call intrc0(npts,olats(n),olons(n),x,y,z,datain,iadj,iend,&
                     ist,odata(n),ierr1)
         if (ierr1 .ne. 0) then
           !print *,n,'warning: ierr = ',ierr1,' in intrc0_n'
           !print *,olats(n), olons(n), npts
           !stop
           ierr = ierr + ierr1
         endif
      enddo
      end subroutine intrpc0_n
      subroutine intrpnn_n(npts,nptso,olats,olons,x,y,z,datain,iadj,&
                           iend,odata,ierr)
      integer, intent(in) :: npts, nptso
      integer, intent(out) :: ierr
      double precision, intent(in), dimension(nptso) :: olats,olons
      double precision, intent(in), dimension(npts) :: datain,x,y,z
      double precision, intent(out), dimension(nptso) :: odata
      integer, intent(in), dimension(npts) :: iend
      integer, intent(in), dimension(6*(npts-1)) :: iadj
      integer n,ierr1,ist
      ist = 1
      ierr = 0
      do n=1,nptso
         call intrnn(npts,olats(n),olons(n),x,y,z,datain,iadj,iend,&
                     ist,odata(n),ierr1)
         if (ierr1 .ne. 0) then
           !print *,n,'warning: ierr = ',ierr1,' in intrc0_n'
           !print *,olats(n), olons(n), npts
           !stop
           ierr = ierr + ierr1
         endif
      enddo
      end subroutine intrpnn_n
