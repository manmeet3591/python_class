      SUBROUTINE ADDNOD (NST,K,X,Y,Z, LIST,LPTR,LEND,
     .                   LNEW, IER)
      INTEGER NST, K, LIST(*), LPTR(*), LEND(K), LNEW, IER
      DOUBLE PRECISION    X(K), Y(K), Z(K)
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/08/99
C
C   This subroutine adds node K to a triangulation of the
C convex hull of nodes 1,...,K-1, producing a triangulation
C of the convex hull of nodes 1,...,K.
C
C   The algorithm consists of the following steps:  node K
C is located relative to the triangulation (TRFIND), its
C index is added to the data structure (INTADD or BDYADD),
C and a sequence of swaps (SWPTST and SWAP) are applied to
C the arcs opposite K so that all arcs incident on node K
C and opposite node K are locally optimal (satisfy the cir-
C cumcircle test).  Thus, if a Delaunay triangulation is
C input, a Delaunay triangulation will result.
C
C
C On input:
C
C       NST = Index of a node at which TRFIND begins its
C             search.  Search time depends on the proximity
C             of this node to K.  If NST < 1, the search is
C             begun at node K-1.
C
C       K = Nodal index (index for X, Y, Z, and LEND) of the
C           new node to be added.  K .GE. 4.
C
C       X,Y,Z = Arrays of length .GE. K containing Car-
C               tesian coordinates of the nodes.
C               (X(I),Y(I),Z(I)) defines node I for
C               I = 1,...,K.
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR,LEND,LNEW = Data structure associated with
C                             the triangulation of nodes 1
C                             to K-1.  The array lengths are
C                             assumed to be large enough to
C                             add node K.  Refer to Subrou-
C                             tine TRMESH.
C
C On output:
C
C       LIST,LPTR,LEND,LNEW = Data structure updated with
C                             the addition of node K as the
C                             last entry unless IER .NE. 0
C                             and IER .NE. -3, in which case
C                             the arrays are not altered.
C
C       IER = Error indicator:
C             IER =  0 if no errors were encountered.
C             IER = -1 if K is outside its valid range
C                      on input.
C             IER = -2 if all nodes (including K) are col-
C                      linear (lie on a common geodesic).
C             IER =  L if nodes L and K coincide for some
C                      L < K.
C
C Modules required by ADDNOD:  BDYADD, COVSPH, INSERT,
C                                INTADD, JRAND, LSTPTR,
C                                STORE, SWAP, SWPTST,
C                                TRFIND
C
C Intrinsic function called by ADDNOD:  ABS
C
C***********************************************************
C
      INTEGER LSTPTR
      INTEGER I1, I2, I3, IO1, IO2, IN1, IST, KK, KM1, L,
     .        LP, LPF, LPO1, LPO1S
      LOGICAL SWPTST
      DOUBLE PRECISION    B1, B2, B3, P(3)
C
C Local parameters:
C
C B1,B2,B3 = Unnormalized barycentric coordinates returned
C              by TRFIND.
C I1,I2,I3 = Vertex indexes of a triangle containing K
C IN1 =      Vertex opposite K:  first neighbor of IO2
C              that precedes IO1.  IN1,IO1,IO2 are in
C              counterclockwise order.
C IO1,IO2 =  Adjacent neighbors of K defining an arc to
C              be tested for a swap
C IST =      Index of node at which TRFIND begins its search
C KK =       Local copy of K
C KM1 =      K-1
C L =        Vertex index (I1, I2, or I3) returned in IER
C              if node K coincides with a vertex
C LP =       LIST pointer
C LPF =      LIST pointer to the first neighbor of K
C LPO1 =     LIST pointer to IO1
C LPO1S =    Saved value of LPO1
C P =        Cartesian coordinates of node K
C
      KK = K
      IF (KK .LT. 4) GO TO 3
C
C Initialization:
C
      KM1 = KK - 1
      IST = NST
      IF (IST .LT. 1) IST = KM1
      P(1) = X(KK)
      P(2) = Y(KK)
      P(3) = Z(KK)
C
C Find a triangle (I1,I2,I3) containing K or the rightmost
C   (I1) and leftmost (I2) visible boundary nodes as viewed
C   from node K.
C
      CALL TRFIND (IST,P,KM1,X,Y,Z,LIST,LPTR,LEND, B1,B2,B3,
     .             I1,I2,I3)
C
C   Test for collinear or duplicate nodes.
C
      IF (I1 .EQ. 0) GO TO 4
      IF (I3 .NE. 0) THEN
        L = I1
        IF (P(1) .EQ. X(L)  .AND.  P(2) .EQ. Y(L)  .AND.
     .      P(3) .EQ. Z(L)) GO TO 5
        L = I2
        IF (P(1) .EQ. X(L)  .AND.  P(2) .EQ. Y(L)  .AND.
     .      P(3) .EQ. Z(L)) GO TO 5
        L = I3
        IF (P(1) .EQ. X(L)  .AND.  P(2) .EQ. Y(L)  .AND.
     .      P(3) .EQ. Z(L)) GO TO 5
        CALL INTADD (KK,I1,I2,I3, LIST,LPTR,LEND,LNEW )
      ELSE
        IF (I1 .NE. I2) THEN
          CALL BDYADD (KK,I1,I2, LIST,LPTR,LEND,LNEW )
        ELSE
          CALL COVSPH (KK,I1, LIST,LPTR,LEND,LNEW )
        ENDIF
      ENDIF
      IER = 0
C
C Initialize variables for optimization of the
C   triangulation.
C
      LP = LEND(KK)
      LPF = LPTR(LP)
      IO2 = LIST(LPF)
      LPO1 = LPTR(LPF)
      IO1 = ABS(LIST(LPO1))
C
C Begin loop:  find the node opposite K.
C
    1 LP = LSTPTR(LEND(IO1),IO2,LIST,LPTR)
        IF (LIST(LP) .LT. 0) GO TO 2
        LP = LPTR(LP)
        IN1 = ABS(LIST(LP))
C
C Swap test:  if a swap occurs, two new arcs are
C             opposite K and must be tested.
C
        LPO1S = LPO1
        IF ( .NOT. SWPTST(IN1,KK,IO1,IO2,X,Y,Z) ) GO TO 2
        CALL SWAP (IN1,KK,IO1,IO2, LIST,LPTR,LEND, LPO1)
        IF (LPO1 .EQ. 0) THEN
C
C   A swap is not possible because KK and IN1 are already
C     adjacent.  This error in SWPTST only occurs in the
C     neutral case and when there are nearly duplicate
C     nodes.
C
          LPO1 = LPO1S
          GO TO 2
        ENDIF
        IO1 = IN1
        GO TO 1
C
C No swap occurred.  Test for termination and reset
C   IO2 and IO1.
C
    2   IF (LPO1 .EQ. LPF  .OR.  LIST(LPO1) .LT. 0) RETURN
        IO2 = IO1
        LPO1 = LPTR(LPO1)
        IO1 = ABS(LIST(LPO1))
        GO TO 1
C
C KK < 4.
C
    3 IER = -1
      RETURN
C
C All nodes are collinear.
C
    4 IER = -2
      RETURN
C
C Nodes L and K coincide.
C
    5 IER = L
      RETURN
      END
      DOUBLE PRECISION FUNCTION AREAS (V1,V2,V3)
      DOUBLE PRECISION V1(3), V2(3), V3(3)
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/18/90
C
C   This function returns the area of a spherical triangle
C on the unit sphere.
C
C
C On input:
C
C       V1,V2,V3 = Arrays of length 3 containing the Carte-
C                  sian coordinates of unit vectors (the
C                  three triangle vertices in any order).
C                  These vectors, if nonzero, are implicitly
C                  scaled to have length 1.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       AREAS = Area of the spherical triangle defined by
C               V1, V2, and V3 in the range 0 to 2*PI (the
C               area of a hemisphere).  AREAS = 0 (or 2*PI)
C               if and only if V1, V2, and V3 lie in (or
C               close to) a plane containing the origin.
C
C Modules required by AREAS:  None
C
C Intrinsic functions called by AREAS:  ACOS, DBLE, REAL,
C                                         SQRT
C
C***********************************************************
C
      DOUBLE PRECISION A1, A2, A3, CA1, CA2, CA3, DV1(3),
     .                 DV2(3), DV3(3), S12, S23, S31,
     .                 U12(3), U23(3), U31(3)
      INTEGER          I
C
C Local parameters:
C
C A1,A2,A3 =    Interior angles of the spherical triangle
C CA1,CA2,CA3 = cos(A1), cos(A2), and cos(A3), respectively
C DV1,DV2,DV3 = Double Precision copies of V1, V2, and V3
C I =           DO-loop index and index for Uij
C S12,S23,S31 = Sum of squared components of U12, U23, U31
C U12,U23,U31 = Unit normal vectors to the planes defined by
C                 pairs of triangle vertices
C
      DO 1 I = 1,3
        DV1(I) = DBLE(V1(I))
        DV2(I) = DBLE(V2(I))
        DV3(I) = DBLE(V3(I))
    1   CONTINUE
C
C Compute cross products Uij = Vi X Vj.
C
      U12(1) = DV1(2)*DV2(3) - DV1(3)*DV2(2)
      U12(2) = DV1(3)*DV2(1) - DV1(1)*DV2(3)
      U12(3) = DV1(1)*DV2(2) - DV1(2)*DV2(1)
C
      U23(1) = DV2(2)*DV3(3) - DV2(3)*DV3(2)
      U23(2) = DV2(3)*DV3(1) - DV2(1)*DV3(3)
      U23(3) = DV2(1)*DV3(2) - DV2(2)*DV3(1)
C
      U31(1) = DV3(2)*DV1(3) - DV3(3)*DV1(2)
      U31(2) = DV3(3)*DV1(1) - DV3(1)*DV1(3)
      U31(3) = DV3(1)*DV1(2) - DV3(2)*DV1(1)
C
C Normalize Uij to unit vectors.
C
      S12 = 0.D0
      S23 = 0.D0
      S31 = 0.D0
      DO 2 I = 1,3
        S12 = S12 + U12(I)*U12(I)
        S23 = S23 + U23(I)*U23(I)
        S31 = S31 + U31(I)*U31(I)
    2   CONTINUE
C
C Test for a degenerate triangle associated with collinear
C   vertices.
C
      IF (S12 .EQ. 0.D0  .OR.  S23 .EQ. 0.D0  .OR.
     .    S31 .EQ. 0.D0) THEN
        AREAS = 0.
        RETURN
      ENDIF
      S12 = SQRT(S12)
      S23 = SQRT(S23)
      S31 = SQRT(S31)
      DO 3 I = 1,3
        U12(I) = U12(I)/S12
        U23(I) = U23(I)/S23
        U31(I) = U31(I)/S31
    3   CONTINUE
C
C Compute interior angles Ai as the dihedral angles between
C   planes:
C           CA1 = cos(A1) = -<U12,U31>
C           CA2 = cos(A2) = -<U23,U12>
C           CA3 = cos(A3) = -<U31,U23>
C
      CA1 = -U12(1)*U31(1)-U12(2)*U31(2)-U12(3)*U31(3)
      CA2 = -U23(1)*U12(1)-U23(2)*U12(2)-U23(3)*U12(3)
      CA3 = -U31(1)*U23(1)-U31(2)*U23(2)-U31(3)*U23(3)
      IF (CA1 .LT. -1.D0) CA1 = -1.D0
      IF (CA1 .GT. 1.D0) CA1 = 1.D0
      IF (CA2 .LT. -1.D0) CA2 = -1.D0
      IF (CA2 .GT. 1.D0) CA2 = 1.D0
      IF (CA3 .LT. -1.D0) CA3 = -1.D0
      IF (CA3 .GT. 1.D0) CA3 = 1.D0
      A1 = ACOS(CA1)
      A2 = ACOS(CA2)
      A3 = ACOS(CA3)
C
C Compute AREAS = A1 + A2 + A3 - PI.
C
      AREAS = REAL(A1 + A2 + A3 - ACOS(-1.D0))
      IF (AREAS .LT. 0.) AREAS = 0.
      RETURN
      END
      SUBROUTINE BDYADD (KK,I1,I2, LIST,LPTR,LEND,LNEW )
      INTEGER KK, I1, I2, LIST(*), LPTR(*), LEND(*), LNEW
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/11/96
C
C   This subroutine adds a boundary node to a triangulation
C of a set of KK-1 points on the unit sphere.  The data
C structure is updated with the insertion of node KK, but no
C optimization is performed.
C
C   This routine is identical to the similarly named routine
C in TRIPACK.
C
C
C On input:
C
C       KK = Index of a node to be connected to the sequence
C            of all visible boundary nodes.  KK .GE. 1 and
C            KK must not be equal to I1 or I2.
C
C       I1 = First (rightmost as viewed from KK) boundary
C            node in the triangulation that is visible from
C            node KK (the line segment KK-I1 intersects no
C            arcs.
C
C       I2 = Last (leftmost) boundary node that is visible
C            from node KK.  I1 and I2 may be determined by
C            Subroutine TRFIND.
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR,LEND,LNEW = Triangulation data structure
C                             created by Subroutine TRMESH.
C                             Nodes I1 and I2 must be in-
C                             cluded in the triangulation.
C
C On output:
C
C       LIST,LPTR,LEND,LNEW = Data structure updated with
C                             the addition of node KK.  Node
C                             KK is connected to I1, I2, and
C                             all boundary nodes in between.
C
C Module required by BDYADD:  INSERT
C
C***********************************************************
C
      INTEGER K, LP, LSAV, N1, N2, NEXT, NSAV
C
C Local parameters:
C
C K =     Local copy of KK
C LP =    LIST pointer
C LSAV =  LIST pointer
C N1,N2 = Local copies of I1 and I2, respectively
C NEXT =  Boundary node visible from K
C NSAV =  Boundary node visible from K
C
      K = KK
      N1 = I1
      N2 = I2
C
C Add K as the last neighbor of N1.
C
      LP = LEND(N1)
      LSAV = LPTR(LP)
      LPTR(LP) = LNEW
      LIST(LNEW) = -K
      LPTR(LNEW) = LSAV
      LEND(N1) = LNEW
      LNEW = LNEW + 1
      NEXT = -LIST(LP)
      LIST(LP) = NEXT
      NSAV = NEXT
C
C Loop on the remaining boundary nodes between N1 and N2,
C   adding K as the first neighbor.
C
    1 LP = LEND(NEXT)
        CALL INSERT (K,LP, LIST,LPTR,LNEW )
        IF (NEXT .EQ. N2) GO TO 2
        NEXT = -LIST(LP)
        LIST(LP) = NEXT
        GO TO 1
C
C Add the boundary nodes between N1 and N2 as neighbors
C   of node K.
C
    2 LSAV = LNEW
      LIST(LNEW) = N1
      LPTR(LNEW) = LNEW + 1
      LNEW = LNEW + 1
      NEXT = NSAV
C
    3 IF (NEXT .EQ. N2) GO TO 4
        LIST(LNEW) = NEXT
        LPTR(LNEW) = LNEW + 1
        LNEW = LNEW + 1
        LP = LEND(NEXT)
        NEXT = LIST(LP)
        GO TO 3
C
    4 LIST(LNEW) = -N2
      LPTR(LNEW) = LSAV
      LEND(K) = LNEW
      LNEW = LNEW + 1
      RETURN
      END
      SUBROUTINE BNODES (N,LIST,LPTR,LEND, NODES,NB,NA,NT)
      INTEGER N, LIST(*), LPTR(*), LEND(N), NODES(*), NB,
     .        NA, NT
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/26/96
C
C   Given a triangulation of N nodes on the unit sphere
C created by Subroutine TRMESH, this subroutine returns an
C array containing the indexes (if any) of the counterclock-
C wise-ordered sequence of boundary nodes -- the nodes on
C the boundary of the convex hull of the set of nodes.  (The
C boundary is empty if the nodes do not lie in a single
C hemisphere.)  The numbers of boundary nodes, arcs, and
C triangles are also returned.
C
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to Subroutine
C                        TRMESH.
C
C The above parameters are not altered by this routine.
C
C       NODES = Integer array of length at least NB
C               (NB .LE. N).
C
C On output:
C
C       NODES = Ordered sequence of boundary node indexes
C               in the range 1 to N (in the first NB loca-
C               tions).
C
C       NB = Number of boundary nodes.
C
C       NA,NT = Number of arcs and triangles, respectively,
C               in the triangulation.
C
C Modules required by BNODES:  None
C
C***********************************************************
C
      INTEGER K, LP, N0, NN, NST
C
C Local parameters:
C
C K =   NODES index
C LP =  LIST pointer
C N0 =  Boundary node to be added to NODES
C NN =  Local copy of N
C NST = First element of nodes (arbitrarily chosen to be
C         the one with smallest index)
C
      NN = N
C
C Search for a boundary node.
C
      DO 1 NST = 1,NN
        LP = LEND(NST)
        IF (LIST(LP) .LT. 0) GO TO 2
    1   CONTINUE
C
C The triangulation contains no boundary nodes.
C
      NB = 0
      NA = 3*(NN-2)
      NT = 2*(NN-2)
      RETURN
C
C NST is the first boundary node encountered.  Initialize
C   for traversal of the boundary.
C
    2 NODES(1) = NST
      K = 1
      N0 = NST
C
C Traverse the boundary in counterclockwise order.
C
    3 LP = LEND(N0)
        LP = LPTR(LP)
        N0 = LIST(LP)
        IF (N0 .EQ. NST) GO TO 4
        K = K + 1
        NODES(K) = N0
        GO TO 3
C
C Store the counts.
C
    4 NB = K
      NT = 2*N - NB - 2
      NA = NT + N - 1
      RETURN
      END
      SUBROUTINE CIRCUM (V1,V2,V3, C,IER)
      INTEGER IER
      DOUBLE PRECISION    V1(3), V2(3), V3(3), C(3)
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/29/95
C
C   This subroutine returns the circumcenter of a spherical
C triangle on the unit sphere:  the point on the sphere sur-
C face that is equally distant from the three triangle
C vertices and lies in the same hemisphere, where distance
C is taken to be arc-length on the sphere surface.
C
C
C On input:
C
C       V1,V2,V3 = Arrays of length 3 containing the Carte-
C                  sian coordinates of the three triangle
C                  vertices (unit vectors) in CCW order.
C
C The above parameters are not altered by this routine.
C
C       C = Array of length 3.
C
C On output:
C
C       C = Cartesian coordinates of the circumcenter unless
C           IER > 0, in which case C is not defined.  C =
C           (V2-V1) X (V3-V1) normalized to a unit vector.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if V1, V2, and V3 lie on a common
C                     line:  (V2-V1) X (V3-V1) = 0.
C             (The vertices are not tested for validity.)
C
C Modules required by CIRCUM:  None
C
C Intrinsic function called by CIRCUM:  SQRT
C
C***********************************************************
C
      INTEGER I
      DOUBLE PRECISION    CNORM, CU(3), E1(3), E2(3)
C
C Local parameters:
C
C CNORM = Norm of CU:  used to compute C
C CU =    Scalar multiple of C:  E1 X E2
C E1,E2 = Edges of the underlying planar triangle:
C           V2-V1 and V3-V1, respectively
C I =     DO-loop index
C
      DO 1 I = 1,3
        E1(I) = V2(I) - V1(I)
        E2(I) = V3(I) - V1(I)
    1   CONTINUE
C
C Compute CU = E1 X E2 and CNORM**2.
C
      CU(1) = E1(2)*E2(3) - E1(3)*E2(2)
      CU(2) = E1(3)*E2(1) - E1(1)*E2(3)
      CU(3) = E1(1)*E2(2) - E1(2)*E2(1)
      CNORM = CU(1)*CU(1) + CU(2)*CU(2) + CU(3)*CU(3)
C
C The vertices lie on a common line if and only if CU is
C   the zero vector.
C
      IF (CNORM .NE. 0.) THEN
C
C   No error:  compute C.
C
        CNORM = SQRT(CNORM)
        DO 2 I = 1,3
          C(I) = CU(I)/CNORM
    2     CONTINUE
        IER = 0
      ELSE
C
C   CU = 0.
C
        IER = 1
      ENDIF
      RETURN
      END
      SUBROUTINE COVSPH (KK,N0, LIST,LPTR,LEND,LNEW )
      INTEGER KK, N0, LIST(*), LPTR(*), LEND(*), LNEW
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/17/96
C
C   This subroutine connects an exterior node KK to all
C boundary nodes of a triangulation of KK-1 points on the
C unit sphere, producing a triangulation that covers the
C sphere.  The data structure is updated with the addition
C of node KK, but no optimization is performed.  All boun-
C dary nodes must be visible from node KK.
C
C
C On input:
C
C       KK = Index of the node to be connected to the set of
C            all boundary nodes.  KK .GE. 4.
C
C       N0 = Index of a boundary node (in the range 1 to
C            KK-1).  N0 may be determined by Subroutine
C            TRFIND.
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR,LEND,LNEW = Triangulation data structure
C                             created by Subroutine TRMESH.
C                             Node N0 must be included in
C                             the triangulation.
C
C On output:
C
C       LIST,LPTR,LEND,LNEW = Data structure updated with
C                             the addition of node KK as the
C                             last entry.  The updated
C                             triangulation contains no
C                             boundary nodes.
C
C Module required by COVSPH:  INSERT
C
C***********************************************************
C
      INTEGER K, LP, LSAV, NEXT, NST
C
C Local parameters:
C
C K =     Local copy of KK
C LP =    LIST pointer
C LSAV =  LIST pointer
C NEXT =  Boundary node visible from K
C NST =   Local copy of N0
C
      K = KK
      NST = N0
C
C Traverse the boundary in clockwise order, inserting K as
C   the first neighbor of each boundary node, and converting
C   the boundary node to an interior node.
C
      NEXT = NST
    1 LP = LEND(NEXT)
        CALL INSERT (K,LP, LIST,LPTR,LNEW )
        NEXT = -LIST(LP)
        LIST(LP) = NEXT
        IF (NEXT .NE. NST) GO TO 1
C
C Traverse the boundary again, adding each node to K's
C   adjacency lst.
C
      LSAV = LNEW
    2 LP = LEND(NEXT)
        LIST(LNEW) = NEXT
        LPTR(LNEW) = LNEW + 1
        LNEW = LNEW + 1
        NEXT = LIST(LP)
        IF (NEXT .NE. NST) GO TO 2
C
      LPTR(LNEW-1) = LSAV
      LEND(K) = LNEW - 1
      RETURN
      END
      SUBROUTINE CRLIST (N,NCOL,X,Y,Z,LIST,LEND, LPTR,LNEW,
     .                   LTRI, LISTC,NB,XC,YC,ZC,RC,IER)
      INTEGER  N, NCOL, LIST(*), LEND(N), LPTR(*), LNEW,
     .         LTRI(6,NCOL), LISTC(*), NB, IER
      DOUBLE PRECISION X(N), Y(N), Z(N), XC(*), YC(*), ZC(*), RC(*)
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/05/98
C
C   Given a Delaunay triangulation of nodes on the surface
C of the unit sphere, this subroutine returns the set of
C triangle circumcenters corresponding to Voronoi vertices,
C along with the circumradii and a lst of triangle indexes
C LISTC stored in one-to-one correspondence with LIST/LPTR
C entries.
C
C   A triangle circumcenter is the point (unit vector) lying
C at the same angular distance from the three vertices and
C contained in the same hemisphere as the vertices.  (Note
C that the negative of a circumcenter is also equidistant
C from the vertices.)  If the triangulation covers the sur-
C face, the Voronoi vertices are the circumcenters of the
C triangles in the Delaunay triangulation.  LPTR, LEND, and
C LNEW are not altered in this case.
C
C   On the other hand, if the nodes are contained in a sin-
C gle hemisphere, the triangulation is implicitly extended
C to the entire surface by adding pseudo-arcs (of length
C greater than 180 degrees) between boundary nodes forming
C pseudo-triangles whose 'circumcenters' are included in the
C lst.  This extension to the triangulation actually con-
C sists of a triangulation of the set of boundary nodes in
C which the swap test is reversed (a non-empty circumcircle
C test).  The negative circumcenters are stored as the
C pseudo-triangle 'circumcenters'.  LISTC, LPTR, LEND, and
C LNEW contain a data structure corresponding to the ex-
C tended triangulation (Voronoi diagram), but LIST is not
C altered in this case.  Thus, if it is necessary to retain
C the original (unextended) triangulation data structure,
C copies of LPTR and LNEW must be saved before calling this
C routine.
C
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C           Note that, if N = 3, there are only two Voronoi
C           vertices separated by 180 degrees, and the
C           Voronoi regions are not well defined.
C
C       NCOL = Number of columns reserved for LTRI.  This
C              must be at least NB-2, where NB is the number
C              of boundary nodes.
C
C       X,Y,Z = Arrays of length N containing the Cartesian
C               coordinates of the nodes (unit vectors).
C
C       LIST = Integer array containing the set of adjacency
C              lsts.  Refer to Subroutine TRMESH.
C
C       LEND = Set of pointers to ends of adjacency lsts.
C              Refer to Subroutine TRMESH.
C
C The above parameters are not altered by this routine.
C
C       LPTR = Array of pointers associated with LIST.  Re-
C              fer to Subroutine TRMESH.
C
C       LNEW = Pointer to the first empty location in LIST
C              and LPTR (lst length plus one).
C
C       LTRI = Integer work space array dimensioned 6 by
C              NCOL, or unused dummy parameter if NB = 0.
C
C       LISTC = Integer array of length at least 3*NT, where
C               NT = 2*N-4 is the number of triangles in the
C               triangulation (after extending it to cover
C               the entire surface if necessary).
C
C       XC,YC,ZC,RC = Arrays of length NT = 2*N-4.
C
C On output:
C
C       LPTR = Array of pointers associated with LISTC:
C              updated for the addition of pseudo-triangles
C              if the original triangulation contains
C              boundary nodes (NB > 0).
C
C       LNEW = Pointer to the first empty location in LISTC
C              and LPTR (lst length plus one).  LNEW is not
C              altered if NB = 0.
C
C       LTRI = Triangle lst whose first NB-2 columns con-
C              tain the indexes of a clockwise-ordered
C              sequence of vertices (first three rows)
C              followed by the LTRI column indexes of the
C              triangles opposite the vertices (or 0
C              denoting the exterior region) in the last
C              three rows.  This array is not generally of
C              any use.
C
C       LISTC = Array containing triangle indexes (indexes
C               to XC, YC, ZC, and RC) stored in 1-1 corres-
C               pondence with LIST/LPTR entries (or entries
C               that would be stored in LIST for the
C               extended triangulation):  the index of tri-
C               angle (N1,N2,N3) is stored in LISTC(K),
C               LISTC(L), and LISTC(M), where LIST(K),
C               LIST(L), and LIST(M) are the indexes of N2
C               as a neighbor of N1, N3 as a neighbor of N2,
C               and N1 as a neighbor of N3.  The Voronoi
C               region associated with a node is defined by
C               the CCW-ordered sequence of circumcenters in
C               one-to-one correspondence with its adjacency
C               lst (in the extended triangulation).
C
C       NB = Number of boundary nodes unless IER = 1.
C
C       XC,YC,ZC = Arrays containing the Cartesian coordi-
C                  nates of the triangle circumcenters
C                  (Voronoi vertices).  XC(I)**2 + YC(I)**2
C                  + ZC(I)**2 = 1.  The first NB-2 entries
C                  correspond to pseudo-triangles if NB > 0.
C
C       RC = Array containing circumradii (the arc lengths
C            or angles between the circumcenters and associ-
C            ated triangle vertices) in 1-1 correspondence
C            with circumcenters.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N < 3.
C             IER = 2 if NCOL < NB-2.
C             IER = 3 if a triangle is degenerate (has ver-
C                     tices lying on a common geodesic).
C
C Modules required by CRLIST:  CIRCUM, LSTPTR, SWPTST
C
C Intrinsic functions called by CRLIST:  ABS, ACOS
C
C***********************************************************
C
      INTEGER LSTPTR
      INTEGER I1, I2, I3, I4, IERR, KT, KT1, KT2, KT11,
     .        KT12, KT21, KT22, LP, LPL, LPN, N0, N1, N2,
     .        N3, N4, NM2, NN, NT
      LOGICAL SWPTST
      LOGICAL SWP
      DOUBLE PRECISION    C(3), T, V1(3), V2(3), V3(3)
C
C Local parameters:
C
C C =         Circumcenter returned by Subroutine CIRCUM
C I1,I2,I3 =  Permutation of (1,2,3):  LTRI row indexes
C I4 =        LTRI row index in the range 1 to 3
C IERR =      Error flag for calls to CIRCUM
C KT =        Triangle index
C KT1,KT2 =   Indexes of a pair of adjacent pseudo-triangles
C KT11,KT12 = Indexes of the pseudo-triangles opposite N1
C               and N2 as vertices of KT1
C KT21,KT22 = Indexes of the pseudo-triangles opposite N1
C               and N2 as vertices of KT2
C LP,LPN =    LIST pointers
C LPL =       LIST pointer of the last neighbor of N1
C N0 =        Index of the first boundary node (initial
C               value of N1) in the loop on boundary nodes
C               used to store the pseudo-triangle indexes
C               in LISTC
C N1,N2,N3 =  Nodal indexes defining a triangle (CCW order)
C               or pseudo-triangle (clockwise order)
C N4 =        Index of the node opposite N2 -> N1
C NM2 =       N-2
C NN =        Local copy of N
C NT =        Number of pseudo-triangles:  NB-2
C SWP =       Logical variable set to TRUE in each optimiza-
C               tion loop (loop on pseudo-arcs) iff a swap
C               is performed
C V1,V2,V3 =  Vertices of triangle KT = (N1,N2,N3) sent to
C               Subroutine CIRCUM
C
      NN = N
      NB = 0
      NT = 0
      IF (NN .LT. 3) GO TO 21
C
C Search for a boundary node N1.
C
      DO 1 N1 = 1,NN
        LP = LEND(N1)
        IF (LIST(LP) .LT. 0) GO TO 2
    1   CONTINUE
C
C The triangulation already covers the sphere.
C
      GO TO 9
C
C There are NB .GE. 3 boundary nodes.  Add NB-2 pseudo-
C   triangles (N1,N2,N3) by connecting N3 to the NB-3
C   boundary nodes to which it is not already adjacent.
C
C   Set N3 and N2 to the first and last neighbors,
C     respectively, of N1.
C
    2 N2 = -LIST(LP)
      LP = LPTR(LP)
      N3 = LIST(LP)
C
C   Loop on boundary arcs N1 -> N2 in clockwise order,
C     storing triangles (N1,N2,N3) in column NT of LTRI
C     along with the indexes of the triangles opposite
C     the vertices.
C
    3 NT = NT + 1
        IF (NT .LE. NCOL) THEN
          LTRI(1,NT) = N1
          LTRI(2,NT) = N2
          LTRI(3,NT) = N3
          LTRI(4,NT) = NT + 1
          LTRI(5,NT) = NT - 1
          LTRI(6,NT) = 0
        ENDIF
        N1 = N2
        LP = LEND(N1)
        N2 = -LIST(LP)
        IF (N2 .NE. N3) GO TO 3
C
      NB = NT + 2
      IF (NCOL .LT. NT) GO TO 22
      LTRI(4,NT) = 0
      IF (NT .EQ. 1) GO TO 7
C
C Optimize the exterior triangulation (set of pseudo-
C   triangles) by applying swaps to the pseudo-arcs N1-N2
C   (pairs of adjacent pseudo-triangles KT1 and KT2 > KT1).
C   The loop on pseudo-arcs is repeated until no swaps are
C   performed.
C
    4 SWP = .FALSE.
      DO 6 KT1 = 1,NT-1
        DO 5 I3 = 1,3
          KT2 = LTRI(I3+3,KT1)
          IF (KT2 .LE. KT1) GO TO 5
C
C   The LTRI row indexes (I1,I2,I3) of triangle KT1 =
C     (N1,N2,N3) are a cyclical permutation of (1,2,3).
C
          IF (I3 .EQ. 1) THEN
            I1 = 2
            I2 = 3
          ELSEIF (I3 .EQ. 2) THEN
            I1 = 3
            I2 = 1
          ELSE
            I1 = 1
            I2 = 2
          ENDIF
          N1 = LTRI(I1,KT1)
          N2 = LTRI(I2,KT1)
          N3 = LTRI(I3,KT1)
C
C   KT2 = (N2,N1,N4) for N4 = LTRI(I,KT2), where
C     LTRI(I+3,KT2) = KT1.
C
          IF (LTRI(4,KT2) .EQ. KT1) THEN
            I4 = 1
          ELSEIF (LTRI(5,KT2) .EQ. KT1) THEN
            I4 = 2
          ELSE
            I4 = 3
          ENDIF
          N4 = LTRI(I4,KT2)
C
C   The empty circumcircle test is reversed for the pseudo-
C     triangles.  The reversal is implicit in the clockwise
C     ordering of the vertices.
C
          IF ( .NOT. SWPTST(N1,N2,N3,N4,X,Y,Z) ) GO TO 5
C
C   Swap arc N1-N2 for N3-N4.  KTij is the triangle opposite
C     Nj as a vertex of KTi.
C
          SWP = .TRUE.
          KT11 = LTRI(I1+3,KT1)
          KT12 = LTRI(I2+3,KT1)
          IF (I4 .EQ. 1) THEN
            I2 = 2
            I1 = 3
          ELSEIF (I4 .EQ. 2) THEN
            I2 = 3
            I1 = 1
          ELSE
            I2 = 1
            I1 = 2
          ENDIF
          KT21 = LTRI(I1+3,KT2)
          KT22 = LTRI(I2+3,KT2)
          LTRI(1,KT1) = N4
          LTRI(2,KT1) = N3
          LTRI(3,KT1) = N1
          LTRI(4,KT1) = KT12
          LTRI(5,KT1) = KT22
          LTRI(6,KT1) = KT2
          LTRI(1,KT2) = N3
          LTRI(2,KT2) = N4
          LTRI(3,KT2) = N2
          LTRI(4,KT2) = KT21
          LTRI(5,KT2) = KT11
          LTRI(6,KT2) = KT1
C
C   Correct the KT11 and KT22 entries that changed.
C
          IF (KT11 .NE. 0) THEN
            I4 = 4
            IF (LTRI(4,KT11) .NE. KT1) THEN
              I4 = 5
              IF (LTRI(5,KT11) .NE. KT1) I4 = 6
            ENDIF
            LTRI(I4,KT11) = KT2
          ENDIF
          IF (KT22 .NE. 0) THEN
            I4 = 4
            IF (LTRI(4,KT22) .NE. KT2) THEN
              I4 = 5
              IF (LTRI(5,KT22) .NE. KT2) I4 = 6
            ENDIF
            LTRI(I4,KT22) = KT1
          ENDIF
    5     CONTINUE
    6   CONTINUE
      IF (SWP) GO TO 4
C
C Compute and store the negative circumcenters and radii of
C   the pseudo-triangles in the first NT positions.
C
    7 DO 8 KT = 1,NT
        N1 = LTRI(1,KT)
        N2 = LTRI(2,KT)
        N3 = LTRI(3,KT)
        V1(1) = X(N1)
        V1(2) = Y(N1)
        V1(3) = Z(N1)
        V2(1) = X(N2)
        V2(2) = Y(N2)
        V2(3) = Z(N2)
        V3(1) = X(N3)
        V3(2) = Y(N3)
        V3(3) = Z(N3)
        CALL CIRCUM (V1,V2,V3, C,IERR)
        IF (IERR .NE. 0) GO TO 23
C
C   Store the negative circumcenter and radius (computed
C     from <V1,C>).
C
        XC(KT) = C(1)
        YC(KT) = C(2)
        ZC(KT) = C(3)
        T = V1(1)*C(1) + V1(2)*C(2) + V1(3)*C(3)
        IF (T .LT. -1.0) T = -1.0
        IF (T .GT. 1.0) T = 1.0
        RC(KT) = ACOS(T)
    8   CONTINUE
C
C Compute and store the circumcenters and radii of the
C   actual triangles in positions KT = NT+1, NT+2, ...
C   Also, store the triangle indexes KT in the appropriate
C   LISTC positions.
C
    9 KT = NT
C
C   Loop on nodes N1.
C
      NM2 = NN - 2
      DO 12 N1 = 1,NM2
        LPL = LEND(N1)
        LP = LPL
        N3 = LIST(LP)
C
C   Loop on adjacent neighbors N2,N3 of N1 for which N2 > N1
C     and N3 > N1.
C
   10   LP = LPTR(LP)
          N2 = N3
          N3 = ABS(LIST(LP))
          IF (N2 .LE. N1  .OR.  N3 .LE. N1) GO TO 11
          KT = KT + 1
C
C   Compute the circumcenter C of triangle KT = (N1,N2,N3).
C
          V1(1) = X(N1)
          V1(2) = Y(N1)
          V1(3) = Z(N1)
          V2(1) = X(N2)
          V2(2) = Y(N2)
          V2(3) = Z(N2)
          V3(1) = X(N3)
          V3(2) = Y(N3)
          V3(3) = Z(N3)
          CALL CIRCUM (V1,V2,V3, C,IERR)
          IF (IERR .NE. 0) GO TO 23
C
C   Store the circumcenter, radius and triangle index.
C
          XC(KT) = C(1)
          YC(KT) = C(2)
          ZC(KT) = C(3)
          T = V1(1)*C(1) + V1(2)*C(2) + V1(3)*C(3)
          IF (T .LT. -1.0) T = -1.0
          IF (T .GT. 1.0) T = 1.0
          RC(KT) = ACOS(T)
C
C   Store KT in LISTC(LPN), where Abs(LIST(LPN)) is the
C     index of N2 as a neighbor of N1, N3 as a neighbor
C     of N2, and N1 as a neighbor of N3.
C
          LPN = LSTPTR(LPL,N2,LIST,LPTR)
          LISTC(LPN) = KT
          LPN = LSTPTR(LEND(N2),N3,LIST,LPTR)
          LISTC(LPN) = KT
          LPN = LSTPTR(LEND(N3),N1,LIST,LPTR)
          LISTC(LPN) = KT
   11     IF (LP .NE. LPL) GO TO 10
   12   CONTINUE
      IF (NT .EQ. 0) GO TO 20
C
C Store the first NT triangle indexes in LISTC.
C
C   Find a boundary triangle KT1 = (N1,N2,N3) with a
C     boundary arc opposite N3.
C
      KT1 = 0
   13 KT1 = KT1 + 1
      IF (LTRI(4,KT1) .EQ. 0) THEN
        I1 = 2
        I2 = 3
        I3 = 1
        GO TO 14
      ELSEIF (LTRI(5,KT1) .EQ. 0) THEN
        I1 = 3
        I2 = 1
        I3 = 2
        GO TO 14
      ELSEIF (LTRI(6,KT1) .EQ. 0) THEN
        I1 = 1
        I2 = 2
        I3 = 3
        GO TO 14
      ENDIF
      GO TO 13
   14 N1 = LTRI(I1,KT1)
      N0 = N1
C
C   Loop on boundary nodes N1 in CCW order, storing the
C     indexes of the clockwise-ordered sequence of triangles
C     that contain N1.  The first triangle overwrites the
C     last neighbor position, and the remaining triangles,
C     if any, are appended to N1's adjacency lst.
C
C   A pointer to the first neighbor of N1 is saved in LPN.
C
   15 LP = LEND(N1)
      LPN = LPTR(LP)
      LISTC(LP) = KT1
C
C   Loop on triangles KT2 containing N1.
C
   16 KT2 = LTRI(I2+3,KT1)
      IF (KT2 .NE. 0) THEN
C
C   Append KT2 to N1's triangle lst.
C
        LPTR(LP) = LNEW
        LP = LNEW
        LISTC(LP) = KT2
        LNEW = LNEW + 1
C
C   Set KT1 to KT2 and update (I1,I2,I3) such that
C     LTRI(I1,KT1) = N1.
C
        KT1 = KT2
        IF (LTRI(1,KT1) .EQ. N1) THEN
          I1 = 1
          I2 = 2
          I3 = 3
        ELSEIF (LTRI(2,KT1) .EQ. N1) THEN
          I1 = 2
          I2 = 3
          I3 = 1
        ELSE
          I1 = 3
          I2 = 1
          I3 = 2
        ENDIF
        GO TO 16
      ENDIF
C
C   Store the saved first-triangle pointer in LPTR(LP), set
C     N1 to the next boundary node, test for termination,
C     and permute the indexes:  the last triangle containing
C     a boundary node is the first triangle containing the
C     next boundary node.
C
      LPTR(LP) = LPN
      N1 = LTRI(I3,KT1)
      IF (N1 .NE. N0) THEN
        I4 = I3
        I3 = I2
        I2 = I1
        I1 = I4
        GO TO 15
      ENDIF
C
C No errors encountered.
C
   20 IER = 0
      RETURN
C
C N < 3.
C
   21 IER = 1
      RETURN
C
C Insufficient space reserved for LTRI.
C
   22 IER = 2
      RETURN
C
C Error flag returned by CIRCUM: KT indexes a null triangle.
C
   23 IER = 3
      RETURN
      END
      SUBROUTINE DELARC (N,IO1,IO2, LIST,LPTR,LEND,
     .                   LNEW, IER)
      INTEGER N, IO1, IO2, LIST(*), LPTR(*), LEND(N), LNEW,
     .        IER
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/17/96
C
C   This subroutine deletes a boundary arc from a triangula-
C tion.  It may be used to remove a null triangle from the
C convex hull boundary.  Note, however, that if the union of
C triangles is rendered nonconvex, Subroutines DELNOD, EDGE,
C and TRFIND (and hence ADDNOD) may fail.  Also, Function
C NEARND should not be called following an arc deletion.
C
C   This routine is identical to the similarly named routine
C in TRIPACK.
C
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 4.
C
C       IO1,IO2 = Indexes (in the range 1 to N) of a pair of
C                 adjacent boundary nodes defining the arc
C                 to be removed.
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR,LEND,LNEW = Triangulation data structure
C                             created by Subroutine TRMESH.
C
C On output:
C
C       LIST,LPTR,LEND,LNEW = Data structure updated with
C                             the removal of arc IO1-IO2
C                             unless IER > 0.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N, IO1, or IO2 is outside its valid
C                     range, or IO1 = IO2.
C             IER = 2 if IO1-IO2 is not a boundary arc.
C             IER = 3 if the node opposite IO1-IO2 is al-
C                     ready a boundary node, and thus IO1
C                     or IO2 has only two neighbors or a
C                     deletion would result in two triangu-
C                     lations sharing a single node.
C             IER = 4 if one of the nodes is a neighbor of
C                     the other, but not vice versa, imply-
C                     ing an invalid triangulation data
C                     structure.
C
C Module required by DELARC:  DELNB, LSTPTR
C
C Intrinsic function called by DELARC:  ABS
C
C***********************************************************
C
      INTEGER LSTPTR
      INTEGER LP, LPH, LPL, N1, N2, N3
C
C Local parameters:
C
C LP =       LIST pointer
C LPH =      LIST pointer or flag returned by DELNB
C LPL =      Pointer to the last neighbor of N1, N2, or N3
C N1,N2,N3 = Nodal indexes of a triangle such that N1->N2
C              is the directed boundary edge associated
C              with IO1-IO2
C
      N1 = IO1
      N2 = IO2
C
C Test for errors, and set N1->N2 to the directed boundary
C   edge associated with IO1-IO2:  (N1,N2,N3) is a triangle
C   for some N3.
C
      IF (N .LT. 4  .OR.  N1 .LT. 1  .OR.  N1 .GT. N  .OR.
     .    N2 .LT. 1  .OR.  N2 .GT. N  .OR.  N1 .EQ. N2) THEN
        IER = 1
        RETURN
      ENDIF
C
      LPL = LEND(N2)
      IF (-LIST(LPL) .NE. N1) THEN
        N1 = N2
        N2 = IO1
        LPL = LEND(N2)
        IF (-LIST(LPL) .NE. N1) THEN
          IER = 2
          RETURN
        ENDIF
      ENDIF
C
C Set N3 to the node opposite N1->N2 (the second neighbor
C   of N1), and test for error 3 (N3 already a boundary
C   node).
C
      LPL = LEND(N1)
      LP = LPTR(LPL)
      LP = LPTR(LP)
      N3 = ABS(LIST(LP))
      LPL = LEND(N3)
      IF (LIST(LPL) .LE. 0) THEN
        IER = 3
        RETURN
      ENDIF
C
C Delete N2 as a neighbor of N1, making N3 the first
C   neighbor, and test for error 4 (N2 not a neighbor
C   of N1).  Note that previously computed pointers may
C   no longer be valid following the call to DELNB.
C
      CALL DELNB (N1,N2,N, LIST,LPTR,LEND,LNEW, LPH)
      IF (LPH .LT. 0) THEN
        IER = 4
        RETURN
      ENDIF
C
C Delete N1 as a neighbor of N2, making N3 the new last
C   neighbor.
C
      CALL DELNB (N2,N1,N, LIST,LPTR,LEND,LNEW, LPH)
C
C Make N3 a boundary node with first neighbor N2 and last
C   neighbor N1.
C
      LP = LSTPTR(LEND(N3),N1,LIST,LPTR)
      LEND(N3) = LP
      LIST(LP) = -N1
C
C No errors encountered.
C
      IER = 0
      RETURN
      END
      SUBROUTINE DELNB (N0,NB,N, LIST,LPTR,LEND,LNEW, LPH)
      INTEGER N0, NB, N, LIST(*), LPTR(*), LEND(N), LNEW,
     .        LPH
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/29/98
C
C   This subroutine deletes a neighbor NB from the adjacency
C lst of node N0 (but N0 is not deleted from the adjacency
C lst of NB) and, if NB is a boundary node, makes N0 a
C boundary node.  For pointer (LIST index) LPH to NB as a
C neighbor of N0, the empty LIST,LPTR location LPH is filled
C in with the values at LNEW-1, pointer LNEW-1 (in LPTR and
C possibly in LEND) is changed to LPH, and LNEW is decremen-
C ted.  This requires a search of LEND and LPTR entailing an
C expected operation count of O(N).
C
C   This routine is identical to the similarly named routine
C in TRIPACK.
C
C
C On input:
C
C       N0,NB = Indexes, in the range 1 to N, of a pair of
C               nodes such that NB is a neighbor of N0.
C               (N0 need not be a neighbor of NB.)
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR,LEND,LNEW = Data structure defining the
C                             triangulation.
C
C On output:
C
C       LIST,LPTR,LEND,LNEW = Data structure updated with
C                             the removal of NB from the ad-
C                             jacency lst of N0 unless
C                             LPH < 0.
C
C       LPH = List pointer to the hole (NB as a neighbor of
C             N0) filled in by the values at LNEW-1 or error
C             indicator:
C             LPH > 0 if no errors were encountered.
C             LPH = -1 if N0, NB, or N is outside its valid
C                      range.
C             LPH = -2 if NB is not a neighbor of N0.
C
C Modules required by DELNB:  None
C
C Intrinsic function called by DELNB:  ABS
C
C***********************************************************
C
      INTEGER I, LNW, LP, LPB, LPL, LPP, NN
C
C Local parameters:
C
C I =   DO-loop index
C LNW = LNEW-1 (output value of LNEW)
C LP =  LIST pointer of the last neighbor of NB
C LPB = Pointer to NB as a neighbor of N0
C LPL = Pointer to the last neighbor of N0
C LPP = Pointer to the neighbor of N0 that precedes NB
C NN =  Local copy of N
C
      NN = N
C
C Test for error 1.
C
      IF (N0 .LT. 1  .OR.  N0 .GT. NN  .OR.  NB .LT. 1  .OR.
     .    NB .GT. NN  .OR.  NN .LT. 3) THEN
        LPH = -1
        RETURN
      ENDIF
C
C   Find pointers to neighbors of N0:
C
C     LPL points to the last neighbor,
C     LPP points to the neighbor NP preceding NB, and
C     LPB points to NB.
C
      LPL = LEND(N0)
      LPP = LPL
      LPB = LPTR(LPP)
    1 IF (LIST(LPB) .EQ. NB) GO TO 2
        LPP = LPB
        LPB = LPTR(LPP)
        IF (LPB .NE. LPL) GO TO 1
C
C   Test for error 2 (NB not found).
C
      IF (ABS(LIST(LPB)) .NE. NB) THEN
        LPH = -2
        RETURN
      ENDIF
C
C   NB is the last neighbor of N0.  Make NP the new last
C     neighbor and, if NB is a boundary node, then make N0
C     a boundary node.
C
      LEND(N0) = LPP
      LP = LEND(NB)
      IF (LIST(LP) .LT. 0) LIST(LPP) = -LIST(LPP)
      GO TO 3
C
C   NB is not the last neighbor of N0.  If NB is a boundary
C     node and N0 is not, then make N0 a boundary node with
C     last neighbor NP.
C
    2 LP = LEND(NB)
      IF (LIST(LP) .LT. 0  .AND.  LIST(LPL) .GT. 0) THEN
        LEND(N0) = LPP
        LIST(LPP) = -LIST(LPP)
      ENDIF
C
C   Update LPTR so that the neighbor following NB now fol-
C     lows NP, and fill in the hole at location LPB.
C
    3 LPTR(LPP) = LPTR(LPB)
      LNW = LNEW-1
      LIST(LPB) = LIST(LNW)
      LPTR(LPB) = LPTR(LNW)
      DO 4 I = NN,1,-1
        IF (LEND(I) .EQ. LNW) THEN
          LEND(I) = LPB
          GO TO 5
        ENDIF
    4   CONTINUE
C
    5 DO 6 I = 1,LNW-1
        IF (LPTR(I) .EQ. LNW) THEN
          LPTR(I) = LPB
        ENDIF
    6   CONTINUE
C
C No errors encountered.
C
      LNEW = LNW
      LPH = LPB
      RETURN
      END
      SUBROUTINE DELNOD (K, N,X,Y,Z,LIST,LPTR,LEND,LNEW,LWK,
     .                   IWK, IER)
      INTEGER K, N, LIST(*), LPTR(*), LEND(*), LNEW, LWK,
     .        IWK(2,*), IER
      DOUBLE PRECISION    X(*), Y(*), Z(*)
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/30/99
C
C   This subroutine deletes node K (along with all arcs
C incident on node K) from a triangulation of N nodes on the
C unit sphere, and inserts arcs as necessary to produce a
C triangulation of the remaining N-1 nodes.  If a Delaunay
C triangulation is input, a Delaunay triangulation will
C result, and thus, DELNOD reverses the effect of a call to
C Subroutine ADDNOD.
C
C
C On input:
C
C       K = Index (for X, Y, and Z) of the node to be
C           deleted.  1 .LE. K .LE. N.
C
C K is not altered by this routine.
C
C       N = Number of nodes in the triangulation on input.
C           N .GE. 4.  Note that N will be decremented
C           following the deletion.
C
C       X,Y,Z = Arrays of length N containing the Cartesian
C               coordinates of the nodes in the triangula-
C               tion.
C
C       LIST,LPTR,LEND,LNEW = Data structure defining the
C                             triangulation.  Refer to Sub-
C                             routine TRMESH.
C
C       LWK = Number of columns reserved for IWK.  LWK must
C             be at least NNB-3, where NNB is the number of
C             neighbors of node K, including an extra
C             pseudo-node if K is a boundary node.
C
C       IWK = Integer work array dimensioned 2 by LWK (or
C             array of length .GE. 2*LWK).
C
C On output:
C
C       N = Number of nodes in the triangulation on output.
C           The input value is decremented unless 1 .LE. IER
C           .LE. 4.
C
C       X,Y,Z = Updated arrays containing nodal coordinates
C               (with elements K+1,...,N+1 shifted up one
C               position, thus overwriting element K) unless
C               1 .LE. IER .LE. 4.
C
C       LIST,LPTR,LEND,LNEW = Updated triangulation data
C                             structure reflecting the dele-
C                             tion unless 1 .LE. IER .LE. 4.
C                             Note that the data structure
C                             may have been altered if IER >
C                             3.
C
C       LWK = Number of IWK columns required unless IER = 1
C             or IER = 3.
C
C       IWK = Indexes of the endpoints of the new arcs added
C             unless LWK = 0 or 1 .LE. IER .LE. 4.  (Arcs
C             are associated with columns, or pairs of
C             adjacent elements if IWK is declared as a
C             singly-subscripted array.)
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if K or N is outside its valid range
C                     or LWK < 0 on input.
C             IER = 2 if more space is required in IWK.
C                     Refer to LWK.
C             IER = 3 if the triangulation data structure is
C                     invalid on input.
C             IER = 4 if K indexes an interior node with
C                     four or more neighbors, none of which
C                     can be swapped out due to collineari-
C                     ty, and K cannot therefore be deleted.
C             IER = 5 if an error flag (other than IER = 1)
C                     was returned by OPTIM.  An error
C                     message is written to the standard
C                     output unit in this case.
C             IER = 6 if error flag 1 was returned by OPTIM.
C                     This is not necessarily an error, but
C                     the arcs may not be optimal.
C
C   Note that the deletion may result in all remaining nodes
C being collinear.  This situation is not flagged.
C
C Modules required by DELNOD:  DELNB, LEFT, LSTPTR, NBCNT,
C                                OPTIM, SWAP, SWPTST
C
C Intrinsic function called by DELNOD:  ABS
C
C***********************************************************
C
      INTEGER LSTPTR, NBCNT
      INTEGER I, IERR, IWL, J, LNW, LP, LP21, LPF, LPH, LPL,
     .        LPL2, LPN, LWKL, N1, N2, NFRST, NIT, NL, NN,
     .        NNB, NR
      LOGICAL LEFT
      LOGICAL BDRY
      DOUBLE PRECISION    X1, X2, XL, XR, Y1, Y2, YL, YR, Z1, Z2, ZL, ZR
C
C Local parameters:
C
C BDRY =    Logical variable with value TRUE iff N1 is a
C             boundary node
C I,J =     DO-loop indexes
C IERR =    Error flag returned by OPTIM
C IWL =     Number of IWK columns containing arcs
C LNW =     Local copy of LNEW
C LP =      LIST pointer
C LP21 =    LIST pointer returned by SWAP
C LPF,LPL = Pointers to the first and last neighbors of N1
C LPH =     Pointer (or flag) returned by DELNB
C LPL2 =    Pointer to the last neighbor of N2
C LPN =     Pointer to a neighbor of N1
C LWKL =    Input value of LWK
C N1 =      Local copy of K
C N2 =      Neighbor of N1
C NFRST =   First neighbor of N1:  LIST(LPF)
C NIT =     Number of iterations in OPTIM
C NR,NL =   Neighbors of N1 preceding (to the right of) and
C             following (to the left of) N2, respectively
C NN =      Number of nodes in the triangulation
C NNB =     Number of neighbors of N1 (including a pseudo-
C             node representing the boundary if N1 is a
C             boundary node)
C X1,Y1,Z1 = Coordinates of N1
C X2,Y2,Z2 = Coordinates of N2
C XL,YL,ZL = Coordinates of NL
C XR,YR,ZR = Coordinates of NR
C
C
C Set N1 to K and NNB to the number of neighbors of N1 (plus
C   one if N1 is a boundary node), and test for errors.  LPF
C   and LPL are LIST indexes of the first and last neighbors
C   of N1, IWL is the number of IWK columns containing arcs,
C   and BDRY is TRUE iff N1 is a boundary node.
C
      N1 = K
      NN = N
      IF (N1 .LT. 1  .OR.  N1 .GT. NN  .OR.  NN .LT. 4  .OR.
     .    LWK .LT. 0) GO TO 21
      LPL = LEND(N1)
      LPF = LPTR(LPL)
      NNB = NBCNT(LPL,LPTR)
      BDRY = LIST(LPL) .LT. 0
      IF (BDRY) NNB = NNB + 1
      IF (NNB .LT. 3) GO TO 23
      LWKL = LWK
      LWK = NNB - 3
      IF (LWKL .LT. LWK) GO TO 22
      IWL = 0
      IF (NNB .EQ. 3) GO TO 3
C
C Initialize for loop on arcs N1-N2 for neighbors N2 of N1,
C   beginning with the second neighbor.  NR and NL are the
C   neighbors preceding and following N2, respectively, and
C   LP indexes NL.  The loop is exited when all possible
C   swaps have been applied to arcs incident on N1.
C
      X1 = X(N1)
      Y1 = Y(N1)
      Z1 = Z(N1)
      NFRST = LIST(LPF)
      NR = NFRST
      XR = X(NR)
      YR = Y(NR)
      ZR = Z(NR)
      LP = LPTR(LPF)
      N2 = LIST(LP)
      X2 = X(N2)
      Y2 = Y(N2)
      Z2 = Z(N2)
      LP = LPTR(LP)
C
C Top of loop:  set NL to the neighbor following N2.
C
    1 NL = ABS(LIST(LP))
      IF (NL .EQ. NFRST  .AND.  BDRY) GO TO 3
      XL = X(NL)
      YL = Y(NL)
      ZL = Z(NL)
C
C   Test for a convex quadrilateral.  To avoid an incorrect
C     test caused by collinearity, use the fact that if N1
C     is a boundary node, then N1 LEFT NR->NL and if N2 is
C     a boundary node, then N2 LEFT NL->NR.
C
      LPL2 = LEND(N2)
      IF ( .NOT. ((BDRY  .OR.  LEFT(XR,YR,ZR,XL,YL,ZL,X1,Y1,
     .      Z1))  .AND.  (LIST(LPL2) .LT. 0  .OR.
     .      LEFT(XL,YL,ZL,XR,YR,ZR,X2,Y2,Z2))) ) THEN
C
C   Nonconvex quadrilateral -- no swap is possible.
C
        NR = N2
        XR = X2
        YR = Y2
        ZR = Z2
        GO TO 2
      ENDIF
C
C   The quadrilateral defined by adjacent triangles
C     (N1,N2,NL) and (N2,N1,NR) is convex.  Swap in
C     NL-NR and store it in IWK unless NL and NR are
C     already adjacent, in which case the swap is not
C     possible.  Indexes larger than N1 must be decremented
C     since N1 will be deleted from X, Y, and Z.
C
      CALL SWAP (NL,NR,N1,N2, LIST,LPTR,LEND, LP21)
      IF (LP21 .EQ. 0) THEN
        NR = N2
        XR = X2
        YR = Y2
        ZR = Z2
        GO TO 2
      ENDIF
      IWL = IWL + 1
      IF (NL .LE. N1) THEN
        IWK(1,IWL) = NL
      ELSE
        IWK(1,IWL) = NL - 1
      ENDIF
      IF (NR .LE. N1) THEN
        IWK(2,IWL) = NR
      ELSE
        IWK(2,IWL) = NR - 1
      ENDIF
C
C   Recompute the LIST indexes and NFRST, and decrement NNB.
C
      LPL = LEND(N1)
      NNB = NNB - 1
      IF (NNB .EQ. 3) GO TO 3
      LPF = LPTR(LPL)
      NFRST = LIST(LPF)
      LP = LSTPTR(LPL,NL,LIST,LPTR)
      IF (NR .EQ. NFRST) GO TO 2
C
C   NR is not the first neighbor of N1.
C     Back up and test N1-NR for a swap again:  Set N2 to
C     NR and NR to the previous neighbor of N1 -- the
C     neighbor of NR which follows N1.  LP21 points to NL
C     as a neighbor of NR.
C
      N2 = NR
      X2 = XR
      Y2 = YR
      Z2 = ZR
      LP21 = LPTR(LP21)
      LP21 = LPTR(LP21)
      NR = ABS(LIST(LP21))
      XR = X(NR)
      YR = Y(NR)
      ZR = Z(NR)
      GO TO 1
C
C   Bottom of loop -- test for termination of loop.
C
    2 IF (N2 .EQ. NFRST) GO TO 3
      N2 = NL
      X2 = XL
      Y2 = YL
      Z2 = ZL
      LP = LPTR(LP)
      GO TO 1
C
C Delete N1 and all its incident arcs.  If N1 is an interior
C   node and either NNB > 3 or NNB = 3 and N2 LEFT NR->NL,
C   then N1 must be separated from its neighbors by a plane
C   containing the origin -- its removal reverses the effect
C   of a call to COVSPH, and all its neighbors become
C   boundary nodes.  This is achieved by treating it as if
C   it were a boundary node (setting BDRY to TRUE, changing
C   a sign in LIST, and incrementing NNB).
C
    3 IF (.NOT. BDRY) THEN
        IF (NNB .GT. 3) THEN
          BDRY = .TRUE.
        ELSE
          LPF = LPTR(LPL)
          NR = LIST(LPF)
          LP = LPTR(LPF)
          N2 = LIST(LP)
          NL = LIST(LPL)
          BDRY = LEFT(X(NR),Y(NR),Z(NR),X(NL),Y(NL),Z(NL),
     .                X(N2),Y(N2),Z(N2))
        ENDIF
        IF (BDRY) THEN
C
C   IF a boundary node already exists, then N1 and its
C     neighbors cannot be converted to boundary nodes.
C     (They must be collinear.)  This is a problem if
C     NNB > 3.
C
          DO 4 I = 1,NN
            IF (LIST(LEND(I)) .LT. 0) THEN
              BDRY = .FALSE.
              GO TO 5
            ENDIF
    4       CONTINUE
          LIST(LPL) = -LIST(LPL)
          NNB = NNB + 1
        ENDIF
      ENDIF
    5 IF (.NOT. BDRY  .AND.  NNB .GT. 3) GO TO 24
C
C Initialize for loop on neighbors.  LPL points to the last
C   neighbor of N1.  LNEW is stored in local variable LNW.
C
      LP = LPL
      LNW = LNEW
C
C Loop on neighbors N2 of N1, beginning with the first.
C
    6 LP = LPTR(LP)
        N2 = ABS(LIST(LP))
        CALL DELNB (N2,N1,N, LIST,LPTR,LEND,LNW, LPH)
        IF (LPH .LT. 0) GO TO 23
C
C   LP and LPL may require alteration.
C
        IF (LPL .EQ. LNW) LPL = LPH
        IF (LP .EQ. LNW) LP = LPH
        IF (LP .NE. LPL) GO TO 6
C
C Delete N1 from X, Y, Z, and LEND, and remove its adjacency
C   lst from LIST and LPTR.  LIST entries (nodal indexes)
C   which are larger than N1 must be decremented.
C
      NN = NN - 1
      IF (N1 .GT. NN) GO TO 9
      DO 7 I = N1,NN
        X(I) = X(I+1)
        Y(I) = Y(I+1)
        Z(I) = Z(I+1)
        LEND(I) = LEND(I+1)
    7   CONTINUE
C
      DO 8 I = 1,LNW-1
        IF (LIST(I) .GT. N1) LIST(I) = LIST(I) - 1
        IF (LIST(I) .LT. -N1) LIST(I) = LIST(I) + 1
    8   CONTINUE
C
C   For LPN = first to last neighbors of N1, delete the
C     preceding neighbor (indexed by LP).
C
C   Each empty LIST,LPTR location LP is filled in with the
C     values at LNW-1, and LNW is decremented.  All pointers
C     (including those in LPTR and LEND) with value LNW-1
C     must be changed to LP.
C
C  LPL points to the last neighbor of N1.
C
    9 IF (BDRY) NNB = NNB - 1
      LPN = LPL
      DO 13 J = 1,NNB
        LNW = LNW - 1
        LP = LPN
        LPN = LPTR(LP)
        LIST(LP) = LIST(LNW)
        LPTR(LP) = LPTR(LNW)
        IF (LPTR(LPN) .EQ. LNW) LPTR(LPN) = LP
        IF (LPN .EQ. LNW) LPN = LP
        DO 10 I = NN,1,-1
          IF (LEND(I) .EQ. LNW) THEN
            LEND(I) = LP
            GO TO 11
          ENDIF
   10     CONTINUE
C
   11   DO 12 I = LNW-1,1,-1
          IF (LPTR(I) .EQ. LNW) LPTR(I) = LP
   12     CONTINUE
   13   CONTINUE
C
C Update N and LNEW, and optimize the patch of triangles
C   containing K (on input) by applying swaps to the arcs
C   in IWK.
C
      N = NN
      LNEW = LNW
      IF (IWL .GT. 0) THEN
        NIT = 4*IWL
        CALL OPTIM (X,Y,Z,IWL, LIST,LPTR,LEND,NIT,IWK, IERR)
        IF (IERR .NE. 0  .AND.  IERR .NE. 1) GO TO 25
        IF (IERR .EQ. 1) GO TO 26
      ENDIF
C
C Successful termination.
C
      IER = 0
      RETURN
C
C Invalid input parameter.
C
   21 IER = 1
      RETURN
C
C Insufficient space reserved for IWK.
C
   22 IER = 2
      RETURN
C
C Invalid triangulation data structure.  NNB < 3 on input or
C   N2 is a neighbor of N1 but N1 is not a neighbor of N2.
C
   23 IER = 3
      RETURN
C
C N1 is interior but NNB could not be reduced to 3.
C
   24 IER = 4
      RETURN
C
C Error flag (other than 1) returned by OPTIM.
C
   25 IER = 5
      WRITE (*,100) NIT, IERR
  100 FORMAT (//5X,'*** Error in OPTIM (called from ',
     .        'DELNOD):  NIT = ',I4,', IER = ',I1,' ***'/)
      RETURN
C
C Error flag 1 returned by OPTIM.
C
   26 IER = 6
      RETURN
      END
      SUBROUTINE EDGE (IN1,IN2,X,Y,Z, LWK,IWK,LIST,LPTR,
     .                 LEND, IER)
      INTEGER IN1, IN2, LWK, IWK(2,*), LIST(*), LPTR(*),
     .        LEND(*), IER
      DOUBLE PRECISION    X(*), Y(*), Z(*)
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/30/98
C
C   Given a triangulation of N nodes and a pair of nodal
C indexes IN1 and IN2, this routine swaps arcs as necessary
C to force IN1 and IN2 to be adjacent.  Only arcs which
C intersect IN1-IN2 are swapped out.  If a Delaunay triangu-
C lation is input, the resulting triangulation is as close
C as possible to a Delaunay triangulation in the sense that
C all arcs other than IN1-IN2 are locally optimal.
C
C   A sequence of calls to EDGE may be used to force the
C presence of a set of edges defining the boundary of a non-
C convex and/or multiply connected region, or to introduce
C barriers into the triangulation.  Note that Subroutine
C GETNP will not necessarily return closest nodes if the
C triangulation has been constrained by a call to EDGE.
C However, this is appropriate in some applications, such
C as triangle-based interpolation on a nonconvex domain.
C
C
C On input:
C
C       IN1,IN2 = Indexes (of X, Y, and Z) in the range 1 to
C                 N defining a pair of nodes to be connected
C                 by an arc.
C
C       X,Y,Z = Arrays of length N containing the Cartesian
C               coordinates of the nodes.
C
C The above parameters are not altered by this routine.
C
C       LWK = Number of columns reserved for IWK.  This must
C             be at least NI -- the number of arcs that
C             intersect IN1-IN2.  (NI is bounded by N-3.)
C
C       IWK = Integer work array of length at least 2*LWK.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to Subroutine
C                        TRMESH.
C
C On output:
C
C       LWK = Number of arcs which intersect IN1-IN2 (but
C             not more than the input value of LWK) unless
C             IER = 1 or IER = 3.  LWK = 0 if and only if
C             IN1 and IN2 were adjacent (or LWK=0) on input.
C
C       IWK = Array containing the indexes of the endpoints
C             of the new arcs other than IN1-IN2 unless
C             IER > 0 or LWK = 0.  New arcs to the left of
C             IN1->IN2 are stored in the first K-1 columns
C             (left portion of IWK), column K contains
C             zeros, and new arcs to the right of IN1->IN2
C             occupy columns K+1,...,LWK.  (K can be deter-
C             mined by searching IWK for the zeros.)
C
C       LIST,LPTR,LEND = Data structure updated if necessary
C                        to reflect the presence of an arc
C                        connecting IN1 and IN2 unless IER >
C                        0.  The data structure has been
C                        altered if IER >= 4.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if IN1 < 1, IN2 < 1, IN1 = IN2,
C                     or LWK < 0 on input.
C             IER = 2 if more space is required in IWK.
C                     Refer to LWK.
C             IER = 3 if IN1 and IN2 could not be connected
C                     due to either an invalid data struc-
C                     ture or collinear nodes (and floating
C                     point error).
C             IER = 4 if an error flag other than IER = 1
C                     was returned by OPTIM.
C             IER = 5 if error flag 1 was returned by OPTIM.
C                     This is not necessarily an error, but
C                     the arcs other than IN1-IN2 may not
C                     be optimal.
C
C   An error message is written to the standard output unit
C in the case of IER = 3 or IER = 4.
C
C Modules required by EDGE:  LEFT, LSTPTR, OPTIM, SWAP,
C                              SWPTST
C
C Intrinsic function called by EDGE:  ABS
C
C***********************************************************
C
      LOGICAL LEFT
      INTEGER I, IERR, IWC, IWCP1, IWEND, IWF, IWL, LFT, LP,
     .        LP21, LPL, N0, N1, N1FRST, N1LST, N2, NEXT,
     .        NIT, NL, NR
      DOUBLE PRECISION    DP12, DP1L, DP1R, DP2L, DP2R, X0, X1, X2, Y0,
     .        Y1, Y2, Z0, Z1, Z2
C
C Local parameters:
C
C DPij =     Dot product <Ni,Nj>
C I =        DO-loop index and column index for IWK
C IERR =     Error flag returned by Subroutine OPTIM
C IWC =      IWK index between IWF and IWL -- NL->NR is
C              stored in IWK(1,IWC)->IWK(2,IWC)
C IWCP1 =    IWC + 1
C IWEND =    Input or output value of LWK
C IWF =      IWK (column) index of the first (leftmost) arc
C              which intersects IN1->IN2
C IWL =      IWK (column) index of the last (rightmost) are
C              which intersects IN1->IN2
C LFT =      Flag used to determine if a swap results in the
C              new arc intersecting IN1-IN2 -- LFT = 0 iff
C              N0 = IN1, LFT = -1 implies N0 LEFT IN1->IN2,
C              and LFT = 1 implies N0 LEFT IN2->IN1
C LP =       List pointer (index for LIST and LPTR)
C LP21 =     Unused parameter returned by SWAP
C LPL =      Pointer to the last neighbor of IN1 or NL
C N0 =       Neighbor of N1 or node opposite NR->NL
C N1,N2 =    Local copies of IN1 and IN2
C N1FRST =   First neighbor of IN1
C N1LST =    (Signed) last neighbor of IN1
C NEXT =     Node opposite NL->NR
C NIT =      Flag or number of iterations employed by OPTIM
C NL,NR =    Endpoints of an arc which intersects IN1-IN2
C              with NL LEFT IN1->IN2
C X0,Y0,Z0 = Coordinates of N0
C X1,Y1,Z1 = Coordinates of IN1
C X2,Y2,Z2 = Coordinates of IN2
C
C
C Store IN1, IN2, and LWK in local variables and test for
C   errors.
C
      N1 = IN1
      N2 = IN2
      IWEND = LWK
      IF (N1 .LT. 1  .OR.  N2 .LT. 1  .OR.  N1 .EQ. N2  .OR.
     .    IWEND .LT. 0) GO TO 31
C
C Test for N2 as a neighbor of N1.  LPL points to the last
C   neighbor of N1.
C
      LPL = LEND(N1)
      N0 = ABS(LIST(LPL))
      LP = LPL
    1 IF (N0 .EQ. N2) GO TO 30
        LP = LPTR(LP)
        N0 = LIST(LP)
        IF (LP .NE. LPL) GO TO 1
C
C Initialize parameters.
C
      IWL = 0
      NIT = 0
C
C Store the coordinates of N1 and N2.
C
    2 X1 = X(N1)
      Y1 = Y(N1)
      Z1 = Z(N1)
      X2 = X(N2)
      Y2 = Y(N2)
      Z2 = Z(N2)
C
C Set NR and NL to adjacent neighbors of N1 such that
C   NR LEFT N2->N1 and NL LEFT N1->N2,
C   (NR Forward N1->N2 or NL Forward N1->N2), and
C   (NR Forward N2->N1 or NL Forward N2->N1).
C
C   Initialization:  Set N1FRST and N1LST to the first and
C     (signed) last neighbors of N1, respectively, and
C     initialize NL to N1FRST.
C
      LPL = LEND(N1)
      N1LST = LIST(LPL)
      LP = LPTR(LPL)
      N1FRST = LIST(LP)
      NL = N1FRST
      IF (N1LST .LT. 0) GO TO 4
C
C   N1 is an interior node.  Set NL to the first candidate
C     for NR (NL LEFT N2->N1).
C
    3 IF (LEFT(X2,Y2,Z2,X1,Y1,Z1,X(NL),Y(NL),Z(NL))) GO TO 4
        LP = LPTR(LP)
        NL = LIST(LP)
        IF (NL .NE. N1FRST) GO TO 3
C
C   All neighbors of N1 are strictly left of N1->N2.
C
      GO TO 5
C
C   NL = LIST(LP) LEFT N2->N1.  Set NR to NL and NL to the
C     following neighbor of N1.
C
    4 NR = NL
        LP = LPTR(LP)
        NL = ABS(LIST(LP))
        IF (LEFT(X1,Y1,Z1,X2,Y2,Z2,X(NL),Y(NL),Z(NL)) ) THEN
C
C   NL LEFT N1->N2 and NR LEFT N2->N1.  The Forward tests
C     are employed to avoid an error associated with
C     collinear nodes.
C
          DP12 = X1*X2 + Y1*Y2 + Z1*Z2
          DP1L = X1*X(NL) + Y1*Y(NL) + Z1*Z(NL)
          DP2L = X2*X(NL) + Y2*Y(NL) + Z2*Z(NL)
          DP1R = X1*X(NR) + Y1*Y(NR) + Z1*Z(NR)
          DP2R = X2*X(NR) + Y2*Y(NR) + Z2*Z(NR)
          IF ( (DP2L-DP12*DP1L .GE. 0.  .OR.
     .          DP2R-DP12*DP1R .GE. 0.)  .AND.
     .         (DP1L-DP12*DP2L .GE. 0.  .OR.
     .          DP1R-DP12*DP2R .GE. 0.) ) GO TO 6
C
C   NL-NR does not intersect N1-N2.  However, there is
C     another candidate for the first arc if NL lies on
C     the line N1-N2.
C
          IF ( .NOT. LEFT(X2,Y2,Z2,X1,Y1,Z1,X(NL),Y(NL),
     .                    Z(NL)) ) GO TO 5
        ENDIF
C
C   Bottom of loop.
C
        IF (NL .NE. N1FRST) GO TO 4
C
C Either the triangulation is invalid or N1-N2 lies on the
C   convex hull boundary and an edge NR->NL (opposite N1 and
C   intersecting N1-N2) was not found due to floating point
C   error.  Try interchanging N1 and N2 -- NIT > 0 iff this
C   has already been done.
C
    5 IF (NIT .GT. 0) GO TO 33
      NIT = 1
      N1 = N2
      N2 = IN1
      GO TO 2
C
C Store the ordered sequence of intersecting edges NL->NR in
C   IWK(1,IWL)->IWK(2,IWL).
C
    6 IWL = IWL + 1
      IF (IWL .GT. IWEND) GO TO 32
      IWK(1,IWL) = NL
      IWK(2,IWL) = NR
C
C   Set NEXT to the neighbor of NL which follows NR.
C
      LPL = LEND(NL)
      LP = LPTR(LPL)
C
C   Find NR as a neighbor of NL.  The search begins with
C     the first neighbor.
C
    7 IF (LIST(LP) .EQ. NR) GO TO 8
        LP = LPTR(LP)
        IF (LP .NE. LPL) GO TO 7
C
C   NR must be the last neighbor, and NL->NR cannot be a
C     boundary edge.
C
      IF (LIST(LP) .NE. NR) GO TO 33
C
C   Set NEXT to the neighbor following NR, and test for
C     termination of the store loop.
C
    8 LP = LPTR(LP)
      NEXT = ABS(LIST(LP))
      IF (NEXT .EQ. N2) GO TO 9
C
C   Set NL or NR to NEXT.
C
      IF ( LEFT(X1,Y1,Z1,X2,Y2,Z2,X(NEXT),Y(NEXT),Z(NEXT)) )
     .    THEN
        NL = NEXT
      ELSE
        NR = NEXT
      ENDIF
      GO TO 6
C
C IWL is the number of arcs which intersect N1-N2.
C   Store LWK.
C
    9 LWK = IWL
      IWEND = IWL
C
C Initialize for edge swapping loop -- all possible swaps
C   are applied (even if the new arc again intersects
C   N1-N2), arcs to the left of N1->N2 are stored in the
C   left portion of IWK, and arcs to the right are stored in
C   the right portion.  IWF and IWL index the first and last
C   intersecting arcs.
C
      IWF = 1
C
C Top of loop -- set N0 to N1 and NL->NR to the first edge.
C   IWC points to the arc currently being processed.  LFT
C   .LE. 0 iff N0 LEFT N1->N2.
C
   10 LFT = 0
      N0 = N1
      X0 = X1
      Y0 = Y1
      Z0 = Z1
      NL = IWK(1,IWF)
      NR = IWK(2,IWF)
      IWC = IWF
C
C   Set NEXT to the node opposite NL->NR unless IWC is the
C     last arc.
C
   11 IF (IWC .EQ. IWL) GO TO 21
      IWCP1 = IWC + 1
      NEXT = IWK(1,IWCP1)
      IF (NEXT .NE. NL) GO TO 16
      NEXT = IWK(2,IWCP1)
C
C   NEXT RIGHT N1->N2 and IWC .LT. IWL.  Test for a possible
C     swap.
C
      IF ( .NOT. LEFT(X0,Y0,Z0,X(NR),Y(NR),Z(NR),X(NEXT),
     .                Y(NEXT),Z(NEXT)) ) GO TO 14
      IF (LFT .GE. 0) GO TO 12
      IF ( .NOT. LEFT(X(NL),Y(NL),Z(NL),X0,Y0,Z0,X(NEXT),
     .                Y(NEXT),Z(NEXT)) ) GO TO 14
C
C   Replace NL->NR with N0->NEXT.
C
      CALL SWAP (NEXT,N0,NL,NR, LIST,LPTR,LEND, LP21)
      IWK(1,IWC) = N0
      IWK(2,IWC) = NEXT
      GO TO 15
C
C   Swap NL-NR for N0-NEXT, shift columns IWC+1,...,IWL to
C     the left, and store N0-NEXT in the right portion of
C     IWK.
C
   12 CALL SWAP (NEXT,N0,NL,NR, LIST,LPTR,LEND, LP21)
      DO 13 I = IWCP1,IWL
        IWK(1,I-1) = IWK(1,I)
        IWK(2,I-1) = IWK(2,I)
   13   CONTINUE
      IWK(1,IWL) = N0
      IWK(2,IWL) = NEXT
      IWL = IWL - 1
      NR = NEXT
      GO TO 11
C
C   A swap is not possible.  Set N0 to NR.
C
   14 N0 = NR
      X0 = X(N0)
      Y0 = Y(N0)
      Z0 = Z(N0)
      LFT = 1
C
C   Advance to the next arc.
C
   15 NR = NEXT
      IWC = IWC + 1
      GO TO 11
C
C   NEXT LEFT N1->N2, NEXT .NE. N2, and IWC .LT. IWL.
C     Test for a possible swap.
C
   16 IF ( .NOT. LEFT(X(NL),Y(NL),Z(NL),X0,Y0,Z0,X(NEXT),
     .                Y(NEXT),Z(NEXT)) ) GO TO 19
      IF (LFT .LE. 0) GO TO 17
      IF ( .NOT. LEFT(X0,Y0,Z0,X(NR),Y(NR),Z(NR),X(NEXT),
     .                Y(NEXT),Z(NEXT)) ) GO TO 19
C
C   Replace NL->NR with NEXT->N0.
C
      CALL SWAP (NEXT,N0,NL,NR, LIST,LPTR,LEND, LP21)
      IWK(1,IWC) = NEXT
      IWK(2,IWC) = N0
      GO TO 20
C
C   Swap NL-NR for N0-NEXT, shift columns IWF,...,IWC-1 to
C     the right, and store N0-NEXT in the left portion of
C     IWK.
C
   17 CALL SWAP (NEXT,N0,NL,NR, LIST,LPTR,LEND, LP21)
      DO 18 I = IWC-1,IWF,-1
        IWK(1,I+1) = IWK(1,I)
        IWK(2,I+1) = IWK(2,I)
   18   CONTINUE
      IWK(1,IWF) = N0
      IWK(2,IWF) = NEXT
      IWF = IWF + 1
      GO TO 20
C
C   A swap is not possible.  Set N0 to NL.
C
   19 N0 = NL
      X0 = X(N0)
      Y0 = Y(N0)
      Z0 = Z(N0)
      LFT = -1
C
C   Advance to the next arc.
C
   20 NL = NEXT
      IWC = IWC + 1
      GO TO 11
C
C   N2 is opposite NL->NR (IWC = IWL).
C
   21 IF (N0 .EQ. N1) GO TO 24
      IF (LFT .LT. 0) GO TO 22
C
C   N0 RIGHT N1->N2.  Test for a possible swap.
C
      IF ( .NOT. LEFT(X0,Y0,Z0,X(NR),Y(NR),Z(NR),X2,Y2,Z2) )
     .  GO TO 10
C
C   Swap NL-NR for N0-N2 and store N0-N2 in the right
C     portion of IWK.
C
      CALL SWAP (N2,N0,NL,NR, LIST,LPTR,LEND, LP21)
      IWK(1,IWL) = N0
      IWK(2,IWL) = N2
      IWL = IWL - 1
      GO TO 10
C
C   N0 LEFT N1->N2.  Test for a possible swap.
C
   22 IF ( .NOT. LEFT(X(NL),Y(NL),Z(NL),X0,Y0,Z0,X2,Y2,Z2) )
     .  GO TO 10
C
C   Swap NL-NR for N0-N2, shift columns IWF,...,IWL-1 to the
C     right, and store N0-N2 in the left portion of IWK.
C
      CALL SWAP (N2,N0,NL,NR, LIST,LPTR,LEND, LP21)
      I = IWL
   23 IWK(1,I) = IWK(1,I-1)
      IWK(2,I) = IWK(2,I-1)
      I = I - 1
      IF (I .GT. IWF) GO TO 23
      IWK(1,IWF) = N0
      IWK(2,IWF) = N2
      IWF = IWF + 1
      GO TO 10
C
C IWF = IWC = IWL.  Swap out the last arc for N1-N2 and
C   store zeros in IWK.
C
   24 CALL SWAP (N2,N1,NL,NR, LIST,LPTR,LEND, LP21)
      IWK(1,IWC) = 0
      IWK(2,IWC) = 0
C
C Optimization procedure --
C
      IER = 0
      IF (IWC .GT. 1) THEN
C
C   Optimize the set of new arcs to the left of IN1->IN2.
C
        NIT = 4*(IWC-1)
        CALL OPTIM (X,Y,Z,IWC-1, LIST,LPTR,LEND,NIT,
     .              IWK, IERR)
        IF (IERR .NE. 0  .AND.  IERR .NE. 1) GO TO 34
        IF (IERR .EQ. 1) IER = 5
      ENDIF
      IF (IWC .LT. IWEND) THEN
C
C   Optimize the set of new arcs to the right of IN1->IN2.
C
        NIT = 4*(IWEND-IWC)
        CALL OPTIM (X,Y,Z,IWEND-IWC, LIST,LPTR,LEND,NIT,
     .              IWK(1,IWC+1), IERR)
        IF (IERR .NE. 0  .AND.  IERR .NE. 1) GO TO 34
        IF (IERR .EQ. 1) GO TO 35
      ENDIF
      IF (IER .EQ. 5) GO TO 35
C
C Successful termination (IER = 0).
C
      RETURN
C
C IN1 and IN2 were adjacent on input.
C
   30 IER = 0
      RETURN
C
C Invalid input parameter.
C
   31 IER = 1
      RETURN
C
C Insufficient space reserved for IWK.
C
   32 IER = 2
      RETURN
C
C Invalid triangulation data structure or collinear nodes
C   on convex hull boundary.
C
   33 IER = 3
      WRITE (*,130) IN1, IN2
  130 FORMAT (//5X,'*** Error in EDGE:  Invalid triangula',
     .        'tion or null triangles on boundary'/
     .        9X,'IN1 =',I4,', IN2=',I4/)
      RETURN
C
C Error flag (other than 1) returned by OPTIM.
C
   34 IER = 4
      WRITE (*,140) NIT, IERR
  140 FORMAT (//5X,'*** Error in OPTIM (called from EDGE):',
     .        '  NIT = ',I4,', IER = ',I1,' ***'/)
      RETURN
C
C Error flag 1 returned by OPTIM.
C
   35 IER = 5
      RETURN
      END
      SUBROUTINE GETNP (X,Y,Z,LIST,LPTR,LEND,L, NPTS, DF,
     .                  IER)
      INTEGER LIST(*), LPTR(*), LEND(*), L, NPTS(L), IER
      DOUBLE PRECISION    X(*), Y(*), Z(*), DF
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/28/98
C
C   Given a Delaunay triangulation of N nodes on the unit
C sphere and an array NPTS containing the indexes of L-1
C nodes ordered by angular distance from NPTS(1), this sub-
C routine sets NPTS(L) to the index of the next node in the
C sequence -- the node, other than NPTS(1),...,NPTS(L-1),
C that is closest to NPTS(1).  Thus, the ordered sequence
C of K closest nodes to N1 (including N1) may be determined
C by K-1 calls to GETNP with NPTS(1) = N1 and L = 2,3,...,K
C for K .GE. 2.
C
C   The algorithm uses the property of a Delaunay triangula-
C tion that the K-th closest node to N1 is a neighbor of one
C of the K-1 closest nodes to N1.
C
C
C On input:
C
C       X,Y,Z = Arrays of length N containing the Cartesian
C               coordinates of the nodes.
C
C       LIST,LPTR,LEND = Triangulation data structure.  Re-
C                        fer to Subroutine TRMESH.
C
C       L = Number of nodes in the sequence on output.  2
C           .LE. L .LE. N.
C
C The above parameters are not altered by this routine.
C
C       NPTS = Array of length .GE. L containing the indexes
C              of the L-1 closest nodes to NPTS(1) in the
C              first L-1 locations.
C
C On output:
C
C       NPTS = Array updated with the index of the L-th
C              closest node to NPTS(1) in position L unless
C              IER = 1.
C
C       DF = Value of an increasing function (negative cos-
C            ine) of the angular distance between NPTS(1)
C            and NPTS(L) unless IER = 1.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if L < 2.
C
C Modules required by GETNP:  None
C
C Intrinsic function called by GETNP:  ABS
C
C***********************************************************
C
      INTEGER I, LM1, LP, LPL, N1, NB, NI, NP
      DOUBLE PRECISION    DNB, DNP, X1, Y1, Z1
C
C Local parameters:
C
C DNB,DNP =  Negative cosines of the angular distances from
C              N1 to NB and to NP, respectively
C I =        NPTS index and DO-loop index
C LM1 =      L-1
C LP =       LIST pointer of a neighbor of NI
C LPL =      Pointer to the last neighbor of NI
C N1 =       NPTS(1)
C NB =       Neighbor of NI and candidate for NP
C NI =       NPTS(I)
C NP =       Candidate for NPTS(L)
C X1,Y1,Z1 = Coordinates of N1
C
      LM1 = L - 1
      IF (LM1 .LT. 1) GO TO 6
      IER = 0
C
C Store N1 = NPTS(1) and mark the elements of NPTS.
C
      N1 = NPTS(1)
      X1 = X(N1)
      Y1 = Y(N1)
      Z1 = Z(N1)
      DO 1 I = 1,LM1
        NI = NPTS(I)
        LEND(NI) = -LEND(NI)
    1   CONTINUE
C
C Candidates for NP = NPTS(L) are the unmarked neighbors
C   of nodes in NPTS.  DNP is initially greater than -cos(PI)
C   (the maximum distance).
C
      DNP = 2.
C
C Loop on nodes NI in NPTS.
C
      DO 4 I = 1,LM1
        NI = NPTS(I)
        LPL = -LEND(NI)
        LP = LPL
C
C Loop on neighbors NB of NI.
C
    2   NB = ABS(LIST(LP))
          IF (LEND(NB) .LT. 0) GO TO 3
C
C NB is an unmarked neighbor of NI.  Replace NP if NB is
C   closer to N1.
C
          DNB = -(X(NB)*X1 + Y(NB)*Y1 + Z(NB)*Z1)
          IF (DNB .GE. DNP) GO TO 3
          NP = NB
          DNP = DNB
    3     LP = LPTR(LP)
          IF (LP .NE. LPL) GO TO 2
    4   CONTINUE
      NPTS(L) = NP
      DF = DNP
C
C Unmark the elements of NPTS.
C
      DO 5 I = 1,LM1
        NI = NPTS(I)
        LEND(NI) = -LEND(NI)
    5   CONTINUE
      RETURN
C
C L is outside its valid range.
C
    6 IER = 1
      RETURN
      END
      SUBROUTINE INSERT (K,LP, LIST,LPTR,LNEW )
      INTEGER K, LP, LIST(*), LPTR(*), LNEW
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/17/96
C
C   This subroutine inserts K as a neighbor of N1 following
C N2, where LP is the LIST pointer of N2 as a neighbor of
C N1.  Note that, if N2 is the last neighbor of N1, K will
C become the first neighbor (even if N1 is a boundary node).
C
C   This routine is identical to the similarly named routine
C in TRIPACK.
C
C
C On input:
C
C       K = Index of the node to be inserted.
C
C       LP = LIST pointer of N2 as a neighbor of N1.
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR,LNEW = Data structure defining the trian-
C                        gulation.  Refer to Subroutine
C                        TRMESH.
C
C On output:
C
C       LIST,LPTR,LNEW = Data structure updated with the
C                        addition of node K.
C
C Modules required by INSERT:  None
C
C***********************************************************
C
      INTEGER LSAV
C
      LSAV = LPTR(LP)
      LPTR(LP) = LNEW
      LIST(LNEW) = K
      LPTR(LNEW) = LSAV
      LNEW = LNEW + 1
      RETURN
      END
      LOGICAL FUNCTION INSIDE (P,LV,XV,YV,ZV,NV,LISTV, IER)
      INTEGER LV, NV, LISTV(NV), IER
      DOUBLE PRECISION    P(3), XV(LV), YV(LV), ZV(LV)
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   12/27/93
C
C   This function locates a point P relative to a polygonal
C region R on the surface of the unit sphere, returning
C INSIDE = TRUE if and only if P is contained in R.  R is
C defined by a cyclically ordered sequence of vertices which
C form a positively-oriented simple closed curve.  Adjacent
C vertices need not be distinct but the curve must not be
C self-intersecting.  Also, while polygon edges are by defi-
C nition restricted to a single hemisphere, R is not so
C restricted.  Its interior is the region to the left as the
C vertices are traversed in order.
C
C   The algorithm consists of selecting a point Q in R and
C then finding all points at which the great circle defined
C by P and Q intersects the boundary of R.  P lies inside R
C if and only if there is an even number of intersection
C points between Q and P.  Q is taken to be a point immedi-
C ately to the left of a directed boundary edge -- the first
C one that results in no consistency-check failures.
C
C   If P is close to the polygon boundary, the problem is
C ill-conditioned and the decision may be incorrect.  Also,
C an incorrect decision may result from a poor choice of Q
C (if, for example, a boundary edge lies on the great cir-
C cle defined by P and Q).  A more reliable result could be
C obtained by a sequence of calls to INSIDE with the ver-
C tices cyclically permuted before each call (to alter the
C choice of Q).
C
C
C On input:
C
C       P = Array of length 3 containing the Cartesian
C           coordinates of the point (unit vector) to be
C           located.
C
C       LV = Length of arrays XV, YV, and ZV.
C
C       XV,YV,ZV = Arrays of length LV containing the Carte-
C                  sian coordinates of unit vectors (points
C                  on the unit sphere).  These values are
C                  not tested for validity.
C
C       NV = Number of vertices in the polygon.  3 .LE. NV
C            .LE. LV.
C
C       LISTV = Array of length NV containing the indexes
C               (for XV, YV, and ZV) of a cyclically-ordered
C               (and CCW-ordered) sequence of vertices that
C               define R.  The last vertex (indexed by
C               LISTV(NV)) is followed by the first (indexed
C               by LISTV(1)).  LISTV entries must be in the
C               range 1 to LV.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       INSIDE = TRUE if and only if P lies inside R unless
C                IER .NE. 0, in which case the value is not
C                altered.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if LV or NV is outside its valid
C                     range.
C             IER = 2 if a LISTV entry is outside its valid
C                     range.
C             IER = 3 if the polygon boundary was found to
C                     be self-intersecting.  This error will
C                     not necessarily be detected.
C             IER = 4 if every choice of Q (one for each
C                     boundary edge) led to failure of some
C                     internal consistency check.  The most
C                     likely cause of this error is invalid
C                     input:  P = (0,0,0), a null or self-
C                     intersecting polygon, etc.
C
C Module required by INSIDE:  INTRSC
C
C Intrinsic function called by INSIDE:  SQRT
C
C***********************************************************
C
      INTEGER I1, I2, IERR, IMX, K, K0, N, NI
      LOGICAL EVEN, LFT1, LFT2, PINR, QINR
      DOUBLE PRECISION    B(3), BP, BQ, CN(3), D, EPS, PN(3), Q(3),
     .        QN(3), QNRM, V1(3), V2(3), VN(3), VNRM
C
C Local parameters:
C
C B =         Intersection point between the boundary and
C               the great circle defined by P and Q
C BP,BQ =     <B,P> and <B,Q>, respectively, maximized over
C               intersection points B that lie between P and
C               Q (on the shorter arc) -- used to find the
C               closest intersection points to P and Q
C CN =        Q X P = normal to the plane of P and Q
C D =         Dot product <B,P> or <B,Q>
C EPS =       Parameter used to define Q as the point whose
C               orthogonal distance to (the midpoint of)
C               boundary edge V1->V2 is approximately EPS/
C               (2*Cos(A/2)), where <V1,V2> = Cos(A).
C EVEN =      TRUE iff an even number of intersection points
C               lie between P and Q (on the shorter arc)
C I1,I2 =     Indexes (LISTV elements) of a pair of adjacent
C               boundary vertices (endpoints of a boundary
C               edge)
C IERR =      Error flag for calls to INTRSC (not tested)
C IMX =       Local copy of LV and maximum value of I1 and
C               I2
C K =         DO-loop index and LISTV index
C K0 =        LISTV index of the first endpoint of the
C               boundary edge used to compute Q
C LFT1,LFT2 = Logical variables associated with I1 and I2 in
C               the boundary traversal:  TRUE iff the vertex
C               is strictly to the left of Q->P (<V,CN> > 0)
C N =         Local copy of NV
C NI =        Number of intersections (between the boundary
C               curve and the great circle P-Q) encountered
C PINR =      TRUE iff P is to the left of the directed
C               boundary edge associated with the closest
C               intersection point to P that lies between P
C               and Q (a left-to-right intersection as
C               viewed from Q), or there is no intersection
C               between P and Q (on the shorter arc)
C PN,QN =     P X CN and CN X Q, respectively:  used to
C               locate intersections B relative to arc Q->P
C Q =         (V1 + V2 + EPS*VN/VNRM)/QNRM, where V1->V2 is
C               the boundary edge indexed by LISTV(K0) ->
C               LISTV(K0+1)
C QINR =      TRUE iff Q is to the left of the directed
C               boundary edge associated with the closest
C               intersection point to Q that lies between P
C               and Q (a right-to-left intersection as
C               viewed from Q), or there is no intersection
C               between P and Q (on the shorter arc)
C QNRM =      Euclidean norm of V1+V2+EPS*VN/VNRM used to
C               compute (normalize) Q
C V1,V2 =     Vertices indexed by I1 and I2 in the boundary
C               traversal
C VN =        V1 X V2, where V1->V2 is the boundary edge
C               indexed by LISTV(K0) -> LISTV(K0+1)
C VNRM =      Euclidean norm of VN
C
      DATA EPS/1.E-3/
C
C Store local parameters, test for error 1, and initialize
C   K0.
C
      IMX = LV
      N = NV
      IF (N .LT. 3  .OR.  N .GT. IMX) GO TO 11
      K0 = 0
      I1 = LISTV(1)
      IF (I1 .LT. 1  .OR.  I1 .GT. IMX) GO TO 12
C
C Increment K0 and set Q to a point immediately to the left
C   of the midpoint of edge V1->V2 = LISTV(K0)->LISTV(K0+1):
C   Q = (V1 + V2 + EPS*VN/VNRM)/QNRM, where VN = V1 X V2.
C
    1 K0 = K0 + 1
      IF (K0 .GT. N) GO TO 14
      I1 = LISTV(K0)
      IF (K0 .LT. N) THEN
        I2 = LISTV(K0+1)
      ELSE
        I2 = LISTV(1)
      ENDIF
      IF (I2 .LT. 1  .OR.  I2 .GT. IMX) GO TO 12
      VN(1) = YV(I1)*ZV(I2) - ZV(I1)*YV(I2)
      VN(2) = ZV(I1)*XV(I2) - XV(I1)*ZV(I2)
      VN(3) = XV(I1)*YV(I2) - YV(I1)*XV(I2)
      VNRM = SQRT(VN(1)*VN(1) + VN(2)*VN(2) + VN(3)*VN(3))
      IF (VNRM .EQ. 0.) GO TO 1
      Q(1) = XV(I1) + XV(I2) + EPS*VN(1)/VNRM
      Q(2) = YV(I1) + YV(I2) + EPS*VN(2)/VNRM
      Q(3) = ZV(I1) + ZV(I2) + EPS*VN(3)/VNRM
      QNRM = SQRT(Q(1)*Q(1) + Q(2)*Q(2) + Q(3)*Q(3))
      Q(1) = Q(1)/QNRM
      Q(2) = Q(2)/QNRM
      Q(3) = Q(3)/QNRM
C
C Compute CN = Q X P, PN = P X CN, and QN = CN X Q.
C
      CN(1) = Q(2)*P(3) - Q(3)*P(2)
      CN(2) = Q(3)*P(1) - Q(1)*P(3)
      CN(3) = Q(1)*P(2) - Q(2)*P(1)
      IF (CN(1) .EQ. 0.  .AND.  CN(2) .EQ. 0.  .AND.
     .    CN(3) .EQ. 0.) GO TO 1
      PN(1) = P(2)*CN(3) - P(3)*CN(2)
      PN(2) = P(3)*CN(1) - P(1)*CN(3)
      PN(3) = P(1)*CN(2) - P(2)*CN(1)
      QN(1) = CN(2)*Q(3) - CN(3)*Q(2)
      QN(2) = CN(3)*Q(1) - CN(1)*Q(3)
      QN(3) = CN(1)*Q(2) - CN(2)*Q(1)
C
C Initialize parameters for the boundary traversal.
C
      NI = 0
      EVEN = .TRUE.
      BP = -2.
      BQ = -2.
      PINR = .TRUE.
      QINR = .TRUE.
      I2 = LISTV(N)
      IF (I2 .LT. 1  .OR.  I2 .GT. IMX) GO TO 12
      LFT2 = CN(1)*XV(I2) + CN(2)*YV(I2) +
     .       CN(3)*ZV(I2) .GT. 0.
C
C Loop on boundary arcs I1->I2.
C
      DO 2 K = 1,N
        I1 = I2
        LFT1 = LFT2
        I2 = LISTV(K)
        IF (I2 .LT. 1  .OR.  I2 .GT. IMX) GO TO 12
        LFT2 = CN(1)*XV(I2) + CN(2)*YV(I2) +
     .         CN(3)*ZV(I2) .GT. 0.
        IF (LFT1 .EQV. LFT2) GO TO 2
C
C   I1 and I2 are on opposite sides of Q->P.  Compute the
C     point of intersection B.
C
        NI = NI + 1
        V1(1) = XV(I1)
        V1(2) = YV(I1)
        V1(3) = ZV(I1)
        V2(1) = XV(I2)
        V2(2) = YV(I2)
        V2(3) = ZV(I2)
        CALL INTRSC (V1,V2,CN, B,IERR)
C
C   B is between Q and P (on the shorter arc) iff
C     B Forward Q->P and B Forward P->Q       iff
C     <B,QN> > 0 and <B,PN> > 0.
C
        IF (B(1)*QN(1) + B(2)*QN(2) + B(3)*QN(3) .GT. 0.
     .      .AND.
     .      B(1)*PN(1) + B(2)*PN(2) + B(3)*PN(3) .GT. 0.)
     .    THEN
C
C   Update EVEN, BQ, QINR, BP, and PINR.
C
          EVEN = .NOT. EVEN
          D = B(1)*Q(1) + B(2)*Q(2) + B(3)*Q(3)
          IF (D .GT. BQ) THEN
            BQ = D
            QINR = LFT2
          ENDIF
          D = B(1)*P(1) + B(2)*P(2) + B(3)*P(3)
          IF (D .GT. BP) THEN
            BP = D
            PINR = LFT1
          ENDIF
        ENDIF
    2   CONTINUE
C
C Test for consistency:  NI must be even and QINR must be
C   TRUE.
C
      IF (NI .NE. 2*(NI/2)  .OR.  .NOT. QINR) GO TO 1
C
C Test for error 3:  different values of PINR and EVEN.
C
      IF (PINR .NEQV. EVEN) GO TO 13
C
C No error encountered.
C
      IER = 0
      INSIDE = EVEN
      RETURN
C
C LV or NV is outside its valid range.
C
   11 IER = 1
      RETURN
C
C A LISTV entry is outside its valid range.
C
   12 IER = 2
      RETURN
C
C The polygon boundary is self-intersecting.
C
   13 IER = 3
      RETURN
C
C Consistency tests failed for all values of Q.
C
   14 IER = 4
      RETURN
      END
      SUBROUTINE INTADD (KK,I1,I2,I3, LIST,LPTR,LEND,LNEW )
      INTEGER KK, I1, I2, I3, LIST(*), LPTR(*), LEND(*),
     .        LNEW
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/17/96
C
C   This subroutine adds an interior node to a triangulation
C of a set of points on the unit sphere.  The data structure
C is updated with the insertion of node KK into the triangle
C whose vertices are I1, I2, and I3.  No optimization of the
C triangulation is performed.
C
C   This routine is identical to the similarly named routine
C in TRIPACK.
C
C
C On input:
C
C       KK = Index of the node to be inserted.  KK .GE. 1
C            and KK must not be equal to I1, I2, or I3.
C
C       I1,I2,I3 = Indexes of the counterclockwise-ordered
C                  sequence of vertices of a triangle which
C                  contains node KK.
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR,LEND,LNEW = Data structure defining the
C                             triangulation.  Refer to Sub-
C                             routine TRMESH.  Triangle
C                             (I1,I2,I3) must be included
C                             in the triangulation.
C
C On output:
C
C       LIST,LPTR,LEND,LNEW = Data structure updated with
C                             the addition of node KK.  KK
C                             will be connected to nodes I1,
C                             I2, and I3.
C
C Modules required by INTADD:  INSERT, LSTPTR
C
C***********************************************************
C
      INTEGER LSTPTR
      INTEGER K, LP, N1, N2, N3
C
C Local parameters:
C
C K =        Local copy of KK
C LP =       LIST pointer
C N1,N2,N3 = Local copies of I1, I2, and I3
C
      K = KK
C
C Initialization.
C
      N1 = I1
      N2 = I2
      N3 = I3
C
C Add K as a neighbor of I1, I2, and I3.
C
      LP = LSTPTR(LEND(N1),N2,LIST,LPTR)
      CALL INSERT (K,LP, LIST,LPTR,LNEW )
      LP = LSTPTR(LEND(N2),N3,LIST,LPTR)
      CALL INSERT (K,LP, LIST,LPTR,LNEW )
      LP = LSTPTR(LEND(N3),N1,LIST,LPTR)
      CALL INSERT (K,LP, LIST,LPTR,LNEW )
C
C Add I1, I2, and I3 as neighbors of K.
C
      LIST(LNEW) = N1
      LIST(LNEW+1) = N2
      LIST(LNEW+2) = N3
      LPTR(LNEW) = LNEW + 1
      LPTR(LNEW+1) = LNEW + 2
      LPTR(LNEW+2) = LNEW
      LEND(K) = LNEW + 2
      LNEW = LNEW + 3
      RETURN
      END
      SUBROUTINE INTRSC (P1,P2,CN, P,IER)
      INTEGER IER
      DOUBLE PRECISION    P1(3), P2(3), CN(3), P(3)
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/19/90
C
C   Given a great circle C and points P1 and P2 defining an
C arc A on the surface of the unit sphere, where A is the
C shorter of the two portions of the great circle C12 assoc-
C iated with P1 and P2, this subroutine returns the point
C of intersection P between C and C12 that is closer to A.
C Thus, if P1 and P2 lie in opposite hemispheres defined by
C C, P is the point of intersection of C with A.
C
C
C On input:
C
C       P1,P2 = Arrays of length 3 containing the Cartesian
C               coordinates of unit vectors.
C
C       CN = Array of length 3 containing the Cartesian
C            coordinates of a nonzero vector which defines C
C            as the intersection of the plane whose normal
C            is CN with the unit sphere.  Thus, if C is to
C            be the great circle defined by P and Q, CN
C            should be P X Q.
C
C The above parameters are not altered by this routine.
C
C       P = Array of length 3.
C
C On output:
C
C       P = Point of intersection defined above unless IER
C           .NE. 0, in which case P is not altered.
C
C       IER = Error indicator.
C             IER = 0 if no errors were encountered.
C             IER = 1 if <CN,P1> = <CN,P2>.  This occurs
C                     iff P1 = P2 or CN = 0 or there are
C                     two intersection points at the same
C                     distance from A.
C             IER = 2 if P2 = -P1 and the definition of A is
C                     therefore ambiguous.
C
C Modules required by INTRSC:  None
C
C Intrinsic function called by INTRSC:  SQRT
C
C***********************************************************
C
      INTEGER I
      DOUBLE PRECISION    D1, D2, PP(3), PPN, T
C
C Local parameters:
C
C D1 =  <CN,P1>
C D2 =  <CN,P2>
C I =   DO-loop index
C PP =  P1 + T*(P2-P1) = Parametric representation of the
C         line defined by P1 and P2
C PPN = Norm of PP
C T =   D1/(D1-D2) = Parameter value chosen so that PP lies
C         in the plane of C
C
      D1 = CN(1)*P1(1) + CN(2)*P1(2) + CN(3)*P1(3)
      D2 = CN(1)*P2(1) + CN(2)*P2(2) + CN(3)*P2(3)
C
      IF (D1 .EQ. D2) THEN
        IER = 1
        RETURN
      ENDIF
C
C Solve for T such that <PP,CN> = 0 and compute PP and PPN.
C
      T = D1/(D1-D2)
      PPN = 0.
      DO 1 I = 1,3
        PP(I) = P1(I) + T*(P2(I)-P1(I))
        PPN = PPN + PP(I)*PP(I)
    1   CONTINUE
C
C PPN = 0 iff PP = 0 iff P2 = -P1 (and T = .5).
C
      IF (PPN .EQ. 0.) THEN
        IER = 2
        RETURN
      ENDIF
      PPN = SQRT(PPN)
C
C Compute P = PP/PPN.
C
      DO 2 I = 1,3
        P(I) = PP(I)/PPN
    2   CONTINUE
      IER = 0
      RETURN
      END
      INTEGER FUNCTION JRAND (N, IX,IY,IZ )
      INTEGER N, IX, IY, IZ
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/28/98
C
C   This function returns a uniformly distributed pseudo-
C random integer in the range 1 to N.
C
C
C On input:
C
C       N = Maximum value to be returned.
C
C N is not altered by this function.
C
C       IX,IY,IZ = Integer seeds initialized to values in
C                  the range 1 to 30,000 before the first
C                  call to JRAND, and not altered between
C                  subsequent calls (unless a sequence of
C                  random numbers is to be repeated by
C                  reinitializing the seeds).
C
C On output:
C
C       IX,IY,IZ = Updated integer seeds.
C
C       JRAND = Random integer in the range 1 to N.
C
C Reference:  B. A. Wichmann and I. D. Hill, "An Efficient
C             and Portable Pseudo-random Number Generator",
C             Applied Statistics, Vol. 31, No. 2, 1982,
C             pp. 188-190.
C
C Modules required by JRAND:  None
C
C Intrinsic functions called by JRAND:  INT, MOD, REAL
C
C***********************************************************
C
      DOUBLE PRECISION U, X
C
C Local parameters:
C
C U = Pseudo-random number uniformly distributed in the
C     interval (0,1).
C X = Pseudo-random number in the range 0 to 3 whose frac-
C       tional part is U.
C
      IX = MOD(171*IX,30269)
      IY = MOD(172*IY,30307)
      IZ = MOD(170*IZ,30323)
      X = (REAL(IX)/30269.) + (REAL(IY)/30307.) +
     .    (REAL(IZ)/30323.)
      U = X - INT(X)
      JRAND = REAL(N)*U + 1.
      RETURN
      END
      LOGICAL FUNCTION LEFT (X1,Y1,Z1,X2,Y2,Z2,X0,Y0,Z0)
      DOUBLE PRECISION X1, Y1, Z1, X2, Y2, Z2, X0, Y0, Z0
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/15/96
C
C   This function determines whether node N0 is in the
C (closed) left hemisphere defined by the plane containing
C N1, N2, and the origin, where left is defined relative to
C an observer at N1 facing N2.
C
C
C On input:
C
C       X1,Y1,Z1 = Coordinates of N1.
C
C       X2,Y2,Z2 = Coordinates of N2.
C
C       X0,Y0,Z0 = Coordinates of N0.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       LEFT = TRUE if and only if N0 is in the closed
C              left hemisphere.
C
C Modules required by LEFT:  None
C
C***********************************************************
C
C LEFT = TRUE iff <N0,N1 X N2> = det(N0,N1,N2) .GE. 0.
C
      LEFT = X0*(Y1*Z2-Y2*Z1) - Y0*(X1*Z2-X2*Z1) +
     .       Z0*(X1*Y2-X2*Y1) .GE. 0.
      RETURN
      END
      INTEGER FUNCTION LSTPTR (LPL,NB,LIST,LPTR)
      INTEGER LPL, NB, LIST(*), LPTR(*)
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/15/96
C
C   This function returns the index (LIST pointer) of NB in
C the adjacency lst for N0, where LPL = LEND(N0).
C
C   This function is identical to the similarly named
C function in TRIPACK.
C
C
C On input:
C
C       LPL = LEND(N0)
C
C       NB = Index of the node whose pointer is to be re-
C            turned.  NB must be connected to N0.
C
C       LIST,LPTR = Data structure defining the triangula-
C                   tion.  Refer to Subroutine TRMESH.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       LSTPTR = Pointer such that LIST(LSTPTR) = NB or
C                LIST(LSTPTR) = -NB, unless NB is not a
C                neighbor of N0, in which case LSTPTR = LPL.
C
C Modules required by LSTPTR:  None
C
C***********************************************************
C
      INTEGER LP, ND
C
C Local parameters:
C
C LP = LIST pointer
C ND = Nodal index
C
      LP = LPTR(LPL)
    1 ND = LIST(LP)
        IF (ND .EQ. NB) GO TO 2
        LP = LPTR(LP)
        IF (LP .NE. LPL) GO TO 1
C
    2 LSTPTR = LP
      RETURN
      END
      INTEGER FUNCTION NBCNT (LPL,LPTR)
      INTEGER LPL, LPTR(*)
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/15/96
C
C   This function returns the number of neighbors of a node
C N0 in a triangulation created by Subroutine TRMESH.
C
C   This function is identical to the similarly named
C function in TRIPACK.
C
C
C On input:
C
C       LPL = LIST pointer to the last neighbor of N0 --
C             LPL = LEND(N0).
C
C       LPTR = Array of pointers associated with LIST.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       NBCNT = Number of neighbors of N0.
C
C Modules required by NBCNT:  None
C
C***********************************************************
C
      INTEGER K, LP
C
C Local parameters:
C
C K =  Counter for computing the number of neighbors
C LP = LIST pointer
C
      LP = LPL
      K = 1
C
    1 LP = LPTR(LP)
        IF (LP .EQ. LPL) GO TO 2
        K = K + 1
        GO TO 1
C
    2 NBCNT = K
      RETURN
      END
      INTEGER FUNCTION NEARND (P,IST,N,X,Y,Z,LIST,LPTR,
     .                         LEND, AL)
      INTEGER IST, N, LIST(*), LPTR(*), LEND(N)
      DOUBLE PRECISION    P(3), X(N), Y(N), Z(N), AL
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/28/98
C
C   Given a point P on the surface of the unit sphere and a
C Delaunay triangulation created by Subroutine TRMESH, this
C function returns the index of the nearest triangulation
C node to P.
C
C   The algorithm consists of implicitly adding P to the
C triangulation, finding the nearest neighbor to P, and
C implicitly deleting P from the triangulation.  Thus, it
C is based on the fact that, if P is a node in a Delaunay
C triangulation, the nearest node to P is a neighbor of P.
C
C
C On input:
C
C       P = Array of length 3 containing the Cartesian coor-
C           dinates of the point P to be located relative to
C           the triangulation.  It is assumed without a test
C           that P(1)**2 + P(2)**2 + P(3)**2 = 1.
C
C       IST = Index of a node at which TRFIND begins the
C             search.  Search time depends on the proximity
C             of this node to P.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y,Z = Arrays of length N containing the Cartesian
C               coordinates of the nodes.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to TRMESH.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       NEARND = Nodal index of the nearest node to P, or 0
C                if N < 3 or the triangulation data struc-
C                ture is invalid.
C
C       AL = Arc length (angular distance in radians) be-
C            tween P and NEARND unless NEARND = 0.
C
C       Note that the number of candidates for NEARND
C       (neighbors of P) is limited to LMAX defined in
C       the PARAMETER statement below.
C
C Modules required by NEARND:  JRAND, LSTPTR, TRFIND, STORE
C
C Intrinsic functions called by NEARND:  ABS, ACOS
C
C***********************************************************
C
      INTEGER   LSTPTR
      INTEGER   LMAX
      PARAMETER (LMAX=25)
      INTEGER   I1, I2, I3, L, LISTP(LMAX), LP, LP1, LP2,
     .          LPL, LPTRP(LMAX), N1, N2, N3, NN, NR, NST
      DOUBLE PRECISION      B1, B2, B3, DS1, DSR, DX1, DX2, DX3, DY1,
     .          DY2, DY3, DZ1, DZ2, DZ3
C
C Local parameters:
C
C B1,B2,B3 =  Unnormalized barycentric coordinates returned
C               by TRFIND
C DS1 =       (Negative cosine of the) distance from P to N1
C DSR =       (Negative cosine of the) distance from P to NR
C DX1,..DZ3 = Components of vectors used by the swap test
C I1,I2,I3 =  Nodal indexes of a triangle containing P, or
C               the rightmost (I1) and leftmost (I2) visible
C               boundary nodes as viewed from P
C L =         Length of LISTP/LPTRP and number of neighbors
C               of P
C LMAX =      Maximum value of L
C LISTP =     Indexes of the neighbors of P
C LPTRP =     Array of pointers in 1-1 correspondence with
C               LISTP elements
C LP =        LIST pointer to a neighbor of N1 and LISTP
C               pointer
C LP1,LP2 =   LISTP indexes (pointers)
C LPL =       Pointer to the last neighbor of N1
C N1 =        Index of a node visible from P
C N2 =        Index of an endpoint of an arc opposite P
C N3 =        Index of the node opposite N1->N2
C NN =        Local copy of N
C NR =        Index of a candidate for the nearest node to P
C NST =       Index of the node at which TRFIND begins the
C               search
C
C
C Store local parameters and test for N invalid.
C
      NN = N
      IF (NN .LT. 3) GO TO 6
      NST = IST
      IF (NST .LT. 1  .OR.  NST .GT. NN) NST = 1
C
C Find a triangle (I1,I2,I3) containing P, or the rightmost
C   (I1) and leftmost (I2) visible boundary nodes as viewed
C   from P.
C
      CALL TRFIND (NST,P,N,X,Y,Z,LIST,LPTR,LEND, B1,B2,B3,
     .             I1,I2,I3)
C
C Test for collinear nodes.
C
      IF (I1 .EQ. 0) GO TO 6
C
C Store the linked lst of 'neighbors' of P in LISTP and
C   LPTRP.  I1 is the first neighbor, and 0 is stored as
C   the last neighbor if P is not contained in a triangle.
C   L is the length of LISTP and LPTRP, and is limited to
C   LMAX.
C
      IF (I3 .NE. 0) THEN
        LISTP(1) = I1
        LPTRP(1) = 2
        LISTP(2) = I2
        LPTRP(2) = 3
        LISTP(3) = I3
        LPTRP(3) = 1
        L = 3
      ELSE
        N1 = I1
        L = 1
        LP1 = 2
        LISTP(L) = N1
        LPTRP(L) = LP1
C
C   Loop on the ordered sequence of visible boundary nodes
C     N1 from I1 to I2.
C
    1   LPL = LEND(N1)
          N1 = -LIST(LPL)
          L = LP1
          LP1 = L+1
          LISTP(L) = N1
          LPTRP(L) = LP1
          IF (N1 .NE. I2  .AND.  LP1 .LT. LMAX) GO TO 1
        L = LP1
        LISTP(L) = 0
        LPTRP(L) = 1
      ENDIF
C
C Initialize variables for a loop on arcs N1-N2 opposite P
C   in which new 'neighbors' are 'swapped' in.  N1 follows
C   N2 as a neighbor of P, and LP1 and LP2 are the LISTP
C   indexes of N1 and N2.
C
      LP2 = 1
      N2 = I1
      LP1 = LPTRP(1)
      N1 = LISTP(LP1)
C
C Begin loop:  find the node N3 opposite N1->N2.
C
    2 LP = LSTPTR(LEND(N1),N2,LIST,LPTR)
        IF (LIST(LP) .LT. 0) GO TO 3
        LP = LPTR(LP)
        N3 = ABS(LIST(LP))
C
C Swap test:  Exit the loop if L = LMAX.
C
        IF (L .EQ. LMAX) GO TO 4
        DX1 = X(N1) - P(1)
        DY1 = Y(N1) - P(2)
        DZ1 = Z(N1) - P(3)
C
        DX2 = X(N2) - P(1)
        DY2 = Y(N2) - P(2)
        DZ2 = Z(N2) - P(3)
C
        DX3 = X(N3) - P(1)
        DY3 = Y(N3) - P(2)
        DZ3 = Z(N3) - P(3)
        IF ( DX3*(DY2*DZ1 - DY1*DZ2) -
     .       DY3*(DX2*DZ1 - DX1*DZ2) +
     .       DZ3*(DX2*DY1 - DX1*DY2) .LE. 0. ) GO TO 3
C
C Swap:  Insert N3 following N2 in the adjacency lst for P.
C        The two new arcs opposite P must be tested.
C
        L = L+1
        LPTRP(LP2) = L
        LISTP(L) = N3
        LPTRP(L) = LP1
        LP1 = L
        N1 = N3
        GO TO 2
C
C No swap:  Advance to the next arc and test for termination
C           on N1 = I1 (LP1 = 1) or N1 followed by 0.
C
    3   IF (LP1 .EQ. 1) GO TO 4
        LP2 = LP1
        N2 = N1
        LP1 = LPTRP(LP1)
        N1 = LISTP(LP1)
        IF (N1 .EQ. 0) GO TO 4
        GO TO 2
C
C Set NR and DSR to the index of the nearest node to P and
C   an increasing function (negative cosine) of its distance
C   from P, respectively.
C
    4 NR = I1
      DSR = -(X(NR)*P(1) + Y(NR)*P(2) + Z(NR)*P(3))
      DO 5 LP = 2,L
        N1 = LISTP(LP)
        IF (N1 .EQ. 0) GO TO 5
        DS1 = -(X(N1)*P(1) + Y(N1)*P(2) + Z(N1)*P(3))
        IF (DS1 .LT. DSR) THEN
          NR = N1
          DSR = DS1
        ENDIF
    5   CONTINUE
      DSR = -DSR
      IF (DSR .GT. 1.0) DSR = 1.0
      AL = ACOS(DSR)
      NEARND = NR
      RETURN
C
C Invalid input.
C
    6 NEARND = 0
      RETURN
      END
      SUBROUTINE OPTIM (X,Y,Z,NA, LIST,LPTR,LEND,NIT,
     .                  IWK, IER)
      INTEGER NA, LIST(*), LPTR(*), LEND(*), NIT, IWK(2,NA),
     .        IER
      DOUBLE PRECISION    X(*), Y(*), Z(*)
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/30/98
C
C   Given a set of NA triangulation arcs, this subroutine
C optimizes the portion of the triangulation consisting of
C the quadrilaterals (pairs of adjacent triangles) which
C have the arcs as diagonals by applying the circumcircle
C test and appropriate swaps to the arcs.
C
C   An iteration consists of applying the swap test and
C swaps to all NA arcs in the order in which they are
C stored.  The iteration is repeated until no swap occurs
C or NIT iterations have been performed.  The bound on the
C number of iterations may be necessary to prevent an
C infinite loop caused by cycling (reversing the effect of a
C previous swap) due to floating point inaccuracy when four
C or more nodes are nearly cocircular.
C
C
C On input:
C
C       X,Y,Z = Arrays containing the nodal coordinates.
C
C       NA = Number of arcs in the set.  NA .GE. 0.
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to Subroutine
C                        TRMESH.
C
C       NIT = Maximum number of iterations to be performed.
C             NIT = 4*NA should be sufficient.  NIT .GE. 1.
C
C       IWK = Integer array dimensioned 2 by NA containing
C             the nodal indexes of the arc endpoints (pairs
C             of endpoints are stored in columns).
C
C On output:
C
C       LIST,LPTR,LEND = Updated triangulation data struc-
C                        ture reflecting the swaps.
C
C       NIT = Number of iterations performed.
C
C       IWK = Endpoint indexes of the new set of arcs
C             reflecting the swaps.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if a swap occurred on the last of
C                     MAXIT iterations, where MAXIT is the
C                     value of NIT on input.  The new set
C                     of arcs is not necessarily optimal
C                     in this case.
C             IER = 2 if NA < 0 or NIT < 1 on input.
C             IER = 3 if IWK(2,I) is not a neighbor of
C                     IWK(1,I) for some I in the range 1
C                     to NA.  A swap may have occurred in
C                     this case.
C             IER = 4 if a zero pointer was returned by
C                     Subroutine SWAP.
C
C Modules required by OPTIM:  LSTPTR, SWAP, SWPTST
C
C Intrinsic function called by OPTIM:  ABS
C
C***********************************************************
C
      INTEGER I, IO1, IO2, ITER, LP, LP21, LPL, LPP, MAXIT,
     .        N1, N2, NNA
      LOGICAL SWPTST
      LOGICAL SWP
C
C Local parameters:
C
C I =       Column index for IWK
C IO1,IO2 = Nodal indexes of the endpoints of an arc in IWK
C ITER =    Iteration count
C LP =      LIST pointer
C LP21 =    Parameter returned by SWAP (not used)
C LPL =     Pointer to the last neighbor of IO1
C LPP =     Pointer to the node preceding IO2 as a neighbor
C             of IO1
C MAXIT =   Input value of NIT
C N1,N2 =   Nodes opposite IO1->IO2 and IO2->IO1,
C             respectively
C NNA =     Local copy of NA
C SWP =     Flag set to TRUE iff a swap occurs in the
C             optimization loop
C
      NNA = NA
      MAXIT = NIT
      IF (NNA .LT. 0  .OR.  MAXIT .LT. 1) GO TO 7
C
C Initialize iteration count ITER and test for NA = 0.
C
      ITER = 0
      IF (NNA .EQ. 0) GO TO 5
C
C Top of loop --
C   SWP = TRUE iff a swap occurred in the current iteration.
C
    1 IF (ITER .EQ. MAXIT) GO TO 6
      ITER = ITER + 1
      SWP = .FALSE.
C
C   Inner loop on arcs IO1-IO2 --
C
      DO 4 I = 1,NNA
        IO1 = IWK(1,I)
        IO2 = IWK(2,I)
C
C   Set N1 and N2 to the nodes opposite IO1->IO2 and
C     IO2->IO1, respectively.  Determine the following:
C
C     LPL = pointer to the last neighbor of IO1,
C     LP = pointer to IO2 as a neighbor of IO1, and
C     LPP = pointer to the node N2 preceding IO2.
C
        LPL = LEND(IO1)
        LPP = LPL
        LP = LPTR(LPP)
    2   IF (LIST(LP) .EQ. IO2) GO TO 3
          LPP = LP
          LP = LPTR(LPP)
          IF (LP .NE. LPL) GO TO 2
C
C   IO2 should be the last neighbor of IO1.  Test for no
C     arc and bypass the swap test if IO1 is a boundary
C     node.
C
        IF (ABS(LIST(LP)) .NE. IO2) GO TO 8
        IF (LIST(LP) .LT. 0) GO TO 4
C
C   Store N1 and N2, or bypass the swap test if IO1 is a
C     boundary node and IO2 is its first neighbor.
C
    3   N2 = LIST(LPP)
        IF (N2 .LT. 0) GO TO 4
        LP = LPTR(LP)
        N1 = ABS(LIST(LP))
C
C   Test IO1-IO2 for a swap, and update IWK if necessary.
C
        IF ( .NOT. SWPTST(N1,N2,IO1,IO2,X,Y,Z) ) GO TO 4
        CALL SWAP (N1,N2,IO1,IO2, LIST,LPTR,LEND, LP21)
        IF (LP21 .EQ. 0) GO TO 9
        SWP = .TRUE.
        IWK(1,I) = N1
        IWK(2,I) = N2
    4   CONTINUE
      IF (SWP) GO TO 1
C
C Successful termination.
C
    5 NIT = ITER
      IER = 0
      RETURN
C
C MAXIT iterations performed without convergence.
C
    6 NIT = MAXIT
      IER = 1
      RETURN
C
C Invalid input parameter.
C
    7 NIT = 0
      IER = 2
      RETURN
C
C IO2 is not a neighbor of IO1.
C
    8 NIT = ITER
      IER = 3
      RETURN
C
C Zero pointer returned by SWAP.
C
    9 NIT = ITER
      IER = 4
      RETURN
      END
      SUBROUTINE SCOORD (PX,PY,PZ, PLAT,PLON,PNRM)
      DOUBLE PRECISION PX, PY, PZ, PLAT, PLON, PNRM
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   08/27/90
C
C   This subroutine converts a point P from Cartesian coor-
C dinates to spherical coordinates.
C
C
C On input:
C
C       PX,PY,PZ = Cartesian coordinates of P.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       PLAT = Latitude of P in the range -PI/2 to PI/2, or
C              0 if PNRM = 0.  PLAT should be scaled by
C              180/PI to obtain the value in degrees.
C
C       PLON = Longitude of P in the range -PI to PI, or 0
C              if P lies on the Z-axis.  PLON should be
C              scaled by 180/PI to obtain the value in
C              degrees.
C
C       PNRM = Magnitude (Euclidean norm) of P.
C
C Modules required by SCOORD:  None
C
C Intrinsic functions called by SCOORD:  ASIN, ATAN2, SQRT
C
C***********************************************************
C
      PNRM = SQRT(PX*PX + PY*PY + PZ*PZ)
      IF (PX .NE. 0.  .OR.  PY .NE. 0.) THEN
        PLON = ATAN2(PY,PX)
      ELSE
        PLON = 0.
      ENDIF
      IF (PNRM .NE. 0.) THEN
        PLAT = ASIN(PZ/PNRM)
      ELSE
        PLAT = 0.
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION STORE (X)
      DOUBLE PRECISION X
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   05/09/92
C
C   This function forces its argument X to be stored in a
C memory location, thus providing a means of determining
C floating point number characteristics (such as the machine
C precision) when it is necessary to avoid computation in
C high precision registers.
C
C
C On input:
C
C       X = Value to be stored.
C
C X is not altered by this function.
C
C On output:
C
C       STORE = Value of X after it has been stored and
C               possibly truncated or rounded to the single
C               precision word length.
C
C Modules required by STORE:  None
C
C***********************************************************
C
      DOUBLE PRECISION Y
      COMMON/STCOM/Y
      Y = X
      STORE = Y
      RETURN
      END
      SUBROUTINE SWAP (IN1,IN2,IO1,IO2, LIST,LPTR,
     .                 LEND, LP21)
      INTEGER IN1, IN2, IO1, IO2, LIST(*), LPTR(*), LEND(*),
     .        LP21
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/22/98
C
C   Given a triangulation of a set of points on the unit
C sphere, this subroutine replaces a diagonal arc in a
C strictly convex quadrilateral (defined by a pair of adja-
C cent triangles) with the other diagonal.  Equivalently, a
C pair of adjacent triangles is replaced by another pair
C having the same union.
C
C
C On input:
C
C       IN1,IN2,IO1,IO2 = Nodal indexes of the vertices of
C                         the quadrilateral.  IO1-IO2 is re-
C                         placed by IN1-IN2.  (IO1,IO2,IN1)
C                         and (IO2,IO1,IN2) must be trian-
C                         gles on input.
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to Subroutine
C                        TRMESH.
C
C On output:
C
C       LIST,LPTR,LEND = Data structure updated with the
C                        swap -- triangles (IO1,IO2,IN1) and
C                        (IO2,IO1,IN2) are replaced by
C                        (IN1,IN2,IO2) and (IN2,IN1,IO1)
C                        unless LP21 = 0.
C
C       LP21 = Index of IN1 as a neighbor of IN2 after the
C              swap is performed unless IN1 and IN2 are
C              adjacent on input, in which case LP21 = 0.
C
C Module required by SWAP:  LSTPTR
C
C Intrinsic function called by SWAP:  ABS
C
C***********************************************************
C
      INTEGER LSTPTR
      INTEGER LP, LPH, LPSAV
C
C Local parameters:
C
C LP,LPH,LPSAV = LIST pointers
C
C
C Test for IN1 and IN2 adjacent.
C
      LP = LSTPTR(LEND(IN1),IN2,LIST,LPTR)
      IF (ABS(LIST(LP)) .EQ. IN2) THEN
        LP21 = 0
        RETURN
      ENDIF
C
C Delete IO2 as a neighbor of IO1.
C
      LP = LSTPTR(LEND(IO1),IN2,LIST,LPTR)
      LPH = LPTR(LP)
      LPTR(LP) = LPTR(LPH)
C
C If IO2 is the last neighbor of IO1, make IN2 the
C   last neighbor.
C
      IF (LEND(IO1) .EQ. LPH) LEND(IO1) = LP
C
C Insert IN2 as a neighbor of IN1 following IO1
C   using the hole created above.
C
      LP = LSTPTR(LEND(IN1),IO1,LIST,LPTR)
      LPSAV = LPTR(LP)
      LPTR(LP) = LPH
      LIST(LPH) = IN2
      LPTR(LPH) = LPSAV
C
C Delete IO1 as a neighbor of IO2.
C
      LP = LSTPTR(LEND(IO2),IN1,LIST,LPTR)
      LPH = LPTR(LP)
      LPTR(LP) = LPTR(LPH)
C
C If IO1 is the last neighbor of IO2, make IN1 the
C   last neighbor.
C
      IF (LEND(IO2) .EQ. LPH) LEND(IO2) = LP
C
C Insert IN1 as a neighbor of IN2 following IO2.
C
      LP = LSTPTR(LEND(IN2),IO2,LIST,LPTR)
      LPSAV = LPTR(LP)
      LPTR(LP) = LPH
      LIST(LPH) = IN1
      LPTR(LPH) = LPSAV
      LP21 = LPH
      RETURN
      END
      LOGICAL FUNCTION SWPTST (N1,N2,N3,N4,X,Y,Z)
      INTEGER N1, N2, N3, N4
      DOUBLE PRECISION    X(*), Y(*), Z(*)
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   03/29/91
C
C   This function decides whether or not to replace a
C diagonal arc in a quadrilateral with the other diagonal.
C The decision will be to swap (SWPTST = TRUE) if and only
C if N4 lies above the plane (in the half-space not contain-
C ing the origin) defined by (N1,N2,N3), or equivalently, if
C the projection of N4 onto this plane is interior to the
C circumcircle of (N1,N2,N3).  The decision will be for no
C swap if the quadrilateral is not strictly convex.
C
C
C On input:
C
C       N1,N2,N3,N4 = Indexes of the four nodes defining the
C                     quadrilateral with N1 adjacent to N2,
C                     and (N1,N2,N3) in counterclockwise
C                     order.  The arc connecting N1 to N2
C                     should be replaced by an arc connec-
C                     ting N3 to N4 if SWPTST = TRUE.  Refer
C                     to Subroutine SWAP.
C
C       X,Y,Z = Arrays of length N containing the Cartesian
C               coordinates of the nodes.  (X(I),Y(I),Z(I))
C               define node I for I = N1, N2, N3, and N4.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       SWPTST = TRUE if and only if the arc connecting N1
C                and N2 should be swapped for an arc con-
C                necting N3 and N4.
C
C Modules required by SWPTST:  None
C
C***********************************************************
C
      DOUBLE PRECISION DX1, DX2, DX3, DY1, DY2, DY3, DZ1, DZ2, DZ3,
     .     X4, Y4, Z4
C
C Local parameters:
C
C DX1,DY1,DZ1 = Coordinates of N4->N1
C DX2,DY2,DZ2 = Coordinates of N4->N2
C DX3,DY3,DZ3 = Coordinates of N4->N3
C X4,Y4,Z4 =    Coordinates of N4
C
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
C
C N4 lies above the plane of (N1,N2,N3) iff N3 lies above
C   the plane of (N2,N1,N4) iff Det(N3-N4,N2-N4,N1-N4) =
C   (N3-N4,N2-N4 X N1-N4) > 0.
C
      SWPTST = DX3*(DY2*DZ1 - DY1*DZ2)
     .        -DY3*(DX2*DZ1 - DX1*DZ2)
     .        +DZ3*(DX2*DY1 - DX1*DY2) .GT. 0.
      RETURN
      END
      SUBROUTINE TRANS (N,RLAT,RLON, X,Y,Z)
      INTEGER N
      DOUBLE PRECISION    RLAT(N), RLON(N), X(N), Y(N), Z(N)
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   04/08/90
C
C   This subroutine transforms spherical coordinates into
C Cartesian coordinates on the unit sphere for input to
C Subroutine TRMESH.  Storage for X and Y may coincide with
C storage for RLAT and RLON if the latter need not be saved.
C
C
C On input:
C
C       N = Number of nodes (points on the unit sphere)
C           whose coordinates are to be transformed.
C
C       RLAT = Array of length N containing latitudinal
C              coordinates of the nodes in radians.
C
C       RLON = Array of length N containing longitudinal
C              coordinates of the nodes in radians.
C
C The above parameters are not altered by this routine.
C
C       X,Y,Z = Arrays of length at least N.
C
C On output:
C
C       X,Y,Z = Cartesian coordinates in the range -1 to 1.
C               X(I)**2 + Y(I)**2 + Z(I)**2 = 1 for I = 1
C               to N.
C
C Modules required by TRANS:  None
C
C Intrinsic functions called by TRANS:  COS, SIN
C
C***********************************************************
C
      INTEGER I, NN
      DOUBLE PRECISION    COSPHI, PHI, THETA
C
C Local parameters:
C
C COSPHI = cos(PHI)
C I =      DO-loop index
C NN =     Local copy of N
C PHI =    Latitude
C THETA =  Longitude
C
      NN = N
      DO 1 I = 1,NN
        PHI = RLAT(I)
        THETA = RLON(I)
        COSPHI = COS(PHI)
        X(I) = COSPHI*COS(THETA)
        Y(I) = COSPHI*SIN(THETA)
        Z(I) = SIN(PHI)
    1   CONTINUE
      RETURN
      END
      SUBROUTINE TRFIND (NST,P,N,X,Y,Z,LIST,LPTR,LEND, B1,
     .                   B2,B3,I1,I2,I3)
      INTEGER NST, N, LIST(*), LPTR(*), LEND(N), I1, I2, I3
      DOUBLE PRECISION    P(3), X(N), Y(N), Z(N), B1, B2, B3
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/30/99
C
C   This subroutine locates a point P relative to a triangu-
C lation created by Subroutine TRMESH.  If P is contained in
C a triangle, the three vertex indexes and barycentric coor-
C dinates are returned.  Otherwise, the indexes of the
C visible boundary nodes are returned.
C
C
C On input:
C
C       NST = Index of a node at which TRFIND begins its
C             search.  Search time depends on the proximity
C             of this node to P.
C
C       P = Array of length 3 containing the x, y, and z
C           coordinates (in that order) of the point P to be
C           located.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y,Z = Arrays of length N containing the Cartesian
C               coordinates of the triangulation nodes (unit
C               vectors).  (X(I),Y(I),Z(I)) defines node I
C               for I = 1 to N.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to Subroutine
C                        TRMESH.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       B1,B2,B3 = Unnormalized barycentric coordinates of
C                  the central projection of P onto the un-
C                  derlying planar triangle if P is in the
C                  convex hull of the nodes.  These parame-
C                  ters are not altered if I1 = 0.
C
C       I1,I2,I3 = Counterclockwise-ordered vertex indexes
C                  of a triangle containing P if P is con-
C                  tained in a triangle.  If P is not in the
C                  convex hull of the nodes, I1 and I2 are
C                  the rightmost and leftmost (boundary)
C                  nodes that are visible from P, and
C                  I3 = 0.  (If all boundary nodes are vis-
C                  ible from P, then I1 and I2 coincide.)
C                  I1 = I2 = I3 = 0 if P and all of the
C                  nodes are coplanar (lie on a common great
C                  circle.
C
C Modules required by TRFIND:  JRAND, LSTPTR, STORE
C
C Intrinsic function called by TRFIND:  ABS
C
C***********************************************************
C
      INTEGER JRAND, LSTPTR
      INTEGER IX, IY, IZ, LP, N0, N1, N1S, N2, N2S, N3, N4,
     .        NEXT, NF, NL
      DOUBLE PRECISION    STORE
      DOUBLE PRECISION    DET, EPS, PTN1, PTN2, Q(3), S12, TOL, XP, YP,
     .        ZP
      DOUBLE PRECISION    X0, X1, X2, Y0, Y1, Y2, Z0, Z1, Z2
C
      SAVE    IX, IY, IZ
      DATA    IX/1/, IY/2/, IZ/3/
C
C Local parameters:
C
C EPS =      Machine precision
C IX,IY,IZ = Integer seeds for JRAND
C LP =       LIST pointer
C N0,N1,N2 = Nodes in counterclockwise order defining a
C              cone (with vertex N0) containing P, or end-
C              points of a boundary edge such that P Right
C              N1->N2
C N1S,N2S =  Initially-determined values of N1 and N2
C N3,N4 =    Nodes opposite N1->N2 and N2->N1, respectively
C NEXT =     Candidate for I1 or I2 when P is exterior
C NF,NL =    First and last neighbors of N0, or first
C              (rightmost) and last (leftmost) nodes
C              visible from P when P is exterior to the
C              triangulation
C PTN1 =     Scalar product <P,N1>
C PTN2 =     Scalar product <P,N2>
C Q =        (N2 X N1) X N2  or  N1 X (N2 X N1) -- used in
C              the boundary traversal when P is exterior
C S12 =      Scalar product <N1,N2>
C TOL =      Tolerance (multiple of EPS) defining an upper
C              bound on the magnitude of a negative bary-
C              centric coordinate (B1 or B2) for P in a
C              triangle -- used to avoid an infinite number
C              of restarts with 0 <= B3 < EPS and B1 < 0 or
C              B2 < 0 but small in magnitude
C XP,YP,ZP = Local variables containing P(1), P(2), and P(3)
C X0,Y0,Z0 = Dummy arguments for DET
C X1,Y1,Z1 = Dummy arguments for DET
C X2,Y2,Z2 = Dummy arguments for DET
C
C Statement function:
C
C DET(X1,...,Z0) .GE. 0 if and only if (X0,Y0,Z0) is in the
C                       (closed) left hemisphere defined by
C                       the plane containing (0,0,0),
C                       (X1,Y1,Z1), and (X2,Y2,Z2), where
C                       left is defined relative to an ob-
C                       server at (X1,Y1,Z1) facing
C                       (X2,Y2,Z2).
C
      DET (X1,Y1,Z1,X2,Y2,Z2,X0,Y0,Z0) = X0*(Y1*Z2-Y2*Z1)
     .     - Y0*(X1*Z2-X2*Z1) + Z0*(X1*Y2-X2*Y1)
C
C Initialize variables.
C
      XP = P(1)
      YP = P(2)
      ZP = P(3)
      N0 = NST
      IF (N0 .LT. 1  .OR.  N0 .GT. N)
     .  N0 = JRAND(N, IX,IY,IZ )
C
C Compute the relative machine precision EPS and TOL.
C
      EPS = 1.E0
    1 EPS = EPS/2.E0
        IF (STORE(EPS+1.E0) .GT. 1.E0) GO TO 1
      EPS = 2.E0*EPS
      TOL = 100.E0*EPS
C
C Set NF and NL to the first and last neighbors of N0, and
C   initialize N1 = NF.
C
    2 LP = LEND(N0)
      NL = LIST(LP)
      LP = LPTR(LP)
      NF = LIST(LP)
      N1 = NF
C
C Find a pair of adjacent neighbors N1,N2 of N0 that define
C   a wedge containing P:  P LEFT N0->N1 and P RIGHT N0->N2.
C
      IF (NL .GT. 0) THEN
C
C   N0 is an interior node.  Find N1.
C
    3   IF ( DET(X(N0),Y(N0),Z(N0),X(N1),Y(N1),Z(N1),
     .           XP,YP,ZP) .LT. 0. ) THEN
          LP = LPTR(LP)
          N1 = LIST(LP)
          IF (N1 .EQ. NL) GO TO 6
          GO TO 3
        ENDIF
      ELSE
C
C   N0 is a boundary node.  Test for P exterior.
C
        NL = -NL
        IF ( DET(X(N0),Y(N0),Z(N0),X(NF),Y(NF),Z(NF),
     .           XP,YP,ZP) .LT. 0. ) THEN
C
C   P is to the right of the boundary edge N0->NF.
C
          N1 = N0
          N2 = NF
          GO TO 9
        ENDIF
        IF ( DET(X(NL),Y(NL),Z(NL),X(N0),Y(N0),Z(N0),
     .           XP,YP,ZP) .LT. 0. ) THEN
C
C   P is to the right of the boundary edge NL->N0.
C
          N1 = NL
          N2 = N0
          GO TO 9
        ENDIF
      ENDIF
C
C P is to the left of arcs N0->N1 and NL->N0.  Set N2 to the
C   next neighbor of N0 (following N1).
C
    4 LP = LPTR(LP)
        N2 = ABS(LIST(LP))
        IF ( DET(X(N0),Y(N0),Z(N0),X(N2),Y(N2),Z(N2),
     .           XP,YP,ZP) .LT. 0. ) GO TO 7
        N1 = N2
        IF (N1 .NE. NL) GO TO 4
      IF ( DET(X(N0),Y(N0),Z(N0),X(NF),Y(NF),Z(NF),
     .         XP,YP,ZP) .LT. 0. ) GO TO 6
C
C P is left of or on arcs N0->NB for all neighbors NB
C   of N0.  Test for P = +/-N0.
C
      IF (STORE(ABS(X(N0)*XP + Y(N0)*YP + Z(N0)*ZP))
     .   .LT. 1.0-4.0*EPS) THEN
C
C   All points are collinear iff P Left NB->N0 for all
C     neighbors NB of N0.  Search the neighbors of N0.
C     Note:  N1 = NL and LP points to NL.
C
    5   IF ( DET(X(N1),Y(N1),Z(N1),X(N0),Y(N0),Z(N0),
     .           XP,YP,ZP) .GE. 0. ) THEN
          LP = LPTR(LP)
          N1 = ABS(LIST(LP))
          IF (N1 .EQ. NL) GO TO 14
          GO TO 5
        ENDIF
      ENDIF
C
C P is to the right of N1->N0, or P = +/-N0.  Set N0 to N1
C   and start over.
C
      N0 = N1
      GO TO 2
C
C P is between arcs N0->N1 and N0->NF.
C
    6 N2 = NF
C
C P is contained in a wedge defined by geodesics N0-N1 and
C   N0-N2, where N1 is adjacent to N2.  Save N1 and N2 to
C   test for cycling.
C
    7 N3 = N0
      N1S = N1
      N2S = N2
C
C Top of edge-hopping loop:
C
    8 B3 = DET(X(N1),Y(N1),Z(N1),X(N2),Y(N2),Z(N2),XP,YP,ZP)
      IF (B3 .LT. 0.) THEN
C
C   Set N4 to the first neighbor of N2 following N1 (the
C     node opposite N2->N1) unless N1->N2 is a boundary arc.
C
        LP = LSTPTR(LEND(N2),N1,LIST,LPTR)
        IF (LIST(LP) .LT. 0) GO TO 9
        LP = LPTR(LP)
        N4 = ABS(LIST(LP))
C
C   Define a new arc N1->N2 which intersects the geodesic
C     N0-P.
C
        IF ( DET(X(N0),Y(N0),Z(N0),X(N4),Y(N4),Z(N4),
     .           XP,YP,ZP) .LT. 0. ) THEN
          N3 = N2
          N2 = N4
          N1S = N1
          IF (N2 .NE. N2S  .AND.  N2 .NE. N0) GO TO 8
        ELSE
          N3 = N1
          N1 = N4
          N2S = N2
          IF (N1 .NE. N1S  .AND.  N1 .NE. N0) GO TO 8
        ENDIF
C
C   The starting node N0 or edge N1-N2 was encountered
C     again, implying a cycle (infinite loop).  Restart
C     with N0 randomly selected.
C
        N0 = JRAND(N, IX,IY,IZ )
        GO TO 2
      ENDIF
C
C P is in (N1,N2,N3) unless N0, N1, N2, and P are collinear
C   or P is close to -N0.
C
      IF (B3 .GE. EPS) THEN
C
C   B3 .NE. 0.
C
        B1 = DET(X(N2),Y(N2),Z(N2),X(N3),Y(N3),Z(N3),
     .           XP,YP,ZP)
        B2 = DET(X(N3),Y(N3),Z(N3),X(N1),Y(N1),Z(N1),
     .           XP,YP,ZP)
        IF (B1 .LT. -TOL  .OR.  B2 .LT. -TOL) THEN
C
C   Restart with N0 randomly selected.
C
          N0 = JRAND(N, IX,IY,IZ )
          GO TO 2
        ENDIF
      ELSE
C
C   B3 = 0 and thus P lies on N1->N2. Compute
C     B1 = Det(P,N2 X N1,N2) and B2 = Det(P,N1,N2 X N1).
C
        B3 = 0.
        S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
        PTN1 = XP*X(N1) + YP*Y(N1) + ZP*Z(N1)
        PTN2 = XP*X(N2) + YP*Y(N2) + ZP*Z(N2)
        B1 = PTN1 - S12*PTN2
        B2 = PTN2 - S12*PTN1
        IF (B1 .LT. -TOL  .OR.  B2 .LT. -TOL) THEN
C
C   Restart with N0 randomly selected.
C
          N0 = JRAND(N, IX,IY,IZ )
          GO TO 2
        ENDIF
      ENDIF
C
C P is in (N1,N2,N3).
C
      I1 = N1
      I2 = N2
      I3 = N3
      IF (B1 .LT. 0.0) B1 = 0.0
      IF (B2 .LT. 0.0) B2 = 0.0
      RETURN
C
C P Right N1->N2, where N1->N2 is a boundary edge.
C   Save N1 and N2, and set NL = 0 to indicate that
C   NL has not yet been found.
C
    9 N1S = N1
      N2S = N2
      NL = 0
C
C           Counterclockwise Boundary Traversal:
C
   10 LP = LEND(N2)
      LP = LPTR(LP)
      NEXT = LIST(LP)
      IF ( DET(X(N2),Y(N2),Z(N2),X(NEXT),Y(NEXT),Z(NEXT),
     .         XP,YP,ZP) .GE. 0. ) THEN
C
C   N2 is the rightmost visible node if P Forward N2->N1
C     or NEXT Forward N2->N1.  Set Q to (N2 X N1) X N2.
C
        S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
        Q(1) = X(N1) - S12*X(N2)
        Q(2) = Y(N1) - S12*Y(N2)
        Q(3) = Z(N1) - S12*Z(N2)
        IF (XP*Q(1) + YP*Q(2) + ZP*Q(3) .GE. 0.) GO TO 11
        IF (X(NEXT)*Q(1) + Y(NEXT)*Q(2) + Z(NEXT)*Q(3)
     .      .GE. 0.) GO TO 11
C
C   N1, N2, NEXT, and P are nearly collinear, and N2 is
C     the leftmost visible node.
C
        NL = N2
      ENDIF
C
C Bottom of counterclockwise loop:
C
      N1 = N2
      N2 = NEXT
      IF (N2 .NE. N1S) GO TO 10
C
C All boundary nodes are visible from P.
C
      I1 = N1S
      I2 = N1S
      I3 = 0
      RETURN
C
C N2 is the rightmost visible node.
C
   11 NF = N2
      IF (NL .EQ. 0) THEN
C
C Restore initial values of N1 and N2, and begin the search
C   for the leftmost visible node.
C
        N2 = N2S
        N1 = N1S
C
C           Clockwise Boundary Traversal:
C
   12   LP = LEND(N1)
        NEXT = -LIST(LP)
        IF ( DET(X(NEXT),Y(NEXT),Z(NEXT),X(N1),Y(N1),Z(N1),
     .           XP,YP,ZP) .GE. 0. ) THEN
C
C   N1 is the leftmost visible node if P or NEXT is
C     forward of N1->N2.  Compute Q = N1 X (N2 X N1).
C
          S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
          Q(1) = X(N2) - S12*X(N1)
          Q(2) = Y(N2) - S12*Y(N1)
          Q(3) = Z(N2) - S12*Z(N1)
          IF (XP*Q(1) + YP*Q(2) + ZP*Q(3) .GE. 0.) GO TO 13
          IF (X(NEXT)*Q(1) + Y(NEXT)*Q(2) + Z(NEXT)*Q(3)
     .        .GE. 0.) GO TO 13
C
C   P, NEXT, N1, and N2 are nearly collinear and N1 is the
C     rightmost visible node.
C
          NF = N1
        ENDIF
C
C Bottom of clockwise loop:
C
        N2 = N1
        N1 = NEXT
        IF (N1 .NE. N1S) GO TO 12
C
C All boundary nodes are visible from P.
C
        I1 = N1
        I2 = N1
        I3 = 0
        RETURN
C
C N1 is the leftmost visible node.
C
   13   NL = N1
      ENDIF
C
C NF and NL have been found.
C
      I1 = NF
      I2 = NL
      I3 = 0
      RETURN
C
C All points are collinear (coplanar).
C
   14 I1 = 0
      I2 = 0
      I3 = 0
      RETURN
      END
      SUBROUTINE TRLIST (N,LIST,LPTR,LEND,NROW, NT,LTRI,IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), NROW, NT,
     .        LTRI(NROW,*), IER
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/20/96
C
C   This subroutine converts a triangulation data structure
C from the linked lst created by Subroutine TRMESH to a
C triangle lst.
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       LIST,LPTR,LEND = Linked lst data structure defin-
C                        ing the triangulation.  Refer to
C                        Subroutine TRMESH.
C
C       NROW = Number of rows (entries per triangle) re-
C              served for the triangle lst LTRI.  The value
C              must be 6 if only the vertex indexes and
C              neighboring triangle indexes are to be
C              stored, or 9 if arc indexes are also to be
C              assigned and stored.  Refer to LTRI.
C
C The above parameters are not altered by this routine.
C
C       LTRI = Integer array of length at least NROW*NT,
C              where NT is at most 2N-4.  (A sufficient
C              length is 12N if NROW=6 or 18N if NROW=9.)
C
C On output:
C
C       NT = Number of triangles in the triangulation unless
C            IER .NE. 0, in which case NT = 0.  NT = 2N-NB-2
C            if NB .GE. 3 or 2N-4 if NB = 0, where NB is the
C            number of boundary nodes.
C
C       LTRI = NROW by NT array whose J-th column contains
C              the vertex nodal indexes (first three rows),
C              neighboring triangle indexes (second three
C              rows), and, if NROW = 9, arc indexes (last
C              three rows) associated with triangle J for
C              J = 1,...,NT.  The vertices are ordered
C              counterclockwise with the first vertex taken
C              to be the one with smallest index.  Thus,
C              LTRI(2,J) and LTRI(3,J) are larger than
C              LTRI(1,J) and index adjacent neighbors of
C              node LTRI(1,J).  For I = 1,2,3, LTRI(I+3,J)
C              and LTRI(I+6,J) index the triangle and arc,
C              respectively, which are opposite (not shared
C              by) node LTRI(I,J), with LTRI(I+3,J) = 0 if
C              LTRI(I+6,J) indexes a boundary arc.  Vertex
C              indexes range from 1 to N, triangle indexes
C              from 0 to NT, and, if included, arc indexes
C              from 1 to NA, where NA = 3N-NB-3 if NB .GE. 3
C              or 3N-6 if NB = 0.  The triangles are or-
C              dered on first (smallest) vertex indexes.
C
C       IER = Error indicator.
C             IER = 0 if no errors were encountered.
C             IER = 1 if N or NROW is outside its valid
C                     range on input.
C             IER = 2 if the triangulation data structure
C                     (LIST,LPTR,LEND) is invalid.  Note,
C                     however, that these arrays are not
C                     completely tested for validity.
C
C Modules required by TRLIST:  None
C
C Intrinsic function called by TRLIST:  ABS
C
C***********************************************************
C
      INTEGER I, I1, I2, I3, ISV, J, KA, KN, KT, LP, LP2,
     .        LPL, LPLN1, N1, N2, N3, NM2
      LOGICAL ARCS
C
C Local parameters:
C
C ARCS =     Logical variable with value TRUE iff are
C              indexes are to be stored
C I,J =      LTRI row indexes (1 to 3) associated with
C              triangles KT and KN, respectively
C I1,I2,I3 = Nodal indexes of triangle KN
C ISV =      Variable used to permute indexes I1,I2,I3
C KA =       Arc index and number of currently stored arcs
C KN =       Index of the triangle that shares arc I1-I2
C              with KT
C KT =       Triangle index and number of currently stored
C              triangles
C LP =       LIST pointer
C LP2 =      Pointer to N2 as a neighbor of N1
C LPL =      Pointer to the last neighbor of I1
C LPLN1 =    Pointer to the last neighbor of N1
C N1,N2,N3 = Nodal indexes of triangle KT
C NM2 =      N-2
C
C
C Test for invalid input parameters.
C
      IF (N .LT. 3  .OR.  (NROW .NE. 6  .AND.  NROW .NE. 9))
     .  GO TO 11
C
C Initialize parameters for loop on triangles KT = (N1,N2,
C   N3), where N1 < N2 and N1 < N3.
C
C   ARCS = TRUE iff arc indexes are to be stored.
C   KA,KT = Numbers of currently stored arcs and triangles.
C   NM2 = Upper bound on candidates for N1.
C
      ARCS = NROW .EQ. 9
      KA = 0
      KT = 0
      NM2 = N-2
C
C Loop on nodes N1.
C
      DO 9 N1 = 1,NM2
C
C Loop on pairs of adjacent neighbors (N2,N3).  LPLN1 points
C   to the last neighbor of N1, and LP2 points to N2.
C
        LPLN1 = LEND(N1)
        LP2 = LPLN1
    1     LP2 = LPTR(LP2)
          N2 = LIST(LP2)
          LP = LPTR(LP2)
          N3 = ABS(LIST(LP))
          IF (N2 .LT. N1  .OR.  N3 .LT. N1) GO TO 8
C
C Add a new triangle KT = (N1,N2,N3).
C
          KT = KT + 1
          LTRI(1,KT) = N1
          LTRI(2,KT) = N2
          LTRI(3,KT) = N3
C
C Loop on triangle sides (I2,I1) with neighboring triangles
C   KN = (I1,I2,I3).
C
          DO 7 I = 1,3
            IF (I .EQ. 1) THEN
              I1 = N3
              I2 = N2
            ELSEIF (I .EQ. 2) THEN
              I1 = N1
              I2 = N3
            ELSE
              I1 = N2
              I2 = N1
            ENDIF
C
C Set I3 to the neighbor of I1 that follows I2 unless
C   I2->I1 is a boundary arc.
C
            LPL = LEND(I1)
            LP = LPTR(LPL)
    2       IF (LIST(LP) .EQ. I2) GO TO 3
              LP = LPTR(LP)
              IF (LP .NE. LPL) GO TO 2
C
C   I2 is the last neighbor of I1 unless the data structure
C     is invalid.  Bypass the search for a neighboring
C     triangle if I2->I1 is a boundary arc.
C
            IF (ABS(LIST(LP)) .NE. I2) GO TO 12
            KN = 0
            IF (LIST(LP) .LT. 0) GO TO 6
C
C   I2->I1 is not a boundary arc, and LP points to I2 as
C     a neighbor of I1.
C
    3       LP = LPTR(LP)
            I3 = ABS(LIST(LP))
C
C Find J such that LTRI(J,KN) = I3 (not used if KN > KT),
C   and permute the vertex indexes of KN so that I1 is
C   smallest.
C
            IF (I1 .LT. I2  .AND.  I1 .LT. I3) THEN
              J = 3
            ELSEIF (I2 .LT. I3) THEN
              J = 2
              ISV = I1
              I1 = I2
              I2 = I3
              I3 = ISV
            ELSE
              J = 1
              ISV = I1
              I1 = I3
              I3 = I2
              I2 = ISV
            ENDIF
C
C Test for KN > KT (triangle index not yet assigned).
C
            IF (I1 .GT. N1) GO TO 7
C
C Find KN, if it exists, by searching the triangle lst in
C   reverse order.
C
            DO 4 KN = KT-1,1,-1
              IF (LTRI(1,KN) .EQ. I1  .AND.  LTRI(2,KN) .EQ.
     .            I2  .AND.  LTRI(3,KN) .EQ. I3) GO TO 5
    4         CONTINUE
            GO TO 7
C
C Store KT as a neighbor of KN.
C
    5       LTRI(J+3,KN) = KT
C
C Store KN as a neighbor of KT, and add a new arc KA.
C
    6       LTRI(I+3,KT) = KN
            IF (ARCS) THEN
              KA = KA + 1
              LTRI(I+6,KT) = KA
              IF (KN .NE. 0) LTRI(J+6,KN) = KA
            ENDIF
    7       CONTINUE
C
C Bottom of loop on triangles.
C
    8     IF (LP2 .NE. LPLN1) GO TO 1
    9     CONTINUE
C
C No errors encountered.
C
      NT = KT
      IER = 0
      RETURN
C
C Invalid input parameter.
C
   11 NT = 0
      IER = 1
      RETURN
C
C Invalid triangulation data structure:  I1 is a neighbor of
C   I2, but I2 is not a neighbor of I1.
C
   12 NT = 0
      IER = 2
      RETURN
      END
      SUBROUTINE TRLPRT (N,X,Y,Z,IFLAG,NROW,NT,LTRI,LOUT)
      INTEGER N, IFLAG, NROW, NT, LTRI(NROW,NT), LOUT
      DOUBLE PRECISION X(N), Y(N), Z(N)
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/02/98
C
C   This subroutine prints the triangle lst created by Sub-
C routine TRLIST and, optionally, the nodal coordinates
C (either latitude and longitude or Cartesian coordinates)
C on logical unit LOUT.  The numbers of boundary nodes,
C triangles, and arcs are also printed.
C
C
C On input:
C
C       N = Number of nodes in the triangulation.
C           3 .LE. N .LE. 9999.
C
C       X,Y,Z = Arrays of length N containing the Cartesian
C               coordinates of the nodes if IFLAG = 0, or
C               (X and Y only) arrays of length N containing
C               longitude and latitude, respectively, if
C               IFLAG > 0, or unused dummy parameters if
C               IFLAG < 0.
C
C       IFLAG = Nodal coordinate option indicator:
C               IFLAG = 0 if X, Y, and Z (assumed to contain
C                         Cartesian coordinates) are to be
C                         printed (to 6 decimal places).
C               IFLAG > 0 if only X and Y (assumed to con-
C                         tain longitude and latitude) are
C                         to be printed (to 6 decimal
C                         places).
C               IFLAG < 0 if only the adjacency lsts are to
C                         be printed.
C
C       NROW = Number of rows (entries per triangle) re-
C              served for the triangle lst LTRI.  The value
C              must be 6 if only the vertex indexes and
C              neighboring triangle indexes are stored, or 9
C              if arc indexes are also stored.
C
C       NT = Number of triangles in the triangulation.
C            1 .LE. NT .LE. 9999.
C
C       LTRI = NROW by NT array whose J-th column contains
C              the vertex nodal indexes (first three rows),
C              neighboring triangle indexes (second three
C              rows), and, if NROW = 9, arc indexes (last
C              three rows) associated with triangle J for
C              J = 1,...,NT.
C
C       LOUT = Logical unit number for output.  If LOUT is
C              not in the range 0 to 99, output is written
C              to unit 6.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C   The triangle lst and nodal coordinates (as specified by
C IFLAG) are written to unit LOUT.
C
C Modules required by TRLPRT:  None
C
C***********************************************************
C
      INTEGER I, K, LUN, NA, NB, NL, NLMAX, NMAX
      DATA    NMAX/9999/,  NLMAX/58/
C
C Local parameters:
C
C I =     DO-loop, nodal index, and row index for LTRI
C K =     DO-loop and triangle index
C LUN =   Logical unit number for output
C NA =    Number of triangulation arcs
C NB =    Number of boundary nodes
C NL =    Number of lines printed on the current page
C NLMAX = Maximum number of print lines per page (except
C           for the last page which may have two addi-
C           tional lines)
C NMAX =  Maximum value of N and NT (4-digit format)
C
      LUN = LOUT
      IF (LUN .LT. 0  .OR.  LUN .GT. 99) LUN = 6
C
C Print a heading and test for invalid input.
C
      WRITE (LUN,100) N
      NL = 3
      IF (N .LT. 3  .OR.  N .GT. NMAX  .OR.
     .    (NROW .NE. 6  .AND.  NROW .NE. 9)  .OR.
     .    NT .LT. 1  .OR.  NT .GT. NMAX) THEN
C
C Print an error message and exit.
C
        WRITE (LUN,110) N, NROW, NT
        RETURN
      ENDIF
      IF (IFLAG .EQ. 0) THEN
C
C Print X, Y, and Z.
C
        WRITE (LUN,101)
        NL = 6
        DO 1 I = 1,N
          IF (NL .GE. NLMAX) THEN
            WRITE (LUN,108)
            NL = 0
          ENDIF
          WRITE (LUN,103) I, X(I), Y(I), Z(I)
          NL = NL + 1
    1     CONTINUE
      ELSEIF (IFLAG .GT. 0) THEN
C
C Print X (longitude) and Y (latitude).
C
        WRITE (LUN,102)
        NL = 6
        DO 2 I = 1,N
          IF (NL .GE. NLMAX) THEN
            WRITE (LUN,108)
            NL = 0
          ENDIF
          WRITE (LUN,104) I, X(I), Y(I)
          NL = NL + 1
    2     CONTINUE
      ENDIF
C
C Print the triangulation LTRI.
C
      IF (NL .GT. NLMAX/2) THEN
        WRITE (LUN,108)
        NL = 0
      ENDIF
      IF (NROW .EQ. 6) THEN
        WRITE (LUN,105)
      ELSE
        WRITE (LUN,106)
      ENDIF
      NL = NL + 5
      DO 3 K = 1,NT
        IF (NL .GE. NLMAX) THEN
          WRITE (LUN,108)
          NL = 0
        ENDIF
        WRITE (LUN,107) K, (LTRI(I,K), I = 1,NROW)
        NL = NL + 1
    3   CONTINUE
C
C Print NB, NA, and NT (boundary nodes, arcs, and
C   triangles).
C
      NB = 2*N - NT - 2
      IF (NB .LT. 3) THEN
        NB = 0
        NA = 3*N - 6
      ELSE
        NA = NT + N - 1
      ENDIF
      WRITE (LUN,109) NB, NA, NT
      RETURN
C
C Print formats:
C
  100 FORMAT (///18X,'STRIPACK (TRLIST) Output,  N = ',I4)
  101 FORMAT (//8X,'Node',10X,'X(Node)',10X,'Y(Node)',10X,
     .        'Z(Node)'//)
  102 FORMAT (//16X,'Node',8X,'Longitude',9X,'Latitude'//)
  103 FORMAT (8X,I4,3E17.6)
  104 FORMAT (16X,I4,2E17.6)
  105 FORMAT (//1X,'Triangle',8X,'Vertices',12X,'Neighbors'/
     .        4X,'KT',7X,'N1',5X,'N2',5X,'N3',4X,'KT1',4X,
     .        'KT2',4X,'KT3'/)
  106 FORMAT (//1X,'Triangle',8X,'Vertices',12X,'Neighbors',
     .        14X,'Arcs'/
     .        4X,'KT',7X,'N1',5X,'N2',5X,'N3',4X,'KT1',4X,
     .        'KT2',4X,'KT3',4X,'KA1',4X,'KA2',4X,'KA3'/)
  107 FORMAT (2X,I4,2X,6(3X,I4),3(2X,I5))
  108 FORMAT (///)
  109 FORMAT (/1X,'NB = ',I4,' Boundary Nodes',5X,
     .        'NA = ',I5,' Arcs',5X,'NT = ',I5,
     .        ' Triangles')
  110 FORMAT (//1X,10X,'*** Invalid Parameter:  N =',I5,
     .        ', NROW =',I5,', NT =',I5,' ***')
      END
      SUBROUTINE TRMESH (N,X,Y,Z,LIST,LPTR,LEND,IER)
      INTEGER N, LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N), LNEW, NEAR(N),
     .        NEXT(N), IER
      DOUBLE PRECISION    X(N), Y(N), Z(N), DIST(N)
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/08/99
C
C   This subroutine creates a Delaunay triangulation of a
C set of N arbitrarily distributed points, referred to as
C nodes, on the surface of the unit sphere.  The Delaunay
C triangulation is defined as a set of (spherical) triangles
C with the following five properties:
C
C  1)  The triangle vertices are nodes.
C  2)  No triangle contains a node other than its vertices.
C  3)  The interiors of the triangles are pairwise disjoint.
C  4)  The union of triangles is the convex hull of the set
C        of nodes (the smallest convex set that contains
C        the nodes).  If the nodes are not contained in a
C        single hemisphere, their convex hull is the en-
C        tire sphere and there are no boundary nodes.
C        Otherwise, there are at least three boundary nodes.
C  5)  The interior of the circumcircle of each triangle
C        contains no node.
C
C The first four properties define a triangulation, and the
C last property results in a triangulation which is as close
C as possible to equiangular in a certain sense and which is
C uniquely defined unless four or more nodes lie in a common
C plane.  This property makes the triangulation well-suited
C for solving closest-point problems and for triangle-based
C interpolation.
C
C   Provided the nodes are randomly ordered, the algorithm
C has expected time complexity O(N*log(N)) for most nodal
C distributions.  Note, however, that the complexity may be
C as high as O(N**2) if, for example, the nodes are ordered
C on increasing latitude.
C
C   Spherical coordinates (latitude and longitude) may be
C converted to Cartesian coordinates by Subroutine TRANS.
C
C   The following is a lst of the software package modules
C which a user may wish to call directly:
C
C  ADDNOD - Updates the triangulation by appending a new
C             node.
C
C  AREAS  - Returns the area of a spherical triangle.
C
C  BNODES - Returns an array containing the indexes of the
C             boundary nodes (if any) in counterclockwise
C             order.  Counts of boundary nodes, triangles,
C             and arcs are also returned.
C
C  CIRCUM - Returns the circumcenter of a spherical trian-
C             gle.
C
C  CRLIST - Returns the set of triangle circumcenters
C             (Voronoi vertices) and circumradii associated
C             with a triangulation.
C
C  DELARC - Deletes a boundary arc from a triangulation.
C
C  DELNOD - Updates the triangulation with a nodal deletion.
C
C  EDGE   - Forces an arbitrary pair of nodes to be connec-
C             ted by an arc in the triangulation.
C
C  GETNP  - Determines the ordered sequence of L closest
C             nodes to a given node, along with the associ-
C             ated distances.
C
C  INSIDE - Locates a point relative to a polygon on the
C             surface of the sphere.
C
C  INTRSC - Returns the point of intersection between a
C             pair of great circle arcs.
C
C  JRAND  - Generates a uniformly distributed pseudo-random
C             integer.
C
C  LEFT   - Locates a point relative to a great circle.
C
C  NEARND - Returns the index of the nearest node to an
C             arbitrary point, along with its squared
C             distance.
C
C  SCOORD - Converts a point from Cartesian coordinates to
C             spherical coordinates.
C
C  STORE  - Forces a value to be stored in main memory so
C             that the precision of floating point numbers
C             in memory locations rather than registers is
C             computed.
C
C  TRANS  - Transforms spherical coordinates into Cartesian
C             coordinates on the unit sphere for input to
C             Subroutine TRMESH.
C
C  TRLIST - Converts the triangulation data structure to a
C             triangle lst more suitable for use in a fin-
C             ite element code.
C
C  TRLPRT - Prints the triangle lst created by Subroutine
C             TRLIST.
C
C  TRMESH - Creates a Delaunay triangulation of a set of
C             nodes.
C
C  TRPLOT - Creates a level-2 Encapsulated Postscript (EPS)
C             file containing a triangulation plot.
C
C  TRPRNT - Prints the triangulation data structure and,
C             optionally, the nodal coordinates.
C
C  VRPLOT - Creates a level-2 Encapsulated Postscript (EPS)
C             file containing a Voronoi diagram plot.
C
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y,Z = Arrays of length N containing the Cartesian
C               coordinates of distinct nodes.  (X(K),Y(K),
C               Z(K)) is referred to as node K, and K is re-
C               ferred to as a nodal index.  It is required
C               that X(K)**2 + Y(K)**2 + Z(K)**2 = 1 for all
C               K.  The first three nodes must not be col-
C               linear (lie on a common great circle).
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR = Arrays of length at least 6N-12.
C
C       LEND = Array of length at least N.
C
C       NEAR,NEXT,DIST = Work space arrays of length at
C                        least N.  The space is used to
C                        efficiently determine the nearest
C                        triangulation node to each un-
C                        processed node for use by ADDNOD.
C
C On output:
C
C       LIST = Set of nodal indexes which, along with LPTR,
C              LEND, and LNEW, define the triangulation as a
C              set of N adjacency lsts -- counterclockwise-
C              ordered sequences of neighboring nodes such
C              that the first and last neighbors of a bound-
C              ary node are boundary nodes (the first neigh-
C              bor of an interior node is arbitrary).  In
C              order to distinguish between interior and
C              boundary nodes, the last neighbor of each
C              boundary node is represented by the negative
C              of its index.
C
C       LPTR = Set of pointers (LIST indexes) in one-to-one
C              correspondence with the elements of LIST.
C              LIST(LPTR(I)) indexes the node which follows
C              LIST(I) in cyclical counterclockwise order
C              (the first neighbor follows the last neigh-
C              bor).
C
C       LEND = Set of pointers to adjacency lsts.  LEND(K)
C              points to the last neighbor of node K for
C              K = 1,...,N.  Thus, LIST(LEND(K)) < 0 if and
C              only if K is a boundary node.
C
C       LNEW = Pointer to the first empty location in LIST
C              and LPTR (lst length plus one).  LIST, LPTR,
C              LEND, and LNEW are not altered if IER < 0,
C              and are incomplete if IER > 0.
C
C       NEAR,NEXT,DIST = Garbage.
C
C       IER = Error indicator:
C             IER =  0 if no errors were encountered.
C             IER = -1 if N < 3 on input.
C             IER = -2 if the first three nodes are
C                      collinear.
C             IER =  L if nodes L and M coincide for some
C                      M > L.  The data structure represents
C                      a triangulation of nodes 1 to M-1 in
C                      this case.
C
C Modules required by TRMESH:  ADDNOD, BDYADD, COVSPH,
C                                INSERT, INTADD, JRAND,
C                                LEFT, LSTPTR, STORE, SWAP,
C                                SWPTST, TRFIND
C
C Intrinsic function called by TRMESH:  ABS
C
C***********************************************************
C
      INTEGER I, I0, J, K, LP, LPL, NEXTI, NN
      LOGICAL LEFT
      DOUBLE PRECISION    D, D1, D2, D3
C
C Local parameters:
C
C D =        (Negative cosine of) distance from node K to
C              node I
C D1,D2,D3 = Distances from node K to nodes 1, 2, and 3,
C              respectively
C I,J =      Nodal indexes
C I0 =       Index of the node preceding I in a sequence of
C              unprocessed nodes:  I = NEXT(I0)
C K =        Index of node to be added and DO-loop index:
C              K > 3
C LP =       LIST index (pointer) of a neighbor of K
C LPL =      Pointer to the last neighbor of K
C NEXTI =    NEXT(I)
C NN =       Local copy of N
C
      NN = N
      IF (NN .LT. 3) THEN
        IER = -1
        RETURN
      ENDIF
C
C Store the first triangle in the linked lst.
C
      IF ( .NOT. LEFT (X(1),Y(1),Z(1),X(2),Y(2),Z(2),
     .                 X(3),Y(3),Z(3)) ) THEN
C
C   The first triangle is (3,2,1) = (2,1,3) = (1,3,2).
C
        LIST(1) = 3
        LPTR(1) = 2
        LIST(2) = -2
        LPTR(2) = 1
        LEND(1) = 2
C
        LIST(3) = 1
        LPTR(3) = 4
        LIST(4) = -3
        LPTR(4) = 3
        LEND(2) = 4
C
        LIST(5) = 2
        LPTR(5) = 6
        LIST(6) = -1
        LPTR(6) = 5
        LEND(3) = 6
C
      ELSEIF ( .NOT. LEFT(X(2),Y(2),Z(2),X(1),Y(1),Z(1),
     .                    X(3),Y(3),Z(3)) )
     .       THEN
C
C   The first triangle is (1,2,3):  3 Strictly Left 1->2,
C     i.e., node 3 lies in the left hemisphere defined by
C     arc 1->2.
C
        LIST(1) = 2
        LPTR(1) = 2
        LIST(2) = -3
        LPTR(2) = 1
        LEND(1) = 2
C
        LIST(3) = 3
        LPTR(3) = 4
        LIST(4) = -1
        LPTR(4) = 3
        LEND(2) = 4
C
        LIST(5) = 1
        LPTR(5) = 6
        LIST(6) = -2
        LPTR(6) = 5
        LEND(3) = 6
C
      ELSE
C
C   The first three nodes are collinear.
C
        IER = -2
        RETURN
      ENDIF
C
C Initialize LNEW and test for N = 3.
C
      LNEW = 7
      IF (NN .EQ. 3) THEN
        IER = 0
        RETURN
      ENDIF
C
C A nearest-node data structure (NEAR, NEXT, and DIST) is
C   used to obtain an expected-time (N*log(N)) incremental
C   algorithm by enabling constant search time for locating
C   each new node in the triangulation.
C
C For each unprocessed node K, NEAR(K) is the index of the
C   triangulation node closest to K (used as the starting
C   point for the search in Subroutine TRFIND) and DIST(K)
C   is an increasing function of the arc length (angular
C   distance) between nodes K and NEAR(K):  -Cos(a) for arc
C   length a.
C
C Since it is necessary to efficiently find the subset of
C   unprocessed nodes associated with each triangulation
C   node J (those that have J as their NEAR entries), the
C   subsets are stored in NEAR and NEXT as follows:  for
C   each node J in the triangulation, I = NEAR(J) is the
C   first unprocessed node in J's set (with I = 0 if the
C   set is empty), L = NEXT(I) (if I > 0) is the second,
C   NEXT(L) (if L > 0) is the third, etc.  The nodes in each
C   set are initially ordered by increasing indexes (which
C   maximizes efficiency) but that ordering is not main-
C   tained as the data structure is updated.
C
C Initialize the data structure for the single triangle.
C
      NEAR(1) = 0
      NEAR(2) = 0
      NEAR(3) = 0
      DO 1 K = NN,4,-1
        D1 = -(X(K)*X(1) + Y(K)*Y(1) + Z(K)*Z(1))
        D2 = -(X(K)*X(2) + Y(K)*Y(2) + Z(K)*Z(2))
        D3 = -(X(K)*X(3) + Y(K)*Y(3) + Z(K)*Z(3))
        IF (D1 .LE. D2  .AND.  D1 .LE. D3) THEN
          NEAR(K) = 1
          DIST(K) = D1
          NEXT(K) = NEAR(1)
          NEAR(1) = K
        ELSEIF (D2 .LE. D1  .AND.  D2 .LE. D3) THEN
          NEAR(K) = 2
          DIST(K) = D2
          NEXT(K) = NEAR(2)
          NEAR(2) = K
        ELSE
          NEAR(K) = 3
          DIST(K) = D3
          NEXT(K) = NEAR(3)
          NEAR(3) = K
        ENDIF
    1   CONTINUE
C
C Add the remaining nodes
C
      DO 6 K = 4,NN
        CALL ADDNOD (NEAR(K),K,X,Y,Z, LIST,LPTR,LEND,
     .               LNEW, IER)
        IF (IER .NE. 0) RETURN
C
C Remove K from the set of unprocessed nodes associated
C   with NEAR(K).
C
        I = NEAR(K)
        IF (NEAR(I) .EQ. K) THEN
          NEAR(I) = NEXT(K)
        ELSE
          I = NEAR(I)
    2     I0 = I
            I = NEXT(I0)
            IF (I .NE. K) GO TO 2
          NEXT(I0) = NEXT(K)
        ENDIF
        NEAR(K) = 0
C
C Loop on neighbors J of node K.
C
        LPL = LEND(K)
        LP = LPL
    3   LP = LPTR(LP)
          J = ABS(LIST(LP))
C
C Loop on elements I in the sequence of unprocessed nodes
C   associated with J:  K is a candidate for replacing J
C   as the nearest triangulation node to I.  The next value
C   of I in the sequence, NEXT(I), must be saved before I
C   is moved because it is altered by adding I to K's set.
C
          I = NEAR(J)
    4     IF (I .EQ. 0) GO TO 5
          NEXTI = NEXT(I)
C
C Test for the distance from I to K less than the distance
C   from I to J.
C
          D = -(X(I)*X(K) + Y(I)*Y(K) + Z(I)*Z(K))
          IF (D .LT. DIST(I)) THEN
C
C Replace J by K as the nearest triangulation node to I:
C   update NEAR(I) and DIST(I), and remove I from J's set
C   of unprocessed nodes and add it to K's set.
C
            NEAR(I) = K
            DIST(I) = D
            IF (I .EQ. NEAR(J)) THEN
              NEAR(J) = NEXTI
            ELSE
              NEXT(I0) = NEXTI
            ENDIF
            NEXT(I) = NEAR(K)
            NEAR(K) = I
          ELSE
            I0 = I
          ENDIF
C
C Bottom of loop on I.
C
          I = NEXTI
          GO TO 4
C
C Bottom of loop on neighbors J.
C
    5     IF (LP .NE. LPL) GO TO 3
    6   CONTINUE
      RETURN
      END
      SUBROUTINE TRPLOT (LUN,PLTSIZ,ELAT,ELON,A,N,X,Y,Z,
     .                   LIST,LPTR,LEND,TITLE,NUMBR, IER)
      CHARACTER*(*) TITLE
      INTEGER   LUN, N, LIST(*), LPTR(*), LEND(N), IER
      LOGICAL   NUMBR
      DOUBLE PRECISION      PLTSIZ, ELAT, ELON, A, X(N), Y(N), Z(N)
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/16/98
C
C   This subroutine creates a level-2 Encapsulated Post-
C script (EPS) file containing a graphical display of a
C triangulation of a set of nodes on the unit sphere.  The
C visible nodes are projected onto the plane that contains
C the origin and has normal defined by a user-specified eye-
C position.  Projections of adjacent (visible) nodes are
C connected by line segments.
C
C
C On input:
C
C       LUN = Logical unit number in the range 0 to 99.
C             The unit should be opened with an appropriate
C             file name before the call to this routine.
C
C       PLTSIZ = Plot size in inches.  A circular window in
C                the projection plane is mapped to a circu-
C                lar viewport with diameter equal to .88*
C                PLTSIZ (leaving room for labels outside the
C                viewport).  The viewport is centered on the
C                8.5 by 11 inch page, and its boundary is
C                drawn.  1.0 .LE. PLTSIZ .LE. 8.5.
C
C       ELAT,ELON = Latitude and longitude (in degrees) of
C                   the center of projection E (the center
C                   of the plot).  The projection plane is
C                   the plane that contains the origin and
C                   has E as unit normal.  In a rotated
C                   coordinate system for which E is the
C                   north pole, the projection plane con-
C                   tains the equator, and only northern
C                   hemisphere nodes are visible (from the
C                   point at infinity in the direction E).
C                   These are projected orthogonally onto
C                   the projection plane (by zeroing the z-
C                   component in the rotated coordinate
C                   system).  ELAT and ELON must be in the
C                   range -90 to 90 and -180 to 180, respec-
C                   tively.
C
C       A = Angular distance in degrees from E to the boun-
C           dary of a circular window against which the
C           triangulation is clipped.  The projected window
C           is a disk of radius r = Sin(A) centered at the
C           origin, and only visible nodes whose projections
C           are within distance r of the origin are included
C           in the plot.  Thus, if A = 90, the plot includes
C           the entire hemisphere centered at E.  0 .LT. A
C           .LE. 90.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y,Z = Arrays of length N containing the Cartesian
C               coordinates of the nodes (unit vectors).
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to Subroutine
C                        TRMESH.
C
C       TITLE = Type CHARACTER variable or constant contain-
C               ing a string to be centered above the plot.
C               The string must be enclosed in parentheses;
C               i.e., the first and last characters must be
C               '(' and ')', respectively, but these are not
C               displayed.  TITLE may have at most 80 char-
C               acters including the parentheses.
C
C       NUMBR = Option indicator:  If NUMBR = TRUE, the
C               nodal indexes are plotted next to the nodes.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if LUN, PLTSIZ, or N is outside its
C                     valid range.
C             IER = 2 if ELAT, ELON, or A is outside its
C                     valid range.
C             IER = 3 if an error was encountered in writing
C                     to unit LUN.
C
C   The values in the data statement below may be altered
C in order to modify various plotting options.
C
C Modules required by TRPLOT:  None
C
C Intrinsic functions called by TRPLOT:  ABS, ATAN, COS,
C                                          NINT, REAL, SIN,
C                                          SQRT
C
C***********************************************************
C
      INTEGER IPX1, IPX2, IPY1, IPY2, IR, LP, LPL, N0, N1
      LOGICAL ANNOT
      DOUBLE PRECISION    CF, CT, EX, EY, EZ, FSIZN, FSIZT, R11, R12,
     .        R21, R22, R23, SF, T, TX, TY, WR, WRS, X0, X1,
     .        Y0, Y1, Z0, Z1
C
      DATA    ANNOT/.TRUE./,  FSIZN/10.0/,  FSIZT/16.0/
C
C Local parameters:
C
C ANNOT =     Logical variable with value TRUE iff the plot
C               is to be annotated with the values of ELAT,
C               ELON, and A
C CF =        Conversion factor for degrees to radians
C CT =        Cos(ELAT)
C EX,EY,EZ =  Cartesian coordinates of the eye-position E
C FSIZN =     Font size in points for labeling nodes with
C               their indexes if NUMBR = TRUE
C FSIZT =     Font size in points for the title (and
C               annotation if ANNOT = TRUE)
C IPX1,IPY1 = X and y coordinates (in points) of the lower
C               left corner of the bounding box or viewport
C               box
C IPX2,IPY2 = X and y coordinates (in points) of the upper
C               right corner of the bounding box or viewport
C               box
C IR =        Half the width (height) of the bounding box or
C               viewport box in points -- viewport radius
C LP =        LIST index (pointer)
C LPL =       Pointer to the last neighbor of N0
C N0 =        Index of a node whose incident arcs are to be
C               drawn
C N1 =        Neighbor of N0
C R11...R23 = Components of the first two rows of a rotation
C               that maps E to the north pole (0,0,1)
C SF =        Scale factor for mapping world coordinates
C               (window coordinates in [-WR,WR] X [-WR,WR])
C               to viewport coordinates in [IPX1,IPX2] X
C               [IPY1,IPY2]
C T =         Temporary variable
C TX,TY =     Translation vector for mapping world coordi-
C               nates to viewport coordinates
C WR =        Window radius r = Sin(A)
C WRS =       WR**2
C X0,Y0,Z0 =  Coordinates of N0 in the rotated coordinate
C               system or label location (X0,Y0)
C X1,Y1,Z1 =  Coordinates of N1 in the rotated coordinate
C               system or intersection of edge N0-N1 with
C               the equator (in the rotated coordinate
C               system)
C
C
C Test for invalid parameters.
C
      IF (LUN .LT. 0  .OR.  LUN .GT. 99  .OR.
     .    PLTSIZ .LT. 1.0  .OR.  PLTSIZ .GT. 8.5  .OR.
     .    N .LT. 3)
     .  GO TO 11
      IF (ABS(ELAT) .GT. 90.0  .OR.  ABS(ELON) .GT. 180.0
     .    .OR.  A .GT. 90.0) GO TO 12
C
C Compute a conversion factor CF for degrees to radians
C   and compute the window radius WR.
C
      CF = ATAN(1.0)/45.0
      WR = SIN(CF*A)
      WRS = WR*WR
C
C Compute the lower left (IPX1,IPY1) and upper right
C   (IPX2,IPY2) corner coordinates of the bounding box.
C   The coordinates, specified in default user space units
C   (points, at 72 points/inch with origin at the lower
C   left corner of the page), are chosen to preserve the
C   square aspect ratio, and to center the plot on the 8.5
C   by 11 inch page.  The center of the page is (306,396),
C   and IR = PLTSIZ/2 in points.
C
      IR = NINT(36.0*PLTSIZ)
      IPX1 = 306 - IR
      IPX2 = 306 + IR
      IPY1 = 396 - IR
      IPY2 = 396 + IR
C
C Output header comments.
C
      WRITE (LUN,100,ERR=13) IPX1, IPY1, IPX2, IPY2
  100 FORMAT ('%!PS-Adobe-3.0 EPSF-3.0'/
     .        '%%BoundingBox:',4I4/
     .        '%%Title:  Triangulation'/
     .        '%%Creator:  STRIPACK'/
     .        '%%EndComments')
C
C Set (IPX1,IPY1) and (IPX2,IPY2) to the corner coordinates
C   of a viewport box obtained by shrinking the bounding box
C   by 12% in each dimension.
C
      IR = NINT(0.88*REAL(IR))
      IPX1 = 306 - IR
      IPX2 = 306 + IR
      IPY1 = 396 - IR
      IPY2 = 396 + IR
C
C Set the line thickness to 2 points, and draw the
C   viewport boundary.
C
      T = 2.0
      WRITE (LUN,110,ERR=13) T
      WRITE (LUN,120,ERR=13) IR
      WRITE (LUN,130,ERR=13)
  110 FORMAT (F12.6,' setlinewidth')
  120 FORMAT ('306 396 ',I3,' 0 360 arc')
  130 FORMAT ('stroke')
C
C Set up an affine mapping from the window box [-WR,WR] X
C   [-WR,WR] to the viewport box.
C
      SF = REAL(IR)/WR
      TX = IPX1 + SF*WR
      TY = IPY1 + SF*WR
      WRITE (LUN,140,ERR=13) TX, TY, SF, SF
  140 FORMAT (2F12.6,' translate'/
     .        2F12.6,' scale')
C
C The line thickness must be changed to reflect the new
C   scaling which is applied to all subsequent output.
C   Set it to 1.0 point.
C
      T = 1.0/SF
      WRITE (LUN,110,ERR=13) T
C
C Save the current graphics state, and set the clip path to
C   the boundary of the window.
C
      WRITE (LUN,150,ERR=13)
      WRITE (LUN,160,ERR=13) WR
      WRITE (LUN,170,ERR=13)
  150 FORMAT ('gsave')
  160 FORMAT ('0 0 ',F12.6,' 0 360 arc')
  170 FORMAT ('clip newpath')
C
C Compute the Cartesian coordinates of E and the components
C   of a rotation R which maps E to the north pole (0,0,1).
C   R is taken to be a rotation about the z-axis (into the
C   yz-plane) followed by a rotation about the x-axis chosen
C   so that the view-up direction is (0,0,1), or (-1,0,0) if
C   E is the north or south pole.
C
C           ( R11  R12  0   )
C       R = ( R21  R22  R23 )
C           ( EX   EY   EZ  )
C
      T = CF*ELON
      CT = COS(CF*ELAT)
      EX = CT*COS(T)
      EY = CT*SIN(T)
      EZ = SIN(CF*ELAT)
      IF (CT .NE. 0.0) THEN
        R11 = -EY/CT
        R12 = EX/CT
      ELSE
        R11 = 0.0
        R12 = 1.0
      ENDIF
      R21 = -EZ*R12
      R22 = EZ*R11
      R23 = CT
C
C Loop on visible nodes N0 that project to points (X0,Y0) in
C   the window.
C
      DO 3 N0 = 1,N
        Z0 = EX*X(N0) + EY*Y(N0) + EZ*Z(N0)
        IF (Z0 .LT. 0.) GO TO 3
        X0 = R11*X(N0) + R12*Y(N0)
        Y0 = R21*X(N0) + R22*Y(N0) + R23*Z(N0)
        IF (X0*X0 + Y0*Y0 .GT. WRS) GO TO 3
        LPL = LEND(N0)
        LP = LPL
C
C Loop on neighbors N1 of N0.  LPL points to the last
C   neighbor of N0.  Copy the components of N1 into P.
C
    1   LP = LPTR(LP)
          N1 = ABS(LIST(LP))
          X1 = R11*X(N1) + R12*Y(N1)
          Y1 = R21*X(N1) + R22*Y(N1) + R23*Z(N1)
          Z1 = EX*X(N1) + EY*Y(N1) + EZ*Z(N1)
          IF (Z1 .LT. 0.) THEN
C
C   N1 is a 'southern hemisphere' point.  Move it to the
C     intersection of edge N0-N1 with the equator so that
C     the edge is clipped properly.  Z1 is implicitly set
C     to 0.
C
            X1 = Z0*X1 - Z1*X0
            Y1 = Z0*Y1 - Z1*Y0
            T = SQRT(X1*X1+Y1*Y1)
            X1 = X1/T
            Y1 = Y1/T
          ENDIF
C
C   If node N1 is in the window and N1 < N0, bypass edge
C     N0->N1 (since edge N1->N0 has already been drawn).
C
          IF ( Z1 .GE. 0.0  .AND.  X1*X1 + Y1*Y1 .LE. WRS
     .         .AND.  N1 .LT. N0 ) GO TO 2
C
C   Add the edge to the path.
C
          WRITE (LUN,180,ERR=13) X0, Y0, X1, Y1
  180     FORMAT (2F12.6,' moveto',2F12.6,' lineto')
C
C Bottom of loops.
C
    2     IF (LP .NE. LPL) GO TO 1
    3   CONTINUE
C
C Paint the path and restore the saved graphics state (with
C   no clip path).
C
      WRITE (LUN,130,ERR=13)
      WRITE (LUN,190,ERR=13)
  190 FORMAT ('grestore')
      IF (NUMBR) THEN
C
C Nodes in the window are to be labeled with their indexes.
C   Convert FSIZN from points to world coordinates, and
C   output the commands to select a font and scale it.
C
        T = FSIZN/SF
        WRITE (LUN,200,ERR=13) T
  200   FORMAT ('/Helvetica findfont'/
     .          F12.6,' scalefont setfont')
C
C Loop on visible nodes N0 that project to points (X0,Y0) in
C   the window.
C
        DO 4 N0 = 1,N
          IF (EX*X(N0) + EY*Y(N0) + EZ*Z(N0) .LT. 0.)
     .      GO TO 4
          X0 = R11*X(N0) + R12*Y(N0)
          Y0 = R21*X(N0) + R22*Y(N0) + R23*Z(N0)
          IF (X0*X0 + Y0*Y0 .GT. WRS) GO TO 4
C
C   Move to (X0,Y0) and draw the label N0.  The first char-
C     acter will will have its lower left corner about one
C     character width to the right of the nodal position.
C
          WRITE (LUN,210,ERR=13) X0, Y0
          WRITE (LUN,220,ERR=13) N0
  210     FORMAT (2F12.6,' moveto')
  220     FORMAT ('(',I3,') show')
    4     CONTINUE
      ENDIF
C
C Convert FSIZT from points to world coordinates, and output
C   the commands to select a font and scale it.
C
      T = FSIZT/SF
      WRITE (LUN,200,ERR=13) T
C
C Display TITLE centered above the plot:
C
      Y0 = WR + 3.0*T
      WRITE (LUN,230,ERR=13) TITLE, Y0
  230 FORMAT (A80/'  stringwidth pop 2 div neg ',F12.6,
     .        ' moveto')
      WRITE (LUN,240,ERR=13) TITLE
  240 FORMAT (A80/'  show')
      IF (ANNOT) THEN
C
C Display the window center and radius below the plot.
C
        X0 = -WR
        Y0 = -WR - 50.0/SF
        WRITE (LUN,210,ERR=13) X0, Y0
        WRITE (LUN,250,ERR=13) ELAT, ELON
        Y0 = Y0 - 2.0*T
        WRITE (LUN,210,ERR=13) X0, Y0
        WRITE (LUN,260,ERR=13) A
  250   FORMAT ('(Window center:  ELAT = ',F7.2,
     .          ',  ELON = ',F8.2,') show')
  260   FORMAT ('(Angular extent:  A = ',F5.2,') show')
      ENDIF
C
C Paint the path and output the showpage command and
C   end-of-file indicator.
C
      WRITE (LUN,270,ERR=13)
  270 FORMAT ('stroke'/
     .        'showpage'/
     .        '%%EOF')
C
C HP's interpreters require a one-byte End-of-PostScript-Job
C   indicator (to eliminate a timeout error message):
C   ASCII 4.
C
      WRITE (LUN,280,ERR=13) CHAR(4)
  280 FORMAT (A1)
C
C No error encountered.
C
      IER = 0
      RETURN
C
C Invalid input parameter LUN, PLTSIZ, or N.
C
   11 IER = 1
      RETURN
C
C Invalid input parameter ELAT, ELON, or A.
C
   12 IER = 2
      RETURN
C
C Error writing to unit LUN.
C
   13 IER = 3
      RETURN
      END
      SUBROUTINE TRPRNT (N,X,Y,Z,IFLAG,LIST,LPTR,LEND,LOUT)
      INTEGER N, IFLAG, LIST(*), LPTR(*), LEND(N), LOUT
      DOUBLE PRECISION    X(N), Y(N), Z(N)
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/25/98
C
C   This subroutine prints the triangulation adjacency lsts
C created by Subroutine TRMESH and, optionally, the nodal
C coordinates (either latitude and longitude or Cartesian
C coordinates) on logical unit LOUT.  The lst of neighbors
C of a boundary node is followed by index 0.  The numbers of
C boundary nodes, triangles, and arcs are also printed.
C
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 3
C           and N .LE. 9999.
C
C       X,Y,Z = Arrays of length N containing the Cartesian
C               coordinates of the nodes if IFLAG = 0, or
C               (X and Y only) arrays of length N containing
C               longitude and latitude, respectively, if
C               IFLAG > 0, or unused dummy parameters if
C               IFLAG < 0.
C
C       IFLAG = Nodal coordinate option indicator:
C               IFLAG = 0 if X, Y, and Z (assumed to contain
C                         Cartesian coordinates) are to be
C                         printed (to 6 decimal places).
C               IFLAG > 0 if only X and Y (assumed to con-
C                         tain longitude and latitude) are
C                         to be printed (to 6 decimal
C                         places).
C               IFLAG < 0 if only the adjacency lsts are to
C                         be printed.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to Subroutine
C                        TRMESH.
C
C       LOUT = Logical unit for output.  If LOUT is not in
C              the range 0 to 99, output is written to
C              logical unit 6.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C   The adjacency lsts and nodal coordinates (as specified
C by IFLAG) are written to unit LOUT.
C
C Modules required by TRPRNT:  None
C
C***********************************************************
C
      INTEGER I, INC, K, LP, LPL, LUN, NA, NABOR(400), NB,
     .        ND, NL, NLMAX, NMAX, NODE, NN, NT
      DATA  NMAX/9999/,  NLMAX/58/
C
C Local parameters:
C
C I =     NABOR index (1 to K)
C INC =   Increment for NL associated with an adjacency lst
C K =     Counter and number of neighbors of NODE
C LP =    LIST pointer of a neighbor of NODE
C LPL =   Pointer to the last neighbor of NODE
C LUN =   Logical unit for output (copy of LOUT)
C NA =    Number of arcs in the triangulation
C NABOR = Array containing the adjacency lst associated
C           with NODE, with zero appended if NODE is a
C           boundary node
C NB =    Number of boundary nodes encountered
C ND =    Index of a neighbor of NODE (or negative index)
C NL =    Number of lines that have been printed on the
C           current page
C NLMAX = Maximum number of print lines per page (except
C           for the last page which may have two addi-
C           tional lines)
C NMAX =  Upper bound on N (allows 4-digit indexes)
C NODE =  Index of a node and DO-loop index (1 to N)
C NN =    Local copy of N
C NT =    Number of triangles in the triangulation
C
      NN = N
      LUN = LOUT
      IF (LUN .LT. 0  .OR.  LUN .GT. 99) LUN = 6
C
C Print a heading and test the range of N.
C
      WRITE (LUN,100) NN
      IF (NN .LT. 3  .OR.  NN .GT. NMAX) THEN
C
C N is outside its valid range.
C
        WRITE (LUN,110)
        RETURN
      ENDIF
C
C Initialize NL (the number of lines printed on the current
C   page) and NB (the number of boundary nodes encountered).
C
      NL = 6
      NB = 0
      IF (IFLAG .LT. 0) THEN
C
C Print LIST only.  K is the number of neighbors of NODE
C   that have been stored in NABOR.
C
        WRITE (LUN,101)
        DO 2 NODE = 1,NN
          LPL = LEND(NODE)
          LP = LPL
          K = 0
C
    1     K = K + 1
            LP = LPTR(LP)
            ND = LIST(LP)
            NABOR(K) = ND
            IF (LP .NE. LPL) GO TO 1
          IF (ND .LE. 0) THEN
C
C   NODE is a boundary node.  Correct the sign of the last
C     neighbor, add 0 to the end of the lst, and increment
C     NB.
C
            NABOR(K) = -ND
            K = K + 1
            NABOR(K) = 0
            NB = NB + 1
          ENDIF
C
C   Increment NL and print the lst of neighbors.
C
          INC = (K-1)/14 + 2
          NL = NL + INC
          IF (NL .GT. NLMAX) THEN
            WRITE (LUN,108)
            NL = INC
          ENDIF
          WRITE (LUN,104) NODE, (NABOR(I), I = 1,K)
          IF (K .NE. 14) WRITE (LUN,107)
    2     CONTINUE
      ELSEIF (IFLAG .GT. 0) THEN
C
C Print X (longitude), Y (latitude), and LIST.
C
        WRITE (LUN,102)
        DO 4 NODE = 1,NN
          LPL = LEND(NODE)
          LP = LPL
          K = 0
C
    3     K = K + 1
            LP = LPTR(LP)
            ND = LIST(LP)
            NABOR(K) = ND
            IF (LP .NE. LPL) GO TO 3
          IF (ND .LE. 0) THEN
C
C   NODE is a boundary node.
C
            NABOR(K) = -ND
            K = K + 1
            NABOR(K) = 0
            NB = NB + 1
          ENDIF
C
C   Increment NL and print X, Y, and NABOR.
C
          INC = (K-1)/8 + 2
          NL = NL + INC
          IF (NL .GT. NLMAX) THEN
            WRITE (LUN,108)
            NL = INC
          ENDIF
          WRITE (LUN,105) NODE, X(NODE), Y(NODE),
     .                    (NABOR(I), I = 1,K)
          IF (K .NE. 8) WRITE (LUN,107)
    4     CONTINUE
      ELSE
C
C Print X, Y, Z, and LIST.
C
        WRITE (LUN,103)
        DO 6 NODE = 1,NN
          LPL = LEND(NODE)
          LP = LPL
          K = 0
C
    5     K = K + 1
            LP = LPTR(LP)
            ND = LIST(LP)
            NABOR(K) = ND
            IF (LP .NE. LPL) GO TO 5
          IF (ND .LE. 0) THEN
C
C   NODE is a boundary node.
C
            NABOR(K) = -ND
            K = K + 1
            NABOR(K) = 0
            NB = NB + 1
          ENDIF
C
C   Increment NL and print X, Y, Z, and NABOR.
C
          INC = (K-1)/5 + 2
          NL = NL + INC
          IF (NL .GT. NLMAX) THEN
            WRITE (LUN,108)
            NL = INC
          ENDIF
          WRITE (LUN,106) NODE, X(NODE), Y(NODE),
     .                    Z(NODE), (NABOR(I), I = 1,K)
          IF (K .NE. 5) WRITE (LUN,107)
    6     CONTINUE
      ENDIF
C
C Print NB, NA, and NT (boundary nodes, arcs, and
C   triangles).
C
      IF (NB .NE. 0) THEN
        NA = 3*NN - NB - 3
        NT = 2*NN - NB - 2
      ELSE
        NA = 3*NN - 6
        NT = 2*NN - 4
      ENDIF
      WRITE (LUN,109) NB, NA, NT
      RETURN
C
C Print formats:
C
  100 FORMAT (///15X,'STRIPACK Triangulation Data ',
     .        'Structure,  N = ',I5//)
  101 FORMAT (1X,'Node',31X,'Neighbors of Node'//)
  102 FORMAT (1X,'Node',5X,'Longitude',6X,'Latitude',
     .        18X,'Neighbors of Node'//)
  103 FORMAT (1X,'Node',5X,'X(Node)',8X,'Y(Node)',8X,
     .        'Z(Node)',11X,'Neighbors of Node'//)
  104 FORMAT (1X,I4,4X,14I5/(1X,8X,14I5))
  105 FORMAT (1X,I4,2E15.6,4X,8I5/(1X,38X,8I5))
  106 FORMAT (1X,I4,3E15.6,4X,5I5/(1X,53X,5I5))
  107 FORMAT (1X)
  108 FORMAT (///)
  109 FORMAT (/1X,'NB = ',I4,' Boundary Nodes',5X,
     .        'NA = ',I5,' Arcs',5X,'NT = ',I5,
     .        ' Triangles')
  110 FORMAT (1X,10X,'*** N is outside its valid',
     .        ' range ***')
      END
      SUBROUTINE VRPLOT (LUN,PLTSIZ,ELAT,ELON,A,N,X,Y,Z,
     .                   NT,LISTC,LPTR,LEND,XC,YC,ZC,TITLE,
     .                   NUMBR, IER)
      CHARACTER*(*) TITLE
      INTEGER LUN, N, NT, LISTC(*), LPTR(*), LEND(N), IER
      LOGICAL NUMBR
      DOUBLE PRECISION    PLTSIZ, ELAT, ELON, A, X(N), Y(N), Z(N),
     .        XC(NT), YC(NT), ZC(NT)
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/16/98
C
C   This subroutine creates a level-2 Encapsulated Post-
C script (EPS) file containing a graphical depiction of a
C Voronoi diagram of a set of nodes on the unit sphere.
C The visible vertices are projected onto the plane that
C contains the origin and has normal defined by a user-
C specified eye-position.  Projections of adjacent (visible)
C Voronoi vertices are connected by line segments.
C
C   The parameters defining the Voronoi diagram may be com-
C puted by Subroutine CRLIST.
C
C
C On input:
C
C       LUN = Logical unit number in the range 0 to 99.
C             The unit should be opened with an appropriate
C             file name before the call to this routine.
C
C       PLTSIZ = Plot size in inches.  A circular window in
C                the projection plane is mapped to a circu-
C                lar viewport with diameter equal to .88*
C                PLTSIZ (leaving room for labels outside the
C                viewport).  The viewport is centered on the
C                8.5 by 11 inch page, and its boundary is
C                drawn.  1.0 .LE. PLTSIZ .LE. 8.5.
C
C       ELAT,ELON = Latitude and longitude (in degrees) of
C                   the center of projection E (the center
C                   of the plot).  The projection plane is
C                   the plane that contains the origin and
C                   has E as unit normal.  In a rotated
C                   coordinate system for which E is the
C                   north pole, the projection plane con-
C                   tains the equator, and only northern
C                   hemisphere points are visible (from the
C                   point at infinity in the direction E).
C                   These are projected orthogonally onto
C                   the projection plane (by zeroing the z-
C                   component in the rotated coordinate
C                   system).  ELAT and ELON must be in the
C                   range -90 to 90 and -180 to 180, respec-
C                   tively.
C
C       A = Angular distance in degrees from E to the boun-
C           dary of a circular window against which the
C           Voronoi diagram is clipped.  The projected win-
C           dow is a disk of radius r = Sin(A) centered at
C           the origin, and only visible vertices whose
C           projections are within distance r of the origin
C           are included in the plot.  Thus, if A = 90, the
C           plot includes the entire hemisphere centered at
C           E.  0 .LT. A .LE. 90.
C
C       N = Number of nodes (Voronoi centers) and Voronoi
C           regions.  N .GE. 3.
C
C       X,Y,Z = Arrays of length N containing the Cartesian
C               coordinates of the nodes (unit vectors).
C
C       NT = Number of Voronoi region vertices (triangles,
C            including those in the extended triangulation
C            if the number of boundary nodes NB is nonzero):
C            NT = 2*N-4.
C
C       LISTC = Array of length 3*NT containing triangle
C               indexes (indexes to XC, YC, and ZC) stored
C               in 1-1 correspondence with LIST/LPTR entries
C               (or entries that would be stored in LIST for
C               the extended triangulation):  the index of
C               triangle (N1,N2,N3) is stored in LISTC(K),
C               LISTC(L), and LISTC(M), where LIST(K),
C               LIST(L), and LIST(M) are the indexes of N2
C               as a neighbor of N1, N3 as a neighbor of N2,
C               and N1 as a neighbor of N3.  The Voronoi
C               region associated with a node is defined by
C               the CCW-ordered sequence of circumcenters in
C               one-to-one correspondence with its adjacency
C               lst (in the extended triangulation).
C
C       LPTR = Array of length 3*NT = 6*N-12 containing a
C              set of pointers (LISTC indexes) in one-to-one
C              correspondence with the elements of LISTC.
C              LISTC(LPTR(I)) indexes the triangle which
C              follows LISTC(I) in cyclical counterclockwise
C              order (the first neighbor follows the last
C              neighbor).
C
C       LEND = Array of length N containing a set of
C              pointers to triangle lsts.  LP = LEND(K)
C              points to a triangle (indexed by LISTC(LP))
C              containing node K for K = 1 to N.
C
C       XC,YC,ZC = Arrays of length NT containing the
C                  Cartesian coordinates of the triangle
C                  circumcenters (Voronoi vertices).
C                  XC(I)**2 + YC(I)**2 + ZC(I)**2 = 1.
C
C       TITLE = Type CHARACTER variable or constant contain-
C               ing a string to be centered above the plot.
C               The string must be enclosed in parentheses;
C               i.e., the first and last characters must be
C               '(' and ')', respectively, but these are not
C               displayed.  TITLE may have at most 80 char-
C               acters including the parentheses.
C
C       NUMBR = Option indicator:  If NUMBR = TRUE, the
C               nodal indexes are plotted at the Voronoi
C               region centers.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if LUN, PLTSIZ, N, or NT is outside
C                     its valid range.
C             IER = 2 if ELAT, ELON, or A is outside its
C                     valid range.
C             IER = 3 if an error was encountered in writing
C                     to unit LUN.
C
C Modules required by VRPLOT:  None
C
C Intrinsic functions called by VRPLOT:  ABS, ATAN, COS,
C                                          NINT, REAL, SIN,
C                                          SQRT
C
C***********************************************************
C
      INTEGER IPX1, IPX2, IPY1, IPY2, IR, KV1, KV2, LP, LPL,
     .        N0
      LOGICAL ANNOT, IN1, IN2
      DOUBLE PRECISION    CF, CT, EX, EY, EZ, FSIZN, FSIZT, R11, R12,
     .        R21, R22, R23, SF, T, TX, TY, WR, WRS, X0, X1,
     .        X2, Y0, Y1, Y2, Z1, Z2
C
      DATA    ANNOT/.TRUE./,  FSIZN/10.0/,  FSIZT/16.0/
C
C Local parameters:
C
C ANNOT =     Logical variable with value TRUE iff the plot
C               is to be annotated with the values of ELAT,
C               ELON, and A
C CF =        Conversion factor for degrees to radians
C CT =        Cos(ELAT)
C EX,EY,EZ =  Cartesian coordinates of the eye-position E
C FSIZN =     Font size in points for labeling nodes with
C               their indexes if NUMBR = TRUE
C FSIZT =     Font size in points for the title (and
C               annotation if ANNOT = TRUE)
C IN1,IN2 =   Logical variables with value TRUE iff the
C               projections of vertices KV1 and KV2, respec-
C               tively, are inside the window
C IPX1,IPY1 = X and y coordinates (in points) of the lower
C               left corner of the bounding box or viewport
C               box
C IPX2,IPY2 = X and y coordinates (in points) of the upper
C               right corner of the bounding box or viewport
C               box
C IR =        Half the width (height) of the bounding box or
C               viewport box in points -- viewport radius
C KV1,KV2 =   Endpoint indexes of a Voronoi edge
C LP =        LIST index (pointer)
C LPL =       Pointer to the last neighbor of N0
C N0 =        Index of a node
C R11...R23 = Components of the first two rows of a rotation
C               that maps E to the north pole (0,0,1)
C SF =        Scale factor for mapping world coordinates
C               (window coordinates in [-WR,WR] X [-WR,WR])
C               to viewport coordinates in [IPX1,IPX2] X
C               [IPY1,IPY2]
C T =         Temporary variable
C TX,TY =     Translation vector for mapping world coordi-
C               nates to viewport coordinates
C WR =        Window radius r = Sin(A)
C WRS =       WR**2
C X0,Y0 =     Projection plane coordinates of node N0 or
C               label location
C X1,Y1,Z1 =  Coordinates of vertex KV1 in the rotated
C               coordinate system
C X2,Y2,Z2 =  Coordinates of vertex KV2 in the rotated
C               coordinate system or intersection of edge
C               KV1-KV2 with the equator (in the rotated
C               coordinate system)
C
C
C Test for invalid parameters.
C
      IF (LUN .LT. 0  .OR.  LUN .GT. 99  .OR.
     .    PLTSIZ .LT. 1.0  .OR.  PLTSIZ .GT. 8.5  .OR.
     .    N .LT. 3  .OR.  NT .NE. 2*N-4)
     .  GO TO 11
      IF (ABS(ELAT) .GT. 90.0  .OR.  ABS(ELON) .GT. 180.0
     .    .OR.  A .GT. 90.0) GO TO 12
C
C Compute a conversion factor CF for degrees to radians
C   and compute the window radius WR.
C
      CF = ATAN(1.0)/45.0
      WR = SIN(CF*A)
      WRS = WR*WR
C
C Compute the lower left (IPX1,IPY1) and upper right
C   (IPX2,IPY2) corner coordinates of the bounding box.
C   The coordinates, specified in default user space units
C   (points, at 72 points/inch with origin at the lower
C   left corner of the page), are chosen to preserve the
C   square aspect ratio, and to center the plot on the 8.5
C   by 11 inch page.  The center of the page is (306,396),
C   and IR = PLTSIZ/2 in points.
C
      IR = NINT(36.0*PLTSIZ)
      IPX1 = 306 - IR
      IPX2 = 306 + IR
      IPY1 = 396 - IR
      IPY2 = 396 + IR
C
C Output header comments.
C
      WRITE (LUN,100,ERR=13) IPX1, IPY1, IPX2, IPY2
  100 FORMAT ('%!PS-Adobe-3.0 EPSF-3.0'/
     .        '%%BoundingBox:',4I4/
     .        '%%Title:  Voronoi diagram'/
     .        '%%Creator:  STRIPACK'/
     .        '%%EndComments')
C
C Set (IPX1,IPY1) and (IPX2,IPY2) to the corner coordinates
C   of a viewport box obtained by shrinking the bounding box
C   by 12% in each dimension.
C
      IR = NINT(0.88*REAL(IR))
      IPX1 = 306 - IR
      IPX2 = 306 + IR
      IPY1 = 396 - IR
      IPY2 = 396 + IR
C
C Set the line thickness to 2 points, and draw the
C   viewport boundary.
C
      T = 2.0
      WRITE (LUN,110,ERR=13) T
      WRITE (LUN,120,ERR=13) IR
      WRITE (LUN,130,ERR=13)
  110 FORMAT (F12.6,' setlinewidth')
  120 FORMAT ('306 396 ',I3,' 0 360 arc')
  130 FORMAT ('stroke')
C
C Set up an affine mapping from the window box [-WR,WR] X
C   [-WR,WR] to the viewport box.
C
      SF = REAL(IR)/WR
      TX = IPX1 + SF*WR
      TY = IPY1 + SF*WR
      WRITE (LUN,140,ERR=13) TX, TY, SF, SF
  140 FORMAT (2F12.6,' translate'/
     .        2F12.6,' scale')
C
C The line thickness must be changed to reflect the new
C   scaling which is applied to all subsequent output.
C   Set it to 1.0 point.
C
      T = 1.0/SF
      WRITE (LUN,110,ERR=13) T
C
C Save the current graphics state, and set the clip path to
C   the boundary of the window.
C
      WRITE (LUN,150,ERR=13)
      WRITE (LUN,160,ERR=13) WR
      WRITE (LUN,170,ERR=13)
  150 FORMAT ('gsave')
  160 FORMAT ('0 0 ',F12.6,' 0 360 arc')
  170 FORMAT ('clip newpath')
C
C Compute the Cartesian coordinates of E and the components
C   of a rotation R which maps E to the north pole (0,0,1).
C   R is taken to be a rotation about the z-axis (into the
C   yz-plane) followed by a rotation about the x-axis chosen
C   so that the view-up direction is (0,0,1), or (-1,0,0) if
C   E is the north or south pole.
C
C           ( R11  R12  0   )
C       R = ( R21  R22  R23 )
C           ( EX   EY   EZ  )
C
      T = CF*ELON
      CT = COS(CF*ELAT)
      EX = CT*COS(T)
      EY = CT*SIN(T)
      EZ = SIN(CF*ELAT)
      IF (CT .NE. 0.0) THEN
        R11 = -EY/CT
        R12 = EX/CT
      ELSE
        R11 = 0.0
        R12 = 1.0
      ENDIF
      R21 = -EZ*R12
      R22 = EZ*R11
      R23 = CT
C
C Loop on nodes (Voronoi centers) N0.
C   LPL indexes the last neighbor of N0.
C
      DO 3 N0 = 1,N
        LPL = LEND(N0)
C
C Set KV2 to the first (and last) vertex index and compute
C   its coordinates (X2,Y2,Z2) in the rotated coordinate
C   system.
C
        KV2 = LISTC(LPL)
        X2 = R11*XC(KV2) + R12*YC(KV2)
        Y2 = R21*XC(KV2) + R22*YC(KV2) + R23*ZC(KV2)
        Z2 = EX*XC(KV2) + EY*YC(KV2) + EZ*ZC(KV2)
C
C   IN2 = TRUE iff KV2 is in the window.
C
        IN2 = Z2 .GE. 0.  .AND.  X2*X2 + Y2*Y2 .LE. WRS
C
C Loop on neighbors N1 of N0.  For each triangulation edge
C   N0-N1, KV1-KV2 is the corresponding Voronoi edge.
C
        LP = LPL
    1   LP = LPTR(LP)
          KV1 = KV2
          X1 = X2
          Y1 = Y2
          Z1 = Z2
          IN1 = IN2
          KV2 = LISTC(LP)
C
C   Compute the new values of (X2,Y2,Z2) and IN2.
C
          X2 = R11*XC(KV2) + R12*YC(KV2)
          Y2 = R21*XC(KV2) + R22*YC(KV2) + R23*ZC(KV2)
          Z2 = EX*XC(KV2) + EY*YC(KV2) + EZ*ZC(KV2)
          IN2 = Z2 .GE. 0.  .AND.  X2*X2 + Y2*Y2 .LE. WRS
C
C Add edge KV1-KV2 to the path iff both endpoints are inside
C   the window and KV2 > KV1, or KV1 is inside and KV2 is
C   outside (so that the edge is drawn only once).
C
          IF (.NOT. IN1  .OR.  (IN2  .AND.  KV2 .LE. KV1))
     .      GO TO 2
          IF (Z2 .LT. 0.) THEN
C
C   KV2 is a 'southern hemisphere' point.  Move it to the
C     intersection of edge KV1-KV2 with the equator so that
C     the edge is clipped properly.  Z2 is implicitly set
C     to 0.
C
            X2 = Z1*X2 - Z2*X1
            Y2 = Z1*Y2 - Z2*Y1
            T = SQRT(X2*X2+Y2*Y2)
            X2 = X2/T
            Y2 = Y2/T
          ENDIF
          WRITE (LUN,180,ERR=13) X1, Y1, X2, Y2
  180     FORMAT (2F12.6,' moveto',2F12.6,' lineto')
C
C Bottom of loops.
C
    2     IF (LP .NE. LPL) GO TO 1
    3   CONTINUE
C
C Paint the path and restore the saved graphics state (with
C   no clip path).
C
      WRITE (LUN,130,ERR=13)
      WRITE (LUN,190,ERR=13)
  190 FORMAT ('grestore')
      IF (NUMBR) THEN
C
C Nodes in the window are to be labeled with their indexes.
C   Convert FSIZN from points to world coordinates, and
C   output the commands to select a font and scale it.
C
        T = FSIZN/SF
        WRITE (LUN,200,ERR=13) T
  200   FORMAT ('/Helvetica findfont'/
     .          F12.6,' scalefont setfont')
C
C Loop on visible nodes N0 that project to points (X0,Y0) in
C   the window.
C
        DO 4 N0 = 1,N
          IF (EX*X(N0) + EY*Y(N0) + EZ*Z(N0) .LT. 0.)
     .      GO TO 4
          X0 = R11*X(N0) + R12*Y(N0)
          Y0 = R21*X(N0) + R22*Y(N0) + R23*Z(N0)
          IF (X0*X0 + Y0*Y0 .GT. WRS) GO TO 4
C
C   Move to (X0,Y0), and draw the label N0 with the origin
C     of the first character at (X0,Y0).
C
          WRITE (LUN,210,ERR=13) X0, Y0
          WRITE (LUN,220,ERR=13) N0
  210     FORMAT (2F12.6,' moveto')
  220     FORMAT ('(',I3,') show')
    4     CONTINUE
      ENDIF
C
C Convert FSIZT from points to world coordinates, and output
C   the commands to select a font and scale it.
C
      T = FSIZT/SF
      WRITE (LUN,200,ERR=13) T
C
C Display TITLE centered above the plot:
C
      Y0 = WR + 3.0*T
      WRITE (LUN,230,ERR=13) TITLE, Y0
  230 FORMAT (A80/'  stringwidth pop 2 div neg ',F12.6,
     .        ' moveto')
      WRITE (LUN,240,ERR=13) TITLE
  240 FORMAT (A80/'  show')
      IF (ANNOT) THEN
C
C Display the window center and radius below the plot.
C
        X0 = -WR
        Y0 = -WR - 50.0/SF
        WRITE (LUN,210,ERR=13) X0, Y0
        WRITE (LUN,250,ERR=13) ELAT, ELON
        Y0 = Y0 - 2.0*T
        WRITE (LUN,210,ERR=13) X0, Y0
        WRITE (LUN,260,ERR=13) A
  250   FORMAT ('(Window center:  ELAT = ',F7.2,
     .          ',  ELON = ',F8.2,') show')
  260   FORMAT ('(Angular extent:  A = ',F5.2,') show')
      ENDIF
C
C Paint the path and output the showpage command and
C   end-of-file indicator.
C
      WRITE (LUN,270,ERR=13)
  270 FORMAT ('stroke'/
     .        'showpage'/
     .        '%%EOF')
C
C HP's interpreters require a one-byte End-of-PostScript-Job
C   indicator (to eliminate a timeout error message):
C   ASCII 4.
C
      WRITE (LUN,280,ERR=13) CHAR(4)
  280 FORMAT (A1)
C
C No error encountered.
C
      IER = 0
      RETURN
C
C Invalid input parameter LUN, PLTSIZ, N, or NT.
C
   11 IER = 1
      RETURN
C
C Invalid input parameter ELAT, ELON, or A.
C
   12 IER = 2
      RETURN
C
C Error writing to unit LUN.
C
   13 IER = 3
      RETURN
      END
      SUBROUTINE APLYR (X,Y,Z,CX,SX,CY,SY, XP,YP,ZP)
      DOUBLE PRECISION X, Y, Z, CX, SX, CY, SY, XP, YP, ZP
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   05/09/92
C
C   This subroutine applies the rotation R defined by Sub-
C routine CONSTR to the unit vector (X Y Z)**T, i,e. (X,Y,Z)
C is rotated to (XP,YP,ZP).  If (XP,YP,ZP) lies in the
C southern hemisphere (ZP < 0), (XP,YP) are set to the
C coordinates of the nearest point of the equator, ZP re-
C maining unchanged.
C
C On input:
C
C       X,Y,Z = Coordinates of a point on the unit sphere.
C
C       CX,SX,CY,SY = Elements of the rotation defined by
C                     Subroutine CONSTR.
C
C Input parameters are not altered except as noted below.
C
C On output:
C
C       XP,YP,ZP = Coordinates of the rotated point on the
C                  sphere unless ZP < 0, in which case
C                  (XP,YP,0) is the closest point of the
C                  equator to the rotated point.  Storage
C                  for XP, YP, and ZP may coincide with
C                  storage for X, Y, and Z, respectively,
C                  if the latter need not be saved.
C
C Modules required by APLYR:  None
C
C Intrinsic function called by APLYR:  SQRT
C
C***********************************************************
C
      DOUBLE PRECISION T
C
C Local parameter:
C
C T = Temporary variable
C
      T = SX*Y + CX*Z
      YP = CX*Y - SX*Z
      ZP = SY*X + CY*T
      XP = CY*X - SY*T
      IF (ZP .GE. 0.) RETURN
C
C Move (XP,YP,ZP) to the equator.
C
      T = SQRT(XP*XP + YP*YP)
      IF (T .EQ. 0.) GO TO 1
      XP = XP/T
      YP = YP/T
      RETURN
C
C Move the south pole to an arbitrary point of the equator.
C
    1 XP = 1.
      YP = 0.
      RETURN
      END
      SUBROUTINE APLYRT (G1P,G2P,CX,SX,CY,SY, G)
      DOUBLE PRECISION G1P, G2P, CX, SX, CY, SY, G(3)
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   05/09/92
C
C   This subroutine applies the inverse (transpose) of the
C rotation defined by Subroutine CONSTR to the vector
C (G1P G2P 0)**T, i.e., the gradient (G1P,G2P,0) in the rot-
C ated coordinate system is mapped to (G1,G2,G3) in the
C original coordinate system.
C
C On input:
C
C       G1P,G2P = X and Y components, respectively, of the
C                 gradient in the rotated coordinate system.
C
C       CX,SX,CY,SY = Elements of the rotation R constructed
C                     by Subroutine CONSTR.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       G = X, Y, and Z components (in that order) of the
C           inverse rotation applied to (G1P,G2P,0) --
C           gradient in the original coordinate system.
C
C Modules required by APLYRT:  None
C
C***********************************************************
C
      DOUBLE PRECISION T
C
C Local parameters:
C
C T = Temporary variable
C
      T = SY*G1P
      G(1) = CY*G1P
      G(2) = CX*G2P - SX*T
      G(3) = -SX*G2P - CX*T
      RETURN
      END
      SUBROUTINE ARCINT (P,P1,P2,F1,F2,G1,G2,SIGMA, F,G,GN)
      DOUBLE PRECISION    P(3), P1(3), P2(3), F1, F2, G1(3), G2(3),
     .        SIGMA, F, G(3), GN
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/21/96
C
C   Given 3 points P, P1, and P2 lying on a common geodesic
C of the unit sphere with P between P1 and P2, along with
C data values and gradients at P1 and P2, this subroutine
C computes an interpolated value F and a gradient vector G
C AT P.  F and the tangential component of G are taken to be
C the value and derivative (with respect to arc-length) of
C a Hermite interpolatory tension spline defined by the end-
C point values and tangential gradient components.  The nor-
C mal component of G is obtained by linear interpolation of
C the normal components of the gradients at P1 and P2.
C
C On input:
C
C       P = Cartesian coordinates of a point lying on the
C           arc defined by P1 and P2.  P(1)**2 + P(2)**2 +
C           P(3)**2 = 1.
C
C       P1,P2 = Coordinates of distinct points on the unit
C               sphere defining an arc with length less than
C               180 degrees.
C
C       F1,F2 = Data values associated with P1 and P2,
C               respectively.
C
C       G1,G2 = Gradient vectors associated with P1 and P2.
C               G1 and G2 are orthogonal to P1 and P2,
C               respectively.
C
C       SIGMA = Tension factor associated with P1-P2.
C
C The above parameters are not altered by this routine.
C
C       G = Array of length 3.
C
C On output:
C
C       F = Interpolated value at P.
C
C       G = Interpolated gradient at P.
C
C       GN = Normal component of G with the direction
C            P1 X P2 taken to be positive.  The extrapola-
C            tion procedure requires this component.
C
C   For each vector V, V(1), V(2), and V(3) contain X, Y,
C and Z components, respectively.
C
C SSRFPACK modules required by ARCINT:  ARCLEN, SNHCSH
C
C Intrinsic functions called by ARCINT:  ABS, EXP, SQRT
C
C***********************************************************
C
      DOUBLE PRECISION    ARCLEN
      INTEGER I, LUN
      DOUBLE PRECISION    A, AL, B1, B2, CM, CMM, CM2, DUMMY, D1, D2, E,
     .        EMS, E1, E2, GT, S, SB1, SB2, SIG, SINH,
     .        SINH2, SM, SM2, TAU1, TAU2, TM, TM1, TM2, TP1,
     .        TP2, TS, UN(3), UNORM
      DATA    LUN/6/
C
C Local parameters:
C
C A =         Angle in radians (arc-length) between P1 and
C               P2
C AL =        Arc-length between P1 and P
C B1,B2 =     Local coordinates of P with respect to P1-P2
C CM,CMM =    Coshm(SIG) and Coshmm(SIG) -- refer to SNHCSH
C CM2 =       Coshm(SB2)
C DUMMY =     Dummy parameter for SNHCSH
C D1,D2 =     Scaled second differences
C E =         CM**2 - SM*Sinh = SIG*SM - 2*CMM (scaled by
C               2*EMS if SIG > .5)
C EMS =       Exp(-SIG)
C E1,E2 =     Exp(-SB1), Exp(-SB2)
C GT =        Tangential component of G -- component in the
C               direction UN X P
C I =         DO-loop index
C LUN =       Logical unit for error messages
C S =         Slope:  (F2-F1)/A
C SB1,SB2 =   SIG*B1, SIG*B2
C SIG =       Abs(SIGMA)
C SINH =      Sinh(SIGMA)
C SINH2 =     Sinh(SB2)
C SM,SM2 =    Sinhm(SIG), Sinhm(SB2)
C TAU1,TAU2 = Tangential derivatives (components of G1,G2)
C               at P1 and P2
C TM =        1-EMS
C TM1,TM2 =   1-E1, 1-E2
C TP1,TP2 =   1+E1, 1+E2
C TS =        TM**2
C UN =        Unit normal to the plane of P, P1, and P2
C UNORM =     Euclidean norm of P1 X P2 -- used to normalize
C               UN
C
C
C Compute unit normal UN.
C
      UN(1) = P1(2)*P2(3) - P1(3)*P2(2)
      UN(2) = P1(3)*P2(1) - P1(1)*P2(3)
      UN(3) = P1(1)*P2(2) - P1(2)*P2(1)
      UNORM = SQRT(UN(1)*UN(1) + UN(2)*UN(2) + UN(3)*UN(3))
      IF (UNORM .EQ. 0.) GO TO 2
C
C Normalize UN.
C
      DO 1 I = 1,3
        UN(I) = UN(I)/UNORM
    1   CONTINUE
C
C Compute tangential derivatives at the endpoints:
C   TAU1 = (G1,UN X P1) = (G1,P2)/UNORM and
C   TAU2 = (G2,UN X P2) = -(G2,P1)/UNORM.
C
      TAU1 = (G1(1)*P2(1) + G1(2)*P2(2) + G1(3)*P2(3))/UNORM
      TAU2 =-(G2(1)*P1(1) + G2(2)*P1(2) + G2(3)*P1(3))/UNORM
C
C Compute arc-lengths A, AL.
C
      A = ARCLEN(P1,P2)
      IF (A .EQ. 0.) GO TO 2
      AL = ARCLEN(P1,P)
C
C Compute local coordinates, slope, and second differences.
C
      B2 = AL/A
      B1 = 1. - B2
      S = (F2-F1)/A
      D1 = S - TAU1
      D2 = TAU2 - S
C
C Test the range of SIGMA.
C
      SIG = ABS(SIGMA)
      IF (SIG .LT. 1.E-9) THEN
C
C Hermite cubic interpolation.
C
        F = F1 + AL*(TAU1 + B2*(D1 + B1*(D1 - D2)))
        GT = TAU1 + B2*(D1 + D2 + 3.*B1*(D1 - D2))
      ELSEIF (SIG .LE. .5) THEN
C
C 0 < SIG .LE. .5.  Use approximations designed to avoid
C   cancellation error in the hyperbolic functions.
C
        SB2 = SIG*B2
        CALL SNHCSH (SIG, SM,CM,CMM)
        CALL SNHCSH (SB2, SM2,CM2,DUMMY)
        SINH = SM + SIG
        SINH2 = SM2 + SB2
        E = SIG*SM - CMM - CMM
        F = F1 + AL*TAU1 + A*((CM*SM2-SM*CM2)*(D1+D2) + SIG*
     .                        (CM*CM2-SINH*SM2)*D1)/(SIG*E)
        GT = TAU1 + ((CM*CM2-SM*SINH2)*(D1+D2) + SIG*
     .               (CM*SINH2-SINH*CM2)*D1)/E
      ELSE
C
C SIG > .5.  Use negative exponentials in order to avoid
C   overflow.  Note that EMS = EXP(-SIG).
C
        SB1 = SIG*B1
        SB2 = SIG - SB1
        E1 = EXP(-SB1)
        E2 = EXP(-SB2)
        EMS = E1*E2
        TM = 1. - EMS
        TS = TM*TM
        TM1 = 1. - E1
        TM2 = 1. - E2
        E = TM*(SIG*(1.+EMS) - TM - TM)
        F = F1 + AL*S + A*(TM*TM1*TM2*(D1+D2) + SIG*
     .                     ((E2*TM1*TM1-B1*TS)*D1 +
     .                      (E1*TM2*TM2-B2*TS)*D2))/(SIG*E)
        TP1 = 1. + E1
        TP2 = 1. + E2
        GT = S + (TM1*(TM*TP2-SIG*E2*TP1)*D1 -
     .            TM2*(TM*TP1-SIG*E1*TP2)*D2)/E
      ENDIF
C
C Compute GN.
C
      GN = B1*(UN(1)*G1(1) + UN(2)*G1(2) + UN(3)*G1(3)) +
     .     B2*(UN(1)*G2(1) + UN(2)*G2(2) + UN(3)*G2(3))
C
C Compute G = GT*(UN X P) + GN*UN.
C
      G(1) = GT*(UN(2)*P(3) - UN(3)*P(2)) + GN*UN(1)
      G(2) = GT*(UN(3)*P(1) - UN(1)*P(3)) + GN*UN(2)
      G(3) = GT*(UN(1)*P(2) - UN(2)*P(1)) + GN*UN(3)
      RETURN
C
C P1 X P2 = 0.  Print an error message and terminate
C   processing.
C
    2 WRITE (LUN,100) (P1(I),I=1,3), (P2(I),I=1,3)
  100 FORMAT ('1','ERROR IN ARCINT -- P1 = ',2(F9.6,',  '),
     .        F9.6/1X,19X,'P2 = ',2(F9.6,',  '),F9.6)
      STOP
      END
      DOUBLE PRECISION FUNCTION ARCLEN (P,Q)
      DOUBLE PRECISION P(3), Q(3)
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   05/09/92
C
C   This function computes the arc-length (angle in radians)
C between a pair of points on the unit sphere.
C
C On input:
C
C       P,Q = Arrays of length 3 containing the X, Y, and Z
C             coordinates (in that order) of points on the
C             unit sphere.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       ARCLEN = Angle in radians between the unit vectors
C                P and Q.  0 .LE. ARCLEN .LE. PI.
C
C Modules required by ARCLEN:  None
C
C Intrinsic functions called by ARCLEN:  ATAN, SQRT
C
C***********************************************************
C
      INTEGER I
      DOUBLE PRECISION    D
C
C Local parameters:
C
C D = Euclidean norm squared of P+Q
C I = DO-loop index
C
      D = 0.
      DO 1 I = 1,3
        D = D + (P(I) + Q(I))**2
    1   CONTINUE
      IF (D .EQ. 0.) THEN
C
C P and Q are separated by 180 degrees.
C
        ARCLEN = 4.*ATAN(1.)
      ELSEIF (D .GE. 4.) THEN
C
C P and Q coincide.
C
        ARCLEN = 0.
      ELSE
        ARCLEN = 2.*ATAN(SQRT((4.-D)/D))
      ENDIF
      RETURN
      END
      SUBROUTINE CONSTR (XK,YK,ZK, CX,SX,CY,SY)
      DOUBLE PRECISION XK, YK, ZK, CX, SX, CY, SY
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   05/09/92
C
C   This subroutine constructs the elements of a 3 by 3
C orthogonal matrix R which rotates a point (XK,YK,ZK) on
C the unit sphere to the north pole, i.e.,
C
C      (XK)     (CY  0 -SY)   (1   0   0)   (XK)     (0)
C  R * (YK)  =  ( 0  1   0) * (0  CX -SX) * (YK)  =  (0)
C      (ZK)     (SY  0  CY)   (0  SX  CX)   (ZK)     (1)
C
C On input:
C
C       XK,YK,ZK = Components of a unit vector to be
C                  rotated to (0,0,1).
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       CX,SX,CY,SY = Elements of R:  CX,SX define a rota-
C                     tion about the X-axis and CY,SY define
C                     a rotation about the Y-axis.
C
C Modules required by CONSTR:  None
C
C Intrinsic function called by CONSTR:  SQRT
C
C***********************************************************
C
      CY = SQRT(YK*YK + ZK*ZK)
      SY = XK
      IF (CY .NE. 0.) THEN
        CX = ZK/CY
        SX = YK/CY
      ELSE
C
C (XK,YK,ZK) lies on the X-axis.
C
        CX = 1.
        SX = 0.
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION FVAL (B1,B2,B3,V1,V2,V3,F1,F2,F3,G1,G2,
     .                    G3,SIG1,SIG2,SIG3)
      DOUBLE PRECISION B1, B2, B3, V1(3), V2(3), V3(3), F1, F2, F3,
     .     G1(3), G2(3), G3(3), SIG1, SIG2, SIG3
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   05/09/92
C
C   Given data values and gradients at the three vertices of
C a spherical triangle containing a point P, this routine
C computes the value of F at P where F interpolates the ver-
C tex data.  Along the triangle sides, the interpolatory
C function F is the Hermite interpolatory tension spline
C defined by the values and tangential gradient components
C at the endpoints, and the gradient component normal to the
C triangle side varies linearly with respect to arc-length
C between the normal gradient components at the endpoints.
C A first-order C-1 blending method is used on the underly-
C ing planar triangle.  Since values and gradients on an arc
C depend only on the vertex data, the method results in C-1
C continuity when used to interpolate over a triangulation.
C
C   The blending method consists of taking F(P) to be a
C weighted sum of the values at PP of three univariate Her-
C mite interpolatory tension splines defined on the line
C segments which join the vertices to the opposite sides and
C pass through PP:  the central projection of P onto the
C underlying planar triangle.  The tension factors for these
C splines are obtained by linear interpolation between the
C pair of tension factors associated with the triangle sides
C which join at the appropriate vertex.
C
C   A tension factor SIGMA associated with a Hermite interp-
C olatory tension spline is a nonnegative parameter which
C determines the curviness of the spline.  SIGMA = 0 results
C in a cubic spline, and the spline approaches the linear
C interpolant as SIGMA increases.
C
C On input:
C
C       B1,B2,B3 = Barycentric coordinates of PP with re-
C                  spect to the (planar) underlying triangle
C                  (V1,V2,V3), where PP is the central
C                  projection of P onto this triangle.
C
C       V1,V2,V3 = Cartesian coordinates of the vertices of
C                  a spherical triangle containing P.  V3
C                  Left V1->V2.
C
C       F1,F2,F3 = Data values associated with the vertices.
C
C       G1,G2,G3 = Gradients associated with the vertices.
C                  Gi is orthogonal to Vi for i = 1,2,3.
C
C       SIG1,SIG2,SIG3 = Tension factors associated with the
C                        triangle sides opposite V1, V2, and
C                        V3, respectively.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       FVAL = Interpolated value at P.
C
C Each vector V above contains X, Y, and Z components in
C   V(1), V(2), and V(3), respectively.
C
C SSRFPACK modules required by FVAL:  ARCINT, ARCLEN, HVAL
C
C Intrinsic function called by FVAL:  SQRT
C
C***********************************************************
C
      DOUBLE PRECISION    HVAL
      INTEGER I
      DOUBLE PRECISION    C1, C2, C3, DS, DUM, DV, F, G(3),
     .        Q1(3), Q2(3), Q3(3), SIG, SUM, S1, S2, S3,
     .        U1(3), U2(3), U3(3), U1N, U2N, U3N, VAL
C
C Local parameters:
C
C C1,C2,C3 =    Coefficients (weight functions) of partial
C                 interpolants.  C1 = 1 on the edge opposite
C                 V1 and C1 = 0 on the other edges.  Simi-
C                 larly for C2 and C3.  C1+C2+C3 = 1.
C DS =          Directional derivative (scaled by distnace)
C                 at U1, U2, or U3:  DS = (G,U1-V1)/U1N =
C                 -(G,V1)/U1N on side opposite V1, where G/
C                 U1N (plus an orthogonal component) is the
C                 projection of G onto the planar triangle
C DUM =         Dummy variable for calls to ARCINT
C DV =          Directional derivatives (scaled by distance)
C                 at a vertex:  D1 = (G1,U1-V1) = (G1,U1)
C F,G =         Value and gradient at Q1 Q2, or Q3 obtained
C                 by interpolation along one of the arcs of
C                 the spherical triangle
C I =           DO-loop index
C Q1,Q2,Q3 =    Central projections of U1, U2, and U3 onto
C                 the sphere and thus lying on an arc of the
C                 spherical triangle
C SIG =         Tension factor for a side-vertex (partial)
C                 interpolant:  obtained by linear interpo-
C                 lation applied to triangle side tensions
C SUM =         Quantity used to normalize C1, C2, and C3
C S1,S2,S3 =    Sums of pairs of barycentric coordinates:
C                 used to compute U1, U2, U3, and SIG
C U1,U2,U3 =    Points on the boundary of the planar trian-
C                 gle and lying on the lines containing PP
C                 and the vertices.  U1 is opposite V1, etc.
C U1N,U2N,U3N = Quantities used to compute Q1, Q2, and Q3
C                 (magnitudes of U1, U2, and U3)
C VAL =         Local variable used to accumulate the con-
C                 tributions to FVAL
C
C
C Compute weight functions C1, C2, and C3.
C
      C1 = B2*B3
      C2 = B3*B1
      C3 = B1*B2
      SUM = C1 + C2 + C3
      IF (SUM .LE. 0.) THEN
C
C P coincides with a vertex.
C
        FVAL = B1*F1 + B2*F2 + B3*F3
        RETURN
      ENDIF
C
C Normalize C1, C2, and C3.
C
      C1 = C1/SUM
      C2 = C2/SUM
      C3 = C3/SUM
C
C Compute (S1,S2,S3), (U1,U2,U3) and (U1N,U2N,U3N).
C
      S1 = B2 + B3
      S2 = B3 + B1
      S3 = B1 + B2
      U1N = 0.
      U2N = 0.
      U3N = 0.
      DO 1 I = 1,3
        U1(I) = (B2*V2(I) + B3*V3(I))/S1
        U2(I) = (B3*V3(I) + B1*V1(I))/S2
        U3(I) = (B1*V1(I) + B2*V2(I))/S3
        U1N = U1N + U1(I)*U1(I)
        U2N = U2N + U2(I)*U2(I)
        U3N = U3N + U3(I)*U3(I)
    1   CONTINUE
C
C Compute Q1, Q2, and Q3.
C
      U1N = SQRT(U1N)
      U2N = SQRT(U2N)
      U3N = SQRT(U3N)
      DO 2 I = 1,3
        Q1(I) = U1(I)/U1N
        Q2(I) = U2(I)/U2N
        Q3(I) = U3(I)/U3N
    2   CONTINUE
C
C Compute interpolated value (VAL) at P by looping on
C   triangle sides.
C
      VAL = 0.
C
C Contribution from side opposite V1:
C
C   Compute value and gradient at Q1 by interpolating
C     between V2 and V3.
C
      CALL ARCINT (Q1,V2,V3,F2,F3,G2,G3,SIG1, F,G,DUM)
C
C   Add in the contribution.
C
      DV = G1(1)*U1(1) + G1(2)*U1(2) + G1(3)*U1(3)
      DS = -(G(1)*V1(1) + G(2)*V1(2) + G(3)*V1(3))/U1N
      SIG = (B2*SIG3 + B3*SIG2)/S1
      VAL = VAL + C1*HVAL(B1,F1,F,DV,DS,SIG)
C
C Contribution from side opposite V2:
C
C   Compute value and gradient at Q2 by interpolating
C     between V3 and V1.
C
      CALL ARCINT (Q2,V3,V1,F3,F1,G3,G1,SIG2, F,G,DUM)
C
C   Add in the contribution.
C
      DV = G2(1)*U2(1) + G2(2)*U2(2) + G2(3)*U2(3)
      DS = -(G(1)*V2(1) + G(2)*V2(2) + G(3)*V2(3))/U2N
      SIG = (B3*SIG1 + B1*SIG3)/S2
      VAL = VAL + C2*HVAL(B2,F2,F,DV,DS,SIG)
C
C Contribution from side opposite V3:
C
C   Compute interpolated value and gradient at Q3
C     by interpolating between V1 and V2.
C
      CALL ARCINT (Q3,V1,V2,F1,F2,G1,G2,SIG3, F,G,DUM)
C
C   Add in the final contribution.
C
      DV = G3(1)*U3(1) + G3(2)*U3(2) + G3(3)*U3(3)
      DS = -(G(1)*V3(1) + G(2)*V3(2) + G(3)*V3(3))/U3N
      SIG = (B1*SIG2 + B2*SIG1)/S3
      FVAL = VAL + C3*HVAL(B3,F3,F,DV,DS,SIG)
      RETURN
      END
      SUBROUTINE GETSIG (N,X,Y,Z,H,LIST,LPTR,LEND,GRAD,
     .                   TOL, SIGMA, DSMAX,IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), IER
      DOUBLE PRECISION    X(N), Y(N), Z(N), H(N), GRAD(3,N), TOL,
     .        SIGMA(*), DSMAX
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/21/96
C
C   Given a triangulation of a set of nodes on the unit
C sphere, along with data values H and gradients GRAD at the
C nodes, this subroutine determines, for each triangulation
C arc, the smallest (nonnegative) tension factor SIGMA such
C that the Hermite interpolatory tension spline H(A), de-
C fined by SIGMA and the endpoint values and directional
C derivatives, preserves local shape properties of the data.
C In order to define the shape properties on an arc, it is
C convenient to map the arc to an interval (A1,A2).  Then,
C denoting the endpoint data values by H1,H2 and the deriva-
C tives (tangential gradient components) by HP1,HP2, and
C letting S = (H2-H1)/(A2-A1), the data properties are
C
C       Monotonicity:  S, HP1, and HP2 are nonnegative or
C                        nonpositive,
C   and
C
C       Convexity:     HP1 .LE. S .LE. HP2  or  HP1 .GE. S
C                        .GE. HP2.
C
C The corresponding properties of H are constant sign of the
C first and second derivatives, respectively.  Note that,
C unless HP1 = S = HP2, infinite tension is required (and H
C is linear on the interval) if S = 0 in the case of mono-
C tonicity, or if HP1 = S or HP2 = S in the case of
C convexity.
C
C   Note that if gradients are to be computed by Subroutine
C GRADG or function values and gradients are computed by
C SMSURF, it may be desirable to alternate those computa-
C tions (which require tension factors) with calls to this
C subroutine.  This iterative procedure should terminate
C with a call to GETSIG in order to ensure that the shape
C properties are preserved, and convergence can be achieved
C (at the cost of optimality) by allowing only increases in
C tension factors (refer to the parameter descriptions for
C SIGMA, DSMAX, and IER).
C
C   Refer to functions SIG0, SIG1, and SIG2 for means of
C selecting minimum tension factors to preserve more general
C properties.
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y,Z = Arrays of length N containing the Cartesian
C               coordinates of the nodes.
C
C       H = Array of length N containing data values at the
C           nodes.  H(I) is associated with (X(I),Y(I),Z(I))
C           for I = 1,...,N.
C
C       LIST,LPTR,LEND = Data structure defining the tri-
C                        angulation.  Refer to STRIPACK
C                        Subroutine TRMESH.
C
C       GRAD = Array dimensioned 3 by N whose columns con-
C              tain gradients at the nodes.  GRAD( ,J) must
C              be orthogonal to node J:  GRAD(1,J)*X(J) +
C              GRAD(2,J)*Y(J) + GRAD(3,J)*Z(J) = 0..  Refer
C              to Subroutines GRADG, GRADL, and SMSURF.
C
C       TOL = Tolerance whose magnitude determines how close
C             each tension factor is to its optimal value
C             when nonzero finite tension is necessary and
C             sufficient to satisfy the constraint --
C             abs(TOL) is an upper bound on the magnitude
C             of the smallest (nonnegative) or largest (non-
C             positive) value of the first or second deriva-
C             tive of H in the interval.  Thus, the con-
C             straint is satisfied, but possibly with more
C             tension than necessary.
C
C The above parameters are not altered by this routine.
C
C       SIGMA = Array of length 2*NA = 6*(N-1)-2*NB, where
C               NA and NB are the numbers of arcs and boun-
C               dary nodes, respectively, containing minimum
C               values of the tension factors.  The tension
C               factors are associated with arcs in one-to-
C               one correspondence with LIST entries.  Note
C               that each arc N1-N2 has two LIST entries and
C               thus, the tension factor is stored in both
C               SIGMA(I) and SIGMA(J) where LIST(I) = N2 (in
C               the adjacency lst for N1) and LIST(J) = N1
C               (in the lst associated with N2).  SIGMA
C               should be set to all zeros if minimal ten-
C               sion is desired, and should be unchanged
C               from a previous call in order to ensure con-
C               vergence of the iterative procedure describ-
C               ed in the header comments.
C
C On output:
C
C       SIGMA = Array containing tension factors for which
C               H(A) preserves the local data properties on
C               each triangulation arc, with the restriction
C               that SIGMA(I) .LE. 85 for all I (unless the
C               input value is larger).  The factors are as
C               small as possible (within the tolerance) but
C               not less than their input values.  If infin-
C               ite tension is required on an arc, the cor-
C               responding factor is SIGMA(I) = 85 (and H
C               is an approximation to the linear inter-
C               polant on the arc), and if neither property
C               is satisfied by the data, then SIGMA(I) = 0
C               (assuming its input value is 0), and thus H
C               is cubic on the arc.
C
C       DSMAX = Maximum increase in a component of SIGMA
C               from its input value.
C
C       IER = Error indicator and information flag:
C             IER = I if no errors were encountered and I
C                     components of SIGMA were altered from
C                     their input values for I .GE. 0.
C             IER = -1 if N < 3.  SIGMA is not altered in
C                      this case.
C             IER = -2 if duplicate nodes were encountered.
C
C STRIPACK modules required by GETSIG:  LSTPTR, STORE
C
C SSRFPACK modules required by GETSIG:  ARCLEN, SNHCSH
C
C Intrinsic functions called by GETSIG:  ABS, EXP, MAX, MIN,
C                                          SIGN, SQRT
C
C***********************************************************
C
      INTEGER LSTPTR
      DOUBLE PRECISION    ARCLEN, STORE
      INTEGER ICNT, LP1, LP2, LPL, LUN, N1, N2, NIT, NM1
      DOUBLE PRECISION    A, AL, C1, C2, COSHM, COSHMM, D0, D1, D1D2,
     .        D1PD2, D2, DMAX, DSIG, DSM, E, EMS, EMS2, F,
     .        F0, FMAX, FNEG, FP, FTOL, P1(3), P2(3), RTOL,
     .        S, S1, S2, SBIG, SCM, SGN, SIG, SIGIN, SINHM,
     .        SSINH, SSM, STOL, T, T0, T1, T2, TM, TP1,
     .        UN(3), UNORM
C
      DATA SBIG/85./,  LUN/-1/
      NM1 = N - 1
      IF (NM1 .LT. 2) GO TO 11
C
C Compute an absolute tolerance FTOL = abs(TOL) and a
C   relative tolerance RTOL = 100*Macheps.
C
      FTOL = ABS(TOL)
      RTOL = 1.
    1 RTOL = RTOL/2.
        IF (STORE(RTOL+1.) .GT. 1.) GO TO 1
      RTOL = RTOL*200.
C
C Print a heading.
C
      IF (LUN .GE. 0) WRITE (LUN,100) N, FTOL
  100 FORMAT ('1',13X,'GETSIG -- N =',I4,', TOL = ',E10.3//)
C
C Initialize change counter ICNT and maximum change DSM for
C   the loop on arcs.
C
      ICNT = 0
      DSM = 0.
C
C Loop on arcs N1-N2 for which N2 > N1.  LPL points to the
C   last neighbor of N1.
C
      DO 10 N1 = 1,NM1
        LPL = LEND(N1)
        LP1 = LPL
C
C   Top of loop on neighbors N2 of N1.
C
    2   LP1 = LPTR(LP1)
        N2 = ABS(LIST(LP1))
        IF (N2 .LE. N1) GO TO 9
C
C Print a message and compute parameters for the arc:
C   nodal coordinates P1 and P2, arc-length AL,
C   UNORM = magnitude of P1 X P2, and
C   SIGIN = input SIGMA value.
C
        IF (LUN .GE. 0) WRITE (LUN,110) N1, N2
  110   FORMAT (/1X,'ARC',I4,' -',I4)
        P1(1) = X(N1)
        P1(2) = Y(N1)
        P1(3) = Z(N1)
        P2(1) = X(N2)
        P2(2) = Y(N2)
        P2(3) = Z(N2)
        AL = ARCLEN(P1,P2)
        UN(1) = P1(2)*P2(3) - P1(3)*P2(2)
        UN(2) = P1(3)*P2(1) - P1(1)*P2(3)
        UN(3) = P1(1)*P2(2) - P1(2)*P2(1)
        UNORM = SQRT(UN(1)*UN(1)+UN(2)*UN(2)+UN(3)*UN(3))
        IF (UNORM .EQ. 0.  .OR.  AL .EQ. 0.) GO TO 12
        SIGIN = SIGMA(LP1)
        IF (SIGIN .GE. SBIG) GO TO 9
C
C Compute scaled directional derivatives S1,S2 at the end-
C   points (for the direction N1->N2), first difference S,
C   and second differences D1,D2.
C
        S1 = AL*(GRAD(1,N1)*P2(1) + GRAD(2,N1)*P2(2) +
     .               GRAD(3,N1)*P2(3))/UNORM
        S2 = -AL*(GRAD(1,N2)*P1(1) + GRAD(2,N2)*P1(2) +
     .            GRAD(3,N2)*P1(3))/UNORM
        S = H(N2) - H(N1)
        D1 = S - S1
        D2 = S2 - S
        D1D2 = D1*D2
C
C Test for infinite tension required to satisfy either
C   property.
C
        SIG = SBIG
        IF ((D1D2 .EQ. 0.  .AND.  S1 .NE. S2)  .OR.
     .      (S .EQ. 0.  .AND.  S1*S2 .GT. 0.)) GO TO 8
C
C Test for SIGMA = 0 sufficient.  The data satisfies convex-
C   ity iff D1D2 .GE. 0, and D1D2 = 0 implies S1 = S = S2.
C
        SIG = 0.
        IF (D1D2 .LT. 0.) GO TO 4
        IF (D1D2 .EQ. 0.) GO TO 8
        T = MAX(D1/D2,D2/D1)
        IF (T .LE. 2.) GO TO 8
        TP1 = T + 1.
C
C Convexity:  find a zero of F(SIG) = SIG*coshm(SIG)/
C   sinhm(SIG) - TP1.
C
C   F(0) = 2-T < 0, F(TP1) .GE. 0, the derivative of F
C     vanishes at SIG = 0, and the second derivative of F is
C     .2 at SIG = 0.  A quadratic approximation is used to
C     obtain a starting point for the Newton method.
C
        SIG = SQRT(10.*T-20.)
        NIT = 0
C
C   Top of loop:
C
    3   IF (SIG .LE. .5) THEN
          CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
          T1 = COSHM/SINHM
          FP = T1 + SIG*(SIG/SINHM - T1*T1 + 1.)
        ELSE
C
C   Scale sinhm and coshm by 2*exp(-SIG) in order to avoid
C     overflow with large SIG.
C
          EMS = EXP(-SIG)
          SSM = 1. - EMS*(EMS+SIG+SIG)
          T1 = (1.-EMS)*(1.-EMS)/SSM
          FP = T1 + SIG*(2.*SIG*EMS/SSM - T1*T1 + 1.)
        ENDIF
C
        F = SIG*T1 - TP1
        IF (LUN .GE. 0) WRITE (LUN,120) SIG, F, FP
  120   FORMAT (1X,'CONVEXITY -- SIG = ',E15.8,
     .          ', F(SIG) = ',E15.8/1X,35X,'FP(SIG) = ',
     .          E15.8)
        NIT = NIT + 1
C
C   Test for convergence.
C
        IF (FP .LE. 0.) GO TO 8
        DSIG = -F/FP
        IF (ABS(DSIG) .LE. RTOL*SIG  .OR.  (F .GE. 0.  .AND.
     .      F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 8
C
C   Update SIG.
C
        SIG = SIG + DSIG
        GO TO 3
C
C Convexity cannot be satisfied.  Monotonicity can be satis-
C   fied iff S1*S .GE. 0 and S2*S .GE. 0 since S .NE. 0.
C
    4   IF (S1*S .LT. 0.  .OR.  S2*S .LT. 0.) GO TO 8
        T0 = 3.*S - S1 - S2
        D0 = T0*T0 - S1*S2
C
C SIGMA = 0 is sufficient for monotonicity iff S*T0 .GE. 0
C   or D0 .LE. 0.
C
        IF (D0 .LE. 0.  .OR.  S*T0 .GE. 0.) GO TO 8
C
C Monotonicity:  find a zero of F(SIG) = sign(S)*HP(R),
C   where HPP(R) = 0 and HP, HPP denote derivatives of H.
C   F has a unique zero, F(0) < 0, and F approaches
C   abs(S) as SIG increases.
C
C   Initialize parameters for the secant method.  The method
C     uses three points:  (SG0,F0), (SIG,F), and
C     (SNEG,FNEG), where SG0 and SNEG are defined implicitly
C     by DSIG = SIG - SG0 and DMAX = SIG - SNEG.
C
        SGN = SIGN(1.D0,S)
        SIG = SBIG
        FMAX = SGN*(SIG*S-S1-S2)/(SIG-2.)
        IF (FMAX .LE. 0.) GO TO 8
        STOL = RTOL*SIG
        F = FMAX
        F0 = SGN*D0/(3.*(D1-D2))
        FNEG = F0
        DSIG = SIG
        DMAX = SIG
        D1PD2 = D1 + D2
        NIT = 0
C
C   Top of loop:  compute the change in SIG by linear
C     interpolation.
C
    5   DSIG = -F*DSIG/(F-F0)
        IF (LUN .GE. 0) WRITE (LUN,130) DSIG
  130   FORMAT (1X,'MONOTONICITY -- DSIG = ',E15.8)
        IF ( ABS(DSIG) .GT. ABS(DMAX)  .OR.
     .       DSIG*DMAX .GT. 0. ) GO TO 7
C
C   Restrict the step-size such that abs(DSIG) .GE. STOL/2.
C     Note that DSIG and DMAX have opposite signs.
C
        IF (ABS(DSIG) .LT. STOL/2.) DSIG = -SIGN(STOL/2.,
     .                              DMAX)
C
C   Update SIG, F0, and F.
C
        SIG = SIG + DSIG
        F0 = F
        IF (SIG .LE. .5) THEN
C
C   Use approximations to the hyperbolic functions designed
C     to avoid cancellation error with small SIG.
C
          CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
          C1 = SIG*COSHM*D2 - SINHM*D1PD2
          C2 = SIG*(SINHM+SIG)*D2 - COSHM*D1PD2
          A = C2 - C1
          E = SIG*SINHM - COSHMM - COSHMM
        ELSE
C
C   Scale sinhm and coshm by 2*exp(-SIG) in order to avoid
C     overflow with large SIG.
C
          EMS = EXP(-SIG)
          EMS2 = EMS + EMS
          TM = 1. - EMS
          SSINH = TM*(1.+EMS)
          SSM = SSINH - SIG*EMS2
          SCM = TM*TM
          C1 = SIG*SCM*D2 - SSM*D1PD2
          C2 = SIG*SSINH*D2 - SCM*D1PD2
C
C   R is in (0,1) and well-defined iff HPP(T1)*HPP(T2) < 0.
C
          F = FMAX
          IF (C1*(SIG*SCM*D1 - SSM*D1PD2) .GE. 0.) GO TO 6
          A = EMS2*(SIG*TM*D2 + (TM-SIG)*D1PD2)
          IF (A*(C2+C1) .LT. 0.) GO TO 6
          E = SIG*SSINH - SCM - SCM
        ENDIF
C
        F = (SGN*(E*S2-C2) + SQRT(A*(C2+C1)))/E
C
C   Update the number of iterations NIT.
C
    6   NIT = NIT + 1
        IF (LUN .GE. 0) WRITE (LUN,140) NIT, SIG, F
  140   FORMAT (1X,11X,I2,' -- SIG = ',E15.8,', F = ',
     .          E15.8)
C
C   Test for convergence.
C
        STOL = RTOL*SIG
        IF (ABS(DMAX) .LE. STOL  .OR.  (F .GE. 0.  .AND.
     .      F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 8
        DMAX = DMAX + DSIG
        IF (F0*F .GT. 0.  .AND.  ABS(F) .GE. ABS(F0))
     .     GO TO 7
        IF (F0*F .LE. 0.) THEN
C
C   F and F0 have opposite signs.  Update (SNEG,FNEG) to
C     (SG0,F0) so that F and FNEG always have opposite
C     signs.  If SIG is closer to SNEG than SG0 and abs(F)
C     < abs(FNEG), then swap (SNEG,FNEG) with (SG0,F0).
C
          T1 = DMAX
          T2 = FNEG
          DMAX = DSIG
          FNEG = F0
          IF ( ABS(DSIG) .GT. ABS(T1)  .AND.
     .         ABS(F) .LT. ABS(T2) ) THEN
C
            DSIG = T1
            F0 = T2
          ENDIF
        ENDIF
        GO TO 5
C
C   Bottom of loop:  F0*F > 0 and the new estimate would
C     be outside of the bracketing interval of length
C     abs(DMAX).  Reset (SG0,F0) to (SNEG,FNEG).
C
    7   DSIG = DMAX
        F0 = FNEG
        GO TO 5
C
C  Update SIGMA, ICNT, and DSM if necessary.
C
    8   SIG = MIN(SIG,SBIG)
        IF (SIG .GT. SIGIN) THEN
          SIGMA(LP1) = SIG
          LP2 = LSTPTR(LEND(N2),N1,LIST,LPTR)
          SIGMA(LP2) = SIG
          ICNT = ICNT + 1
          DSM = MAX(DSM,SIG-SIGIN)
        ENDIF
C
C Bottom of loop on neighbors N2 of N1.
C
    9   IF (LP1 .NE. LPL) GO TO 2
   10   CONTINUE
C
C No errors encountered.
C
      DSMAX = DSM
      IER = ICNT
      RETURN
C
C N < 3.
C
   11 DSMAX = 0.
      IER = -1
      RETURN
C
C Nodes N1 and N2 coincide.
C
   12 DSMAX = DSM
      IER = -2
      RETURN
      END
      SUBROUTINE GIVENS ( A,B, C,S)
      DOUBLE PRECISION A, B, C, S
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   05/09/92
C
C   This subroutine constructs the Givens plane rotation,
C
C           ( C  S)
C       G = (     ) , where C*C + S*S = 1,
C           (-S  C)
C
C which zeros the second component of the vector (A,B)**T
C (transposed).  Subroutine ROTATE may be called to apply
C the transformation to a 2 by N matrix.
C
C   This routine is identical to Subroutine SROTG from the
C LINPACK BLAS (Basic Linear Algebra Subroutines).
C
C On input:
C
C       A,B = Components of the vector defining the rota-
C             tion.  These are overwritten by values R
C             and Z (described below) which define C and S.
C
C On output:
C
C       A = Signed Euclidean norm R of the input vector:
C           R = +/-SQRT(A*A + B*B)
C
C       B = Value Z such that:
C             C = SQRT(1-Z*Z) and S=Z if ABS(Z) .LE. 1, and
C             C = 1/Z and S = SQRT(1-C*C) if ABS(Z) > 1.
C
C       C = +/-(A/R) or 1 if R = 0.
C
C       S = +/-(B/R) or 0 if R = 0.
C
C Modules required by GIVENS:  None
C
C Intrinsic functions called by GIVENS:  ABS, SQRT
C
C***********************************************************
C
      DOUBLE PRECISION AA, BB, R, U, V
C
C Local parameters:
C
C AA,BB = Local copies of A and B
C R =     C*A + S*B = +/-SQRT(A*A+B*B)
C U,V =   Variables used to scale A and B for computing R
C
      AA = A
      BB = B
      IF (ABS(AA) .GT. ABS(BB)) THEN
C
C ABS(A) > ABS(B).
C
        U = AA + AA
        V = BB/U
        R = SQRT(.25 + V*V) * U
        C = AA/R
        S = V * (C + C)
C
C Note that R has the sign of A, C > 0, and S has
C   SIGN(A)*SIGN(B).
C
        B = S
        A = R
      ELSEIF (BB .NE. 0.) THEN
C
C ABS(A) .LE. ABS(B).
C
        U = BB + BB
        V = AA/U
C
C Store R in A.
C
        A = SQRT(.25 + V*V) * U
        S = BB/A
        C = V * (S + S)
C
C Note that R has the sign of B, S > 0, and C has
C   SIGN(A)*SIGN(B).
C
        B = 1.
        IF (C .NE. 0.) B = 1./C
      ELSE
C
C A = B = 0.
C
        C = 1.
        S = 0.
      ENDIF
      RETURN
      END
      SUBROUTINE GRADG (N,X,Y,Z,F,LIST,LPTR,LEND,IFLGS,
     .                  SIGMA, NIT,DGMAX,GRAD, IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), IFLGS, NIT, IER
      DOUBLE PRECISION    X(N), Y(N), Z(N), F(N), SIGMA(*), DGMAX,
     .        GRAD(3,N)
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/24/96
C
C   Given a triangulation of N nodes on the unit sphere with
C data values F at the nodes and tension factors SIGMA asso-
C ciated with the arcs, this routine uses a global method
C to compute estimated gradients at the nodes.  The method
C consists of minimizing a quadratic functional Q(G) over
C the N-vector G of gradients, where Q approximates the
C linearized curvature of the restriction to arcs of the
C interpolatory function F defined by Function FVAL.  The
C restriction of F to an arc of the triangulation is the
C Hermite interpolatory tension spline defined by the data
C values and tangential gradient components at the endpoints
C of the arc.  Letting D1F(A) and D2F(A) denote first and
C second derivatives of F with respect to a parameter A var-
C ying along a triangulation arc, Q is the sum of integrals
C over the arcs of D2F(A)**2 + ((SIGMA/L)*(D1F(A)-S))**2,
C where L denotes arc-length, SIGMA is the appropriate ten-
C sion factor, and S is the slope of the linear function of
C A which interpolates the values of F at the endpoints of
C the arc.
C
C   Since the gradient at node K lies in the plane tangent
C to the sphere surface at K, it is effectively defined by
C only two components -- its X and Y components in the coor-
C dinate system obtained by rotating K to the north pole.
C Thus, the minimization problem corresponds to an order-2N
C symmetric positive-definite sparse linear system which is
C solved by a block Gauss-Seidel method with 2 by 2 blocks.
C
C   An alternative method, Subroutine GRADL, computes a
C local approximation to the gradient at a single node and,
C although less efficient when all gradients are needed, was
C found to be generally more accurate (in the case of uni-
C form zero tension) when the nodal distribution is very
C dense, varies greatly, or does not cover the sphere.
C GRADG, on the other hand, was found to be slightly more
C accurate on a uniform distribution of 514 nodes.
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y,Z = Arrays of length N containing Cartesian
C               coordinates of the nodes.  X(I)**2 + Y(I)**2
C               + Z(I)**2 = 1 for I = 1,...,N.
C
C       F = Array of length N containing data values at the
C           nodes.  F(I) is associated with (X(I),Y(I),Z(I))
C           for I = 1,...,N.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to STRIPACK
C                        Subroutine TRMESH.
C
C       IFLGS = Tension factor option:
C               IFLGS .LE. 0 if a single uniform tension
C                            factor is to be used.
C               IFLGS .GE. 1 if variable tension is desired.
C
C       SIGMA = Uniform tension factor (IFLGS .LE. 0), or
C               array containing tension factors associated
C               with arcs in one-to-one correspondence with
C               LIST entries (IFLGS .GE. 1).  Refer to Sub-
C               programs GETSIG, SIG0, SIG1, and SIG2.
C
C The above parameters are not altered by this routine.
C
C       NIT = Maximum number of Gauss-Seidel iterations to
C             be applied.  This maximum will likely be a-
C             chieved if DGMAX is smaller than the machine
C             precision.  NIT .GE. 0.
C
C       DGMAX = Nonnegative convergence criterion.  The
C               method is terminated when the maximum change
C               in a gradient between iterations is at most
C               DGMAX.  The change in a gradient is taken to
C               be the Euclidean norm of the difference (in
C               the rotated coordinate system) relative to 1
C               plus the norm of the old gradient value.
C
C       GRAD = 3 by N array whose columns contain initial
C              solution estimates (zero vectors are suffici-
C              ent).  GRAD(I,J) contains component I of the
C              gradient at node J for I = 1,2,3 (X,Y,Z) and
C              J = 1,...,N.  GRAD( ,J) must be orthogonal to
C              node J -- GRAD(1,J)*X(J) + GRAD(2,J)*Y(J) +
C              GRAD(3,J)*Z(J) = 0.
C
C On output:
C
C       NIT = Number of Gauss-Seidel iterations employed.
C
C       DGMAX = Maximum change in a gradient at the last
C               iteration.
C
C       GRAD = Estimated gradients.  See the description
C              under input parameters.  GRAD is not changed
C              if IER = -1.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered and the
C                     convergence criterion was achieved.
C             IER = 1 if no errors were encountered but con-
C                     vergence was not achieved within NIT
C                     iterations.
C             IER = -1 if N or DGMAX is outside its valid
C                      range or NIT .LT. 0 on input.
C             IER = -2 if all nodes are collinear or the
C                      triangulation is invalid.
C             IER = -3 if duplicate nodes were encountered.
C
C SSRFPACK modules required by GRADG:  APLYRT, CONSTR,
C                                        GRCOEF, SNHCSH
C
C Intrinsic functions called by GRADG:  ATAN, MAX, SQRT
C
C***********************************************************
C
      INTEGER IFL, ITER, J, K, LPJ, LPL, MAXIT, NN
      DOUBLE PRECISION    ALFA, A11, A12, A22, CX, CY, D, DEN, DET,
     .        DGK(3), DGMX, DG1, DG2, FK, G1, G2, G3, R1,
     .        R2, SD, SIG, SINAL, SX, SY, T, TOL, XK, YK,
     .        ZK, XJ, YJ, ZJ, XS, YS
C
C Local parameters:
C
C ALFA =        Arc-length between nodes K and J
C A11,A12,A22 = Matrix components of the 2 by 2 block A*DG
C                 = R where A is symmetric, (DG1,DG2,0) is
C                 the change in the gradient at K, and R is
C                 the residual
C CX,CY =       Components of a rotation mapping K to the
C                 north pole (0,0,1)
C D =           Function of SIG computed by GRCOEF -- factor
C                 in the order-2 system
C DEN =         ALFA*SINAL**2 -- factor in the 2 by 2 system
C DET =         Determinant of the order-2 matrix
C DGK =         Change in GRAD( ,K) from the previous esti-
C                 mate in the original coordinate system
C DGMX =        Maximum change in a gradient between itera-
C                 tions
C DG1,DG2 =     Solution of the 2 by 2 system -- first 2
C                 components of DGK in the rotated coordi-
C                 nate system
C FK =          Data value F(K)
C G1,G2,G3 =    Components of GRAD( ,K)
C IFL =         Local copy of IFLGS
C ITER =        Number of iterations used
C J =           Neighbor of K
C K =           DO-loop and node index
C LPJ =         LIST pointer of node J as a neighbor of K
C LPL =         Pointer to the last neighbor of K
C MAXIT =       Input value of NIT
C NN =          Local copy of N
C R1,R2 =       Components of the residual -- derivatives of
C                 Q with respect to the components of the
C                 gradient at node K
C SD =          Function of SIG computed by GRCOEF -- factor
C                 in the order-2 system
C SIG =         Tension factor associated with ARC K-J
C SINAL =       SIN(ALFA) -- magnitude of the vector cross
C                 product between nodes K and J
C SX,SY =       Components of a rotation mapping K to the
C                 north pole (0,0,1)
C T =           Temporary storage for factors in the system
C                 components
C TOL =         Local copy of DGMAX
C XK,YK,ZK =    Coordinates of node K -- X(K), Y(K), Z(K)
C XJ,YJ,ZJ =    Coordinates of node J in the rotated coor-
C                 dinate system
C XS,YS =       XJ**2, YJ**2
C
      NN = N
      IFL = IFLGS
      MAXIT = NIT
      TOL = DGMAX
C
C Test for errors in input, and initialize iteration count,
C   tension factor, and output value of DGMAX.
C
      IF (NN .LT. 3  .OR.  MAXIT .LT. 0  .OR.  TOL .LT. 0.)
     .   GO TO 11
      ITER = 0
      SIG = SIGMA(1)
      DGMX = 0.
C
C Top of iteration loop.
C
    1 IF (ITER .EQ. MAXIT) GO TO 4
      DGMX = 0.
C
C Loop on nodes.
C
      DO 3 K = 1,NN
        XK = X(K)
        YK = Y(K)
        ZK = Z(K)
        FK = F(K)
        G1 = GRAD(1,K)
        G2 = GRAD(2,K)
        G3 = GRAD(3,K)
C
C   Construct the rotation mapping node K to the north pole.
C
        CALL CONSTR (XK,YK,ZK, CX,SX,CY,SY)
C
C   Initialize components of the 2 by 2 system for the
C     change (DG1,DG2,0) in the K-th solution components
C     (symmetric matrix in A and residual in R).
C
        A11 = 0.
        A12 = 0.
        A22 = 0.
        R1 = 0.
        R2 = 0.
C
C   Loop on neighbors J of node K.
C
        LPL = LEND(K)
        LPJ = LPL
    2   LPJ = LPTR(LPJ)
          J = ABS(LIST(LPJ))
C
C   Compute the coordinates of J in the rotated system.
C
          T = SX*Y(J) + CX*Z(J)
          YJ = CX*Y(J) - SX*Z(J)
          ZJ = SY*X(J) + CY*T
          XJ = CY*X(J) - SY*T
C
C   Compute arc-length ALFA between J and K, SINAL =
C     SIN(ALFA), and DEN = ALFA*SIN(ALFA)**2.
C
          ALFA = 2.*ATAN(SQRT((1.-ZJ)/(1.+ZJ)))
          XS = XJ*XJ
          YS = YJ*YJ
          SINAL = SQRT(XS+YS)
          DEN = ALFA*(XS+YS)
C
C   Test for coincident nodes and compute functions of SIG:
C     D = SIG*(SIG*COSHM-SINHM)/E and SD = SIG*SINHM/E for
C     E = SIG*SINH-2*COSHM.
C
          IF (DEN .EQ. 0.) GO TO 13
          IF (IFL .GE. 1) SIG = SIGMA(LPJ)
          CALL GRCOEF (SIG, D,SD)
C
C   Update the system components for node J.
C
          T = D/DEN
          A11 = A11 + T*XS
          A12 = A12 + T*XJ*YJ
          A22 = A22 + T*YS
          T = (D+SD)*(FK-F(J))/(ALFA*ALFA*SINAL) +
     .        ( D*(G1*X(J) + G2*Y(J) + G3*Z(J)) -
     .          SD*(GRAD(1,J)*XK + GRAD(2,J)*YK +
     .                 GRAD(3,J)*ZK) )/DEN
          R1 = R1 - T*XJ
          R2 = R2 - T*YJ
C
C   Bottom of loop on neighbors.
C
          IF (LPJ .NE. LPL) GO TO 2
C
C   Solve the 2 by 2 system and update DGMAX.
C
        DET = A11*A22 - A12*A12
        IF (DET .EQ. 0.  .OR.  A11 .EQ. 0.) GO TO 12
        DG2 = (A11*R2 - A12*R1)/DET
        DG1 = (R1 - A12*DG2)/A11
        DGMX = MAX(DGMX,SQRT(DG1*DG1+DG2*DG2)/
     .             (1.+SQRT(G1*G1+G2*G2+G3*G3)))
C
C   Rotate (DG1,DG2,0) back to the original coordinate
C     system and update GRAD( ,K).
C
        CALL APLYRT (DG1,DG2,CX,SX,CY,SY, DGK)
        GRAD(1,K) = G1 + DGK(1)
        GRAD(2,K) = G2 + DGK(2)
        GRAD(3,K) = G3 + DGK(3)
    3   CONTINUE
C
C   Increment ITER and test for convergence.
C
      ITER = ITER + 1
      IF (DGMX .GT. TOL) GO TO 1
C
C The method converged.
C
      NIT = ITER
      DGMAX = DGMX
      IER = 0
      RETURN
C
C The method failed to converge within NIT iterations.
C
    4 DGMAX = DGMX
      IER = 1
      RETURN
C
C Invalid input parameter.
C
   11 NIT = 0
      DGMAX = 0.
      IER = -1
      RETURN
C
C Node K and its neighbors are collinear.
C
   12 NIT = 0
      DGMAX = DGMX
      IER = -2
      RETURN
C
C Nodes K and J coincide.
C
   13 NIT = 0
      DGMAX = DGMX
      IER = -3
      RETURN
      END
      SUBROUTINE GRADL (N,K,X,Y,Z,W,LIST,LPTR,LEND, G,IER)
      INTEGER N, K, LIST(*), LPTR(*), LEND(N), IER
      DOUBLE PRECISION    X(N), Y(N), Z(N), W(N), G(3)
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/24/96
C
C   Given a triangulation of a set of nodes on the unit
C sphere with their associated data values W, this routine
C estimates a gradient vector at node K as follows:  the
C coordinate system is rotated so that K becomes the north
C pole, node K and a set of nearby nodes are projected
C orthogonally onto the X-Y plane (in the new coordinate
C system), a quadratic is fitted in a weighted least squares
C sense to the data values at the projected nodes such that
C the value (associated with K) at (0,0) is interpolated, X
C and Y-partial derivative estimates DX and DY are computed
C by differentiating the quadratic at (0,0), and the esti-
C mated gradient G is obtained by rotating (DX,DY,0) back to
C the original coordinate system.  Note that G lies in the
C plane tangent to the sphere at node K, i.e., G is orthogo-
C nal to the unit vector represented by node K.  A Marquardt
C stabilization factor is used if necessary to ensure a
C well-conditioned least squares system, and a unique solu-
C tion exists unless the nodes are collinear.
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 7.
C
C       K = Node at which the gradient is sought.  1 .LE. K
C           .LE. N.
C
C       X,Y,Z = Arrays containing the Cartesian coordinates
C               of the nodes.
C
C       W = Array containing the data values at the nodes.
C           W(I) is associated with (X(I),Y(I),Z(I)) for
C           I = 1,...,N.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to STRIPACK
C                        Subroutine TRMESH.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       G = X, Y, and Z components (in that order) of the
C           estimated gradient at node K unless IER < 0.
C
C       IER = Error indicator:
C             IER .GE. 6 if no errors were encountered.
C                        IER contains the number of nodes
C                        (including K) used in the least
C                        squares fit.
C             IER = -1 if N or K is outside its valid range.
C             IER = -2 if the least squares system has no
C                      unique solution due to duplicate or
C                      collinear nodes.
C
C STRIPACK module required by GRADL:  GETNP
C
C SSRFPACK modules required by GRADL:  APLYR, APLYRT,
C                                        CONSTR, GIVENS,
C                                        ROTATE, SETUP
C
C Intrinsic functions called by GRADL:  ABS, MIN, REAL, SQRT
C
C***********************************************************
C
      INTEGER   LMN, LMX
      PARAMETER (LMN=10,  LMX=30)
      INTEGER I, IERR, IM1, IP1, J, JP1, KK, L, LM1, LMAX,
     .        LMIN, LNP, NN, NP, NPTS(LMX)
      DOUBLE PRECISION    A(6,6), AV, AVSQ, C, CX, CY, DF, DMIN, DTOL,
     .        DX, DY, RF, RIN, RTOL, S, SF, SUM, SX, SY,
     .        WK, WT, XP, YP, ZP
C
      DATA    RTOL/1.E-6/, DTOL/.01/, SF/1./
C
C Local parameters:
C
C A =         Transpose of the (upper triangle of the) aug-
C               mented regression matrix
C AV =        Root-mean-square distance (in the rotated
C               coordinate system) between the origin and
C               the nodes (other than K) in the least
C               squares fit.  The first 3 columns of A**T
C               are scaled by 1/AVSQ, the next 2 by 1/AV.
C AVSQ =      AV*AV:  accumulated in SUM
C C,S =       Components of the plane rotation used to
C               triangularize the regression matrix
C CX,SX =     Components of a plane rotation about the X-
C               axis which, together with CY and SY, define
C               a mapping from node K to the north pole
C               (0,0,1)
C CY,SY =     Components of a plane rotation about the Y-
C               axis
C DF =        Negative Z component (in the rotated coordi-
C               nate system) of an element NP of NPTS --
C               increasing function of the angular distance
C               between K and NP.  DF lies in the interval
C               (-1,1).
C DMIN =      Minimum of the magnitudes of the diagonal
C               elements of the triangularized regression
C               matrix
C DTOL =      Tolerance for detecting an ill-conditioned
C               system (DMIN is required to be at least
C               DTOL)
C DX,DY =     X and Y components of the estimated gradient
C               in the rotated coordinate system
C I,J =       Loop indexes
C IERR =      Error flag for calls to GETNP (not checked)
C IM1,IP1 =   I-1, I+1
C JP1 =       J+1
C KK =        Local copy of K
C L =         Number of columns of A**T to which a rotation
C               is applied
C LM1 =       LMIN-1
C LMIN,LMAX = Min(LMN,N), Min(LMX,N)
C LMN,LMX =   Minimum and maximum values of LNP for N
C               sufficiently large.  In most cases LMN-1
C               nodes are used in the fit.  7 .LE. LMN .LE.
C               LMX.
C LNP =       Length of NPTS or LMAX+1
C NN =        Local copy of N
C NP =        Element of NPTS to be added to the system
C NPTS =      Array containing the indexes of a sequence of
C               nodes ordered by angular distance from K.
C               NPTS(1)=K and the first LNP-1 elements of
C               NPTS are used in the least squares fit.
C               unless LNP = LMAX+1, NPTS(LNP) determines R
C               (see RIN).
C RF =        Value of DF associated with NPTS(LNP) unless
C               LNP = LMAX+1 (see RIN)
C RIN =       Inverse of a radius of influence R which
C               enters into WT:  R = 1+RF unless all ele-
C               ments of NPTS are used in the fit (LNP =
C               LMAX+1), in which case R is the distance
C               function associated with some point more
C               distant from K than NPTS(LMAX)
C RTOL =      Tolerance for determining LNP (and hence R):
C               if the increase in DF between two successive
C               elements of NPTS is less than RTOL, they are
C               treated as being the same distance from node
C               K and an additional node is added
C SF =        Marquardt stabilization factor used to damp
C               out the first 3 solution components (second
C               partials of the quadratic) when the system
C               is ill-conditioned.  Increasing SF results
C               in more damping (a more nearly linear fit).
C SUM =       Sum of squared Euclidean distances (in the
C               rotated coordinate system) between the
C               origin and the nodes used in the least
C               squares fit
C WK =        W(K) -- data value at node K
C WT =        Weight for the equation coreesponding to NP:
C               WT = (R-D)/(R*D) = 1/D - RIN, where D = 1-ZP
C               is associated with NP
C XP,YP,ZP =  Coordinates of NP in the rotated coordinate
C               system unless ZP < 0, in which case
C               (XP,YP,0) lies on the equator
C
      NN = N
      KK = K
      WK = W(KK)
C
C Check for errors and initialize LMIN, LMAX.
C
      IF (NN .LT. 7  .OR.  KK .LT. 1  .OR.  KK .GT. NN)
     .   GO TO 13
      LMIN = MIN(LMN,NN)
      LMAX = MIN(LMX,NN)
C
C Compute NPTS, LNP, AVSQ, AV, and R.
C   Set NPTS to the closest LMIN-1 nodes to K.  DF contains
C   the negative Z component (in the rotated coordinate
C   system) of the new node on return from GETNP.
C
      SUM = 0.
      NPTS(1) = KK
      LM1 = LMIN - 1
      DO 1 LNP = 2,LM1
        CALL GETNP (X,Y,Z,LIST,LPTR,LEND,LNP, NPTS, DF,IERR)
        SUM = SUM + 1. - DF*DF
    1   CONTINUE
C
C   Add additional nodes to NPTS until the increase in
C     R = 1+RF is at least RTOL.
C
      DO 2 LNP = LMIN,LMAX
        CALL GETNP (X,Y,Z,LIST,LPTR,LEND,LNP, NPTS, RF,IERR)
        IF (RF-DF .GE. RTOL) GO TO 3
        SUM = SUM + 1. - RF*RF
    2   CONTINUE
C
C   Use all LMAX nodes in the least squares fit.  R is
C     arbitrarily increased by 5 percent.
C
      RF = 1.05*RF + .05
      LNP = LMAX + 1
C
C   There are LNP-2 equations corresponding to nodes
C     NPTS(2),...,NPTS(LNP-1).
C
    3 AVSQ = SUM/REAL(LNP-2)
      AV = SQRT(AVSQ)
      RIN = 1./(1.+RF)
C
C Construct the rotation.
C
      CALL CONSTR (X(KK),Y(KK),Z(KK), CX,SX,CY,SY)
C
C Set up the first 5 equations of the augmented regression
C   matrix (transposed) as the columns of A, and zero out
C   the lower triangle (upper triangle of A) with Givens
C   rotations.
C
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
          CALL ROTATE (L,C,S, A(JP1,J),A(JP1,I) )
    4     CONTINUE
    5   CONTINUE
C
C Add the additional equations to the system using
C   the last column of A.  I .LE. LNP.
C
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
          CALL ROTATE (L,C,S, A(JP1,J),A(JP1,6) )
    7     CONTINUE
        I = I + 1
        GO TO 6
C
C Test the system for ill-conditioning.
C
    8 DMIN = MIN( ABS(A(1,1)),ABS(A(2,2)),ABS(A(3,3)),
     .            ABS(A(4,4)),ABS(A(5,5)) )
      IF (DMIN .GE. DTOL) GO TO 12
      IF (LNP .LE. LMAX) THEN
C
C Add another node to the system and increase R.
C   I = LNP.
C
        LNP = LNP + 1
        IF (LNP .LE. LMAX) CALL GETNP (X,Y,Z,LIST,LPTR,LEND,
     .                                 LNP,NPTS, RF,IERR)
        RIN = 1./(1.05*(1.+RF))
        GO TO 6
      ENDIF
C
C Stabilize the system by damping second partials.  Add
C   multiples of the first three unit vectors to the first
C   three equations.
C
      DO 11 I = 1,3
        A(I,6) = SF
        IP1 = I + 1
        DO 9 J = IP1,6
          A(J,6) = 0.
    9     CONTINUE
        DO 10 J = I,5
          JP1 = J + 1
          L = 6 - J
          CALL GIVENS ( A(J,J),A(J,6), C,S)
          CALL ROTATE (L,C,S, A(JP1,J),A(JP1,6) )
   10     CONTINUE
   11   CONTINUE
C
C Test the linear portion of the stabilized system for
C   ill-conditioning.
C
      DMIN = MIN( ABS(A(4,4)),ABS(A(5,5)) )
      IF (DMIN .LT. DTOL) GO TO 14
C
C Solve the 2 by 2 triangular system for the estimated
C   partial derivatives.
C
   12 DY = A(6,5)/A(5,5)
      DX = (A(6,4) - A(5,4)*DY)/A(4,4)/AV
      DY = DY/AV
C
C Rotate the gradient (DX,DY,0) back into the original
C   coordinate system.
C
      CALL APLYRT (DX,DY,CX,SX,CY,SY, G)
      IER = LNP - 1
      RETURN
C
C N or K is outside its valid range.
C
   13 IER = -1
      RETURN
C
C No unique solution due to collinear nodes.
C
   14 IER = -2
      RETURN
      END
      SUBROUTINE GRCOEF (SIGMA, D,SD)
      DOUBLE PRECISION SIGMA, D, SD
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/21/96
C
C   This subroutine computes factors involved in the linear
C systems solved by Subroutines GRADG and SMSGS.
C
C On input:
C
C       SIGMA = Nonnegative tension factor associated with a
C               triangulation arc.
C
C SIGMA is not altered by this routine.
C
C On output:
C
C       D = Diagonal factor.  D = SIG*(SIG*Coshm(SIG) -
C           Sinhm(SIG))/E where E = SIG*Sinh(SIG) - 2*
C           Coshm(SIG).  D > 0, and D = 4 at SIG = 0.
C
C       SD = Off-diagonal factor.  SD = SIG*Sinhm(SIG)/E.
C            SD > 0, and SD = 2 at SIG = 0.
C
C SSRFPACK module required by GRCOEF:  SNHCSH
C
C Intrinsic function called by GRCOEF:  EXP
C
C***********************************************************
C
      DOUBLE PRECISION COSHM, COSHMM, E, EMS, SCM, SIG, SINHM, SSINH,
     .     SSM
      SIG = SIGMA
      IF (SIG .LT. 1.E-9) THEN
C
C Cubic function:
C
        D = 4.
        SD = 2.
      ELSEIF (SIG .LE. .5) THEN
C
C 0 < SIG .LE. .5.
C
C Use approximations designed to avoid cancellation error
C   in the hyperbolic functions when SIGMA is small.
C
        CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
        E = SIG*SINHM - COSHMM - COSHMM
        D = SIG*(SIG*COSHM-SINHM)/E
        SD = SIG*SINHM/E
      ELSE
C
C SIG > .5.
C
C Scale SINHM, COSHM, and E by 2*EXP(-SIG) in order to
C   avoid overflow when SIGMA is large.
C
        EMS = EXP(-SIG)
        SSINH = 1. - EMS*EMS
        SSM = SSINH - 2.*SIG*EMS
        SCM = (1.-EMS)*(1.-EMS)
        E = SIG*SSINH - SCM - SCM
        D = SIG*(SIG*SCM-SSM)/E
        SD = SIG*SSM/E
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION HVAL (B,H1,H2,HP1,HP2,SIGMA)
      DOUBLE PRECISION B, H1, H2, HP1, HP2, SIGMA
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/21/96
C
C   Given a line segment P1-P2 containing a point P, along
C with values and derivatives at the endpoints, this func-
C tion returns the value H(P), where H is the Hermite inter-
C polatory tension spline defined by the endpoint data.
C
C On input:
C
C       B = Local coordinate of P with respect to P1-P2:
C           P = B*P1 + (1-B)*P2, and thus B = d(P,P2)/
C           d(P1,P2), where d(P1,P2) is the distance between
C           P1 and P2.  B < 0 or B > 1 results in extrapola-
C           tion.
C
C       H1,H2 = Values interpolated at P1 and P2, respec-
C               tively.
C
C       HP1,HP2 = Products of d(P1,P2) with first order der-
C                 ivatives at P1 and P2, respectively.  HP1
C                 may, for example, be the scalar product of
C                 P2-P1 with a gradient at P1.
C
C       SIGMA = Nonnegative tension factor associated with
C               the spline.  SIGMA = 0 corresponds to a
C               cubic spline, and H approaches the linear
C               interpolant of H1 and H2 as SIGMA increases.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       HVAL = Interpolated value H(P).
C
C SSRFPACK module required by HVAL:  SNHCSH
C
C Intrinsic functions called by HVAL:  ABS, EXP
C
C***********************************************************
C
      DOUBLE PRECISION B1, B2, CM, CM2, CMM, D1, D2, DUMMY, E, E1, E2,
     .     EMS, S, SB1, SB2, SIG, SM, SM2, TM, TM1, TM2, TS
      B1 = B
      B2 = 1. - B1
C
C Compute slope S and second differences D1 and D2 scaled
C   by the separation between P1 and P2.
C
      S = H2 - H1
      D1 = S - HP1
      D2 = HP2 - S
C
C Test the range of SIGMA.
C
      SIG = ABS(SIGMA)
      IF (SIG .LT. 1.E-9) THEN
C
C Hermite cubic interpolation:
C
        HVAL = H1 + B2*(HP1 + B2*(D1 + B1*(D1 - D2)))
      ELSEIF (SIG .LE. .5) THEN
C
C 0 < SIG .LE. .5.  Use approximations designed to avoid
C   cancellation error in the hyperbolic functions.
C
        SB2 = SIG*B2
        CALL SNHCSH (SIG, SM,CM,CMM)
        CALL SNHCSH (SB2, SM2,CM2,DUMMY)
        E = SIG*SM - CMM - CMM
        HVAL = H1 + B2*HP1 + ((CM*SM2-SM*CM2)*(D1+D2) + SIG*
     .                     (CM*CM2-(SM+SIG)*SM2)*D1)/(SIG*E)
      ELSE
C
C SIG > .5.  Use negative exponentials in order to avoid
C   overflow.  Note that EMS = EXP(-SIG).
C
        SB1 = SIG*B1
        SB2 = SIG - SB1
        E1 = EXP(-SB1)
        E2 = EXP(-SB2)
        EMS = E1*E2
        TM = 1. - EMS
        TS = TM*TM
        TM1 = 1. - E1
        TM2 = 1. - E2
        E = TM*(SIG*(1.+EMS) - TM - TM)
        HVAL = H1 + B2*S + (TM*TM1*TM2*(D1+D2) + SIG*
     .                      ((E2*TM1*TM1-B1*TS)*D1 +
     .                       (E1*TM2*TM2-B2*TS)*D2))/(SIG*E)
      ENDIF
      RETURN
      END
      SUBROUTINE INTRNN (N,PLAT,PLON,X,Y,Z,W,LIST,LPTR,LEND, 
     .                   IST,PW,IER)
      INTEGER N, LIST(*),LPTR(*), LEND(N), IST, IER
      DOUBLE PRECISION  PLAT, PLON, X(N), Y(N), Z(N), W(N), PW
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
      DOUBLE PRECISION    P(3), P1(3), P2(3), P3(3), B1,B2,B3, 
     .        DIST(3), ARCLEN
      IF (N .LT. 3  .OR.  IST .LT. 1  .OR.  IST .GT. N)
     .    GO TO 11
C
C Transform (PLAT,PLON) to Cartesian coordinates.
C
      P(1) = COS(PLAT)*COS(PLON)
      P(2) = COS(PLAT)*SIN(PLON)
      P(3) = SIN(PLAT)
C
C Find the vertex indexes of a triangle containing P.
C
      CALL TRFIND(IST,P,N,X,Y,Z,LIST,LPTR,LEND, B1,B2,B3,
     .            I1,I2,I3)
      IF (I1 .EQ. 0) GO TO 12
      IST = I1
      IF (I3 .LE. 0) GO TO 13
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
C
C N or IST is outside its valid range.
C
   11 IER = -1
      RETURN
C
C Collinear nodes.
C
   12 IER = -2
      RETURN
C
   13 IER = 1
      RETURN
      END
      SUBROUTINE INTRC0 (N,PLAT,PLON,X,Y,Z,W,LIST,LPTR,
     .                   LEND, IST, PW,IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), IST, IER
      DOUBLE PRECISION    PLAT, PLON, X(N), Y(N), Z(N), W(N), PW
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/24/96
C
C   Given a triangulation of a set of nodes on the unit
C sphere, along with data values at the nodes, this sub-
C routine computes the value at a point P of a continuous
C function which interpolates the data values.  The interp-
C olatory function is linear on each underlying triangle
C (planar triangle with the same vertices as a spherical
C triangle).  If P is not contained in a triangle, an ex-
C trapolated value is taken to be the interpolated value at
C the nearest point of the triangulation boundary.
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       PLAT,PLON = Latitude and longitude of P in radians.
C
C       X,Y,Z = Arrays containing Cartesian coordinates of
C               the nodes.
C
C       W = Array containing data values at the nodes.  W(I)
C           is associated with (X(I),Y(I),Z(I)) for I =
C           1,...,N.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to STRIPACK
C                        Subroutine TRMESH.
C
C       IST = Index of the starting node in the search for a
C             triangle containing P.  1 .LE. IST .LE. N.
C             The output value of IST from a previous call
C             may be a good choice.
C
C Input parameters other than IST are not altered by this
C   routine.
C
C On output:
C
C       IST = Index of one of the vertices of the triangle
C             containing P (or nearest P) unless IER = -1
C             or IER = -2.
C
C       PW = Value of the interpolatory function at P if
C            IER .GE. 0.
C
C       IER = Error indicator:
C             IER = 0 if interpolation was performed
C                     successfully.
C             IER = 1 if extrapolation was performed
C                     successfully.
C             IER = -1 if N < 3 or IST is outside its valid
C                      range.
C             IER = -2 if the nodes are collinear.
C             IER = -3 if P is not in a triangle and the
C                      angle between P and the nearest boun-
C                      dary point is at least 90 degrees.
C
C STRIPACK modules required by INTRC0:  JRAND, LSTPTR,
C                                         STORE, TRFIND
C
C Intrinsic functions called by INTRC0:  COS, SIN
C
C***********************************************************
C
      INTEGER I1, I2, I3, LP, N1, N2
      DOUBLE PRECISION    B1, B2, B3, P(3), PTN1, PTN2, S12, SUM
C
C Local parameters:
C
C B1,B2,B3 = Barycentric coordinates of the central projec-
C              tion of P onto the underlying planar trian-
C              gle, or (B1 and B2) projection of Q onto the
C              underlying line segment N1-N2 when P is
C              exterior.  Unnormalized coordinates are
C              computed by TRFIND when P is in a triangle.
C I1,I2,I3 = Vertex indexes returned by TRFIND
C LP =       LIST pointer to N1 as a neighbor of N2 or N2
C              as a neighbor of N1
C N1,N2 =    Endpoints of a boundary arc which is visible
C              from P when P is not contained in a triangle
C P =        Cartesian coordinates of P
C PTN1 =     Scalar product (P,N1)
C PTN2 =     Scalar product (P,N2)
C S12 =      Scalar product (N1,N2)
C SUM =      Quantity used to normalize the barycentric
C              coordinates
C
      IF (N .LT. 3  .OR.  IST .LT. 1  .OR.  IST .GT. N)
     .    GO TO 11
C
C Transform (PLAT,PLON) to Cartesian coordinates.
C
      P(1) = COS(PLAT)*COS(PLON)
      P(2) = COS(PLAT)*SIN(PLON)
      P(3) = SIN(PLAT)
C
C Find the vertex indexes of a triangle containing P.
C
      CALL TRFIND(IST,P,N,X,Y,Z,LIST,LPTR,LEND, B1,B2,B3,
     .            I1,I2,I3)
      IF (I1 .EQ. 0) GO TO 12
      IST = I1
      IF (I3 .NE. 0) THEN
C
C P is contained in the triangle (I1,I2,I3).  Normalize the
C   barycentric coordinates.
C
        SUM = B1 + B2 + B3
        B1 = B1/SUM
        B2 = B2/SUM
        B3 = B3/SUM
        PW = B1*W(I1) + B2*W(I2) + B3*W(I3)
        IER = 0
        RETURN
      ENDIF
C
C P is exterior to the triangulation, and I1 and I2 are
C   boundary nodes which are visible from P.  Set PW to the
C   interpolated value at Q, where Q is the closest boundary
C   point to P.
C
C Traverse the boundary starting from the rightmost visible
C   node I1.
C
      N1 = I1
      PTN1 = P(1)*X(N1) + P(2)*Y(N1) + P(3)*Z(N1)
      IF (I1 .NE. I2) GO TO 2
C
C All boundary nodes are visible from P.  Find a boundary
C   arc N1->N2 such that P Left (N2 X N1)->N1.
C
C Counterclockwise boundary traversal:
C   Set N2 to the first neighbor of N1.
C
    1 LP = LEND(N1)
        LP = LPTR(LP)
        N2 = LIST(LP)
C
C Compute inner products (N1,N2) and (P,N2), and compute
C   B2 = DET(P,N1,N2 X N1).
C
        S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
        PTN2 = P(1)*X(N2) + P(2)*Y(N2) + P(3)*Z(N2)
        B2 = PTN2 - S12*PTN1
        IF (B2 .LE. 0.) GO TO 2
C
C P Right (N2 X N1)->N1 -- Iterate.
C
        N1 = N2
        I1 = N1
        PTN1 = PTN2
        GO TO 1
C
C P Left (N2 X N1)->N1, where N2 is the first neighbor of P1.
C   Clockwise boundary traversal:
C
    2 N2 = N1
        PTN2 = PTN1
C
C Set N1 to the last neighbor of N2 and test for
C   termination.
C
        LP = LEND(N2)
        N1 = -LIST(LP)
        IF (N1 .EQ. I1) GO TO 13
C
C Compute inner products (N1,N2) and (P,N1).
C
        S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
        PTN1 = P(1)*X(N1) + P(2)*Y(N1) + P(3)*Z(N1)
C
C Compute B2 = DET(P,N1,N2 X N1) = DET(Q,N1,N2 X N1)*(P,Q).
C
        B2 = PTN2 - S12*PTN1
        IF (B2 .LE. 0.) GO TO 2
C
C Compute B1 = DET(P,N2 X N1,N2) = DET(Q,N2 X N1,N2)*(P,Q).
C
      B1 = PTN1 - S12*PTN2
      IF (B1 .LE. 0.) THEN
C
C Q = N2.
C
        PW = W(N2)
      ELSE
C
C P Strictly Left (N2 X N1)->N2 and P Strictly Left
C   N1->(N2 X N1).  Thus Q lies on the interior of N1->N2.
C   Normalize the coordinates and compute PW.
C
        SUM = B1 + B2
        PW = (B1*W(N1) + B2*W(N2))/SUM
      ENDIF
      IER = 1
      RETURN
C
C N or IST is outside its valid range.
C
   11 IER = -1
      RETURN
C
C Collinear nodes.
C
   12 IER = -2
      RETURN
C
C The angular distance between P and the closest boundary
C   point to P is at least 90 degrees.
C
   13 IER = -3
      RETURN
      END
      SUBROUTINE INTRC1 (N,PLAT,PLON,X,Y,Z,F,LIST,LPTR,LEND,
     .                   IFLGS,SIGMA,IFLGG,GRAD, IST, FP,
     .                   IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), IFLGS, IFLGG,
     .        IST, IER
      DOUBLE PRECISION    PLAT, PLON, X(N), Y(N), Z(N), F(N), SIGMA(*),
     .        GRAD(3,N), FP
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/25/96
C
C   Given a triangulation of a set of nodes on the unit
C sphere, along with data values and gradients at the nodes,
C this routine computes a value F(P), where F interpolates
C the nodal data and is once-continuously differentiable
C over the convex hull of the nodes.  Refer to Function FVAL
C for further details.  If P is not contained in a triangle,
C an extrapolated value is computed by extending F beyond
C the boundary in a continuous fashion.
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 3
C           and N .GE. 7 if IFLGG .LE. 0.
C
C       PLAT,PLON = Latitude and longitude in radians of the
C                   point P at which F is to be evaluated.
C
C       X,Y,Z = Arrays of length N containing Cartesian
C               coordinates of the nodes.
C
C       F = Array of length N containing values of F at the
C           nodes:  F(I) = F(X(I),Y(I),Z(I)).
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to STRIPACK
C                        Subroutine TRMESH.
C
C       IFLGS = Tension factor option:
C               IFLGS .LE. 0 if a single uniform tension
C                            factor is to be used.
C               IFLGS .GE. 1 if variable tension is desired.
C
C       SIGMA = Uniform tension factor (IFLGS .LE. 0), or
C               array containing tension factors associated
C               with arcs in one-to-one correspondence with
C               LIST entries (IFLGS .GE. 1).  Refer to Sub-
C               programs FVAL, GETSIG, SIG0, SIG1, and SIG2.
C
C       IFLGG = Gradient option:
C               IFLGG .LE. 0 if INTRC1 is to provide grad-
C                            ient estimates as needed (from
C                            GRADL).
C               IFLGG .GE. 1 if gradients are user-provided
C                            in GRAD.  This is more effici-
C                            ent if INTRC1 is to be called
C                            several times.
C
C       GRAD = 3 by N array whose I-th column contains
C              an estimated gradient at node I if IFLGG .GE.
C              1, or unused dummy parameter if IFLGG .LE. 0.
C              Refer to Subroutines GRADL and GRADG.
C
C       IST = Index of the starting node in the search for a
C             triangle containing P.  The output value of
C             IST from a previous call may be a good choice.
C             1 .LE. IST .LE. N.
C
C Input parameters other than IST are not altered by this
C   routine.
C
C On output:
C
C       IST = Index of one of the vertices of the triangle
C             containing P (or a boundary node if P is not
C             contained in a triangle) unless IER = -1 or
C             IER = -2.
C
C       FP = Value of F at P unless IER < 0, in which case
C            FP is not defined.
C
C       IER = Error indicator and information flag:
C             IER = 0 if no errors were encountered and P is
C                     contained in a triangle.
C             IER = 1 if no errors were encountered and
C                     extrapolation was required.
C             IER = -1 if N or IST is outside its valid
C                      range.
C             IER = -2 if the nodes are collinear.
C             IER = -3 if the angular distance between P and
C                      the nearest point of the triangula-
C                      tion is at least 90 degrees.
C
C STRIPACK modules required by INTRC1: JRAND, LSTPTR, STORE,
C                                        TRFIND
C                    (and optionally)  GETNP if IFLGG .LE. 0
C
C SSRFPACK modules required by INTRC1:  ARCINT, ARCLEN,
C                                         FVAL, HVAL, SNHCSH
C              (and if IFLGG .LE. 0)  APLYR, APLYRT, CONSTR,
C                                       GIVENS, GRADL,
C                                       ROTATE, SETUP
C
C Intrinsic functions called by INTRC1:  COS, SIN, SQRT
C
C***********************************************************
C
      INTEGER LSTPTR
      DOUBLE PRECISION    ARCLEN, FVAL
      INTEGER I, IERR, I1, I2, I3, LP, N1, N2, NN
      DOUBLE PRECISION    A, B1, B2, B3, DUM(3), FQ, GQ(3), GQN, G1(3),
     .        G2(3), G3(3), P(3), P1(3), P2(3), P3(3), PTGQ,
     .        PTN1, PTN2, Q(3), QNORM, S1, S2, S3, S12, SUM
C
C Local parameters:
C
C A =        Angular separation between P and Q
C B1,B2,B3 = Barycentric coordinates of the central projec-
C              tion of P onto the underlying planar triangle,
C              or (B1 and B2) projection of Q onto the
C              underlying line segment N1-N2 when P is
C              exterior.  Unnormalized coordinates are
C              computed by TRFIND when P is in a triangle.
C DUM =      Dummy parameter for ARCINT
C FQ,GQ =    Interpolated value and gradient at Q
C GQN =      Negative of the component of GQ in the direction
C              Q->P
C G1,G2,G3 = Gradients at I1, I2, and I3, or (G1 and G2) at
C              N1 and N2
C I =        DO-loop index
C IERR =     Error flag for calls to GRADL
C I1,I2,I3 = Vertex indexes returned by TRFIND
C LP =       LIST pointer
C N1,N2 =    Indexes of the endpoints of a boundary arc when
C              P is exterior (not contained in a triangle)
C NN =       Local copy of N
C P =        Cartesian coordinates of P
C P1,P2,P3 = Cartesian coordinates of the vertices I1, I2,
C              and I3, or (P1 and P2) coordinates of N1 and
C              N2 if P is exterior
C PTGQ =     Scalar product (P,GQ) -- factor of the component
C              of GQ in the direction Q->P
C PTN1 =     Scalar product (P,N1) -- factor of B1 and B2
C PTN2 =     Scalar product (P,N2) -- factor of B1 and B2
C Q =        Closest boundary point to P when P is exterior
C QNORM =    Factor used to normalize Q
C S1,S2,S3 = Tension factors associated with the triangle
C              sides opposite I1, I2, and I3, or (S1) the
C              boundary arc N1-N2
C S12 =      Scalar product (N1,N2) -- factor of B1 and B2
C SUM =      Quantity used to normalize the barycentric
C              coordinates
C
      NN = N
      IF (NN .LT. 3  .OR.  (IFLGG .LE. 0  .AND.  NN .LT. 7)
     .    .OR.  IST .LT. 1  .OR.  IST .GT. NN) GO TO 11
C
C Transform (PLAT,PLON) to Cartesian coordinates.
C
      P(1) = COS(PLAT)*COS(PLON)
      P(2) = COS(PLAT)*SIN(PLON)
      P(3) = SIN(PLAT)
C
C Locate P with respect to the triangulation.
C
      CALL TRFIND (IST,P,NN,X,Y,Z,LIST,LPTR,LEND, B1,B2,B3,
     .             I1,I2,I3)
      IF (I1 .EQ. 0) GO TO 12
      IST = I1
      IF (I3 .NE. 0) THEN
C
C P is contained in the triangle (I1,I2,I3).  Store the
C   vertex coordinates, gradients, and tension factors in
C   local variables.
C
        P1(1) = X(I1)
        P1(2) = Y(I1)
        P1(3) = Z(I1)
        P2(1) = X(I2)
        P2(2) = Y(I2)
        P2(3) = Z(I2)
        P3(1) = X(I3)
        P3(2) = Y(I3)
        P3(3) = Z(I3)
        IF (IFLGG .GT. 0) THEN
C
C   Gradients are user-provided.
C
          DO 1 I = 1,3
            G1(I) = GRAD(I,I1)
            G2(I) = GRAD(I,I2)
            G3(I) = GRAD(I,I3)
    1       CONTINUE
        ELSE
C
C   Compute gradient estimates at the vertices.
C
          CALL GRADL (NN,I1,X,Y,Z,F,LIST,LPTR,LEND, G1,IERR)
          IF (IERR .LT. 0) GO TO 12
          CALL GRADL (NN,I2,X,Y,Z,F,LIST,LPTR,LEND, G2,IERR)
          IF (IERR .LT. 0) GO TO 12
          CALL GRADL (NN,I3,X,Y,Z,F,LIST,LPTR,LEND, G3,IERR)
          IF (IERR .LT. 0) GO TO 12
        ENDIF
C
        IF (IFLGS .GT. 0) THEN
C
C   Variable tension:
C
          LP = LSTPTR(LEND(I2),I3,LIST,LPTR)
          S1 = SIGMA(LP)
          LP = LSTPTR(LEND(I3),I1,LIST,LPTR)
          S2 = SIGMA(LP)
          LP = LSTPTR(LEND(I1),I2,LIST,LPTR)
          S3 = SIGMA(LP)
        ELSE
C
C   Uniform tension:
C
          S1 = SIGMA(1)
          S2 = S1
          S3 = S1
        ENDIF
C
C Normalize the coordinates.
C
        SUM = B1 + B2 + B3
        B1 = B1/SUM
        B2 = B2/SUM
        B3 = B3/SUM
        FP = FVAL(B1,B2,B3,P1,P2,P3,F(I1),F(I2),F(I3),G1,
     .            G2,G3,S1,S2,S3)
        IER = 0
        RETURN
      ENDIF
C
C P is exterior to the triangulation, and I1 and I2 are
C   boundary nodes which are visible from P.  Extrapolate to
C   P by linear (with respect to arc-length) interpolation
C   of the value and directional derivative (gradient comp-
C   onent in the direction Q->P) of the interpolatory
C   surface at Q where Q is the closest boundary point to P.
C
C Determine Q by traversing the boundary starting from I1.
C
      N1 = I1
      PTN1 = P(1)*X(N1) + P(2)*Y(N1) + P(3)*Z(N1)
      IF (I1 .NE. I2) GO TO 3
C
C All boundary nodes are visible from P.  Find a boundary
C   arc N1->N2 such that P Left (N2 X N1)->N1.
C
C Counterclockwise boundary traversal:
C   Set N2 to the first neighbor of N1.
C
    2 LP = LEND(N1)
        LP = LPTR(LP)
        N2 = LIST(LP)
C
C Compute inner products (N1,N2) and (P,N2), and compute
C   B2 = Det(P,N1,N2 X N1).
C
        S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
        PTN2 = P(1)*X(N2) + P(2)*Y(N2) + P(3)*Z(N2)
        B2 = PTN2 - S12*PTN1
        IF (B2 .LE. 0.) GO TO 3
C
C P Right (N2 X N1)->N1:  iterate.
C
        N1 = N2
        I1 = N1
        PTN1 = PTN2
        GO TO 2
C
C P Left (N2 X N1)->N1 where N2 is the first neighbor of N1.
C   Clockwise boundary traversal:
C
    3 N2 = N1
        PTN2 = PTN1
C
C Set N1 to the last neighbor of N2 and test for
C   termination.
C
        LP = LEND(N2)
        N1 = -LIST(LP)
        IF (N1 .EQ. I1) GO TO 13
C
C Compute inner products (N1,N2) and (P,N1).
C
        S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
        PTN1 = P(1)*X(N1) + P(2)*Y(N1) + P(3)*Z(N1)
C
C Compute B2 = Det(P,N1,N2 X N1) = Det(Q,N1,N2 X N1)*(P,Q).
C
        B2 = PTN2 - S12*PTN1
        IF (B2 .LE. 0.) GO TO 3
C
C Compute B1 = Det(P,N2 X N1,N2) = Det(Q,N2 X N1,N2)*(P,Q).
C
      B1 = PTN1 - S12*PTN2
      IF (B1 .LE. 0.) THEN
C
C Q = N2.  Store value, coordinates, and gradient at Q.
C
        FQ = F(N2)
        Q(1) = X(N2)
        Q(2) = Y(N2)
        Q(3) = Z(N2)
        IF (IFLGG .GT. 0) THEN
          DO 4 I = 1,3
            GQ(I) = GRAD(I,N2)
    4       CONTINUE
        ELSE
          CALL GRADL (NN,N2,X,Y,Z,F,LIST,LPTR,LEND, GQ,IERR)
          IF (IERR .LT. 0) GO TO 12
        ENDIF
C
C Extrapolate to P:  FP = FQ + A*(GQ,Q X (PXQ)/SIN(A)),
C   where A is the angular separation between Q and P,
C   and Sin(A) is the magnitude of P X Q.
C
        A = ARCLEN(Q,P)
        PTGQ = P(1)*GQ(1) + P(2)*GQ(2) + P(3)*GQ(3)
        FP = FQ
        IF (A .NE. 0.) FP = FP + PTGQ*A/SIN(A)
        IER = 1
        RETURN
      ENDIF
C
C P Strictly Left (N2 X N1)->N2 and P Strictly Left
C   N1->(N2 X N1).  Thus Q lies on the interior of N1->N2.
C   Store coordinates of N1 and N2 in local variables.
C
      P1(1) = X(N1)
      P1(2) = Y(N1)
      P1(3) = Z(N1)
      P2(1) = X(N2)
      P2(2) = Y(N2)
      P2(3) = Z(N2)
C
C Compute the central projection of Q onto P2-P1 and
C   normalize to obtain Q.
C
      QNORM = 0.
      DO 5 I = 1,3
        Q(I) = B1*P1(I) + B2*P2(I)
        QNORM = QNORM + Q(I)*Q(I)
    5   CONTINUE
      QNORM = SQRT(QNORM)
      DO 6 I = 1,3
        Q(I) = Q(I)/QNORM
    6   CONTINUE
C
C Store gradients at N1 and N2 and tension factor S1.
C
      IF (IFLGG .GT. 0) THEN
        DO 7 I = 1,3
          G1(I) = GRAD(I,N1)
          G2(I) = GRAD(I,N2)
    7     CONTINUE
      ELSE
        CALL GRADL (NN,N1,X,Y,Z,F,LIST,LPTR,LEND, G1,IERR)
        IF (IERR .LT. 0) GO TO 12
        CALL GRADL (NN,N2,X,Y,Z,F,LIST,LPTR,LEND, G2,IERR)
        IF (IERR .LT. 0) GO TO 12
      ENDIF
C
      IF (IFLGS .LE. 0) S1 = SIGMA(1)
      IF (IFLGS .GE. 1) S1 = SIGMA(LP)
C
C Compute an interpolated value and normal gradient
C   component at Q.
C
      CALL ARCINT (Q,P1,P2,F(N1),F(N2),G1,G2,S1, FQ,DUM,GQN)
C
C Extrapolate to P:  the normal gradient component GQN is
C   the negative of the component in the direction Q->P.
C
      FP = FQ - GQN*ARCLEN(Q,P)
      IER = 1
      RETURN
C
C N or IST is outside its valid range.
C
   11 IER = -1
      RETURN
C
C Collinear nodes encountered.
C
   12 IER = -2
      RETURN
C
C The distance between P and the closest boundary point
C   is at least 90 degrees.
C
   13 IER = -3
      RETURN
      END
      SUBROUTINE ROTATE (N,C,S, X,Y )
      INTEGER N
      DOUBLE PRECISION    C, S, X(N), Y(N)
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/01/88
C
C                                                ( C  S)
C   This subroutine applies the Givens rotation  (     )  to
C                                                (-S  C)
C                    (X(1) ... X(N))
C the 2 by N matrix  (             ) .
C                    (Y(1) ... Y(N))
C
C   This routine is identical to Subroutine SROT from the
C LINPACK BLAS (Basic Linear Algebra Subroutines).
C
C On input:
C
C       N = Number of columns to be rotated.
C
C       C,S = Elements of the Givens rotation.  Refer to
C             Subroutine GIVENS.
C
C The above parameters are not altered by this routine.
C
C       X,Y = Arrays of length .GE. N containing the compo-
C             nents of the vectors to be rotated.
C
C On output:
C
C       X,Y = Arrays containing the rotated vectors (not
C             altered if N < 1).
C
C Modules required by ROTATE:  None
C
C***********************************************************
C
      INTEGER I
      DOUBLE PRECISION    XI, YI
C
      DO 1 I = 1,N
        XI = X(I)
        YI = Y(I)
        X(I) = C*XI + S*YI
        Y(I) = -S*XI + C*YI
    1   CONTINUE
      RETURN
      END
      SUBROUTINE SETUP (XI,YI,WI,WK,S1,S2,WT, ROW)
      DOUBLE PRECISION XI, YI, WI, WK, S1, S2, WT, ROW(6)
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   05/09/92
C
C   This subroutine sets up the I-th row of an augmented
C regression matrix for a weighted least squares fit of a
C quadratic function Q(X,Y) to a set of data values Wi,
C where Q(0,0) = Wk.  The first 3 columns (quadratic terms)
C are scaled by 1/S2 and the fourth and fifth columns (lin-
C ear terms) are scaled by 1/S1.
C
C On input:
C
C       XI,YI = Coordinates of node I.
C
C       WI = Data value at node I.
C
C       WK = Data value interpolated by Q at the origin.
C
C       S1,S2 = Inverse scale factors.
C
C       WT = Weight factor corresponding to the I-th
C            equation.
C
C       ROW = Array of length 6.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       ROW = Array containing a row of the augmented re-
C             gression matrix.
C
C Modules required by SETUP:  None
C
C***********************************************************
C
      DOUBLE PRECISION W1, W2
C
C Local parameters:
C
C W1 = Weighted scale factor for the linear terms
C W2 = Weighted scale factor for the quadratic terms
C
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
      SUBROUTINE SGPRNT (N,LUNIT,LIST,LPTR,LEND,SIGMA)
      INTEGER N, LUNIT, LIST(*), LPTR(*), LEND(N)
      DOUBLE PRECISION    SIGMA(*)
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/21/98
C
C   Given a triangulation of a set of nodes on the unit
C sphere, along with an array of tension factors associated
C with the triangulation arcs, this subroutine prints the
C lst of arcs (with tension factors) ordered by endpoint
C nodal indexes.  An arc is identified with its smaller
C endpoint index:  N1-N2, where N1 < N2.
C
C   This routine is identical to the similarly named routine
C in SRFPACK.
C
C
C On input:
C
C       N = Number of nodes in the triangulation.  3 .LE. N
C           .LE. 9999.
C
C       LUNIT = Logical unit for output.  0 .LE. LUNIT .LE.
C               99.  Output is printed on unit 6 if LUNIT is
C               outside its valid range.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to STRIPACK
C                        Subroutine TRMESH.
C
C       SIGMA = Array of length 2*NA = 6*(N-1)-2*NB, where
C               NA and NB are the numbers of arcs and boun-
C               dary nodes, respectively, containing tension
C               factors associated with arcs in one-to-one
C               correspondence with LIST entries.  Note that
C               each arc N1-N2 has two LIST entries and
C               thus, SIGMA(I) and SIGMA(J) should be iden-
C               tical, where LIST(I) = N2 (in the adjacency
C               lst for N1) and LIST(J) = N1 (in the lst
C               associated with N2).  Both SIGMA(I) and
C               SIGMA(J) are printed if they are not iden-
C               tical.
C
C None of the parameters are altered by this routine.
C
C STRIPACK module required by SGPRNT:  LSTPTR
C
C Intrinsic function called by SGPRNT:  ABS
C
C***********************************************************
C
      INTEGER LSTPTR
      INTEGER LP1, LP2, LPL, LUN, N1, N2, NA, NAT, NB, NE,
     .        NL, NLMAX, NM1, NMAX
      LOGICAL ERROR
      DOUBLE PRECISION    SIG
      DATA NMAX/9999/,  NLMAX/58/
C
      LUN = LUNIT
      IF (LUN .LT. 0  .OR.  LUN .GT. 99) LUN = 6
C
C Print a heading, test for invalid N, and initialize coun-
C   ters:
C
C NL = Number of lines printed on the current page
C NA = Number of arcs encountered
C NE = Number of errors in SIGMA encountered
C NB = Number of boundary nodes encountered
C
      WRITE (LUN,100) N
      IF (N .LT. 3  .OR.  N .GT. NMAX) GO TO 4
      NL = 6
      NA = 0
      NE = 0
      NB = 0
C
C Outer loop on nodes N1.  LPL points to the last neighbor
C   of N1.
C
      NM1 = N - 1
      DO 3 N1 = 1,NM1
        LPL = LEND(N1)
        IF (LIST(LPL) .LT. 0) NB = NB + 1
        LP1 = LPL
C
C Inner loop on neighbors N2 of N1 such that N1 < N2.
C
    1   LP1 = LPTR(LP1)
          N2 = ABS(LIST(LP1))
          IF (N2 .LT. N1) GO TO 2
          NA = NA + 1
          SIG = SIGMA(LP1)
C
C   Test for an invalid SIGMA entry.
C
          LP2 = LSTPTR (LEND(N2),N1,LIST,LPTR)
          ERROR = SIGMA(LP2) .NE. SIG
          IF (ERROR) NE = NE + 1
C
C   Print a line and update the counters.
C
          IF (.NOT. ERROR) WRITE (LUN,110) N1, N2, SIG
          IF (ERROR) WRITE (LUN,120) N1, N2, SIG, SIGMA(LP2)
          NL = NL + 1
          IF (NL .GE. NLMAX) THEN
            WRITE (LUN,130)
            NL = 1
          ENDIF
C
C Bottom of loop on neighbors N2 of N1.
C
    2     IF (LP1 .NE. LPL) GO TO 1
    3   CONTINUE
      LPL = LEND(N)
      IF (LIST(LPL) .LT. 0) NB = NB + 1
C
C Test for errors in SIGMA.
C
      IF (NE .GT. 0) WRITE (LUN,200) NE
C
C Print NA and test for an invalid triangulation.
C
      WRITE (LUN,140) NA
      IF (NB .NE. 0) THEN
        NAT = 3*NM1 - NB
      ELSE
        NAT = 3*N - 6
      ENDIF
      IF (NAT .NE. NA) WRITE (LUN,210) NAT
      RETURN
C
C N is outside its valid range.
C
    4 WRITE (LUN,220) NMAX
      RETURN
C
C Print formats:
C
  100 FORMAT ('1',14X,'TENSION FACTORS,  N =',I5,
     .        ' NODES'//1X,18X,'N1',5X,'N2',8X,'TENSION'//)
  110 FORMAT (1X,16X,I4,3X,I4,5X,F12.8)
  120 FORMAT (1X,16X,I4,3X,I4,5X,F12.8,3X,F12.8,' *')
  130 FORMAT ('1')
  140 FORMAT (//1X,10X,'NA =',I5,' ARCS')
C
C Error messages:
C
  200 FORMAT (//1X,10X,'*',I5,' ERRORS IN SIGMA')
  210 FORMAT (/1X,10X,'*** ERROR IN TRIANGULATION -- ',
     .        '3N-NB-3 = ',I5,' ***')
  220 FORMAT (1X,10X,'*** N IS OUT OF RANGE -- NMAX = ',
     .        I4,' ***')
      END
      DOUBLE PRECISION FUNCTION SIG0 (N1,N2,N,X,Y,Z,H,LIST,LPTR,LEND,
     .                    GRAD,IFLGB,HBND,TOL,
     .                    IFLGS, SIGMA, IER)
      INTEGER N1, N2, N, LIST(*), LPTR(*), LEND(N), IFLGB,
     .        IFLGS, IER
      DOUBLE PRECISION    X(N), Y(N), Z(N), H(N), GRAD(3,N), HBND, TOL,
     .        SIGMA(*)
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/21/96
C
C   Given a triangulation of a set of nodes on the unit
C sphere, along with data values H and gradients GRAD at the
C nodes, this function determines the smallest tension fac-
C tor SIG0 such that the Hermite interpolatory tension
C spline H(A), defined by SIG0 and the endpoint values and
C directional derivatives associated with an arc N1-N2, is
C bounded (either above or below) by HBND for all A in
C (A1,A2), where (A1,A2) denotes an interval corresponding
C to the arc and A is the arc-length.
C
C On input:
C
C       N1,N2 = Nodal indexes of the endpoints of an arc for
C               which the tension factor is to be computed.
C               The indexes must be distinct and lie in the
C               range 1 to N, and if IFLGS .GE. 1, they must
C               correspond to adjacent nodes in the triangu-
C               lation.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y,Z = Arrays of length N containing coordinates of
C               the nodes.  X(I)**2 + Y(I)**2 + Z(I)**2 = 1.
C
C       H = Array of length N containing data values at the
C           nodes.  H(I) is associated with (X(I),Y(I),Z(I))
C           for I = 1 to N.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to STRIPACK
C                        Subroutine TRMESH.
C
C       GRAD = Array dimensioned 3 by N whose columns con-
C              tain gradients at the nodes.  GRAD( ,J) must
C              be orthogonal to node J:  GRAD(1,J)*X(J) +
C              GRAD(2,J)*Y(J) + GRAD(3,J)*Z(J) = 0.  Refer
C              to Subroutines GRADG, GRADL, and SMSURF.
C
C       IFLGB = Bound option indicator:
C               IFLGB = -1 if HBND is a lower bound on H.
C               IFLGB = 1 if HBND is an upper bound on H.
C
C       HBND = Bound on H.  HBND .LE. min(H1,H2) if IFLGB =
C              -1 and HBND .GE. max(H1,H2) if IFLGB = 1,
C              where H1 and H2 are the data values at the
C              endpoints of the arc N1-N2.
C
C       TOL = Tolerance whose magnitude determines how close
C             SIG0 is to its optimal value when nonzero
C             finite tension is necessary and sufficient to
C             satisfy the constraint.  For a lower bound,
C             SIG0 is chosen so that HBND .LE. HMIN .LE.
C             HBND + abs(TOL), where HMIN is the minimum
C             value of H on the arc, and for an upper bound,
C             the maximum of H satisfies HBND - abs(TOL)
C             .LE. HMAX .LE. HBND.  Thus, the constraint is
C             satisfied but possibly with more tension than
C             necessary.
C
C       IFLGS = Tension array option indicator:
C               IFLGS .LE. 0 if SIGMA is not to be used.
C               IFLGS .GE. 1 if SIGMA is to be updated by
C                            storing SIG0 in the appropriate
C                            locations.
C
C The above parameters are not altered by this function.
C
C       SIGMA = Dummy parameter (IFLGS .LE. 0) or array con-
C               taining tension factors associated with arcs
C               in one-to-one correspondence with LIST
C               entries (IFLGS .GE. 1).  Refer to Subroutine
C               GETSIG.
C
C On output:
C
C       SIGMA = Tension factor array updated with the new
C               value if and only if IFLGS .GE. 1 and IER
C               .GE. 0.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered and the
C                     constraint can be satisfied with fin-
C                     ite tension.
C             IER = 1 if no errors were encountered but in-
C                     finite tension is required to satisfy
C                     the constraint (e.g., IFLGB = -1, HBND
C                     = H(A1), and the directional deriva-
C                     tive of H at A1 is negative).
C             IER = -1 if N1, N2, N, or IFLGB is outside its
C                      valid range.
C             IER = -2 if nodes N1 and N2 coincide or IFLGS
C                      .GE. 1 and the nodes are not adja-
C                      cent.
C             IER = -3 if HBND is outside its valid range.
C
C       SIG0 = Minimum tension factor defined above unless
C              IER < 0, in which case SIG0 = -1.  If IER
C              = 1, SIG0 is set to 85, resulting in an
C              approximation to the linear interpolant of
C              the endpoint values.
C
C STRIPACK module required by SIG0:  STORE
C
C SSRFPACK modules required by SIG0:  ARCLEN, SNHCSH
C
C Intrinsic functions called by SIG0:  ABS, EXP, LOG, MAX,
C                                        MIN, REAL, SIGN,
C                                        SQRT
C
C***********************************************************
C
      DOUBLE PRECISION    ARCLEN, STORE
      INTEGER LP1, LP2, LPL, LUN, NIT
      DOUBLE PRECISION    A, A0, AA, AL, B, B0, BND, C, C1, C2, COSHM,
     .        COSHMM, D, D0, D1PD2, D2, DMAX, DSIG, E, EMS,
     .        F, F0, FMAX, FNEG, FTOL, H1, H2, P1(3), P2(3),
     .        R, RF, RSIG, RTOL, S, S1, S2, SBIG, SCM, SIG,
     .        SINHM, SNEG, SSINH, SSM, STOL, T, T0, T1, T2,
     .        TM, UN(3), UNORM
C
      DATA SBIG/85./,  LUN/-1/
      RF = REAL(IFLGB)
      BND = HBND
C
C Print a heading.
C
      IF (LUN .GE. 0  .AND.  RF .LT. 0.) WRITE (LUN,100) N1,
     .                                   N2, BND
      IF (LUN .GE. 0  .AND.  RF .GT. 0.) WRITE (LUN,110) N1,
     .                                   N2, BND
  100 FORMAT (//1X,'SIG0 -- N1 =',I4,', N2 =',I4,
     .        ', LOWER BOUND = ',E15.8)
  110 FORMAT (//1X,'SIG0 -- N1 =',I4,', N2 =',I4,
     .        ', UPPER BOUND = ',E15.8)
C
C Test for errors and store local parameters.
C
      IER = -1
      IF (MIN(N1,N2) .LT. 1  .OR.  N1 .EQ. N2  .OR.
     .    MAX(N1,N2,3) .GT. N  .OR.  ABS(RF) .NE. 1.)
     .   GO TO 11
      IER = -2
      IF (IFLGS .GT. 0) THEN
C
C   Set LP1 and LP2 to the pointers to N2 as a neighbor of
C     N1 and N1 as a neighbor of N2, respectively.
C
        LPL = LEND(N1)
        LP1 = LPTR(LPL)
    1   IF (LIST(LP1) .EQ. N2) GO TO 2
          LP1 = LPTR(LP1)
          IF (LP1 .NE. LPL) GO TO 1
        IF (ABS(LIST(LP1)) .NE. N2) GO TO 11
C
    2   LPL = LEND(N2)
        LP2 = LPTR(LPL)
    3   IF (LIST(LP2) .EQ. N1) GO TO 4
          LP2 = LPTR(LP2)
          IF (LP2 .NE. LPL) GO TO 3
        IF (ABS(LIST(LP2)) .NE. N1) GO TO 11
      ENDIF
C
C Store nodal coordinates P1 and P2, compute arc-length AL
C   and unit normal UN = (P1 X P2)/UNORM, and test for
C   coincident nodes.
C
    4 P1(1) = X(N1)
      P1(2) = Y(N1)
      P1(3) = Z(N1)
      P2(1) = X(N2)
      P2(2) = Y(N2)
      P2(3) = Z(N2)
      AL = ARCLEN(P1,P2)
      UN(1) = P1(2)*P2(3) - P1(3)*P2(2)
      UN(2) = P1(3)*P2(1) - P1(1)*P2(3)
      UN(3) = P1(1)*P2(2) - P1(2)*P2(1)
      UNORM = SQRT(UN(1)*UN(1) + UN(2)*UN(2) + UN(3)*UN(3))
      IF (UNORM .EQ. 0.  .OR.  AL .EQ. 0.) GO TO 11
C
C Store endpoint data values and test for valid constraint.
C
      H1 = H(N1)
      H2 = H(N2)
      IER = -3
      IF ((RF .LT. 0.  .AND.  MIN(H1,H2) .LT. BND)  .OR.
     .    (RF .GT. 0.  .AND.  BND .LT. MAX(H1,H2)))
     .   GO TO 11
C
C Compute scaled directional derivatives S1,S2 at the end-
C   points (for the direction N1->N2) and test for infinite
C   tension required.
C
      S1 = AL*(GRAD(1,N1)*P2(1) + GRAD(2,N1)*P2(2) +
     .         GRAD(3,N1)*P2(3))/UNORM
      S2 = -AL*(GRAD(1,N2)*P1(1) + GRAD(2,N2)*P1(2) +
     .          GRAD(3,N2)*P1(3))/UNORM
      IER = 1
      SIG = SBIG
      IF ((H1 .EQ. BND  .AND.  RF*S1 .GT. 0.)  .OR.
     .    (H2 .EQ. BND  .AND.  RF*S2 .LT. 0.)) GO TO 10
C
C Test for SIG = 0 sufficient.
C
      IER = 0
      SIG = 0.
      IF (RF*S1 .LE. 0.  .AND.  RF*S2 .GE. 0.) GO TO 10
C
C   Compute first difference S and coefficients A0 and B0
C     of the Hermite cubic interpolant H0(A) = H2 - (S2*R +
C     B0*R**2 + (A0/3)*R**3), where R(A) = (A2-A)/AL.
C
      S = H2 - H1
      T0 = 3.*S - S1 - S2
      A0 = 3.*(S-T0)
      B0 = T0 - S2
      D0 = T0*T0 - S1*S2
C
C   H0 has local extrema in (A1,A2) iff S1*S2 < 0 or
C     (T0*(S1+S2) < 0 and D0 .GE. 0).
C
      IF (S1*S2 .GE. 0.  .AND.  (T0*(S1+S2) .GE. 0.  .OR.
     .    D0 .LT. 0.)) GO TO 10
      IF (A0 .EQ. 0.) THEN
C
C   H0 is quadratic and has an extremum at R = -S2/(2*B0).
C     H0(R) = H2 + S2**2/(4*B0).  Note that A0 = 0 implies
C     2*B0 = S1-S2, and S1*S2 < 0 implies B0 .NE. 0.
C     Also, the extremum is a min iff HBND is a lower bound.
C
        F0 = (BND - H2 - S2*S2/(4.*B0))*RF
      ELSE
C
C   A0 .NE. 0 and H0 has extrema at R = (-B0 +/- SQRT(D0))/
C     A0 = S2/(-B0 -/+ SQRT(D0)), where the negative root
C     corresponds to a min.  The expression for R is chosen
C     to avoid cancellation error.  H0(R) = H2 + (S2*B0 +
C     2*D0*R)/(3*A0).
C
        T = -B0 - SIGN(SQRT(D0),B0)
        R = T/A0
        IF (RF*B0 .GT. 0.) R = S2/T
        F0 = (BND - H2 - (S2*B0+2.*D0*R)/(3.*A0))*RF
      ENDIF
C
C   F0 .GE. 0 iff SIG = 0 is sufficient to satisfy the
C     constraint.
C
      IF (F0 .GE. 0.) GO TO 10
C
C Find a zero of F(SIG) = (BND-H(R))*RF where the derivative
C   of H, HP, vanishes at R.  F is a nondecreasing function,
C   F(0) < 0, and F = FMAX for SIG sufficiently large.
C
C Initialize parameters for the secant method.  The method
C   uses three points:  (SG0,F0), (SIG,F), and (SNEG,FNEG),
C   where SG0 and SNEG are defined implicitly by DSIG = SIG
C   - SG0 and DMAX = SIG - SNEG.  SG0 is initially zero and
C   SNEG is initialized to a sufficiently large value that
C   FNEG > 0.  This value is used only if the initial value
C   of F is negative.
C
      FMAX = MAX(1.E-3,MIN(ABS(H1-BND),ABS(H2-BND)))
      T = MAX(ABS(H1-BND),ABS(H2-BND))
      SIG = MAX(ABS(S1),ABS(S2))/T
      DMAX = SIG*(1.-T/FMAX)
      SNEG = SIG - DMAX
      IF (LUN .GE. 0) WRITE (LUN,120) SIG, SNEG, F0, FMAX
  120 FORMAT (1X,8X,'SIG = ',E15.8,', SNEG = ',E15.8/
     .        1X,9X,'F0 = ',E15.8,', FMAX = ',E15.8/)
      DSIG = SIG
      FNEG = FMAX
      D2 = S2 - S
      D1PD2 = S2 - S1
      NIT = 0
C
C Compute an absolute tolerance FTOL = abs(TOL) and a
C   relative tolerance RTOL = 100*Macheps.
C
      FTOL = ABS(TOL)
      RTOL = 1.
    5 RTOL = RTOL/2.
        IF (STORE(RTOL+1.) .GT. 1.) GO TO 5
      RTOL = RTOL*200.
C
C Top of loop:  compute F.
C
    6 EMS = EXP(-SIG)
      IF (SIG .LE. .5) THEN
C
C   Use approximations designed to avoid cancellation error
C     (associated with small SIG) in the modified hyperbolic
C     functions.
C
        CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
        C1 = SIG*COSHM*D2 - SINHM*D1PD2
        C2 = SIG*(SINHM+SIG)*D2 - COSHM*D1PD2
        A = C2 - C1
        AA = A/EMS
        E = SIG*SINHM - COSHMM - COSHMM
      ELSE
C
C   Scale SINHM and COSHM by 2*exp(-SIG) in order to avoid
C     overflow.
C
        TM = 1. - EMS
        SSINH = TM*(1.+EMS)
        SSM = SSINH - 2.*SIG*EMS
        SCM = TM*TM
        C1 = SIG*SCM*D2 - SSM*D1PD2
        C2 = SIG*SSINH*D2 - SCM*D1PD2
        AA = 2.*(SIG*TM*D2 + (TM-SIG)*D1PD2)
        A = EMS*AA
        E = SIG*SSINH - SCM - SCM
      ENDIF
C
C   HP(R) = (S2 - (C1*sinh(SIG*R) - C2*coshm(SIG*R))/E)/DT
C     = 0 for ESR = (-B +/- sqrt(D))/A = C/(-B -/+ sqrt(D))
C     where ESR = exp(SIG*R), A = C2-C1, D = B**2 - A*C, and
C     B and C are defined below.
C
      B = E*S2 - C2
      C = C2 + C1
      D = B*B - A*C
      F = 0.
      IF (AA*C .EQ. 0.  .AND.  B .EQ. 0.) GO TO 7
      F = FMAX
      IF (D .LT. 0.) GO TO 7
      T1 = SQRT(D)
      T = -B - SIGN(T1,B)
      RSIG = 0.
      IF (RF*B .LT. 0.  .AND.  AA .NE. 0.) THEN
        IF (T/AA .GT. 0.) RSIG = SIG + LOG(T/AA)
      ENDIF
      IF ((RF*B .GT. 0.  .OR.  AA .EQ. 0.)  .AND.
     .    C/T .GT. 0.) RSIG = LOG(C/T)
      IF ((RSIG .LE. 0.  .OR.  RSIG .GE. SIG)  .AND.
     .    B .NE. 0.) GO TO 7
C
C   H(R) = H2 - (B*SIG*R + C1 + RF*sqrt(D))/(SIG*E).
C
      F = (BND - H2 + (B*RSIG+C1+RF*T1)/(SIG*E))*RF
C
C   Update the number of iterations NIT.
C
    7 NIT = NIT + 1
      IF (LUN .GE. 0) WRITE (LUN,130) NIT, SIG, F
  130 FORMAT (1X,3X,I2,' -- SIG = ',E15.8,', F = ',
     .        E15.8)
      IF (F0*F .LT. 0.) THEN
C
C   F0*F < 0.  Update (SNEG,FNEG) to (SG0,F0) so that F and
C     FNEG always have opposite signs.  If SIG is closer to
C     SNEG than SG0, then swap (SNEG,FNEG) with (SG0,F0).
C
        T1 = DMAX
        T2 = FNEG
        DMAX = DSIG
        FNEG = F0
        IF (ABS(DSIG) .GT. ABS(T1)) THEN
C
          DSIG = T1
          F0 = T2
        ENDIF
      ENDIF
C
C   Test for convergence.
C
      STOL = RTOL*SIG
      IF (ABS(DMAX) .LE. STOL  .OR.  (F .GE. 0.  .AND.
     .    F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 10
C
C   Test for F0 = F = FMAX or F < 0 on the first iteration.
C
      IF (F0 .NE. F  .AND.  (NIT .GT. 1  .OR.  F .GT. 0.))
     .   GO TO 9
C
C   F*F0 > 0 and either the new estimate would be outside
C     of the bracketing interval of length abs(DMAX) or
C     F < 0 on the first iteration.  Reset (SG0,F0) to
C     (SNEG,FNEG).
C
    8 DSIG = DMAX
      F0 = FNEG
C
C   Compute the change in SIG by linear interpolation
C     between (SG0,F0) and (SIG,F).
C
    9 DSIG = -F*DSIG/(F-F0)
      IF (LUN .GE. 0) WRITE (LUN,140) DSIG
  140 FORMAT (1X,8X,'DSIG = ',E15.8)
      IF ( ABS(DSIG) .GT. ABS(DMAX)  .OR.
     .     DSIG*DMAX .GT. 0. ) GO TO 8
C
C   Restrict the step-size such that abs(DSIG) .GE. STOL/2.
C     Note that DSIG and DMAX have opposite signs.
C
      IF (ABS(DSIG) .LT. STOL/2.) DSIG = -SIGN(STOL/2.,DMAX)
C
C   Bottom of loop:  Update SIG, DMAX, and F0.
C
      SIG = SIG + DSIG
      DMAX = DMAX + DSIG
      F0 = F
      GO TO 6
C
C No errors encountered.
C
   10 SIG0 = SIG
      IF (IFLGS .LE. 0) RETURN
      SIGMA(LP1) = SIG
      SIGMA(LP2) = SIG
      RETURN
C
C Error termination.
C
   11 SIG0 = -1.
      RETURN
      END
      DOUBLE PRECISION FUNCTION SIG1 (N1,N2,N,X,Y,Z,H,LIST,LPTR,LEND,
     .                    GRAD,IFLGB,HPBND,TOL,IFLGS, SIGMA, IER)
      INTEGER N1, N2, N, LIST(*), LPTR(*), LEND(N), IFLGB,
     .        IFLGS, IER
      DOUBLE PRECISION    X(N), Y(N), Z(N), H(N), GRAD(3,N), HPBND, TOL,
     .        SIGMA(*)
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/21/96
C
C   Given a triangulation of a set of nodes on the unit
C sphere, along with data values H and gradients GRAD at the
C nodes, this function determines the smallest tension fac-
C tor SIG1 such that the first derivative HP(A) of the
C Hermite interpolatory tension spline H(A), defined by SIG1
C and the endpoint values and directional derivatives asso-
C ciated with an arc N1-N2, is bounded (either above or
C below) by HPBND for all A in (A1,A2), where (A1,A2) de-
C notes an interval corresponding to the arc and A denotes
C arc-length.
C
C On input:
C
C       N1,N2 = Nodal indexes of the endpoints of an arc for
C               which the tension factor is to be computed.
C               The indexes must be distinct and lie in the
C               range 1 to N, and if IFLGS .GE. 1, they must
C               correspond to adjacent nodes in the triangu-
C               lation.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y,Z = Arrays of length N containing coordinates of
C               the nodes.  X(I)**2 + Y(I)**2 + Z(I)**2 = 1.
C
C       H = Array of length N containing data values at the
C           nodes.  H(I) is associated with (X(I),Y(I),Z(I))
C           for I = 1 to N.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to STRIPACK
C                        Subroutine TRMESH.
C
C       GRAD = Array dimensioned 3 by N whose columns con-
C              gradients at the nodes.  GRAD( ,J) must be
C              orthogonal to node J:  GRAD(1,J)*X(J) +
C              GRAD(2,J)*Y(J) + GRAD(3,J)*Z(J) = 0.  Refer
C              to Subroutines GRADG, GRADL, and SMSURF.
C
C       IFLGB = Bound option indicator:
C               IFLGB = -1 if HPBND is a lower bound on HP.
C               IFLGB = 1 if HPBND is an upper bound on HP.
C
C       HPBND = Bound on HP.  HPBND .LE. min(HP1,HP2,S) if
C               IFLGB = -1 and HPBND .GE. max(HP1,HP2,S) if
C               IFLGB = 1, where HP1 and HP2 are the direc-
C               tional derivatives at the endpoints of the
C               arc N1-N2, and S is the slope of the linear
C               interpolant of the endpoint data values.
C
C       TOL = Tolerance whose magnitude determines how close
C             SIG1 is to its optimal value when nonzero
C             finite tension is necessary and sufficient to
C             satisfy the constraint.  For a lower bound,
C             SIG1 is chosen so that HPBND .LE. HPMIN .LE.
C             HPBND + abs(TOL), where HPMIN is the minimum
C             value of HP on the arc.  For an upper bound,
C             the maximum of HP satisfies HPBND - abs(TOL)
C             .LE. HPMAX .LE. HPBND.  Thus, the constraint
C             is satisfied but possibly with more tension
C             than necessary.
C
C       IFLGS = Tension array option indicator:
C               IFLGS .LE. 0 if SIGMA is not to be used.
C               IFLGS .GE. 1 if SIGMA is to be updated by
C                            storing SIG1 in the appropriate
C                            locations.
C
C The above parameters are not altered by this function.
C
C       SIGMA = Dummy parameter (IFLGS .LE. 0) or array
C               containing tension factors associated with
C               arcs in one-to-one correspondence with LIST
C               entries (IFLGS .GE. 1).  Refer to Subroutine
C               GETSIG.
C
C On output:
C
C       SIGMA = Tension factor array updated with the new
C               value if and only if IFLGS .GE. 1 and IER
C               .GE. 0.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered and the
C                     constraint can be satisfied with fin-
C                     ite tension.
C             IER = 1 if no errors were encountered but in-
C                     finite tension is required to satisfy
C                     the constraint (e.g., IFLGB = -1,
C                     HPBND = S, and HP1 > S).
C             IER = -1 if N1, N2, N, or IFLGB is outside its
C                      valid range.
C             IER = -2 if nodes N1 and N2 coincide or IFLGS
C                      .GE. 1 and the nodes are not adja-
C                      cent.
C             IER = -3 if HPBND is outside its valid range.
C
C       SIG1 = Minimum tension factor defined above unless
C              IER < 0, in which case SIG1 = -1.  If IER
C              = 1, SIG1 is set to 85, resulting in an
C              approximation to the linear interpolant of
C              the endpoint values.
C
C STRIPACK module required by SIG1:  STORE
C
C SSRFPACK modules required by SIG1:  ARCLEN, SNHCSH
C
C Intrinsic functions called by SIG1:   ABS, EXP, MAX, MIN,
C                                         REAL, SIGN, SQRT
C
C***********************************************************
C
      DOUBLE PRECISION    ARCLEN, STORE
      INTEGER LP1, LP2, LPL, LUN, NIT
      DOUBLE PRECISION    A, A0, AL, B0, BND, C0, C1, C2, COSHM, COSHMM,
     .        D0, D1, D1PD2, D2, DMAX, DSIG, E, EMS, EMS2,
     .        F, F0, FMAX, FNEG, FTOL, P1(3), P2(3), RF,
     .        RTOL, S, S1, S2, SBIG, SIG, SINH, SINHM, STOL,
     .        T0, T1, T2, TM, UN(3), UNORM
C
      DATA SBIG/85./,  LUN/-1/
      RF = REAL(IFLGB)
      BND = HPBND
C
C Print a heading.
C
      IF (LUN .GE. 0  .AND.  RF .LT. 0.) WRITE (LUN,100) N1,
     .                                   N2, BND
      IF (LUN .GE. 0  .AND.  RF .GT. 0.) WRITE (LUN,110) N1,
     .                                   N2, BND
  100 FORMAT (//1X,'SIG1 -- N1 =',I4,', N2 =',I4,
     .        ', LOWER BOUND = ',E15.8)
  110 FORMAT (//1X,'SIG1 -- N1 =',I4,', N2 =',I4,
     .        ', UPPER BOUND = ',E15.8)
C
C Test for errors and store local parameters.
C
      IER = -1
      IF (MIN(N1,N2) .LT. 1  .OR.  N1 .EQ. N2  .OR.
     .    MAX(N1,N2,3) .GT. N  .OR.  ABS(RF) .NE. 1.)
     .   GO TO 11
      IER = -2
      IF (IFLGS .GT. 0) THEN
C
C   Set LP1 and LP2 to the pointers to N2 as a neighbor of
C     N1 and N1 as a neighbor of N2, respectively.
C
        LPL = LEND(N1)
        LP1 = LPTR(LPL)
    1   IF (LIST(LP1) .EQ. N2) GO TO 2
          LP1 = LPTR(LP1)
          IF (LP1 .NE. LPL) GO TO 1
        IF (ABS(LIST(LP1)) .NE. N2) GO TO 11
C
    2   LPL = LEND(N2)
        LP2 = LPTR(LPL)
    3   IF (LIST(LP2) .EQ. N1) GO TO 4
          LP2 = LPTR(LP2)
          IF (LP2 .NE. LPL) GO TO 3
        IF (ABS(LIST(LP2)) .NE. N1) GO TO 11
      ENDIF
C
C Store nodal coordinates P1 and P2, compute arc-length AL
C   and unit normal UN = (P1 X P2)/UNORM, and test for
C   coincident nodes.
C
    4 P1(1) = X(N1)
      P1(2) = Y(N1)
      P1(3) = Z(N1)
      P2(1) = X(N2)
      P2(2) = Y(N2)
      P2(3) = Z(N2)
      AL = ARCLEN(P1,P2)
      UN(1) = P1(2)*P2(3) - P1(3)*P2(2)
      UN(2) = P1(3)*P2(1) - P1(1)*P2(3)
      UN(3) = P1(1)*P2(2) - P1(2)*P2(1)
      UNORM = SQRT(UN(1)*UN(1) + UN(2)*UN(2) + UN(3)*UN(3))
      IF (UNORM .EQ. 0.  .OR.  AL .EQ. 0.) GO TO 11
C
C Compute first difference S and scaled directional deriva-
C   tives S1,S2 at the endpoints (for the direction N1->N2).
C
      S = H(N2) - H(N1)
      S1 = AL*(GRAD(1,N1)*P2(1) + GRAD(2,N1)*P2(2) +
     .         GRAD(3,N1)*P2(3))/UNORM
      S2 = -AL*(GRAD(1,N2)*P1(1) + GRAD(2,N2)*P1(2) +
     .          GRAD(3,N2)*P1(3))/UNORM
C
C Test for a valid constraint.
C
      IER = -3
      IF ((RF .LT. 0.  .AND.  MIN(S1,S2,S) .LT. BND)  .OR.
     .    (RF .GT. 0.  .AND.  BND .LT. MAX(S1,S2,S)))
     .   GO TO 11
C
C Test for infinite tension required.
C
      IER = 1
      SIG = SBIG
      IF (S .EQ. BND  .AND.  (S1 .NE. S  .OR.  S2 .NE. S))
     .   GO TO 10
C
C Test for SIG = 0 sufficient.  The Hermite cubic interpo-
C   lant H0 has derivative HP0(T) = (S2 + 2*B0*R + A0*R**2)/
C   AL, where R = (T2-T)/AL.
C
      IER = 0
      SIG = 0.
      T0 = 3.*S - S1 - S2
      B0 = T0 - S2
      C0 = T0 - S1
      A0 = -B0 - C0
C
C   HP0(R) has an extremum (at R = -B0/A0) in (0,1) iff
C     B0*C0 > 0 and the third derivative of H0 has the
C     sign of A0.
C
      IF (B0*C0 .LE. 0.  .OR.  A0*RF .GT. 0.) GO TO 10
C
C   A0*RF < 0 and HP0(R) = -D0/(DT*A0) at R = -B0/A0.
C
      D0 = T0*T0 - S1*S2
      F0 = (BND + D0/(A0*AL))*RF
      IF (F0 .GE. 0.) GO TO 10
C
C Find a zero of F(SIG) = (BND-HP(R))*RF, where HP has an
C   extremum at R.  F has a unique zero, F(0) = F0 < 0, and
C   F = (BND-S)*RF > 0 for SIG sufficiently large.
C
C Initialize parameters for the secant method.  The method
C   uses three points:  (SG0,F0), (SIG,F), and (SNEG,FNEG),
C   where SG0 and SNEG are defined implicitly by DSIG = SIG
C   - SG0 and DMAX = SIG - SNEG.  SG0 is initially zero and
C   SIG is initialized to the zero of (BND - (SIG*S-S1-S2)/
C   (AL*(SIG-2.)))*RF -- a value for which F(SIG) .GE. 0 and
C   F(SIG) = 0 for SIG sufficiently large that 2*SIG is in-
C   significant relative to exp(SIG).
C
      FMAX = (BND-S/AL)*RF
      SIG = 2. - A0/(3.*(AL*BND-S))
      IF (LUN .GE. 0) WRITE (LUN,120) F0, FMAX, SIG
  120 FORMAT (1X,9X,'F0 = ',E15.8,', FMAX = ',E15.8/
     .        1X,8X,'SIG = ',E15.8/)
      IF (STORE(SIG*EXP(-SIG)+.5) .EQ. .5) GO TO 10
      DSIG = SIG
      DMAX = -2.*SIG
      FNEG = FMAX
      D1 = S - S1
      D2 = S2 - S
      D1PD2 = D1 + D2
      NIT = 0
C
C Compute an absolute tolerance FTOL = abs(TOL), and a
C   relative tolerance RTOL = 100*Macheps.
C
      FTOL = ABS(TOL)
      RTOL = 1.
    5 RTOL = RTOL/2.
        IF (STORE(RTOL+1.) .GT. 1.) GO TO 5
      RTOL = RTOL*200.
C
C Top of loop:  compute F.
C
    6 IF (SIG .LE. .5) THEN
C
C   Use approximations designed to avoid cancellation
C     error (associated with small SIG) in the modified
C     hyperbolic functions.
C
        CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
        C1 = SIG*COSHM*D2 - SINHM*D1PD2
        C2 = SIG*(SINHM+SIG)*D2 - COSHM*D1PD2
        A = C2 - C1
        E = SIG*SINHM - COSHMM - COSHMM
      ELSE
C
C   Scale SINHM and COSHM by 2*exp(-SIG) in order to avoid
C     overflow.
C
        EMS = EXP(-SIG)
        EMS2 = EMS + EMS
        TM = 1. - EMS
        SINH = TM*(1.+EMS)
        SINHM = SINH - SIG*EMS2
        COSHM = TM*TM
        C1 = SIG*COSHM*D2 - SINHM*D1PD2
        C2 = SIG*SINH*D2 - COSHM*D1PD2
        A = EMS2*(SIG*TM*D2 + (TM-SIG)*D1PD2)
        E = SIG*SINH - COSHM - COSHM
      ENDIF
C
C   The second derivative HPP of H(R) has a zero at exp(SIG*
C     R) = SQRT((C2+C1)/A) and R is in (0,1) and well-
C     defined iff HPP(T1)*HPP(T2) < 0.
C
      F = FMAX
      T1 = A*(C2+C1)
      IF (T1 .GE. 0.) THEN
        IF (C1*(SIG*COSHM*D1 - SINHM*D1PD2) .LT. 0.) THEN
C
C   HP(R) = (B+SIGN(A)*SQRT(A*C))/(AL*E) at the critical
C     value of R, where A = C2-C1, B = E*S2-C2, and C = C2 +
C     C1.  Note that RF*A < 0.
C
          F = (BND - (E*S2-C2 - RF*SQRT(T1))/(AL*E))*RF
        ENDIF
      ENDIF
C
C   Update the number of iterations NIT.
C
      NIT = NIT + 1
      IF (LUN .GE. 0) WRITE (LUN,130) NIT, SIG, F
  130 FORMAT (1X,3X,I2,' -- SIG = ',E15.8,', F = ',
     .        E15.8)
      IF (F0*F .LT. 0.) THEN
C
C   F0*F < 0.  Update (SNEG,FNEG) to (SG0,F0) so that F
C     and FNEG always have opposite signs.  If SIG is closer
C     to SNEG than SG0 and abs(F) < abs(FNEG), then swap
C     (SNEG,FNEG) with (SG0,F0).
C
        T1 = DMAX
        T2 = FNEG
        DMAX = DSIG
        FNEG = F0
        IF ( ABS(DSIG) .GT. ABS(T1)  .AND.
     .       ABS(F) .LT. ABS(T2) ) THEN
C
          DSIG = T1
          F0 = T2
        ENDIF
      ENDIF
C
C   Test for convergence.
C
      STOL = RTOL*SIG
      IF (ABS(DMAX) .LE. STOL  .OR.  (F .GE. 0.  .AND.
     .    F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 10
      IF (F0*F .LT. 0.  .OR.  ABS(F) .LT. ABS(F0)) GO TO 8
C
C   F*F0 > 0 and the new estimate would be outside of the
C     bracketing interval of length abs(DMAX).  Reset
C     (SG0,F0) to (SNEG,FNEG).
C
    7 DSIG = DMAX
      F0 = FNEG
C
C   Compute the change in SIG by linear interpolation
C     between (SG0,F0) and (SIG,F).
C
    8 DSIG = -F*DSIG/(F-F0)
      IF (LUN .GE. 0) WRITE (LUN,140) DSIG
  140 FORMAT (1X,8X,'DSIG = ',E15.8)
      IF ( ABS(DSIG) .GT. ABS(DMAX)  .OR.
     .     DSIG*DMAX .GT. 0. ) GO TO 7
C
C   Restrict the step-size such that abs(DSIG) .GE. STOL/2.
C     Note that DSIG and DMAX have opposite signs.
C
      IF (ABS(DSIG) .LT. STOL/2.) DSIG = -SIGN(STOL/2.,DMAX)
C
C   Bottom of loop:  update SIG, DMAX, and F0.
C
      SIG = SIG + DSIG
      DMAX = DMAX + DSIG
      F0 = F
      GO TO 6
C
C No errors encountered.
C
   10 SIG1 = SIG
      IF (IFLGS .LE. 0) RETURN
      SIGMA(LP1) = SIG
      SIGMA(LP2) = SIG
      RETURN
C
C Error termination.
C
   11 SIG1 = -1.
      RETURN
      END
      DOUBLE PRECISION FUNCTION SIG2 (N1,N2,N,X,Y,Z,H,LIST,LPTR,LEND,
     .                    GRAD,TOL,IFLGS, SIGMA, IER)
      INTEGER N1, N2, N, LIST(*), LPTR(*), LEND(N), IFLGS,
     .        IER
      DOUBLE PRECISION    X(N), Y(N), Z(N), H(N), GRAD(3,N), TOL,
     .        SIGMA(*)
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/25/96
C
C   Given a triangulation of a set of nodes on the unit
C sphere, along with data values H and gradients GRAD at the
C nodes, this function determines the smallest tension fac-
C tor SIG2 such that the Hermite interpolatory tension
C spline H(A), defined by SIG2 and the endpoint values and
C directional derivatives associated with an arc N1-N2,
C preserves convexity (or concavity) of the data:
C
C   HP1 .LE. S .LE. HP2 implies HPP(A) .GE. 0, and
C   HP1 .GE. S .GE. HP2 implies HPP(A) .LE. 0
C
C for all A in the open interval (A1,A2) corresponding to
C the arc, where HP1 and HP2 are the derivative values of H
C at the endpoints, S is the slope of the linear interpolant
C of the endpoint data values, HPP denotes the second deriv-
C ative of H, and A is arc-length.  Note, however, that
C infinite tension is required if HP1 = S or HP2 = S (unless
C HP1 = HP2 = S).
C
C On input:
C
C       N1,N2 = Nodal indexes of the endpoints of an arc for
C               which the tension factor is to be computed.
C               The indexes must be distinct and lie in the
C               range 1 to N, and if IFLGS .GE. 1, they must
C               correspond to adjacent nodes in the triangu-
C               lation.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y,Z = Arrays of length N containing coordinates of
C               the nodes.  X(I)**2 + Y(I)**2 + Z(I)**2 = 1.
C
C       H = Array of length N containing data values at the
C           nodes.  H(I) is associated with (X(I),Y(I),Z(I))
C           for I = 1 to N.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to STRIPACK
C                        Subroutine TRMESH.
C
C       GRAD = Array dimensioned 3 by N whose columns con-
C              gradients at the nodes.  GRAD( ,J) must be
C              orthogonal to node J:  GRAD(1,J)*X(J) +
C              GRAD(2,J)*Y(J) + GRAD(3,J)*Z(J) = 0.  Refer
C              to Subroutines GRADG, GRADL, and SMSURF.
C
C       TOL = Tolerance whose magnitude determines how close
C             SIG2 is to its optimal value when nonzero
C             finite tension is necessary and sufficient to
C             satisfy convexity or concavity.  In the case
C             convexity, SIG2 is chosen so that 0 .LE.
C             HPPMIN .LE. abs(TOL), where HPPMIN is the
C             minimum value of HPP on the arc.  In the case
C             of concavity, the maximum value of HPP satis-
C             fies -abs(TOL) .LE. HPPMAX .LE. 0.  Thus, the
C             constraint is satisfied but possibly with more
C             tension than necessary.
C
C       IFLGS = Tension array option indicator:
C               IFLGS .LE. 0 if SIGMA is not to be used.
C               IFLGS .GE. 1 if SIGMA is to be updated by
C                            storing SIG2 in the appropriate
C                            locations.
C
C The above parameters are not altered by this function.
C
C       SIGMA = Dummy parameter (IFLGS .LE. 0) or array
C               containing tension factors associated with
C               arcs in one-to-one correspondence with LIST
C               entries (IFLGS .GE. 1).  Refer to Subroutine
C               GETSIG.
C
C On output:
C
C       SIGMA = Tension factor array updated with the new
C               value if and only if IFLGS .GE. 1 and IER
C               .GE. 0.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered and fin-
C                     ite tension is sufficient to satisfy
C                     convexity (or concavity).
C             IER = 1 if no errors were encountered but in-
C                     finite tension is required to satisfy
C                     convexity.
C             IER = 2 if the data does not satisfy convexity
C                     or concavity.
C             IER = -1 if N1, N2, or N is outside its valid
C                      range.
C             IER = -2 if nodes N1 and N2 coincide or IFLGS
C                      .GE. 1 and the nodes are not adja-
C                      cent.
C
C       SIG2 = Minimum tension factor defined above unless
C              IER < 0, in which case SIG2 = -1.  If IER
C              = 1, SIG2 is set to 85, resulting in an
C              approximation to the linear interpolant of
C              the endpoint values.  If IER = 2, SIG2 = 0,
C              resulting in the Hermite cubic interpolant.
C
C STRIPACK module required by SIG2:  STORE
C
C SSRFPACK modules required by SIG2:  ARCLEN, SNHCSH
C
C Intrinsic functions called by SIG2:  ABS, EXP, MAX, MIN,
C                                        SQRT
C
C***********************************************************
C
      DOUBLE PRECISION    ARCLEN, STORE
      INTEGER LP1, LP2, LPL, LUN, NIT
      DOUBLE PRECISION    AL, COSHM, D1, D1D2, D2, DSIG, DUMMY, EMS, F,
     .        FP, FTOL, P1(3), P2(3), RTOL, S, SBIG, SIG,
     .        SINHM, SSM, T, T1, TP1, UN(3), UNORM
C
      DATA SBIG/85./,  LUN/-1/
C
C Print a heading.
C
      IF (LUN .GE. 0) WRITE (LUN,100) N1, N2
  100 FORMAT (//1X,'SIG2 -- N1 =',I4,', N2 =',I4)
C
C Test for errors and set local parameters.
C
      IER = -1
      IF (MIN(N1,N2) .LT. 1  .OR.  N1 .EQ. N2  .OR.
     .    MAX(N1,N2,3) .GT. N) GO TO 11
      IER = -2
      IF (IFLGS .GT. 0) THEN
C
C   Set LP1 and LP2 to the pointers to N2 as a neighbor of
C     N1 and N1 as a neighbor of N2, respectively.
C
        LPL = LEND(N1)
        LP1 = LPTR(LPL)
    1   IF (LIST(LP1) .EQ. N2) GO TO 2
          LP1 = LPTR(LP1)
          IF (LP1 .NE. LPL) GO TO 1
        IF (ABS(LIST(LP1)) .NE. N2) GO TO 11
C
    2   LPL = LEND(N2)
        LP2 = LPTR(LPL)
    3   IF (LIST(LP2) .EQ. N1) GO TO 4
          LP2 = LPTR(LP2)
          IF (LP2 .NE. LPL) GO TO 3
        IF (ABS(LIST(LP2)) .NE. N1) GO TO 11
      ENDIF
C
C Store nodal coordinates P1 and P2, compute arc-length AL
C   and unit normal UN = (P1 X P2)/UNORM, and test for
C   coincident nodes.
C
    4 P1(1) = X(N1)
      P1(2) = Y(N1)
      P1(3) = Z(N1)
      P2(1) = X(N2)
      P2(2) = Y(N2)
      P2(3) = Z(N2)
      AL = ARCLEN(P1,P2)
      UN(1) = P1(2)*P2(3) - P1(3)*P2(2)
      UN(2) = P1(3)*P2(1) - P1(1)*P2(3)
      UN(3) = P1(1)*P2(2) - P1(2)*P2(1)
      UNORM = SQRT(UN(1)*UN(1) + UN(2)*UN(2) + UN(3)*UN(3))
      IF (UNORM .EQ. 0.  .OR.  AL .EQ. 0.) GO TO 11
C
C Compute first and second differences and test for infinite
C   tension required.
C
      S = H(N2) - H(N1)
      D1 = S - AL*(GRAD(1,N1)*P2(1) + GRAD(2,N1)*P2(2) +
     .             GRAD(3,N1)*P2(3))/UNORM
      D2 = -AL*(GRAD(1,N2)*P1(1) + GRAD(2,N2)*P1(2) +
     .          GRAD(3,N2)*P1(3))/UNORM - S
      D1D2 = D1*D2
      IER = 1
      SIG = SBIG
      IF (D1D2 .EQ. 0.  .AND.  D1 .NE. D2) GO TO 10
C
C Test for a valid constraint.
C
      IER = 2
      SIG = 0.
      IF (D1D2 .LT. 0.) GO TO 10
C
C Test for SIG = 0 sufficient.
C
      IER = 0
      IF (D1D2 .EQ. 0.) GO TO 10
      T = MAX(D1/D2,D2/D1)
      IF (T .LE. 2.) GO TO 10
C
C Find a zero of F(SIG) = SIG*COSHM(SIG)/SINHM(SIG) - (T+1).
C   Since the derivative of F vanishes at the origin, a
C   quadratic approximation is used to obtain an initial
C   estimate for the Newton method.
C
      TP1 = T + 1.
      SIG = SQRT(10.*T-20.)
      NIT = 0
C
C   Compute an absolute tolerance FTOL = abs(TOL) and a
C     relative tolerance RTOL = 100*Macheps.
C
      FTOL = ABS(TOL)
      RTOL = 1.
    5 RTOL = RTOL/2.
        IF (STORE(RTOL+1.) .GT. 1.) GO TO 5
      RTOL = RTOL*200.
C
C Top of loop:  evaluate F and its derivative FP.
C
    6 IF (SIG .LE. .5) THEN
C
C   Use approximations designed to avoid cancellation error
C     in the hyperbolic functions.
C
        CALL SNHCSH (SIG, SINHM,COSHM,DUMMY)
        T1 = COSHM/SINHM
        FP = T1 + SIG*(SIG/SINHM - T1*T1 + 1.)
      ELSE
C
C   Scale SINHM and COSHM by 2*exp(-SIG) in order to avoid
C     overflow.
C
        EMS = EXP(-SIG)
        SSM = 1. - EMS*(EMS+SIG+SIG)
        T1 = (1.-EMS)*(1.-EMS)/SSM
        FP = T1 + SIG*(2.*SIG*EMS/SSM - T1*T1 + 1.)
      ENDIF
C
      F = SIG*T1 - TP1
C
C   Update the number of iterations NIT.
C
      NIT = NIT + 1
      IF (LUN .GE. 0) WRITE (LUN,110) NIT, SIG, F, FP
  110 FORMAT (1X,3X,I2,' -- SIG = ',E15.8,', F = ',
     .        E15.8/1X,31X,'FP = ',E15.8)
C
C   Test for convergence.
C
      IF (FP .LE. 0.) GO TO 10
      DSIG = -F/FP
      IF (ABS(DSIG) .LE. RTOL*SIG  .OR.  (F .GE. 0.  .AND.
     .    F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 10
C
C   Bottom of loop:  update SIG.
C
      SIG = SIG + DSIG
      GO TO 6
C
C No errors encountered.
C
   10 SIG2 = SIG
      IF (IFLGS .LE. 0) RETURN
      SIGMA(LP1) = SIG
      SIGMA(LP2) = SIG
      RETURN
C
C Error termination.
C
   11 SIG2 = -1.
      RETURN
      END
      SUBROUTINE SMSGS (N,X,Y,Z,U,LIST,LPTR,LEND,IFLGS,
     .                  SIGMA,W,P, NIT,DFMAX,F,GRAD, IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), IFLGS, NIT, IER
      DOUBLE PRECISION    X(N), Y(N), Z(N), U(N), SIGMA(*), W(N), P,
     .        DFMAX, F(N), GRAD(3,N)
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/25/96
C
C   This subroutine solves the symmetric positive definite
C linear system associated with minimizing the quadratic
C functional Q(F,FX,FY,FZ) described in Subroutine SMSURF.
C Since the gradient at node K lies in the plane tangent to
C the sphere surface at K, it is effectively defined by only
C two components -- its X and Y components in the coordinate
C system obtained by rotating K to the north pole.  Thus,
C the minimization problem corresponds to an order-3N system
C which is solved by the block Gauss-Seidel method with 3 by
C 3 blocks.
C
C On input:
C
C       N,X,Y,Z,U,LIST,LPTR,LEND,IFLGS,SIGMA,W = Parameters
C           as described in Subroutine SMSURF.
C
C       P = Positive smoothing parameter defining Q.
C
C The above parameters are not altered by this routine.
C
C       NIT = Maximum number of iterations to be used.  This
C             maximum will likely be achieved if DFMAX is
C             smaller than the machine precision.  NIT .GE.
C             0.
C
C       DFMAX = Nonnegative convergence criterion.  The
C               method is terminated when the maximum
C               change in a solution F-component between
C               iterations is at most DFMAX.  The change in
C               a component is taken to be the absolute
C               difference relative to 1 plus the old value.
C
C       F = Initial estimate of the first N solution compo-
C           nents.
C
C       GRAD = 3 by N array containing initial estimates of
C              the last 3N solution components (the gradi-
C              ent with FX, FY, and FZ in rows 1, 2, and 3,
C              respectively).
C
C On output:
C
C       NIT = Number of Gauss-Seidel iterations employed.
C
C       DFMAX = Maximum relative change in a solution F-
C               component at the last iteration.
C
C       F = First N solution components -- function values
C           at the nodes.
C
C       GRAD = Last 3N solution components -- gradients at
C              the nodes.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered and the
C                     convergence criterion was achieved.
C             IER = 1 if no errors were encountered but con-
C                     vergence was not achieved within NIT
C                     iterations.
C             IER = -1 if N, P, NIT, or DFMAX is outside its
C                      valid range on input.  F and GRAD are
C                      not altered in this case.
C             IER = -2 if all nodes are collinear or the
C                      triangulation is invalid.
C             IER = -3 if duplicate nodes were encountered.
C
C SSRFPACK modules required by SMSGS:  APLYRT, CONSTR,
C                                        GRCOEF, SNHCSH
C
C Intrinsic functions called by SMSGS:  ABS, ATAN, MAX, SQRT
C
C***********************************************************
C
      INTEGER IFL, ITER, ITMAX, J, K, LPJ, LPL, NN
      DOUBLE PRECISION    ALFA, ALFSQ, C11, C12, C13, C22, C23, C33,
     .        CC22, CC23, CC33, CX, CY, DEN1, DEN2, DET, DF,
     .        DFMX, DGK(3), DGX, DGY, FK, G1, G2, G3, GJK,
     .        GKJ, PP, R1, R2, R3, RR2, RR3, SIG, SINAL, SX,
     .        SY, T, T1, T2, T3, T4, T5, T6, TOL, XJ, XK,
     .        XS, YJ, YK, YS, ZJ, ZK
C
      NN = N
      IFL = IFLGS
      PP = P
      ITMAX = NIT
      TOL = DFMAX
C
C Test for errors in input and initialize iteration count,
C   tension factor, and output value of DFMAX.
C
      IF (NN .LT. 3  .OR.  PP .LE. 0.  .OR.  ITMAX .LT. 0
     .    .OR.  TOL .LT. 0.) GO TO 5
      ITER = 0
      SIG = SIGMA(1)
      DFMX = 0.
C
C Top of iteration loop.
C
    1 IF (ITER .EQ. ITMAX) GO TO 4
      DFMX = 0.
C
C   Loop on nodes.
C
      DO 3 K = 1,NN
        XK = X(K)
        YK = Y(K)
        ZK = Z(K)
        FK = F(K)
        G1 = GRAD(1,K)
        G2 = GRAD(2,K)
        G3 = GRAD(3,K)
C
C   Construct the rotation mapping node K to the north pole.
C
        CALL CONSTR (XK,YK,ZK, CX,SX,CY,SY)
C
C   Initialize components of the order-3 system for the
C     change (DF,DGX,DGY) in the K-th solution components.
C
        C11 = PP*W(K)
        C12 = 0.
        C13 = 0.
        C22 = 0.
        C23 = 0.
        C33 = 0.
        R1 = C11*(U(K)-FK)
        R2 = 0.
        R3 = 0.
C
C   Loop on neighbors J of node K.
C
        LPL = LEND(K)
        LPJ = LPL
    2   LPJ = LPTR(LPJ)
          J = ABS(LIST(LPJ))
C
C   Compute the coordinates of J in the rotated system.
C
          T = SX*Y(J) + CX*Z(J)
          YJ = CX*Y(J) - SX*Z(J)
          ZJ = SY*X(J) + CY*T
          XJ = CY*X(J) - SY*T
C
C   Compute arc-length ALFA between K and J, ALFSQ = ALFA*
C     ALFA, SINAL = SIN(ALFA), DEN1 = ALFA*SIN(ALFA)**2, and
C     DEN2 = ALFSQ*SINAL.
C
          ALFA = 2.*ATAN(SQRT((1.-ZJ)/(1.+ZJ)))
          ALFSQ = ALFA*ALFA
          XS = XJ*XJ
          YS = YJ*YJ
          SINAL = SQRT(XS+YS)
          DEN1 = ALFA*(XS+YS)
          DEN2 = ALFSQ*SINAL
C
C   Test for coincident nodes and compute functions of SIG:
C     T1 = SIG*SIG*COSHM/E, T2 = SIG*SINHM/E, and T3 = SIG*
C     (SIG*COSHM-SINHM)/E for E = SIG*SINH - 2*COSHM.
C
          IF (DEN1 .EQ. 0.) GO TO 7
          IF (IFL .GE. 1) SIG = SIGMA(LPJ)
          CALL GRCOEF (SIG, T3,T2)
          T1 = T2 + T3
C
C   Update system components for node J.
C
          T4 = 2.*T1/(ALFA*ALFSQ)
          T5 = T1/DEN2
          T6 = T3/DEN1
          C11 = C11 + T4
          C12 = C12 + T5*XJ
          C13 = C13 + T5*YJ
          C22 = C22 + T6*XS
          C23 = C23 + T6*XJ*YJ
          C33 = C33 + T6*YS
          GKJ = G1*X(J) + G2*Y(J) + G3*Z(J)
          GJK = GRAD(1,J)*XK + GRAD(2,J)*YK + GRAD(3,J)*ZK
          R1 = R1 + T4*(F(J)-FK) + T5*(GJK-GKJ)
          T = T5*(F(J)-FK) - T6*GKJ + T2*GJK/DEN1
          R2 = R2 + T*XJ
          R3 = R3 + T*YJ
C
C   Bottom of loop on neighbors.
C
          IF (LPJ .NE. LPL) GO TO 2
C
C   Solve the system associated with the K-th block.
C
        CC22 = C11*C22 - C12*C12
        CC23 = C11*C23 - C12*C13
        CC33 = C11*C33 - C13*C13
        RR2 = C11*R2 - C12*R1
        RR3 = C11*R3 - C13*R1
        DET = CC22*CC33 - CC23*CC23
        IF (DET .EQ. 0.  .OR.  CC22 .EQ. 0.  .OR.
     .      C11 .EQ. 0.) GO TO 6
        DGY = (CC22*RR3 - CC23*RR2)/DET
        DGX = (RR2 - CC23*DGY)/CC22
        DF = (R1 - C12*DGX - C13*DGY)/C11
C
C   Rotate (DGX,DGY,0) back to the original coordinate
C     system, and update GRAD( ,K), F(K), and DFMX.
C
        CALL APLYRT (DGX,DGY,CX,SX,CY,SY, DGK)
        GRAD(1,K) = G1 + DGK(1)
        GRAD(2,K) = G2 + DGK(2)
        GRAD(3,K) = G3 + DGK(3)
        F(K) = FK + DF
        DFMX = MAX(DFMX,ABS(DF)/(1.+ABS(FK)))
    3   CONTINUE
C
C   Increment ITER and test for convergence.
C
      ITER = ITER + 1
      IF (DFMX .GT. TOL) GO TO 1
C
C The method converged.
C
      NIT = ITER
      DFMAX = DFMX
      IER = 0
      RETURN
C
C The method failed to converge within NIT iterations.
C
    4 DFMAX = DFMX
      IER = 1
      RETURN
C
C Invalid input parameter.
C
    5 NIT = 0
      DFMAX = 0.
      IER = -1
      RETURN
C
C Node K and its neighbors are collinear.
C
    6 NIT = 0
      DFMAX = DFMX
      IER = -2
      RETURN
C
C Nodes J and K coincide.
C
    7 NIT = 0
      DFMAX = DFMX
      IER = -3
      RETURN
      END
      SUBROUTINE SMSURF (N,X,Y,Z,U,LIST,LPTR,LEND,IFLGS,
     .                   SIGMA,W,SM,SMTOL,GSTOL,LPRNT, F,
     .                   GRAD,IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), IFLGS, LPRNT,
     .        IER
      DOUBLE PRECISION    X(N), Y(N), Z(N), U(N), SIGMA(*), W(N), SM,
     .        SMTOL, GSTOL, F(N), GRAD(3,N)
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/21/98
C
C   Given a triangulation of N nodes on the unit sphere with
C data values U at the nodes and tension factors SIGMA
C associated with the arcs, this routine determines a set of
C nodal function values F and gradients GRAD = (FX,FY,FZ)
C such that a quadratic functional Q1(F,GRAD) is minimized
C subject to the constraint Q2(F) .LE. SM for Q2(F) =
C (U-F)**T*W*(U-F), where W is a diagonal matrix of positive
C weights.  The functional Q1 is an approximation to the
C linearized curvature over the triangulation of a C-1 fun-
C ction F(V), V a unit vector, which interpolates the nodal
C values and gradients.  Subroutines INTRC1 and UNIF may be
C called to evaluate F at arbitrary points.
C
C   The smoothing procedure is an extension of the method
C for cubic spline smoothing due to C. Reinsch -- Numer.
C Math., 10 (1967) and 16 (1971).  Refer to Function FVAL
C for a further description of the interpolant F.  Letting
C D1F(T) and D2F(T) denote first and second derivatives of F
C with respect to a parameter T varying along a triangula-
C tion arc, Q1 is the sum of integrals over the arcs of
C D2F(T)**2 + ((SIGMA/L)*(D1F(T)-S))**2 where L denotes arc-
C length, SIGMA is the appropriate tension factor, and S is
C the slope of the linear function of T which interpolates
C the values of F at the endpoints of the arc.  Introducing
C a smoothing parameter P, and assuming the constraint is
C active, the problem is equivalent to minimizing Q(P,F,
C GRAD) = Q1(F,GRAD) + P*(Q2(F)-SM).  The secant method is
C used to find a zero of G(P) = 1/SQRT(Q2) - 1/SQRT(SM)
C where F(P) satisfies the order-3N symmetric positive def-
C inite linear system obtained by setting the gradient of Q
C (treated as a function of F and GRAD with GRAD tangent to
C the sphere surface) to zero.  The linear system is solved
C by the block Gauss-Seidel method (refer to SMSGS).
C
C   Note that the method can also be used to select grad-
C ients for the interpolation problem (F = U, SM = 0, and P
C infinite).  This is achieved by a call to Subroutine
C GRADG.
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y,Z = Arrays of length N containing Cartesian
C               coordinates of the nodes.
C
C       U = Array of length N containing data values at the
C           nodes.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to STRIPACK
C                        Subroutine TRMESH.
C
C       IFLGS = Tension factor option:
C               IFLGS .LE. 0 if a single uniform tension
C                            factor is to be used.
C               IFLGS .GE. 1 if variable tension is desired.
C
C       SIGMA = Uniform tension factor (IFLGS .LE. 0), or
C               array containing tension factors associated
C               with arcs in one-to-one correspondence with
C               LIST entries (IFLGS .GE. 1).  Refer to Sub-
C               programs GETSIG, SIG0, SIG1, and SIG2.
C
C       W = Array of length N containing positive weights
C           associated with the data values.  The recommend-
C           ed value of W(I) is 1/DU**2 where DU is the
C           standard deviation associated with U(I).  DU**2
C           is the expected value of the squared error in
C           the measurement of U(I).  (The mean error is
C           assumed to be zero.)
C
C       SM = Positive parameter specifying an upper bound on
C            Q2(F).  Note that F is constant (and Q2(F)
C            is minimized) if SM is sufficiently large that
C            the constraint is not active.  It is recommend-
C            ed that SM satisfy N-SQRT(2N) .LE. SM .LE. N+
C            SQRT(2N).
C
C       SMTOL = Parameter in the open interval (0,1) speci-
C               fying the relative error allowed in satisfy-
C               ing the constraint -- the constraint is
C               assumed to be satisfied if SM*(1-SMTOL) .LE.
C               Q2 .LE. SM*(1+SMTOL).  A reasonable value
C               for SMTOL is SQRT(2/N).
C
C       GSTOL = Nonnegative tolerance defining the conver-
C               gence criterion for the Gauss-Seidel method.
C               Refer to parameter DFMAX in Subroutine
C               SMSGS.  A recommended value is .05*DU**2,
C               where DU is an average standard deviation
C               in the data values.
C
C       LPRNT = Logical unit on which diagnostic messages
C               are printed, or negative integer specifying
C               no diagnostics.  For each secant iteration,
C               the following values are printed:  P, G(P),
C               NIT, DFMAX, and DP, where NIT denotes the
C               number of Gauss-Seidel iterations used in
C               the computation of G, DFMAX denotes the max-
C               imum relative change in a solution component
C               in the last Gauss-Seidel iteration, and DP
C               is the change in P computed by linear inter-
C               polation between the current point (P,G) and
C               a previous point.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       F = Array of length N containing nodal function val-
C           ues unless IER < 0.
C
C       GRAD = 3 by N array whose columns contain gradients
C              of F at the nodes unless IER < 0.
C
C       IER = Error indicator and information flag:
C             IER = 0 if no errors were encountered and the
C                     constraint is active -- Q2(F) is ap-
C                     proximately equal to SM.
C             IER = 1 if no errors were encountered but the
C                     constraint is not active -- F and GRAD
C                     are the values and gradients of a con-
C                     stant function which minimizes Q2(F),
C                     and Q1 = 0.
C             IER = 2 if the constraint could not be satis-
C                     fied to within SMTOL due to
C                     ill-conditioned linear systems.
C             IER = -1 if N, W, SM, SMTOL, or GSTOL is out-
C                      side its valid range on input.
C             IER = -2 if all nodes are collinear or the
C                      triangulation is invalid.
C             IER = -3 if duplicate nodes were encountered.
C
C SSRFPACK modules required by SMSURF:  APLYRT, CONSTR,
C                                         GRCOEF, SMSGS,
C                                         SNHCSH
C
C Intrinsic functions called by SMSURF:  ABS, SQRT
C
C***********************************************************
C
      INTEGER I, IERR, ITER, ITMAX, LUN, NIT, NITMAX, NN
      DOUBLE PRECISION    C, DFMAX, DMAX, DP, G, G0, GNEG, P, Q2, Q2MAX,
     .        Q2MIN, S, SUMW, TOL, WI
C
C Local parameters:
C
C ITMAX = Maximum number of secant iterations.
C LUN = Local copy of LPRNT.
C NITMAX = Maximum number of Gauss-Seidel iterations for
C          each secant iteration.
C NN = Local copy of N.
C TOL = Local copy of GSTOL.
C
      DATA ITMAX/50/,  NITMAX/40/
C
      NN = N
      TOL = GSTOL
      LUN = LPRNT
      IF (LUN .GT. 99) LUN = -1
C
C Test for errors and initialize F to the weighted least
C   squares fit of a constant function to the data.
C
      IER = -1
      IF (NN .LT. 3  .OR.  SM .LE. 0.  .OR.  SMTOL .LE. 0.
     .    .OR.  SMTOL .GE. 1.  .OR.  TOL .LE. 0.) RETURN
      C = 0.
      SUMW = 0.
      DO 1 I = 1,NN
        WI = W(I)
        IF (WI .LE. 0.) RETURN
        C = C + WI*U(I)
        SUMW = SUMW + WI
    1   CONTINUE
      C = C/SUMW
C
C Compute nodal values and gradients, and accumulate Q2 =
C   (U-F)**T*W*(U-F).
C
      Q2 = 0.
      DO 2 I = 1,NN
        F(I) = C
        GRAD(1,I) = 0.
        GRAD(2,I) = 0.
        GRAD(3,I) = 0.
        Q2 = Q2 + W(I)*(U(I)-F(I))**2
    2   CONTINUE
C
C Compute bounds on Q2 defined by SMTOL, and test for the
C   constraint satisfied by the constant fit.
C
      Q2MIN = SM*(1.-SMTOL)
      Q2MAX = SM*(1.+SMTOL)
      IF (Q2 .LE. Q2MAX) THEN
C
C The constraint is satisfied by a constant function.
C
        IER = 1
        IF (LUN .GE. 0) WRITE (LUN,100)
  100   FORMAT (///1X,'SMSURF -- THE CONSTRAINT IS NOT ',
     .          'ACTIVE AND THE FITTING FCN IS CONSTANT.')
        RETURN
      ENDIF
C
C Compute G0 = G(0) and print a heading.
C
      IER = 0
      S = 1./SQRT(SM)
      G0 = 1./SQRT(Q2) - S
      IF (LUN .GE. 0) WRITE (LUN,110) SM, TOL, NITMAX, G0
  110 FORMAT (///1X,'SMSURF -- SM = ',E10.4,', GSTOL = ',
     .        E7.1,', NITMAX = ',I2,', G(0) = ',E15.8)
C
C G(P) is strictly increasing and concave, and G(0) .LT. 0.
C   Initialize parameters for the secant method.  The method
C   uses three points -- (P0,G0), (P,G), and (PNEG,GNEG)
C   where P0 and PNEG are defined implicitly by DP = P - P0
C   and DMAX = P - PNEG.
C
      P = 10.*SM
      DP = P
      DMAX = 0.
      ITER = 0
C
C Top of loop -- compute G.
C
    3 NIT = NITMAX
      DFMAX = TOL
      CALL SMSGS (NN,X,Y,Z,U,LIST,LPTR,LEND,IFLGS,SIGMA,W,
     .            P, NIT,DFMAX,F,GRAD, IERR)
      IF (IERR .LT. 0) IER = IERR
C
C   IERR = -1 in SMSGS could be caused by P = 0 as a result
C     of inaccurate solutions to ill-conditioned systems.
C
      IF (IERR .EQ. -1) IER = 2
      IF (IERR .LT. 0) RETURN
      Q2 = 0.
      DO 4 I = 1,NN
        Q2 = Q2 + W(I)*(U(I)-F(I))**2
    4   CONTINUE
      G = 1./SQRT(Q2) - S
      ITER = ITER + 1
      IF (LUN .GE. 0) WRITE (LUN,120) ITER, P, G, NIT, DFMAX
  120 FORMAT (/1X,I2,' -- P = ',E15.8,', G = ',E15.8,
     .        ', NIT = ',I2,', DFMAX = ',E12.6)
C
C   Test for convergence.
C
      IF (Q2MIN .LE. Q2  .AND.  Q2 .LE. Q2MAX) RETURN
      IF (ITER .GE. ITMAX) THEN
        IER = 2
        RETURN
      ENDIF
      IF (DMAX .EQ. 0.  .AND.  G .LE. 0.) THEN
C
C   Increase P until G(P) > 0.
C
        P = 10.*P
        DP = P
        GO TO 3
      ENDIF
C
C   A bracketing interval [P0,P] has been found.
C
      IF (G0*G .LE. 0.) THEN
C
C   G0*G < 0.  Update (PNEG,GNEG) to (P0,G0) so that G
C     and GNEG always have opposite signs.
C
        DMAX = DP
        GNEG = G0
      ENDIF
C
C   Compute the change in P by linear interpolation between
C     (P0,G0) and (P,G).
C
    5 DP = -G*DP/(G-G0)
      IF (LUN .GE. 0) WRITE (LUN,130) DP
  130 FORMAT (1X,5X,'DP = ',E15.8)
      IF (ABS(DP) .GT. ABS(DMAX)) THEN
C
C   G0*G .GT. 0 and the new estimate would be outside of the
C     bracketing interval of length ABS(DMAX).  Reset
C     (P0,G0) to (PNEG,GNEG).
C
        DP = DMAX
        G0 = GNEG
        GO TO 5
      ENDIF
C
C   Bottom of loop -- update P, DMAX, and G0.
C
      P = P + DP
      DMAX = DMAX + DP
      G0 = G
      GO TO 3
      END
      SUBROUTINE SNHCSH (X, SINHM,COSHM,COSHMM)
      DOUBLE PRECISION X, SINHM, COSHM, COSHMM
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   03/18/90
C
C   This subroutine computes approximations to the modified
C hyperbolic functions defined below with relative error
C bounded by 4.7E-12 for a floating point number system with
C sufficient precision.  For IEEE standard single precision,
C the relative error is less than 1.E-5 for all x.
C
C   Note that the 13-digit constants in the data statements
C below may not be acceptable to all compilers.
C
C On input:
C
C       X = Point at which the functions are to be
C           evaluated.
C
C X is not altered by this routine.
C
C On output:
C
C       SINHM = sinh(X) - X.
C
C       COSHM = cosh(X) - 1.
C
C       COSHMM = cosh(X) - 1 - X*X/2.
C
C Modules required by SNHCSH:  None
C
C Intrinsic functions called by SNHCSH:  ABS, EXP
C
C***********************************************************
C
      DOUBLE PRECISION AX, C1, C2, C3, C4, EXPX, F, XC, XS, XSD2, XSD4
C
      DATA C1/.1666666666659E0/,
     .     C2/.8333333431546E-2/,
     .     C3/.1984107350948E-3/,
     .     C4/.2768286868175E-5/
      AX = ABS(X)
      XS = AX*AX
      IF (AX .LE. .5) THEN
C
C Approximations for small X:
C
        XC = X*XS
        SINHM = XC*(((C4*XS+C3)*XS+C2)*XS+C1)
        XSD4 = .25*XS
        XSD2 = XSD4 + XSD4
        F = (((C4*XSD4+C3)*XSD4+C2)*XSD4+C1)*XSD4
        COSHMM = XSD2*F*(F+2.)
        COSHM = COSHMM + XSD2
      ELSE
C
C Approximations for large X:
C
        EXPX = EXP(AX)
        SINHM = -(((1./EXPX+AX)+AX)-EXPX)/2.
        IF (X .LT. 0.) SINHM = -SINHM
        COSHM = ((1./EXPX-2.)+EXPX)/2.
        COSHMM = COSHM - XS/2.
      ENDIF
      RETURN
      END
      SUBROUTINE UNIF (N,X,Y,Z,F,LIST,LPTR,LEND,IFLGS,SIGMA,
     .                 NROW,NI,NJ,PLAT,PLON,IFLGG, GRAD, FF,
     .                 IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), IFLGS, NROW, NI,
     .        NJ, IFLGG, IER
      DOUBLE PRECISION    X(N), Y(N), Z(N), F(N), SIGMA(*), PLAT(NI),
     .        PLON(NJ), GRAD(3,N), FF(NROW,NJ)
C
C***********************************************************
C
C                                              From SSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/25/96
C
C   Given a Delaunay triangulation of a set of nodes on the
C unit sphere, along with data values and tension factors
C associated with the triangulation arcs, this routine
C interpolates the data values to a uniform grid for such
C applications as contouring.  The interpolant is once con-
C tinuously differentiable.  Extrapolation is performed at
C grid points exterior to the triangulation when the nodes
C do not cover the entire sphere.
C
C On input:
C
C       N = Number of nodes.  N .GE. 3 and N .GE. 7 if
C           IFLAG .NE. 1.
C
C       X,Y,Z = Arrays containing Cartesian coordinates of
C               the nodes.  X(I)**2 + Y(I)**2 + Z(I)**2 = 1
C               for I = 1 to N.
C
C       F = Array containing data values.  F(I) is associ-
C           ated with (X(I),Y(I),Z(I)).
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to STRIPACK
C                        Subroutine TRMESH.
C
C       IFLGS = Tension factor option:
C               IFLGS .LE. 0 if a single uniform tension
C                            factor is to be used.
C               IFLGS .GE. 1 if variable tension is desired.
C
C       SIGMA = Uniform tension factor (IFLGS .LE. 0), or
C               array containing tension factors associated
C               with arcs in one-to-one correspondence with
C               LIST entries (IFLGS .GE. 1).  Refer to Sub-
C               programs GETSIG, SIG0, SIG1, and SIG2.
C
C       NROW = Number of rows in the dimension statement of
C              FF.
C
C       NI,NJ = Number of rows and columns in the uniform
C               grid.  1 .LE. NI .LE. NROW and 1 .LE. NJ.
C
C       PLAT,PLON = Arrays of length NI and NJ, respective-
C                   ly, containing the latitudes and
C                   longitudes of the grid lines.
C
C       IFLGG = Option indicator:
C               IFLGG = 0 if gradient estimates at the ver-
C                         tices of a triangle are to be
C                         recomputed for each grid point in
C                         the triangle and not saved.
C               IFLGG = 1 if gradient estimates are input in
C                         GRAD.
C               IFLGG = 2 if gradient estimates are to be
C                         computed once for each node (by
C                         GRADL) and saved in GRAD.
C
C The above parameters are not altered by this routine.
C
C       GRAD = 3 by N array whose columns contain the X, Y,
C              and Z components (in that order) of the grad-
C              ients at the nodes if IFLGG = 1, array of
C              sufficient size if IFLGG = 2, or dummy para-
C              meter if IFLGG = 0.
C
C Gradient estimates may be computed by Subroutines GRADL or
C   GRADG if IFLGG = 1.
C
C       FF = NROW by NCOL array with NROW .GE. NI and NCOL
C            .GE. NJ.
C
C On output:
C
C       GRAD = Array containing estimated gradients as de-
C              fined above if IFLGG = 2 and IER .GE. 0.
C              GRAD is not altered if IFLGG .NE. 2.
C
C       FF = Interpolated values at the grid points if IER
C            .GE. 0.  FF(I,J) = F(PLAT(I),PLON(J)) for I =
C            1,...,NI and J = 1,...,NJ.
C
C       IER = Error indicator:
C             IER = K if no errors were encountered and K
C                     grid points required extrapolation for
C                     K .GE. 0.
C             IER = -1 if N, NI, NJ, or IFLGG is outside its
C                      valid range.
C             IER = -2 if the nodes are collinear.
C             IER = -3 if extrapolation failed due to the
C                      uniform grid extending too far beyond
C                      the triangulation boundary.
C
C STRIPACK modules required by UNIF:  GETNP, JRAND, LSTPTR,
C                                       STORE, TRFIND
C
C SSRFPACK modules required by UNIF:  APLYR, APLYRT, ARCINT,
C                                       ARCLEN, CONSTR,
C                                       FVAL, GIVENS, GRADL,
C                                       HVAL, INTRC1,
C                                       ROTATE, SETUP,
C                                       SNHCSH
C
C***********************************************************
C
      INTEGER I, J, IERR, IFL, IST, NEX, NN, NST, NX, NY
      DATA    NST/1/
C
C Local parameters:
C
C I,J =   DO-loop indexes
C IERR =  Error flag for calls to GRADL and INTRC1
C IFL =   Local copy of IFLGG
C IST =   Parameter for INTRC1
C NEX =   Number of grid points exterior to the triangula-
C           tion boundary (number of extrapolated values)
C NN =    Local copy of N
C NST =   Initial value for IST
C NX,NY = Local copies of NI and NJ
C
      NN = N
      NX = NI
      NY = NJ
      IFL = IFLGG
      IF (NX .LT. 1  .OR.  NX .GT. NROW  .OR.  NY .LT. 1
     .   .OR.  IFL .LT. 0  .OR.  IFL .GT. 2) GO TO 4
      IST = NST
      IF (IFL .EQ. 2) THEN
C
C Compute gradient estimates at the nodes.
C
        DO 1 I = 1,NN
          CALL GRADL (NN,I,X,Y,Z,F,LIST,LPTR,
     .                LEND, GRAD(1,I),IERR)
          IF (IERR .LT. 0) GO TO 5
    1     CONTINUE
        IFL = 1
      ENDIF
C
C Compute uniform grid points and interpolated values.
C
      NEX = 0
      DO 3 J = 1,NY
        DO 2 I = 1,NX
          CALL INTRC1 (NN,PLAT(I),PLON(J),X,Y,Z,F,LIST,LPTR,
     .                 LEND,IFLGS,SIGMA,IFL,
     .                 GRAD, IST, FF(I,J),IERR)
          IF (IERR .LT. 0) GO TO 5
          NEX = NEX + IERR
    2     CONTINUE
    3   CONTINUE
      IER = NEX
      RETURN
C
C NI, NJ, or IFLGG is outside its valid range.
C
    4 IER = -1
      RETURN
C
C Error in GRADL or INTRC1.
C
    5 IER = IERR
      RETURN
      END
      subroutine intrpc1_n(npts,nptso,olats,olons,x,y,z,datain,lst,
     .                     lptr,lend,odata,ierr)
      integer, intent(in) :: npts, nptso
      integer, intent(out) :: ierr
      double precision, intent(in), dimension(nptso) :: olats,olons
      double precision, intent(in), dimension(npts) :: datain,x,y,z
      double precision, intent(out), dimension(nptso) :: odata
      double precision, dimension(3,nptso) :: grad
      double precision, dimension(npts) :: sigma
      integer, intent(in), dimension(npts) :: lend
      integer, intent(in), dimension(6*(npts-2)) :: lst,lptr
      integer n,ierr1,ist
      ist = 1
      ierr = 0
      sigma = 0
      do n=1,nptso
         call intrc1(npts,olats(n),olons(n),x,y,z,datain,lst,lptr,
     .               lend,-1,sigma,0,grad,ist,odata(n),ierr1)
         if (ierr1 .ne. 0) then
           !print *,n,'warning: ierr = ',ierr1,' in intrc1_n'
           !print *,olats(n), olons(n), npts
           !stop
           ierr = ierr + ierr1
         endif
      enddo
      end subroutine intrpc1_n
      subroutine intrpc0_n(npts,nptso,olats,olons,x,y,z,datain,lst,
     .                     lptr,lend,odata,ierr)
      integer, intent(in) :: npts, nptso
      integer, intent(out) :: ierr
      double precision, intent(in), dimension(nptso) :: olats,olons
      double precision, intent(in), dimension(npts) :: datain,x,y,z
      double precision, intent(out), dimension(nptso) :: odata
      integer, intent(in), dimension(npts) :: lend
      integer, intent(in), dimension(6*(npts-2)) :: lst,lptr
      integer n,ierr1,ist
      ist = 1
      ierr = 0
      do n=1,nptso
         call intrc0(npts,olats(n),olons(n),x,y,z,datain,lst,lptr,
     .               lend,ist,odata(n),ierr1)
         if (ierr1 .ne. 0) then
           !print *,n,'warning: ierr = ',ierr1,' in intrc0_n'
           !print *,olats(n), olons(n), npts
           !stop
           ierr = ierr + ierr1
         endif
      enddo
      end subroutine intrpc0_n
      subroutine intrpnn_n(npts,nptso,olats,olons,x,y,z,datain,lst,
     .                     lptr,lend,odata,ierr)
      integer, intent(in) :: npts, nptso
      integer, intent(out) :: ierr
      double precision, intent(in), dimension(nptso) :: olats,olons
      double precision, intent(in), dimension(npts) :: datain,x,y,z
      double precision, intent(out), dimension(nptso) :: odata
      integer, intent(in), dimension(npts) :: lend
      integer, intent(in), dimension(6*(npts-2)) :: lst,lptr
      integer n,ierr1,ist
      ist = 1
      ierr = 0
      do n=1,nptso
         call intrnn(npts,olats(n),olons(n),x,y,z,datain,lst,lptr,
     .               lend,ist,odata(n),ierr1)
         if (ierr1 .ne. 0) then
           !print *,n,'warning: ierr = ',ierr1,' in intrc0_n'
           !print *,olats(n), olons(n), npts
           !stop
           ierr = ierr + ierr1
         endif
      enddo
      end subroutine intrpnn_n
