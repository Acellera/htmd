      SUBROUTINE TSTARC(SCALE, DENS, 
     1 NATOMS, COORD, IAN, NESP, POTPT, PTYPE, ITYPE,
     2 TESTARC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION COORD(3,*)
      DOUBLE PRECISION POTPT(3,*)
      INTEGER          PTYPE(*)
      DOUBLE PRECISION SCALE(*)
      DOUBLE PRECISION DENS(*)
      INTEGER          TESTARC(*)
      INTEGER NESP
      INCLUDE 'SIZES'
      INTEGER          IATOM,JATOM,KATOM, ISKIP
      DOUBLE PRECISION RADIUS,ANGLE, RSKIP

      COMMON /WORK1/  PAD1(2*MESP), RAD(MESP), IAS(MESP*2)
      COMMON /POTESP/ XC,YC,ZC,ESPNUC,ESPELE
      COMMON /THRESH/ RRIONS, RRATMS
      DIMENSION VANREAD(NUMATM)
      DIMENSION VANDER(100)
      DIMENSION CON(3,10000)
      DIMENSION IAN(NATOMS)
      DIMENSION INBR(2000),CNBR(3,2000),RNBR(2000)
      LOGICAL SNBR(2000),MNBR(2000)
      DIMENSION CI(3), IELDAT(56), TEMP(3)
      DIMENSION CW(3,2)
      LOGICAL SI
      LOGICAL COLLID
      DOUBLE PRECISION A11,A21,A31, A12,A22,A32, A13,A23,A33
      DOUBLE PRECISION XP(3),YP(3),ZP(3), JK(3), LENGTH
C
C     DATA FOR VANDER VALL RADII
C
      CHARACTER MARKER*3, MARKSS*3, MYNAM*3, IELDAT*4
      DATA VANDER/1.20D0,1.20D0,1.37D0,1.45D0,1.45D0,1.50D0,1.50D0,
     1            1.40D0,1.35D0,1.30D0,1.57D0,1.36D0,1.24D0,1.17D0,
     2            1.80D0,1.75D0,1.70D0,17*0.0D0,2.3D0,65*0.0D0/
      DATA MARKER/'A  '/,MARKSS/'SS0'/,MYNAM/'UC '/
C
      DATA IELDAT/'  BQ','  H ','  HE','  LI','  BE','  B ',
     1            '  C ','  N ','  O ','  F ','  NE','  NA',
     2            '  MG','  AL','  SI','  P ','  S ','  CL',
     3            '  AR','  K ','  CA','  SC','  TI','  V ',
     4            '  CR','  MN','  FE','  CO','  NI','  CU',
     5            '  ZN','  GA','  GE','  AS','  SE','  BR',
     6            '  KR','  RB','  SR','   Y','  ZR','  NB',
     7            '  MO','  TC','  RU','  RH','  PD','  AG',
     8            '  CD','  IN','  SN','  SB','  TE','   I',
     9            '   X','  CS'/

      IF (DENS(1) .LT. 0.02d0) RETURN
      PI = 4.D0*ATAN(1.D0) / 180.0D0
C
C    J--->I (XP).... P
C     \
C (YP) K
C
      IATOM = TESTARC(1)
      JATOM = TESTARC(2)
      KATOM = TESTARC(3)
      ANGLE = FLOAT(TESTARC(4)) * PI
C obtain XP axis from JI vector and normilize XP vector
      LENGTH = SQRT(DIST2(COORD(1,IATOM),COORD(1,JATOM)))
      XP(1) = (COORD(1,IATOM) - COORD(1,JATOM)) / LENGTH
      XP(2) = (COORD(2,IATOM) - COORD(2,JATOM)) / LENGTH
      XP(3) = (COORD(3,IATOM) - COORD(3,JATOM)) / LENGTH
C normilize JK vector prior collinearity test
      LENGTH = SQRT(DIST2(COORD(1,KATOM),COORD(1,JATOM)))
      JK(1) = (COORD(1,KATOM) - COORD(1,JATOM)) / LENGTH
      JK(2) = (COORD(2,KATOM) - COORD(2,JATOM)) / LENGTH
      JK(3) = (COORD(3,KATOM) - COORD(3,JATOM)) / LENGTH
C test for collinear vectors
      IF(ABS(XP(1)*JK(1)+XP(2)*JK(2)+XP(3)*JK(3)) .GT. 0.99) THEN
         WRITE(*,*) 'Error: reference atoms are on stright line ',
     $              IATOM,JATOM,KATOM
         STOP
      ENDIF
C obtain ZP axis as ZP = XP x JK and normalize it
      ZP(1) = XP(2)*JK(3) - XP(3)*JK(2)
      ZP(2) = XP(3)*JK(1) - XP(1)*JK(3)
      ZP(3) = XP(1)*JK(2) - XP(2)*JK(1)
      LENGTH = SQRT(ZP(1)*ZP(1) + ZP(2)*ZP(2) + ZP(3)*ZP(3))
      ZP(1) = ZP(1) / LENGTH
      ZP(2) = ZP(2) / LENGTH
      ZP(3) = ZP(3) / LENGTH
C obtain YP axis as YP = ZP x XP
      YP(1) = ZP(2)*XP(3) - ZP(3)*XP(2)
      YP(2) = ZP(3)*XP(1) - ZP(1)*XP(3)
      YP(3) = ZP(1)*XP(2) - ZP(2)*XP(1)
C rotate the YP axis in YP:ZP plane on ANGLE degrees
      YP(1) = YP(1)*COS(ANGLE) + ZP(1)*SIN(ANGLE)
      YP(2) = YP(2)*COS(ANGLE) + ZP(2)*SIN(ANGLE)
      YP(3) = YP(3)*COS(ANGLE) + ZP(3)*SIN(ANGLE)
C get new ZP axis as ZP = XP x YP
      ZP(1) = XP(2)*YP(3) - XP(3)*YP(2)
      ZP(2) = XP(3)*YP(1) - XP(1)*YP(3)
      ZP(3) = XP(1)*YP(2) - XP(2)*YP(1)
C Cosines between two coordinate systems: X,Y,Z and XP,YP,ZP
      A11 = XP(1)
      A21 = YP(1)
      A31 = ZP(1)
        A12 = XP(2)
        A22 = YP(2)
        A32 = ZP(2)
          A13 = XP(3)
          A23 = YP(3)
          A33 = ZP(3)
C radius of the arc
      RADIUS = VANDER(IAN(IATOM)) * SCALE(IATOM)

C determine density of points to be generated
      RSKIP=20.0D0/DENS(IATOM)
      ISKIP=RSKIP
      IF(FLOAT(ISKIP)+0.5D0.LT.RSKIP) ISKIP=ISKIP+1
C generate points on the arc drawn around I-atom; the arc plane is parallel to XY plane
      DO I=0,359,ISKIP
C generate a point in the local frame
        ANGLE = FLOAT(I) * PI
        TEMP(1) = COS(ANGLE) * RADIUS
        TEMP(2) = SIN(ANGLE) * RADIUS
        TEMP(3) = 0.0D0
C apply rotation
        JK(1) = TEMP(1)*A11 + TEMP(2)*A21 + TEMP(3)*A31
        JK(2) = TEMP(1)*A12 + TEMP(2)*A22 + TEMP(3)*A32
        JK(3) = TEMP(1)*A13 + TEMP(2)*A23 + TEMP(3)*A33
C translation
        TEMP(1) = JK(1) + COORD(1,IATOM)
        TEMP(2) = JK(2) + COORD(2,IATOM)
        TEMP(3) = JK(3) + COORD(3,IATOM)
C check distance to atoms
        DO K=1,NATOMS
          IF(K.NE.IATOM.AND.DIST2(TEMP,COORD(1,K)).LT.RRATMS) GOTO 100
        ENDDO
C save the point
        CALL ADDPT(POTPT,NESP,TEMP,PTYPE,ITYPE)
  100   CONTINUE
      ENDDO       
      RETURN
      END


      SUBROUTINE SURFAC(SCALE, DENS, 
     1 NATOMS, COORD, IAN, NESP, POTPT, PTYPE, ITYPE,
     2 Exclude,Ini_Flag)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C***********************************************************************
C
C      THIS SUBROUTINE CALCULATES THE MOLECULAR SURFACE OF A MOLECULE
C      GIVEN THE COORDINATES OF ITS ATOMS.  VAN DER WAALS' RADII FOR
C      THE ATOMS AND THE PROBE RADIUS MUST ALSO BE SPECIFIED.
C
C      ON INPUT    SCALE = INITIAL VAN DER WAALS' SCALE FACTOR
C                  DENS  = DENSITY OF POINTS PER UNIT AREA
C
C      THIS SUBROUTINE WAS LIFTED FROM MICHAEL CONNOLLY'S SURFACE
C      PROGRAM FOR UCSF GRAPHICS SYSTEM BY U.CHANDRA SINGH AND
C      P.A.KOLLMAN AND MODIFIED FOR USE IN QUEST. K.M.MERZ
C      ADAPTED AND CLEANED UP THIS PROGRAM FOR USE IN AMPAC/MOPAC
C      IN FEB. 1989 AT UCSF.
C
C***********************************************************************
      
      DOUBLE PRECISION COORD(3,*)
      DOUBLE PRECISION POTPT(3,*)
      INTEGER          PTYPE(*)
      DOUBLE PRECISION SCALE(*)
      DOUBLE PRECISION DENS(*)
      INTEGER          Exclude(*)
      INTEGER NESP
	LOGICAL Ini_Flag
      INCLUDE 'SIZES'

      COMMON /WORK1/  PAD1(2*MESP), RAD(MESP), IAS(MESP*2)
      COMMON /POTESP/ XC,YC,ZC,ESPNUC,ESPELE
      COMMON /THRESH/ RRIONS, RRATMS

      DIMENSION VANREAD(NATOMS)
      DIMENSION VANDER(100)
      DIMENSION CON(3,10000)
      DIMENSION IAN(NATOMS)
C
C     NEIGHBOR ARRAYS
C
C     THIS SAME DIMENSION FOR THE MAXIMUM NUMBER OF NEIGHBORS
C     IS USED TO DIMENSION ARRAYS IN THE LOGICAL FUNCTION COLLID
C
      DIMENSION INBR(2000),CNBR(3,2000),RNBR(2000)
      LOGICAL SNBR(2000),MNBR(2000)
C
C     ARRAYS FOR ALL ATOMS
C
C     IATOM, JATOM AND KATOM COORDINATES
C
      DIMENSION CI(3), IELDAT(56), TEMP(3)
C
C     GEOMETRIC CONSTRUCTION VECTORS
C
      DIMENSION CW(3,2)
C
C     LOGICAL VARIABLES
C
      LOGICAL SI
C
C     LOGICAL FUNCTIONS
C
      LOGICAL COLLID
C
C     DATA FOR VANDER VALL RADII
C
      CHARACTER MARKER*3, MARKSS*3, MYNAM*3, IELDAT*4
      DATA VANDER/1.20D0,1.20D0,1.37D0,1.45D0,1.45D0,1.50D0,1.50D0,
     1            1.40D0,1.35D0,1.30D0,1.57D0,1.36D0,1.24D0,1.17D0,
     2            1.80D0,1.75D0,1.70D0,17*0.0D0,2.3D0,65*0.0D0/
      DATA MARKER/'A  '/,MARKSS/'SS0'/,MYNAM/'UC '/
C
      DATA IELDAT/'  BQ','  H ','  HE','  LI','  BE','  B ',
     1            '  C ','  N ','  O ','  F ','  NE','  NA',
     2            '  MG','  AL','  SI','  P ','  S ','  CL',
     3            '  AR','  K ','  CA','  SC','  TI','  V ',
     4            '  CR','  MN','  FE','  CO','  NI','  CU',
     5            '  ZN','  GA','  GE','  AS','  SE','  BR',
     6            '  KR','  RB','  SR','   Y','  ZR','  NB',
     7            '  MO','  TC','  RU','  RH','  PD','  AG',
     8            '  CD','  IN','  SN','  SB','  TE','   I',
     9            '   X','  CS'/

C     check if van der waals radii file is saved in current file
C     format is one radius per line in same order as atoms
      LOGICAL :: VANDER_EXISTS
      
      INQUIRE(FILE="vander.txt", EXIST=VANDER_EXISTS)
      IF (VANDER_EXISTS .eqv. .TRUE.) THEN
         OPEN(UNIT=91, FILE='vander.txt', 
     1        FORM='FORMATTED', STATUS='OLD')
         DO 20 I=1,NATOMS
            READ(UNIT=91,FMT=*) VANREAD(I)
            write(*,*) I, VANREAD(I)
 20      CONTINUE
         close(91)
      ENDIF

      IF (DENS(1) .LT. 0.02d0) RETURN

      PI=4.D0*ATAN(1.D0)
C     INSERT VAN DER WAAL RADII FOR ZINC
      VANDER(30)=1.00D0
C
C     CONVERT INTERNAL TO CARTESIAN COORDINATES
C
C     ONLY VAN DER WAALS' TYPE SURFACE IS GENERATED
C
      DO 30 I=1,NATOMS
         RAD(I) = VANDER(IAN(I))*SCALE(I)
C     replace with values that were read in if present
         IF(VANDER_EXISTS .EQV. .TRUE.) THEN
            RAD(I) = VANREAD(I) * SCALE(I)
         ENDIF
         IF (RAD(I) .LT. 0.01D0) THEN
            WRITE(IW,'(T2,''VAN DER WAALS'''' RADIUS FOR ATOM '',I3,
     1         '' IS ZERO)'' )') I
            STOP
         ENDIF
         IAS(I) = 2
   30 CONTINUE
C
C     BIG LOOP FOR EACH ATOM
C
      DO 110 IATOM = 1, NATOMS
         IF (IAS(IATOM) .EQ. 0) GO TO 110
C
C     TRANSFER VALUES FROM LARGE ARRAYS TO IATOM VARIABLES
C
         RI = RAD(IATOM)
         SI = IAS(IATOM) .EQ. 2
         DO 40 K = 1,3
            CI(K) = COORD(K,IATOM)
   40    CONTINUE
C
C     GATHER THE NEIGHBORING ATOMS OF IATOM
C
         NNBR = 0
         DO 60 JATOM = 1, NATOMS
            IF (IATOM .EQ. JATOM .OR. IAS(JATOM) .EQ. 0) GO TO 60
            D2 = DIST2(CI,COORD(1,JATOM))
            IF (D2 .GE. (2*RW+RI+RAD(JATOM)) ** 2) GO TO 60
C
C     WE HAVE A NEW NEIGHBOR
C     TRANSFER ATOM COORDINATES, RADIUS AND SURFACE REQUEST NUMBER
C
            NNBR = NNBR + 1
            IF (NNBR .GT. 2000)THEN
               WRITE(IW,'(''ERROR'',2X,''TOO MANY NEIGHBORS:'',I5)')NNBR
               STOP
            ENDIF
            INBR(NNBR) = JATOM
            DO 50 K = 1,3
               CNBR(K,NNBR) = COORD(K,JATOM)
   50       CONTINUE
            RNBR(NNBR) = RAD(JATOM)
            SNBR(NNBR) = IAS(JATOM) .EQ. 2
   60    CONTINUE
C
C     CONTACT SURFACE
C
         IF (.NOT. SI) GO TO 110
         NCON = (4 * PI * RI ** 2) * DENS(IATOM)
         IF (NCON .GT. 1000) NCON = 1000
C         IF(NCON .LT. 50) NCON = 50
C
C     THIS CALL MAY DECREASE NCON SOMEWHAT
C
         IF ( NCON .EQ. 0) THEN
            WRITE(IW,'(T2,''VECTOR LENGTH OF ZERO IN SURFAC'')')
            STOP
         ENDIF
         CALL GENUN(CON,NCON)
C         CALL GENUN(CON(1,INN+1),NCON)
         
         AREA = (4 * PI * RI ** 2) / NCON
C
C     CONTACT PROBE PLACEMENT LOOP
C
         DO 100 I = 1, NCON
            DO 70 K = 1,3
               CW(K,1) = CI(K) + (RI + RW) * CON(K,I)
   70       CONTINUE
C
C     CHECK FOR COLLISION WITH NEIGHBORING ATOMS
C
            IF (COLLID(CW(1,1),RW,CNBR,
     1      RNBR,MNBR,NNBR,1,JNBR,KNBR)) GO TO 100
            DO 80 KK=1,3
               TEMP(KK) =CI(KK)+RI*CON(KK,I)
   80       CONTINUE
* check proximity to existing points
            if(ITYPE .ne. 10) then
C              this is a perturbation point
               DO K = 2, Exclude(1)+1
C                 PERTURBATION point - ATOM
                  IF(DIST2(TEMP,COORD(1,Exclude(K))).LT.RRATMS) GOTO 100
               ENDDO
               DO K = 1, NESP
                  IF(PTYPE(K).ne.10) then
C                   PERTURBATION point - PERTURBATION point
                    IF(DIST2(TEMP,POTPT(1,K)) .LT. RRIONS) GOTO 100
                  ENDIF
               ENDDO
            else
C              this is a grid point
               DO K = 2, Exclude(1)+1
C                 PERTURBATION point - "Exclude(K)" atom
                  IF(DIST2(TEMP,COORD(1,Exclude(K))).LT.RRATMS) GOTO 100
               ENDDO
               IF(.NOT. Ini_Flag) then
                  DO K = 1, NESP
                     IF(PTYPE(K).eq.10) then
C                       ESP grid point - ESP grid point
                        IF(DIST2(TEMP,POTPT(1,K)) .LT. RRIONS) GOTO 100
                     ENDIF
                  ENDDO
               ENDIF
            endif
C
C     PROXIMITY TEST IS PASSED; STORE POINT IN POTPT AND INCREMENT NESP VARIABLE
C
            CALL ADDPT(POTPT,NESP,TEMP,PTYPE,ITYPE)

  100    CONTINUE
  110 CONTINUE
      RETURN
      END


C****************************************************************
      LOGICAL FUNCTION COLLID(CW,RW,CNBR,RNBR,MNBR,NNBR,ISHAPE,
     1JNBR,KNBR)
C****************************************************************
C
C     COLLISION CHECK OF PROBE WITH NEIGHBORING ATOMS
C     USED BY SURFAC ONLY.
C
C****************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CW(3)
      DIMENSION CNBR(3,2000)
      DIMENSION RNBR(2000)
      LOGICAL MNBR(2000)
      IF (NNBR .LE. 0) GO TO 20
C
C     CHECK WHETHER PROBE IS TOO CLOSE TO ANY NEIGHBOR
C
      DO 10 I = 1, NNBR
         IF (ISHAPE .GT. 1 .AND. I .EQ. JNBR) GO TO 10
         IF (ISHAPE .EQ. 3 .AND. (I .EQ. KNBR .OR. .NOT. MNBR(I)))
     1   GO TO 10
         SUMRAD = RW + RNBR(I)
         VECT1 = DABS(CW(1) - CNBR(1,I))
         IF (VECT1 .GE. SUMRAD) GO TO 10
         VECT2 = DABS(CW(2) - CNBR(2,I))
         IF (VECT2 .GE. SUMRAD) GO TO 10
         VECT3 = DABS(CW(3) - CNBR(3,I))
         IF (VECT3 .GE. SUMRAD) GO TO 10
         SR2 = SUMRAD ** 2
         DD2 = VECT1 ** 2 + VECT2 ** 2 + VECT3 ** 2
         IF (DD2 .LT. SR2) GO TO 30
   10 CONTINUE
   20 CONTINUE
      COLLID = .FALSE.
      GO TO 40
   30 CONTINUE
      COLLID = .TRUE.
   40 CONTINUE
      RETURN
      END


C
C     DETERMINE DISTANCE BETWEEN COORDINATE POINTS
C
      DOUBLE PRECISION FUNCTION DIST2(A,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(3),B(3)
      DIST2 = (A(1)-B(1))**2 + (A(2)-B(2))**2 + (A(3)-B(3))**2
      RETURN
      END


      SUBROUTINE GENUN(U,N)
C****************************************************************
C
C     GENERATE UNIT VECTORS OVER SPHERE. USED BY SURFAC ONLY.
C
C****************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION :: U(3,10000)
C      DIMENSION U(3,10000)
      INTEGER :: I, J, NU
      
      DO I=1,10000
         U(1,I)=0.0
         U(2,I)=0.0
         U(3,I)=0.0
      ENDDO

      PI=4.D0*ATAN(1.D0)
      NEQUAT = SQRT(N * PI)
      NVERT = NEQUAT/2
      NU = 0
      write(*,*) 'genun'
      write(*,*) N
      DO 20 I = 1,NVERT+1
         FI = (PI * (I-1)) / NVERT
         Z = COS(FI)
         XY = SIN(FI)
         NHOR = NEQUAT * XY
         IF (NHOR .LT. 1) NHOR = 1
         DO 10 J = 1,NHOR
            FJ = (2.D0 * PI * (J-1)) / NHOR
            X = DCOS(FJ) * XY
            Y = DSIN(FJ) * XY
            IF (NU .GE. N) GO TO 30
            NU = NU + 1
            write(*,*) NU
            U(1,NU) = X
            U(2,NU) = Y
            U(3,NU) = Z
   10    CONTINUE
   20 CONTINUE
   30 CONTINUE
      N = NU
      RETURN
      END




      SUBROUTINE PBONDS(SCALE, DENS, 
     1 NATOMS, COORD, IAN, NESP, POTPT, PTYPE, ITYPE, GRID, NGRID,
     2 Exclude)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DOUBLE PRECISION COORD(3,*)
      DOUBLE PRECISION POTPT(3,*), GRID(3,NGRID)
      INTEGER          PTYPE(*)
      DOUBLE PRECISION SCALE(*)
      DOUBLE PRECISION DENS(*)
      INTEGER          Exclude(*)
      INTEGER NESP
      INCLUDE 'SIZES'

      COMMON /DEBUG/ NBONDS,IPOINT,IAT,JAT,DEVMAX
      COMMON /THRESH/ RRIONS, RRATMS

      DIMENSION RAD(100)
      DIMENSION VANDER(100)
      DIMENSION IAN(NATOMS)
C
      DIMENSION CI(3), IELDAT(56), TEMP(3), CNEAR(3)
C
C     DATA FOR VANDER VALL RADII
C
      CHARACTER MARKER*3, MARKSS*3, MYNAM*3, IELDAT*4
      DATA VANDER/1.20D0,1.20D0,1.37D0,1.45D0,1.45D0,1.50D0,1.50D0,
     1            1.40D0,1.35D0,1.30D0,1.57D0,1.36D0,1.24D0,1.17D0,
     2            1.80D0,1.75D0,1.70D0,17*0.0D0,2.3D0,65*0.0D0/
      DATA MARKER/'A  '/,MARKSS/'SS0'/,MYNAM/'UC '/
C
      DATA IELDAT/'  BQ','  H ','  HE','  LI','  BE','  B ',
     1            '  C ','  N ','  O ','  F ','  NE','  NA',
     2            '  MG','  AL','  SI','  P ','  S ','  CL',
     3            '  AR','  K ','  CA','  SC','  TI','  V ',
     4            '  CR','  MN','  FE','  CO','  NI','  CU',
     5            '  ZN','  GA','  GE','  AS','  SE','  BR',
     6            '  KR','  RB','  SR','   Y','  ZR','  NB',
     7            '  MO','  TC','  RU','  RH','  PD','  AG',
     8            '  CD','  IN','  SN','  SB','  TE','   I',
     9            '   X','  CS'/
      PI=4.D0*ATAN(1.D0)
      NBONDS = 0
      NESPOLD = NESP
C
      DO 30 I=1,NATOMS
         RAD(I) = VANDER(IAN(I))*SCALE(I)
         IF (RAD(I) .LT. 0.01D0) THEN
            WRITE(IW,'(T2,''VAN DER WAALS'''' RADIUS FOR ATOM '',I3,
     1         '' IS ZERO)'' )')
            STOP
         ENDIF
   30 CONTINUE
C
C     BIG LOOP FOR EACH ATOM
C
      DO 110 IATOM = 1, NATOMS
C
C     Prepare IATOM data for fast access
C
         RI = RAD(IATOM)
         DO 40 K = 1,3
            CI(K) = COORD(K,IATOM)
   40    CONTINUE
C
C     GATHER THE NEIGHBORING ATOMS OF IATOM
C
         RADI = VANDER(IAN(IATOM))
         DO 60 JATOM = 1, NATOMS
C  Skip excluded atoms
            DO K = 2, Exclude(1)+1
               IF(JATOM.EQ.Exclude(K)) GOTO 60
            ENDDO
C
            IF (IATOM .EQ. JATOM ) GO TO 60
            D2 = SQRT(DIST2(CI,COORD(1,JATOM)))
            IF( D2*2.0d0 .LE. RADI+VANDER(IAN(JATOM))+0.8d0) THEN
C
C  Set new point along chemical bond
               RIJ = (D2 + RAD(JATOM)) / D2
               TEMP(1) = CI(1) + (COORD(1,JATOM) - CI(1)) * RIJ
               TEMP(2) = CI(2) + (COORD(2,JATOM) - CI(2)) * RIJ
               TEMP(3) = CI(3) + (COORD(3,JATOM) - CI(3)) * RIJ
C
               DO K = 1, NATOMS
C                 PERTURBATION point - "Exclude(K)" atom
                  IF(DIST2(TEMP,COORD(1,K)).LT.RRATMS) THEN
                    GOTO 60
                  ENDIF
               ENDDO
C
C              WE HAVE A BOND
               NBONDS = NBONDS + 1
C
               CALL ONSURF(CI,TEMP,COORD,NATOMS,GRID,NGRID,
     1                        POTPT,NESP,PTYPE,ITYPE,
     2                        IATOM,JATOM)
            ENDIF
   60    CONTINUE
  110 CONTINUE
C
      NBONDS = NBONDS / 2
      RETURN
      END


      SUBROUTINE PANGLE(SCALE, DENS, 
     1 NATOMS, COORD, IAN, NESP, POTPT, PTYPE, ITYPE, GRID, NGRID)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DOUBLE PRECISION COORD(3,*)
      DOUBLE PRECISION POTPT(3,*), GRID(3,NGRID)
      INTEGER          PTYPE(*)
      DOUBLE PRECISION SCALE(*)
      DOUBLE PRECISION DENS(*)
      INTEGER NESP
      INCLUDE 'SIZES'

      COMMON /DEBUG/ NBONDS,IPOINT,IAT,JAT,DEVMAX
      COMMON /THRESH/ RRIONS, RRATMS

      DIMENSION RAD(100)
      DIMENSION VANDER(100)
      DIMENSION IAN(NATOMS)
C
      DIMENSION CI(3), IELDAT(56), TEMP(3), CNEAR(3)
C
C     DATA FOR VANDER VALL RADII
C
      CHARACTER MARKER*3, MARKSS*3, MYNAM*3, IELDAT*4
      DATA VANDER/1.20D0,1.20D0,1.37D0,1.45D0,1.45D0,1.50D0,1.50D0,
     1            1.40D0,1.35D0,1.30D0,1.57D0,1.36D0,1.24D0,1.17D0,
     2            1.80D0,1.75D0,1.70D0,17*0.0D0,2.3D0,65*0.0D0/
      DATA MARKER/'A  '/,MARKSS/'SS0'/,MYNAM/'UC '/
C
      DATA IELDAT/'  BQ','  H ','  HE','  LI','  BE','  B ',
     1            '  C ','  N ','  O ','  F ','  NE','  NA',
     2            '  MG','  AL','  SI','  P ','  S ','  CL',
     3            '  AR','  K ','  CA','  SC','  TI','  V ',
     4            '  CR','  MN','  FE','  CO','  NI','  CU',
     5            '  ZN','  GA','  GE','  AS','  SE','  BR',
     6            '  KR','  RB','  SR','   Y','  ZR','  NB',
     7            '  MO','  TC','  RU','  RH','  PD','  AG',
     8            '  CD','  IN','  SN','  SB','  TE','   I',
     9            '   X','  CS'/
      PI=4.D0*ATAN(1.D0)
      NESPOLD = NESP
C
      DO 30 I=1,NATOMS
         RAD(I) = VANDER(IAN(I))*SCALE(I)
         IF (RAD(I) .LT. 0.01D0) THEN
            WRITE(IW,'(T2,''VAN DER WAALS'''' RADIUS FOR ATOM '',I3,
     1         '' IS ZERO)'' )')
            STOP
         ENDIF
   30 CONTINUE
C
C     BIG LOOP FOR EACH ATOM
C
      DO 110 IATOM = 1, NATOMS
C
C Ignore H-atoms
C
         IF(IAN(IATOM).EQ.1) GOTO 110
C
C     Prepare IATOM data for fast access
C
         RI = RAD(IATOM)
         DO 40 K = 1,3
            CI(K) = COORD(K,IATOM)
   40    CONTINUE
C
C     GATHER THE NEIGHBORING ATOMS OF IATOM
C
         RADI = VANDER(IAN(IATOM))
         NBONDS = 0
         DO 60 JATOM = 1, NATOMS
            IF (IATOM .EQ. JATOM ) GO TO 60
            D2 = SQRT(DIST2(CI,COORD(1,JATOM)))
            IF( D2*2.0d0 .LE. RADI+VANDER(IAN(JATOM))+0.8d0) THEN
C
C              WE HAVE A BOND
C
               NBONDS = NBONDS + 1
               IF(NBONDS.EQ.1) IB = JATOM
               IF(NBONDS.EQ.2) IC = JATOM
               IF(NBONDS.GE.3) GOTO 110
            ENDIF
   60    CONTINUE
C
         IF(NBONDS.EQ.2) THEN
C
C           We have got IATOM having two neighbours only
C
            AB = SQRT(DIST2(CI,COORD(1,IB)))
            AC = SQRT(DIST2(CI,COORD(1,IC)))
C Vector AB
            BX = (COORD(1,IB) - CI(1)) / AB
            BY = (COORD(2,IB) - CI(2)) / AB
            BZ = (COORD(3,IB) - CI(3)) / AB
C Vector AC
            CX = (COORD(1,IC) - CI(1)) / AC
            CY = (COORD(2,IC) - CI(2)) / AC
            CZ = (COORD(3,IC) - CI(3)) / AC
C
C Check collinearity of AB and AC vectors
C
            COLLIN = ABS(BX*CX + BY*CY + BZ*CZ)
            IF(COLLIN .LT. 0.0001d0) THEN
               write(*,"(A,2I3,A,2I3,A)") ' PROGRAM ERROR: VECTORS ',
     *            IATOM,IB,' AND ',IATOM,IC, ' ARE COLINEAR'
               stop
            ENDIF
C
C Vector A = AB + AC
            AX = BX + CX
            AY = BY + CY
            AZ = BZ + CZ
            ALEN = SQRT(AX*AX + AY*AY + AZ*AZ)
C
            CNEAR(1) = CI(1) - (BX + CX) * RI / ALEN
            CNEAR(2) = CI(2) - (BY + CY) * RI / ALEN
            CNEAR(3) = CI(3) - (BZ + CZ) * RI / ALEN
* add point
            CALL ONSURF(CI,CNEAR,COORD,NATOMS,GRID,NGRID,
     1                        POTPT,NESP,PTYPE,ITYPE,
     2                        IATOM,IATOM)
C
C Add points out of plane
C
            DX = BY*CZ - BZ*CY
            DY = BZ*CX - BX*CZ
            DZ = BX*CY - BY*CX
C Z-direction: Positive
            CNEAR(1) = CI(1) + DX * RI
            CNEAR(2) = CI(2) + DY * RI
            CNEAR(3) = CI(3) + DZ * RI 
* add point
            CALL ONSURF(CI,CNEAR,COORD,NATOMS,GRID,NGRID,
     1                        POTPT,NESP,PTYPE,ITYPE,
     2                        IATOM,IATOM)
C Z-direction: Negative
            CNEAR(1) = CI(1) - DX * RI
            CNEAR(2) = CI(2) - DY * RI
            CNEAR(3) = CI(3) - DZ * RI 
* add point
            CALL ONSURF(CI,CNEAR,COORD,NATOMS,GRID,NGRID,
     1                        POTPT,NESP,PTYPE,ITYPE,
     2                        IATOM,IATOM)
C Diagonal: Positive direction
            PX = DX - AX
            PY = DY - AY
            PZ = DZ - AZ
            PLEN = SQRT(PX*PX + PY*PY + PZ*PZ)
            CNEAR(1) = CI(1) + PX * RI / PLEN
            CNEAR(2) = CI(2) + PY * RI / PLEN
            CNEAR(3) = CI(3) + PZ * RI / PLEN
* add point
            CALL ONSURF(CI,CNEAR,COORD,NATOMS,GRID,NGRID,
     1                        POTPT,NESP,PTYPE,ITYPE,
     2                        IATOM,IATOM)
C Diagonal: Negative direction
            PX = DX + AX
            PY = DY + AY
            PZ = DZ + AZ
            PLEN = SQRT(PX*PX + PY*PY + PZ*PZ)
            CNEAR(1) = CI(1) - PX * RI / PLEN
            CNEAR(2) = CI(2) - PY * RI / PLEN
            CNEAR(3) = CI(3) - PZ * RI / PLEN
* add point
            CALL ONSURF(CI,CNEAR,COORD,NATOMS,GRID,NGRID,
     1                        POTPT,NESP,PTYPE,ITYPE,
     2                        IATOM,IATOM)
         ENDIF
  110 CONTINUE
C
      RETURN
      END


      SUBROUTINE ADDPT(POOL,NESP,POINT,PTYPE,ITYPE)
      DOUBLE PRECISION POOL(3,*), POINT(3)
      INTEGER          NESP, PTYPE(*), ITYPE
      INCLUDE 'SIZES'
      NESP = NESP + 1
      IF(NESP.LE.MESP) THEN
         POOL(1,NESP) = POINT(1)
         POOL(2,NESP) = POINT(2)
         POOL(3,NESP) = POINT(3)
         PTYPE(NESP)   = ITYPE
      ELSE
         WRITE(*,*) " INSUFFICIENT <MESP> STORAGE = ", MESP
         STOP
      ENDIF
      RETURN
      END


C
C  Search point CPOINT = CI + CVEC along vector CVEC 
C     on (+-)WIDTH% of its length untill the closest 
C     point on the Connolly surface (GRID) is identified
C  Take the closest point from the grid surface so the points 
C     generated are always on the Connolly surface
C  Add point to POTPT; increment NESP
C
      SUBROUTINE ONSURF(CI,CVEC,COORD,NATOMS,GRID,NGRID,
     1                        POTPT,NESP,PTYPE,ITYPE,
     2                        IATOM,JATOM)
C
C  CI       Cartesian coordinates of IATOM
C  CVEC     Cartesian coordinates of vector
C              CI + CVEC give coordinate of the search point
C  GRID     Cartesian coordinates of grid surface points
C
C  IATOM and JATOM are used for debugging only
C
      IMPLICIT NONE
      INTEGER          NATOMS, NGRID
      DOUBLE PRECISION CI(3), CVEC(3), COORD(3,NATOMS), GRID(3,NGRID),
     1                 POTPT(3,*)
      INTEGER          NESP, PTYPE(*), ITYPE, IATOM, JATOM
C
C Common variables:
C
      DOUBLE PRECISION RRIONS, RRATMS
      COMMON /THRESH/ RRIONS, RRATMS
      INTEGER          NBONDS,IPOINT,IAT,JAT
      DOUBLE PRECISION DEVMAX
      COMMON /DEBUG/ NBONDS,IPOINT,IAT,JAT,DEVMAX
C
C  Functions:
C
      DOUBLE PRECISION DIST2
C
C  Local variables:
C
C  NTRY     Number of search trials
C  WIDTH    Half width of search, %; will search both directions (+-)
C  VLEN     Vector length multiplier
C  VINCR    Length increment
C  RMIN     Shortest distance between grid and the search point
C  DGRID    Distance to grid point
C  INCR,K   Loop index
C  TEMP     Coordinates of the search point
C  CPOINT   Search point
C  CVEC     Radius-vector to point CVEC
C
      INTEGER          NTRY, INCR, K
      DOUBLE PRECISION WIDTH, VLEN, RMIN, DGRID, TEMP(3), CPOINT(3),
     1                 VINCR, RI
C
C Debugging variables
C
C  DEVIAT   Angle between initial and obtained vectors (numerical error)
C  VX,VY,VZ Resulting vector
      DOUBLE PRECISION VX,VY,VZ,DEVIAT
C
      DATA NTRY /50/, WIDTH /0.3d0/
C
C Check initial distance to atoms and eliminate bad points from begin
*      CALL ADDPT(POTPT,NESP,CVEC,PTYPE,ITYPE)
      DO K = 1, NATOMS
         IF(DIST2(CVEC,COORD(1,K)) .LT. RRATMS) RETURN
      ENDDO
C
C Find nearest grid surface point 
C
C Set increment value for NTRY/2 steps 
      VINCR = WIDTH *2.0d0 / FLOAT(NTRY)
C Set initial value of VLEN for the search
      VLEN = 1.0d0 - WIDTH
C Convert CVEC to a vector along line IATOM->CVEC
      CVEC(1) = CVEC(1) - CI(1)
      CVEC(2) = CVEC(2) - CI(2)
      CVEC(3) = CVEC(3) - CI(3)
C
      RMIN  = 100.0d0
      DO 10 INCR = 1, NTRY
         VLEN = VLEN + VINCR
         TEMP(1) = CI(1) + CVEC(1) * VLEN
         TEMP(2) = CI(2) + CVEC(2) * VLEN
         TEMP(3) = CI(3) + CVEC(3) * VLEN
*         CALL ADDPT(POTPT,NESP,TEMP,PTYPE,ITYPE)
C Check distance to atoms; remove points, which are too close to atoms
         DO K = 1, NATOMS
            IF(DIST2(TEMP,COORD(1,K)) .LT. RRATMS) GOTO 10
         ENDDO
C Find nearest grid surface point 
         DO K = 1, NGRID
            DGRID = DIST2(TEMP,GRID(1,K))
            IF(DGRID .LT. RMIN) THEN
               CPOINT(1) = GRID(1,K)
               CPOINT(2) = GRID(2,K)
               CPOINT(3) = GRID(3,K)
               RMIN = DGRID
            ENDIF
         ENDDO
   10 CONTINUE 
C
      IF(RMIN .LT. 0.1d0) THEN
C Check distance to other ions
         DO K = 1, NESP
            IF(PTYPE(K).NE.10 .AND. 
     1           DIST2(CPOINT,POTPT(1,K)) .LT. RRIONS) RETURN
         ENDDO
C All tests are passed; Add point to the pool, now
         CALL ADDPT(POTPT,NESP,CPOINT,PTYPE,ITYPE)
C
C DEBUGING: check angle between vectors CVEC and (CPOINT-CI)
C
         VX = CPOINT(1) - CI(1)
         VY = CPOINT(2) - CI(2)
         VZ = CPOINT(3) - CI(3)
         DEVIAT = ABS( CVEC(1)*VX + CVEC(2)*VY + CVEC(3)*VZ ) / 
     1      ( SQRT(CVEC(1)*CVEC(1)+CVEC(2)*CVEC(2)+CVEC(3)*CVEC(3)) * 
     2        SQRT(VX*VX+VY*VY+VZ*VZ) )
         IF(DEVIAT .LT. DEVMAX) THEN
            IPOINT = NESP
            IAT    = IATOM
            JAT    = JATOM
            DEVMAX = DEVIAT
         ENDIF
      ENDIF
C
      RETURN
      END

