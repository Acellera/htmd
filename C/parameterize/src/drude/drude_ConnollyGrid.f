CGRID: Victor Anisimov, 2004-2005 
C Perturbation ion and grid generation program placing ions and grid points on Connolly surfaces
C
      program ConnollyGrid
      INCLUDE 'SIZES'

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      parameter (nSurfMax = 128)
      integer NSURF
      integer nesp, natoms
      integer nat(NUMATM)
      integer i, j, k, l, l1, l2, m, n, iPoint, iAtom
      integer  nPoints(nSurfMax)
      integer  nDiffElems
      integer  DiffElems(100)

      double precision  POTPT(3,MESP), SPNT(3,MESP)
      INTEGER           PTYPE(MESP), TMP(MESP)
      double precision  COORD(3,NUMATM)
      double precision  SCALE(NUMATM), DENS(NUMATM), SDENS(NUMATM)
      double precision  VANDER(100)

      COMMON /DEBUG/ NBONDS,IPOINT,IAT,JAT,DEVMAX
      COMMON /THRESH/ RRIONS, RRATMS

      character*256 JOBNAM, CHARMM
      character*256 EXENAME, str, STRING
      character*80  GAUFILE(200), SURFSTR(100)
      logical       Connolly(nSurfMax), Bonds(nSurfMax)
      logical       Angles(nSurfMax), Tests(nSurfMax)
      integer       Exclude(NUMATM,nSurfMax), Include(NUMATM,nSurfMax)
      integer       Testarc(4,nSurfMax)
      INTEGER       NPTSURF(nSurfMax)
      LOGICAL       OPENFL

      double precision surf_param(nSurfMax), dens_param(nSurfMax)
      double precision RATMS(nSurfMax), RIONS(nSurfMax)
      integer  natread
      character*4   pLabel
      logical       GenGauFile
      CHARACTER*4   RESIDUE
      INTEGER       SUBSTR, RGLINE, STRLEN
      double precision extra(3)
      logical       first

      double precision READFLOAT
      integer          READINT

      integer          h

      DATA VANDER/1.20D0,                                     1.20D0,
     *       1.37D0,1.45D0,1.45D0,1.50D0,1.50D0,1.40D0,1.35D0,1.30D0,
     *       1.57D0,1.36D0,1.24D0,1.17D0,1.80D0,1.75D0,1.70D0,1.0000,
     *       16*0.0D0,2.3D0,65*0.0D0/

      character*2 ELDAT(56)
      DATA ELDAT / 'H ',                              'HE',
     *             'LI','BE','B ','C ','N ','O ','F ','NE',
     *             'NA','MG','AL','SI','P ','S ','CL','AR',
     *             'K ','CA',
     * 'SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN',
     *                       'GA','GE','AS','SE','BR','KR',  
     *             'RB','SR',
     * 'Y ','ZR','NB','MO','TC','RU','RH','PD','AG','CD', 
     *                       'IN','SN','SB','TE','I ','XE', 
     *             'CS','BA'/

      DEVMAX = 1.0d0
      IAT = 0
      JAT = 0
      PI = 3.14159265358d0

      NARGS = IARGC()
      CALL GETARG(1, JOBNAM)
      CALL GETARG(0, EXENAME)
      GenGauFile = NARGS.GT.1

      do j = 1, nSurfMax
         do i = 1, NumAtm
            Exclude(i,j) = 0
            Include(i,j) = 0
         enddo
      enddo

*Read configuration (cfg) file
      IF(OPENFL(JOBNAM,10,0)) GOTO 40
      read(10,'(A80)',ERR=40, END=40) STRING
      i = INDEX(STRING, "=")
      read(STRING(i+1:),*,err=40) NSURF
      do i = 1, NSURF
        read(10,'(A256)') str
        STRING = str
        h = index(STRING,"#")
        if(h.eq.0) h = STRLEN(STRING)
        CALL UPCASE(STRING)

        read(STRING,*) surf_param(i), dens_param(i)
C Read distance thresholds in Angsroms
        RIONS(I)    = READFLOAT(STRING, 3, 1.0d0)
        RATMS(I)    = READFLOAT(STRING, 4, 2.0d0)
C Read point types
        Angles(i)   = INDEX(STRING(1:h),'A') .GT. 0 
        Bonds(i)    = INDEX(STRING(1:h),'B') .GT. 0 
        Connolly(i) = INDEX(STRING(1:h),'C') .GT. 0 
        Tests(i)    = INDEX(STRING(1:h),'T') .GT. 0 
        if( Tests(i) ) then
C Read arc instructions; 1,2,3 atoms, and rotation Angle
           Testarc(1,i) = READINT(STRING, 6, 0)
           Testarc(2,i) = READINT(STRING, 7, 0)
           Testarc(3,i) = READINT(STRING, 8, 0)
           Testarc(4,i) = READINT(STRING, 9, 0)
        elseif( INDEX(STRING(i:h),'E') .GT. 0 ) then
C Read excluded atoms
           j = 2
           k = 6
           n = 0
    5      continue
              Exclude(j,i) = READINT(STRING, k, 0)
              if(Exclude(j,i) .gt. 0) then
                 j = j + 1
                 k = k + 1
                 n = n + 1
                 if(n.lt.NumAtm) goto 5
              else
C                total number of excluded atoms
                 if(n.gt.0) Exclude(1,i) = n
              endif
        elseif( INDEX(STRING(i:h),'I') .GT. 0 ) then
C Read included atoms; the exclusion list will be built later out of the inclusion list
           j = 2
           k = 6
           n = 0
    6      continue
              Include(j,i) = READINT(STRING, k, 0)
              if(Include(j,i) .gt. 0) then
                 j = j + 1
                 k = k + 1
                 n = n + 1
                 if(n.lt.NumAtm) goto 6
              else
C                total number of included atoms
                 if(n.gt.0) Include(1,i) = n
              endif
        endif
        SURFSTR(I)  = str
      enddo
      CLOSE(10)

*Open Gaussian file on input
      IF(OPENFL(JOBNAM,11,0)) GOTO 20

*Read Gaussian keywords
      ISTR = 0
      do While( RGLINE(11,GAUFILE,ISTR,100) .GT. 0 )
      enddo
 
*Read commentaries
      I = RGLINE(11,GAUFILE,ISTR,100)
C Read residue name from the first comment line
      CALL CPDATA(GAUFILE(ISTR),RESIDUE)
      IF(RESIDUE .EQ. ' ') THEN
         WRITE(*,*) ' NO CHARMM NAME IN COMMENT LINE OF GAUSSIAN FILE'
         STOP
      ENDIF
      do While( RGLINE(11,GAUFILE,ISTR,100) .GT. 0 )
      enddo

*Charge and multiplicity
      I = RGLINE(11,GAUFILE,ISTR,100)
      natoms = 0
*Memorize the starting line number for atomic coordinates
      ICOORD = ISTR + 1
*Atomic coordinates
      do While( RGLINE(11,GAUFILE,ISTR,100) .GT. 0 )
        STRING = GAUFILE(ISTR)
        natoms = natoms + 1
        nat(natoms) = NATREAD(STRING,80,ELDAT,56)
        if(nat(natoms) .eq. 0) then
           write(*,'(a/a)') 'Error: unknown atom label',GAUFILE(ISTR)
           stop
        endif
        read(STRING,*,err=10) (COORD(k,natoms),k=1,3) 
      end do
   10 continue
      write(*,*) natoms
*Check end of the Gaussian file
      IF(INDEX(GAUFILE(ISTR),".") .GT. 0) THEN
C        This is a coordinate line; add blank line after it
         ISTR = ISTR + 1
         GAUFILE(ISTR) = ' '
      ENDIF
      CLOSE(11)

*Convert Inclusion list into Exclusion one
      do i = 1, NSURF
         if(Include(1,i).gt.0) then
            j = 2
            do k = 1, natoms
               do l = 2, Include(1,i)+1
                  if(k.eq.Include(l,i)) goto 15
               enddo
               Exclude(j,i) = k
               j = j + 1
   15          continue
            enddo
            Exclude(1,i) = j - 2
         endif
c        write(*,*) i,":",(Exclude(k,i),k=2,Exclude(1,i)+1)
      enddo

*Start creating ions and grid points
      first = .false.
      NESP  = 0
      NIONS = 0
      do i = 1, NSURF
        do iAtom = 1, nAtoms
          SCALE(iAtom) = surf_param(i)
          DENS(iAtom)  = dens_param(i)
          SDENS(iAtom) = 10.0d0
        enddo
        RRIONS = RIONS(I) * RIONS(I)
        RRATMS = RATMS(I) * RATMS(I)

        KIONS = 0
        if((Bonds(i) .or. Angles(i)) .and..not.Tests(i)) then
          NSPNT = 0
          NESPOLD = NESP
          CALL SURFAC(SCALE,SDENS,natoms,COORD,nat,nspnt,SPNT,TMP,10,
     *                Exclude(1,i),.TRUE.)
          IF(Angles(i)) CALL PANGLE(SCALE,DENS,natoms,COORD,nat,NESP,
     *                          POTPT,PTYPE,2,SPNT,nspnt)
          IF(Bonds(i))  CALL PBONDS(SCALE,DENS,natoms,COORD,nat,NESP,
     *                          POTPT,PTYPE,1,SPNT,nspnt,Exclude(1,i))
C begin hack; add ions to the first surface 
          if(first) then
             extra(1)=1.52D0
             extra(2)=2.08D0
             extra(3)=0.00D0
             call addpt(POTPT,nesp,extra,PTYPE,1)
             extra(1)=-1.52D0
             extra(2)=2.08D0
             extra(3)=0.00D0
             call addpt(POTPT,nesp,extra,PTYPE,1)
             extra(1)=1.82D0
             extra(2)=2.20D0
             extra(3)=0.00D0
             call addpt(POTPT,nesp,extra,PTYPE,1)
             extra(1)=-1.82D0
             extra(2)=2.20D0
             extra(3)=0.00D0
             call addpt(POTPT,nesp,extra,PTYPE,1)
             first = .false.
          endif
C end hack
          KIONS = NESP - NESPOLD
          NIONS = NIONS + NESP - NESPOLD
        endif
        NESPOLD = NESP
        if(Connolly(i) .AND. Bonds(i)) then
          CALL SURFAC(SCALE,DENS,natoms,COORD,nat,NESP,POTPT,PTYPE,5,
     *                Exclude(1,i),.FALSE.)
          NIONS = NIONS + NESP - NESPOLD
        else if((Bonds(i).or.Angles(i)) .AND. Tests(i)) then
          CALL TSTARC(SCALE,DENS,natoms,COORD,nat,NESP,POTPT,PTYPE,5,
     *                Testarc(1,i))
          NIONS = NIONS + NESP - NESPOLD
        else if(Connolly(i) .AND. Tests(i)) then
          CALL TSTARC(SCALE,DENS,natoms,COORD,nat,NESP,POTPT,PTYPE,10,
     *                Testarc(1,i))
        else if(Connolly(i) .AND..NOT.Bonds(i) .AND..NOT.Tests(i)) then
          CALL SURFAC(SCALE,DENS,natoms,COORD,nat,NESP,POTPT,PTYPE,10,
     *                Exclude(1,i),.FALSE.)
        endif
        KESP = NESP - NESPOLD
        NPTSURF(I) = KIONS + KESP
        IF( Bonds(i) .or. Angles(i)) THEN
          write(*,'(A,I2,1X,A1,1X,3(A,I4,2X),A,2I3,A,F6.2)') 
     1         'Surface: ',i, CHAR(64+i),
     2         'pbonds:',   KIONS,
     3         'pconnolly:',KESP,
     4         'nBonds:',   NBONDS,
     5         'Bond:', IAT,JAT, ' Dev:',ACOS(DEVMAX)*180.0d0/PI
          DEVMAX = 1.0d0
          IAT = 0
          JAT = 0
        ELSE
           write(*,'(A,I2,1X,A1,1X,2(A,I4,2X))') 
     1         'Surface: ',i, CHAR(64+i),
     2         'pbonds:',   KIONS,
     3         'pconnolly:',KESP
        ENDIF
      ENDDO

*count number of ions
      NIONS = 0
      do i = 1, NESP
         IF(PTYPE(I).NE.10) NIONS = NIONS + 1
      enddo

*write MESP (pgrid.xyz) for QM calculations
      IF(NARGS.GT.1) THEN
        IF(OPENFL(JOBNAM,12,0)) GOTO 30
        iResidue = 1
        jResidue = 1
        DO iPoint=1, nesp
          if(PTYPE(iPoint) .eq. 10) then 
            write(12, '(3F10.5)') (POTPT(i,iPoint), i=1,3)
            iResidue = iResidue + 1
          endif
        ENDDO
        CLOSE(12)
      ENDIF

*CHARMM visualization
      IF(OPENFL(JOBNAM,20,0)) GOTO 30
      CALL WRTITLE(20,SURFSTR,NSURF)
      WRITE(20,'(I5)') NESP + nAtoms
      CALL WRCHARMM(20,'MOL ',GAUFILE(ICOORD),nAtoms, COORD)
      CHARMM(1:4) = ' '
      iStart = 1
      iEnd   = 0
      iResidue = 1
      DO i=1, NSURF
        IF(Bonds(i) .or. Angles(i)) THEN
          CHARMM(1:2) = 'CH'
        ELSE 
          CHARMM(1:2) = 'ES'
        ENDIF 
        CHARMM(3:3) = CHAR(64+i)
        iEnd = iEnd + NPTSURF(I)
        DO iPoint = iStart, iEnd
          iResidue = iResidue + 1
          IF(PTYPE(iPoint) .EQ. 1) THEN
             pLabel = "B   "
          ELSEIF(PTYPE(iPoint) .EQ. 5) THEN
             pLabel = "F   "
          ELSE
             pLabel = "H   "
          ENDIF
          M = nAtoms + iPoint
          CALL WCHSTR(20,CHARMM,pLabel,POTPT(1,iPoint),M,iResidue)
        ENDDO
        iStart = iEnd + 1
      ENDDO
      CLOSE(20)

* CHARMM gpos.N.crd files
      IF(NARGS.GT.1) THEN
        ICAL = 0
        M = nAtoms + 1
        IF(OPENFL(JOBNAM,21,ICAL)) GOTO 30
        DO iPoint=1, NESP
          IF(PTYPE(iPoint) .NE. 10) THEN 

* write cumulative file
            do iAtom=1,natoms
               write(21,'(3F10.5)') (COORD(k,iAtom),k=1,3)
            enddo
            write(21,'(3F10.5)') (POTPT(k,iPoint),k=1,3)

            ICAL = ICAL + 1
            IF(ICAL.EQ.1) THEN
* do not write individual .crd files but one
              IF(OPENFL(JOBNAM,22,ICAL)) GOTO 30
              CALL WRTITLE(22,SURFSTR,NSURF)
              WRITE(22,'(I5)') M
              CALL WRCHARMM(22,RESIDUE,GAUFILE(ICOORD),nAtoms,COORD) 
              CALL WCHSTR(22,'CAL ','CAL ',POTPT(1,iPoint),M,2)
              CLOSE(22)
            ENDIF
          ENDIF
        ENDDO
        CLOSE(21)
      ENDIF

* CHARMM pgrid.crd file
*      IF(NARGS.GT.1) THEN
*        IF(OPENFL(JOBNAM,25,0)) GOTO 30
*        CALL WRTITLE(25,SURFSTR,NSURF)
*        WRITE(25,'(I5)') nAtoms+NESP-NIONS
*        CALL WRCHARMM(25,RESIDUE,GAUFILE(ICOORD),nAtoms,COORD) 
*        M = nAtoms
*        iResidue = 1
*        DO iPoint=1, NESP
*          IF(PTYPE(iPoint) .EQ. 10) THEN 
*            M = M + 1
*            iResidue = iResidue + 1
*            CALL WCHSTR(25,'DUM ','DUM ',POTPT(1,iPoint),M,iResidue)
*          ENDIF
*        ENDDO
*        CLOSE(25)
*      ENDIF
      
*generate Gaussian files
      IF(GenGauFile) THEN
         ION = 0
         IPTR = INDEX(JOBNAM,'.')
            IF(OPENFL(JOBNAM,29,ION)) GOTO 30
            CALL WRITEGAU(29,GAUFILE,ISTR,.TRUE.)
            write(29,"(A)") '@pgrid.xyz /N'
            write(29,*) ' '
C
* individual .gjf files
*            IF(OPENFL(JOBNAM,30,ION)) GOTO 30
*            CALL WRITEGAU(30,GAUFILE,ISTR,.TRUE.)
*            write(30,"(A)") '@'//JOBNAM(:IPTR)//'pgrid.xyz /N'
*            write(30,*) ' '
*            close(30)
         do i = 1, NESP
            IF(PTYPE(I).NE.10) THEN
               write(29,"(A9)") '--Link1--'
               CALL WRITEGAU(29,GAUFILE,ISTR,.FALSE.)
               WRITE(29,'(3F10.5,3x,"0.5")') (POTPT(K,I),K=1,3)
               write(29,*) ' '
               write(29,"(A)") '@pgrid.xyz /N'
               write(29,*) ' '
C
               ION = ION + 1
* individual .gjf files
*               IF(OPENFL(JOBNAM,30,ION)) GOTO 30
*               CALL WRITEGAU(30,GAUFILE,ISTR,.FALSE.)
*               WRITE(30,'(3F10.5,3x,"0.5")') (POTPT(K,I),K=1,3)
*               WRITE(30,*) ''
*               WRITE(30,"(A)") '@'//JOBNAM(:iptr)//'pgrid.xyz /N'
*               WRITE(30,*) ''
*               CLOSE(30)
            ENDIF
         ENDDO
         WRITE(*,*) ION, " Gaussian jobs were created"
         CLOSE(29)
      ENDIF

      write(iw,"(2(i4,a))") NESP-NIONS,' ESP points ',NIONS,' Ions'
      stop
   20 write(iw,'(a)') " Cannot open (or read) input file."
      stop
   30 write(iw,'(a)') " Cannot write to CHARMM file."
      stop
   40 write(iw,'(a)') " Cannot open file with surfaces parameters."
      stop
      end 


* convert all characters to upper case 
      SUBROUTINE UPCASE(STRING)
      CHARACTER*(*) STRING
      CHARACTER*1   CH
      INTEGER       NSIZE, I, ICH
      INTEGER       STRLEN
      NSIZE = STRLEN(STRING)
      DO I=1,NSIZE
         ICH = ICHAR(STRING(I:I))
         IF(ICH.GE.97 .AND. ICH.LE.122) STRING(I:I)=CHAR(ICH-32)
      ENDDO
      RETURN
      END


* read atom number from string
      INTEGER FUNCTION NATREAD(STRING,NSIZE,TABLE,NEL)
      CHARACTER*(*) STRING
      CHARACTER*2   TABLE(NEL)
      INTEGER       I, IPTR, ICH
      INTEGER       SUBSTR
C
      NATREAD=0
      CALL UPCASE(STRING)
C     Find first non-space character
      IPTR = SUBSTR(STRING,1)
      IF(IPTR .EQ. 0) RETURN
C     Find space (32) or number (48-57); leave two characters of atom; 
C                                                      clean the rest
      DO I=IPTR+1,NSIZE
         ICH = ICHAR(STRING(I:I))
         IF(ICH .EQ. 32) THEN
            GOTO 20
         ELSEIF((ICH.GE.48 .AND. ICH.LE.57) .OR. I.GT.IPTR+1) THEN
            STRING(I:I) = ' '
         ENDIF
      ENDDO
   20 CONTINUE
C     Find atom symbol
      DO I=1,NEL
         IF(STRING(IPTR:IPTR+1) .EQ. TABLE(I)) THEN
            NATREAD = I
            STRING(IPTR:IPTR+1) = ' '
            RETURN
         ENDIF
      ENDDO
      RETURN
      END


      LOGICAL FUNCTION OPENFL(JOBNAM,ITYPE,NUMBER)
      CHARACTER*(*) JOBNAM
      CHARACTER*256 EXT
      INTEGER       I, ITYPE, NUMBER
      CHARACTER*3   STR
      OPENFL = .TRUE.
      I = INDEX(JOBNAM,".")
      write (*,*) 'openfl'
      write (*,*) ITYPE, NUMBER
      IF(ITYPE.EQ.10) THEN
C CONFIGURATION FILE
         OPEN(UNIT=ITYPE, FILE='ConnollyGrid.cfg', 
     1        STATUS='OLD', ERR=10)
         OPENFL = .FALSE.
         RETURN
      ELSEIF(ITYPE.EQ.11) THEN
C GAUSSIAN TEMPLATE FILE
C         EXT = JOBNAM(:I)//'gjf'
         OPEN(UNIT=ITYPE, FILE=JOBNAM, STATUS='OLD', ERR=10)
         OPENFL = .FALSE.
         RETURN
      ELSEIF(ITYPE.EQ.12) THEN
         EXT = 'pgrid.xyz'
C         EXT = JOBNAM(:I)//'pgrid.xyz' 
      ELSEIF(ITYPE.EQ.20) THEN
C VISUALIZATION OF IONS AND GRID POINTS
         EXT = JOBNAM(:I)//'crd'
      ELSEIF(ITYPE.EQ.21) THEN
C QPOS CHARMM CUMULATIVE FILE
         CALL NUMSTR(STR,NUMBER,IN)
C         EXT = JOBNAM(:I)//'qpos.all'
         EXT = 'qpos.all'
      ELSEIF(ITYPE.EQ.22) THEN
C QPOS CHARMM FILE
         CALL NUMSTR(STR,NUMBER,IN)
         EXT = JOBNAM(:I)//'qpos.'//STR(1:IN)//'.crd'
      ELSEIF(ITYPE.EQ.25) THEN
C GRIP POINTS IN CHARMM FORMAT
         EXT = JOBNAM(:I)//'pgrid.crd'
      ELSEIF(ITYPE.EQ.29)THEN
C OUTPUT GJF FILE
         CALL NUMSTR(STR,NUMBER,IN)
         EXT = JOBNAM(:I)//'gjf.all' 
      ELSEIF(ITYPE.EQ.30)THEN
C OUTPUT GJF FILE
         CALL NUMSTR(STR,NUMBER,IN)
         EXT = JOBNAM(:I)//STR(1:IN)//'.gjf' 
      ELSE
C RETURN ERROR FOR UNKNOWN CASE
         RETURN
      ENDIF
      OPEN(UNIT=ITYPE, FILE=EXT, STATUS='UNKNOWN', ERR=10)
      CLOSE(ITYPE,STATUS='DELETE')
      OPEN(UNIT=ITYPE, FILE=EXT, STATUS='NEW', ERR=10)
      OPENFL = .FALSE.
   10 RETURN
      END


C Convert number to string with leading zeros
      SUBROUTINE NUMSTR(STR,NUMBER,IN)
      CHARACTER*3 STR
      INTEGER     NUMBER
C
      IF(NUMBER.GT.999) THEN
         WRITE(*,*) ' NUMBER OF FILES TO BE GENERATED ',NUMBER
         WRITE(*,*) ' EXCEEDS 999'
         STOP
      ENDIF
      IF(NUMBER .LE. 999 .AND. NUMBER .GT. 99) THEN
         WRITE(STR,'(I3)') NUMBER
         IN = 3
      ELSEIF (NUMBER .LE. 99 .AND. NUMBER .GT. 9) THEN
         WRITE(STR(1:2),'(I2)') NUMBER
         IN = 2
      ELSEIF (NUMBER .LE. 9 .AND. NUMBER .GE. 0) THEN
         WRITE(STR(1:1),'(I1)') NUMBER
         IN = 1
      ENDIF
C
C      IN = 3
C      STR = '000'
C      IF(NUMBER .LE. 999 .AND. NUMBER .GT. 99) THEN
C         WRITE(STR,'(I3)') NUMBER
C      ELSEIF (NUMBER .LE. 99 .AND. NUMBER .GT. 9) THEN
C         WRITE(STR(2:3),'(I2)') NUMBER
C      ELSEIF (NUMBER .LE. 9 .AND. NUMBER .GE. 0) THEN
C         WRITE(STR(3:3),'(I1)') NUMBER
C      ELSE
C         WRITE(*,*) ' NUMBER OF FILES TO BE GENERATED ',NUMBER
C         WRITE(*,*) ' EXCEEDS 999'
C         STOP
C      ENDIF
C
      RETURN
      END


C Find Nth substring and return its index
      INTEGER FUNCTION SUBSTR(STRING,N)
      CHARACTER*(*) STRING
      INTEGER       NSIZE, N
      INTEGER       STRLEN
C
      INTEGER       I, ISUBS ,ISHIFT
C
      NSIZE = STRLEN(STRING)
      SUBSTR = 0
      ISUBS  = 0
C Loop over substrings
      I = 1
   10 IF(STRING(I:I) .NE. ' ') THEN
C Count substrings
         ISUBS = ISUBS + 1
         IF(ISUBS .EQ. N) THEN
            SUBSTR = I
            RETURN
         ENDIF
C Find first blank character
   20    I = I + 1
         IF(I.LE.NSIZE .AND. STRING(I:I).NE.' ') GOTO 20
      ENDIF
C Skip blank characters
      I = I + 1
      IF(I .LE. NSIZE) GOTO 10
      RETURN
      END


C Write coordinates in CHARMM format
      SUBROUTINE WRCHARMM(IU,RESIDUE,GAUFILE,nAtoms,COORD)
      INTEGER          IU, nAtoms
      DOUBLE PRECISION COORD(3,nAtoms)
      CHARACTER*80     GAUFILE(nAtoms)
      CHARACTER*4      RESIDUE, ALABEL
C
      CHARACTER*80 STRING
      INTEGER      iAtom, LSTART, I
      INTEGER      SUBSTR, IPTR
C
      CALL UPCASE(RESIDUE)
      DO iAtom = 1, nAtoms
* extract atom label
        STRING = GAUFILE(iAtom)
* CHARMM name is expected after exclamation "!" mark
        IPTR = INDEX(STRING,"!")
        IF(IPTR.GT.0) THEN
* check space character between "!" and CHARMM name; read 6th argument
           LSTART = SUBSTR(STRING,6)
* if failed, set the pointer to next character after "!" 
           IF(LSTART .LE. 0) LSTART = IPTR + 1
* finally, check that CHARMM name is found 
           IF(STRING(LSTART:LSTART).EQ.' ') LSTART = 0
        ELSE
* read first substring
           LSTART = SUBSTR(STRING,1)
        ENDIF
        IF(LSTART.LE.0) THEN
          WRITE(*,*) ' WRONG CHARMM ATOM LABELS'
          STOP
        ENDIF
        ALABEL = ' '
        I = 1
* make left alignement
   10   CONTINUE
        IF(STRING(LSTART:LSTART).NE.' ' .AND. I.LE.4) THEN
           ALABEL(I:I) = STRING(LSTART:LSTART)
           I = I + 1
           LSTART = LSTART + 1
           GOTO 10
        ENDIF
        CALL WCHSTR(IU,RESIDUE,ALABEL,COORD(1,iAtom),iAtom,1)
        ICOORD = ICOORD + 1
      ENDDO
      RETURN
      END


C Write single CHARMM string
      SUBROUTINE WCHSTR(IU,RESIDUE,NAME,COORD,M,iResidue)
      INTEGER          IU, M, iResidue
      CHARACTER*4      RESIDUE
      CHARACTER*4      NAME
      DOUBLE PRECISION COORD(3)
C
      INTEGER          I
C
      WRITE(IU,'(2i5,1x,a4,1x,a4,3F10.5,1x,a4,1x,i4,f10.5)') 
     * M,iResidue,RESIDUE,NAME,(COORD(I),I=1,3),
     * RESIDUE, iResidue, 0.
      RETURN
      END


C Write CHARMM title
      SUBROUTINE WRTITLE(IU,SURFSTR,NSURF)
      INTEGER       IU,NSURF
      CHARACTER*80  SURFSTR(NSURF)
      INTEGER       I
      CHARACTER*80  FNAME
      WRITE(IU,"('* NSURF=',I2)") NSURF
      DO I = 1, NSURF
         WRITE(IU,"('* ',A)") SURFSTR(I)
      ENDDO
* Write file name
      INQUIRE(UNIT=IU,NAME=FNAME)
      I = INDEX(FNAME,' ') - 1
      WRITE(IU,"('* ',A)") FNAME(:I)
      WRITE(IU,"('* ')")
      RETURN
      END


C Write Gaussian Input file
      SUBROUTINE WRITEGAU(IU,GAUFILE,ISTR,FIRST)
      INTEGER   IU, ISTR
      CHARACTER*80 GAUFILE(200), STRING
      LOGICAL   FIRST
      INTEGER   I,L
      INTEGER   STRLEN
      LOGICAL   BEGIN
C
      BEGIN = .TRUE.
      WRITE(IU,"(A)") "%chk=esp.chk"
C      WRITE(IU,"(A)") "%mem=1000MB"
C      WRITE(IU,"(A)") "%nproc=4"
      DO I = 1, ISTR
         STRING = GAUFILE(I)
         L = INDEX(STRING,'!') - 1
         IF(L.LE.0) L = STRLEN(STRING)
         WRITE(IU,"(A)") STRING(:L)
         IF( BEGIN .AND. INDEX(STRING,'#') .GT. 0 ) THEN
C   Add the extra string only once
            BEGIN = .FALSE.
            IF(FIRST) THEN
               WRITE(IU,"(A)") "pop=nbo"
            ELSE
               WRITE(IU,"(A)") "charge guess=read"
            ENDIF
         ENDIF
      ENDDO
C
      RETURN
      END


C Find string length
      INTEGER FUNCTION STRLEN(STRING)
      CHARACTER*(*) STRING
      INTEGER   I
      STRLEN = 0
      DO I = LEN(STRING), 1, -1
         IF(STRING(I:I) .NE. ' ') THEN
            STRLEN = I
            RETURN
         ENDIF
      ENDDO
      RETURN
      END


C Read float variable
      DOUBLE PRECISION FUNCTION READFLOAT(STRING,IARG,DEF)
      CHARACTER*(*)     STRING
      INTEGER           IARG
      DOUBLE PRECISION  DEF
      LOGICAL           ISFLOAT
C
      INTEGER           I, SUBSTR
C
      READFLOAT = DEF
      I = SUBSTR(STRING,IARG)
      IF(I.GT.0 .AND. ISFLOAT(STRING(I:I)) ) THEN
         READ(STRING(I:),*,ERR=10, END=10) READFLOAT
      ENDIF
   10 CONTINUE
C
      RETURN
      END


C Read integer variable
      INTEGER FUNCTION READINT(STRING,IARG,IDEF)
      CHARACTER*(*)     STRING
      INTEGER           IARG
      INTEGER           IDEF
      LOGICAL           ISINT
C
      INTEGER           I, SUBSTR
C
      READINT = IDEF
      I = SUBSTR(STRING,IARG)

      IF( I.GT.0 ) THEN
	   IF( ISINT(STRING(I:I)) ) THEN
            READ(STRING(I:),*,ERR=10, END=10) READINT
  	   ENDIF
      ENDIF
   10 CONTINUE
C
      RETURN
      END


C Check float type character
      LOGICAL FUNCTION ISFLOAT(CH)
      CHARACTER*1  CH
      IF(CH.EQ.'.' .OR. CH.EQ.'-' .OR. CH.EQ.'+' .OR.
     +          (ICHAR(CH).GE.48 .AND. ICHAR(CH).LE.57) ) THEN
         ISFLOAT = .TRUE.
      ELSE
         ISFLOAT = .FALSE.
      ENDIF
      RETURN
      END


C Check integer type character
      LOGICAL FUNCTION ISINT(CH)
      CHARACTER*1  CH
      IF(ICHAR(CH).GE.48 .AND. ICHAR(CH).LE.57) THEN
         ISINT = .TRUE.
      ELSE
         ISINT = .FALSE.
      ENDIF
      RETURN
      END


      INTEGER FUNCTION RGLINE(IU,GAUFILE,ISTR,LIMIT) 
      INTEGER       IU, ISTR, LIMIT
      CHARACTER*80  GAUFILE(200)
      CHARACTER*80  STRING
C
      READ(IU,'(A80)',ERR=10, END=10) STRING
      ISTR = ISTR + 1
      IF(ISTR.GT.LIMIT) THEN
         WRITE(*,*) 'Too many atoms in the Gaussian file'
         STOP
      ENDIF
      GAUFILE(ISTR) = STRING
      IF (STRING .EQ. ' ') THEN
         RGLINE = 0
      ELSE
         RGLINE = ISTR
      ENDIF
      RETURN
   10 CONTINUE
      RGLINE = -1
C
      RETURN
      END


C Copy "data" from first string to second one ignoring leading blanks 
C                                      up to the receiver string size
      SUBROUTINE CPDATA(STRFROM,STRTO)
      CHARACTER*(*) STRFROM, STRTO
      INTEGER       LFROM,   LTO, ISTART, I
      INTEGER       STRLEN, SUBSTR
      LOGICAL       BLANK
C
      ISTART = SUBSTR(STRFROM,1)
      IF(ISTART .GT. 0) THEN
         ISTART = ISTART - 1
         LFROM = STRLEN(STRFROM) - ISTART
         LTO   = STRLEN(STRTO)
         IF(LFROM .GT. LTO) THEN
* Copy data untill first blank is found
            
            BLANK = .FALSE.
            DO I = 1, LTO
               IF(BLANK) THEN
                  STRTO(I:I) = ' '
               ELSE
                  STRTO(I:I) = STRFROM(ISTART+I:ISTART+I)
                  IF(STRTO(I:I).EQ.' ') BLANK = .TRUE.
               ENDIF
            ENDDO
         ELSE
* Copy data and add blanks
            DO I = 1, LFROM
               STRTO(I:I) = STRFROM(ISTART+I:ISTART+I)
            ENDDO
            DO I = LFROM+1, LTO
               STRTO(I:I) = ' '
            ENDDO
         ENDIF
         IF(STRTO(LTO:LTO).EQ.'_' .OR. STRTO(LTO:LTO).EQ.'-') 
     $       STRTO(LTO:LTO) = ' '
      ELSE
         STRTO = ' '
      ENDIF
C
      RETURN
      END

