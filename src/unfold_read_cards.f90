!
! Copyright (C) 2013 Pietro Bonfa' and
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
MODULE unfold_read_cards_module
   !---------------------------------------------------------------------------
   !
   ! ...  This module handles the reading of cards from standard input
   ! ...  Written by Carlo Cavazzoni and modified for "path" implementation
   ! ...  by Carlo Sbraccia
   !
   USE kinds,     ONLY : DP
   USE parser,    ONLY : parse_unit,field_count, read_line
   USE io_global, ONLY : meta_ionode
   !
   USE unfold_input_parameters
   !
   IMPLICIT NONE
   !
   SAVE
   !
   PRIVATE
   !
   PUBLIC :: unfold_read_cards
   !
   ! ... end of module-scope declarations
   !
   !  ----------------------------------------------
   !
CONTAINS
   !
   ! ... Read CARDS ....
   !
   ! ... subroutines
   !
   !----------------------------------------------------------------------
   !
   !----------------------------------------------------------------------
   SUBROUTINE unfold_read_cards(unit)
      !----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: unit
      !
      CHARACTER(len=256)         :: input_line
      CHARACTER(len=80)          :: card
      CHARACTER(len=1), EXTERNAL :: capital
      LOGICAL                    :: tend
      INTEGER                    :: i
      !
      !
      parse_unit = unit
      !
100   CALL read_line( input_line, end_of_file=tend )
      !
      IF( tend ) GOTO 120
      IF( input_line == ' ' .or. input_line(1:1) == '#' ) GOTO 100
      !
      READ (input_line, *) card
      !
      DO i = 1, len_trim( input_line )
         input_line( i : i ) = capital( input_line( i : i ) )
      ENDDO
      !
      IF( trim(card) =='TRMAT' ) THEN
         !
         CALL card_trmat( input_line )
         !
      ELSEIF( trim(card) =='UNKPTS' ) THEN
         !
         CALL card_unkpts( input_line )
         !
      ELSEIF( trim(card) =='CELL_PARAMETERS' ) THEN
         !
         CALL card_cell_parameters( input_line )
         !
      ELSE
         !
         IF ( meta_ionode ) &
            WRITE( 0,'(A)') 'Warning: card '//trim(input_line)//' ignored'
         !
      ENDIF
      !
      ! ... END OF LOOP ... !
      !
      GOTO 100
      !
120      CONTINUE
      !
      RETURN
      !
   END SUBROUTINE unfold_read_cards
   !
   !------------------------------------------------------------------------
   !    BEGIN manual
   !----------------------------------------------------------------------
   !
   ! TRMAT
   !
   !   specify transformation matrix from base cell to supercell
   !
   ! Syntax:
   !
   !    TRMAT
   !      HT(1,1) HT(1,2) HT(1,3)
   !      HT(2,1) HT(2,2) HT(2,3)
   !      HT(3,1) HT(3,2) HT(3,3)
   !
   ! Example:
   !
   ! TRMAT
   !     2.0    0.0    0.0
   !     0.0    2.0    0.0
   !     0.0    0.0    3.0
   !
   ! Where:
   !
   !      HT(i,j) (real)  supercell dimensions ( in units of base cell ),
   !
   !
   !
   !----------------------------------------------------------------------
   !    END manual
   !------------------------------------------------------------------------
   !
   SUBROUTINE card_trmat( input_line )
      !
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: input_line
      INTEGER            :: i, j
      !
      !
      IF ( tunfold ) THEN
         CALL errore( ' card_trmat', ' two occurrences', 2 )
      ENDIF
      !
      DO i = 1, 3
         CALL read_line( input_line )
         READ(input_line,*) ( trmat( i, j ), j = 1, 3 )
      ENDDO
      !
      tunfold  = .true.
      !
      RETURN
      !
   END SUBROUTINE card_trmat
   !
   !
   !------------------------------------------------------------------------
   !    BEGIN manual
   !----------------------------------------------------------------------
   !
   ! UNKPTS
   !
   !   use the specified set of k points for unfolded band structure
   !   k points are referred to the SUPERCELL BZ
   !
   ! Syntax:
   !
   !   UNKPTS (mesh_option)
   !     n
   !     xk(1,1) xk(2,1) xk(3,1) wk(1)
   !     ...     ...     ...     ...
   !     xk(1,n) xk(2,n) xk(3,n) wk(n)
   !
   ! Example:
   !
   ! UNKPTS
   !   10
   !    0.1250000  0.1250000  0.1250000   1.00
   !    0.1250000  0.1250000  0.3750000   3.00
   !    0.1250000  0.1250000  0.6250000   3.00
   !    0.1250000  0.1250000  0.8750000   3.00
   !    0.1250000  0.3750000  0.3750000   3.00
   !    0.1250000  0.3750000  0.6250000   6.00
   !    0.1250000  0.3750000  0.8750000   6.00
   !    0.1250000  0.6250000  0.6250000   3.00
   !    0.3750000  0.3750000  0.3750000   1.00
   !    0.3750000  0.3750000  0.6250000   3.00
   !
   ! Where:
   !
   !   mesh_option == crystal    k points mesh is given in stdin in scaled
   !                             units
   !   mesh_option == tpiba      k points mesh is given in stdin in units
   !                             of ( 2 PI / alat )
   !   mesh_option == gamma      only gamma point is used
   !   mesh_option == tpiba_b    as tpiba but the weights gives the
   !                             number of points between this point
   !                             and the next
   !   mesh_option == crystal_b  as crystal but the weights gives the
   !                             number of points between this point and
   !                             the next
   !
   !   n       ( integer )  number of k points
   !   xk(:,i) ( real )     coordinates of i-th k point
   !   wk(i)   ( real )     weights of i-th k point
   !
   !----------------------------------------------------------------------
   !    END manual
   !------------------------------------------------------------------------
   !
   SUBROUTINE card_unkpts( input_line )
      !
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: input_line
      INTEGER            :: i, j
      INTEGER            :: unnkaux
      INTEGER, ALLOCATABLE :: unwkaux(:)
      REAL(DP), ALLOCATABLE :: unxkaux(:,:)
      REAL(DP) :: undelta, unwk0
      LOGICAL, EXTERNAL  :: matches
      LOGICAL            :: tend,terr
      LOGICAL            :: kband = .false.
      !
      !

      IF ( tunkpts_inp ) THEN
         CALL errore( ' card_unkpts ', ' two occurrences', 2 )
      ENDIF
      !
      IF ( matches( "CRYSTAL", input_line ) ) THEN
         !  input k-points are in crystal (reciprocal lattice) axis
         unk_points = 'crystal'
         IF ( matches( "_B", input_line ) ) kband=.true.
      ELSEIF ( matches( "TPIBA", input_line ) ) THEN
         !  input k-points are in 2pi/a units
         unk_points = 'tpiba'
         IF ( matches( "_B", input_line ) ) kband=.true.
      ELSEIF ( matches( "GAMMA", input_line ) ) THEN
         !  Only Gamma (k=0) is used
         unk_points = 'gamma'
      ELSE
         !  by default, input k-points are in 2pi/a units
         unk_points = 'tpiba'
      ENDIF
      !

      IF ( ( unk_points == 'tpiba' ) .or. ( unk_points == 'crystal' ) ) THEN
         !
         ! ... input k-points are in 2pi/a units
         !
         CALL read_line( input_line, end_of_file = tend, error = terr )
         IF (tend) GOTO 101
         IF (terr) GOTO 201
         READ(input_line, *, END=101, ERR=201) unnkstot
         !
         IF (.NOT. kband) THEN
            ALLOCATE ( unxk(3, unnkstot), unwk(unnkstot) )
            DO i = 1, unnkstot
               CALL read_line( input_line, end_of_file = tend, error = terr )
               IF (tend) GOTO 101
               IF (terr) GOTO 201
               READ(input_line,*, END=101, ERR=201) unxk(1,i),unxk(2,i),unxk(3,i),unwk(i)
            ENDDO
         ELSE
            unnkaux=unnkstot
            ALLOCATE(unxkaux(3,unnkstot), unwkaux(unnkstot))
            DO i = 1, unnkstot
               CALL read_line( input_line, end_of_file = tend, error = terr )
               IF (tend) GOTO 101
               IF (terr) GOTO 201
               READ(input_line,*, END=101, ERR=201) unxkaux(1,i), unxkaux(2,i), &
                                                  unxkaux(3,i), unwk0
               unwkaux(i) = NINT ( unwk0 ) ! beware: wkaux is integer
            ENDDO
            ! Count k-points first
            unnkstot=0
            DO i=1,unnkaux-1
               IF ( unwkaux(i) > 0 ) THEN
                  unnkstot=unnkstot+unwkaux(i)
               ELSEIF ( unwkaux(i) == 0 ) THEN
                  unnkstot=unnkstot+1
               ELSE
                  CALL errore ('card_unkpts', 'wrong number of points',i)
               ENDIF
            ENDDO
            unnkstot=unnkstot+1
            ALLOCATE ( unxk(3,unnkstot), unwk(unnkstot) )
            ! Now fill the points
            unnkstot=0
            DO i=1,unnkaux-1
               IF (unwkaux(i)>0) THEN
                  undelta=1.0_DP/unwkaux(i)
                  DO j=0,unwkaux(i)-1
                     unnkstot=unnkstot+1
                     unxk(:,unnkstot)=unxkaux(:,i)+undelta*j*(unxkaux(:,i+1)-unxkaux(:,i))
                     unwk(unnkstot)=1.0_DP
                  ENDDO
               ELSEIF (unwkaux(i)==0) THEN
                  unnkstot=unnkstot+1
                  unxk(:,unnkstot)=unxkaux(:,i)
                  unwk(unnkstot)=1.0_DP
               ELSE
                  CALL errore ('card_unkpts', 'wrong number of points',i)
               ENDIF
            ENDDO
            unnkstot=unnkstot+1
            unxk(:,unnkstot)=unxkaux(:,unnkaux)
            unwk(unnkstot)=1.0_DP
            DEALLOCATE(unxkaux)
            DEALLOCATE(unwkaux)
         ENDIF
         !
      ELSEIF ( unk_points == 'gamma' ) THEN
         !
         unnkstot = 1
         ALLOCATE ( unxk(3,1), unwk(1) )
         unxk(:,1) = 0.0_DP
         unwk(1) = 1.0_DP
         !
      ENDIF
      !
      tunkpts_inp = .true.
      !
      RETURN
101     CALL errore ('card_unkpts', ' end of file while reading ' &
            & // trim(unk_points) // ' k points', 1)
201     CALL errore ('card_unkpts', ' error while reading ' &
            & // trim(unk_points) // ' k points', 1)
      !
   END SUBROUTINE card_unkpts
   !
   !
   !------------------------------------------------------------------------
   !    BEGIN manual
   !----------------------------------------------------------------------
   !
   ! CELL_PARAMETERS
   !
   !   use the specified cell dimensions
   !
   ! Syntax:
   !
   !    CELL_PARAMETERS (cell_option)
   !      HT(1,1) HT(1,2) HT(1,3)
   !      HT(2,1) HT(2,2) HT(2,3)
   !      HT(3,1) HT(3,2) HT(3,3)
   !
   !   cell_option == alat      lattice vectors in units of alat
   !   cell_option == bohr      lattice vectors in Bohr
   !   cell_option == angstrom  lattice vectors in Angstrom
   !
   ! Example:
   !
   ! CELL_PARAMETERS
   !    24.50644311    0.00004215   -0.14717844
   !    -0.00211522    8.12850030    1.70624903
   !     0.16447787    0.74511792   23.07395418
   !
   ! Where:
   !
   !      HT(i,j) (real)  cell dimensions ( in a.u. ),
   !                      note the relation with lattice vectors:
   !                      HT(1,:) = A1, HT(2,:) = A2, HT(3,:) = A3
   !
   !----------------------------------------------------------------------
   !    END manual
   !------------------------------------------------------------------------
   !
   SUBROUTINE card_cell_parameters( input_line )
      !
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: input_line
      INTEGER            :: i, j
      LOGICAL, EXTERNAL  :: matches
      !
      !
      IF ( tcell ) THEN
         CALL errore( ' card_cell_parameters ', ' two occurrences', 2 )
      ENDIF
      !
      IF ( matches( "BOHR", input_line ) ) THEN
         cell_units = 'bohr'
      ELSEIF ( matches( "ANGSTROM", input_line ) ) THEN
         cell_units = 'angstrom'
      ELSEIF ( matches( "ALAT", input_line ) ) THEN
         cell_units = 'alat'
      ELSE
         cell_units = 'none'
         CALL infomsg( 'read_cards ', &
            & 'DEPRECATED: no units specified in CELL_PARAMETERS card' )
         ! Cell parameters are set in cell_base_init
      ENDIF
      !
      DO i = 1, 3
         CALL read_line( input_line )
         READ(input_line,*) ( rd_ht( i, j ), j = 1, 3 )
      ENDDO
      !
      trd_ht = .true.
      tcell  = .true.
      !
      RETURN
      !
   END SUBROUTINE card_cell_parameters
END MODULE unfold_read_cards_module
