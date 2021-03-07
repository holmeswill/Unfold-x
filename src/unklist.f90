!
! Copyright (C) 2013 Pietro Bonfa' and
! This file is distributed under the terms of the
! GNU General Public License.
! See http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM unklist
  !-----------------------------------------------------------------------
  !
  !
  USE io_global, ONLY: stdout
  USE mp_global,  ONLY : mp_startup, mp_global_end
  USE kinds,     ONLY: DP
  USE matrix_inversion, ONLY : invmat
  USE read_cards_module, ONLY : read_cards
  USE unfold_input_parameters, ONLY : trmat, unxk, unnkstot, unk_points, &
                                        trd_ht, rd_ht, cell_units
  USE unfold_data, ONLY : get_klist
  USE unfold_read_cards_module, ONLY : unfold_read_cards
  USE cell_base,     ONLY : at, bg, alat, cell_base_init

  IMPLICIT NONE
  ! variables for cell structure
  INTEGER :: ibrav = 14
      ! index of the the Bravais lattice
      ! Note: in variable cell CP molecular dynamics, usually one does
      !       not want to put constraints on the cell symmetries, thus
      !       ibrav = 14 is used

  REAL(DP) :: celldm(6) = 0.0_DP
      ! dimensions of the cell (lattice parameters and angles)
  REAL(DP) :: a = 0.0_DP
  REAL(DP) :: c = 0.0_DP
  REAL(DP) :: b = 0.0_DP
  REAL(DP) :: cosab = 0.0_DP
  REAL(DP) :: cosac = 0.0_DP
  REAL(DP) :: cosbc = 0.0_DP
  !
  ! Variables for new k points
  INTEGER  :: nrpkp
  REAL(DP), ALLOCATABLE :: rpkp (:,:) !reduced kpoints
  REAL(DP), ALLOCATABLE :: G (:,:) !corresponding G vectors
  !
  ! and new direct and reciprocal space
  !
  REAL(DP) :: nat(3,3) = RESHAPE( (/ 0.0_DP /), (/ 3, 3 /), (/ 0.0_DP /) )
  REAL(DP) :: nbg(3,3) = RESHAPE( (/ 0.0_DP /), (/ 3, 3 /), (/ 0.0_DP /) )
  !  inverse of trmat
  REAL(DP) :: itrmat(3,3)
  !  determinant of the trmat matrix
  REAL(DP) :: det
  !
  ! single kpoint shifted to the small BZ
  REAL(DP) :: pkp(3)
  !
  ! counters
  INTEGER :: ipkp, i, j
  INTEGER :: ipol, apol, ios
  !
  ! input variables
  LOGICAL :: uniq = .true.
  LOGICAL :: symreduce = .false. !not implemented
  LOGICAL :: printg = .false.
  LOGICAL :: cart = .true.       !not implemented
  !
  !
  NAMELIST / input_unklist / printg, cart, uniq, ibrav, celldm, a, b, c, cosab, cosac, cosbc

#ifdef __MPI
  CALL mp_startup ( )
#endif

  CALL input_from_file ( )
  ! initialise data from input
  ios = 0
  !
  !     reading the namelist inputpp
  !
  READ (5, input_unklist, iostat = ios)
  !
  IF ( ios /= 0) CALL errore ('unklist', 'reading input_unklist namelist', abs(ios))
  !
  CALL unfold_read_cards(5)
  !
  !
  ! ... set up atomic positions and crystal lattice
  !
  call cell_base_init ( ibrav, celldm, a, b, c, cosab, cosac, cosbc, &
                        trd_ht, rd_ht, cell_units )
  !
  ! Init unfold data
  !

  !=============================================
  ! Invert trmat to get reciprocal space vectors and print summary
  CALL invmat(3,trmat,itrmat,det)
  !=============================================

  ! Print Det
  WRITE( stdout, '(5X, "The determinant of TRMAT is: ", f5.1)') det
      !Original zones:
  WRITE( stdout, '(5X, &
          &     "Supercell crystal axes: (cart. coord. in units of SuperCell_alat)",/, &
          &       3(15x,"a(",i1,") = (",3f11.6," )  ",/ ) )')  (apol,  &
          (at (ipol, apol) , ipol = 1, 3) , apol = 1, 3)

  WRITE( stdout, '(5x, &
      &   "Supercell reciprocal axes: (cart. coord. in units 2 pi/SC_alat)",/, &
      &            3(15x,"b(",i1,") = (",3f10.6," )  ",/ ) )')  (apol,&
      &  (bg (ipol, apol) , ipol = 1, 3) , apol = 1, 3)



  DO i=1,3
      DO j=1,3
          nat(i,j) = itrmat(j,1) * at(i,1) + itrmat(j,2) * at(i,2) + itrmat(j,3) * at(i,3)
      ENDDO
  ENDDO

  CALL recips( nat(1,1), nat(1,2), nat(1,3), nbg(1,1), nbg(1,2), nbg(1,3) )

  !New zones:
  WRITE( stdout, '(5X, &
      &     "New crystal axes: (cart. coord. in units of SC_alat)",/, &
      &       3(15x,"a(",i1,") = (",3f11.6," )  ",/ ) )')  (apol,  &
      (nat (ipol, apol) , ipol = 1, 3) , apol = 1, 3)

  WRITE( stdout, '(5x, &
      &   "New reciprocal axes: (cart. coord. in units 2 pi/SC_alat)",/, &
      &            3(15x,"b(",i1,") = (",3f10.6," )  ",/ ) )')  (apol,&
      &  (nbg (ipol, apol) , ipol = 1, 3) , apol = 1, 3)



  !=============================================
  ! Get points in small BZ
  !
  ! Transform to tpiba units if needed
  IF ( unk_points == 'crystal' .or. unk_points == 'crystal_b' ) THEN
      CALL cryst_to_cart(unnkstot, unxk, nbg, 1)
  ELSE
      DO ipkp = 1,unnkstot
          ! transform vectors to new BZ
          DO i = 1, 3
              pkp(i) = trmat(i,1)*unxk(1,ipkp) + trmat(i,2)*unxk(2,ipkp) + &
                       trmat(i,3)*unxk(3,ipkp)
          ENDDO
          unxk(:,ipkp) = pkp
      ENDDO
  ENDIF
  !

  ALLOCATE (rpkp(3, unnkstot))
  ALLOCATE (G(3, unnkstot))

  CALL get_klist(nrpkp,rpkp,G,symreduce,uniq,.false.)

  !
  ! Output points
  write (stdout,'(3X, "List of K points for supercell band calculation (NOT reduced by Symmetry)")')
  write (stdout,'(3X, "K_POINTS tpiba")')
  write (stdout,'(3X, i3)') nrpkp
  DO ipkp = 1,nrpkp
      write (stdout,'(5x,3f12.7," ",i4)') rpkp(:,ipkp), ipkp
      IF ( printg ) write (stdout,'(5x, "G: ", 3f12.7)') G(:,ipkp)
  ENDDO
  !
  ! To crrystal
  CALL cryst_to_cart(nrpkp, rpkp, at, -1)
  write (stdout,'(3X, "K_POINTS crystal")')
  write (stdout,'(3X, i3)') nrpkp
  DO ipkp = 1,nrpkp
      write (stdout,'(5x,3f12.7," ",i4)') rpkp(:,ipkp), ipkp
      IF ( printg ) write (stdout,'(5x, "G: ", 3f12.7)') G(:,ipkp)
  ENDDO

  DEALLOCATE (rpkp,G)
#ifdef __MPI
  CALL mp_global_end()
#endif
END PROGRAM unklist

