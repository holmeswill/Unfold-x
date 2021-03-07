!
! Copyright (C) 2013 Pietro Bonfa' and
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!------------------------------------------------------------------------------!
  MODULE unfold_data
!------------------------------------------------------------------------------!

    USE kinds, ONLY : DP
    USE constants, ONLY : pi, bohr_radius_angs
    USE io_global, ONLY : stdout
    USE cell_base, ONLY : at, bg
    USE matrix_inversion, ONLY : invmat
    USE unfold_input_parameters, ONLY : trmat, unxk, unnkstot, unk_points
!
    IMPLICIT NONE
    SAVE
    !
    PRIVATE
    !
    PUBLIC :: nat, nbg, itrmat, det, xk, nkstot, k_points
    !
    !  alat: lattice parameter - often used to scale quantities, or
    !  in combination to other parameters/constants to define new units
    !REAL(DP) :: alat = 0.0_DP
    ! omega: volume of the simulation cell
    !REAl(DP) :: omega = 0.0_DP
    ! tpiba: 2 PI/alat, tpiba2=tpiba^2
    !REAL(DP) :: tpiba  = 0.0_DP, tpiba2 = 0.0_DP
    !  new direct and reciprocal lattice primitive vectors
    !  at(:,i) are the lattice vectors of the simulation cell, a_i,
    !          in alat units: a_i(:) = at(:,i)/alat
    !  bg(:,i) are the reciprocal lattice vectors, b_i,
    !          in tpiba=2pi/alat units: b_i(:) = bg(:,i)/tpiba
    REAL(DP) :: nat(3,3) = RESHAPE( (/ 0.0_DP /), (/ 3, 3 /), (/ 0.0_DP /) )
    REAL(DP) :: nbg(3,3) = RESHAPE( (/ 0.0_DP /), (/ 3, 3 /), (/ 0.0_DP /) )
    !  inverse of trmat
    REAL(DP) :: itrmat(3,3)
    !  determinant of the trmat matrix
    REAL(DP) :: det
    !  kpoints in big BZ
    REAL(DP), ALLOCATABLE :: xk(:,:)
    INTEGER :: nkstot = 0
    CHARACTER(len=80) :: k_points = 'gamma'
    !
    PUBLIC :: unfold_data_init, get_klist, deallocate_unfold_data
!------------------------------------------------------------------------------!
  CONTAINS
!------------------------------------------------------------------------------!
!
  SUBROUTINE unfold_data_init( )
    !
    ! ... initialize cell_base module variables, set up crystal lattice
    !

    IMPLICIT NONE
    !
    INTEGER :: ipkp, i, j, ipol, apol
    REAL(DP) :: temp(3)
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
    ! first save input points given for the big BZ
    ! the will be used for producing kpath segments in output
    !
    ALLOCATE (xk(3, unnkstot))
    xk = unxk
    nkstot = unnkstot
    k_points = unk_points
    !
    ! now transform to tpiba units if needed
    IF ( unk_points == 'crystal' .or. unk_points == 'crystal_b' ) THEN
        CALL cryst_to_cart(unnkstot, unxk, nbg, 1)
    ELSE
        !
        ! Now transform the points using trmat
        !
        DO ipkp = 1,unnkstot
            ! transform vectors to new BZ
            DO i = 1, 3
                temp(i) = trmat(i,1)*unxk(1,ipkp) + trmat(i,2)*unxk(2,ipkp) + &
                            trmat(i,3)*unxk(3,ipkp)
            ENDDO
            unxk(:,ipkp) = temp
        ENDDO
    ENDIF
    !
  RETURN
  !
  END SUBROUTINE unfold_data_init
!
  SUBROUTINE get_klist(nrpkp,rpkp,G,symreduce,uniq,verbose)
    !
    ! This subroutine is used to bring the kpoints back to the IBZ of
    ! the supercell (the small BZ).
    ! ...
    ! ... nrpkp     : number of reduced k points (out)
    ! ... rpkp      : kpoints (possibly) reduced by symmetry (out)
    ! ... G         : G vector providing the shift into the small BZ
    !                (and possibly changed by symmetry) (out)
    ! ... symreduce : whether to use symmetry to reduce the number of
    !                 k-points in to produce the band structure (in)
    ! ... uniq      : if true, return only inequivalent kpoints
    !                 (used only for printing)
    ! ... verbose   : makes the code rambling
    !
    USE kinds,     ONLY : DP
    !USE io_global, ONLY: stdout
    !USE input_parameters, ONLY : unxk, unnkstot, trmat
    USE cell_base, ONLY : bg, at

    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: nrpkp
    REAL(DP), INTENT(OUT) :: rpkp(3,unnkstot)
    REAL(DP), INTENT(OUT) :: G(3,unnkstot)
    LOGICAL, INTENT(IN) :: uniq
    LOGICAL, INTENT(IN) :: symreduce
    LOGICAL, INTENT(IN) :: verbose

    REAL(DP) :: pkp(3)
    INTEGER :: nG(3)
    INTEGER :: ipkp, i, j, a, iuniq

    REAL(DP), ALLOCATABLE :: auxrpkp (:,:)
    !
    ALLOCATE (auxrpkp(3, unnkstot))
    !

    DO ipkp = 1,unnkstot
        ! go to crystal coordinates
        pkp = unxk(:,ipkp)
        CALL cryst_to_cart(1,pkp,at,-1)
        !DO i = 1, 3
        !   pkp(i) = at(1,i)*unxk(1,ipkp) + at(2,i)*unxk(2,ipkp) + &
        !               at(3,i)*unxk(3,ipkp)
        !ENDDO
        nG(:) = 0 ! variable where g vectors are stored

        ! reduce to IBZ
        !
        CALL vecfbz(1.d-8,bg,pkp,nG)
        auxrpkp(:,ipkp) = pkp
        G(:,ipkp) = nG
    ENDDO
    !
    IF ( symreduce ) THEN
        CALL get_kequiv ( unnkstot, auxrpkp, nrpkp, rpkp, G, uniq, verbose)
        !go to cartesian coordinates
        CALL cryst_to_cart(nrpkp,rpkp,bg,1)

    ELSE
        ! back to 2pi/a coordinates
        iuniq = 0 ! counter for points
        DO ipkp = 1,unnkstot
            !if only unique kpoints are needed check if this kpoint is already in the list.
            IF ( uniq .eqv. .true. ) THEN
                DO a= 1, ipkp-1
                    IF ( (ABS(auxrpkp(1,ipkp) - auxrpkp(1,a)) < 1.d-8) .and. &
                         (ABS(auxrpkp(2,ipkp) - auxrpkp(2,a)) < 1.d-8) .and. &
                         (ABS(auxrpkp(3,ipkp) - auxrpkp(3,a)) < 1.d-8) )  EXIT
                ENDDO

                IF ( a < (ipkp-1) ) THEN
                    ! write (*,*) 'Skipping ', auxrpkp(:,ipkp)
                    CYCLE
                ENDIF
            ENDIF
            !
            iuniq = iuniq + 1
            rpkp(:,iuniq) = auxrpkp(:,ipkp)
            !
        ENDDO
        nrpkp = iuniq
        !go to cartesian coordinates
        CALL cryst_to_cart(nrpkp,rpkp,bg,1)
    ENDIF

    DEALLOCATE (auxrpkp)

  END SUBROUTINE get_klist
  !------------------------------------------------------------------------------!
  !
  SUBROUTINE get_kequiv ( nkr, xkg, nks, xk , G, uniq, verbose) !uniq: flag if you want only unique list
  !-----------------------------------------------------------------------
  !
  !  This subroutine checks if two kpoints are related by the symmetries of the system.
  !  Returns a list of equivalent kpoints and new G vectors obtained with
  !  the inverse symmetry which relates the two kpoints (so that we do not rotate the
  !  wfc but the direction of relevant g vectors in the unfold procedure).
  !
    USE kinds, ONLY: DP
    USE symm_base,          ONLY : s, t_rev, irt, nrot, nsym, invsym, nosym, &
                                   d1,d2,d3, time_reversal, sname, set_sym_bl, &
                                   find_sym, inverse_s, no_t_rev, invs
    !USE cell_base, ONLY : bg
    USE io_global, ONLY: stdout

    IMPLICIT NONE
    !
    INTEGER, INTENT(in):: nkr
    real(DP), INTENT(in) ::xkg(3,nkr)
    INTEGER, INTENT(out) :: nks
    real(DP), INTENT(out):: xk(3,nkr)
    real(DP), INTENT(inout):: G(3,nkr)
    LOGICAL , INTENT(in):: uniq
    LOGICAL , INTENT(in):: verbose
    ! LOCAL:
    real(DP), PARAMETER :: eps=1.0d-5
    real(DP) :: xkr(3), sg(3), fact

    !real(DP), ALLOCATABLE:: xkg(:,:), wkk(:)

    INTEGER :: i, ns, n, nk
    INTEGER, ALLOCATABLE :: equiv(:)
    LOGICAL :: in_the_list = .false.
    !
    ALLOCATE (equiv( nkr))
    !
    !  equiv(nk) =nk : k-point nk is not equivalent to any previous k-point
    !  equiv(nk)!=nk : k-point nk is equivalent to k-point equiv(nk)
    !
    DO nk=1,nkr
       equiv(nk)=nk
    ENDDO
    !
    DO nk=1,nkr
    !  check if this k-point has already been found equivalent to another

        IF (equiv(nk) == nk) THEN
            !  check if there are equivalent k-point to this in the list
            !  (excepted those previously found to be equivalent to another)
            !  check both k and -k
            DO ns=1,nrot

                DO i=1,3
                    xkr(i) = s(i,1,ns) * xkg(1,nk) &
                            + s(i,2,ns) * xkg(2,nk) &
                            + s(i,3,ns) * xkg(3,nk)
                    xkr(i) = xkr(i) - nint( xkr(i) )
                ENDDO

                IF(t_rev(ns)==1) xkr = -xkr

                DO n=nk+1,nkr
                  in_the_list = abs(xkr(1) - xkg(1,n) )<=eps .and. &
                                  abs(xkr(2) - xkg(2,n))<=eps .and. &
                                  abs(xkr(3) - xkg(3,n))<=eps

                  !write (stdout,*)   in_the_list, nk, xkr
                  IF (in_the_list) THEN

                      IF (n>nk .and. equiv(n)==n) THEN
                          IF ( verbose ) write (stdout, '("ROT - K point ", I3, " EQUIVALENT TO: ", I3)')  n, nk
                          sg(:) = s(:,1,invs(ns)) * G(1,n) + &
                          s(:,2,invs(ns)) * G(2,n) + &
                          s(:,3,invs(ns)) * G(3,n)
                          IF ( verbose ) write (stdout, '("Initial G: ", 3f10.3)') G(:,n)
                          IF ( verbose ) write (stdout, '("Rotated G: ", 3f10.3)') sg
                          G(:,n) = sg
                          equiv(n) = nk
                          EXIT
                      !ELSE
                      !    IF (equiv(n)/=n .or. n>nk ) CALL errore('kpoint_grid', &
                      !    'something wrong in the checking algorithm',1)
                      ENDIF
                  ENDIF
                ENDDO


                IF ( time_reversal ) THEN
                    !
                    DO n=nk+1,nkr
                          in_the_list = abs(-xkr(1) - xkg(1,n) )<=eps .and. &
                                      abs(-xkr(2) - xkg(2,n))<=eps .and. &
                                      abs(-xkr(3) - xkg(3,n))<=eps
                          !write (stdout,*)   in_the_list , nk, -xkr
                          IF (in_the_list) THEN

                              IF (n>nk .and. equiv(n)==n) THEN
                                  IF ( verbose ) write (stdout, '("TR - K point ", I3, " EQUIVALENT TO: ", I3)')  n, nk
                                  sg(:) = s(:,1,invs(ns)) * G(1,n) + &
                                      s(:,2,invs(ns)) * G(2,n) + &
                                      s(:,3,invs(ns)) * G(3,n)
                                  IF ( verbose ) write (stdout, '("Initial G: ", 3f10.3)') G(:,n)
                                  IF ( verbose ) write (stdout, '("Rotated G: ", 3f10.3)') sg
                                  G(:,n) = sg
                                  equiv(n) = nk
                                  EXIT
                              !ELSE
                              !write (stdout,*)   n , equiv(n), nk
                              !IF (equiv(n)/=n.or.n<nk) CALL errore('kpoint_grid', &
                              !'something wrong in the checking algorithm',2)
                              ENDIF
                          ENDIF
                   ENDDO
                ENDIF
                !
            ENDDO
        ENDIF
    ENDDO

    !  counts irreducible points and order them

    nks=0
    fact=0.0d0
    DO nk=1,nkr
       IF (equiv(nk)==nk .and. uniq .eqv. .true.) THEN
          nks=nks+1
          IF (nks>nkr) CALL errore('kequiv','too many k-points',1)
          !wk(nks) = wkk(nk)
          !fact    = fact+wk(nks)
          !  bring back into to the first BZ
          DO i=1,3
             xk(i,nks) = xkg(i,nk) !-nint(xkg(i,nk))
          ENDDO
       ELSEIF (uniq .eqv. .false.) THEN
          nks=nks+1
          !DO i=1,3
          xk(:,nk) = xkg(:,equiv(nk)) !-nint(xkg(i,nk))
          !write (stdout,('(3f5.10)')) xk(:,nk)
          !ENDDO
       ENDIF
    ENDDO
    !
    !
    DEALLOCATE(equiv)

    RETURN
  END SUBROUTINE get_kequiv

  SUBROUTINE deallocate_unfold_data ( )
      IF ( allocated( xk ) ) DEALLOCATE( xk )
  END SUBROUTINE deallocate_unfold_data

!
!------------------------------------------------------------------------------!
END MODULE unfold_data
!------------------------------------------------------------------------------!
