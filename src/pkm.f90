!
! Copyright (C) 2013 Pietro Bonfa'
! This file is distributed under the terms of the
! GNU General Public License.
! See http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE pkm(itrmat, pkp, pG, ispin, any_uspp, verbose, pkml, ens, foundk)
    !----------------------------------------------------------------------------
    !
    ! ... Calculates P_Km coefficients (See Eq. 15 - PRB 85 085201)
    ! ...
    ! ... itrmat: transformation matrix for the two BZ (Eq. 2 - article above)
    ! ... pkp  : K point shifted from the large BZ to the small BZ
    ! ... pG   : G vector providing the shift into the small BZ
    ! ... pkml : array containing the coefficients for each band of each spin
    ! ... ens  : energies of each band at kpoint pkp
    !
    USE kinds, ONLY: DP
    USE cell_base, ONLY: at
    USE constants, ONLY: rytoev
    USE gvect, ONLY: g
    USE klist , ONLY: nks, xk, ngk, igk_k
    USE lsda_mod, ONLY : lsda, isk
    USE wvfct, ONLY: npw, npwx, nbnd, et
    USE control_flags, ONLY : gamma_only
    USE uspp, ONLY: nkb, vkb
    USE io_global, ONLY: stdout
    USE io_files, ONLY: nwordwfc, iunwfc
    USE wavefunctions, ONLY : evc
    USE mp_global, ONLY: intra_pool_comm
    USE mp, ONLY: mp_sum
    USE becmod, ONLY: becp, calbec, allocate_bec_type, deallocate_bec_type
    USE basic_algebra_routines, ONLY : dot_product_
    IMPLICIT NONE
    REAL(DP), INTENT(in) :: itrmat(3,3)
    REAL(DP), INTENT(in) :: pkp(3)
    REAL(DP), INTENT(in) :: pG(3)
    INTEGER,  INTENT(in) :: ispin
    LOGICAL,  INTENT(in) :: verbose
    LOGICAL,  INTENT(in) :: any_uspp
    REAL(DP), INTENT(out) :: pkml(nbnd)
    REAL(DP), INTENT(out) :: ens(nbnd)
    LOGICAL,  INTENT(out) :: foundk
    !
    INTEGER :: ibnd, j, ig, ik, ikk
    !
    REAL(DP) :: P_Km
    REAL (DP) :: gdpmx,gdpmy,gdpmz
    REAL (DP) :: cG(3) ! crystal g vectors
    COMPLEX(DP), ALLOCATABLE :: aux(:,:) ! for USPP
    !
    gdpmx = 0.d0
    gdpmy = 0.d0
    gdpmz = 0.d0
    P_Km  = 0.d0
    pkml(:)=0
    !
    IF ( verbose ) write (stdout,'(5x,"Using K = ",3f12.7,", spin=",i1)') pkp, ispin
    IF ( verbose ) write (stdout,'(5x,"Using G = ",3f12.7,", spin=",i1)') pG, ispin
    !
    ! Here we seek the k-vector pkp among the kvectors
    ! stored by pwscf during the previous band structure calculation.
    foundk = .false.
    !
    DO ik = 1, nks

        IF ( lsda ) THEN
           IF ( isk(ik) .ne. ispin ) CYCLE
        ENDIF

        IF ((abs(pkp(1) - xk (1, ik))<1.d-5) .and. &
            (abs(pkp(2) - xk (2, ik))<1.d-5) .and. &
            (abs(pkp(3) - xk (3, ik))<1.d-5)) THEN
            ikk = ik
            foundk = .true.
            EXIT
        ENDIF
    ENDDO

    !if k vector was not found exit!
    IF ( .not. foundk ) THEN
        ens(:) = 0
        pkml(:) = 0
        RETURN
    ENDIF

    CALL allocate_bec_type ( nkb, nbnd, becp ) ! Is this needed at every loop?

    ! read wfc
    CALL davcio (evc, 2*nwordwfc, iunwfc, ikk, - 1)

    ! If needed allocate resources for USPP
    IF ( any_uspp ) THEN
        !write (stdout,*) 'npwx', npwx, 'npw', npw
        ALLOCATE( aux( npwx, nbnd ) )
        ! This is needed when mixing US with other Pseudo
        CALL init_us_1
        !
        IF ( verbose ) write (stdout,*) '---> init_us'
        CALL init_us_2 (ngk(ikk), igk_k(1, ikk), xk (1, ikk), vkb)
        IF ( verbose ) write (stdout,*) '---> done init_us'
        !
        IF ( verbose ) write (stdout,*) '---> calbec and s_psi'
        CALL calbec ( npw, vkb, evc, becp )
        CALL s_psi ( npwx, npw, nbnd, evc, aux )
        IF ( verbose ) write (stdout,*) '---> done calbec and s_psi'

    ENDIF



    DO ibnd = 1, nbnd
            ens(ibnd) = et(ibnd,ikk)*rytoev
            !
            ! calculate the P_Km
            !
            DO j = 1, ngk(ikk)

                ! Transform pG vector to cartesian coordinates.
                ! The code below is equivalent to: cG = g(:,igk_k(j,ikk)) ; CALL cryst_to_cart(1, cG, bg, 1)
                DO ig = 1, 3
                    cG(ig) = at(1,ig)*g(1,igk_k(j,ikk)) + at(2,ig)*g(2,igk_k(j,ikk)) + &
                       at(3,ig)*g(3,igk_k(j,ikk))
                ENDDO
                !
                ! multiply by the inverse of the transformation matrix to move to big BZ
                gdpmx = dot_product_(itrmat(1,:), cG - real(pG) )
                gdpmy = dot_product_(itrmat(2,:), cG - real(pG) )
                gdpmz = dot_product_(itrmat(3,:), cG - real(pG) )
                !
                ! Check if this is a g vector of the big BZ (if not it has a fractional coordinate)
                IF ((abs(gdpmx - anint(gdpmx)) < 1.d-8) .AND. &
                    & (abs(gdpmy - anint(gdpmy)) < 1.d-8) .AND. &
                    & (abs(gdpmz - anint(gdpmz)) < 1.d-8)) THEN
                    !
                    IF ( any_uspp ) THEN
                        IF(gamma_only)THEN
                            P_Km = P_Km + conjg(evc(j,ibnd)) * aux(j,ibnd)
                        ELSE
                            P_Km = P_Km + conjg(evc(j,ibnd)) * aux(j,ibnd)
                        ENDIF
                    ELSE
                        IF(gamma_only)THEN
                            P_Km = P_Km +  2*conjg(evc(j,ibnd)) * evc(j,ibnd)
                        ELSE
                            P_Km = P_Km + conjg(evc(j,ibnd)) * evc(j,ibnd)
                        ENDIF
                    ENDIF

                ENDIF
            ENDDO
            !
            CALL mp_sum(P_Km, intra_pool_comm)
            !
            IF ( verbose ) write (stdout, '(8x,"PK_m = ",f12.7,", en=",f12.7)')  P_Km,  ens(ibnd)
            pkml(ibnd) = P_Km
            P_Km = 0.d0
    ENDDO

    CALL deallocate_bec_type (becp)
    IF ( any_uspp ) DEALLOCATE( aux )

END SUBROUTINE pkm


