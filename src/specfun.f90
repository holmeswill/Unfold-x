!
! Copyright (C) 2013 Pietro Bonfa'
! This file is distributed under the terms of the
! GNU General Public License.
! See http://www.gnu.org/copyleft/gpl.txt .
!

FUNCTION specfun(ens, pkml, E, w) RESULT(r)
    !
    ! This subroutine calculates the specral function
    ! from P_Km coefficients
    !

    USE kinds,     ONLY: DP
    USE wvfct,     ONLY: nbnd
    USE constants, ONLY : pi

    IMPLICIT NONE
    REAL(DP), INTENT(in)  :: ens(nbnd)
    REAL(DP), INTENT(in)  :: pkml(nbnd)
    REAL(DP), INTENT(in)  :: E
    REAL(DP), INTENT(in)  :: w
    REAL(DP) :: r
    INTEGER :: ibnd
    r = 0
    DO ibnd = 1, nbnd
        IF (ens(ibnd) > (E - 6*w) .and. ens(ibnd) < (E + 6*w) ) THEN
            r = r + pkml(ibnd)*(1.d0/(w*sqrt(pi)))*exp(-(((E-ens(ibnd))**2)/(w**2))) ! pkml(ibnd)
        ENDIF
    ENDDO

END FUNCTION specfun

