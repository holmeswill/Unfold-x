!
! Copyright (C) 2013 Pietro Bonfa'
! This file is distributed under the terms of the
! GNU General Public License.
! See http://www.gnu.org/copyleft/gpl.txt .
!

SUBROUTINE wxml(ispin, ipkp, nbnd, pkml, ens, filout)
    !
    ! This routine write P_Km coefficients to XML files
    !
    USE kinds, ONLY: DP
    USE unfold_data, ONLY : xk
    USE xmltools,    ONLY : xml_openfile, xml_closefile, xmlw_writetag,&
                            xmlw_opentag, xmlw_closetag, add_attr
    implicit none

    INTEGER, INTENT(in) :: ispin, ipkp, nbnd
    REAL(DP), INTENT(in) :: pkml(nbnd)
    REAL(DP), INTENT(in) :: ens(nbnd)
    CHARACTER(len=256), INTENT(in) :: filout
    !
    CHARACTER(len=256) :: dirname, filename
    CHARACTER(len=9)   :: auxfname
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    !
    INTEGER :: iunout
    !
    dirname = TRIM( filout ) // '.save/'
    !
    ! ... create the k-points FILENAMES
    !
    WRITE(auxfname,'("eigenval",I1)') ispin
    filename = TRIM(dirname) // auxfname // TRIM(int_to_char(ipkp)) // '.xml'
    !
    iunout = xml_openfile ( filename )
    IF ( iunout == -1 ) RETURN
    !
    CALL xmlw_opentag ( "K-POINT" )
    !
    CALL xmlw_writetag("K-POINT_COORDS", xk(:,ipkp)) !CALL iotk_write_dat( iunout, "K-POINT_COORDS", xk(:,ipkp), COLUMNS=3 )
    !
    CALL add_attr ( "nbnd", nbnd )
    CALL add_attr ( "ik", ipkp )
    CALL add_attr ( "ispin", ispin )
    !
    CALL add_attr ( "UNITS", "eV" )
    !
    CALL xmlw_writetag( "EIGENVALUES", ens(:) )
    !
    !
    CALL xmlw_writetag( "PROJECTIONS", pkml(:) )
    !
    CALL xmlw_closetag ( )
    !
    CALL xml_closefile ( )
    !
END SUBROUTINE wxml

