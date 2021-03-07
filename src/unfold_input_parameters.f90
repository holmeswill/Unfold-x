!
! Copyright (C) 2013 Pietro Bonfa' and
! Copyright (C) 2002-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!=----------------------------------------------------------------------------=!
!
MODULE unfold_input_parameters
!
!=----------------------------------------------------------------------------=!
!
!
!=----------------------------------------------------------------------------=!
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : nsx
  !
  IMPLICIT NONE
  !
  SAVE
  !
!=----------------------------------------------------------------------------=!
! BEGIN manual
!
!
! * DESCRIPTION OF THE INPUT FILE
!  (to be given as standard input)
!

        !
! ----------------------------------------------------------------------

!    TRMAT - UNFOLD
       REAL(DP)  :: trmat(3,3) = 0.0_DP
       LOGICAL   :: tunfold = .false.

!
!    UNKPTS - UNFOLD
!
! ...   k-points inputs
        LOGICAL :: tunkpts_inp = .false.
        REAL(DP), ALLOCATABLE :: unxk(:,:), unwk(:)
        INTEGER :: unnkstot = 0, unnk1 = 0, unnk2 = 0, unnk3 = 0, unk1 = 0, unk2 = 0, unk3 = 0
        CHARACTER(len=80) :: unk_points = 'gamma'
          ! unk_points = 'crystal' | 'tpiba' | 'gamma'*
          ! unk_points = 'crystal_b' | 'tpiba_b'
          ! select the unk points mesh
          ! 'crystal'    k points mesh is given in stdin in scaled units
          ! 'tpiba'      k points mesh is given in stdin in units of ( 2 PI / alat )
          ! 'gamma'      only gamma point is used ( default in CPMD simulation )
          ! _b means that a band input is given. The weights is a integer
          !  number that gives the number of points between the present point
          !  and the next. The weight of the last point is not used.
  !
  !    CELL_PARAMETERS - UNFOLD
  !
  LOGICAL   :: tcell = .false.
  CHARACTER(10) :: cell_units = 'bohr'
  LOGICAL:: trd_ht = .false.
  REAL(DP) :: rd_ht (3,3)
  !
PUBLIC :: deallocate_unfold_kpoints

CONTAINS
  !
  SUBROUTINE deallocate_unfold_kpoints()
    !
    IF ( allocated( unxk ) ) DEALLOCATE( unxk )
    IF ( allocated( unwk ) ) DEALLOCATE( unwk )
    !
    RETURN
    !
  END SUBROUTINE deallocate_unfold_kpoints
  !
!=----------------------------------------------------------------------------=!
!
END MODULE unfold_input_parameters
!
!=----------------------------------------------------------------------------=!
