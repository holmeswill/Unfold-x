!
! Copyright (C) 2013 Pietro Bonfa'
! This file is distributed under the terms of the
! GNU General Public License.
! See http://www.gnu.org/copyleft/gpl.txt .
!

SUBROUTINE wbin(ispin, ipkp, nbnd, pkml, ens, filout)
    !
    ! This routine write P_Km coefficients to binary files
    !
    USE kinds, ONLY: DP
    implicit none

    INTEGER, INTENT(in) :: ispin, ipkp, nbnd
    REAL(DP), INTENT(in) :: pkml(nbnd)
    REAL(DP), INTENT(in) :: ens(nbnd)
    CHARACTER(len=256), INTENT(in) :: filout
    !
    REAL(DP), ALLOCATABLE :: outmat(:,:)
    CHARACTER(len=256) :: auxfilout
    !real, array dimension
    integer :: irec

    write(auxfilout,'(a,"_K",I3.3,"_S",I1)') trim(filout), ipkp, ispin

    ALLOCATE(outmat(2,nbnd))
    outmat(1,:) = pkml
    outmat(2,:) = ens
! ... initialize array ...

    inquire( iolength=irec ) outmat
    open( 36, file=auxfilout, form='unformatted', access='direct', recl=irec )
    write( 36, rec=1 ) outmat
    close( 36, status='keep' )
    DEALLOCATE(outmat)

END SUBROUTINE wbin

SUBROUTINE wbin2(m, n, spr,E,unxk,ispin,filout,kpathunits)
    !
    ! This routine writes the spectral function to binary file(s)
    !
    USE kinds, ONLY: DP
    USE unfold_data, ONLY: xk, nbg, k_points, nkstot
    USE lsda_mod, ONLY: nspin
    implicit none

    INTEGER, INTENT(in) :: n
    INTEGER, INTENT(in) :: m

    REAL(DP), INTENT(in) :: spr(m, n)
    REAL(DP), INTENT(in) :: E(m)
    REAL(DP), INTENT(in) :: unxk(3,n)
    INTEGER, INTENT(in) :: ispin
    CHARACTER(len=256), INTENT(in) :: filout
    CHARACTER(len=25), INTENT(in) :: kpathunits
    CHARACTER(len=256) :: auxfilout
    !
    ! FOR K PATH CALC
    REAL(DP) :: dt = 0
    REAL(DP) :: vl(3)
    REAL(DP) :: r = 0
    REAL(DP), ALLOCATABLE :: kpathlen(:)
    INTEGER :: switch = 0 !whether coordinates should be changed or not
    !

    ! For Gnuplot Matrix allocation and writing
    INTEGER :: irec, stat
    REAL(DP), ALLOCATABLE :: outmat(:)
    !
    INTEGER :: i, j, ind, ipkp


    ! get filename
    IF (nspin > 1) THEN
        write(auxfilout,'(I1)') ispin
        auxfilout = trim(filout)//auxfilout
    ELSE
        auxfilout = filout
    ENDIF


    !=============================================
    ! CALCULATE KPATH LENGTH
    !
    IF ( (k_points == 'crystal' .or. k_points == 'crystal_b') .and. &
         (kpathunits == 'tpiba') ) switch = 1
    IF ( (k_points == 'tpiba' .or. k_points == 'tpiba_b') .and. &
         (kpathunits == 'crystal') ) switch = -1
    !
    ! change coordinates if needed
    ! switch = 1 -> go to cartesian cartesian coordinates
    ! switch =-1 -> go to crystal cartesian coordinates
    !
    ! only if spin == 1 otherwise it's applied twice!
    IF ( ( switch /= 0 ) .and. ( ispin == 1 ) ) CALL cryst_to_cart (nkstot, xk, nbg, switch)
    !
    !
    ALLOCATE (kpathlen(n))
    dt =  0.d0
    ! calculate segments length
    DO ipkp = 1,n-1
       kpathlen(ipkp) = dt
       vl(:)=xk(:,ipkp+1)-xk(:,ipkp)
       r=sqrt(vl(1)**2+vl(2)**2+vl(3)**2)
       dt=dt+r
    ENDDO
    kpathlen(n) = dt
    !=============================================


    ALLOCATE(outmat(3*n*m))
    ind = 1
    DO i=1,n
        DO j=1,m
            outmat(ind)=kpathlen(i)
            ind = ind + 1
            outmat(ind)=E(j)
            ind = ind + 1
            outmat(ind)=spr(j,i)
            ind = ind + 1
        ENDDO
    ENDDO

    !delete file if it's already present
    open(unit=36, iostat=stat, file=auxfilout, status='old')
    if (stat.eq.0) close(36, status='delete')


    inquire( iolength=irec ) outmat
    open( 36, file=auxfilout, form='unformatted', access='direct', recl=irec )
    write( 36, rec=1 ) outmat
    close( 36, status='keep' )
    DEALLOCATE(outmat, kpathlen)
END SUBROUTINE wbin2
