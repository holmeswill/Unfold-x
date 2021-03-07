!
! Copyright (C) 2013 Pietro Bonfa' and
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License.
! See http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
PROGRAM unfold
  !-----------------------------------------------------------------------
  !
  !    Program for unfolding band structures.
  !    The two basic steps are:
  !    1) produce a list of k-points for 'bands' in the supercell
  !    2) run the program to unfold the band structure obtained
  !       with step 1.
  !
  !    DESCRIPTION of the INPUT : see file Doc/INPUT_UN.*
  !
  USE io_global,  ONLY : ionode
  USE mp_global,  ONLY : mp_startup
  USE environment,ONLY : environment_start
  !
  IMPLICIT NONE
  !
  ! initialise environment
  !
#ifdef __MPI
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'UNFOLD' )
  !
  IF ( ionode )  CALL input_from_file ( )
  !
  CALL run()
  !
  CALL stop_unfold()
  !
END PROGRAM unfold
!
!-----------------------------------------------------------------------
SUBROUTINE run ()
  !-----------------------------------------------------------------------
  !
  !    This subroutine reads the data for the output file produced by pw.x
  !    extracts and calculates the Spectral Function and
  !    writes it to a binary file for plotting.
  !
  !    DESCRIPTION of the INPUT: see README file
  !
  USE read_cards_module, ONLY : read_cards
  USE io_global, ONLY: stdout
  USE unfold_input_parameters, ONLY : trmat, unxk, unnkstot, unk_points
  USE kinds,     ONLY : DP
  USE lsda_mod, ONLY: lsda, nspin
  USE noncollin_module, ONLY : noncolin
  !USE ener,      ONLY : ef
  USE ions_base, ONLY : nat, ntyp=>nsp, ityp, tau
  !USE uspp,            ONLY : nkb, vkb
  USE uspp_param,      ONLY : upf
  USE klist,     ONLY : nks
  USE io_files,  ONLY : tmp_dir, prefix
  USE io_global, ONLY : ionode, ionode_id
  USE mp_global,     ONLY : nproc_pool, nproc_file, nproc_pool_file
  !USE noncollin_module, ONLY : i_cons
  !USE paw_variables, ONLY : okpaw
  USE mp,        ONLY : mp_bcast, mp_sum
  USE mp_world,         ONLY : world_comm
  USE mp_pools,         ONLY : npool, me_pool, root_pool, inter_pool_comm
  !USE constants, ONLY : rytoev
  USE wvfct, ONLY: nbnd, et
  USE io_files, ONLY : create_directory
  USE unfold_data, ONLY : unfold_data_init, get_klist, itrmat
  USE unfold_read_cards_module, ONLY : unfold_read_cards

  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  !
  REAL(DP), ALLOCATABLE :: pkml(:) !al pkml
  REAL(DP), ALLOCATABLE :: ens(:) !energies
  REAL(DP), ALLOCATABLE :: rspf(:,:,:) !response function
  REAL(DP), ALLOCATABLE :: rpkp (:,:) !reduced kpoints
  REAL(DP), ALLOCATABLE :: nG (:,:) !reduced kpoints
  REAL(DP), ALLOCATABLE :: E(:)
  INTEGER,  ALLOCATABLE :: kdone(:,:)
  !
  REAL(DP), EXTERNAL :: specfun

  INTEGER :: ios
  INTEGER :: nrpkp
  INTEGER :: ipkp, ispin !, npkp
  LOGICAL :: kfound

  LOGICAL :: any_uspp

!!! Input parameters

  CHARACTER(len=256) :: outdir   ! directory for temporary files
  CHARACTER(len=256) :: filout

  CHARACTER(len=256) :: dirname
  CHARACTER(len=25) :: kpathunit
  REAL(DP) :: Emin = 4.0d0
  REAL(DP) :: Emax = 16.0d0
  REAL(DP) :: DeltaE = 0.1d0
  REAL(DP) :: w = 0.25d0
  REAL(DP) :: E_unset=1000000.d0
  LOGICAL :: nscfklist = .false.
  LOGICAL :: symreduce = .true.
  LOGICAL :: verbose = .false.
  LOGICAL :: write_pkm = .false. !! WARNING: this will write a lot of files!
  INTEGER :: ispfp = 0
  INTEGER :: nspfp

  NAMELIST / inputun / outdir, prefix, Emin, Emax, DeltaE, w, filout, nscfklist, symreduce, verbose, write_pkm, kpathunit



  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  IF ( trim( filout ) == ' ' ) filout = './out.dat'
  !
  ios = 0
  !
  IF ( ionode )  THEN
     !
     !     reading the namelist inputpp
     !
     READ (5, inputun, iostat = ios)
     !
     tmp_dir = trimcheck ( outdir )
     !
  ENDIF
  !
  CALL unfold_read_cards(5)


  CALL mp_bcast (ios, ionode_id, world_comm )
  !
  IF ( ios /= 0) CALL errore ('unfold', 'reading inputun namelist', abs(ios))

  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id , world_comm )
  CALL mp_bcast( prefix, ionode_id , world_comm )
  CALL mp_bcast( trmat, ionode_id , world_comm )
  CALL mp_bcast( Emin, ionode_id , world_comm )
  CALL mp_bcast( Emax, ionode_id , world_comm )
  CALL mp_bcast( DeltaE, ionode_id , world_comm )
  CALL mp_bcast( w, ionode_id , world_comm )
  CALL mp_bcast( filout, ionode_id , world_comm )
  CALL mp_bcast( nscfklist, ionode_id , world_comm )
  CALL mp_bcast( symreduce, ionode_id , world_comm )
  CALL mp_bcast( write_pkm , ionode_id , world_comm )
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file ( )
  !
  IF (noncolin) CALL errore('unfold',&
       'Non-collinear not implemented',1)
  !
  IF (write_pkm .and. (npool > 1)) CALL errore('unfold',&
       'Cannot save PKM files when using pool parallelism',1)

  CALL openfil_pp ( )

  !=============================================
  ! Initialize data (Brillouin zones, kpoints)
  CALL unfold_data_init()
  !=============================================

  !=============================================
  !Check USPP
  any_uspp = any(upf(1:ntyp)%tvanp)
  !
  if(any_uspp) then
     ! WARNING - I'm not expert enought to guarantee that USPS will work
     write (stdout,'(5X, "WARNING: USPP blindly implemented!")')
  end if
  !=============================================

  !=============================================
  !NSCFKPOINTS (deprecated)
  !
  !If user wants k point, just print the kpoint list and exit
   IF ( nscfklist .and. ionode ) THEN
        CALL print_klist(unnkstot,symreduce,verbose)
        RETURN
   ELSEIF ( nscfklist ) THEN
        RETURN
   ENDIF
  !=============================================

  !=============================================
  ! CALCULATE ENERGIES
  !
  ! find min and max energy for plot (band extrema if not set)
  !
  IF ( Emin == -E_unset ) THEN
      Emin = MINVAL ( et(1, 1:nks) )
      IF ( w > 0.0_dp ) Emin = Emin - 3.0_dp * w
  END IF
  IF ( Emax  == E_unset ) THEN
      Emax = MINVAL ( et(nbnd, 1:nks) )
      IF ( w > 0.0_dp ) Emax = Emax + 3.0_dp * w
  END IF
  !
  ! calculate number of points
  nspfp = nint ( (Emax - Emin) / DeltaE+0.500001d0)
  !
  ! prepare vector of energies
  ALLOCATE(E(nspfp))
  DO ispfp= 1, nspfp
      E(ispfp) = Emin + (ispfp - 1) * DeltaE
  ENDDO
  !=============================================

  !=============================================
  ! beeing nice to users is good :)
  write (stdout,'(5X, "Plotting command for Gnuplot")')
  write (stdout,'(5X, "set pm3d map interpolate 2,2")')
  write (stdout,'(5X, "splot ''",a,"'' binary record=(",i4,",-1) format=''%double'' u 1:2:3")') trim(filout), nspfp
  !=============================================


  !=============================================
  !!! GET REEUCED POINTS
  !
  ! Here the basecell k points (given in input) are mapped to the
  ! supercell BZ. Symmetry is eventually used to reduce the number
  ! of kpoints needed.
  !
  ALLOCATE (rpkp(3,unnkstot))
  ALLOCATE (nG(3,unnkstot))
  !
  CALL get_klist(nrpkp,rpkp,nG,symreduce,.false.,verbose)
  !
  !=============================================

  !=============================================
  !!! CREATE DIRECTORY FOR PKML FILES (if needed)
  !
  IF ( write_pkm ) THEN
    dirname = TRIM( filout ) // '.save'
    !
    ! ... create the main saving directory
    !
    CALL create_directory( dirname )
  ENDIF
  !=============================================


  !=============================================
  !!! Calculate P_Km coefficients
  ALLOCATE (pkml(nbnd))           ! P_Km List : P_Km for each m band
  ALLOCATE (ens(nbnd))            ! EigeNvalueS
  ALLOCATE (rspf(nspfp,unnkstot,nspin)) ! ResultsSPectralFunction(Skpoints, Number of SPectral Function Points)
  ALLOCATE (kdone(nrpkp,nspin))

  rspf = 0.d0
  kdone(:,:) = 0

  DO ispin = 1, nspin
    ens  = 0.d0
    !
    DO ipkp = 1,nrpkp
        write (stdout,'(5x,"Doing K = ",3f12.7,", spin=",i1)') unxk(:,ipkp), ispin
        ! Get P_Km for point ipkp. Output to pkml, also eigenvalues are written to ens.
        CALL pkm (itrmat, rpkp(:,ipkp), nG(:,ipkp), ispin, any_uspp, verbose, pkml, ens, kfound)
        !
        IF (.not. kfound) THEN
            IF (npool == 1) write (stdout,*) '---> WARNING: couldn t find corresponding vector for' , rpkp(:,ipkp)
            IF (npool == 1) write (stdout,*) '---> Giving up with this point'
            CYCLE
        ENDIF
        !
        IF ( write_pkm .and. ionode ) CALL wxml(ispin, ipkp, nbnd, pkml, ens, filout)
        !
        DO ispfp= 1, nspfp !loop over energies from Emin to Emax (Number of SPectral Function Points)
            rspf(ispfp,ipkp,ispin) = specfun(ens, pkml, E(ispfp), w)
        ENDDO
        !
        kdone(ipkp, ispin) = 1
    ENDDO
    !
  ENDDO

  !=============================================
  !!! CHECK AND COLLECT DATA
  !
  ! Check missing points
  CALL mp_sum( kdone, inter_pool_comm )
  IF ( ANY(kdone == 0) ) &
  WRITE (stdout,*) '---> WARNING: couldn t find corresponding vector for' ,  COUNT(kdone==0), ' k points'
  !
  ! Some points may be computed twice when pool parallelism is active
  IF ( ANY(kdone > 1) ) THEN
      WRITE (stdout,*) '---> INFO: ' ,  COUNT(kdone>1), ' k points were calculated multiple times.'
      WRITE (stdout,*) '---> INFO: this will not affect the results, only the performance of the code, slightly.'
      WRITE (stdout,*) '---> INFO: Sorry about that.'
      DO ispin = 1, nspin
          DO ipkp = 1,nrpkp
              IF ( kdone(ipkp, ispin) > 1 ) rspf(:,ipkp,ispin) = rspf(:,ipkp,ispin) / DBLE(  kdone(ipkp, ispin) )
          ENDDO
      ENDDO
  ENDIF
  !
  ! add contributions from different pools
  CALL mp_sum( rspf, inter_pool_comm )
  !
  ! Write data to disk
  DO ispin = 1, nspin
      IF ( ionode ) CALL wbin2(nspfp, unnkstot, rspf(:,:,ispin), E, unxk, ispin, filout, kpathunit)
  ENDDO
  !=============================================
  write (stdout,'(5x,"Done!")')
  !
  DEALLOCATE (pkml, ens, rspf, rpkp, nG, E, kdone)
  !
END SUBROUTINE run

SUBROUTINE print_klist(unnkstot,symreduce, verbose)
  USE io_global, ONLY: stdout
  USE kinds,     ONLY: DP
  USE unfold_data, ONLY : get_klist

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: unnkstot
  LOGICAL, INTENT(IN) :: symreduce
  LOGICAL, INTENT(IN) :: verbose

  REAL(DP), ALLOCATABLE :: rpkp (:,:) !reduced kpoints
  REAL(DP), ALLOCATABLE :: nG (:,:) !reduced kpoints

  INTEGER :: nrpkp
  INTEGER :: ipkp !, npkp
          ALLOCATE (rpkp(3,unnkstot))
          ALLOCATE (nG(3,unnkstot))
          IF (symreduce) THEN
            write (stdout,'(3X, "List of K points for supercell band calculation (Reduced by Symmetry)")')
          ELSE
            write (stdout,'(3X, "List of K points for supercell band calculation (NOT reduced by Symmetry)")')
          ENDIF
          CALL get_klist(nrpkp,rpkp,nG,symreduce,.true.,verbose)
          write (stdout,'(3X, "K_POINTS tpiba")')
          write (stdout,'(3X, i3)') nrpkp
          DO ipkp = 1,nrpkp
              write (stdout,'(5x,3f12.7," ",i4)') rpkp(:,ipkp), ipkp
          ENDDO
          DEALLOCATE (rpkp)
          DEALLOCATE (nG)
END SUBROUTINE print_klist


SUBROUTINE stop_unfold
  !--------------------------------------------------------------------
  !
  ! Synchronize processes before stopping.
  !
  USE io_files, ONLY: iunwfc
  USE mp_global, ONLY: mp_global_end
  USE parallel_include
  USE unfold_data , ONLY: deallocate_unfold_data
  IMPLICIT NONE
#ifdef __MPI

  INTEGER :: info
  LOGICAL :: op

  INQUIRE ( iunwfc, opened = op )

  IF ( op ) THEN
     CLOSE (unit = iunwfc, status = 'keep')
  ENDIF

  CALL mp_global_end()

#endif
  CALL deallocate_unfold_data()
  STOP
END SUBROUTINE stop_unfold



