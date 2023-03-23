module time_constraint

  use utils
  use main_data
  contains 

  SUBROUTINE CalcTimeConstraintsForElems()
  
    IMPLICIT NONE
  
  ! evaluate time constraint
    CALL CalcCourantConstraintForElems()
  
  ! check hydro constraint
    CALL CalcHydroConstraintForElems()
  
  END SUBROUTINE CalcTimeConstraintsForElems

  SUBROUTINE CalcCourantConstraintForElems()
  
    IMPLICIT NONE
    REAL(KIND=8)    :: dtcourant
    INTEGER(KIND=4) :: COURANT_ELEM
  
    REAL(KIND=8) :: qqc, qqc2, dtf
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: dtcourant_per_thread
    INTEGER(KIND=4) :: length, threads, i, indx, thread_num
    INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: courant_elem_per_thread
  
    dtcourant    = 1.0e+20_RLK
    COURANT_ELEM = -1
  
    qqc = domain%m_qqc
    length = domain%m_numElem
  
    qqc2 = (64.0_RLK) * qqc * qqc
  
  ! For sequential run (non OpenMP, replace line above with next line.  
    threads = 1_4
  
    ALLOCATE(dtcourant_per_thread(0:threads-1))
    ALLOCATE(courant_elem_per_thread(0:threads-1))
  
    DO i = 0, threads-1
      courant_elem_per_thread(i) = -1
      dtcourant_per_thread(i) =  (1.0e+20_RLK)
    ENDDO
  
  
    DO i = 0, length-1
      indx = domain%m_matElemlist(i)
  
      dtf = domain%m_ss(indx) * domain%m_ss(indx)
  
      IF ( domain%m_vdov(indx) < (0.0_RLK) ) THEN
  
        dtf = dtf + qqc2 * domain%m_arealg(indx) * domain%m_arealg(indx)  &
                  * domain%m_vdov(indx) * domain%m_vdov(indx)
      ENDIF
  
      dtf = SQRT(dtf)
  
      dtf = domain%m_arealg(indx) / dtf
  
  !   determine minimum timestep with its corresponding elem
      IF (domain%m_vdov(indx) /= (0.0_RLK)) THEN
  
  !     For sequential run (non OpenMP, replace line above with next line.
        thread_num = (0_4)
  
  
        IF ( dtf < dtcourant_per_thread(thread_num) ) THEN
  
          dtcourant_per_thread(thread_num) = dtf
          courant_elem_per_thread(thread_num) = indx
        ENDIF
      ENDIF
  
    ENDDO
  
    DO i = 0, threads-1
      IF(dtcourant_per_thread(i) < dtcourant) THEN
        dtcourant = dtcourant_per_thread(i)
        courant_elem =  courant_elem_per_thread(i)
      ENDIF
    ENDDO
  
  
  ! Don't try to register a time constraint if none of the elements
  ! were active
    IF (courant_elem /= -1) THEN
      domain%m_dtcourant = dtcourant
    ENDIF
  
    DEALLOCATE(dtcourant_per_thread)
    DEALLOCATE(courant_elem_per_thread)
  
    RETURN
  
  END SUBROUTINE CalcCourantConstraintForElems
  
  SUBROUTINE CalcHydroConstraintForElems()
  
    IMPLICIT NONE 
  
    REAL(KIND=8) :: dthydro
    REAL(KIND=8) :: dvovmax, dtdvov
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dthydro_per_thread
    INTEGER(KIND=4) :: hydro_elem
    INTEGER(KIND=4) :: threads, i, length, indx, thread_num
  
    INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: hydro_elem_per_thread
  
    dthydro = 1.0e+20_RLK
    hydro_elem = -1
    dvovmax = domain%m_dvovmax
    length = domain%m_numElem
  
  ! For sequential run (non OpenMP, replace line above with next line.
    threads = (1_4)
  
    ALLOCATE(dthydro_per_thread(0:threads-1))
    ALLOCATE(hydro_elem_per_thread(0:threads-1))
  
    DO i = 0, threads-1
      hydro_elem_per_thread(i) = hydro_elem
      dthydro_per_thread(i) = dthydro
    ENDDO
  
  
    DO i = 0, length-1
      indx = domain%m_matElemlist(i) ;
  
      IF (domain%m_vdov(indx) /= (0.0_RLK)) THEN
        dtdvov = dvovmax / (ABS(domain%m_vdov(indx))+(1.e-20_RLK))
  
  !     For sequential run (non OpenMP, replace line above with next line.
        thread_num = (0)
  
        IF ( dthydro_per_thread(thread_num) > dtdvov ) THEN
          dthydro_per_thread(thread_num) = dtdvov
          hydro_elem_per_thread(thread_num) = indx
        ENDIF
      ENDIF
    ENDDO
  
  
    DO i = 0, threads-1
      IF (dthydro_per_thread(i) < dthydro) THEN
        dthydro = dthydro_per_thread(i)
        hydro_elem =  hydro_elem_per_thread(i)
      ENDIF
    ENDDO
  
    IF (hydro_elem /= -1) THEN
      domain%m_dthydro = dthydro
    ENDIF
  
    RETURN
  
  END SUBROUTINE CalcHydroConstraintForElems

end module time_constraint
