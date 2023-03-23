module main_data

  INTEGER(KIND=4), PARAMETER :: VolumeError = -1
  INTEGER(KIND=4), PARAMETER :: QStopError  = -2

  INTEGER(KIND=4), PARAMETER :: RLK = 8

  TYPE domain_type
  
    !****************
    ! Implementation 
    !****************
    
    ! Node-centered 
    
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_x   ! coordinates 
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_y 
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_z 
    
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_xd  ! velocities 
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_yd 
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_zd 
    
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_xdd  ! accelerations 
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_ydd 
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_zdd 
    
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_fx   ! forces 
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_fy 
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_fz 
    
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_nodalMass   ! mass 
    
    INTEGER,     DIMENSION(:), ALLOCATABLE ::m_symmX   ! symmetry plane nodesets 
    INTEGER,     DIMENSION(:), ALLOCATABLE ::m_symmY 
    INTEGER,     DIMENSION(:), ALLOCATABLE ::m_symmZ 
    
    INTEGER,     DIMENSION(:), ALLOCATABLE ::m_nodeElemCount 
    INTEGER,     DIMENSION(:), ALLOCATABLE ::m_nodeElemStart 
    !   INTEGER,     DIMENSION(:), ALLOCATABLE ::m_nodeElemList 
    INTEGER,     DIMENSION(:), ALLOCATABLE ::m_nodeElemCornerList 
    
    ! Element-centered 
    
    INTEGER,     DIMENSION(:), ALLOCATABLE :: m_matElemlist   ! material indexset 
    INTEGER,     DIMENSION(:), POINTER     :: m_nodelist => NULL()  ! elemToNode connectivity 
    
    INTEGER,     DIMENSION(:), ALLOCATABLE :: m_lxim   ! element connectivity across each face 
    INTEGER,     DIMENSION(:), ALLOCATABLE :: m_lxip 
    INTEGER,     DIMENSION(:), ALLOCATABLE :: m_letam 
    INTEGER,     DIMENSION(:), ALLOCATABLE :: m_letap 
    INTEGER,     DIMENSION(:), ALLOCATABLE :: m_lzetam 
    INTEGER,     DIMENSION(:), ALLOCATABLE :: m_lzetap 
    
    INTEGER,     DIMENSION(:), ALLOCATABLE :: m_elemBC   ! symmetry/free-surface flags for each elem face 
    
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_dxx   ! principal strains -- temporary 
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_dyy 
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_dzz 
    
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_delv_xi     ! velocity gradient -- temporary 
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_delv_eta 
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_delv_zeta 
    
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_delx_xi     ! coordinate gradient -- temporary 
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_delx_eta 
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_delx_zeta 
    
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_e    ! energy 
    
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_p    ! pressure 
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_q    ! q 
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_ql   ! linear term for q 
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_qq   ! quadratic term for q 
    
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_v      ! relative volume 
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_volo   ! reference volume 
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_vnew   ! new relative volume -- temporary 
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_delv   ! m_vnew - m_v 
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_vdov   ! volume derivative over volume 
    
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_arealg   ! characteristic length of an element 
    
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_ss       ! "sound speed" 
    
    REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: m_elemMass   ! mass 
  
    ! Parameters 
    
    REAL(KIND=8)       ::   m_dtfixed            ! fixed time increment 
    REAL(KIND=8)       ::   m_time               ! current time 
    REAL(KIND=8)       ::   m_deltatime          ! variable time increment 
    REAL(KIND=8)       ::   m_deltatimemultlb 
    REAL(KIND=8)       ::   m_deltatimemultub 
    REAL(KIND=8)       ::   m_stoptime           ! end time for simulation 
    
    REAL(KIND=8)       ::   m_u_cut              ! velocity tolerance 
    REAL(KIND=8)       ::   m_hgcoef             ! hourglass control 
    REAL(KIND=8)       ::   m_qstop              ! excessive q indicator 
    REAL(KIND=8)       ::   m_monoq_max_slope 
    REAL(KIND=8)       ::   m_monoq_limiter_mult 
    REAL(KIND=8)       ::   m_e_cut              ! energy tolerance 
    REAL(KIND=8)       ::   m_p_cut              ! pressure tolerance 
    REAL(KIND=8)       ::   m_ss4o3 
    REAL(KIND=8)       ::   m_q_cut              ! q tolerance 
    REAL(KIND=8)       ::   m_v_cut              ! relative volume tolerance 
    REAL(KIND=8)       ::   m_qlc_monoq          ! linear term coef for q 
    REAL(KIND=8)       ::   m_qqc_monoq          ! quadratic term coef for q 
    REAL(KIND=8)       ::   m_qqc 
    REAL(KIND=8)       ::   m_eosvmax 
    REAL(KIND=8)       ::   m_eosvmin 
    REAL(KIND=8)       ::   m_pmin               ! pressure floor 
    REAL(KIND=8)       ::   m_emin               ! energy floor 
    REAL(KIND=8)       ::   m_dvovmax            ! maximum allowable volume change 
    REAL(KIND=8)       ::   m_refdens            ! reference density 
    
    REAL(KIND=8)       ::   m_dtcourant          ! courant constraint 
    REAL(KIND=8)       ::   m_dthydro            ! volume change constraint 
    REAL(KIND=8)       ::   m_dtmax              ! maximum allowable time increment 
    
    INTEGER ::   m_cycle              ! iteration count for simulation 
    
    INTEGER ::   m_sizeX            ! X,Y,Z extent of this block 
    INTEGER ::   m_sizeY 
    INTEGER ::   m_sizeZ 
    
    INTEGER ::   m_numElem          ! Elements/Nodes in this domain 
    INTEGER ::   m_numNode 
  
    
  END TYPE domain_type

  TYPE(domain_type) :: domain

  ! Stuff needed for boundary conditions
  ! 2 BCs on each of 6 hexahedral faces (12 bits)
  INTEGER, PARAMETER :: XI_M        = z'003' ! 0x003
  INTEGER, PARAMETER :: XI_M_SYMM   = z'001' ! 0x001
  INTEGER, PARAMETER :: XI_M_FREE   = z'002' ! 0x002
  
  INTEGER, PARAMETER :: XI_P        = z'00c' ! 0x00c
  INTEGER, PARAMETER :: XI_P_SYMM   = z'004' ! 0x004
  INTEGER, PARAMETER :: XI_P_FREE   = z'008' ! 0x008
  
  INTEGER, PARAMETER :: ETA_M       = z'030' ! 0x030
  INTEGER, PARAMETER :: ETA_M_SYMM  = z'010' ! 0x010
  INTEGER, PARAMETER :: ETA_M_FREE  = z'020' ! 0x020
  
  INTEGER, PARAMETER :: ETA_P       = z'0c0' ! 0x0c0
  INTEGER, PARAMETER :: ETA_P_SYMM  = z'040' ! 0x040
  INTEGER, PARAMETER :: ETA_P_FREE  = z'080' ! 0x080
  
  INTEGER, PARAMETER :: ZETA_M      = z'300' ! 0x300
  INTEGER, PARAMETER :: ZETA_M_SYMM = z'100' ! 0x100
  INTEGER, PARAMETER :: ZETA_M_FREE = z'200' ! 0x200
  
  INTEGER, PARAMETER :: ZETA_P      = z'c00' ! 0xc00
  INTEGER, PARAMETER :: ZETA_P_SYMM = z'400' ! 0x400
  INTEGER, PARAMETER :: ZETA_P_FREE = z'800' ! 0x800

end module main_data
