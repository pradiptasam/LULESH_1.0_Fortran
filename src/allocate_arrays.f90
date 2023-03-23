module allocate_arrays 

  use main_data
  use utils
  IMPLICIT NONE 

  contains 

  SUBROUTINE AllocateNodalPersistent(size)
  
    IMPLICIT NONE 
    INTEGER :: size
  
    ALLOCATE(domain%m_x(0:size-1))
    ALLOCATE(domain%m_y(0:size-1)) 
    ALLOCATE(domain%m_z(0:size-1)) 
    
    ALLOCATE(domain%m_xd(0:size-1)) 
    ALLOCATE(domain%m_yd(0:size-1)) 
    ALLOCATE(domain%m_zd(0:size-1)) 
    domain%m_xd = 0.0_RLK
    domain%m_yd = 0.0_RLK
    domain%m_zd = 0.0_RLK
  
    ALLOCATE(domain%m_xdd(0:size-1)) 
    ALLOCATE(domain%m_ydd(0:size-1)) 
    ALLOCATE(domain%m_zdd(0:size-1)) 
    domain%m_xdd = 0.0_RLK
    domain%m_ydd = 0.0_RLK
    domain%m_zdd = 0.0_RLK
    
    ALLOCATE(domain%m_fx(0:size-1)) 
    ALLOCATE(domain%m_fy(0:size-1)) 
    ALLOCATE(domain%m_fz(0:size-1)) 
    
    ALLOCATE(domain%m_nodalMass(0:size-1))
    domain%m_nodalMass = 0.0_RLK
  
  END SUBROUTINE AllocateNodalPersistent
  
  
  
  
  SUBROUTINE AllocateElemPersistent(size)
    IMPLICIT NONE
    INTEGER :: size 
    
    ALLOCATE(domain%m_matElemlist(0:size-1)) 
    ALLOCATE(domain%m_nodelist(0:8*size-1)) 
    
    ALLOCATE(domain%m_lxim(0:size-1)) 
    ALLOCATE(domain%m_lxip(0:size-1)) 
    ALLOCATE(domain%m_letam(0:size-1)) 
    ALLOCATE(domain%m_letap(0:size-1)) 
    ALLOCATE(domain%m_lzetam(0:size-1)) 
    ALLOCATE(domain%m_lzetap(0:size-1)) 
      
    ALLOCATE(domain%m_elemBC(0:size-1)) 
    
    ALLOCATE(domain%m_e(0:size-1)) 
    domain%m_e = 0.0_RLK
    ALLOCATE(domain%m_p(0:size-1)) 
    domain%m_p = 0.0_RLK
    ALLOCATE(domain%m_q(0:size-1)) 
    ALLOCATE(domain%m_ql(0:size-1)) 
    ALLOCATE(domain%m_qq(0:size-1)) 
    
    ALLOCATE(domain%m_v(0:size-1)) 
    domain%m_v = 1.0_RLK
    ALLOCATE(domain%m_volo(0:size-1)) 
    ALLOCATE(domain%m_delv(0:size-1)) 
    ALLOCATE(domain%m_vdov(0:size-1)) 
    
    ALLOCATE(domain%m_arealg(0:size-1)) 
    
    ALLOCATE(domain%m_ss(0:size-1)) 
    
    ALLOCATE(domain%m_elemMass(0:size-1)) 
    
  END SUBROUTINE AllocateElemPersistent
  
  
  
  
  
  !Temporaries should not be initialized in bulk but
  !this is a runnable placeholder for now
  SUBROUTINE AllocateElemTemporary( size)
    IMPLICIT NONE 
    INTEGER :: size
    
    ALLOCATE(domain%m_dxx(0:size-1))
    ALLOCATE(domain%m_dyy(0:size-1))
    ALLOCATE(domain%m_dzz(0:size-1))
    
    ALLOCATE(domain%m_delv_xi(0:size-1))
    ALLOCATE(domain%m_delv_eta(0:size-1))
    ALLOCATE(domain%m_delv_zeta(0:size-1))
    
    ALLOCATE(domain%m_delx_xi(0:size-1))
    ALLOCATE(domain%m_delx_eta(0:size-1))
    ALLOCATE(domain%m_delx_zeta(0:size-1))
    
    ALLOCATE(domain%m_vnew(0:size-1))
    
  END SUBROUTINE AllocateElemTemporary
  
  
  SUBROUTINE AllocateNodesets( size)
    IMPLICIT NONE 
    INTEGER :: size
  
    ALLOCATE(domain%m_symmX(0:size-1))
    ALLOCATE(domain%m_symmY(0:size-1))
    ALLOCATE(domain%m_symmZ(0:size-1))
  
  END SUBROUTINE AllocateNodesets
  
  
  
  
  
  SUBROUTINE AllocateNodeElemIndexes()
    IMPLICIT NONE 
  
    INTEGER :: m
    INTEGER :: i,j,k
    INTEGER :: offset
    INTEGER :: clSize, clv
    INTEGER :: numElem 
    INTEGER :: numNode
    INTEGER :: nodelist_len
  
    numElem = domain%m_numElem
    numNode = domain%m_numNode
  
    ! set up node-centered indexing of elements
    ALLOCATE(domain%m_nodeElemCount(0:numNode-1))
    !  m_nodeElemCount.resize(numNode);
  
    domain%m_nodeElemCount=0
  
    DO i=0, SIZE(domain%m_nodelist)-1  
       domain%m_nodeElemCount(domain%m_nodelist(i))=   &
            domain%m_nodeElemCount(domain%m_nodelist(i))+1
    END DO
  
    ALLOCATE(domain%m_nodeElemStart(0:numNode-1))  
  
    domain%m_nodeElemStart=0
  
    DO i=1,numNode-1
       domain%m_nodeElemStart(i) = domain%m_nodeElemStart(i-1) + domain%m_nodeElemCount(i-1)
    END DO
  
    ALLOCATE(domain%m_nodeElemCornerList(0:(domain%m_nodeElemStart(numNode-1) +  &
                                            domain%m_nodeElemCount(numNode-1))))
  
    domain%m_nodeElemCount=0  
  
  
    !CHECK THIS LOOP - THINK IT'S OK ????
  
    DO i=0, SIZE(domain%m_nodelist)-1
          m=domain%m_nodelist(i)
          offset = domain%m_nodeElemStart(m)+domain%m_nodeElemCount(m)
          domain%m_nodeElemCornerList(offset) = i
          domain%m_nodeElemCount(m) = domain%m_nodeElemCount(m) + 1
    END DO
  
    clSize = SIZE(domain%m_nodeElemCornerList)
    DO i=0, clSize-1
       clv=domain%m_nodeElemCornerList(i)
       IF ((clv.LT.0).OR.(clv.GT.numElem*8))THEN
          PRINT*,"AllocateNodeElemIndexes(): nodeElemCornerList entry out of range!"
          CALL luabort(1)
       END IF
    END DO
  
  
  END SUBROUTINE AllocateNodeElemIndexes

end module allocate_arrays 
