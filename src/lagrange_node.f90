module lagrange_node

  use utils
  use main_data
  use calc_force

  contains

  SUBROUTINE LagrangeNodal()
  
    IMPLICIT NONE 
    REAL(KIND=8) :: delt
    REAL(KIND=8) :: u_cut
  
    delt  = domain%m_deltatime
    u_cut = domain%m_u_cut
  
  ! time of boundary condition evaluation is beginning of step for force and
  ! acceleration boundary conditions.
    CALL CalcForceForNodes()
  
    CALL CalcAccelerationForNodes()
  
    CALL ApplyAccelerationBoundaryConditionsForNodes()
  
    CALL CalcVelocityForNodes( delt, u_cut )
  
    CALL CalcPositionForNodes( delt )
  
  END SUBROUTINE LagrangeNodal
  
  SUBROUTINE CalcForceForNodes()
  
    IMPLICIT NONE 
    INTEGER(KIND=4) :: numNode
    INTEGER(KIND=4) :: i
  
    numNode = domain%m_numNode
    DO i=0, numNode-1
      domain%m_fx(i) = 0.0_RLK
      domain%m_fy(i) = 0.0_RLK
      domain%m_fz(i) = 0.0_RLK
    ENDDO
  
  ! Calcforce calls partial, force, hourq
    CALL CalcVolumeForceForElems()
  
  ! Calculate Nodal Forces at domain boundaries
  ! problem->commSBN->Transfer(CommSBN::forces)
  
  END SUBROUTINE CalcForceForNodes

  SUBROUTINE CalcAccelerationForNodes()
  
    IMPLICIT NONE 
    INTEGER(KIND=4) :: numNode
    INTEGER(KIND=4) :: i
  
    numNode = domain%m_numNode
    DO i=0, numNode-1
      domain%m_xdd(i) = domain%m_fx(i) / domain%m_nodalMass(i)
      domain%m_ydd(i) = domain%m_fy(i) / domain%m_nodalMass(i)
      domain%m_zdd(i) = domain%m_fz(i) / domain%m_nodalMass(i)
    ENDDO
  
  END SUBROUTINE CalcAccelerationForNodes
  
  SUBROUTINE ApplyAccelerationBoundaryConditionsForNodes()
  
    IMPLICIT NONE 
    INTEGER(KIND=4) :: numNodeBC
    INTEGER(KIND=4) :: i
  
    numNodeBC = (domain%m_sizeX+1)*(domain%m_sizeX+1)
  
    DO i=0, numNodeBC-1
      domain%m_xdd(domain%m_symmX(i)) = 0.0_RLK
    ENDDO
    DO i=0, numNodeBC-1
      domain%m_ydd(domain%m_symmY(i)) = 0.0_RLK
    ENDDO
    DO i=0, numNodeBC-1
      domain%m_zdd(domain%m_symmZ(i)) = 0.0_RLK
    ENDDO
  
  END SUBROUTINE ApplyAccelerationBoundaryConditionsForNodes
  
  SUBROUTINE CalcVelocityForNodes( dt, u_cut)
  
    IMPLICIT NONE 
    REAL(KIND=8)    :: dt, u_cut
    INTEGER(KIND=4) :: numNode
    INTEGER(KIND=4) :: i
    REAL(KIND=8)    :: xdtmp, ydtmp, zdtmp
  
    numNode = domain%m_numNode
  
    DO i = 0, numNode-1
  
      xdtmp = domain%m_xd(i) + domain%m_xdd(i) * dt
      if( ABS(xdtmp) < u_cut ) xdtmp = 0.0_RLK
      domain%m_xd(i) = xdtmp
  
      ydtmp = domain%m_yd(i) + domain%m_ydd(i) * dt
      if( ABS(ydtmp) < u_cut ) ydtmp = 0.0_RLK
      domain%m_yd(i) = ydtmp
  
      zdtmp = domain%m_zd(i) + domain%m_zdd(i) * dt
      if( ABS(zdtmp) < u_cut ) zdtmp = 0.0_RLK
      domain%m_zd(i) = zdtmp
    ENDDO
  
  
  END SUBROUTINE CalcVelocityForNodes
  
  SUBROUTINE CalcPositionForNodes(dt)
  
    IMPLICIT NONE 
    REAL(KIND=8)    :: dt
    INTEGER(KIND=4) :: numNode
    INTEGER(KIND=4) :: i
  
    numNode = domain%m_numNode
  
    DO i = 0, numNode-1
      domain%m_x(i) = domain%m_x(i) + domain%m_xd(i) * dt
      domain%m_y(i) = domain%m_y(i) + domain%m_yd(i) * dt
      domain%m_z(i) = domain%m_z(i) + domain%m_zd(i) * dt
    ENDDO
  
  END SUBROUTINE CalcPositionForNodes
  
end module lagrange_node
