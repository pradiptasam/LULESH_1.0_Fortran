module lagrange_elem

  use main_data
  use utils
  contains

  SUBROUTINE LagrangeElements()
  
    IMPLICIT NONE 
    REAL(KIND=8) :: deltatime
  
    deltatime = domain%m_deltatime
  
    CALL CalcLagrangeElements(deltatime)
  
  ! Calculate Q.  (Monotonic q option requires communication)
    CALL CalcQForElems()
  
    CALL ApplyMaterialPropertiesForElems()
  
    CALL UpdateVolumesForElems()
  
  END SUBROUTINE LagrangeElements
  
  SUBROUTINE CalcLagrangeElements( deltatime)
  
    IMPLICIT NONE 
    REAL(KIND=8)    :: deltatime
    REAL(KIND=8)    :: vdov, vdovthird
    INTEGER(KIND=4) :: numElem, k
  
    numElem = domain%m_numElem
    IF (numElem > 0) THEN
      CALL CalcKinematicsForElems(numElem, deltatime)
  
  !   element loop to do some stuff not included in the elemlib function.
  
      DO k=0, numElem-1
  !     calc strain rate and apply as constraint (only done in FB element)
        vdov = domain%m_dxx(k) + domain%m_dyy(k) + domain%m_dzz(k)
        vdovthird = vdov/(3.0_RLK)
  
  !     make the rate of deformation tensor deviatoric
        domain%m_vdov(k) = vdov
        domain%m_dxx(k) = domain%m_dxx(k) - vdovthird
        domain%m_dyy(k) = domain%m_dyy(k) - vdovthird
        domain%m_dzz(k) = domain%m_dzz(k) - vdovthird
  
  !     See if any volumes are negative, and take appropriate action.
        IF (domain%m_vnew(k) <= (0.0_RLK)) THEN
          call luabort(VolumeError)
        ENDIF
      ENDDO
    ENDIF
  
  END SUBROUTINE CalcLagrangeElements
  
  SUBROUTINE CalcQForElems()
  
    IMPLICIT NONE 
    REAL(KIND=8)    :: qstop
    INTEGER(KIND=4) :: numElem, idx, i
  
    qstop = domain%m_qstop
    numElem = domain%m_numElem
  !
  ! MONOTONIC Q option
  !
  
  ! Calculate velocity gradients
    CALL CalcMonotonicQGradientsForElems()
  
  ! Transfer veloctiy gradients in the first order elements
  ! problem->commElements->Transfer(CommElements::monoQ)
    CALL CalcMonotonicQForElems()
  
  ! Don't allow excessive artificial viscosity
    IF (numElem /= 0) THEN
      idx = -1
      DO i = 0, numElem-1
        IF ( domain%m_q(i) > qstop ) THEN
          idx = i
          EXIT
        ENDIF
      ENDDO
  
      IF (idx >= 0) THEN
        CALL luabort(QStopError)
      ENDIF
    ENDIF
  
  END SUBROUTINE CalcQForElems
  
  SUBROUTINE ApplyMaterialPropertiesForElems()
  
    IMPLICIT NONE 
    REAL(KIND=8)    :: eosvmin, eosvmax, vc
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: vnewc
    INTEGER(KIND=4) :: length, zn
    integer :: i
  
    length = domain%m_numElem
  
    IF (length /= 0) THEN
  !   Expose all of the variables needed for material evaluation
      eosvmin = domain%m_eosvmin
      eosvmax = domain%m_eosvmax
      ALLOCATE(vnewc(0:length-1))
  
      DO i = 0, length-1
        zn = domain%m_matElemlist(i)
        vnewc(i) = domain%m_vnew(zn)
      ENDDO
  
      IF (eosvmin /= (0.0_RLK)) THEN
        DO i = 0, length-1
          IF (vnewc(i) < eosvmin) vnewc(i) = eosvmin
        ENDDO
      ENDIF
  
      IF (eosvmax /= (0.0_RLK)) THEN
        DO i = 0, length-1
          IF (vnewc(i) > eosvmax) vnewc(i) = eosvmax
        ENDDO
      ENDIF
  
      DO i = 0, length-1
        zn = domain%m_matElemlist(i)
        vc = domain%m_v(zn)
        IF (eosvmin /= (0.0_RLK)) THEN
          IF (vc < eosvmin) vc = eosvmin
        ENDIF
        IF (eosvmax /= (0.0_RLK)) THEN
          IF (vc > eosvmax) vc = eosvmax
        ENDIF
        IF (vc <= 0.0_RLK) THEN
          CALL luabort(VolumeError)
        ENDIF
      ENDDO
  
      CALL EvalEOSForElems(vnewc, length)
  
      DEALLOCATE(vnewc)
  
    ENDIF
  
  END SUBROUTINE ApplyMaterialPropertiesForElems
  
  SUBROUTINE UpdateVolumesForElems()
  
    IMPLICIT NONE 
    ReAL(KIND=8)    :: v_cut, tmpV
    INTEGER(KIND=4) :: numElem, i
  
    numElem = domain%m_numElem
  
    IF (numElem /= 0) THEN
      v_cut = domain%m_v_cut
  
      DO i = 0, numElem - 1
        tmpV = domain%m_vnew(i)
  
        IF ( ABS(tmpV - (1.0_RLK)) < v_cut ) THEN
          tmpV = (1.0_RLK)
        ENDIF
        domain%m_v(i) = tmpV
      ENDDO
    ENDIF
  
  END SUBROUTINE UpdateVolumesForElems
  
  SUBROUTINE CalcKinematicsForElems( numElem, dt )
  
    IMPLICIT NONE
    INTEGER      :: numElem
    INTEGER      :: k, lnode, gnode, j
    REAL(KIND=8) :: dt
  
    REAL(KIND=8), DIMENSION(0:7,0:2) :: B  ! shape function derivatives
    REAL(KIND=8), DIMENSION(0:5):: D
    REAL(KIND=8), DIMENSION(0:7) :: x_local
    REAL(KIND=8), DIMENSION(0:7) :: y_local
    REAL(KIND=8), DIMENSION(0:7) :: z_local
    REAL(KIND=8), DIMENSION(0:7) :: xd_local
    REAL(KIND=8), DIMENSION(0:7) :: yd_local
    REAL(KIND=8), DIMENSION(0:7) :: zd_local
    REAL(KIND=8) :: detJ, volume, relativeVolume,dt2
    INTEGER(KIND=4), DIMENSION(:), POINTER :: elemToNode => NULL()
  
    detJ = 0.0_RLK
  
  ! loop over all elements
    DO k = 0, numElem-1
      elemToNode => domain%m_nodelist(k*8:)
  
  !   get nodal coordinates from global arrays and copy into local arrays
      DO lnode=0, 7
        gnode = elemToNode(lnode+1)
        x_local(lnode) = domain%m_x(gnode)
        y_local(lnode) = domain%m_y(gnode)
        z_local(lnode) = domain%m_z(gnode)
      ENDDO
  
  !   volume calculations
      volume = CalcElemVolume(x_local, y_local, z_local )
      relativeVolume = volume / domain%m_volo(k)
      domain%m_vnew(k) = relativeVolume
      domain%m_delv(k) = relativeVolume - domain%m_v(k)
  
  !   set characteristic length
      domain%m_arealg(k) = CalcElemCharacteristicLength(x_local, y_local,  &
                                                        z_local, volume)
  
  !   get nodal velocities from global array and copy into local arrays.
      DO lnode=0, 7
        gnode = elemToNode(lnode+1);
        xd_local(lnode) = domain%m_xd(gnode)
        yd_local(lnode) = domain%m_yd(gnode)
        zd_local(lnode) = domain%m_zd(gnode)
      ENDDO
      
      dt2 = (0.5_RLK) * dt
      DO j=0, 7
        x_local(j) = x_local(j) - dt2 * xd_local(j)
        y_local(j) = y_local(j) - dt2 * yd_local(j)
        z_local(j) = z_local(j) - dt2 * zd_local(j)
      ENDDO
      
      CALL CalcElemShapeFunctionDerivatives( x_local, y_local, z_local,  &
                                             B, detJ )
      
      CALL CalcElemVelocityGrandient( xd_local, yd_local, zd_local,  &
                                      B, detJ, D )
      
  !   put velocity gradient quantities into their global arrays.
      domain%m_dxx(k) = D(0);
      domain%m_dyy(k) = D(1);
      domain%m_dzz(k) = D(2);
    ENDDO
  
  END SUBROUTINE CalcKinematicsForElems
  
  SUBROUTINE  CalcMonotonicQGradientsForElems()
  
    IMPLICIT NONE 
    REAL(KIND=8), PARAMETER :: ptiny = 1.e-36_RLK
    REAL(KIND=8)            :: ax,ay,az,dxv,dyv,dzv
    REAL(KIND=8)            :: x0,x1,x2,x3,x4,x5,x6,x7
    REAL(KIND=8)            :: y0,y1,y2,y3,y4,y5,y6,y7
    REAL(KIND=8)            :: z0,z1,z2,z3,z4,z5,z6,z7
    REAL(KIND=8)            :: xv0,xv1,xv2,xv3,xv4,xv5,xv6,xv7
    REAL(KIND=8)            :: yv0,yv1,yv2,yv3,yv4,yv5,yv6,yv7
    REAL(KIND=8)            :: zv0,zv1,zv2,zv3,zv4,zv5,zv6,zv7
    REAL(KIND=8)            :: vol, norm
    REAL(KIND=8)            :: dxi,dxj,dxk
    REAL(KIND=8)            :: dyi,dyj,dyk
    REAL(KIND=8)            :: dzi,dzj,dzk
    INTEGER(KIND=4)         :: numElem, i
    INTEGER(KIND=4)         :: n0,n1,n2,n3,n4,n5,n6,n7
    INTEGER(KIND=4), DIMENSION(:), POINTER :: elemToNode => NULL()
  
    numElem = domain%m_numElem
  
    DO i=0, numElem-1
      
      elemToNode => domain%m_nodelist(i*8:)
      n0 = elemToNode(1)
      n1 = elemToNode(2)
      n2 = elemToNode(3)
      n3 = elemToNode(4)
      n4 = elemToNode(5)
      n5 = elemToNode(6)
      n6 = elemToNode(7)
      n7 = elemToNode(8)
  
      x0 = domain%m_x(n0)
      x1 = domain%m_x(n1)
      x2 = domain%m_x(n2)
      x3 = domain%m_x(n3)
      x4 = domain%m_x(n4)
      x5 = domain%m_x(n5)
      x6 = domain%m_x(n6)
      x7 = domain%m_x(n7)
      
      y0 = domain%m_y(n0)
      y1 = domain%m_y(n1)
      y2 = domain%m_y(n2)
      y3 = domain%m_y(n3)
      y4 = domain%m_y(n4)
      y5 = domain%m_y(n5)
      y6 = domain%m_y(n6)
      y7 = domain%m_y(n7)
      
      z0 = domain%m_z(n0)
      z1 = domain%m_z(n1)
      z2 = domain%m_z(n2)
      z3 = domain%m_z(n3)
      z4 = domain%m_z(n4)
      z5 = domain%m_z(n5)
      z6 = domain%m_z(n6)
      z7 = domain%m_z(n7)
      
      xv0 = domain%m_xd(n0)
      xv1 = domain%m_xd(n1)
      xv2 = domain%m_xd(n2)
      xv3 = domain%m_xd(n3)
      xv4 = domain%m_xd(n4)
      xv5 = domain%m_xd(n5)
      xv6 = domain%m_xd(n6)
      xv7 = domain%m_xd(n7)
      
      yv0 = domain%m_yd(n0)
      yv1 = domain%m_yd(n1)
      yv2 = domain%m_yd(n2)
      yv3 = domain%m_yd(n3)
      yv4 = domain%m_yd(n4)
      yv5 = domain%m_yd(n5)
      yv6 = domain%m_yd(n6)
      yv7 = domain%m_yd(n7)
      
      zv0 = domain%m_zd(n0)
      zv1 = domain%m_zd(n1)
      zv2 = domain%m_zd(n2)
      zv3 = domain%m_zd(n3)
      zv4 = domain%m_zd(n4)
      zv5 = domain%m_zd(n5)
      zv6 = domain%m_zd(n6)
      zv7 = domain%m_zd(n7)
      
      vol = domain%m_volo(i)*domain%m_vnew(i)
      norm = (1.0_RLK) / ( vol + ptiny )
      
      dxj = (-0.25_RLK)*(SUM4(x0,x1,x5,x4) - SUM4(x3,x2,x6,x7))
      dyj = (-0.25_RLK)*(SUM4(y0,y1,y5,y4) - SUM4(y3,y2,y6,y7))
      dzj = (-0.25_RLK)*(SUM4(z0,z1,z5,z4) - SUM4(z3,z2,z6,z7))
      
      dxi = ( 0.25_RLK)*(SUM4(x1,x2,x6,x5) - SUM4(x0,x3,x7,x4))
      dyi = ( 0.25_RLK)*(SUM4(y1,y2,y6,y5) - SUM4(y0,y3,y7,y4))
      dzi = ( 0.25_RLK)*(SUM4(z1,z2,z6,z5) - SUM4(z0,z3,z7,z4))
      
      dxk = ( 0.25_RLK)*(SUM4(x4,x5,x6,x7) - SUM4(x0,x1,x2,x3))
      dyk = ( 0.25_RLK)*(SUM4(y4,y5,y6,y7) - SUM4(y0,y1,y2,y3))
      dzk = ( 0.25_RLK)*(SUM4(z4,z5,z6,z7) - SUM4(z0,z1,z2,z3))
      
  !   find delvk and delxk ( i cross j )
      
      ax = dyi*dzj - dzi*dyj
      ay = dzi*dxj - dxi*dzj
      az = dxi*dyj - dyi*dxj
      
      domain%m_delx_zeta(i) = vol / SQRT(ax*ax + ay*ay + az*az + ptiny)
      
      ax = ax * norm
      ay = ay * norm
      az = az * norm
      
      dxv = (0.25_RLK)*(SUM4(xv4,xv5,xv6,xv7) - SUM4(xv0,xv1,xv2,xv3))
      dyv = (0.25_RLK)*(SUM4(yv4,yv5,yv6,yv7) - SUM4(yv0,yv1,yv2,yv3))
      dzv = (0.25_RLK)*(SUM4(zv4,zv5,zv6,zv7) - SUM4(zv0,zv1,zv2,zv3))
      
      domain%m_delv_zeta(i) = ax*dxv + ay*dyv + az*dzv
      
  !   find delxi and delvi ( j cross k )
      
      ax = dyj*dzk - dzj*dyk ;
      ay = dzj*dxk - dxj*dzk ;
      az = dxj*dyk - dyj*dxk ;
      
      domain%m_delx_xi(i) = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;
      
      ax = ax * norm
      ay = ay * norm
      az = az * norm
      
      dxv = (0.25_RLK)*(SUM4(xv1,xv2,xv6,xv5) - SUM4(xv0,xv3,xv7,xv4)) ;
      dyv = (0.25_RLK)*(SUM4(yv1,yv2,yv6,yv5) - SUM4(yv0,yv3,yv7,yv4)) ;
      dzv = (0.25_RLK)*(SUM4(zv1,zv2,zv6,zv5) - SUM4(zv0,zv3,zv7,zv4)) ;
      
      domain%m_delv_xi(i) = ax*dxv + ay*dyv + az*dzv ;
      
  !   find delxj and delvj ( k cross i )
      
      ax = dyk*dzi - dzk*dyi ;
      ay = dzk*dxi - dxk*dzi ;
      az = dxk*dyi - dyk*dxi ;
      
      domain%m_delx_eta(i) = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;
      
      ax = ax * norm
      ay = ay * norm
      az = az * norm
      
      dxv = (-0.25_RLK)*(SUM4(xv0,xv1,xv5,xv4) - SUM4(xv3,xv2,xv6,xv7)) ;
      dyv = (-0.25_RLK)*(SUM4(yv0,yv1,yv5,yv4) - SUM4(yv3,yv2,yv6,yv7)) ;
      dzv = (-0.25_RLK)*(SUM4(zv0,zv1,zv5,zv4) - SUM4(zv3,zv2,zv6,zv7)) ;
      
      domain%m_delv_eta(i) = ax*dxv + ay*dyv + az*dzv ;
    ENDDO
  
  END SUBROUTINE CalcMonotonicQGradientsForElems
  
  SUBROUTINE CalcMonotonicQForElems()
  
    IMPLICIT NONE 
    REAL(KIND=8), PARAMETER :: ptiny = 1.e-36_RLK
    REAL(KIND=8) :: monoq_max_slope
    REAL(KIND=8) :: monoq_limiter_mult
    REAL(KIND=8) :: qlc_monoq
    REAL(KIND=8) :: qqc_monoq
    INTEGER(KIND=4) :: elength     ! the elementset length
  !
  ! initialize parameters
  !
    monoq_max_slope    = domain%m_monoq_max_slope
    monoq_limiter_mult = domain%m_monoq_limiter_mult
  
  !
  ! calculate the monotonic q for pure regions
  !
    elength = domain%m_numElem
    IF (elength > 0) THEN
      qlc_monoq = domain%m_qlc_monoq
      qqc_monoq = domain%m_qqc_monoq
      CALL CalcMonotonicQRegionForElems( qlc_monoq, qqc_monoq, &
                                         monoq_limiter_mult,   &
                                         monoq_max_slope,      &
                                         ptiny, elength )
    ENDIF
  
  END SUBROUTINE CalcMonotonicQForElems
  
  SUBROUTINE EvalEOSForElems(vnewc, length)
  
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(0:) :: vnewc
    INTEGER :: length
  
    REAL(KIND=8) :: e_cut, p_cut, ss4o3, q_cut
    REAL(KIND=8) :: eosvmax, eosvmin, pmin, emin, rho0
    REAL(KIND=8) :: vchalf
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: e_old, &
                  delvc, p_old, q_old, compression,   &
                  compHalfStep, qq, ql, work, p_new,  &
                  e_new, q_new, bvc, pbvc
    INTEGER      :: i, zidx
  
    e_cut = domain%m_e_cut
    p_cut = domain%m_p_cut
    ss4o3 = domain%m_ss4o3
    q_cut = domain%m_q_cut
  
    eosvmax = domain%m_eosvmax
    eosvmin = domain%m_eosvmin
    pmin    = domain%m_pmin
    emin    = domain%m_emin
    rho0    = domain%m_refdens
  
    ALLOCATE(e_old(0:length-1))
    ALLOCATE(delvc(0:length-1))
    ALLOCATE(p_old(0:length-1))
    ALLOCATE(q_old(0:length-1))
    ALLOCATE(compression(0:length-1))
    ALLOCATE(compHalfStep(0:length-1))
    ALLOCATE(qq(0:length-1))
    ALLOCATE(ql(0:length-1))
    ALLOCATE(work(0:length-1))
    ALLOCATE(p_new(0:length-1))
    ALLOCATE(e_new(0:length-1))
    ALLOCATE(q_new(0:length-1))
    ALLOCATE(bvc(0:length-1))
    ALLOCATE(pbvc(0:length-1))
  
  ! compress data, minimal set
    DO i = 0, length-1
      zidx = domain%m_matElemlist(i) ;
      e_old(i) = domain%m_e(zidx)
    ENDDO
    
    DO i = 0, length-1
      zidx = domain%m_matElemlist(i)
      delvc(i) = domain%m_delv(zidx)
    ENDDO
    
    DO i = 0, length-1
      zidx = domain%m_matElemlist(i)
      p_old(i) = domain%m_p(zidx)
    ENDDO
    
    DO i = 0, length-1
      zidx = domain%m_matElemlist(i)
      q_old(i) = domain%m_q(zidx)
    ENDDO
    
    DO i = 0, length-1
      compression(i) = (1.0_RLK) / vnewc(i) - (1.0_RLK)
      vchalf = vnewc(i) - delvc(i) * (0.5_RLK)
      compHalfStep(i) = (1.0_RLK) / vchalf - (1.0_RLK)
    ENDDO
    
  ! Check for v > eosvmax or v < eosvmin
    IF ( eosvmin /= (0.0_RLK) ) THEN
      DO i = 0, length-1
        IF (vnewc(i) <= eosvmin) THEN  ! impossible due to calling func?
          compHalfStep(i) = compression(i)
        ENDIF
      ENDDO
    ENDIF
    IF ( eosvmax /= (0.0_RLK) ) THEN
      DO i = 0, length-1
        IF (vnewc(i) >= eosvmax) THEN ! impossible due to calling func? 
          p_old(i)        = (0.0_RLK)
          compression(i)  = (0.0_RLK)
          compHalfStep(i) = (0.0_RLK)
        ENDIF
      ENDDO
    ENDIF
  
    DO i = 0, length-1
      zidx = domain%m_matElemlist(i)
      qq(i) = domain%m_qq(zidx)
      ql(i) = domain%m_ql(zidx)
      work(i) = (0.0_RLK)
    ENDDO
  
    CALL CalcEnergyForElems(p_new, e_new, q_new, bvc, pbvc,          &
                            p_old, e_old,  q_old, compression,       &
                            compHalfStep, vnewc, work,  delvc, pmin, &
                            p_cut, e_cut, q_cut, emin,               &
                            qq, ql, rho0, eosvmax, length)
  
  
    DO i = 0, length-1
      zidx = domain%m_matElemlist(i)
      domain%m_p(zidx) = p_new(i)
    ENDDO
  
    DO i = 0, length-1
      zidx = domain%m_matElemlist(i)
      domain%m_e(zidx) = e_new(i)
    ENDDO
  
    DO i = 0, length-1
      zidx = domain%m_matElemlist(i)
      domain%m_q(zidx) = q_new(i)
    ENDDO
  
  
    CALL CalcSoundSpeedForElems(vnewc, rho0, e_new, p_new,  &
                                pbvc, bvc, ss4o3, length)
  
    DEALLOCATE(pbvc)
    DEALLOCATE(bvc)
    DEALLOCATE(q_new)
    DEALLOCATE(e_new)
    DEALLOCATE(p_new)
    DEALLOCATE(work)
    DEALLOCATE(ql)
    DEALLOCATE(qq)
    DEALLOCATE(compHalfStep)
    DEALLOCATE(compression)
    DEALLOCATE(q_old)
    DEALLOCATE(p_old)
    DEALLOCATE(delvc)
    DEALLOCATE(e_old)
  
  
  END SUBROUTINE EvalEOSForElems

  SUBROUTINE CalcElemVelocityGrandient( xvel, yvel, zvel, &
                                        b, detJ, d )
  
    IMPLICIT NONE 
  
    REAL(KIND=8), DIMENSION(0:7),     INTENT(IN)  :: xvel, yvel, zvel
    REAL(KIND=8), DIMENSION(0:7,0:2), INTENT(IN)  :: b   ![3,8]
    REAL(KIND=8),                     INTENT(IN)  :: detJ
    REAL(KIND=8), DIMENSION(0:5),     INTENT(OUT) :: d
  
    REAL(KIND=8) :: dyddx, dxddy, dzddx, dxddz, dzddy, dyddz
    REAL(KIND=8) :: inv_detJ
    REAL(KIND=8), DIMENSION(0:7) :: pfx
    REAL(KIND=8), DIMENSION(0:7) :: pfy
    REAL(KIND=8), DIMENSION(0:7) :: pfz
  
    inv_detJ = (1.0_RLK) / detJ
    pfx = b(:,0)
    pfy = b(:,1)
    pfz = b(:,2)
  
    d(0) = inv_detJ * ( pfx(0) * (xvel(0)-xvel(6))   &
                      + pfx(1) * (xvel(1)-xvel(7))   &
                      + pfx(2) * (xvel(2)-xvel(4))   &
                      + pfx(3) * (xvel(3)-xvel(5)) )
    
    d(1) = inv_detJ * ( pfy(0) * (yvel(0)-yvel(6))   &
                      + pfy(1) * (yvel(1)-yvel(7))   &
                      + pfy(2) * (yvel(2)-yvel(4))   &
                      + pfy(3) * (yvel(3)-yvel(5)) )
    
    d(2) = inv_detJ * ( pfz(0) * (zvel(0)-zvel(6))   &
                      + pfz(1) * (zvel(1)-zvel(7))   &
                      + pfz(2) * (zvel(2)-zvel(4))   &
                      + pfz(3) * (zvel(3)-zvel(5)) )
  
    dyddx = inv_detJ * ( pfx(0) * (yvel(0)-yvel(6))  &
                       + pfx(1) * (yvel(1)-yvel(7))  &
                       + pfx(2) * (yvel(2)-yvel(4))  &
                       + pfx(3) * (yvel(3)-yvel(5)) )
    
    dxddy = inv_detJ * ( pfy(0) * (xvel(0)-xvel(6))  &
                       + pfy(1) * (xvel(1)-xvel(7))  &
                       + pfy(2) * (xvel(2)-xvel(4))  &
                       + pfy(3) * (xvel(3)-xvel(5)) )
    
    dzddx = inv_detJ * ( pfx(0) * (zvel(0)-zvel(6))  &
                       + pfx(1) * (zvel(1)-zvel(7))  &
                       + pfx(2) * (zvel(2)-zvel(4))  &
                       + pfx(3) * (zvel(3)-zvel(5)) )
    
    dxddz = inv_detJ * ( pfz(0) * (xvel(0)-xvel(6))  &
                       + pfz(1) * (xvel(1)-xvel(7))  &
                       + pfz(2) * (xvel(2)-xvel(4))  &
                       + pfz(3) * (xvel(3)-xvel(5)) )
    
    dzddy = inv_detJ * ( pfy(0) * (zvel(0)-zvel(6))  &
                       + pfy(1) * (zvel(1)-zvel(7))  &
                       + pfy(2) * (zvel(2)-zvel(4))  &
                       + pfy(3) * (zvel(3)-zvel(5)) )
    
    dyddz = inv_detJ * ( pfz(0) * (yvel(0)-yvel(6))  &
                       + pfz(1) * (yvel(1)-yvel(7))  &
                       + pfz(2) * (yvel(2)-yvel(4))  &
                       + pfz(3) * (yvel(3)-yvel(5)) )
  
    d(5) = (0.5_RLK) * ( dxddy + dyddx )
    d(4) = (0.5_RLK) * ( dxddz + dzddx )
    d(3) = (0.5_RLK) * ( dzddy + dyddz )
  
  END SUBROUTINE CalcElemVelocityGrandient
  
  SUBROUTINE CalcMonotonicQRegionForElems( qlc_monoq, qqc_monoq,                &
                                           monoq_limiter_mult, monoq_max_slope, &
                                           ptiny,                               &
                                           elength                              ) 
  
    IMPLICIT NONE 
    REAL(KIND=8) :: qlc_monoq,  qqc_monoq
    REAL(KIND=8) :: monoq_limiter_mult,  monoq_max_slope
    REAL(KIND=8) :: ptiny
    INTEGER(KIND=4) :: elength     ! the elementset length
    INTEGER(KIND=4) :: ielem, i, bcMask
    REAL(KIND=8) :: qlin, qquad, phixi, phieta, phizeta, delvm, delvp
    REAL(KIND=8) :: norm, delvxxi, delvxeta, delvxzeta, rho
  
  
    DO ielem = 0, elength-1
      i = domain%m_matElemlist(ielem)
      bcMask = domain%m_elemBC(i)
  
  !   phixi
      norm = (1.0_RLK) / ( domain%m_delv_xi(i) + ptiny )
  
      SELECT CASE(IAND(bcMask, XI_M))
        CASE (0)
          delvm = domain%m_delv_xi(domain%m_lxim(i))
        CASE (XI_M_SYMM)
          delvm = domain%m_delv_xi(i)
        CASE (XI_M_FREE)
          delvm = (0.0_RLK)
        CASE DEFAULT
  !       ERROR
      END SELECT
  
      SELECT CASE(IAND(bcMask, XI_P))
        CASE (0)
          delvp = domain%m_delv_xi(domain%m_lxip(i))
        CASE (XI_P_SYMM)
          delvp = domain%m_delv_xi(i)
        CASE (XI_P_FREE)
          delvp = (0.0_RLK)
        CASE DEFAULT
  !       ERROR 
      END SELECT
  
      delvm = delvm * norm
      delvp = delvp * norm
  
      phixi = (0.5_RLK) * ( delvm + delvp )
      
      delvm = delvm * monoq_limiter_mult
      delvp = delvp * monoq_limiter_mult
      
      if ( delvm < phixi ) phixi = delvm
      if ( delvp < phixi ) phixi = delvp
      if ( phixi < 0.0_RLK ) phixi = (0.0_RLK)
      if ( phixi > monoq_max_slope) phixi = monoq_max_slope
      
      
  !   phieta
      norm = (1.0_RLK) / ( domain%m_delv_eta(i) + ptiny )
      
      SELECT CASE(IAND(bcMask, ETA_M))
        CASE (0)
          delvm = domain%m_delv_eta(domain%m_letam(i))
        CASE (ETA_M_SYMM)
          delvm = domain%m_delv_eta(i)
        CASE (ETA_M_FREE)
          delvm = 0.0_RLK
        CASE DEFAULT
  !       ERROR
      END SELECT
      SELECT CASE(IAND(bcMask, ETA_P))
        CASE (0)
          delvp = domain%m_delv_eta(domain%m_letap(i))
        CASE (ETA_P_SYMM)
          delvp = domain%m_delv_eta(i)
        CASE (ETA_P_FREE)
          delvp = (0.0_RLK)
        CASE DEFAULT
  !       ERROR
      END SELECT
      
      delvm = delvm * norm
      delvp = delvp * norm
      
      phieta = (0.5_RLK) * ( delvm + delvp )
      
      delvm = delvm * monoq_limiter_mult
      delvp = delvp * monoq_limiter_mult
      
      if ( delvm  < phieta ) phieta = delvm
      if ( delvp  < phieta ) phieta = delvp
      if ( phieta < (0.0_RLK)) phieta = (0.0_RLK)
      if ( phieta > monoq_max_slope)  phieta = monoq_max_slope
      
  !   phizeta
      norm = (1.0_RLK) / ( domain%m_delv_zeta(i) + ptiny ) ;
      
      SELECT CASE(IAND(bcMask, ZETA_M))
        CASE (0)
          delvm = domain%m_delv_zeta(domain%m_lzetam(i))
        CASE (ZETA_M_SYMM)
          delvm = domain%m_delv_zeta(i)
        CASE (ZETA_M_FREE)
          delvm = (0.0_RLK)
        CASE DEFAULT
  !       ERROR
      END SELECT
      SELECT CASE(IAND(bcMask, ZETA_P))
        CASE (0)
          delvp = domain%m_delv_zeta(domain%m_lzetap(i))
        CASE (ZETA_P_SYMM)
          delvp = domain%m_delv_zeta(i)
        CASE (ZETA_P_FREE)
          delvp = (0.0_RLK)
        CASE DEFAULT
  !       ERROR
      END SELECT
  
      delvm = delvm * norm
      delvp = delvp * norm
  
      phizeta = (0.5_RLK) * ( delvm + delvp )
  
      delvm = delvm * monoq_limiter_mult
      delvp = delvp * monoq_limiter_mult
  
      IF ( delvm   < phizeta ) phizeta = delvm
      IF ( delvp   < phizeta ) phizeta = delvp
      IF ( phizeta < (0.0_RLK) ) phizeta = (0.0_RLK)
      IF ( phizeta > monoq_max_slope  ) phizeta = monoq_max_slope
  
  !   Remove length scale
  
      IF ( domain%m_vdov(i) > (0.0_RLK) ) THEN
        qlin  = (0.0_RLK)
        qquad = (0.0_RLK)
      ELSE
        delvxxi   = domain%m_delv_xi(i)   * domain%m_delx_xi(i)
        delvxeta  = domain%m_delv_eta(i)  * domain%m_delx_eta(i)
        delvxzeta = domain%m_delv_zeta(i) * domain%m_delx_zeta(i)
  
        IF ( delvxxi   > (0.0_RLK) ) delvxxi   = (0.0_RLK)
        IF ( delvxeta  > (0.0_RLK) ) delvxeta  = (0.0_RLK)
        IF ( delvxzeta > (0.0_RLK) ) delvxzeta = (0.0_RLK)
  
        rho = domain%m_elemMass(i) / (domain%m_volo(i) * domain%m_vnew(i))
  
        qlin = -qlc_monoq * rho *                      &
               (  delvxxi   * ((1.0_RLK) - phixi)  +     &
                  delvxeta  * ((1.0_RLK) - phieta) +     &
                  delvxzeta * ((1.0_RLK) - phizeta)  )
  
        qquad = qqc_monoq * rho *                                       &
               (  delvxxi*delvxxi     * ((1.0_RLK) - phixi*phixi)   +     &
                  delvxeta*delvxeta   * ((1.0_RLK) - phieta*phieta) +     &
                  delvxzeta*delvxzeta * ((1.0_RLK) - phizeta*phizeta)  )
      ENDIF
  
      domain%m_qq(i) = qquad
      domain%m_ql(i) = qlin
    ENDDO
  
  END SUBROUTINE CalcMonotonicQRegionForElems
  
  SUBROUTINE  CalcEnergyForElems( p_new,  e_new,  q_new,          &
                                  bvc,  pbvc,                     &
                                  p_old,  e_old,  q_old,          &
                                  compression,  compHalfStep,     &
                                  vnewc,  work,  delvc,  pmin,    &
                                  p_cut,   e_cut,  q_cut,  emin,  &
                                  qq,  ql,                        &
                                  rho0,                           &
                                  eosvmax,                        &
                                  length                          ) 
  
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(0:) :: p_new, e_new, q_new
    REAL(KIND=8), DIMENSION(0:) :: bvc,  pbvc
    REAL(KIND=8), DIMENSION(0:) :: p_old, e_old, q_old
    REAL(KIND=8), DIMENSION(0:) :: compression, compHalfStep
    REAL(KIND=8), DIMENSION(0:) :: vnewc, work, delvc
    REAL(KIND=8)    :: pmin, p_cut,  e_cut, q_cut, emin
    REAL(KIND=8), DIMENSION(0:) :: qq, ql
    REAL(KIND=8)    :: rho0
    REAL(KIND=8)    :: eosvmax
    INTEGER(KIND=4) :: length
  
    INTEGER(KIND=4) :: i
    REAL(KIND=8)    :: vhalf, ssc, q_tilde
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: pHalfStep
    REAL(KIND=8), PARAMETER :: TINY1 = 0.111111e-36_RLK
    REAL(KIND=8), PARAMETER :: TINY3 = 0.333333e-18_RLK
    REAL(KIND=8), PARAMETER :: SIXTH = (1.0_RLK) / (6.0_RLK)
  
  
    ALLOCATE(pHalfStep(0:length-1))
  
    DO i = 0, length-1
      e_new(i) = e_old(i) - (0.5_RLK) * delvc(i) * (p_old(i) + q_old(i))  &
               + (0.5_RLK) * work(i)
  
      IF (e_new(i)  < emin ) THEN
        e_new(i) = emin
      ENDIF
    ENDDO
  
    CALL CalcPressureForElems(pHalfStep, bvc, pbvc, e_new, compHalfStep,  &
                              vnewc, pmin, p_cut, eosvmax, length)
  
    DO i = 0, length-1
      vhalf = (1.0_RLK) / ((1.0_RLK) + compHalfStep(i))
      
      IF ( delvc(i) > (0.0_RLK) ) THEN
  !      q_new(i) /* = qq(i) = ql(i) */ = Real_t(0.) ;
        q_new(i) = (0.0_RLK)
      ELSE
        ssc = ( pbvc(i) * e_new(i)   &
            + vhalf * vhalf * bvc(i) * pHalfStep(i) ) / rho0
  
        IF ( ssc <= TINY1 ) THEN
          ssc = TINY3
        ELSE
          ssc = SQRT(ssc)
        ENDIF
  
        q_new(i) = (ssc*ql(i) + qq(i))
      ENDIF
  
      e_new(i) = e_new(i) + (0.5_RLK) * delvc(i)   &
         * (  (3.0_RLK)*(p_old(i)     + q_old(i))  &
            - (4.0_RLK)*(pHalfStep(i) + q_new(i)))
    ENDDO
    
    DO i = 0, length-1
      e_new(i) = e_new(i) + (0.5_RLK) * work(i)
  
      IF (ABS(e_new(i)) < e_cut) THEN
        e_new(i) = (0.0_RLK)
      ENDIF
      IF (e_new(i)  < emin ) THEN
        e_new(i) = emin
      ENDIF
    ENDDO
  
    CALL CalcPressureForElems(p_new, bvc, pbvc, e_new, compression,  &
                              vnewc, pmin, p_cut, eosvmax, length)
  
    DO i = 0, length-1
      IF (delvc(i) > (0.0_RLK)) THEN
        q_tilde = (0.0_RLK)
      ELSE
        ssc = ( pbvc(i) * e_new(i)     &
            + vnewc(i) * vnewc(i) * bvc(i) * p_new(i) ) / rho0
  
        IF ( ssc <= TINY1 ) THEN
          ssc = TINY3
        ELSE
          ssc = SQRT(ssc)
        ENDIF
  
        q_tilde = (ssc*ql(i) + qq(i))
      ENDIF
  
      e_new(i) = e_new(i) - (  (7.0_RLK)*(p_old(i)     + q_old(i))   &
                          -    (8.0_RLK)*(pHalfStep(i) + q_new(i))   &
                          + (p_new(i) + q_tilde)) * delvc(i)*SIXTH
  
      IF (ABS(e_new(i)) < e_cut) THEN
        e_new(i) = (0.0_RLK)
      ENDIF
      IF ( e_new(i)  < emin ) THEN
        e_new(i) = emin
      ENDIF
    ENDDO
  
    CALL CalcPressureForElems(p_new, bvc, pbvc, e_new, compression,  &
                              vnewc, pmin, p_cut, eosvmax, length)
  
    DO i = 0, length-1
  
      IF ( delvc(i) <= (0.0_RLK) ) THEN
        ssc = ( pbvc(i) * e_new(i)  &
            + vnewc(i) * vnewc(i) * bvc(i) * p_new(i) ) / rho0
  
        IF ( ssc <= TINY1 ) THEN
          ssc = TINY3
        ELSE
          ssc = SQRT(ssc)
        ENDIF
  
        q_new(i) = (ssc*ql(i) + qq(i))
  
        if (ABS(q_new(i)) < q_cut) q_new(i) = (0.0_RLK)
      ENDIF
    ENDDO
  
    DEALLOCATE(pHalfStep)
  
  END SUBROUTINE CalcEnergyForElems
  
  SUBROUTINE CalcSoundSpeedForElems(vnewc,  rho0, enewc, &
                                    pnewc, pbvc,         &
                                    bvc, ss4o3, nz       )
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(0:) :: vnewc, enewc
    REAL(KIND=8), DIMENSION(0:) :: pnewc, pbvc
    REAL(KIND=8), DIMENSION(0:) :: bvc
    REAL(KIND=8) :: rho0
    REAL(KIND=8) :: ss4o3
    INTEGER      :: nz
    REAL(KIND=8), PARAMETER :: TINY1 = 0.111111e-36_RLK
    REAL(KIND=8), PARAMETER :: TINY3 = 0.333333e-18_RLK
    REAL(KIND=8) :: ssTmp
    INTEGER      :: i, iz
  
    DO i = 0, nz - 1
      iz = domain%m_matElemlist(i)
      ssTmp = (pbvc(i) * enewc(i) + vnewc(i) * vnewc(i) *  &
                           bvc(i) * pnewc(i)) / rho0
      IF (ssTmp <= TINY1) THEN
        ssTmp = TINY3
      ELSE
        ssTmp = SQRT(ssTmp)
      ENDIF
      domain%m_ss(iz) = ssTmp
    ENDDO
  
  END SUBROUTINE CalcSoundSpeedForElems

  SUBROUTINE CalcPressureForElems( p_new, bvc,         &
                                   pbvc, e_old,        &
                                   compression, vnewc, &
                                   pmin,               &
                                   p_cut,eosvmax,      &
                                   length              )
  
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(0:) :: p_new, bvc, pbvc, e_old
    REAL(KIND=8), DIMENSION(0:) ::  compression
    REAL(KIND=8), DIMENSION(0:) :: vnewc
    REAL(KIND=8)    :: pmin
    REAL(KIND=8)    :: p_cut
    REAL(KIND=8)    :: eosvmax
    INTEGER(KIND=4) :: length 
  
    INTEGER(KIND=4) :: i
    REAL(KIND=8), PARAMETER :: c1s = (2.0_RLK)/(3.0_RLK)
  
    DO i = 0, length-1
      bvc(i) = c1s * (compression(i) + (1.0_RLK))
      pbvc(i) = c1s
    ENDDO
  
    DO i = 0, length-1
      p_new(i) = bvc(i) * e_old(i)
  
      IF (ABS(p_new(i)) < p_cut) THEN
        p_new(i) = (0.0_RLK)
      ENDIF
  
      IF ( vnewc(i) >= eosvmax ) THEN  ! impossible condition here?
        p_new(i) = (0.0_RLK)
      ENDIF
  
      IF (p_new(i) < pmin) THEN
        p_new(i) = pmin
      ENDIF
    ENDDO
  
  END SUBROUTINE CalcPressureForElems

end module lagrange_elem
