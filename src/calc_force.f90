module calc_force

  use utils
  use main_data

  contains

  SUBROUTINE CalcVolumeForceForElems()
  
    IMPLICIT NONE
    INTEGER(KIND=4) :: numElem
    INTEGER(KIND=4) :: k
    REAL(KIND=8) :: hgcoef
    REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: sigxx
    REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: sigyy
    REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: sigzz
    REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: determ
  
  
    numElem = domain%m_numElem
    IF (numElem /= 0) THEN
      hgcoef = domain%m_hgcoef
      ALLOCATE(sigxx(0:numElem-1))
      ALLOCATE(sigyy(0:numElem-1))
      ALLOCATE(sigzz(0:numElem-1))
      ALLOCATE(determ(0:numElem-1))
  
  !   Sum contributions to total stress tensor
      CALL InitStressTermsForElems(numElem, sigxx, sigyy, sigzz)
  
  !   call elemlib stress integration loop to produce nodal forces from
  !   material stresses.
      CALL IntegrateStressForElems( numElem, sigxx, sigyy, sigzz, determ)
  
  !   check for negative element volume and abort if found
      DO k=0, numElem-1
         IF (determ(k) <= 0.0_RLK) THEN
           CALL luabort(VolumeError)
         ENDIF
      ENDDO
  
      CALL CalcHourglassControlForElems(determ, hgcoef)
  
      DEALLOCATE(determ)
      DEALLOCATE(sigzz)
      DEALLOCATE(sigyy)
      DEALLOCATE(sigxx)
    ENDIF
  
  END SUBROUTINE CalcVolumeForceForElems
  
  SUBROUTINE InitStressTermsForElems( numElem,  sigxx, sigyy, sigzz)
  
    IMPLICIT NONE
  
    INTEGER         :: numElem
    REAL(KIND=8), DIMENSION(0:) :: sigxx
    REAL(KIND=8), DIMENSION(0:) :: sigyy
    REAL(KIND=8), DIMENSION(0:) :: sigzz
    INTEGER(KIND=4) :: ii
  
    DO ii = 0, numElem-1
      sigxx(ii) =  - domain%m_p(ii) - domain%m_q(ii)
      sigyy(ii) =  - domain%m_p(ii) - domain%m_q(ii)
      sigzz(ii) =  - domain%m_p(ii) - domain%m_q(ii)
    ENDDO
  
  END SUBROUTINE InitStressTermsForElems

  SUBROUTINE IntegrateStressForElems(numElem, sigxx, sigyy, sigzz, determ)
    IMPLICIT NONE
  
    INTEGER      :: numElem
    REAL(KIND=8),DIMENSION(0:) :: sigxx, sigyy, sigzz
    REAL(KIND=8),DIMENSION(0:), INTENT(INOUT) :: determ  
  
    REAL(KIND=8) :: fx, fy, fz
    REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: fx_elem, fy_elem, fz_elem
    REAL(KIND=8),DIMENSION(0:7,0:2) :: B   ! shape function derivatives
    REAL(KIND=8),DIMENSION(0:7)   :: x_local
    REAL(KIND=8),DIMENSION(0:7)   :: y_local
    REAL(KIND=8),DIMENSION(0:7)   :: z_local
    INTEGER(KIND=4), DIMENSION(:), POINTER :: elemNodes => NULL()
    INTEGER      :: lnode, gnode, count, start, elem, kk, i
    INTEGER      :: numNode, numElem8
  
    numElem8 = numElem * 8
    ALLOCATE(fx_elem(0:numElem8-1))
    ALLOCATE(fy_elem(0:numElem8-1))
    ALLOCATE(fz_elem(0:numElem8-1))
    
  ! loop over all elements
    DO kk=0, numElem-1
      elemNodes => domain%m_nodelist(kk*8:)
  
  !   get nodal coordinates from global arrays and copy into local arrays.
      DO lnode=0, 7
        gnode = elemNodes(lnode+1)
        x_local(lnode) = domain%m_x(gnode)
        y_local(lnode) = domain%m_y(gnode)
        z_local(lnode) = domain%m_z(gnode)
      ENDDO
  
  !   Volume calculation involves extra work for numerical consistency.
      CALL CalcElemShapeFunctionDerivatives(x_local, y_local, z_local, &
                                            B, determ(kk))
  
      CALL CalcElemNodeNormals( B(:,0) , B(:,1), B(:,2), x_local, y_local, z_local )
  
      CALL SumElemStressesToNodeForces( B, sigxx(kk), sigyy(kk), sigzz(kk),  &
                                        fx_elem(kk*8), fy_elem(kk*8), fz_elem(kk*8) )
  
! #if 0
! !   copy nodal force contributions to global force arrray.
!     DO lnode=0, 7
!       node = elemNodes(lnode+1)
!       domain%m_fx(gnode) = domain%m_fx(gnode) + fx_local(lnode)
!       domain%m_fy(gnode) = domain%m_fy(gnode) + fy_local(lnode)
!       domain%m_fz(gnode) = domain%m_fz(gnode) + fz_local(lnode)
!     ENDDO
! #endif
    ENDDO
  
    numNode = domain%m_numNode
  
    DO gnode=0, numNode-1
      count = domain%m_nodeElemCount(gnode)
      start = domain%m_nodeElemStart(gnode)
      fx = (0.0_RLK)
      fy = (0.0_RLK)
      fz = (0.0_RLK)
      DO i=0, count-1
        elem = domain%m_nodeElemCornerList(start+i)
        fx = fx + fx_elem(elem)
        fy = fy + fy_elem(elem)
        fz = fz + fz_elem(elem)
      ENDDO
      domain%m_fx(gnode) = fx
      domain%m_fy(gnode) = fy
      domain%m_fz(gnode) = fz
    ENDDO
  
    DEALLOCATE(fz_elem)
    DEALLOCATE(fy_elem)
    DEALLOCATE(fx_elem)
  
  END SUBROUTINE IntegrateStressForElems

  SUBROUTINE CalcHourglassControlForElems(determ, hgcoef)
  
    IMPLICIT NONE
    REAL(KIND=8),DIMENSION(0:) :: determ
    REAL(KIND=8) :: hgcoef
  
    REAL(KIND=8),DIMENSION(:), ALLOCATABLE :: dvdx
    REAL(KIND=8),DIMENSION(:), ALLOCATABLE :: dvdy
    REAL(KIND=8),DIMENSION(:), ALLOCATABLE :: dvdz
    REAL(KIND=8),DIMENSION(:), ALLOCATABLE :: x8n
    REAL(KIND=8),DIMENSION(:), ALLOCATABLE :: y8n
    REAL(KIND=8),DIMENSION(:), ALLOCATABLE :: z8n
    REAL(KIND=8),DIMENSION(0:7) :: x1, y1, z1
    REAL(KIND=8),DIMENSION(0:7) :: pfx, pfy, pfz
    INTEGER(KIND=4) :: numElem, numElem8, i, ii, jj
    INTEGER(KIND=4), DIMENSION(:), POINTER :: elemToNode => NULL()
  
    numElem = domain%m_numElem
    numElem8 = numElem * 8
    ALLOCATE(dvdx(0:numElem8-1))
    ALLOCATE(dvdy(0:numElem8-1))
    ALLOCATE(dvdz(0:numElem8-1))
    ALLOCATE(x8n(0:numElem8-1))
    ALLOCATE(y8n(0:numElem8-1))
    ALLOCATE(z8n(0:numElem8-1))
  
  ! start loop over elements
    DO i=0, numElem-1
      elemToNode => domain%m_nodelist(i*8:)
      CALL CollectDomainNodesToElemNodes(elemToNode, x1, y1, z1)
      CALL CalcElemVolumeDerivative(pfx, pfy, pfz, x1, y1, z1)
  
  !   load into temporary storage for FB Hour Glass control
      DO ii=0, 7
        jj=8*i+ii
  
        dvdx(jj) = pfx(ii)
        dvdy(jj) = pfy(ii)
        dvdz(jj) = pfz(ii)
        x8n(jj)  = x1(ii)
        y8n(jj)  = y1(ii)
        z8n(jj)  = z1(ii)
      ENDDO
  
      determ(i) = domain%m_volo(i) * domain%m_v(i)
  
  !   Do a check for negative volumes
      IF ( domain%m_v(i) <= (0.0_RLK) ) THEN
        CALL luabort(VolumeError)
      ENDIF
    ENDDO
  
    IF ( hgcoef > (0.0_RLK) ) THEN
      CALL CalcFBHourglassForceForElems(determ,x8n,y8n,z8n,dvdx,dvdy,dvdz,hgcoef)
    ENDIF
  
    DEALLOCATE(z8n)
    DEALLOCATE(y8n)
    DEALLOCATE(x8n)
    DEALLOCATE(dvdz)
    DEALLOCATE(dvdy)
    DEALLOCATE(dvdx)
  
    RETURN
  
  END SUBROUTINE CalcHourglassControlForElems
  
  SUBROUTINE CollectDomainNodesToElemNodes(elemToNode, elemX, elemY, elemZ)
  
    IMPLICIT NONE 
  
    INTEGER, DIMENSION(:), POINTER :: elemToNode
    REAL(KIND=8),DIMENSION(0:7)    :: elemX, elemY, elemZ
  
    INTEGER(KIND=4) :: nd0i, nd1i, nd2i, nd3i
    INTEGER(KIND=4) :: nd4i, nd5i, nd6i, nd7i
  
    nd0i = elemToNode(1)
    nd1i = elemToNode(2)
    nd2i = elemToNode(3)
    nd3i = elemToNode(4)
    nd4i = elemToNode(5)
    nd5i = elemToNode(6)
    nd6i = elemToNode(7)
    nd7i = elemToNode(8)
  
    elemX(0) = domain%m_x(nd0i)
    elemX(1) = domain%m_x(nd1i)
    elemX(2) = domain%m_x(nd2i)
    elemX(3) = domain%m_x(nd3i)
    elemX(4) = domain%m_x(nd4i)
    elemX(5) = domain%m_x(nd5i)
    elemX(6) = domain%m_x(nd6i)
    elemX(7) = domain%m_x(nd7i)
  
    elemY(0) = domain%m_y(nd0i)
    elemY(1) = domain%m_y(nd1i)
    elemY(2) = domain%m_y(nd2i)
    elemY(3) = domain%m_y(nd3i)
    elemY(4) = domain%m_y(nd4i)
    elemY(5) = domain%m_y(nd5i)
    elemY(6) = domain%m_y(nd6i)
    elemY(7) = domain%m_y(nd7i)
  
    elemZ(0) = domain%m_z(nd0i)
    elemZ(1) = domain%m_z(nd1i)
    elemZ(2) = domain%m_z(nd2i)
    elemZ(3) = domain%m_z(nd3i)
    elemZ(4) = domain%m_z(nd4i)
    elemZ(5) = domain%m_z(nd5i)
    elemZ(6) = domain%m_z(nd6i)
    elemZ(7) = domain%m_z(nd7i)
  
  END SUBROUTINE CollectDomainNodesToElemNodes
  
  SUBROUTINE CalcElemVolumeDerivative(dvdx,dvdy,dvdz, x, y, z)
  
  
    IMPLICIT NONE
    REAL(KIND=8),DIMENSION(0:7) :: dvdx, dvdy, dvdz
    REAL(KIND=8),DIMENSION(0:7) :: x, y, z
  
    CALL VoluDer(x(1), x(2), x(3), x(4), x(5), x(7),  &
                 y(1), y(2), y(3), y(4), y(5), y(7),  &
                 z(1), z(2), z(3), z(4), z(5), z(7),  &
                 dvdx(0), dvdy(0), dvdz(0))
    CALL VoluDer(x(0), x(1), x(2), x(7), x(4), x(6),  &
                 y(0), y(1), y(2), y(7), y(4), y(6),  &
                 z(0), z(1), z(2), z(7), z(4), z(6),  &
                 dvdx(3), dvdy(3), dvdz(3))
    CALL VoluDer(x(3), x(0), x(1), x(6), x(7), x(5),  &
                 y(3), y(0), y(1), y(6), y(7), y(5),  &
                 z(3), z(0), z(1), z(6), z(7), z(5),  &
                 dvdx(2), dvdy(2), dvdz(2))
    CALL VoluDer(x(2), x(3), x(0), x(5), x(6), x(4),  &
                 y(2), y(3), y(0), y(5), y(6), y(4),  &
                 z(2), z(3), z(0), z(5), z(6), z(4),  &
                 dvdx(1), dvdy(1), dvdz(1))
    CALL VoluDer(x(7), x(6), x(5), x(0), x(3), x(1),  &
                 y(7), y(6), y(5), y(0), y(3), y(1),  &
                 z(7), z(6), z(5), z(0), z(3), z(1),  &
                 dvdx(4), dvdy(4), dvdz(4))
    CALL VoluDer(x(4), x(7), x(6), x(1), x(0), x(2),  &
                 y(4), y(7), y(6), y(1), y(0), y(2),  &
                 z(4), z(7), z(6), z(1), z(0), z(2),  &
                 dvdx(5), dvdy(5), dvdz(5))
    CALL VoluDer(x(5), x(4), x(7), x(2), x(1), x(3),  &
                 y(5), y(4), y(7), y(2), y(1), y(3),  &
                 z(5), z(4), z(7), z(2), z(1), z(3),  &
                 dvdx(6), dvdy(6), dvdz(6))
    CALL VoluDer(x(6), x(5), x(4), x(3), x(2), x(0),  &
                 y(6), y(5), y(4), y(3), y(2), y(0),  &
                 z(6), z(5), z(4), z(3), z(2), z(0),  &
                 dvdx(7), dvdy(7), dvdz(7))
  
  END SUBROUTINE CalcElemVolumeDerivative
  
  SUBROUTINE CalcFBHourglassForceForElems(determ,           &
                                          x8n, y8n, z8n,    &
                                          dvdx, dvdy, dvdz, &
                                          hourg             )
  
  ! *************************************************
  ! *
  ! *     FUNCTION: Calculates the Flanagan-Belytschko anti-hourglass
  ! *               force.
  ! *
  ! *************************************************
  
  
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(0:) :: determ
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: x8n, y8n, z8n
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dvdx, dvdy, dvdz
    REAL(KIND=8) :: hourg
  
    REAL(KIND=8) :: coefficient, volinv, ss1, mass1, volume13
    REAL(KIND=8) :: hourmodx, hourmody, hourmodz
    REAL(KIND=8) :: fx, fy, fz
    REAL(KIND=8), DIMENSION(0:7) :: hgfx, hgfy, hgfz
    REAL(KIND=8), DIMENSION(0:7) :: xd1, yd1, zd1
    REAL(KIND=8), DIMENSION(0:3) :: hourgam0, hourgam1, hourgam2, hourgam3
    REAL(KIND=8), DIMENSION(0:3) :: hourgam4, hourgam5, hourgam6, hourgam7
    REAL(KIND=8), DIMENSION(0:7,0:3) :: gamma
    REAL(KIND=8), DIMENSION(:), POINTER :: fx_elem, fy_elem, fz_elem
    REAL(KIND=8), DIMENSION(:), POINTER :: fx_local, fy_local, fz_local
    INTEGER(KIND=4) :: numElem, numElem8, i, i2, i3, i1
    INTEGER(KIND=4) :: numNode, gnode, elem, count, start
    INTEGER(KIND=4) :: n0si2, n1si2, n2si2, n3si2, n4si2, n5si2, n6si2, n7si2
    INTEGER(KIND=4), DIMENSION(:), POINTER :: elemToNode => NULL()
  
    NULLIFY(fx_local, fy_local, fz_local)
    NULLIFY(fx_elem, fy_elem, fz_elem)
    numElem = domain%m_numElem
    numElem8 = numElem * 8
    ALLOCATE(fx_elem(0:numElem8-1))
    ALLOCATE(fy_elem(0:numElem8-1))
    ALLOCATE(fz_elem(0:numElem8-1))
  
    gamma(0,0) = ( 1.0_RLK)
    gamma(1,0) = ( 1.0_RLK)
    gamma(2,0) = (-1.0_RLK)
    gamma(3,0) = (-1.0_RLK)
    gamma(4,0) = (-1.0_RLK)
    gamma(5,0) = (-1.0_RLK)
    gamma(6,0) = ( 1.0_RLK)
    gamma(7,0) = ( 1.0_RLK)
    gamma(0,1) = ( 1.0_RLK)
    gamma(1,1) = (-1.0_RLK)
    gamma(2,1) = (-1.0_RLK)
    gamma(3,1) = ( 1.0_RLK)
    gamma(4,1) = (-1.0_RLK)
    gamma(5,1) = ( 1.0_RLK)
    gamma(6,1) = ( 1.0_RLK)
    gamma(7,1) = (-1.0_RLK)
    gamma(0,2) = ( 1.0_RLK)
    gamma(1,2) = (-1.0_RLK)
    gamma(2,2) = ( 1.0_RLK)
    gamma(3,2) = (-1.0_RLK)
    gamma(4,2) = ( 1.0_RLK)
    gamma(5,2) = (-1.0_RLK)
    gamma(6,2) = ( 1.0_RLK)
    gamma(7,2) = (-1.0_RLK)
    gamma(0,3) = (-1.0_RLK)
    gamma(1,3) = ( 1.0_RLK)
    gamma(2,3) = (-1.0_RLK)
    gamma(3,3) = ( 1.0_RLK)
    gamma(4,3) = ( 1.0_RLK)
    gamma(5,3) = (-1.0_RLK)
    gamma(6,3) = ( 1.0_RLK)
    gamma(7,3) = (-1.0_RLK)
    
  ! *************************************************
  ! compute the hourglass modes
    
    
    DO i2=0, numElem-1
  
      elemToNode => domain%m_nodelist(i2*8:)
  
      i3=8*i2
      volinv= (1.0_RLK)/determ(i2)
  
      DO i1=0, 3
  
        hourmodx =                                             &
          x8n(i3)   * gamma(0,i1) + x8n(i3+1) * gamma(1,i1) +  &
          x8n(i3+2) * gamma(2,i1) + x8n(i3+3) * gamma(3,i1) +  &
          x8n(i3+4) * gamma(4,i1) + x8n(i3+5) * gamma(5,i1) +  &
          x8n(i3+6) * gamma(6,i1) + x8n(i3+7) * gamma(7,i1)
  
        hourmody =                                             &
          y8n(i3)   * gamma(0,i1) + y8n(i3+1) * gamma(1,i1) +  &
          y8n(i3+2) * gamma(2,i1) + y8n(i3+3) * gamma(3,i1) +  &
          y8n(i3+4) * gamma(4,i1) + y8n(i3+5) * gamma(5,i1) +  &
          y8n(i3+6) * gamma(6,i1) + y8n(i3+7) * gamma(7,i1)
  
        hourmodz =                                             &
          z8n(i3)   * gamma(0,i1) + z8n(i3+1) * gamma(1,i1) +  &
          z8n(i3+2) * gamma(2,i1) + z8n(i3+3) * gamma(3,i1) +  &
          z8n(i3+4) * gamma(4,i1) + z8n(i3+5) * gamma(5,i1) +  &
          z8n(i3+6) * gamma(6,i1) + z8n(i3+7) * gamma(7,i1)
  
        hourgam0(i1) = gamma(0,i1) -  volinv*(dvdx(i3  ) * hourmodx +  &
                      dvdy(i3  ) * hourmody + dvdz(i3  ) * hourmodz )
  
        hourgam1(i1) = gamma(1,i1) -  volinv*(dvdx(i3+1) * hourmodx +  &
                      dvdy(i3+1) * hourmody + dvdz(i3+1) * hourmodz )
  
        hourgam2(i1) = gamma(2,i1) -  volinv*(dvdx(i3+2) * hourmodx +  &
                      dvdy(i3+2) * hourmody + dvdz(i3+2) * hourmodz )
  
        hourgam3(i1) = gamma(3,i1) -  volinv*(dvdx(i3+3) * hourmodx +  &
                      dvdy(i3+3) * hourmody + dvdz(i3+3) * hourmodz )
  
        hourgam4(i1) = gamma(4,i1) -  volinv*(dvdx(i3+4) * hourmodx +  &
                      dvdy(i3+4) * hourmody + dvdz(i3+4) * hourmodz )
  
        hourgam5(i1) = gamma(5,i1) -  volinv*(dvdx(i3+5) * hourmodx +  &
                      dvdy(i3+5) * hourmody + dvdz(i3+5) * hourmodz )
  
        hourgam6(i1) = gamma(6,i1) -  volinv*(dvdx(i3+6) * hourmodx +  &
                      dvdy(i3+6) * hourmody + dvdz(i3+6) * hourmodz )
  
        hourgam7(i1) = gamma(7,i1) -  volinv*(dvdx(i3+7) * hourmodx +  &
                      dvdy(i3+7) * hourmody + dvdz(i3+7) * hourmodz )
  
      ENDDO
  
  !   compute forces
  !   store forces into h arrays (force arrays)
  
      ss1=domain%m_ss(i2)
      mass1=domain%m_elemMass(i2)
      volume13=CBRT(determ(i2))
  
      n0si2 = elemToNode(1)
      n1si2 = elemToNode(2)
      n2si2 = elemToNode(3)
      n3si2 = elemToNode(4)
      n4si2 = elemToNode(5)
      n5si2 = elemToNode(6)
      n6si2 = elemToNode(7)
      n7si2 = elemToNode(8)
  
      xd1(0) = domain%m_xd(n0si2)
      xd1(1) = domain%m_xd(n1si2)
      xd1(2) = domain%m_xd(n2si2)
      xd1(3) = domain%m_xd(n3si2)
      xd1(4) = domain%m_xd(n4si2)
      xd1(5) = domain%m_xd(n5si2)
      xd1(6) = domain%m_xd(n6si2)
      xd1(7) = domain%m_xd(n7si2)
  
      yd1(0) = domain%m_yd(n0si2)
      yd1(1) = domain%m_yd(n1si2)
      yd1(2) = domain%m_yd(n2si2)
      yd1(3) = domain%m_yd(n3si2)
      yd1(4) = domain%m_yd(n4si2)
      yd1(5) = domain%m_yd(n5si2)
      yd1(6) = domain%m_yd(n6si2)
      yd1(7) = domain%m_yd(n7si2)
  
      zd1(0) = domain%m_zd(n0si2)
      zd1(1) = domain%m_zd(n1si2)
      zd1(2) = domain%m_zd(n2si2)
      zd1(3) = domain%m_zd(n3si2)
      zd1(4) = domain%m_zd(n4si2)
      zd1(5) = domain%m_zd(n5si2)
      zd1(6) = domain%m_zd(n6si2)
      zd1(7) = domain%m_zd(n7si2)
  
      coefficient = - hourg * (0.01_RLK) * ss1 * mass1 / volume13
  
      CALL CalcElemFBHourglassForce(xd1,yd1,zd1,                          &
                                    hourgam0,hourgam1,hourgam2,hourgam3,  &
                                    hourgam4,hourgam5,hourgam6,hourgam7,  &
                                    coefficient, hgfx, hgfy, hgfz)
  
      fx_local(0:) => fx_elem(i3:)
      fx_local(0) = hgfx(0)
      fx_local(1) = hgfx(1)
      fx_local(2) = hgfx(2)
      fx_local(3) = hgfx(3)
      fx_local(4) = hgfx(4)
      fx_local(5) = hgfx(5)
      fx_local(6) = hgfx(6)
      fx_local(7) = hgfx(7)
  
      fy_local(0:) => fy_elem(i3:)
      fy_local(0) = hgfy(0)
      fy_local(1) = hgfy(1)
      fy_local(2) = hgfy(2)
      fy_local(3) = hgfy(3)
      fy_local(4) = hgfy(4)
      fy_local(5) = hgfy(5)
      fy_local(6) = hgfy(6)
      fy_local(7) = hgfy(7)
  
      fz_local(0:) => fz_elem(i3:)
      fz_local(0) = hgfz(0)
      fz_local(1) = hgfz(1)
      fz_local(2) = hgfz(2)
      fz_local(3) = hgfz(3)
      fz_local(4) = hgfz(4)
      fz_local(5) = hgfz(5)
      fz_local(6) = hgfz(6)
      fz_local(7) = hgfz(7)
  
! #if 0
!     domain%m_fx(n0si2) = domain%m_fx(n0si2) + hgfx(0)
!     domain%m_fy(n0si2) = domain%m_fy(n0si2) + hgfy(0)
!     domain%m_fz(n0si2) = domain%m_fz(n0si2) + hgfz(0)
! 
!     domain%m_fx(n1si2) = domain%m_fx(n1si2) + hgfx(1)
!     domain%m_fy(n1si2) = domain%m_fy(n1si2) + hgfy(1)
!     domain%m_fz(n1si2) = domain%m_fz(n1si2) + hgfz(1)
! 
!     domain%m_fx(n2si2) = domain%m_fx(n2si2) + hgfx(2)
!     domain%m_fy(n2si2) = domain%m_fy(n2si2) + hgfy(2)
!     domain%m_fz(n2si2) = domain%m_fz(n2si2) + hgfz(2)
! 
!     domain%m_fx(n3si2) = domain%m_fx(n3si2) + hgfx(3)
!     domain%m_fy(n3si2) = domain%m_fy(n3si2) + hgfy(3)
!     domain%m_fz(n3si2) = domain%m_fz(n3si2) + hgfz(3)
! 
!     domain%m_fx(n4si2) = domain%m_fx(n4si2) + hgfx(4)
!     domain%m_fy(n4si2) = domain%m_fy(n4si2) + hgfy(4)
!     domain%m_fz(n4si2) = domain%m_fz(n4si2) + hgfz(4)
! 
!     domain%m_fx(n5si2) = domain%m_fx(n5si2) + hgfx(5)
!     domain%m_fy(n5si2) = domain%m_fy(n5si2) + hgfy(5)
!     domain%m_fz(n5si2) = domain%m_fz(n5si2) + hgfz(5)
! 
!     domain%m_fx(n6si2) = domain%m_fx(n6si2) + hgfx(6)
!     domain%m_fy(n6si2) = domain%m_fy(n6si2) + hgfy(6)
!     domain%m_fz(n6si2) = domain%m_fz(n6si2) + hgfz(6)
! 
!     domain%m_fx(n7si2) = domain%m_fx(n7si2) + hgfx(7)
!     domain%m_fy(n7si2) = domain%m_fy(n7si2) + hgfy(7)
!     domain%m_fz(n7si2) = domain%m_fz(n7si2) + hgfz(7)
! #endif
    ENDDO
  
    numNode = domain%m_numNode
  
    DO gnode=0, numNode-1
  
      count = domain%m_nodeElemCount(gnode)
      start = domain%m_nodeElemStart(gnode)
      fx = (0.0_RLK)
      fy = (0.0_RLK)
      fz = (0.0_RLK)
      DO i = 0, count-1
        elem = domain%m_nodeElemCornerList(start+i)
        fx = fx + fx_elem(elem)
        fy = fy + fy_elem(elem)
        fz = fz + fz_elem(elem)
      ENDDO
      domain%m_fx(gnode) = domain%m_fx(gnode) + fx
      domain%m_fy(gnode) = domain%m_fy(gnode) + fy
      domain%m_fz(gnode) = domain%m_fz(gnode) + fz
    ENDDO
  
    DEALLOCATE(fz_elem)
    DEALLOCATE(fy_elem)
    DEALLOCATE(fx_elem)
  
  
  END SUBROUTINE CalcFBHourglassForceForElems
  
  SUBROUTINE CalcElemNodeNormals( pfx,pfy, pfz, x, y, z  )
  
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(0:) :: pfx,pfy,pfz
    REAL(KIND=8), DIMENSION(0:) :: x, y, z 
    integer ::  i 
  
    DO i = 0, 7
      pfx(i) = 0.0_RLK
      pfy(i) = 0.0_RLK
      pfz(i) = 0.0_RLK
    ENDDO
  
  ! evaluate face one: nodes 0, 1, 2, 3
    CALL SumElemFaceNormal(pfx(0), pfy(0), pfz(0),              &
                           pfx(1), pfy(1), pfz(1),              &
                           pfx(2), pfy(2), pfz(2),              &
                           pfx(3), pfy(3), pfz(3),              &
                           x(0), y(0), z(0), x(1), y(1), z(1),  &
                           x(2), y(2), z(2), x(3), y(3), z(3))
  ! evaluate face two: nodes 0, 4, 5, 1 */
    CALL SumElemFaceNormal(pfx(0), pfy(0), pfz(0),              &
                           pfx(4), pfy(4), pfz(4),              &
                           pfx(5), pfy(5), pfz(5),              &
                           pfx(1), pfy(1), pfz(1),              &
                           x(0), y(0), z(0), x(4), y(4), z(4),  &
                           x(5), y(5), z(5), x(1), y(1), z(1))
  ! evaluate face three: nodes 1, 5, 6, 2 */
    CALL SumElemFaceNormal(pfx(1), pfy(1), pfz(1),              &
                           pfx(5), pfy(5), pfz(5),              &
                           pfx(6), pfy(6), pfz(6),              &
                           pfx(2), pfy(2), pfz(2),              &
                           x(1), y(1), z(1), x(5), y(5), z(5),  &
                           x(6), y(6), z(6), x(2), y(2), z(2))
  ! evaluate face four: nodes 2, 6, 7, 3 */
    CALL SumElemFaceNormal(pfx(2), pfy(2), pfz(2),              &
                           pfx(6), pfy(6), pfz(6),              &
                           pfx(7), pfy(7), pfz(7),              &
                           pfx(3), pfy(3), pfz(3),              &
                           x(2), y(2), z(2), x(6), y(6), z(6),  &
                           x(7), y(7), z(7), x(3), y(3), z(3))
  ! evaluate face five: nodes 3, 7, 4, 0 */
    CALL SumElemFaceNormal(pfx(3), pfy(3), pfz(3),              &
                           pfx(7), pfy(7), pfz(7),              &
                           pfx(4), pfy(4), pfz(4),              &
                           pfx(0), pfy(0), pfz(0),              &
                           x(3), y(3), z(3), x(7), y(7), z(7),  &
                           x(4), y(4), z(4), x(0), y(0), z(0))
  ! evaluate face six: nodes 4, 7, 6, 5 */
    CALL SumElemFaceNormal(pfx(4), pfy(4), pfz(4),              &
                           pfx(7), pfy(7), pfz(7),              &
                           pfx(6), pfy(6), pfz(6),              &
                           pfx(5), pfy(5), pfz(5),              &
                           x(4), y(4), z(4), x(7), y(7), z(7),  &
                           x(6), y(6), z(6), x(5), y(5), z(5))
  
  END SUBROUTINE CalcElemNodeNormals

  SUBROUTINE SumElemStressesToNodeForces(B, stress_xx, stress_yy, stress_zz,  fx,  fy,  fz)
  
    IMPLICIT NONE 
    REAL(KIND=8) ,DIMENSION(0:7,0:2) :: B ! alloc 2nd dim to 8 or 0:7
    REAL(KIND=8) :: stress_xx, stress_yy, stress_zz
    REAL(KIND=8), DIMENSION(0:7) ::  fx,  fy,  fz 
  
    REAL(KIND=8) :: pfx0, pfx1, pfx2, pfx3, pfx4, pfx5, pfx6, pfx7
    REAL(KIND=8) :: pfy0, pfy1, pfy2, pfy3, pfy4, pfy5, pfy6, pfy7
    REAL(KIND=8) :: pfz0, pfz1, pfz2, pfz3, pfz4, pfz5, pfz6, pfz7
  
    pfx0 = B(0,0)
    pfx1 = B(1,0)
    pfx2 = B(2,0)
    pfx3 = B(3,0)
    pfx4 = B(4,0)
    pfx5 = B(5,0)
    pfx6 = B(6,0)
    pfx7 = B(7,0)
  
    pfy0 = B(0,1)
    pfy1 = B(1,1)
    pfy2 = B(2,1)
    pfy3 = B(3,1)
    pfy4 = B(4,1)
    pfy5 = B(5,1)
    pfy6 = B(6,1)
    pfy7 = B(7,1)
  
    pfz0 = B(0,2)
    pfz1 = B(1,2)
    pfz2 = B(2,2)
    pfz3 = B(3,2)
    pfz4 = B(4,2)
    pfz5 = B(5,2)
    pfz6 = B(6,2)
    pfz7 = B(7,2)
  
    fx(0) = -( stress_xx * pfx0 )
    fx(1) = -( stress_xx * pfx1 )
    fx(2) = -( stress_xx * pfx2 )
    fx(3) = -( stress_xx * pfx3 )
    fx(4) = -( stress_xx * pfx4 )
    fx(5) = -( stress_xx * pfx5 )
    fx(6) = -( stress_xx * pfx6 )
    fx(7) = -( stress_xx * pfx7 )
    
    fy(0) = -( stress_yy * pfy0  )
    fy(1) = -( stress_yy * pfy1  )
    fy(2) = -( stress_yy * pfy2  )
    fy(3) = -( stress_yy * pfy3  )
    fy(4) = -( stress_yy * pfy4  )
    fy(5) = -( stress_yy * pfy5  )
    fy(6) = -( stress_yy * pfy6  )
    fy(7) = -( stress_yy * pfy7  )
    
    fz(0) = -( stress_zz * pfz0 )
    fz(1) = -( stress_zz * pfz1 )
    fz(2) = -( stress_zz * pfz2 )
    fz(3) = -( stress_zz * pfz3 )
    fz(4) = -( stress_zz * pfz4 )
    fz(5) = -( stress_zz * pfz5 )
    fz(6) = -( stress_zz * pfz6 )
    fz(7) = -( stress_zz * pfz7 )
  
  END SUBROUTINE SumElemStressesToNodeForces

  SUBROUTINE SumElemFaceNormal(normalX0, normalY0, normalZ0, &
                               normalX1, normalY1, normalZ1, &
                               normalX2, normalY2, normalZ2, &
                               normalX3, normalY3, normalZ3, &
                                x0,  y0,  z0,    &
                                x1,  y1,  z1,    &
                                x2,  y2,  z2,    &
                                x3,  y3,  z3     )
  
    IMPLICIT NONE
  
    REAL(KIND=8) :: normalX0,normalY0,normalZ0
    REAL(KIND=8) :: normalX1,normalY1,normalZ1
    REAL(KIND=8) :: normalX2,normalY2,normalZ2
    REAL(KIND=8) :: normalX3,normalY3,normalZ3
    REAL(KIND=8) :: x0,y0,z0
    REAL(KIND=8) :: x1,y1,z1
    REAL(KIND=8) :: x2,y2,z2  
    REAL(KIND=8) :: x3,y3,z3
   
    REAL(KIND=8) :: bisectX0
    REAL(KIND=8) :: bisectY0
    REAL(KIND=8) :: bisectZ0
    REAL(KIND=8) :: bisectX1
    REAL(KIND=8) :: bisectY1
    REAL(KIND=8) :: bisectZ1
    REAL(KIND=8) :: areaX
    REAL(KIND=8) :: areaY
    REAL(KIND=8) :: areaZ
    REAL(KIND=8), PARAMETER :: RHALF = 0.5_RLK
    REAL(KIND=8), PARAMETER :: RQTR  = 0.25_RLK
  
    bisectX0 = RHALF * (x3 + x2 - x1 - x0)
    bisectY0 = RHALF * (y3 + y2 - y1 - y0)
    bisectZ0 = RHALF * (z3 + z2 - z1 - z0)
    bisectX1 = RHALF * (x2 + x1 - x3 - x0)
    bisectY1 = RHALF * (y2 + y1 - y3 - y0)
    bisectZ1 = RHALF * (z2 + z1 - z3 - z0)
    areaX = RQTR * (bisectY0 * bisectZ1 - bisectZ0 * bisectY1)
    areaY = RQTR * (bisectZ0 * bisectX1 - bisectX0 * bisectZ1)
    areaZ = RQTR * (bisectX0 * bisectY1 - bisectY0 * bisectX1)
  
    normalX0 = normalX0 + areaX
    normalX1 = normalX1 + areaX
    normalX2 = normalX2 + areaX
    normalX3 = normalX3 + areaX
  
    normalY0 = normalY0 + areaY
    normalY1 = normalY1 + areaY
    normalY2 = normalY2 + areaY
    normalY3 = normalY3 + areaY
  
    normalZ0 = normalZ0 + areaZ
    normalZ1 = normalZ1 + areaZ
    normalZ2 = normalZ2 + areaZ
    normalZ3 = normalZ3 + areaZ
  
  END SUBROUTINE SumElemFaceNormal

  SUBROUTINE CalcElemFBHourglassForce(xd, yd, zd, &
                                      hourgam0, hourgam1, &
                                      hourgam2, hourgam3, &
                                      hourgam4, hourgam5, &
                                      hourgam6, hourgam7, &
                                      coefficient, hgfx,  &
                                      hgfy, hgfz          )
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(0:7) :: xd,yd,zd
    REAL(KIND=8), DIMENSION(0:3) :: hourgam0, hourgam1,  &
                                    hourgam2, hourgam3,  &
                                    hourgam4, hourgam5,  &
                                    hourgam6, hourgam7
    REAL(KIND=8) :: coefficient
    REAL(KIND=8),DIMENSION(0:7) :: hgfx,hgfy,hgfz
    REAL(KIND=8) :: h00,h01,h02,h03
  
    INTEGER(KIND=4),PARAMETER :: i00 = 0_4
    INTEGER(KIND=4),PARAMETER :: i01 = 1_4
    INTEGER(KIND=4),PARAMETER :: i02 = 2_4
    INTEGER(KIND=4),PARAMETER :: i03 = 3_4
  
    h00 =                                             &
      hourgam0(i00) * xd(0) + hourgam1(i00) * xd(1) + &
      hourgam2(i00) * xd(2) + hourgam3(i00) * xd(3) + &
      hourgam4(i00) * xd(4) + hourgam5(i00) * xd(5) + &
      hourgam6(i00) * xd(6) + hourgam7(i00) * xd(7)
    
    h01 =                                             &
      hourgam0(i01) * xd(0) + hourgam1(i01) * xd(1) + &
      hourgam2(i01) * xd(2) + hourgam3(i01) * xd(3) + &
      hourgam4(i01) * xd(4) + hourgam5(i01) * xd(5) + &
      hourgam6(i01) * xd(6) + hourgam7(i01) * xd(7)
    
    h02 =                                             &
      hourgam0(i02) * xd(0) + hourgam1(i02) * xd(1) + &
      hourgam2(i02) * xd(2) + hourgam3(i02) * xd(3) + &
      hourgam4(i02) * xd(4) + hourgam5(i02) * xd(5) + &
      hourgam6(i02) * xd(6) + hourgam7(i02) * xd(7)
    
    h03 =                                             &
      hourgam0(i03) * xd(0) + hourgam1(i03) * xd(1) + &
      hourgam2(i03) * xd(2) + hourgam3(i03) * xd(3) + &
      hourgam4(i03) * xd(4) + hourgam5(i03) * xd(5) + &
      hourgam6(i03) * xd(6) + hourgam7(i03) * xd(7)
    
    hgfx(0) = coefficient *                       &
     (hourgam0(i00) * h00 + hourgam0(i01) * h01 + &
      hourgam0(i02) * h02 + hourgam0(i03) * h03)
    
    hgfx(1) = coefficient *                       &
     (hourgam1(i00) * h00 + hourgam1(i01) * h01 + &
      hourgam1(i02) * h02 + hourgam1(i03) * h03)
    
    hgfx(2) = coefficient *                       &
     (hourgam2(i00) * h00 + hourgam2(i01) * h01 + &
      hourgam2(i02) * h02 + hourgam2(i03) * h03)
    
    hgfx(3) = coefficient *                       &
     (hourgam3(i00) * h00 + hourgam3(i01) * h01 + &
      hourgam3(i02) * h02 + hourgam3(i03) * h03)
    
    hgfx(4) = coefficient *                       &
     (hourgam4(i00) * h00 + hourgam4(i01) * h01 + &
      hourgam4(i02) * h02 + hourgam4(i03) * h03)
    
    hgfx(5) = coefficient *                       &
     (hourgam5(i00) * h00 + hourgam5(i01) * h01 + &
      hourgam5(i02) * h02 + hourgam5(i03) * h03)
    
    hgfx(6) = coefficient *                       &
     (hourgam6(i00) * h00 + hourgam6(i01) * h01 + &
      hourgam6(i02) * h02 + hourgam6(i03) * h03)
    
    hgfx(7) = coefficient *                       &
     (hourgam7(i00) * h00 + hourgam7(i01) * h01 + &
      hourgam7(i02) * h02 + hourgam7(i03) * h03)
    
    h00 =                                             &
      hourgam0(i00) * yd(0) + hourgam1(i00) * yd(1) + &
      hourgam2(i00) * yd(2) + hourgam3(i00) * yd(3) + &
      hourgam4(i00) * yd(4) + hourgam5(i00) * yd(5) + &
      hourgam6(i00) * yd(6) + hourgam7(i00) * yd(7)
    
    h01 =                                             &
      hourgam0(i01) * yd(0) + hourgam1(i01) * yd(1) + &
      hourgam2(i01) * yd(2) + hourgam3(i01) * yd(3) + &
      hourgam4(i01) * yd(4) + hourgam5(i01) * yd(5) + &
      hourgam6(i01) * yd(6) + hourgam7(i01) * yd(7)
    
    h02 =                                            &
      hourgam0(i02) * yd(0) + hourgam1(i02) * yd(1)+ &
      hourgam2(i02) * yd(2) + hourgam3(i02) * yd(3)+ &
      hourgam4(i02) * yd(4) + hourgam5(i02) * yd(5)+ &
      hourgam6(i02) * yd(6) + hourgam7(i02) * yd(7)
    
    h03 =                                             &
      hourgam0(i03) * yd(0) + hourgam1(i03) * yd(1) + &
      hourgam2(i03) * yd(2) + hourgam3(i03) * yd(3) + &
      hourgam4(i03) * yd(4) + hourgam5(i03) * yd(5) + &
      hourgam6(i03) * yd(6) + hourgam7(i03) * yd(7)
    
    
    hgfy(0) = coefficient *                       &
     (hourgam0(i00) * h00 + hourgam0(i01) * h01 + &
      hourgam0(i02) * h02 + hourgam0(i03) * h03)
    
    hgfy(1) = coefficient *                       &
     (hourgam1(i00) * h00 + hourgam1(i01) * h01 + &
      hourgam1(i02) * h02 + hourgam1(i03) * h03)
    
    hgfy(2) = coefficient *                       &
     (hourgam2(i00) * h00 + hourgam2(i01) * h01 + &
      hourgam2(i02) * h02 + hourgam2(i03) * h03)
    
    hgfy(3) = coefficient *                       &
     (hourgam3(i00) * h00 + hourgam3(i01) * h01 + &
      hourgam3(i02) * h02 + hourgam3(i03) * h03)
    
    hgfy(4) = coefficient *                       &
     (hourgam4(i00) * h00 + hourgam4(i01) * h01 + &
      hourgam4(i02) * h02 + hourgam4(i03) * h03)
    
    hgfy(5) = coefficient *                       &
     (hourgam5(i00) * h00 + hourgam5(i01) * h01 + &
      hourgam5(i02) * h02 + hourgam5(i03) * h03)
    
    hgfy(6) = coefficient *                       &
     (hourgam6(i00) * h00 + hourgam6(i01) * h01 + &
      hourgam6(i02) * h02 + hourgam6(i03) * h03)
    
    hgfy(7) = coefficient *                       &
     (hourgam7(i00) * h00 + hourgam7(i01) * h01 + &
      hourgam7(i02) * h02 + hourgam7(i03) * h03)
    
    h00 =                                              &
      hourgam0(i00) * zd(0) + hourgam1(i00) * zd(1) +  &
      hourgam2(i00) * zd(2) + hourgam3(i00) * zd(3) +  &
      hourgam4(i00) * zd(4) + hourgam5(i00) * zd(5) +  &
      hourgam6(i00) * zd(6) + hourgam7(i00) * zd(7)
    
    h01 =                                              &
      hourgam0(i01) * zd(0) + hourgam1(i01) * zd(1) +  &
      hourgam2(i01) * zd(2) + hourgam3(i01) * zd(3) +  &
      hourgam4(i01) * zd(4) + hourgam5(i01) * zd(5) +  &
      hourgam6(i01) * zd(6) + hourgam7(i01) * zd(7)
    
    h02 =                                              &
      hourgam0(i02) * zd(0) + hourgam1(i02) * zd(1)+   &
      hourgam2(i02) * zd(2) + hourgam3(i02) * zd(3)+   &
      hourgam4(i02) * zd(4) + hourgam5(i02) * zd(5)+   &
      hourgam6(i02) * zd(6) + hourgam7(i02) * zd(7)
    
    h03 =                                              &
      hourgam0(i03) * zd(0) + hourgam1(i03) * zd(1) +  &
      hourgam2(i03) * zd(2) + hourgam3(i03) * zd(3) +  &
      hourgam4(i03) * zd(4) + hourgam5(i03) * zd(5) +  &
      hourgam6(i03) * zd(6) + hourgam7(i03) * zd(7)
    
    
    hgfz(0) = coefficient *                        &
     (hourgam0(i00) * h00 + hourgam0(i01) * h01 +  &
      hourgam0(i02) * h02 + hourgam0(i03) * h03)
    
    hgfz(1) = coefficient *                        &
     (hourgam1(i00) * h00 + hourgam1(i01) * h01 +  &
      hourgam1(i02) * h02 + hourgam1(i03) * h03)
    
    hgfz(2) = coefficient *                        &
     (hourgam2(i00) * h00 + hourgam2(i01) * h01 +  &
      hourgam2(i02) * h02 + hourgam2(i03) * h03)
    
    hgfz(3) = coefficient *                        &
     (hourgam3(i00) * h00 + hourgam3(i01) * h01 +  &
      hourgam3(i02) * h02 + hourgam3(i03) * h03)
    
    hgfz(4) = coefficient *                        &
     (hourgam4(i00) * h00 + hourgam4(i01) * h01 +  &
      hourgam4(i02) * h02 + hourgam4(i03) * h03)
    
    hgfz(5) = coefficient *                        &
     (hourgam5(i00) * h00 + hourgam5(i01) * h01 +  &
      hourgam5(i02) * h02 + hourgam5(i03) * h03)
    
    hgfz(6) = coefficient *                        &
     (hourgam6(i00) * h00 + hourgam6(i01) * h01 +  &
      hourgam6(i02) * h02 + hourgam6(i03) * h03)
    
    hgfz(7) = coefficient *                        &
     (hourgam7(i00) * h00 + hourgam7(i01) * h01 +  &
      hourgam7(i02) * h02 + hourgam7(i03) * h03)
  
  END SUBROUTINE CalcElemFBHourglassForce

  SUBROUTINE VoluDer(x0, x1, x2,      &
                     x3, x4, x5,      &
                     y0, y1, y2,      &
                     y3, y4, y5,      &
                     z0, z1, z2,      &
                     z3, z4, z5,      &
                     dvdx, dvdy, dvdz )
    IMPLICIT NONE
    REAL(KIND=8) :: x0, x1, x2, x3, x4, x5
    REAL(KIND=8) :: y0, y1, y2, y3, y4, y5
    REAL(KIND=8) :: z0, z1, z2, z3, z4, z5
    REAL(KIND=8) :: dvdx, dvdy, dvdz
  
    REAL(KIND=8), PARAMETER :: twelfth = 1.0_RLK / 12.0_RLK
  
    dvdx =                                              &
      (y1 + y2) * (z0 + z1) - (y0 + y1) * (z1 + z2) +   &
      (y0 + y4) * (z3 + z4) - (y3 + y4) * (z0 + z4) -   &
      (y2 + y5) * (z3 + z5) + (y3 + y5) * (z2 + z5)
  
    dvdy =                                              &
      - (x1 + x2) * (z0 + z1) + (x0 + x1) * (z1 + z2) - &
      (x0 + x4) * (z3 + z4) + (x3 + x4) * (z0 + z4) +   &
      (x2 + x5) * (z3 + z5) - (x3 + x5) * (z2 + z5)
  
    dvdz =                                              &
      - (y1 + y2) * (x0 + x1) + (y0 + y1) * (x1 + x2) - &
      (y0 + y4) * (x3 + x4) + (y3 + y4) * (x0 + x4) +   &
      (y2 + y5) * (x3 + x5) - (y3 + y5) * (x2 + x5)
  
    dvdx = dvdx * twelfth
    dvdy = dvdy * twelfth
    dvdz = dvdz * twelfth
  
  END SUBROUTINE VoluDer

end module calc_force
