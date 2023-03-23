module utils
  use main_data, only: RLK
  contains
  SUBROUTINE luabort(errcode)
    IMPLICIT NONE
    INTEGER(KIND=4) :: errcode
    WRITE(6,*) "ERROR CODE: ", errcode
    WRITE(6,*) "ABORTING"
    STOP
  END SUBROUTINE luabort
  
  REAL(KIND=8) FUNCTION SUM4(a, b, c, d)
  
    IMPLICIT NONE 
    REAL(KIND=8) :: a, b, c, d
  
    SUM4 = a + b + c + d
  
    RETURN
  
  END FUNCTION SUM4
  
  REAL(KIND=8) FUNCTION TRIPLE_PRODUCT(x1, y1, z1, x2, y2, z2, x3, y3, z3)
  
    REAL(KIND=8) :: x1, y1, z1, x2, y2, z2, x3, y3, z3
  
    TRIPLE_PRODUCT = ((x1)*((y2)*(z3) - (z2)*(y3)) + (x2)*((z1)*(y3)  &
                    - (y1)*(z3)) + (x3)*((y1)*(z2) - (z1)*(y2)))
  
    RETURN
  
  END FUNCTION TRIPLE_PRODUCT
  
  REAL(KIND=8) FUNCTION CalcElemVolume( x, y, z )
  
    IMPLICIT NONE
    REAL(KIND=8) ,DIMENSION(0:7) :: x, y, z
  
    REAL(KIND=8)  :: volume=0.0_RLK
  
  !!$REAL(KIND=8)  :: x0, x1, x2, x3, x4, x5, x6, x7
  !!$REAL(KIND=8)  :: y0, y1, y2, y3, y4, y5, y6, y7
  !!$REAL(KIND=8)  :: z0, z1, z2, z3, z4, z5, z6, z7
  
    REAL(KIND=8) :: twelveth = (1.0_RLK)/(12.0_RLK)
  
    REAL(KIND=8) :: dx61
    REAL(KIND=8) :: dy61
    REAL(KIND=8) :: dz61
  
    REAL(KIND=8) :: dx70
    REAL(KIND=8) :: dy70
    REAL(KIND=8) :: dz70
    
    REAL(KIND=8) :: dx63
    REAL(KIND=8) :: dy63 
    REAL(KIND=8) :: dz63
  
    REAL(KIND=8) :: dx20 
    REAL(KIND=8) :: dy20
    REAL(KIND=8) :: dz20
  
    REAL(KIND=8) :: dx50 
    REAL(KIND=8) :: dy50
    REAL(KIND=8) :: dz50
  
    REAL(KIND=8) :: dx64
    REAL(KIND=8) :: dy64
    REAL(KIND=8) :: dz64
  
    REAL(KIND=8) :: dx31
    REAL(KIND=8) :: dy31 
    REAL(KIND=8) :: dz31
  
    REAL(KIND=8) :: dx72
    REAL(KIND=8) :: dy72
    REAL(KIND=8) :: dz72
  
    REAL(KIND=8) :: dx43 
    REAL(KIND=8) :: dy43
    REAL(KIND=8) :: dz43
  
    REAL(KIND=8) :: dx57
    REAL(KIND=8) :: dy57
    REAL(KIND=8) :: dz57
  
    REAL(KIND=8) :: dx14
    REAL(KIND=8) :: dy14 
    REAL(KIND=8) :: dz14
  
    REAL(KIND=8) :: dx25
    REAL(KIND=8) :: dy25 
    REAL(KIND=8) :: dz25
  
  
    dx61 = x(6) - x(1)
    dy61 = y(6) - y(1)
    dz61 = z(6) - z(1)
  
    dx70 = x(7) - x(0)
    dy70 = y(7) - y(0)
    dz70 = z(7) - z(0)
  
    dx63 = x(6) - x(3)
    dy63 = y(6) - y(3)
    dz63 = z(6) - z(3)
  
    dx20 = x(2) - x(0)
    dy20 = y(2) - y(0)
    dz20 = z(2) - z(0)
  
    dx50 = x(5) - x(0)
    dy50 = y(5) - y(0)
    dz50 = z(5) - z(0)
  
    dx64 = x(6) - x(4)
    dy64 = y(6) - y(4)
    dz64 = z(6) - z(4)
  
    dx31 = x(3) - x(1)
    dy31 = y(3) - y(1)
    dz31 = z(3) - z(1)
  
    dx72 = x(7) - x(2)
    dy72 = y(7) - y(2)
    dz72 = z(7) - z(2)
  
    dx43 = x(4) - x(3)
    dy43 = y(4) - y(3)
    dz43 = z(4) - z(3)
  
    dx57 = x(5) - x(7)
    dy57 = y(5) - y(7)
    dz57 = z(5) - z(7)
  
    dx14 = x(1) - x(4)
    dy14 = y(1) - y(4)
    dz14 = z(1) - z(4)
  
    dx25 = x(2) - x(5)
    dy25 = y(2) - y(5)
    dz25 = z(2) - z(5)
  
    volume =  TRIPLE_PRODUCT(dx31 + dx72, dx63, dx20,   &
                             dy31 + dy72, dy63, dy20,   &
                             dz31 + dz72, dz63, dz20) + &
              TRIPLE_PRODUCT(dx43 + dx57, dx64, dx70,   &
                             dy43 + dy57, dy64, dy70,   &
                             dz43 + dz57, dz64, dz70) + &
              TRIPLE_PRODUCT(dx14 + dx25, dx61, dx50,   &
                             dy14 + dy25, dy61, dy50,   &
                             dz14 + dz25, dz61, dz50)
  
    volume = volume*twelveth
  
  
    CalcElemVolume=volume
  
  
  END FUNCTION CalcElemVolume
  
  FUNCTION AreaFace( x0, x1, x2, x3,  &
                     y0, y1, y2, y3,  &
                     z0, z1, z2, z3  ) RESULT(area)
  
  
    IMPLICIT NONE
    REAL(KIND=8)  :: x0, x1, x2, x3
    REAL(KIND=8)  :: y0, y1, y2, y3
    REAL(KIND=8)  :: z0, z1, z2, z3
  
    REAL(KIND=8) :: fx, fy, fz
    REAL(KIND=8) :: gx, gy, gz
    REAL(KIND=8) :: area
  
    fx = (x2 - x0) - (x3 - x1)
    fy = (y2 - y0) - (y3 - y1)
    fz = (z2 - z0) - (z3 - z1)
    gx = (x2 - x0) + (x3 - x1)
    gy = (y2 - y0) + (y3 - y1)
    gz = (z2 - z0) + (z3 - z1)
  
    area =                             &
      (fx * fx + fy * fy + fz * fz) *  &
      (gx * gx + gy * gy + gz * gz) -  &
      (fx * gx + fy * gy + fz * gz) *  &
      (fx * gx + fy * gy + fz * gz)
  
    RETURN
  
  END FUNCTION AreaFace
  
  REAL(KIND=8) FUNCTION CBRT(dat)
  
    IMPLICIT NONE
    REAL(KIND=8) :: dat
  
    CBRT = dat**(1.0_RLK/3.0_RLK)
  
  END FUNCTION CBRT
  
  FUNCTION CalcElemCharacteristicLength( x, y, z, volume) RESULT(charLength)
  
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(0:7) :: x, y, z
    REAL(KIND=8) :: volume
    REAL(KIND=8) :: a
    REAL(KIND=8) :: charLength
  
    charLength = 0.0_RLK
  
    a = AreaFace(x(0),x(1),x(2),x(3),  &
                 y(0),y(1),y(2),y(3),  &
                 z(0),z(1),z(2),z(3))
    charLength = MAX(a,charLength)
    
    a = AreaFace(x(4),x(5),x(6),x(7),  &
                 y(4),y(5),y(6),y(7),  &
                 z(4),z(5),z(6),z(7))
    charLength = MAX(a,charLength)
    
    a = AreaFace(x(0),x(1),x(5),x(4),  &
                 y(0),y(1),y(5),y(4),  &
                 z(0),z(1),z(5),z(4))
    charLength = MAX(a,charLength)
  
    a = AreaFace(x(1),x(2),x(6),x(5),  &
                 y(1),y(2),y(6),y(5),  &
                 z(1),z(2),z(6),z(5))
    charLength = MAX(a,charLength)
    
    a = AreaFace(x(2),x(3),x(7),x(6),  &
                 y(2),y(3),y(7),y(6),  &
                 z(2),z(3),z(7),z(6))
    charLength = MAX(a,charLength)
    
    a = AreaFace(x(3),x(0),x(4),x(7),  &
                 y(3),y(0),y(4),y(7),  &
                 z(3),z(0),z(4),z(7))
    charLength = MAX(a,charLength)
    
    charLength = (4.0_RLK) * volume / SQRT(charLength);
  
    RETURN
  
  END FUNCTION CalcElemCharacteristicLength
  
  SUBROUTINE CalcElemShapeFunctionDerivatives( x, y, z,   &
                                               b,         &
                                               el_volume   )
    IMPLICIT NONE 
  
    REAL(KIND=8), DIMENSION(0:7)  :: x, y, z
    REAL(KIND=8), DIMENSION(0:7,0:2) :: b ! alloc 2nd dim to 8 or 0:7
    REAL(KIND=8), INTENT(INOUT) :: el_volume
    REAL(KIND=8)  :: x0,x1,x2,x3,x4,x5,x6,x7
    REAL(KIND=8)  :: y0,y1,y2,y3,y4,y5,y6,y7
    REAL(KIND=8)  :: z0,z1,z2,z3,z4,z5,z6,z7
  
    REAL(KIND=8)  :: fjxxi, fjxet, fjxze
    REAL(KIND=8)  :: fjyxi, fjyet, fjyze
    REAL(KIND=8)  :: fjzxi, fjzet, fjzze
    REAL(KIND=8)  :: cjxxi, cjxet, cjxze
    REAL(KIND=8)  :: cjyxi, cjyet, cjyze
    REAL(KIND=8)  :: cjzxi, cjzet, cjzze
  
    x0 = x(0)
    x1 = x(1)
    x2 = x(2)
    x3 = x(3)
    x4 = x(4)
    x5 = x(5)
    x6 = x(6)
    x7 = x(7)
  
    y0 = y(0)
    y1 = y(1)
    y2 = y(2)
    y3 = y(3)
    y4 = y(4)
    y5 = y(5)
    y6 = y(6)
    y7 = y(7)
  
    z0 = z(0)
    z1 = z(1)
    z2 = z(2)
    z3 = z(3)
    z4 = z(4)
    z5 = z(5)
    z6 = z(6)
    z7 = z(7)
  
    fjxxi = .125 * ( (x6-x0) + (x5-x3) - (x7-x1) - (x4-x2) )
    fjxet = .125 * ( (x6-x0) - (x5-x3) + (x7-x1) - (x4-x2) )
    fjxze = .125 * ( (x6-x0) + (x5-x3) + (x7-x1) + (x4-x2) )
  
    fjyxi = .125 * ( (y6-y0) + (y5-y3) - (y7-y1) - (y4-y2) )
    fjyet = .125 * ( (y6-y0) - (y5-y3) + (y7-y1) - (y4-y2) )
    fjyze = .125 * ( (y6-y0) + (y5-y3) + (y7-y1) + (y4-y2) )
  
    fjzxi = .125 * ( (z6-z0) + (z5-z3) - (z7-z1) - (z4-z2) )
    fjzet = .125 * ( (z6-z0) - (z5-z3) + (z7-z1) - (z4-z2) )
    fjzze = .125 * ( (z6-z0) + (z5-z3) + (z7-z1) + (z4-z2) )
  
  ! compute cofactors
    cjxxi =    (fjyet * fjzze) - (fjzet * fjyze)
    cjxet =  - (fjyxi * fjzze) + (fjzxi * fjyze)
    cjxze =    (fjyxi * fjzet) - (fjzxi * fjyet)
  
    cjyxi =  - (fjxet * fjzze) + (fjzet * fjxze)
    cjyet =    (fjxxi * fjzze) - (fjzxi * fjxze)
    cjyze =  - (fjxxi * fjzet) + (fjzxi * fjxet)
  
    cjzxi =    (fjxet * fjyze) - (fjyet * fjxze)
    cjzet =  - (fjxxi * fjyze) + (fjyxi * fjxze)
    cjzze =    (fjxxi * fjyet) - (fjyxi * fjxet)
  
  ! calculate partials :
  !     this need only be done for l = 0,1,2,3   since , by symmetry ,
  !     (6,7,4,5) = - (0,1,2,3) .
    b(0,0) =   -  cjxxi  -  cjxet  -  cjxze
    b(1,0) =      cjxxi  -  cjxet  -  cjxze
    b(2,0) =      cjxxi  +  cjxet  -  cjxze
    b(3,0) =   -  cjxxi  +  cjxet  -  cjxze
    b(4,0) = -b(2,0)
    b(5,0) = -b(3,0)
    b(6,0) = -b(0,0)
    b(7,0) = -b(1,0)
  
    b(0,1) =   -  cjyxi  -  cjyet  -  cjyze
    b(1,1) =      cjyxi  -  cjyet  -  cjyze
    b(2,1) =      cjyxi  +  cjyet  -  cjyze
    b(3,1) =   -  cjyxi  +  cjyet  -  cjyze
    b(4,1) = -b(2,1)
    b(5,1) = -b(3,1)
    b(6,1) = -b(0,1)
    b(7,1) = -b(1,1)
  
    b(0,2) =   -  cjzxi  -  cjzet  -  cjzze
    b(1,2) =      cjzxi  -  cjzet  -  cjzze
    b(2,2) =      cjzxi  +  cjzet  -  cjzze
    b(3,2) =   -  cjzxi  +  cjzet  -  cjzze
    b(4,2) = -b(2,2)
    b(5,2) = -b(2,3)
    b(6,2) = -b(2,0)
    b(7,2) = -b(2,1)
  
  ! calculate jacobian determinant (volume)
    el_volume = 8.0_RLK * ( fjxet * cjxet + fjyet * cjyet + fjzet * cjzet)
  
  END SUBROUTINE CalcElemShapeFunctionDerivatives

end module utils
