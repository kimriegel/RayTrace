      PatchNox=4
      allocate(xlimit(PatchNox))
      allocate(Nx(patchNox-1))
      xlimit(1)=0.0
      xlimit(2)=10.0
      xlimit(3)=64.4322
      xlimit(4)=100.0

      Nx(1)=10
      Nx(2)=10
      Nx(3)=10

      qx=3.0
      
      PatchNoy=4
      allocate(ylimit(PatchNoy))
      allocate(Ny(PatchNoy-1))

      ylimit(1)=0.0
      ylimit(2)=10.0
      ylimit(3)=46.9316
      ylimit(4)=60.0

      Ny(1)=10
      Ny(2)=10
      Ny(3)=10

      qy=3.0

      PatchNoz=2
      allocate(zlimit(PatchNoz))
      allocate(Nz(PatchNoz-1))

      zlimit(1)=0.0
      zlimit(2)=8.2423

      Nz(1)=10

      qz=3.0


      PatchNo=Nx(1)*Ny(1)+Nx(2)*Ny(1)+Nx(3)*Ny(1)+Nx(1)*Ny(2)+Nx(2)*
     *     Ny(2)+Nx(3)*Ny(2)+Nx(1)*Ny(3)+Nx(2)*Ny(3)+Nx(3)*Ny(3)+2*Nx(2)
     *     *Nz(1)+2*Ny(2)*Nz(1)

      allocate(patcharray(PatchNo,sizeffttwo,10))
      allocate(formfactors(PatchNo,PatchNo,3))

      patcharray=0.0
C     This creates a patches for the Building. 
      increment=1
      slope=0
      slope1=0
      tempsize=Nx(1)*Ny(1)
      allocate(ddx1(Nx(1)))
      allocate(ddy1(Ny(1)))
      allocate(patcharraytemp(Nx(1)*Ny(1),6))
      b=ylimit(2)
      b1=ylimit(1)
      CALL PATCHESSHORT(xlimit(1),xlimit(2),Nx(1),qx,ddx1)
      CALL PATCHESSHORT(ylimit(1),ylimit(2),Ny(1),qy,ddy1)
      CALL CREATEPATCHARRAY(ddx1,ddy1,Nx(1),Ny(1),xlimit(1),xlimit(2)
     *     ,ylimit(1),ylimit(2),zlimit(1),zlimit(1),patcharraytemp,slope
     *     ,b,slope1,b1,count,FaceNormals(5,1:3))
      DO 101 Q=increment,increment+tempsize-1
         DO 102 W=1,sizeffttwo
         patcharray(Q,W,1)=patcharraytemp(Q-increment+1,1)
         patcharray(Q,W,2)=patcharraytemp(Q-increment+1,2)
         patcharray(Q,W,3)=patcharraytemp(Q-increment+1,3)
         patcharray(Q,W,4)=patcharraytemp(Q-increment+1,4)
         patcharray(Q,W,5)=patcharraytemp(Q-increment+1,5)
         patcharray(Q,W,6)=patcharraytemp(Q-increment+1,6)
         patcharray(Q,W,7)=0.0
         patcharray(Q,W,8)=0.0
         patcharray(Q,W,9)=1.0
         patcharray(Q,W,10)=5.0
 102     CONTINUE
C         print*,patcharray(Q,1,1),patcharray(Q,1,2),patcharray(Q,1,3)
C         print*,patcharraytemp(Q-increment+1,1),patcharraytemp
C     *        (Q-increment+1,2),patcharraytemp(Q-increment+1,3)
 101    CONTINUE
      increment=Q
      tempsize=Nx(2)*Ny(1)
      deallocate(ddx1)
      deallocate(ddy1)
      deallocate(patcharraytemp)

      allocate(ddx1(Nx(2)))
      allocate(ddy1(Ny(1)))
      allocate(patcharraytemp(Nx(2)*Ny(1),6))
      b=ylimit(2)
      b1=ylimit(1)
      CALL PATCHESSHORT(xlimit(2),xlimit(3),Nx(2),qx,ddx1)
      CALL PATCHESSHORT(ylimit(1),ylimit(2),Ny(1),qy,ddy1)
      CALL CREATEPATCHARRAY(ddx1,ddy1,Nx(2),Ny(1),xlimit(2),xlimit(3),
     *     ylimit(1),ylimit(2),zlimit(1),zlimit(1),patcharraytemp,slope,
     *     b,slope1,b1,count,FaceNormals(5,1:3))
      DO 103 Q=increment,increment+tempsize-1
         DO 104 W=1,sizeffttwo
         patcharray(Q,W,1)=patcharraytemp(Q-increment+1,1)
         patcharray(Q,W,2)=patcharraytemp(Q-increment+1,2)
         patcharray(Q,W,3)=patcharraytemp(Q-increment+1,3)
         patcharray(Q,W,4)=patcharraytemp(Q-increment+1,4)
         patcharray(Q,W,5)=patcharraytemp(Q-increment+1,5)
         patcharray(Q,W,6)=patcharraytemp(Q-increment+1,6)
         patcharray(Q,W,7)=0.0
         patcharray(Q,W,8)=0.0
         patcharray(Q,W,9)=2.0
         patcharray(Q,W,10)=5.0
 104     CONTINUE
 103  CONTINUE
      deallocate(ddx1)
      deallocate(ddy1)
      increment=Q
      tempsize=Nx(3)*Ny(1)
      deallocate(patcharraytemp)
      allocate(ddx1(Nx(3)))
      allocate(ddy1(Ny(1)))
      allocate(patcharraytemp(Nx(3)*Ny(1),6))
      b=ylimit(2)
      b1=ylimit(1)
      CALL PATCHESSHORT(xlimit(3),xlimit(4),Nx(3),qx,ddx1)
      CALL PATCHESSHORT(ylimit(1),ylimit(2),Ny(1),qy,ddy1)
      CALL CREATEPATCHARRAY(ddx1,ddy1,Nx(3),Ny(1),xlimit(3),xlimit(4),
     *     ylimit(1),ylimit(2),zlimit(1),zlimit(1),patcharraytemp,slope,
     *     b,slope1,b1,count,FaceNormals(5,1:3))
      DO 105 Q=increment,increment+tempsize-1
         DO 106 W=1,sizeffttwo
            patcharray(Q,W,1)=patcharraytemp(Q-increment+1,1)
            patcharray(Q,W,2)=patcharraytemp(Q-increment+1,2)
            patcharray(Q,W,3)=patcharraytemp(Q-increment+1,3)
            patcharray(Q,W,4)=patcharraytemp(Q-increment+1,4)
            patcharray(Q,W,5)=patcharraytemp(Q-increment+1,5)
            patcharray(Q,W,6)=patcharraytemp(Q-increment+1,6)
            patcharray(Q,W,7)=0.0
            patcharray(Q,W,8)=0.0
            patcharray(Q,W,9)=3.0
            patcharray(Q,W,10)=5.0
 106     CONTINUE
 105  CONTINUE
      deallocate(ddx1)
      deallocate(ddy1)
      increment=Q
      tempsize=Nx(1)*Ny(2)
      deallocate(patcharraytemp)
      allocate(ddx1(Nx(1)))
      allocate(ddy1(Ny(2)))
      allocate(patcharraytemp(Nx(1)*Ny(2),6))
      b=ylimit(3)
      b1=ylimit(2)
      CALL PATCHESSHORT(xlimit(1),xlimit(2),Nx(1),qx,ddx1)
      CALL PATCHESSHORT(ylimit(2),ylimit(3),Ny(2),qy,ddy1)
      CALL CREATEPATCHARRAY(ddx1,ddy1,Nx(1),Ny(2),xlimit(1),xlimit(2),
     *     ylimit(2),ylimit(3),zlimit(1),zlimit(1),patcharraytemp,slope,
     *     b,slope1,b1,count,FaceNormals(5,1:3))
      DO 107 Q=increment,increment+tempsize-1
         DO 108 W=1,sizeffttwo
         patcharray(Q,W,1)=patcharraytemp(Q-increment+1,1)
         patcharray(Q,W,2)=patcharraytemp(Q-increment+1,2)
         patcharray(Q,W,3)=patcharraytemp(Q-increment+1,3)
         patcharray(Q,W,4)=patcharraytemp(Q-increment+1,4)
         patcharray(Q,W,5)=patcharraytemp(Q-increment+1,5)
         patcharray(Q,W,6)=patcharraytemp(Q-increment+1,6)
         patcharray(Q,W,7)=0.0
         patcharray(Q,W,8)=0.0
         patcharray(Q,W,9)=4.0
         patcharray(Q,W,10)=5.0
 108     CONTINUE
 107  CONTINUE
      deallocate(ddx1)
      deallocate(ddy1)
      increment=Q
      tempsize=Nx(3)*Ny(2)
      deallocate(patcharraytemp)
      allocate(ddx1(Nx(3)))
      allocate(ddy1(Ny(2)))
      allocate(patcharraytemp(Nx(3)*Ny(2),6))
      b=ylimit(3)
      b1=ylimit(2)
      CALL PATCHESSHORT(xlimit(3),xlimit(4),Nx(3),qx,ddx1)
      CALL PATCHESSHORT(ylimit(2),ylimit(3),Ny(2),qy,ddy1)
      CALL CREATEPATCHARRAY(ddx1,ddy1,Nx(3),Ny(2),xlimit(3),xlimit(4)
     *     ,ylimit(2),ylimit(3),zlimit(1),0.0,patcharraytemp,slope,b,
     *     slope1,b1,count,FaceNormals(5,1:3))
      DO 109 Q=increment,increment+tempsize-1
         DO 110 W=1,sizeffttwo
         patcharray(Q,W,1)=patcharraytemp(Q-increment+1,1)
         patcharray(Q,W,2)=patcharraytemp(Q-increment+1,2)
         patcharray(Q,W,3)=patcharraytemp(Q-increment+1,3)
         patcharray(Q,W,4)=patcharraytemp(Q-increment+1,4)
         patcharray(Q,W,5)=patcharraytemp(Q-increment+1,5)
         patcharray(Q,W,6)=patcharraytemp(Q-increment+1,6)
         patcharray(Q,W,7)=0.0
         patcharray(Q,W,8)=0.0
         patcharray(Q,W,9)=5.0
         patcharray(Q,W,10)=5.0
 110  CONTINUE
 109  CONTINUE
      deallocate(ddx1)
      deallocate(ddy1)     
      increment=Q
      tempsize=Nx(1)*Ny(3)
      deallocate(patcharraytemp)
      allocate(patcharraytemp(Nx(1)*Ny(3),6))
      b=ylimit(4)
      b1=ylimit(3)
      allocate(ddx1(Nx(1)))
      allocate(ddy1(Ny(3)))
      CALL PATCHESSHORT(xlimit(1),xlimit(2),Nx(1),qx,ddx1)
      CALL PATCHESSHORT(ylimit(3),ylimit(4),Ny(3),qy,ddy1)
      CALL CREATEPATCHARRAY(ddx1,ddy1,Nx(1),Ny(3),xlimit(1),xlimit(2),
     *     ylimit(3),ylimit(4),zlimit(1),zlimit(1),patcharraytemp,slope
     *     ,b,slope1,b1,count,FaceNormals(5,1:3))
      DO 111 Q=increment,increment+tempsize-1
         DO 112 W=1,sizeffttwo
         patcharray(Q,W,1)=patcharraytemp(Q-increment+1,1)
         patcharray(Q,W,2)=patcharraytemp(Q-increment+1,2)
         patcharray(Q,W,3)=patcharraytemp(Q-increment+1,3)
         patcharray(Q,W,4)=patcharraytemp(Q-increment+1,4)
         patcharray(Q,W,5)=patcharraytemp(Q-increment+1,5)
         patcharray(Q,W,6)=patcharraytemp(Q-increment+1,6)
         patcharray(Q,W,7)=0.0
         patcharray(Q,W,8)=0.0
         patcharray(Q,W,9)=6.0
         patcharray(Q,W,10)=5.0
 112  CONTINUE
 111  CONTINUE
      deallocate(ddx1)
      deallocate(ddy1)
      increment=Q
      tempsize=Nx(2)*Ny(3)
      deallocate(patcharraytemp)
      allocate(ddx1(Nx(2)))
      allocate(ddy1(Ny(3)))
      allocate(patcharraytemp(Nx(2)*Ny(3),6))
      b=ylimit(4)
      b1=ylimit(3)
      CALL PATCHESSHORT(xlimit(2),xlimit(3),Nx(2),qx,ddx1)
      CALL PATCHESSHORT(ylimit(3),ylimit(4),Ny(3),qy,ddy1)
      CALL CREATEPATCHARRAY(ddx1,ddy1,Nx(2),Ny(3),xlimit(2),xlimit(3),
     *     ylimit(3),ylimit(4),zlimit(1),zlimit(1),patcharraytemp,slope,
     *     b,slope1,b1,count,FaceNormals(5,1:3))
      DO 113 Q=increment,increment+tempsize-1
         DO 114 W=1,sizeffttwo
         patcharray(Q,W,1)=patcharraytemp(Q-increment+1,1)
         patcharray(Q,W,2)=patcharraytemp(Q-increment+1,2)
         patcharray(Q,W,3)=patcharraytemp(Q-increment+1,3)
         patcharray(Q,W,4)=patcharraytemp(Q-increment+1,4)
         patcharray(Q,W,5)=patcharraytemp(Q-increment+1,5)
         patcharray(Q,W,6)=patcharraytemp(Q-increment+1,6)
         patcharray(Q,W,7)=0.0
         patcharray(Q,W,8)=0.0
         patcharray(Q,W,9)=7.0
         patcharray(Q,W,10)=5.0
 114  CONTINUE
 113  CONTINUE
      deallocate(ddx1)
      deallocate(ddy1)
      increment=Q
      tempsize=Nx(3)*Ny(3)
      deallocate(patcharraytemp)
      allocate(ddx1(Nx(3)))
      allocate(ddy1(Ny(3)))
      allocate(patcharraytemp(Nx(3)*Ny(3),6))
      b=ylimit(4)
      b1=ylimit(3)
      CALL PATCHESSHORT(xlimit(3),xlimit(4),Nx(3),qx,ddx1)
      CALL PATCHESSHORT(ylimit(3),ylimit(4),Ny(3),qy,ddy1)
      CALL CREATEPATCHARRAY(ddx1,ddy1,Nx(3),Ny(3),xlimit(3),xlimit(4),
     *     ylimit(3),ylimit(4),zlimit(1),zlimit(1),patcharraytemp,slope,
     *     b,slope1,b1,count,FaceNormals(5,1:3))
      DO 115 Q=increment,increment+tempsize-1
         DO 116 W=1,sizeffttwo
         patcharray(Q,W,1)=patcharraytemp(Q-increment+1,1)
         patcharray(Q,W,2)=patcharraytemp(Q-increment+1,2)
         patcharray(Q,W,3)=patcharraytemp(Q-increment+1,3)
         patcharray(Q,W,4)=patcharraytemp(Q-increment+1,4)
         patcharray(Q,W,5)=patcharraytemp(Q-increment+1,5)
         patcharray(Q,W,6)=patcharraytemp(Q-increment+1,6)
         patcharray(Q,W,7)=0.0
         patcharray(Q,W,8)=0.0
         patcharray(Q,W,9)=8.0
         patcharray(Q,W,10)=5.0
 116  CONTINUE
 115  CONTINUE
      deallocate(ddx1)
      deallocate(ddy1)
      increment=Q
      tempsize=Nx(2)*Nz(1)
      deallocate(patcharraytemp)
      allocate(ddx1(Nx(2)))
      allocate(ddz1(Nz(1)))
      allocate(patcharraytemp(Nx(2)*Nz(1),6))
      b=zlimit(2)
      b1=zlimit(1)
      CALL PATCHESSHORT(xlimit(2),xlimit(3),Nx(2),qx,ddx1)
      CALL PATCHESSHORT(zlimit(1),zlimit(2),Nz(1),qz,ddz1)
      CALL CREATEPATCHARRAY(ddx1,ddz1,Nx(2),Nz(1),xlimit(2),xlimit(3),
     *     ylimit(2),ylimit(2),zlimit(1),zlimit(2),patcharraytemp,slope
     *     ,b,slope1,b1,count,FaceNormals(4,1:3))
      DO 117 Q=increment,increment+tempsize-1
         DO 118 W=1,sizeffttwo
         patcharray(Q,W,1)=patcharraytemp(Q-increment+1,1)
         patcharray(Q,W,2)=patcharraytemp(Q-increment+1,2)
         patcharray(Q,W,3)=patcharraytemp(Q-increment+1,3)
         patcharray(Q,W,4)=patcharraytemp(Q-increment+1,4)
         patcharray(Q,W,5)=patcharraytemp(Q-increment+1,5)
         patcharray(Q,W,6)=patcharraytemp(Q-increment+1,6)
         patcharray(Q,W,7)=0.0
         patcharray(Q,W,8)=0.0
         patcharray(Q,W,9)=9.0
         patcharray(Q,W,10)=4.0
 118  CONTINUE
 117  CONTINUE
      deallocate(ddx1)
      deallocate(ddz1)
      increment=Q
      tempsize=Nx(2)*Nz(1)
      deallocate(patcharraytemp)
      allocate(ddx1(Nx(2)))
      allocate(ddz1(Nz(1)))
      allocate(patcharraytemp(Nx(2)*Nz(1),6))
      b=zlimit(2)
      b1=zlimit(1)
      CALL PATCHESSHORT(xlimit(2),xlimit(3),Nx(2),qx,ddx1)
      CALL PATCHESSHORT(zlimit(1),zlimit(2),Nz(1),qz,ddz1)
      CALL CREATEPATCHARRAY(ddx1,ddz1,Nx(2),Nz(1),xlimit(2),xlimit(3),
     *     ylimit(3),ylimit(3),zlimit(1),zlimit(2),patcharraytemp,slope,
     *     b,slope1,b1,count,FaceNormals(2,1:3))
      DO 119 Q=increment,increment+tempsize-1
         DO 120 W=1,sizeffttwo
         patcharray(Q,W,1)=patcharraytemp(Q-increment+1,1)
         patcharray(Q,W,2)=patcharraytemp(Q-increment+1,2)
         patcharray(Q,W,3)=patcharraytemp(Q-increment+1,3)
         patcharray(Q,W,4)=patcharraytemp(Q-increment+1,4)
         patcharray(Q,W,5)=patcharraytemp(Q-increment+1,5)
         patcharray(Q,W,6)=patcharraytemp(Q-increment+1,6)
         patcharray(Q,W,7)=0.0
         patcharray(Q,W,8)=0.0
         patcharray(Q,W,9)=10.0
         patcharray(Q,W,10)=2.0
 120  CONTINUE
 119  CONTINUE
      deallocate(ddx1)
      deallocate(ddz1)
      increment=Q
      tempsize=Ny(2)*Nz(1)
      deallocate(patcharraytemp)
      allocate(ddy1(Ny(2)))
      allocate(ddz1(Nz(1)))
      allocate(patcharraytemp(Ny(2)*Nz(1),6))
      b=zlimit(2)
      b1=zlimit(1)
      CALL PATCHESSHORT(ylimit(2),ylimit(3),Ny(2),qy,ddy1)
      CALL PATCHESSHORT(zlimit(1),zlimit(2),Nz(1),qz,ddz1)
      CALL CREATEPATCHARRAY(ddy1,ddz1,Ny(2),Nz(1),xlimit(2),xlimit(2),
     *     ylimit(2),ylimit(3),zlimit(1),zlimit(2),patcharraytemp,slope
     *     ,b,slope1,b1,count,FaceNormals(1,1:3))
      DO 121 Q=increment,increment+tempsize-1
         DO 122 W=1,sizeffttwo
         patcharray(Q,W,1)=patcharraytemp(Q-increment+1,1)
         patcharray(Q,W,2)=patcharraytemp(Q-increment+1,2)
         patcharray(Q,W,3)=patcharraytemp(Q-increment+1,3)
         patcharray(Q,W,4)=patcharraytemp(Q-increment+1,4)
         patcharray(Q,W,5)=patcharraytemp(Q-increment+1,5)
         patcharray(Q,W,6)=patcharraytemp(Q-increment+1,6)
         patcharray(Q,W,7)=0.0
         patcharray(Q,W,8)=0.0
         patcharray(Q,W,9)=11.0
         patcharray(Q,W,10)=1.0
 122  CONTINUE
 121  CONTINUE
      deallocate(ddy1)
      deallocate(ddz1)
      increment=Q
      tempsize=Ny(2)*Nz(1)
      deallocate(patcharraytemp)
      allocate(ddy1(Ny(2)))
      allocate(ddz1(Nz(1)))
      allocate(patcharraytemp(Ny(2)*Nz(1),6))
      b=zlimit(2)
      b1=zlimit(1)
      CALL PATCHESSHORT(ylimit(2),ylimit(3),Ny(2),qy,ddy1)
      CALL PATCHESSHORT(zlimit(1),zlimit(2),Nz(1),qz,ddz1)
      CALL CREATEPATCHARRAY(ddy1,ddz1,Ny(2),Nz(1),xlimit(3),xlimit(3),
     *     ylimit(2),ylimit(3),zlimit(1),zlimit(2),patcharraytemp,slope,
     *     b,slope1,b1,count,FaceNormals(3,1:3))
      DO 123 Q=increment,increment+tempsize-1
         DO 124 W=1,sizeffttwo
         patcharray(Q,W,1)=patcharraytemp(Q-increment+1,1)
         patcharray(Q,W,2)=patcharraytemp(Q-increment+1,2)
         patcharray(Q,W,3)=patcharraytemp(Q-increment+1,3)
         patcharray(Q,W,4)=patcharraytemp(Q-increment+1,4)
         patcharray(Q,W,5)=patcharraytemp(Q-increment+1,5)
         patcharray(Q,W,6)=patcharraytemp(Q-increment+1,6)
         patcharray(Q,W,7)=0.0
         patcharray(Q,W,8)=0.0
         patcharray(Q,W,9)=12.0
         patcharray(Q,W,10)=3.0
 124  CONTINUE
 123  CONTINUE
      deallocate(ddy1)
      deallocate(ddz1)
      increment=Q
      tempsize=Ny(2)*Nx(2)
      deallocate(patcharraytemp)
      allocate(ddx1(Nx(2)))
      allocate(ddy1(Ny(2)))
      allocate(patcharraytemp(Ny(2)*Nx(2),6))
      b=ylimit(3)
      b1=ylimit(2)
      CALL PATCHESSHORT(xlimit(2),xlimit(3),Nx(2),qx,ddx1)
      CALL PATCHESSHORT(ylimit(2),ylimit(3),Ny(2),qy,ddy1)
      CALL CREATEPATCHARRAY(ddx1,ddy1,Nx(2),Ny(2),xlimit(2),xlimit(3),
     *     ylimit(2),ylimit(3),zlimit(2),zlimit(2),patcharraytemp,slope,
     *     b,slope1,b1,count,FaceNormals(5,1:3))
      DO 125 Q=increment,increment+tempsize-1
         DO 126 W=1,sizeffttwo
         patcharray(Q,W,1)=patcharraytemp(Q-increment+1,1)
         patcharray(Q,W,2)=patcharraytemp(Q-increment+1,2)
         patcharray(Q,W,3)=patcharraytemp(Q-increment+1,3)
         patcharray(Q,W,4)=patcharraytemp(Q-increment+1,4)
         patcharray(Q,W,5)=patcharraytemp(Q-increment+1,5)
         patcharray(Q,W,6)=patcharraytemp(Q-increment+1,6)
         patcharray(Q,W,7)=0.0
         patcharray(Q,W,8)=0.0
         patcharray(Q,W,9)=13.0
         patcharray(Q,W,10)=5.0
 126  CONTINUE
 125  CONTINUE
      deallocate(ddx1)
      deallocate(ddy1)
C      DO 129 Q=1,PatchNo
C         print*, patcharray(Q,1,1), patcharray(Q,1,2),patcharray(Q,1,3)
C 129     CONTINUE
      deallocate(patcharraytemp)
      DO 127 Q=1,PatchNo
         DO 128 W=1, PatchNo
C            print*, patcharray(Q,1,9), patcharray(W,1,9)
C            print*,patcharray(Q,1,1),patcharray(Q,1,2),patcharray(Q,1,3)
C            print*,patcharray(W,1,1),patcharray(W,1,2),patcharray(W,1,3)
            if (patcharray(W,1,9).eq.1.0.or.patcharray(W,1,9).eq.2.0
     *           .or.patcharray(W,1,9).eq.3.0.or.patcharray(W,1,9).eq.
     *           4.0.or.patcharray(W,1,9).eq.5.0.or.patcharray(W,1,9)
     *           .eq.6.0.or.patcharray(W,1,9).eq.7.0.or.
     *           patcharray(W,1,9).eq.8.0)then
               formfactors(Q,W,2)=1.0
            endif
            if (patcharray(W,1,9).eq.9.0.or.patcharray(W,1,9).eq.10.0
     *           .or.patcharray(W,1,9).eq.11.0.or.patcharray(W,1,9).eq.
     *           12.0.or.patcharray(W,1,9).eq.13.0)
     *           then 
               formfactors(Q,W,2)=2.0
            endif
            if(patcharray(Q,1,9).eq.1.0.and.(patcharray(W,1,9).eq.
     *           8.0.or.patcharray(W,1,9).eq.9.0.or.patcharray(W,1,9)
     *           .eq.12.0.or.patcharray(W,1,9).eq.13.0.or.
     *           patcharray(W,1,9).eq.14.0.or.patcharray(W,1,9).eq.
     *           15.0))then
               CALL PERPFORMFACTOR(patcharray,PatchNo,sizeffttwo,
     *              formfactors,Q,W,PI,FaceNormals,FaceNormalNo)
            elseif(patcharray(Q,1,9).eq.2.0.and.(patcharray(W,1,9)
     *              .eq.12.0))then
               CALL PERPFORMFACTOR(patcharray,PatchNo,sizeffttwo,
     *              formfactors,Q,W,PI,FaceNormals,FaceNormalNo)
            elseif(patcharray(Q,1,9).eq.3.0.and.(patcharray(W,1,9)
     *              .eq.13.0.or.patcharray(W,1,9).eq.14.0))then
               CALL PERPFORMFACTOR(patcharray,PatchNo,sizeffttwo,
     *              formfactors,Q,W,PI,FaceNormals,FaceNormalNo)
            elseif(patcharray(Q,1,9).eq.4.0.and.(patcharray(W,1,9)
     *              .eq.15.0))then
               CALL PERPFORMFACTOR(patcharray,PatchNo,sizeffttwo,
     *              formfactors,Q,W,PI,FaceNormals,FaceNormalNo)
            elseif(patcharray(Q,1,9).eq.5.0.and.(patcharray(W,1,9)
     *              .eq.10.0.or.patcharray(W,1,9).eq.11.0.or.
     *              patcharray(W,1,9).eq.12.0.or.patcharray(W,1,9).eq.
     *              13.0.or.patcharray(W,1,9).eq.14.0.or.
     *              patcharray(W,1,9).eq.15.0))then
               CALL PERPFORMFACTOR(patcharray,PatchNo,sizeffttwo,
     *              formfactors,Q,W,PI,FaceNormals,FaceNormalNo)
            elseif(patcharray(Q,1,9).eq.8.0.and.(patcharray(W,1,9)
     *              .eq.1.0))then
               CALL PERPFORMFACTOR(patcharray,PatchNo,sizeffttwo,
     *              formfactors,Q,W,PI,FaceNormals,FaceNormalNo)
            elseif(patcharray(Q,1,9).eq.9.0.and.(patcharray(W,1,9)
     *              .eq.1.0))then
               CALL PERPFORMFACTOR(patcharray,PatchNo,sizeffttwo,
     *              formfactors,Q,W,PI,FaceNormals,FaceNormalNo)
            elseif(patcharray(Q,1,9).eq.10.0.and.(patcharray(W,1,9)
     *              .eq.5.0))then
               CALL PERPFORMFACTOR(patcharray,PatchNo,sizeffttwo,
     *              formfactors,Q,W,PI,FaceNormals,FaceNormalNo)
            elseif(patcharray(Q,1,9).eq.11.0.and.(patcharray(W,1,9)
     *              .eq.5.0))then
               CALL PERPFORMFACTOR(patcharray,PatchNo,sizeffttwo,
     *              formfactors,Q,W,PI,FaceNormals,FaceNormalNo)
            elseif(patcharray(Q,1,9).eq.12.0.and.(patcharray(W,1,9)
     *              .eq.1.0.or.patcharray(W,1,9).eq.2.0.or.
     *              patcharray(W,1,9).eq.5.0))then
               CALL PERPFORMFACTOR(patcharray,PatchNo,sizeffttwo,
     *              formfactors,Q,W,PI,FaceNormals,FaceNormalNo)
            elseif(patcharray(Q,1,9).eq.13.0.and.(patcharray(W,1,9)
     *              .eq.1.0.or.patcharray(W,1,9).eq.3.0.or.
     *              patcharray(W,1,9).eq.5.0.or.patcharray(W,1,9).eq.
     *              14.0))then
               CALL PERPFORMFACTOR(patcharray,PatchNo,sizeffttwo,
     *              formfactors,Q,W,PI,FaceNormals,FaceNormalNo)
            elseif(patcharray(Q,1,9).eq.14.0.and.(patcharray(W,1,9)
     *              .eq.1.0.or.patcharray(W,1,9).eq.3.0.or.
     *              patcharray(W,1,9).eq.5.0.or.patcharray(W,1,9).eq.
     *              13.0))then
               CALL PERPFORMFACTOR(patcharray,PatchNo,sizeffttwo,
     *              formfactors,Q,W,PI,FaceNormals,FaceNormalNo)
            elseif(patcharray(Q,1,9).eq.15.0.and.(patcharray(W,1,9)
     *              .eq.1.0.or.patcharray(W,1,9).eq.4.0.or.
     *              patcharray(W,1,9).eq.5.0))then
               CALL PERPFORMFACTOR(patcharray,PatchNo,sizeffttwo,
     *              formfactors,Q,W,PI,FaceNormals,FaceNormalNo)
            else              
               formfactors(Q,W,1)=0.0
               formfactors(Q,W,3)=0.0
            endif
C           print*,Patcharray(Q,1,9),patcharray(W,1,9),formfactors(Q,W,1)
 128     CONTINUE
 127  CONTINUE
      print*, 'Calculated form factors'

