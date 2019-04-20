C $Revision: 1.0 $ : $Date:$ : $Author: kim

C***------------------------------------------------------------***
C***-------------------------- PATCHES -------------------------***
C***------------------------------------------------------------***

      SUBROUTINE PATCHESSHORT(min,max,N,q,dd)

C     This function creates patch lengths using a geometric sum.

      integer N,m
      double precision dd(N), min, max, W,k(N),q
      W=max-min
C     For box 1 we mesh the x direction
      DO 1 m=1, N
         k(m)=W/2.0*(1-q)/(1-q**(N/2))
         if (m.le.(N/2))then
            dd(m)=k(m)*q**(m-1)
         elseif(m.le.N.and.m.gt.(N/2))then
            dd(m)=k(m)*q**(N-m)
         endif
 1    CONTINUE

      END

      SUBROUTINE CREATEPATCHARRAY(ddm,ddl,Nm,Nl,x1,x2
     *     ,y1,y2,z1,z2,patcharray,slope,b,slope1,b1,count,normal)

C     This function creates an array of patches for each plane

      integer Nm,Nl,count
      double precision x1,y1,x2,y2,z1,z2,normal(3)
      real temp1(3),temp2(3),temp3(3),vec1(3),vec2(3)
      real vec3(3),slope,b,ddz,ddy,ddx,zcenter,xcenter,ycenter
      real slope1,b1
      double precision patcharray(Nm*Nl,6),x,y,z,ddm(Nm),ddl(Nl),d
      d=-x1*normal(1)-y1*normal(2)-z1*normal(3)
      if (z1.eq.z2)then
         count=1
         print*, 'ZPLANE',Nl, Nm, x1, y1, z1
         DO 2 l=1, Nl
            DO 3 m=1, Nm
               x=x1-0.5*ddm(m)+SUM(ddm(1:m))
               y=y1-0.5*ddl(l)+SUM(ddl(1:l))
               z=(-d-normal(1)*x-normal(2)*y)/normal(3)
               if(m.eq.1)then
                  zcenter=z
               endif
               if(l.eq.1)then
                  ddz=z-z1
               else
                  ddz=0.5*(z-zcenter)
               endif
               if(y.gt.slope*x+b.or.y.lt.slope1*x+b1)then
                  GOTO 3
               endif

               patcharray(count,1)=x 
               patcharray(count,2)=y
               patcharray(count,3)=z
               patcharray(count,4)=ddm(m)
               patcharray(count,5)=ddl(l)
               patcharray(count,6)=ddz
               count=count+1
              
 3          CONTINUE
 2       CONTINUE
      elseif(y1.eq.y2)then
         count=1
         print*, 'YPLANE',Nl, Nm, x1, y1, z1
         DO 4 l=1, Nl
            DO 5 m=1, Nm
               x=x1-0.5*ddm(m)+SUM(ddm(1:m))
               z=z1-.5*ddl(l)+SUM(ddl(1:l))
               y=(-d-normal(1)*x-normal(3)*z)/normal(2)
               if(m.eq.1)then
                  ycenter=y
               endif
               if(l.eq.1)then
                  ddy=y-y1
               else
                  ddy=0.5*(y-ycenter)
               endif
               if(z.gt.slope*x+b.or.z.lt.slope1*x+b1)then
                  GOTO 5
               endif
               patcharray(count,1)=x
               patcharray(count,2)=y
               patcharray(count,3)=z
               patcharray(count,4)=ddm(m)
               patcharray(count,5)=ddy
               patcharray(count,6)=ddl(l)

               count=count+1
 5          CONTINUE
 4       CONTINUE
      elseif(x1.eq.x2)then
         count=1
         print*, 'XPLANE',Nl, Nm, x1, y1, z1
         DO 6 l=1, Nl
            DO 7 m=1, Nm

               y=y1-.5*ddm(m)+SUM(ddm(1:m))
               z=z1-.5*ddl(l)+SUM(ddl(1:l))
               x=(-d-normal(3)*z-normal(2)*y)/normal(1)
               if(m.eq.1)then
                  xcenter=x
               endif
               if(l.eq.1)then
                  ddx=x-x1
               else
                  ddx=0.5*(x-xcenter)
               endif
               if(z.gt.slope*y+b.or.z.lt.slope1*y+b1)then
                  GOTO 7
               endif
               patcharray(count,1)=x
               patcharray(count,2)=y
               patcharray(count,3)=z
               patcharray(count,4)=ddx
               patcharray(count,5)=ddm(m)
               patcharray(count,6)=ddl(l)
               count=count+1
 7          CONTINUE
 6       CONTINUE
      endif
      END
      
      SUBROUTINE PERPFORMFACTOR(patcharray,PatchNo,sizeffttwo,
     *     formfactors,Q,W,PI,FaceNormals,FaceNormalNo)
      
C     This function creates form factors between two patches

      Integer PatchNo, sizeffttwo,Q,W,I,FaceNormalNo
      real patcharray(PatchNo,sizeffttwo,10)
      double precision FaceNormals(FaceNormalNo,3)
      real S1(3),costheta1,costheta2,length1,length2,vec1(3),vec2(3)
      real S2(3),S1length,S2length,area,minimum
      real formfactors(PatchNo,PatchNo,3),dlnlm,nu
      real PI
      dlnlm=sqrt((patcharray(Q,1,1)-patcharray(W,1,1))**2+
     *           (patcharray(Q,1,2)-patcharray(W,1,2))**2+
     *           (patcharray(Q,1,3)-patcharray(W,1,3))**2)
      formfactors(Q,W,3)=dlnlm
      vec1=FaceNormals(int(patcharray(Q,1,10)),1:3)
      vec2=FaceNormals(int(patcharray(W,1,10)),1:3)
      length1=sqrt(vec1(1)**2+vec1(2)**2+vec1(3)**2)
      length2=sqrt(vec2(1)**2+vec2(2)**2+vec2(3)**2)
      S1(1)=patcharray(W,1,1)-patcharray(Q,1,1)
      S1(2)=patcharray(W,1,2)-patcharray(Q,1,2)
      S1(3)=patcharray(W,1,3)-patcharray(Q,1,3)
      S1length=sqrt(S1(1)**2+S1(2)**2+S1(3)**2)
      S2(1)=patcharray(Q,1,1)-patcharray(W,1,1)
      S2(2)=patcharray(Q,1,2)-patcharray(W,1,2)
      S2(3)=patcharray(Q,1,3)-patcharray(W,1,3)
      S2length=sqrt(S2(1)**2+S2(2)**2+S2(3)**2)
      costheta1=(vec1(1)*S1(1)+vec1(2)*S1(2)+vec1(3)*S1(3))/(length1*
     *     S1length)
      costheta2=(vec2(1)*S2(1)+vec2(2)*S2(2)+vec2(3)*S2(3))/(length2*
     *     S2length)
      minimum=min(patcharray(W,1,4),patcharray(W,1,5),patcharray(W,1,6))
      if(minimum.eq.patcharray(W,1,4))then
         area=patcharray(W,1,5)*patcharray(W,1,6)
      elseif(minimum.eq.patcharray(W,1,5))then
         area=patcharray(W,1,4)*patcharray(W,1,6)
      elseif(minimum.eq.patcharray(W,1,6))then
         area=patcharray(W,1,5)*patcharray(W,1,4)
      endif
      formfactors(Q,W,1)=(costheta1*costheta2)*area/(PI*S1length**2)
      END
