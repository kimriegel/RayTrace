C $Revision: 1.6 $ : $Date: 2010/02/25 15:48:42 $ : $Author: kim

C***--------------------------------------------------------------***
C***------------------------ ABSORPTION --------------------------***
C***--------------------------------------------------------------***

      real FUNCTION ABSORPTION(ps,freq,hr,Temp)
C This function computes the air absorption for a given frequency, 
C ambient pressure, relative humidity and temperature.

C Define all variables and reference values
      real ps0, ps, freq, hr, Temp, T0, T01, F, FrN, FrO
      ps0=1.0
      hr=20.0
      T0=293.15
      T01=273.16
      F=freq/ps

C Compute all relevant parameters
      psat=ps0*10**(-6.8346*(T01/Temp)**1.261+4.6151)
      h=ps0*(hr/ps)*(psat/ps0)
      FrN=1/ps0*(T0/Temp)**(1/2)*(9+280*h*
     * exp(-4.17*((T0/Temp)**(1/3)-1)))
      FrO=1/ps0*(24+4.04*10**4*h*((.02+h)/(.391+h)));
      term1=0.01275*(exp(-2239.1/Temp)/(FrO+F**2/FrO));
      term2=0.1068*(exp(-3352/Temp)/(FrN+F**2/FrN));
      ABSORPTION=ps0*F**2*((1.84*10**(-11.0)*(Temp/T0)**(0.5)*ps0)+
     * (Temp/T0)**(-5.0/2.0)*(term1+term2))
      RETURN
      END

C***--------------------------------------------------------------***
C***---------------------- TIMERECONSTRUCT -----------------------***
C***--------------------------------------------------------------***

      SUBROUTINE TIMERECONSTRUCT(sizefft,timearray,arraysize,temparray,
     *     timetemparray)
      
C This Function computes the timesignal from a given fft.  It writes the
C time signal to an array.

C Define all variables

      INTEGER sizefft,D,arraysize,W
      DOUBLE COMPLEX tempfft(sizefft/2+1)
      real timearray(sizefft)
      COMPLEX XJ
      INTEGER*8 invplan
      double precision timesignal(sizefft)
      real temparray(arraysize,sizefft/2,6)
      real timetemparray(arraysize,sizefft,5)

      INCLUDE 'fftw3.f'
      XJ=(0,1)
 
      print*, 'timereconstruct has been called'
      DO 35 D=1, arraysize
         DO 36 W=1, sizefft
            timetemparray(D,W,1)=temparray(D,1,1)
            timetemparray(D,W,2)=temparray(D,1,2)
            timetemparray(D,W,3)=temparray(D,1,3)
            timetemparray(D,W,4)=timearray(W)
            timetemparray(D,W,5)=0.0
 36      CONTINUE
 35   CONTINUE
      print*, 'timetemparray has been initialized'
C Create the complex array to feed into the inverse fft function
      DO 37 D=1, arraysize
         if (temparray(D,1,5).eq.0.0)then
            DO 40 W=1,sizefft
               timetemparray(D,W,5)=0.0
 40         CONTINUE
         else
         DO 38 W=1,sizefft/2+1
            if (W.eq.1) then
               tempfft(W)=cmplx(0.0)
            else
            tempfft(W)=cmplx(abs(temparray(D,W-1,5))*exp(XJ*
     *           temparray(D,W-1,6)))
            endif
 38      CONTINUE
C            print*, tempfft
         print*, 'created temparray'
C use fftw to compute the inverse fft.
         call dfftw_plan_dft_c2r_1d(invplan,sizefft,tempfft,
     *        timesignal, FFTW_ESTIMATE)
C           print*, timesignal
         call dfftw_execute(invplan, tempfft, timesignal)
         call dfftw_destroy_plan(invplan)
         print*, 'created time signature'
         DO 39 W=1,sizefft         
            timetemparray(D,W,5)=timesignal(W)
 39      CONTINUE
         endif
 37   CONTINUE
      return
      end

C***--------------------------------------------------------------***
C***------------------------ RECEIVERHIT -------------------------***
C***--------------------------------------------------------------***

      SUBROUTINE RECEIVERHITFUNC(sizefft,outputarray,arraysize,
     *     temparray)
      
C This Function adds the pressures from a ray when it hits the receiver.

C Define all variables

      INTEGER sizefft,D,arraysize,W
      real outputarray(sizefft/2,7)
      real temparray(arraysize,sizefft/2,6)
      double complex temp1, temp2,temp3
      COMPLEX XJ
      XJ=(0,1)
      print*, 'everything seems to initiate'
C Add new pressures to existing pressures in temparray 
C First Look for the correct location.
      DO 40 D=1, arraysize
C         print*,'outputarray2', outputarray(1,2)
C         print*,'outputarray3', outputarray(1,3)
C         print*,'outputarray4', outputarray(1,4)
C         print*,'temparray1',temparray(D,1,1)
C         print*,'temparray2',temparray(D,1,2)
C         print*,'temparray3',temparray(D,1,3)
         if (outputarray(1,2).eq.temparray(D,1,1).and.
     *        outputarray(1,3).eq.temparray(D,1,2).and.
     *        outputarray(1,4).eq.temparray(D,1,3))then
            print*, 'first if statement passed'
C If the location is the same, loop through the frequency and add current
C values with new values.
            DO 41 W=1, sizefft/2
               temp1=cmplx(abs(temparray(D,W,5))*exp(XJ*
     *              temparray(D,W,6)))
C               print*,temp1
C               if (W.eq.1)print*, 'temp1 fine'
               temp2=cmplx(abs(outputarray(W,5))*exp(XJ*
     *              outputarray(W,6)))
C               if (W.ge.(sizefft/2)-20) then
C                  print*, temp2
C               endif
C               if (W.eq.1)print*, 'temp2 fine'
               temp3=temp1+temp2
C               if (W.eq.1)print*, 'temp3 fine'
               temparray(D,W,5)=abs(temp3)
C               if (W.eq.1)print*, 'temparray 5 fine'
               temparray(D,W,6)=ATAN2(imagpart(temp3),realpart(temp3))
C               if (W.eq.1)print*, 'temparray 6 fine'
C                  print*, temparray(1,W,5)
 41         CONTINUE
C            print*, outputarray(:,6)
         endif
 40   CONTINUE
      print*, 'got through end'
      return
      end

C***--------------------------------------------------------------***
C***--------------------------- HEADER ---------------------------***
C***--------------------------------------------------------------***

      INTEGER FUNCTION Header(fileid)

C     this function prints the header for the tecplot data

      integer fileid
      Write(fileid,*) 'TITLE = "Pressure at earlevel"'
      Write(fileid,*) 'VARIABLES = "X[m]" "Y[m]" "Z[m]" "P[Pa]"'
      Write(fileid,*) 'TEXT'
      Write(fileid,*) 'CS=FRAME'
      Write(fileid,*) 'X=71.9660948264,Y=82.9866270431'
      Write(fileid,*) 'C=BLACK'
      Write(fileid,*) 'S=LOCAL'
      Write(fileid,*) 'HU=POINT'
      Write(fileid,*) 'LS=1 AN=MIDCENTER'
      Write(fileid,*) 'BX=Filled BXM=60 LT=0.1 BXO=BLACK BXF=WHITE'
      Write(fileid,*) 'F=HELV'
      Write(fileid,*) 'H=20 A=0'
      Write(fileid,*) 'MFC=""'
      Write(fileid,*) 'CLIPPING=CLIPTOVIEWPORT'
      Write(fileid,*) 'T="Time = &(SOLUTIONTIME%4f)"'  
      Header=0
      return
      end

C***--------------------------------------------------------------***
C***------------------------ TIMEHEADER --------------------------***
C***--------------------------------------------------------------***
      
      INTEGER FUNCTION TimeHeader(fileid,time,sizex,sizey,sizez,
     *     planename)

C     this function prints the header between each time step.

      integer fileid,sizex,sizey,sizez
      real time
      CHARACTER*20 planename
      Write(fileid,*) 'ZONE',' T="',planename,'"'
      Write(fileid,*) 'STRANDID=1, SOLUTIONTIME=',time
      Write(fileid,*) 'I=',sizex,'J=',sizey,'K=',sizez,
     * 'ZONETYPE=Ordered'
      Write(fileid,*) 'DATAPACKING=POINT'
      Write(fileid,*) 'DT=(SINGLE SINGLE SINGLE SINGLE )'
      header=0
      return
      end

C***--------------------------------------------------------------***
C***-------------------------- GRID ------------------------------***
C***--------------------------------------------------------------***

      SUBROUTINE Grid(radius,A,B,C,D,xmin,ymin,zmin,xmax,ymax,zmax,
     * receiverarray, arraysize,sizex,sizey,sizez,step)
      
C     This function creates an equally spaced grid of size step 
C     apart for a receiver plane array.

      real radius,A,B,C,D,xmin,ymin,zmin,xmax,ymax,zmax,step
      integer i,j,sizex,sizey,sizez
      integer count,arraysize
      real receiverarray(arraysize,3),s
      s=step/radius
      if(xmin.eq.xmax)then
         count=1
         DO 1 i=1,int((zmax-zmin)/step),1
            Do 2 j=1,int((ymax-ymin)/(step)),1
               receiverarray(count,1)=(D-B*(ymin+(s*j+(1-s))*radius)-
     *             C*(zmin+(s*i+(1-s))*radius))/A
               receiverarray(count,2)=ymin+(s*j+(1-s))*radius
               receiverarray(count,3)=zmin+(s*i+(1-s))*radius
               count=count+1
 2          Continue      
 1       CONTINUE
         sizex=int((ymax-ymin)/(step))
         sizey=int((zmax-zmin)/step)
         sizez=1
      endif
      if(ymin.eq.ymax)then
         count=1
         DO 3 i=1,int((xmax-xmin)/(step)),1
            Do 4 j=1,int((zmax-zmin)/(step)),1
               receiverarray(count,1)=xmin+(s*i+(1-s))*radius
               receiverarray(count,2)=(D-A*(xmin+(s*i+(1-s))*radius)-
     *             C*(zmin+(s*j+(1-s))*radius))/B
               receiverarray(count,3)=zmin+(s*j+(1-s))*radius
               count=count+1
 4             Continue      
 3       CONTINUE
         sizex=int((zmax-zmin)/step)
         sizey=int((xmax-xmin)/(step))
         sizez=1
      endif
      if(zmin.eq.zmax)then
         count=1         
         DO 5 i=1,int((xmax-xmin)/(step)),1
            Do 6 j=1,int((ymax-ymin)/(step)),1          
               receiverarray(count,1)=xmin+(s*i+(1-s))*radius   
               receiverarray(count,2)=ymin+(s*j+(1-s))*radius
               receiverarray(count,3)=(D-A*(xmin+(s*i+(1-s))*radius)-
     *          B*(ymin+(s*j+(1-s))*radius))/C
               count=count+1
 6          Continue      
 5       CONTINUE
         sizex=int((ymax-ymin)/(step))
         sizey=int((xmax-xmin)/step)
         sizez=1
      endif
      Return
      END

      SUBROUTINE InitialGrid(radius,A,B,C,D,theta,phi,xmin,ymin,zmin,
     * xmax,ymax,zmax,receiverarray, arraysize,sizex,sizey,sizez)
      
C     This function creates an equally spaced grid of size step apart

      real xmin,ymin,zmin,xmax,ymax,zmax
      integer i,j,sizex,sizey,sizez
      integer count,arraysize
      real receiverarray(arraysize,3)
      real A,B,C,D
      real xspace, yspace,zspace
      real theta, phi,radius
      real lengthx, lengthy, lengthz
      real vectorx(3), vectory(3), vectorz(3), initial(3)
      yspace=radius*abs(cos(phi))
      zspace=radius*abs(sin(theta))
      if(xmin.eq.xmax)then
         count=1      
         DO 1 i=1,int((zmax-zmin)/zspace),1
            Do 2 j=1,int((ymax-ymin)/(yspace)),1
               receiverarray(count,1)=(D-B*(ymin+j*yspace)-
     *             C*(zmin+i*zspace))/A
               receiverarray(count,2)=ymin+j*yspace
               receiverarray(count,3)=zmin+i*zspace
               count=count+1
 2          Continue      
 1       CONTINUE 
         sizex=int((ymax-ymin)/(yspace))
         sizey=int((zmax-zmin)/zspace)
         sizez=1
      endif
      if(ymin.eq.ymax)then
         count=1
         DO 3 i=1,int((xmax-xmin)/(xspace)),1
            Do 4 j=1,int((zmax-zmin)/(zspace)),1
               receiverarray(count,1)=xmin+i*xspace
               receiverarray(count,2)=(D-A*(xmin+i*xspace)-
     *             C*(zmin+j*zspace))/B
               receiverarray(count,3)=zmin+j*zspace
               count=count+1
 4             Continue      
 3       CONTINUE
         sizex=int((zmax-zmin)/zspace)
         sizey=int((xmax-xmin)/(xspace))
         sizez=1
      endif
      if(zmin.eq.zmax)then
         count=1
         DO 5 i=1,int((xmax-xmin)/(xspace)),1
            Do 6 j=1,int((ymax-ymin)/(yspace)),1          
               receiverarray(count,1)=xmin+i*xspace  
               receiverarray(count,2)=ymin+j*yspace
               receiverarray(count,3)=(D-A*(xmin+i*xspace)-
     *          B*(ymin+j*yspace))/C
               count=count+1
 6          Continue   
 5       CONTINUE
         sizex=int((xmax-xmin)/(xspace))
         sizey=int((ymax-ymin)/yspace)
         sizez=1
      endif
      Return
      END

C***--------------------------------------------------------------***
C***------------------------ SPHERECHECK -------------------------***
C***--------------------------------------------------------------***

      SUBROUTINE SPHERECHECK(Sc,Sr2,F,veci,dx)

C     This function performs a check whether a ray hits a sphere.  If
C     it does hit the function returns the distance to the sphere
      real Sr2,dx,dx0,dx1
      real F(3),veci(3),Sc(3), OC(3), L2oc,tca, t2hc, dir(3)
      real A,B,C,D
      real HUGE
      HUGE=1000000.0
      OC(1)=Sc(1)-veci(1)
      OC(2)=Sc(2)-veci(2)
      OC(3)=Sc(3)-veci(3)
      L2OC=dot_product(OC,OC)
      tca=dot_product(OC,F)
      t2hc=Sr2-L2OC+tca**2
      if(L2oc.lt.Sr2) then
        dx=HUGE
      elseif(tca.lt.0.0)then
         dx=HUGE
      elseif(t2hc.lt.0.0)then
         dx=HUGE
      else
         dx=tca-sqrt(t2hc)
      endif
      RETURN
      END
C***--------------------------------------------------------------***
C***---------------------------- CROSS ---------------------------***
C***--------------------------------------------------------------***

      SUBROUTINE CROSS(A, B, normal)

C     This function calculates a cross product of A and B and returns
C     normal

      real A(3), B(3), normal(3)
      real length
      normal(1)=A(2)*B(3)-A(3)*B(2)
      normal(2)=A(3)*B(1)-A(1)*B(3)
      normal(3)=A(1)*B(2)-A(2)*B(1)
      length=sqrt(normal(1)**2.0+normal(2)**2+normal(3)**2)
      if (length.ne.0.0)then
      normal=normal/length
      endif
      end
C***--------------------------------------------------------------***
C***--------------------------- POLYGON --------------------------***
C***--------------------------------------------------------------***

      SUBROUTINE POLYGON(Vecip1,F,Q,size,Number,PointNumbers,PolyArray
     *     ,BuildingPoints,normal,FaceNormalNo,FaceNormals,dxbuilding,
     *     behind)
      INTEGER Q,size,Number,PointNumbers,P,NC,SH,NSH,behind,FaceNormalNo
      real Vecip1(3),normal(3),d,F(3),HUGE,odd
      real PolyArray(Number,size+1),t,Vd,V0,G(size,2)
      real BuildingPoints(PointNumbers,3),maximum
      real intersection(3),dxbuilding,tempA,tempB,tempC
      double precision FaceNormals(FaceNormalNo,3)
      HUGE=1000000.0
      NC=0
      behind=0
      normal(1)=FaceNormals(int(PolyArray(Q,1)),1)
      normal(2)=FaceNormals(int(PolyArray(Q,1)),2)
      normal(3)=FaceNormals(int(PolyArray(Q,1)),3)
      d=-dot_product(normal,BuildingPoints(int(PolyArray(Q,2)),1:3))
      Vd=dot_product(normal,F)
      if (Vd.ge.0.0)then
         dxbuilding=HUGE
         GOTO 7
      endif
      V0=-(dot_product(normal,Vecip1)+d)
      t=V0/Vd
      if(t.lt.0.0)then
         dxbuilding=HUGE
         behind=1
         GOTO 7
      endif
      intersection(1)=Vecip1(1)+F(1)*t
      intersection(2)=Vecip1(2)+F(2)*t
      intersection(3)=Vecip1(3)+F(3)*t
      maximum=max(abs(normal(1)),abs(normal(2)),abs(normal(3)))
      if(maximum.eq.abs(normal(1)))then
         DO 8 P=1,size
            G(P,1:2)=(/intersection(2)-BuildingPoints(int(PolyArray(Q,
     *           1+P)),2),intersection(3)-BuildingPoints(int(PolyArray(Q
     *           ,1+P)),3)/)
 8       CONTINUE
      elseif(maximum.eq.abs(normal(2)))then
         DO 9 P=1,size
            G(P,1:2)=(/intersection(1)-BuildingPoints(int(PolyArray(Q,
     *           1+P)),1),intersection(3)-BuildingPoints(int(PolyArray(Q
     *           ,1+P)),3)/)
 9       CONTINUE
      elseif(maximum.eq.abs(normal(3)))then
         DO 10 P=1,size
            G(P,1:2)=(/intersection(1)-BuildingPoints(int(PolyArray(Q,
     *           1+P)),1),intersection(2)-BuildingPoints(int(PolyArray(Q
     *           ,1+P)),2)/)
 10      CONTINUE
      endif
      DO 11 P=1,size
         if(P.eq.size)then
            if(G(P,2).lt.0.0)then
               SH=-1
            else
               SH=1
            endif
            if(G(1,2).lt.0.0)then
               NSH=-1
            else
               NSH=1
            endif
         else
            if(G(P,2).lt.0.0)then
               SH=-1
            else
               SH=1
            endif
            if(G(P+1,2).lt.0.0)then
               NSH=-1
            else
               NSH=1
            endif
         endif
         if(SH.ne.NSH)then
            if(P.eq.size)then
               if(G(P,1).gt.0.0.and.G(1,1).gt.0.0)then
                  NC=NC+1
               elseif(G(P,1).gt.0.0.or.G(1,1).gt.0.0)then
                  IF((G(P,1)-(G(P,2)*(G(P+1,1)-G(P,1))/(G(P+1,2)-G(P,2))
     *                 )).GT.0.0)THEN
                     NC=NC+1
                  endif
               endif
            else
               if(G(P,1).gt.0.0.and.G(P+1,1).gt.0.0)then
                  NC=NC+1
               elseif(G(P,1).gt.0.0.or.G(P+1,1).gt.0.0)then
                  IF((G(P,1)-(G(P,2)*(G(P+1,1)-G(P,1))/(G(P+1,2)-G(P,2))
     *                 )).GT.0.0)THEN
                     NC=NC+1
                  endif
               endif
            endif
         endif         
 11   CONTINUE
      odd=MOD(NC,2)
      if(odd.eq.0)then
         dxbuilding=HUGE
      elseif(odd.eq.1)then
         dxbuilding=t
      endif
 7    CONTINUE
      END
C***--------------------------------------------------------------***
C***------------------------- INSIDECHECK ------------------------***
C***--------------------------------------------------------------***

      SUBROUTINE INSIDECHECK(point,BuildingPoints,PointNumbers,F,
     *     TriangleArray,SquareArray,TriangleNumber,SquareNumber,
     *     TriangleSequence,SquareSequence,Triangles,Squares,
     *     PolyBuilding,P,inside,FaceNormalNo,FaceNormals)
      integer PointNumbers,TriangleNumber,SquareNumber,Triangles,Squares
      integer P,behind,inside,PolyBuilding,I,FaceNormalNo
      real point(3),F(3),normal(3),dxnear
      real BuildingPoints(PointNumbers,3)
      real TriangleArray(TriangleNumber,6)
      real SquareArray(SquareNumber,7)
      integer SquareSequence(PolyBuilding,Squares)
      Integer TriangleSequence(PolyBuilding,Triangles)
      double precision FaceNormals(FaceNormalNo,3)
      inside=1
      DO 9 I=1, Triangles
         call Polygon(point,F,TriangleSequence(P,I),3,TriangleNumber,
     *        PointNumbers,TriangleArray,BuildingPoints,normal,
     *        FaceNormalNo,FaceNormals,dxnear,behind)
         if(behind.eq.0)then
            inside=0
            GOTO 11
         endif
 9    CONTINUE
      Do 10 I=1,Squares
         call Polygon(point,F,SquareSequence(P,I),4,SquareNumber,
     *        PointNumbers,SquareArray,BuildingPoints,normal,
     *        FaceNormalNo,FaceNormals,dxnear,behind)
         if(behind.eq.0)then
            inside=0
            GOTO 11
         endif
 10   CONTINUE
 11   Continue
      END
      
C***--------------------------------------------------------------***
C***---------------------------- PLANE ---------------------------***
C***--------------------------------------------------------------***

      SUBROUTINE PLANE(Vecip1, B1, B2, planehit,nbox)

C     This function calculates the normal at the hitpoint of a box.

      real nbox(3),Point1(3)
      real Vecip1(3), B1(3), B2(3),Point3(3), Point2(3)
      integer planehit
      if (planehit.eq.1) then
         if(Vecip1(1).eq.B1(1)) then
            Point2=(/B1(1),B1(2),B2(3)/)  
            Point3=(/B1(1),B2(2),B1(3)/) 
            call CROSS((Point2-B1),(Point3-B1),nbox)
         elseif(Vecip1(1).eq.B2(1))then
            Point1=(/B2(1),B1(2),B1(3)/)
            Point2=(/B2(1),B1(2),B2(3)/)
            Point3=(/B2(1),B2(2),B1(3)/)
            call CROSS((Point3-Point1),(Point2-Point1),nbox)
         endif
      endif
      if (planehit.eq.2) then
         if(Vecip1(2).eq.B1(2)) then
            Point2=(/B2(1),B1(2),B1(3)/)  
            Point3=(/B1(1), B1(2), B2(3)/) 
            call CROSS((Point2-B1),(Point3-B1),nbox)
         elseif(Vecip1(2).eq.B2(2))then
            Point1=(/B1(1),B2(2),B1(1)/)
            Point2=(/B1(1),B2(2),B2(3)/)
            Point3=(/B2(1),B2(2),B1(3)/)
            call CROSS((Point2-Point1),(Point3-Point1),nbox)
         endif
      endif
      if(planehit.eq.3) then
         if(Vecip1(3).eq.B1(3)) then
            Point2=(/B2(1),B1(2),B1(3)/)  
            Point3=(/B1(1), B2(2), B1(3)/) 
            call CROSS((Point3-B1),(Point2-B1),nbox)
         elseif(Vecip1(3).eq.B2(3))then
            Point2=(/B1(1),B2(2),B2(3)/)
            Point3=(/B2(1),B1(2),B2(3)/)
            call CROSS((Point2-B2),(Point3-B2),nbox)
         endif
      endif
      end

C***--------------------------------------------------------------***
C***--------------------------- BOX ------------------------------***
C***--------------------------------------------------------------***

      SUBROUTINE BOX(B1,B2,Vecip1,F,dxnear, dxfar, hit,planehit)

C     This function checks to see if the ray hits a box.  It determines which
C     plane the ray hits

      real T1X, T2X, T1Y, T2Y, T1Z, T2Z, tmp
      real B1(3), B2(3),tempF(3)
      real Vecip1(3), F(3), HUGE,dxnear, dxfar
      INTEGER hit, planehit, tmphit
      hit=5
      HUGE=1000000.0
      dxnear=-HUGE
      dxfar=HUGE
      tempF=F
C      print*, tempF
      if ((F(1).eq.0.0) .or. (F(2).eq.0.0) .or. (F(3).eq.0.0)) then
         if (F(1).eq.0.0) then
            if((vecip1(1).lt.B1(1)) .or. (vecip1(1).gt.B2(1)))then
               hit=0
               dxnear=HUGE
               GO TO 100
            else
            endif
         endif
         if (F(2).eq.0.0) then
            if((vecip1(2).lt.B1(2)) .or. (vecip1(2).gt.B2(2)))then
               hit=0
               dxnear=HUGE
               GO To 100
            endif
         endif
         if (F(3).eq.0.0) then
C            print*, 'Vecip1 inside',vecip1(3)
            if((vecip1(3).lt.B1(3)) .or. (vecip1(3).gt.B2(3)))then
C               print*, 'does this happen?'
               hit=0
               dxnear=HUGE
               Go to 100
            Endif
         endif 
      endif
C      print*, 'hit',hit
      if (hit.ne.0) then
         if(F(1).eq.0)tempF(1)=1.0
         if(F(2).eq.0)tempF(2)=1.0
         if(F(3).eq.0)tempF(3)=1.0
C     print*, B1(1),B2(1), Vecip1(1), F(1)
         if(F(1).ne.0.0)then
            T1X=(B1(1)-Vecip1(1))/tempF(1)
            T2X=(B2(1)-Vecip1(1))/tempF(1)
C      print*, 'T1X,T2x',T1X,T2X
C     print*, T1X, T2X
            If (T1X.gt.T2X) then
               tmp=T1X
               T1X=T2X
               T2X=tmp
            endif
            if(T1X.gt.dxnear)dxnear=T1X
            if(T2X.lt.dxfar)dxfar=T2X
            if (dxnear.GT.dxfar) then
               hit=0
               dxnear=huge
               go to 100
            elseif (dxfar.lt.0.0) then
               hit=0
               dxnear=huge
               goto 100
            endif
         endif
         if(F(2).ne.0.0)then
            T1Y=(B1(2)-Vecip1(2))/tempF(2)
            T2Y=(B2(2)-Vecip1(2))/tempF(2)
C            print*, 'T1Z, T2Z',T1Y,T2Y
            If (T1Y.GT.T2Y) then
               tmp=T1Y
               T1Y=T2Y
               T2Y=tmp
            endif
            if (T1Y.GT.dxnear) dxnear=T1Y
            if (T2Y.LT.dxfar) dxfar=T2Y
            if (dxnear.GT.dxfar) then
               hit=0
               dxnear=huge
               goto 100
            elseif (dxfar.LT.0.0) then
               hit=0
               dxnear=huge
               goto 100
            endif
         endif
C         print*, B1(3), Vecip1(3), tempF(3)
C         print*, B2(3), Vecip1(3), tempF(3)
         if(F(3).ne.0.0)then
            T1Z=(B1(3)-Vecip1(3))/tempF(3)
            T2Z=(B2(3)-Vecip1(3))/tempF(3)
C       print*, 'T1Z,T2Z',T1Z,T2Z
            If (T1Z.GT.T2Z) then
               tmp=T1Z
               T1Z=T2Z
               T2Z=tmp
            endif
            if (T1Z.GT.dxnear) dxnear=T1Z
            if (T2Z.LT.dxfar) dxfar=T2Z
            if (dxnear.GT.dxfar) then
               hit=0
               dxnear=huge
               goto 100
            elseif (dxfar.LT.0) then
               hit=0
               dxnear=huge
               goto 100
            elseif (dxnear.le.0)then
               hit=0
               dxnear=HUGE
               goto 100
            endif
         endif
         if (hit.ne.0) then 
            if (dxnear.LT.dxfar) then
               hit=1
               if (dxnear.EQ.T1X) planehit=1
               if (dxnear.EQ.T1Y) planehit=2
               if (dxnear.EQ.T1Z) planehit=3
            endif
         endif
      endif
 100  CONTINUE
      return 
      end

      SUBROUTINE BOX2(B1,B2,Vecip1,F,dxnear, dxfar, hit,planehit)

C     This function checks to see if the ray hits a box.  It determines which
C     plane the ray hits

      real T1X, T2X, T1Y, T2Y, T1Z, T2Z, tmp
      real B1(3), B2(3),tempF(3)
      real Vecip1(3), F(3), HUGE,dxnear, dxfar
      INTEGER hit, planehit, tmphit
      hit=5
      HUGE=1000000.0
      dxnear=-HUGE
      dxfar=HUGE
      tempF=F
      if((F(1).eq.0.0) .or. (F(2).eq.0.0) .or. (F(3).eq.0.0)) then
         if (F(1).eq.0.0) then
            if((vecip1(1).lt.B1(1)) .or. (vecip1(1).gt.B2(1)))then
               hit=0
               dxnear=HUGE
               GO TO 100
            endif
         endif
         if (F(2).eq.0.0) then
            if((vecip1(2).lt.B1(2)) .or. (vecip1(2).gt.B2(2)))then
               hit=0
               dxnear=HUGE
               GO To 100
            endif
         endif
         if (F(3).eq.0.0) then
            if((vecip1(3).lt.B1(3)) .or. (vecip1(3).gt.B2(3)))then
               hit=0
               dxnear=HUGE
               Go to 100
            endif
         endif 
      endif
C      print*, 'hit',hit
      if (hit.ne.0) then
C         print*, B1(1),B2(1), Vecip1(1), F(1)
C         if(F(1).ne.0.0)then
            T1X=(B1(1)-Vecip1(1))/tempF(1)
            T2X=(B2(1)-Vecip1(1))/tempF(1)
C            print*, 'T1X,T2x',T1X,T2X
C     print*, T1X, T2X
            If (T1X.gt.T2X) then
               tmp=T1X
               T1X=T2X
               T2X=tmp
            endif
            if(T1X.gt.dxnear)dxnear=T1X
            if(T2X.lt.dxfar)dxfar=T2X
            if (dxnear.GT.dxfar) then
               hit=0
               dxnear=huge
               go to 100
            elseif (dxfar.lt.0.0) then
               hit=0
               dxnear=huge
               goto 100
            endif
C         endif
C         if(F(2).ne.0.0)then
            T1Y=(B1(2)-Vecip1(2))/tempF(2)
            T2Y=(B2(2)-Vecip1(2))/tempF(2)
C            print*, 'T1Y,T2Y',T1Y,T2Y
            If (T1Y.GT.T2Y) then
               tmp=T1Y
               T1Y=T2Y
               T2Y=tmp
            endif
            if (T1Y.GT.dxnear) dxnear=T1Y
            if (T2Y.LT.dxfar) dxfar=T2Y
            if (dxnear.GT.dxfar) then
               hit=0
               dxnear=huge
               goto 100
            elseif (dxfar.LT.0.0) then
               hit=0
               dxnear=huge
               goto 100
            endif
C         endif
C         if(F(3).ne.0.0)then
C     print*, B1(3), Vecip1(3), tempF(3)
C     print*, B2(3), Vecip1(3), tempF(3)
            T1Z=(B1(3)-Vecip1(3))/tempF(3)
            T2Z=(B2(3)-Vecip1(3))/tempF(3)
C            print*, 'T1Z,T2Z',T1Z,T2Z
            If (T1Z.GT.T2Z) then
               tmp=T1Z
               T1Z=T2Z
               T2Z=tmp
            endif
            if (T1Z.GT.dxnear) dxnear=T1Z
            if (T2Z.LT.dxfar) dxfar=T2Z
            if (dxnear.GT.dxfar) then
               hit=0
               dxnear=huge
               goto 100
            elseif (dxfar.LT.0) then
               hit=0
               dxnear=huge
               goto 100
            elseif (dxnear.lt.-0.005)then
               hit=0
               dxnear=HUGE
               goto 100
            endif
C         endif  

         if (hit.ne.0) then 
            if (dxnear.LT.dxfar) then
               hit=1
               if (dxnear.EQ.T1X) planehit=1
               if (dxnear.EQ.T1Y) planehit=2
               if (dxnear.EQ.T1Z) planehit=3
            endif
         endif
      endif
 100  CONTINUE
      return 
      end

      SUBROUTINE ROTATION(axis, angle, rotationmatrix)

      real axis(3)
      real angle
      real rotationmatrix(3,3)
      rotationmatrix(1,1)=axis(1)**2+(1-axis(1)**2)*cos(angle)
      rotationmatrix(1,2)=axis(1)*axis(2)*(1-cos(angle))+axis(3)*
     *     sin(angle)
      rotationmatrix(1,3)=axis(1)*axis(3)*(1-cos(angle))-axis(2)*
     *     sin(angle)                
      rotationmatrix(2,1)=axis(1)*axis(2)*(1-cos(angle))-axis(3)*
     *     sin(angle)                 
      rotationmatrix(2,2)=axis(2)**2+(1-axis(2)**2)*cos(angle)              
      rotationmatrix(2,3)=axis(2)*axis(3)*(1-cos(angle))+axis(1)*
     *     sin(angle)
      rotationmatrix(3,1)=axis(1)*axis(3)*(1-cos(angle))+axis(2)*
     *     sin(angle)
      rotationmatrix(3,2)=axis(2)*axis(3)*(1-cos(angle))-axis(1)*
     *     sin(angle)
      rotationmatrix(3,3)=axis(3)**2+(1-axis(3)**2)*cos(angle)

      RETURN
      END
