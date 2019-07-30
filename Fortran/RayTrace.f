      PROGRAM RayTrace
C     Kimberly Lefkowitz created this program to propagate sonic booms around
C     large structures, and to graduate.  It is a ray tracing model that
C     will include specular and diffuse reflections.  It will print out the 
C     sound field at ear height, at relevent microphone locations, and at
C     the building walls.  It will read in the fft of a sonic boom signiture.
 
C $Revision: 1.10 $ : $Date: 2010/02/25 15:48:12 $ : $Author: kim
C     Initialize all Variables

C     Initialize all parameter variables

      INTEGER sum,absorbplanes
      INTEGER tmp4
      real soundspeed, lambda,freq,length,ps, Temp, hr
      real Fs, tempalphaground(8)
      real PI,HUGE,m,percentdiffuse,tmpsum
      PARAMETER (PI=3.14159265358979323846,HUGE=1000000.0)
      COMPLEX XJ

C     Initialize all geometry parameters

      real slope, b,slope1,b1
      INTEGER boxnumber,TriangleNumber,SquareNumber,PointNumbers
      INteGEr Triangles, Squares, PolyBuilding,FaceNormalNo,behind

C     Initialize phase propagation variables

      real phasefinal

C     Initialize variables for complex absorptions

      integer complexabsorption
      real height1, height2, height3

C     Initialize Amplitude propagation variables
      real ampfinal,normalization,normal(3)

C     Initialize Distance and Direction propagation variables

      real Vinitial(3),Vecip1(3),veci(3),r(3),F(3),h
      real Finitial(3),dx

C     Initialize iterators and counters

      INTEGER increment
      INTEGER IMAX,I,Q,P,ray,RAYMAX,W,S,K,count,j,D

C     Initialize receiver variables

      real radius,radius2,receiverpoint(3),receiverpoint2(3)
      real xmin,ymin,zmin,xmax,ymax,zmax
      real dxreceiver,tempreceiver,receivercheck
      real receiverA,receiverB,receiverC,receiverD
      real lastreceiver(3)
      real lastreceiver2(3),checkdirection(3),OC(3),OCLength
      Integer receiverhit,doublehit,hitcount,arraysize,arraysize1
      integer arraysize2,arraysize3,arraysize4,arraysize5,arraysize6
      Integer arraysize7,planenum
      real xspace, yspace, zspace

C     Initialize Ground Variables

      real nground(3),groundheight,GROUNDABC(3),GROUNDD
      double precision GROUNDVD,GROUNDVO,dxground1, GROUNDN(3)
      double precision Ftemp(3),vecitemp(3),vecip1temp(3)
      real dxground
      integer groundhit

C     Initialize intial array variables

      real PLANEABC(4),yinitial, xinitial,zinitial,area,corr
      real theta,phi,ninitial,xiinitial,zetainitial
      real boomspacing

C     Initialize initial signal variables

      INTEGER sizefft
      INTEGER*8 plan

C     Initialize Misc Variables

      real temp1,tmp1,tmp2,tmp3,tmp,n2,twopi,twopih
      integer sizeffttwo
      real twopidx,PIRlm2,dot1
      double complex temp2, temp3, temp4,dot
      CHARACTER*20 FILENAME
      CHARACTER*30 INPUTFILE, OUTPUTFILE

C     Intilize Building Variables

      real nbuilding(3),nbox(3), dxbuilding, dxnear, dxfar
      integer buildinghit,hit,planehit,whichbox

C     Initialize output variables

      real timestep,time1,time,timelength, time2
      integer sizex,sizey,sizez
      integer sizex1,sizey1,sizez1,sizex2,sizey2,sizez2
      integer sizex3,sizey3,sizez3,sizex4,sizey4,sizez4
      integer sizex5,sizey5,sizez5,sizex6,sizey6,sizez6
      integer sizex7,sizey7,sizez7
      character*20 planename1, planename2,planename3, planename4
      character*20 planename5, planename6,planename7

C     Initialize Radiosity Variables

      integer radiosity, PatchNo, tempsize,KMAX,Npatch
      real diffusion,cosgamma,cosdeltagamma,nu,Rlm,patcharea
      real cosxilm,Patchlength,RadFinitial(3),radF(3),diffusionground
      double complex Ek,vec3(3)
      double precision qx,qz,qy
      real minx, maxx, miny, maxy, minz, maxz
      integer PatchNox,PatchNoy,PatchNoz
      real qx1,qx2,qx3,qx4,qy1,qy2,qy3,qy4,qz1,dlnlm
      real dl, dm, dlprime, dnprime,ddm,knu,alpha

C     Initialize all dynamic Arrays

      integer, allocatable::Nx(:)
      integer, allocatable::Ny(:)
      integer, allocatable::Nz(:)
      real, allocatable::airabsorb(:)
      double precision, allocatable::xlimit(:)
      double precision, allocatable::ylimit(:)
      double precision, allocatable::zlimit(:)
      real, allocatable::ampinitial(:)
      real, allocatable::phaseinitial(:)
      real, allocatable::alphabuilding(:,:)
      real, allocatable::timearray(:)
      real, allocatable::inputarray(:,:)
      double complex, allocatable::outputsignal(:)
      double precision, allocatable::inputsignal(:)
      real, allocatable::outputarray1(:,:)
      real, allocatable::dhoutputarray1(:,:)
      real, allocatable::boomarray(:,:)
      real, allocatable::receiverarray(:,:)
      real, allocatable::receiverarray1(:,:)
      real, allocatable::receiverarray2(:,:)
      real, allocatable::receiverarray3(:,:)
      real, allocatable::receiverarray4(:,:)
      real, allocatable::receiverarray5(:,:)
      real, allocatable::receiverarray6(:,:)
      real, allocatable::receiverarray7(:,:)
      real, allocatable::temparray(:,:,:)
      real, allocatable::timetemparray(:,:,:)
      double precision, allocatable::ddx1(:)
      double precision, allocatable::ddy1(:)
      double precision, allocatable::ddz1(:)
      real, allocatable::patcharray(:,:,:)
      real, allocatable::patcharray1(:,:,:)
      double precision, allocatable::patcharraytemp(:,:)
      real, allocatable::formfactors(:,:,:)
      double complex, allocatable::Gk(:,:)
      double complex, allocatable::Gkminus1(:,:)
      real, allocatable::alphaground(:)
      real, allocatable::alphanothing(:)
      real, allocatable::boxarraynear(:,:)
      real, allocatable::boxarrayfar(:,:)
      real, allocatable::BuildingPoints(:,:)
      real, allocatable::TriangleArray(:,:)
      real, allocatable::SquareArray(:,:)
      Integer, allocatable::TriangleSequence(:,:)
      Integer, allocatable::SquareSequence(:,:)
      double precision, allocatable::FaceNormals(:,:)
      real, allocatable::tempalphabuilding(:,:)

C     Assign initial values to variables.  These will change depending on 
C     circumstances.
C     Must initialize a few constants in order to utilize fftw

      INCLUDE 'fftw3.f'
C     Include Parameter file 
      INCLUDE 'Parameterfile.f'

C     Include Structure Geometry

C      INCLUDE 'NoGeometry.f'
      INCLUDE 'BuildingGeometry.f'
C      INCLUDE 'UrbanCanyonH3L35.f'
C      INCLUDE 'UrbanCanyonH3L35Deg90.f'
C      INCLUDE 'UrbanCanyonH3L7.f'
C      INCLUDE 'UrbanCanyonH3L7Deg90.f'
C      INCLUDE 'UrbanCanyonH3L14.f'
C      INCLUDE 'UrbanCanyonH3L14Deg90.f'
C      INCLUDE 'UrbanCanyonH6L35.f'
C      INCLUDE 'UrbanCanyonH6L35Deg90.f'
C      INCLUDE 'UrbanCanyonH6L7.f'
C      INCLUDE 'UrbanCanyonH6L7Deg90.f'
C      INCLUDE 'UrbanCanyonH6L14.f'
C      INCLUDE 'UrbanCanyonH6L14Deg90.f'
C      INCLUDE 'UrbanCanyonH12L35.f'
C      INCLUDE 'UrbanCanyonH12L35Deg90.f'
C      INCLUDE 'UrbanCanyonH12L7.f'
C      INCLUDE 'UrbanCanyonH12L7Deg90.f'
C      INCLUDE 'UrbanCanyonH12L14.f'
C      INCLUDE 'UrbanCanyonH12L14Deg90.f'
C      INCLUDE 'UrbanCanyonH24L35.f'
C      INCLUDE 'UrbanCanyonH24L35Deg90.f'
C      INCLUDE 'UrbanCanyonH24L7.f'
C      INCLUDE 'UrbanCanyonH24L7Deg90.f'
C      INCLUDE 'UrbanCanyonH24L14.f'
C      INCLUDE 'UrbanCanyonH24L14Deg90.f'
C      INCLUDE 'BuildingGeometryVal.f'
C      INCLUDE 'BuildingGeometryMuse.f'
C      INCLUDE 'NASABuildingGeometry.f'
C      INCLUDE 'BuildingGeoMuseComp.f'
C      INCLUDE 'UrbanCanyonG2H3L55.f'
C      INCLUDE 'UrbanCanyonG2H3L7.f'
C      INCLUDE 'UrbanCanyonG2H3L14.f'
C      INCLUDE 'UrbanCanyonG2H6L55.f'
C      INCLUDE 'UrbanCanyonG2H6L7.f'
C      INCLUDE 'UrbanCanyonG2H6L14.f'
C      INCLUDE 'UrbanCanyonG2H12L55.f'
C      INCLUDE 'UrbanCanyonG2H12L7.f'
C      INCLUDE 'UrbanCanyonG2H12L14.f'
C      INCLUDE 'UrbanCanyonG2H24L55.f'
C      INCLUDE 'UrbanCanyonG2H24L7.f'
C      INCLUDE 'UrbanCanyonG2H24L14.f'

C     Initialize counters and calculations that will be done repetitively
      XJ=(0.0,1.0)
      radius2=radius**2
      twopi=2.0*PI
      S=1
      K=0
      raysum=0
      OPEN(UNIT=5,file=INPUTFILE)
C      print*, INPUTFILE

C     Count the number of elements in the input file

 11   Read(5,*,END=12)temp1
      K=K+1
      Go TO 11

C allocate the correct size to the signal and fft arrays

 12   allocate(inputsignal(K))
      allocate(outputsignal(K/2+1))
      allocate(inputarray(K/2,3))
      close(5)
      
C     Read in the input signal
      OPEN(UNIT=5,file=INPUTFILE)
      DO 30 W=1,K
         Read(5,*,END=5)inputsignal(W)
         S=S+1
 30   CONTINUE

C  Take the fft of the input signal with fftw
 5    sizefft=K
      sizeffttwo=sizefft/2.0
      call dfftw_plan_dft_r2c_1d(plan,sizefft, inputsignal, outputsignal
     *     , FFTW_ESTIMATE)
      call dfftw_execute_dft_r2c(plan, inputsignal, outputsignal)
      call dfftw_destroy_plan(plan)
      allocate(timearray(sizefft))
      allocate(ampinitial(sizeffttwo))
      allocate(phaseinitial(sizeffttwo))

C     Create initial signal
      allocate(airabsorb(sizeffttwo))
C      print*, outputsignal
      DO 13 K=1, sizeffttwo
         inputarray(K,1)=(K)*Fs/2*1/(sizeffttwo)
         inputarray(K,2)=abs(outputsignal(K+1)/sizefft)
         inputarray(K,3)=ATAN2(imagpart(outputsignal(K+1)/sizefft),
     *        realpart(outputsignal(K+1)/sizefft))
         airabsorb(K)=ABSORPTION(ps,inputarray(K,1),hr,Temp)
 13   CONTINUE
C      print*, airabsorb
      DO 14 K=1, sizefft
         timearray(K)=(K-1)*1/Fs
 14   CONTINUE
      deallocate(outputsignal)
      deallocate(inputsignal)
C     Set initial values

      Vinitial=(/xinitial,yinitial,zinitial/)
      xiinitial=COS(phi)*sin(theta)
      ninitial=SIN(phi)*sin(theta)
      zetainitial=cos(theta)
      length=sqrt(xiinitial*xiinitial+ninitial*ninitial+zetainitial*
     *     zetainitial)
      Finitial=(/xiinitial,ninitial,zetainitial/)
C      print*, Finitial
      tmp=(Finitial(1)*Vinitial(1)+Finitial(2)*Vinitial(2)+Finitial(3)*
     *     Vinitial(3))
      PLANEABC=(/Finitial(1),Finitial(2),Finitial(3),tmp/)
C     Create initial boom array

      yspace=boomspacing*abs(cos(phi))
      zspace=boomspacing*abs(sin(theta))
      if (xmin.eq.xmax) then
         RAYMAX=int((ymax-ymin)/yspace)*int((zmax-zmin)/zspace)
      elseif(ymin.eq.ymax)then
         RAYMAX=int((xmax-xmin)/xspace)*int((zmax-zmin)/zspace)
      elseif(zmin.eq.zmax) then
         RAYMAX=int((ymax-ymin)/yspace)*int((xmax-xmin)/xspace)
      endif
      allocate(boomarray(RAYMAX,3))
      PRINT*, RAYMAX
      call InitialGrid(boomspacing,PLANEABC(1),PLANEABC(2),PLANEABC(3),
     *     PLANEABC(4),theta,phi,xmin,ymin,zmin,xmax,ymax,zmax,boomarray
     *     ,RAYMAX,sizex,sizey,sizez)
C      DO 500 j=1, RAYMAX
C         print*, boomarray(j,1),boomarray(j,2),boomarray(j,3)
C 500     CONTINUE
C      print*, 'created boom array'

C     Create a receiver array, include a receiver file. 

      allocate(alphanothing(sizeffttwo))
      alphanothing=0.0
C      INCLUDE 'ReceiverEarlevel.f'
C      INCLUDE 'ReceiverCrossSection.f'
C      INCLUDE 'ReceiverFrontWall.f'
C      INCLUDE 'ReceiverRoof.f'
C      INCLUDE 'ReceiverPointSourceVal.f'
      INCLUDE 'ReceiverPointSource.f'
C      INCLUDE 'ReceiverGroundH3W35.f'
C      INCLUDE 'ReceiverGroundH3W35Deg90.f'
C      INCLUDE 'ReceiverGroundH3W7.f'
C      INCLUDE 'ReceiverGroundH3W7Deg90.f'
C      INCLUDE 'ReceiverGroundH3W14.f'
C      INCLUDE 'ReceiverGroundH3W14Deg90.f'
C      INCLUDE 'ReceiverGroundH6W35.f'
C      INCLUDE 'ReceiverGroundH6W35Deg90.f'
C      INCLUDE 'ReceiverGroundH6W7.f'
C      INCLUDE 'ReceiverGroundH6W7Deg90.f'
C      INCLUDE 'ReceiverGroundH6W14.f'
C      INCLUDE 'ReceiverGroundH6W14Deg90.f'
C      INCLUDE 'ReceiverGroundH12W35.f'
C      INCLUDE 'ReceiverGroundH12W35Deg90.f'
C      INCLUDE 'ReceiverGroundH12W7.f'
C      INCLUDE 'ReceiverGroundH12W7Deg90.f'
C      INCLUDE 'ReceiverGroundH12W14.f'
C      INCLUDE 'ReceiverGroundH12W14Deg90.f'
C      INCLUDE 'ReceiverGroundH24W35.f'
C      INCLUDE 'ReceiverGroundH24W35Deg90.f'
C      INCLUDE 'ReceiverGroundH24W7.f'
C      INCLUDE 'ReceiverGroundH24W7Deg90.f'
C      INCLUDE 'ReceiverGroundH24W14.f'
C      INCLUDE 'ReceiverGroundH24W14Deg90.f'
C      INCLUDE 'ReceiverPointSourceMuse.f'
C      INCLUDE 'ReceiverGroundG2H3W55.f'
C      INCLUDE 'ReceiverGroundG2H3W7.f'
C      INCLUDE 'ReceiverGroundG2H3W14.f'
C      INCLUDE 'ReceiverGroundG2H6W55.f'
C      INCLUDE 'ReceiverGroundG2H6W7.f'
C      INCLUDE 'ReceiverGroundG2H6W14.f'
C      INCLUDE 'ReceiverGroundG2H12W55.f'
C      INCLUDE 'ReceiverGroundG2H12W7.f'
C      INCLUDE 'ReceiverGroundG2H12W14.f'
C      INCLUDE 'ReceiverGroundG2H24W55.f'
C      INCLUDE 'ReceiverGroundG2H24W7.f'
C      INCLUDE 'ReceiverGroundG2H24W14.f'
      sum=0

C     deallocate temparary receiver arrays
      deallocate(receiverarray1)
      if (planenum.ge.2) deallocate(receiverarray2)
      if (planenum.ge.3) deallocate(receiverarray3)
      if (planenum.ge.4) deallocate(receiverarray4)
      if (planenum.ge.5) deallocate(receiverarray5)
      if (planenum.ge.6) deallocate(receiverarray6)
      if (planenum.ge.7) deallocate(receiverarray7)

C     initialize normalization factor 
      normalization=(PI*radius2)/(boomspacing*boomspacing)

      allocate(temparray(arraysize,sizeffttwo,6))
      allocate(timetemparray(arraysize,sizefft,5))
      DO 18 D=1, arraysize
         DO 15 W=1, sizeffttwo
            temparray(D,W,1)=receiverarray(D,1)
            temparray(D,W,2)=receiverarray(D,2)
            temparray(D,W,3)=receiverarray(D,3)
            temparray(D,W,4)=inputarray(W,1)
            temparray(D,W,5)=0.0 
            temparray(D,W,6)=0.0 
 15      CONTINUE
 18   CONTINUE
C            print*,temparray(1,sizeffttwo,4)
C            print*,temparray(1,sizeffttwo+1,4)


C     Define ground plane

      groundheight=0.000000000
      GROUNDABC=(/0.000000000,0.000000000,1.00000000/)
      GROUNDD=-groundheight
      nground=(/0.0,0.0,1.0/)
      allocate(alphaground(sizeffttwo))
C     Allocate absorption coefficients for each surface for each frequency

      DO 17 D=1, sizeffttwo
         if(inputarray(D,1).ge.0.0.or.inputarray(D,1).lt.88.0)then
            alphaground(D)=tempalphaground(1)
         elseif(inputarray(D,1).ge.88.0.or.inputarray(D,1).lt.177.0)then
            alphaground(D)=tempalphaground(2)
         elseif(inputarray(D,1).ge.177.0.or.inputarray(D,1).lt.355.0)
     *           then
            alphaground(D)=tempalphaground(3)
         elseif(inputarray(D,1).ge.355.0.or.inputarray(D,1).lt.710.0)
     *           then
            alphaground(D)=tempalphaground(4)
         elseif(inputarray(D,1).ge.710.0.or.inputarray(D,1).lt.1420.0)
     *           then
            alphaground(D)=tempalphaground(5)
         elseif(inputarray(D,1).ge.1420.0.or.inputarray(D,1).lt.2840.0)
     *           then
            alphaground(D)=tempalphaground(6)
         elseif(inputarray(D,1).ge.2840.0.or.inputarray(D,1).lt.5680.0)
     *           then
            alphaground(D)=tempalphaground(7)
         elseif(inputarray(D,1).ge.5680.0.or.inputarray(D,1).lt.
     *           inputarray(sizeffttwo,1))then
            alphaground(D)=tempalphaground(8)
         endif
 17   CONTINUE
      allocate(alphabuilding(absorbplanes,sizeffttwo))
      DO 9 W=1,absorbplanes
         DO 8 D=1, sizeffttwo
            if(inputarray(D,1).ge.0.0.or.inputarray(D,1).lt.88.0)then
               alphabuilding(W,D)=tempalphabuilding(W,1)
            elseif(inputarray(D,1).ge.88.0.or.inputarray(D,1).lt.177.0)
     *              then
               alphabuilding(W,D)=tempalphabuilding(W,2)
            elseif(inputarray(D,1).ge.177.0.or.inputarray(D,1).lt.355.0)
     *              then
               alphabuilding(W,D)=tempalphabuilding(W,3)
            elseif(inputarray(D,1).ge.355.0.or.inputarray(D,1).lt.710.0)
     *              then
               alphabuilding(W,D)=tempalphabuilding(W,4)
            elseif(inputarray(D,1).ge.710.0.or.inputarray(D,1).lt.
     *              1420.0)then
               alphabuilding(W,D)=tempalphabuilding(W,5)
            elseif(inputarray(D,1).ge.1420.0.or.inputarray(D,1).lt.
     *              2840.0)then
               alphabuilding(W,D)=tempalphabuilding(W,6)
            elseif(inputarray(D,1).ge.2840.0.or.inputarray(D,1).lt.
     *              5680.0)then
               alphabuilding(W,D)=tempalphabuilding(W,7)
            elseif(inputarray(D,1).ge.5680.0.or.inputarray(D,1).lt.
     *              inputarray(sizeffttwo,1))then
               alphabuilding(W,D)=tempalphabuilding(W,8)
            endif
 8       CONTINUE
 9    CONTINUE

C     Mesh the patches for the environment.  Include patching file. 
      if(radiosity.eq.1)then
C         INCLUDE 'SingleBuildingGeometry.f'
C         INCLUDE 'SingleBuildingGeometryVal.f'
C         INCLUDE 'NASAEMBuilding.f'
C         INCLUDE 'NASAMuseBuilding.f'
C         INCLUDE 'NASAMuseBuildComplex.f'
C         INCLUDE 'UrbanCanyonH3L35Geo.f'
C         INCLUDE 'UrbanCanyonH3L35GeoDeg90.f'
C         INCLUDE 'UrbanCanyonH3L7Geo.f' 
C         INCLUDE 'UrbanCanyonH3L7GeoDeg90.f'
C         INCLUDE 'UrbanCanyonH3L14Geo.f'
C         INCLUDE 'UrbanCanyonH3L14GeoDeg90.f'
C         INCLUDE 'UrbanCanyonH6L35Geo.f'
C         INCLUDE 'UrbanCanyonH6L35GeoDeg90.f'
C         INCLUDE 'UrbanCanyonH6L7Geo.f'
C         INCLUDE 'UrbanCanyonH6L7GeoDeg90.f'
C         INCLUDE 'UrbanCanyonH6L14Geo.f'
C         INCLUDE 'UrbanCanyonH6L14GeoDeg90.f'
C         INCLUDE 'UrbanCanyonH12L35Geo.f'
C         INCLUDE 'UrbanCanyonH12L35GeoDeg90.f'
C         INCLUDE 'UrbanCanyonH12L7Geo.f'
C         INCLUDE 'UrbanCanyonH12L7GeoDeg90.f'
C         INCLUDE 'UrbanCanyonH12L14Geo.f'
C         INCLUDE 'UrbanCanyonH12L14GeoDeg90.f'
C         INCLUDE 'UrbanCanyonH24L35Geo.f'
C         INCLUDE 'UrbanCanyonH24L35GeoDeg90.f'
C         INCLUDE 'UrbanCanyonH24L7Geo.f'
C         INCLUDE 'UrbanCanyonH24L7GeoDeg90.f'
C         INCLUDE 'UrbanCanyonH24L14Geo.f'
c         INCLUDE 'UrbanCanyonH24L14GeoDeg90.f'
C         INCLUDE 'UrbanCanyonG2H3L55Geo.f'
C         INCLUDE 'UrbanCanyonG2H3L7Geo.f' 
C         INCLUDE 'UrbanCanyonG2H3L14Geo.f'
C         INCLUDE 'UrbanCanyonG2H6L55Geo.f'
C         INCLUDE 'UrbanCanyonG2H6L7Geo.f'
C         INCLUDE 'UrbanCanyonG2H6L14Geo.f'
C         INCLUDE 'UrbanCanyonG2H12L55Geo.f'
C         INCLUDE 'UrbanCanyonG2H12L7Geo.f'
C         INCLUDE 'UrbanCanyonG2H12L14Geo.f'
C         INCLUDE 'UrbanCanyonG2H24L55Geo.f'
C         INCLUDE 'UrbanCanyonG2H24L7Geo.f'
C         INCLUDE 'UrbanCanyonG2H24L14Geo.f'
         diffusion=percentdiffuse
         diffusionground=0.0
      else
         diffusion=0.0
         diffusionground=0.0
      endif
      count=0
C      print*, 'normalization',normalization
C     Loop through the intial ray locations
      DO 40 ray=1,RAYMAX,1
C         ray = 607
         hitcount=0
         tmpsum=0.0
         doublehit=0
         DO 24 W=1, sizeffttwo
            ampinitial(W)=inputarray(W,2)/normalization
            phaseinitial(W)=inputarray(W,3)
 24      CONTINUE
C            print*,ampinitial(sizeffttwo-4:sizeffttwo)
         Vinitial=(/BOOMARRAY(ray,1),BOOMARRAY(ray,2),
     *        BOOMARRAY(ray,3)/)
         if (h.lt.2*radius)then 
C            print*, 'h is less than 2r'
            Call abort
         endif
         F=Finitial

         veci=Vinitial
C     Making small steps along the ray path.  For each step we should return, 
C     location, phase and amplitude
         DO 10 I=1,IMAX,1
            dxreceiver=HUGE
C     Find the closest sphere and store that as the distance
C            print*, veci
C            print*, 'F: ',F 
C            print*, 'dx: ', dx
            DO 16 Q=1,arraysize,1 
               CALL SPHERECHECK(receiverarray(Q,1:3),
     *              radius2,F,veci,tempreceiver)
               if(receiverhit.ge.1) then
                  if(lastreceiver(1).eq.receiverarray(Q,1).and.
     *                 lastreceiver(2).eq.receiverarray(Q,2).and.
     *                 lastreceiver(3).eq.receiverarray(Q,3))then
                     tempreceiver=HUGE
                  endif
                  if(F(1).eq.checkdirection(1).and.F(2).eq.
     *                 checkdirection(2).and.F(3).eq.
     *                 checkdirection(3))then
                     OC(1)=receiverarray(Q,1)-veci(1)
                     OC(2)=receiverarray(Q,2)-veci(2)
                     OC(3)=receiverarray(Q,3)-veci(3)
                     OCLength=OC(1)*OC(1)+OC(2)*OC(2)+OC(3)*OC(3)
                     if(OCLength.lt.radius2)then
                        tempreceiver=HUGE
                     endif
                  endif
               endif
               if(receiverhit.ge.2)then
                  if(lastreceiver2(1).eq.receiverarray(Q,1).and.
     *                 lastreceiver2(2).eq.receiverarray(Q,2).and.
     *                 lastreceiver2(3).eq.receiverarray(Q,3))then
                     tempreceiver=HUGE
                  endif
               endif
               if (tempreceiver.lt.dxreceiver)then      
                  dxreceiver=tempreceiver
                  receiverpoint(1)=receiverarray(Q,1)
                  receiverpoint(2)=receiverarray(Q,2)
                  receiverpoint(3)=receiverarray(Q,3)
               elseif (tempreceiver.eq.dxreceiver.and.
     *                 tempreceiver.ne.HUGE)then
                  receivercheck=tempreceiver
C                  print*, 'receivercheck',receivercheck
                  if(receiverarray(Q,1).eq.receiverpoint(1).and.
     *                 receiverarray(Q,2).eq.receiverpoint(2).and.
     *                 receiverarray(Q,3).eq.receiverpoint(3))then
                     doublehit=0
                  else
                     receiverpoint2(1)=receiverarray(Q,1)
                     receiverpoint2(2)=receiverarray(Q,2)
                     receiverpoint2(3)=receiverarray(Q,3)
                     doublehit=1
                  endif
               endif
 16         CONTINUE
C     Check Intersection with ground plane
            GROUNDN=GroundABC
            Ftemp=F
            vecitemp=veci
            GROUNDVD=GROUNDn(1)*Ftemp(1)+GROUNDN(2)*Ftemp(2)+GROUNDN(3)*
     *           Ftemp(3)
            if (groundhit.eq.1) then
               dxground=huge
            elseif (GROUNDVD.ne.0.0) then
               GROUNDVO=((GROUNDn(1)*vecitemp(1)+GROUNDn(2)*vecitemp(2)+
     *              GROUNDn(3)*vecitemp(3))+GROUNDD)
               dxground1=(-1.0D0)*GROUNDVO*(1.0D0)/GROUNDVD
               dxground=dxground1
               Vecip1=veci+dxground*F
               Vecip1temp=vecitemp+dxground1*Ftemp
               tmp=(GROUNDabc(1)*Vecip1(1)+GROUNDabc(2)*Vecip1(2)+
     *              GROUNDabc(3)*Vecip1(3)+GROUNDD)
               tmp=(GROUNDn(1)*Vecip1temp(1)+GROUNDn(2)*Vecip1temp(2)+
     *              GROUNDn(3)*Vecip1temp(3)+GROUNDD)
               if (dxground.lt.0.0) dxground=HUGE
            else
               dxground=huge
            endif
C     Check intersection with building
            dxbuilding=HUGE
C            if(buildinghit.eq.1) then
C               dxbuilding=huge
C            else
               hit=0
               planehit=0
C     Check intersection with Boxes
               DO 28 Q=1,boxnumber,1 
                  call BOX(Boxarraynear(Q,1:3), Boxarrayfar(Q,1:3)
     *                 ,Veci,F,dxnear, dxfar, hit, planehit)
                  if (dxnear.lt.dxbuilding)then
                     dxbuilding=dxnear
C                     print*, 'dxbuilding',dxnear
                     Vecip1=veci+dxbuilding*F
                     whichbox=Q
                     Call PLANE(Vecip1, boxarraynear(whichbox,1:3), 
     *                    boxarrayfar(whichbox,1:3), planehit, nbox)
                  endif
C                  print*, 'dxbuilding',dxbuilding
 28            CONTINUE
C     Check intersection with Triangles
               if(TriangleNumber.gt.0)then
                  DO 32 Q=1,TriangleNumber,1 
                     call Polygon(veci,F,Q,3,TriangleNumber,PointNumbers
     *                    ,Trianglearray,BuildingPoints,normal,
     *                    FaceNormalNo,FaceNormals,dxnear,behind)
                     if (dxnear.lt.dxbuilding)then
                        dxbuilding=dxnear
                        nbox=normal
                        whichbox=Q
                     endif
 32               CONTINUE
               endif
C     Check intersection with Squares
               if(SquareNumber.gt.0)then
                  DO 33 Q=1,SquareNumber,1 
                     call Polygon(veci,F,Q,4,SquareNumber,PointNumbers,
     *                    SquareArray,BuildingPoints,normal,FaceNormalNo
     *                    ,FaceNormals,dxnear,behind)
                     if (dxnear.lt.dxbuilding)then
                        dxbuilding=dxnear
                        nbox=normal
                        whichbox=Q
                     endif
 33               CONTINUE
               endif
C            endif
            buildinghit=0
            receiverhit=0
            groundhit=0
C     Check to see if ray hits within step size
C            print*, dxreceiver, dxground, dxbuilding
            if (dxreceiver.lt.h.or.dxground.lt.h.or.dxbuilding.lt.h)
     *           then
               dx=MIN(dxreceiver,dxground,dxbuilding)
               tmpsum=tmpsum+dx
C     if the ray hits a receiver, store in an array.  If the ray hits twice
C     Create two arrays to store in.
               if (dx.eq.dxreceiver) then
                  sum=sum+1
                  Vecip1=veci+dx*F
                  veci=Vecip1
C                  print*, 'hit receiver at step ',I
                  receiverhit=1
                  checkdirection=F
                  
                  if(doublehit.eq.1)then
                     receiverhit=2
                  endif
                  hitcount=hitcount+1
C                  print*, 'hit receiver',sum,tmpsum,receiverpoint
                  DO 20 W=1, sizeffttwo
                     m=airabsorb(W)
                     lambda=soundspeed/inputarray(W,1)
                     phasefinal=phaseinitial(W)-(twopi*dx)/lambda   
                     ampfinal=ampinitial(W)*(1-alphanothing(W))*
     *                    exp(-m*dx)
                     ampinitial(W)=ampfinal
                     phaseinitial(W)=mod(phasefinal,twopi)
                     if (phaseinitial(W).GT.PI) then
                        phaseinitial(W)=phaseinitial(W)-twopi
                     endif
                     if(doublehit.eq.1)then
                        if(W.eq.1)allocate(outputarray1(
     *                       sizeffttwo,6))
                        if(W.eq.1)allocate(dhoutputarray1(
     *                       sizeffttwo,6))   
                        outputarray1(W,1)=inputarray(W,1)
                        outputarray1(W,2)=receiverpoint(1)
                        outputarray1(W,3)=receiverpoint(2)
                        outputarray1(W,4)=receiverpoint(3)
                        outputarray1(W,5)=ampinitial(W)/2.0
                        outputarray1(W,6)=phaseinitial(W)
                        dhoutputarray1(W,1)=inputarray(W,1)
                        dhoutputarray1(W,2)=receiverpoint2(1)
                        dhoutputarray1(W,3)=receiverpoint2(2)
                        dhoutputarray1(W,4)=receiverpoint2(3)
                        dhoutputarray1(W,5)=ampinitial(W)/2.0
                        dhoutputarray1(W,6)=phaseinitial(W)
                        lastreceiver(1)=receiverpoint(1)
                        lastreceiver(2)=receiverpoint(2)
                        lastreceiver(3)=receiverpoint(3)
                        lastreceiver2(1)=receiverpoint2(1)
                        lastreceiver2(2)=receiverpoint2(2)
                        lastreceiver2(3)=receiverpoint2(3)
                     else
C                        print*, 'this happens outputarray'
                        if(W.eq.1)allocate(outputarray1(sizeffttwo,6))
                        outputarray1(W,1)=inputarray(W,1)
                        outputarray1(W,2)=receiverpoint(1)
                        outputarray1(W,3)=receiverpoint(2)
                        outputarray1(W,4)=receiverpoint(3)
                        outputarray1(W,5)=ampinitial(W)
                        outputarray1(W,6)=phaseinitial(W)
                        lastreceiver(1)=receiverpoint(1)
                        lastreceiver(2)=receiverpoint(2)
                        lastreceiver(3)=receiverpoint(3)
                     endif                 
 20               CONTINUE
C                   print*,ampinitial(sizeffttwo-4:sizeffttwo)

C                  print*, 'assigned to outputarray1'
C                  print*, outputarray1(1,2), outputarray1(1,3), 
C                  print*,'magnitude', temparray(2,:,5)

C     *                 outputarray1(1,4)
                  Call receiverHITFUNC(sizefft,outputarray1,
     *                 arraysize,temparray)
C                  print*, 'receiverHITFUNC completed'
                  if (doublehit.eq.1)then
                     call receiverHITFUNC(sizefft,dhoutputarray1,
     *                    arraysize,temparray)
                     count=count+1
                  endif
                  count=count+1
                  deallocate(outputarray1)
                  if(doublehit.eq.1)then
                     deallocate(dhoutputarray1)
                  endif
C                  print*, 'got to the end of receiver hit'
               endif
C               print*,'direction', temparray(2,:,5)
C     If the ray hits the ground then bounce off the ground and continue
               if (abs(dx-dxground).lt.10.0**(-13.0)) then
                  Vecip1=veci+dxground*F
                  tmp=(GROUNDabc(1)*Vecip1(1)+GROUNDabc(2)*Vecip1(2)+
     *                 GROUNDabc(3)*Vecip1(3)+GROUNDD)
                  if(tmp.ne.GROUNDD) Vecip1(3)=0.0
C                  print*,'hit ground at step ', I
                  veci=Vecip1
                  dot1=(F(1)*nground(1)+F(2)*nground(2)+F(3)*nground(3))
                  n2=(nground(1)*nground(1)+nground(2)*nground(2)+
     *                 nground(3)*nground(3))
                  r=F-2.0*(dot1/n2)*nground
                  length=sqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
                  F=(/r(1),r(2),r(3)/)
                  groundhit=1
                  twopidx=twopi*dxground
C     Loop through all the frequencies
                  DO 21 W=1, sizeffttwo
                        m=airabsorb(W)
                        lambda=soundspeed/inputarray(W,1)
                        phasefinal=phaseinitial(W)-(twopidx)/lambda
                        ampfinal=ampinitial(W)*(1.0-alphaground(W))
     *                       *(1.0-diffusionground)*exp(-m*dxground)
                        phaseinitial(W)=mod(phasefinal,twopi)
                        if (phaseinitial(W).GT.PI) then
                           phaseinitial(W)=phaseinitial(W)-twopi
                        endif
                        if(radiosity.eq.1.and.(diffusionground.ne.0.0))
     *                       then
                           DO 25 Q=1,PatchNo
                              if (formfactors(1,Q,2).eq.1)then
                                 if(veci(1).le.(patcharray(Q,W,1)+0.5*
     *                                patcharray(Q,W,4)).and.veci(1).ge.
     *                                (patcharray(Q,W,1)-0.5*patcharray
     *                                (Q,W,4)))then
                                    if(veci(2).le.(patcharray(Q,W,2)+0.5
     *                                   *patcharray(Q,W,5)).and.veci(2)
     *                                   .ge.(patcharray(Q,W,2)-0.5*
     *                                   patcharray(Q,W,5)))then
                                       if(veci(3).le.(patcharray(Q,W,3)+
     *                                      0.5*patcharray(Q,W,6)).and.
     *                                      veci(3).ge.(patcharray(Q,W,3
     *                                      )-0.5*patcharray(Q,W,6)))
     *                                      then
                                          temp2=cmplx(abs(patcharray
     *                                         (Q,W,7))*exp(XJ*
     *                                         patcharray(Q,W,8)))
                                          temp3=cmplx(abs(ampinitial(W)*
     *                                         (1.0-alphaground(W))*
     *                                         diffusionground*exp(-m*
     *                                         dxground))*exp(XJ*
     *                                         phasefinal))
                                          temp4=temp2+temp3
                                          patcharray(Q,W,7)=abs(temp4)
                                          patcharray(Q,W,8)=ATAN2(
     *                                         imagpart(temp4),realpart
     *                                         (temp4))
                                          GOTO 26
                                       endif
                                    endif
                                 endif
 26                              CONTINUE
                              endif
 25                        CONTINUE
                        endif   
                     ampinitial(W)=ampfinal                           
 21               CONTINUE
C                print*,ampinitial(sizeffttwo-4:sizeffttwo)

               endif
C               print*,phaseinitial
C               print*, 'dxground: ', dxground
C               print*, 'alphaground: ', alphaground
C               print*, 'diffusionground: ', diffusionground

               
C     if the ray hits the building then change the direction and continue
C               print*, 'dx: ',dx 
C               print*, 'dxbuilding: ', dxbuilding
               if (dx.eq.dxbuilding) then
                  Vecip1=veci+dx*F
                  veci=Vecip1
C                  print*, 'hit building at step ', I
                  n2=(nbox(1)*nbox(1)+nbox(2)*nbox(2)+nbox(3)*nbox(3))
                  nbuilding=nbox/sqrt(n2)
                  dot1=(F(1)*nbuilding(1)+F(2)*nbuilding(2)+F(3)*
     *                 nbuilding(3))
                  r=F-2.0*(dot1/n2)*nbuilding
                  length=sqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
                  F=(/r(1),r(2),r(3)/)
                  buildinghit=1
                  twopidx=twopi*dx
                  DO 22 W=1, sizeffttwo
                        if(complexabsorption.eq.1)then
                           if (absorbplanes.eq.2)then
                           if(veci(3).gt.0.0.and.veci(3).le.height1)then
                              alpha=alphabuilding(1,W)
                           elseif(veci(3).gt.height1.and.veci(3).le.
     *                             height2)then
                              alpha=alphabuilding(2,W)
                           endif
                           endif
                           if(absorbplanes.eq.3)then
                           if(veci(3).gt.height2.and.veci(3).le.
     *                             height3)then
                              alpha=alphabuilding(3,W)
                           endif
                           endif
                           if(absorbplanes.eq.4)then
                            if(veci(3).gt.height3)then
                               alpha=alphabuilding(4,W)
                            endif
                           endif
                        else
                           alpha=alphabuilding(1,W)
                        endif
                        m=airabsorb(W)
                        lambda=soundspeed/inputarray(W,1)                
                        phasefinal=phaseinitial(W)-(twopidx)/lambda
                        ampfinal=ampinitial(W)*(1.0-alpha)*
     *                       (1.0-diffusion)*exp(-m*dx)
                        phaseinitial(W)=mod(phasefinal,twopi)
C                        if (W.le.5) then
C                        print*,'amps: ','(1.0-', alpha,
C     *                   ')*(1.0-', diffusion,')*exp(-', m,'*',
C     *                    dx
C                        endif
                        if (phaseinitial(W).GT.PI) then
                           phaseinitial(W)=phaseinitial(W)-twopi
                        endif
                        if(radiosity.eq.1)then
C     Loop through all patches if radiosity is turned on.
                           DO 29 Q=1,PatchNo
                              if (formfactors(1,Q,2).eq.2.0)then
                                 if((veci(1).le.(patcharray(Q,W,1)+0.5*
     *                                patcharray(Q,W,4))).and.(veci(1)
     *                                .ge.(patcharray(Q,W,1)-0.5*
     *                                patcharray(Q,W,4))))then
                                    if((veci(2).le.(patcharray(Q,W,2)+
     *                                   0.5*patcharray(Q,W,5))).and.
     *                                   (veci(2).ge.(patcharray(Q,W,2)-
     *                                   0.5*patcharray(Q,W,5))))then
                                       if((veci(3).le.(patcharray(Q,W,3)
     *                                      +0.5*patcharray(Q,W,6)))
     *                                      .and.(veci(3).ge.(patcharray
     *                                      (Q,W,3)-0.5*patcharray(Q,W,6
     *                                      ))))then
                                          temp2=cmplx(abs(patcharray
     *                                         (Q,W,7))*exp(XJ*
     *                                         patcharray(Q,W,8)))
                                          temp3=abs(ampinitial(W)*(1.0-
     *                                         alpha)*diffusion*exp
     *                                         (-m*dx))*exp(XJ*
     *                                         phaseinitial(W))
                                          temp4=temp2+temp3
                                          patcharray(Q,W,7)=abs(temp4)
                                          patcharray(Q,W,8)=ATAN2(
     *                                         imagpart(temp4),realpart(
     *                                         temp4))
                                          GOTO 27
                                       endif
                                    endif
                                 endif
 27                              CONTINUE
                              endif
 29                     CONTINUE
                     endif 
                  ampinitial(W)=ampfinal
 22            CONTINUE
C             print*,ampinitial(sizeffttwo-4:sizeffttwo)

            endif
C     If there was no interaction with buildings then proceed with one step. 
         else
            tmpsum=tmpsum+h
            Vecip1=veci+(h)*F
C            print*, 'no hit'
C            print*, phaseinitial(:5)
            veci=Vecip1
            twopih=twopi*h
            DO 23 W=1, sizeffttwo
C     Loop through all frequencies. 
               m=airabsorb(W)
               lambda=soundspeed/inputarray(W,1)
               phasefinal=phaseinitial(W)-(twopih)/lambda
               ampfinal=ampinitial(W)*(1-alphanothing(W))*
     *              exp(-m*h)
               ampinitial(W)=ampfinal  
               phaseinitial(W)=mod(phasefinal,twopi)
C               if (W.ge.(sizeffttwo-5)) then
C                  print*,'(1.0-',alphanothing(w),')*(1.0-',
CC     *             diffusion,
C     *             ')*exp(',-m,'*',h,')'
C               endif
               if (phaseinitial(W).GT.PI) then
                  phaseinitial(W)=phaseinitial(W)-twopi
               endif
 23         CONTINUE
C             print*,ampinitial(sizeffttwo-4:sizeffttwo)

         endif
 10   CONTINUE
      print*, 'finished ray', ray
 40   CONTINUE
C     Once all rays are complete.  Deallocate all arrays that are no longer needed
      deallocate(boomarray)
      deallocate(receiverarray)
      deallocate(ampinitial)
      deallocate(phaseinitial)
      allocate(Gk(PatchNo,sizeffttwo))
      allocate(Gkminus1(PatchNo,sizeffttwo))
      if (radiosity.eq.1)then
C     If radiosity is turned on then do the energy exchange. 
         KMAX=3
         Npatch=10
         DO 53 K=1,KMAX
            DO 55 D=1,PatchNo
               DO 54 W=1,sizeffttwo
                  if (K.eq.1)then
                     Gk(D,W)=cmplx(abs(patcharray(D,W,7)
     *                    )*exp(XJ*patcharray(D,W,8)))
                  else
                     Gkminus1(D,W)=Gk(D,W)
                     DO 56 I=1,PatchNo
                        if (I.eq.1)then
                           Gk(D,W)=0.0
                        else
                           if (formfactors(D,I,2).eq.1.0)then 
                              alpha=alphaground(W)
                           elseif (formfactors(D,I,2).eq.2.0)then
                              if( complexabsorption.eq.1)then
                                 if(absorbplanes.eq.2)then
                                 if(patcharray(I,1,3).gt.0.0.and.
     *                                patcharray(I,1,3).le.height1)then
                                    alpha=alphabuilding(1,W)
                                 elseif(patcharray(I,1,3).gt.height1
     *                                   .and.patcharray(I,1,3).le.
     *                                   height2)then
                                    alpha=alphabuilding(2,W)
                                 endif
                                 endif
                                 if(absorbplanes.eq.3)then
                                 if(patcharray(I,1,3).gt.height2
     *                                   .and.patcharray(I,1,3).le.
     *                                   height3)then
                                    alpha=alphabuilding(3,W)
                                 endif
                                 endif 
                                 if(absorbplanes.eq.4)then
                                 if(patcharray(I,1,3).gt.height3)
     *                                   then
                                    alpha=alphabuilding(4,W)
                                 endif
                                 endif
                              else
                                 alpha=alphabuilding(1,W)
                              endif
                           endif
                           m=airabsorb(W)
                           temp2=(1-alpha)*exp(-m*formfactors(D,I,3))*
     *                          formfactors(D,I,1)*Gkminus1(D,W)*exp(
     *                          -XJ*twopi*inputarray(W,1)*formfactors
     *                          (D,I,3)/soundspeed)
                           Gk(D,W)=Gk(D,W)+temp2
                        endif
 56                  CONTINUE
                  endif
 54            CONTINUE
C               print*, 'finished patch', D, 'of',PatchNo
 55         CONTINUE
C            print*, arraysize,PatchNo,sizefft
C     Do energy exchange with other receivers
            DO 50 D=1,arraysize
               DO 51 Q=1,PatchNo
                  Rlm=0.0
                  DO 58 I=1,Npatch
                     DO 59 J=1,Npatch
                        DO 60 S=1,Npatch
                           tmp1=((patcharray(Q,1,1)-.5*patcharray
     *                          (Q,1,4)+(patcharray(Q,1,4)/Npatch)
     *                          *(I-.5))-temparray(D,1,1))*((patcharray(
     *                          Q,1,1)-.5*patcharray(Q,1,4)+(patcharray(
     *                          Q,1,4)/Npatch)*(I-.5))-temparray(D,1,1))
                           tmp2=((patcharray(Q,1,2)-.5*patcharray
     *                          (Q,1,5)+(patcharray(Q,1,5)/Npatch)
     *                          *(J-.5))-temparray(D,1,2))*((patcharray(
     *                          Q,1,2)-.5*patcharray(Q,1,5)+(patcharray(
     *                          Q,1,5)/Npatch)*(J-.5))-temparray(D,1,2))
                           tmp3=((patcharray(Q,1,3)-.5*patcharray
     *                          (Q,1,6)+(patcharray(Q,1,6)/Npatch)
     *                          *(S-.5))-temparray(D,1,3))*((patcharray(
     *                          Q,1,3)-.5*patcharray(Q,1,6)+(patcharray(
     *                          Q,1,6)/Npatch)*(S-.5))-temparray(D,1,3))
                           Rlm=Rlm+1.0/(NPatch*Npatch*Npatch)*sqrt(tmp1+
     *                          tmp2+tmp3)
 60                     Continue
 59                  CONTINUE
 58               CONTINUE
                  Patchlength=sqrt(((patcharray(Q,1,1)-temparray(D,1,1))
     *                 *(patcharray(Q,1,1)-temparray(D,1,1)))+((patchar
     *                 ray(Q,1,2)-temparray(D,1,2))*(patcharray(Q,1,2)
     *                 -temparray(D,1,2)))+((patcharray(Q,1,3)-temparray
     *                (D,1,3))*(patcharray(Q,1,3)-temparray
     *                (D,1,3))))
                  Finitial(1)=(patcharray(Q,1,1)-temparray(D,1,1))
                  Finitial(2)=(patcharray(Q,1,2)-temparray(D,1,2))
                  Finitial(3)=(patcharray(Q,1,3)-temparray(D,1,3))
                  F=(/Finitial(1)/Rlm,Finitial(2)/Rlm,
     *                 Finitial(3)/Rlm/)
                  dxbuilding=HUGE
C     Check to see that the receiver is visible by patches. 
                  DO 57 I=1,boxnumber,1 
                     call BOX(boxarraynear(I,1:3),
     *                    Boxarrayfar(I,1:3),temparray(D,1,1:3),
     *                    F,dxnear,dxfar,hit, planehit)
                     if (dxnear.lt.dxbuilding)then
                        dxbuilding=dxnear
                     endif
                     if ((temparray(D,1,1).gt.boxarraynear(I,1)
     *                    .and.temparray(D,1,2).gt.boxarraynear(I,2)
     *                    .and.temparray(D,1,3).gt.boxarraynear(I,3))
     *                    .and.(temparray(D,1,1).lt.boxarrayfar(I,1)
     *                    .and.temparray(D,1,2).lt.boxarrayfar(I,2)
     *                    .and.temparray(D,1,3).lt.boxarrayfar(I,3)))
     *                    then
                        dxbuilding=-1.0
                     endif
 57               CONTINUE
                  if(TriangleNumber.gt.0)then
                     DO 67 I=1,TriangleNumber,1 
                        call Polygon(temparray(D,1,1:3),F,I,3,
     *                       TriangleNumber,PointNumbers,Trianglearray,
     *                       BuildingPoints,normal,FaceNormalNo,
     *                       FaceNormals,dxnear,behind)
                        if (dxnear.lt.dxbuilding)then
                           dxbuilding=dxnear
                        endif
 67                  CONTINUE
                  endif
                  if(SquareNumber.gt.0)then
                     DO 68 I=1,SquareNumber,1 
                        call Polygon(temparray(D,1,1:3),F,I,4,
     *                       SquareNumber,PointNumbers,SquareArray,
     *                       BuildingPoints,normal,FaceNormalNo,
     *                       FaceNormals,dxnear,behind)
                        if (dxnear.lt.dxbuilding)then
                           dxbuilding=dxnear
                        endif
 68                  CONTINUE
                  endif
                  if(PolyBuilding.gt.0)then
                     DO 69 I=1,PolyBuilding
C     Check that the receivers are not inside the building
                        CALL INSIDECHECK(temparray(D,1,1:3),
     *                       BuildingPoints,PointNumbers,F,Trianglearray
     *                       ,SquareArray,TriangleNumber,SquareNumber,
     *                       TriangleSequence,SquareSequence,Triangles,
     *                       Squares,PolyBuilding,I,inside,FaceNormalNo,
     *                       FaceNormals)
                        if(inside.eq.1)then
                           dxbuilding=-1.0
                        endif
 69                  CONTINUE
                  endif
                  vec3=FaceNormals(int(patcharray(Q,1,10)),1:3)
                  PIRlm2=PI*Rlm*Rlm
                  DO 52 W=1,sizeffttwo
                     if (Gk(Q,W).ne.0.0)then   
                        if(dxbuilding.lt.Patchlength)then
                           Ek=0.0
                        elseif(dxbuilding.ge.Patchlength)then
                           length=sqrt(vec3(1)*vec3(1)+
     *                          vec3(2)*vec3(2)+vec3(3)*vec3(3))
                           cosxilm=(vec3(1)*F(1)+vec3(2)*F(2)+vec3(3)*
     *                          F(3))/(Rlm*length)
                           m=airabsorb(W)
                           Ek=exp(-m*Rlm)*((cosxilm*Gk(Q,W)*exp(-XJ*
     *                          twopi*inputarray(W,1)*Rlm/soundspeed))/
     *                          (PIRlm2))
                        endif
                        temp2=cmplx(abs(temparray(D,W,5))*exp(XJ*
     *                       temparray(D,W,6)))
                        temp3=Ek+temp2
                        temparray(D,W,5)=ABS(temp3)
                        temparray(D,W,6)=ATAN2(imagpart(temp3)
     *                       ,realpart(temp3))
                     endif
 52               CONTINUE
 51            CONTINUE
C               print*, 'finished receiver', D, 'of', arraysize
 50         CONTINUE
 53      CONTINUE
      endif
C     Reconstruct the time signal
      CALL TIMERECONSTRUCT(sizefft, timearray, arraysize, temparray, 
     *     timetemparray)

C      print*,'magnitude: ',temparray(2,:,5)
C      print*,'direction: ',temparray(2,:,6)

C     Write out time signatures for each receiver. 
      OPEN(UNIT=20,file=OUTPUTFILE,status='new')
      true=Header(20)
      if(planenum.ge.1)then
         DO 61 W=1, sizefft
            true=TimeHeader(20,timetemparray(1,W,4),
     *           sizex1,sizey1,sizez1,planename1)
            DO 62 D=1, arraysize1
               write(20,*) timetemparray(D,W,1:3),timetemparray(D,W,5)
 62         CONTINUE
C            print*, 'finished time',timetemparray(1,W,4)
 61      CONTINUE
      endif
      if(planenum.ge.2)then
         DO 63 W=1, sizefft
            true=TimeHeader(20,timetemparray(1,W,4),
     *           sizex2,sizey2,sizez2,planename2)
            DO 64 D=arraysize1+1, arraysize1+arraysize2
               write(20,*) timetemparray(D,W,1:3),timetemparray(D,W,5)
 64         CONTINUE
 63      CONTINUE
      endif
      if(planenum.ge.3)then
         DO 65 W=1, sizefft
            true=TimeHeader(20,timetemparray(1,W,4),
     *           sizex3,sizey3,sizez3,planename3)
            DO 66 D=arraysize1+arraysize2+1,arraysize1+arraysize2+
     *           arraysize3
               write(20,*) timetemparray(D,W,1:3),timetemparray(D,W,5)
 66         CONTINUE
 65      CONTINUE
      endif
      if(planenum.ge.4)then
         DO 70 W=1, sizefft
            true=TimeHeader(20,timetemparray(1,W,4),
     *           sizex4,sizey4,sizez4,planename4)
            DO 71 D=arraysize1+arraysize2+arraysize3+1, arraysize1+
     *           arraysize2+arraysize3+arraysize4
               write(20,*) timetemparray(D,W,1:3),timetemparray(D,W,5)
 71         CONTINUE
 70      CONTINUE
      endif
      if(planenum.ge.5)then
         DO 72 W=1, sizefft
            true=TimeHeader(20,timetemparray(1,W,4),
     *           sizex5,sizey5,sizez5,planename5)
            DO 73 D=arraysize1+arraysize2+arraysize3+arraysize4+1, 
     *           arraysize1+arraysize2+arraysize3+arraysize4+arraysize5
               write(20,*) timetemparray(D,W,1:3),timetemparray(D,W,5)
 73         CONTINUE
 72      CONTINUE
      endif
      if(planenum.ge.6)then
         DO 74 W=1, sizefft
            true=TimeHeader(20,timetemparray(1,W,4),
     *           sizex6,sizey6,sizez6,planename6)
            DO 75 D=arraysize1+arraysize2+arraysize3+arraysize4+
     *           arraysize5+1,arraysize1+arraysize2+arraysize3+
     *           arraysize4+arraysize5+arraysize6
               write(20,*) timetemparray(D,W,1:3),timetemparray(D,W,5)
 75         CONTINUE
 74      CONTINUE
      endif
      if(planenum.ge.7)then
         DO 76 W=1, sizefft
            true=TimeHeader(20,timetemparray(1,W,4),
     *           sizex7,sizey7,sizez7,planename7)
            DO 77 D=arraysize1+arraysize2+arraysize3+arraysize4+
     *           arraysize5+arraysize6+1,arraysize1+arraysize2+
     *           arraysize3+arraysize4+arraysize5+arraysize6+arraysize7
               write(20,*) timetemparray(D,W,1:3),timetemparray(D,W,5)
 77         CONTINUE
 76      CONTINUE
      endif
      close(20)
      END
