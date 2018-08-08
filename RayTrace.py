# RayTrace

# Kimberly Lefkowitz created this program to propagate sonic booms around
# large structures, and to graduate. It is a ray tracing model that 
# will include specular and diffuse reflections. It will print out the
# sound field at ear height, at relevant microphone locations, and at 
# the building walls. It will read in the fft of a sonic boom signature.

# Dr. Riegel, William Costa, and Kory Seaton porting program from Fortran to python


# linking other python files with initialized variables and functions
import Parameterfile as PF
import BuildingGeometry as BG
import numpy as np
import Functions as F
# port and import receiver file

# Initialize counters 
PI=3.1415965358979323846
XJ=(0.0,1.0)
#radius2 = radius**2
twopi= 2.0*PI
S=1
K=0
raysum=0

# Count the Number of elements in the input file
with open(PF.INPUTFILE) as file:
    K=len(file.read())

# Allocate the correct size to the signal and fft arrays

inputsignal=np.zeros(K)
outputsignal=np.zeros(K/2+1)
# not working, fix it -Will
inputarray=np.zeros(((K/2),3))

# Read in the input signal
# i don't know what this does

with open(PF.INPUTFILE) as file:
    for W in range(1,K):
        file.read()

#take the fft of the input signal with fftw

sizefft = K
print(sizefft)
sizeffttwo=sizefft/2.0
np.fft.ifft(tempfft,sizefft)
timearray=np.zeros(sizefft)
ampinitial=np.zeros(sizeffttwo)
phaseinitial=np.zeros(sizeffttwo)

#       Create initial signal 

airabsorb = np.zeros(sizeffttwo)
for K in range(1,sizeffttwo):
    inputarray[K,1]=[K]*Fs/2*1/(sizeffttwo)
    inputarray[K,2]=abs(outputsignal(K+1)/sizefft)
    inputarray[K,3]=ATAN2(np.imag(outputsignal(K+1)/sizefft),realpart(outputsignal(K+1)/sizefft))
    airabsorb[K]=ABSORPTION(ps,inputarray[K,1],hr,Temp)
    print(airabsorb[K])
for K in range(1,sizefft):
    timearray[K]= F.ABSORPTION(ps,inputarray[K,1],hr,Temp)

#close output signal and input #somehow

#       Set initial values
Vinitial=[xinitial,yinitial,zinitial]
xiinitial=m.cos(phi)*m.sin(theta)
ninitial=m.sin(phi)*m.sin(theta)
zetainitial=m.cos(theta)
length=m.sqrt(xiinitial*xiinitial+ninitial*ninitial+zetainitial*zetainitial)
Finitial=[xiinitial,ninitial,zetainitial]
# print everything 
print(Finitial)
tmp=(Finitial(1)*Vinitial(1)+Finitial(2)*Vinitial(2)+Finitial(3)*Vinitial(3))
PLANEABC=[Finitial[1],Finitial[2],Finitial[3],tmp]
print(Vinitial)
print(xiinitial)
print(ninitial)
print(zetainitial)
print(length)
print(tmp)
print(PLANEABC)
#Our validation check 

#       Create initial boom array

yspace=boomspacing*abs(m.cos(phi))
print(yspace)
zspace=boomspacing*abs(m.sin(theta))
print(zspace)
if (xmin == xmax):
   RAYMAX=int((ymax-ymin)/yspace)*int((zmax-zmin)/zspace)
elif(ymin == ymax):
   RAYMAX=int((xmax-xmin)/xspace)*int((zmax-zmin)/zspace)
elif(zmin == zmax):
   RAYMAX=int((ymax-ymin)/yspace)*int((xmax-xmin)/xspace)
boomarray = np.zeros((RAYMAX,3))
#^ Looks like it will give same syntax bug as earlier, find np func that does similar to zeros
#but has space 
#czech it now -8/5/18
#it twerks
print(RAYMAX)
f.InitialGrid(boomspacing,PLANEABC[1],PLANEABC[2],PLANEABC[3],PLANEABC[4],theta,phi,xmin,ymin,zmin,xmax,ymax,zmax,boomarray,RAYMAX,sizex,sizey,sizez)

#     Create a receiver array, include a receiver file. 

alphanothing = np.zeros(sizeffttwo)
alphanothing = 0.0

import ReceiverPointSource as RPS 
sum = 0

#############################################################
#      deallocate(receiverarray1)
#      if (planenum.ge.2) deallocate(receiverarray2)
#      if (planenum.ge.3) deallocate(receiverarray3)
#      if (planenum.ge.4) deallocate(receiverarray4)
#      if (planenum.ge.5) deallocate(receiverarray5)
#      if (planenum.ge.6) deallocate(receiverarray6)
#      if (planenum.ge.7) deallocate(receiverarray7)
#############################################################
#Find function to deallocate 

#############################################################
#Check
np.empty(receiverarray1)
if planenum >=2 :
    np.empty(receiverarray2)
if planenum >=3 :
    np.empty(receiverarray3)
if planenum >=4 :
    np.empty(receiverarray4)
if planenum >=5 :
    np.empty(receiverarray5)
if planenum >=6 :
    np.empty(receiverarray6)
if planenum >=7 :
    np.empty(receiverarray7)
#############################################################

#       Initialize normalization factor 
normalization=(PI*radius2)/(boomspacing*boomspacing)

temparray=np.zeros(arraysize,sizeffttwo,6)
timetemparray=np.zeros(arraysize,sizeefft,5)
for D in range(1,arraysize):
    for W in range(1,sizeffttwo):
        temparray[D,W,1]=receiverarray[D,1]
        temparray[D,W,2]=receiverarray[D,2]
        temparray[D,W,3]=receiverarray[D,3]
        temparray[D,W,4]=inputarray[W,1]
        temparray[D,W,5]=0.0 
        temparray[D,W,6]=0.0 

#       Define ground plane

groundheight=0.000000000
GROUNDABC=[0.000000000,0.000000000,1.00000000]
GROUNDD=-groundheight
nground=[0.0,0.0,1.0]
alphaground=np.zeros(sizeffttwo)
#     Allocate absorption coefficients for each surface for each frequency

for D in range(1,sizeffttwo):
    if inputarray[D,1] >= 0.0 or inputarray[D,1] < 88.0 :
        alphaground[D]=tempalphaground[1]
    elif inputarray[D,1] >= 88.0 or inputarray[D,1] < 177.0 :
        alphaground[D]=tempalphaground[2]
    elif inputarray[D,1] >= 177.0 or inputarray[D,1] < 355.0 :
        alphaground[D]=tempalphaground[3]
    elif inputarray[D,1] >= 355.0 or inputarray[D,1] < 710.0 :
        alphaground[D]=tempalphaground[4]
    elif inputarray[D,1] >= 710.0 or inputarray[D,1] < 1420.0 :
        alphaground[D]=tempalphaground[5]
    elif inputarray[D,1] >= 1420.0 or inputarray[D,1] < 2840.0 :
        alphaground[D]=tempalphaground[6]
    elif inputarray[D,1] >= 2840.0 or inputarray[D,1] < 5680.0 :
        alphaground[D]=tempalphaground[7]
    elif inputarray[D,1] >= 5680.0 or inputarray[D,1] < inputarray[sizeffttwo,1]:
        alphaground[D]=tempalphaground[8]

alphabuilding = np.zeros((absorbplanes,sizeffttwo))

for W in range(1,absorbplanes):
    for D in range(1,absorbplanes):
        if  inputarray[D,1] >= 0.0 or inputarray[D,1] < 88.0:
            alphabuilding[W,D]=tempalphabuilding[W,1]
        elif inputarray[D,1] >= 88.0 or inputarray[D,1] < 177.0:
            alphabuilding[W,D] = tempalphabuilding [W,2]


        elif inputarray[D,1] >= 177.0 or inputarray[D,1] < 355.0 :
            alphabuilding[W,D] = tempalphabuilding[W,3]
        elif inputarray[D,1] >= 355.0 or inputarray[D,1] < 710.0 :
            alphabuilding[W,D] = tempalphabuilding[W,4]
        elif inputarray[D,1] >= 710.0 or inputarray[D,1] < 1420.0 :
            alphabuilding[W,D] = tempalphabuilding[W,5]
        elif inputarray[D,1] >= 1420.0 or inputarray[D,1] < 2840.0 :
            alphabuilding[W,D] = tempalphabuilding[W,6]
        elif inputarray[D,1] >= 2840.0 or inputarray[D,1] < 5680.0 :
            alphabuilding[W,D] = tempalphabuilding[W,7]
        elif inputarray[D,1] >= 5680.0 or inputarray[D,1] < inputarray[sizeffttwo,1] :
            alphabuilding[W,D] = tempalphabuilding[W,8]

#        Mesh the patches for the environment.  Include patching file. 
if radiosity == 1 :
    import SingleBuildingGeometry
#    import SingleBuildingGeometryVal
#    import NASAEMBuilding
#    import NASAMuseBuilding
#    import NASAMuseBuildComplex
#    import UrbanCanyonH3L35Geo
#    import UrbanCanyonH3L35GeoDeg90
#    import UrbanCanyonH3L7Geo.
#    import UrbanCanyonH3L7GeoDeg90
#    import UrbanCanyonH3L14Geo
#    import UrbanCanyonH3L14GeoDeg90
#    import UrbanCanyonH6L35Geo
#    import UrbanCanyonH6L35GeoDeg90
#    import UrbanCanyonH6L7Geo
#    import UrbanCanyonH6L7GeoDeg90
#    import UrbanCanyonH6L14Geo
#    import UrbanCanyonH6L14GeoDeg90
#    import UrbanCanyonH12L35Geo
#    import UrbanCanyonH12L35GeoDeg90
#    import UrbanCanyonH12L7Geo
#    import UrbanCanyonH12L7GeoDeg90
#    import UrbanCanyonH12L14Geo
#    import UrbanCanyonH12L14GeoDeg90
#    import UrbanCanyonH24L35Geo
#    import UrbanCanyonH24L35GeoDeg90
#    import UrbanCanyonH24L7Geo
#    import UrbanCanyonH24L7GeoDeg90
#    import UrbanCanyonH24L14Geo
#    import UrbanCanyonH24L14GeoDeg90
#    import UrbanCanyonG2H3L55Geo
#    import UrbanCanyonG2H3L7Geo
#    import UrbanCanyonG2H3L14Geo
#    import UrbanCanyonG2H6L55Geo
#    import UrbanCanyonG2H6L7Geo
#    import UrbanCanyonG2H6L14Geo
#    import UrbanCanyonG2H12L55Geo
#    import UrbanCanyonG2H12L7Geo
#    import UrbanCanyonG2H12L14Geo
#    import UrbanCanyonG2H24L55Geo
#    import UrbanCanyonG2H24L7Geo
#    import UrbanCanyonG2H24L14Geo
    diffusion=percentdiffuse
    diffusionground=0.0
else:
    diffusion=0.0
    diffusionground=0.0

count=0
#C      print('normalization',normalization)
#    Loop through the intial ray locations
##############################################################


#     Loop through the intial ray locations
ray=0
while ray <= RAYMAX:
      hitcount=0
      tmpsum=0.0
      doublehit=0
      W=0
      while W < sizeffttwo:
            ampinitial(W)=inputarray(W,1)/normalization
            phaseinitial(W)=inputarray(W,2)
      Vinitial=([BOOMARRAY(ray,0),BOOMARRAY(ray,1),BOOMARRAY(ray,2)])
      if (h < (2*radius)): 
            print('h is less than 2r')
            break
      F=Finitial
      veci=Vinitial
# Making small steps along the ray path.  For each step we should return, 
# location, phase and amplitude
      I=0
      while I < IMAX:
            dxreceiver=HUGE
# Find the closest sphere and store that as the distance
            print(veci)
            Q=0
            while Q < arraysize:
                  tempreceiver=SPHERECHECK(receiverarray[Q,0:2],radius2,F,veci)
                  if (receiverhit > 1):
                        if(lastreceiver[0]==receiverarray[Q,0] and lastreceiver[1]==receiverarray[Q,1] and lastreceiver[2] == receiverarray[Q,2]):
                              tempreceiver=HUGE
                        if(F(0) == checkdirection[0] and F[1] == checkdirection[1] and F[2] == checkdirection[2]):
                              OC(0)=receiverarray[Q,0]-veci[0]
                              OC(1)=receiverarray[Q,1]-veci[1]
                              OC(2)=receiverarray[Q,2]-veci[2] 
                              OCLength=OC[0]*OC[0]+OC[1]*OC[1]+OC[2]*OCp[2]
                              if(OCLength < radius2):
                                    tempreceiver=HUGE
                  if(receiverhit >= 2):
                        if(lastreceiver2[0]== receiverarray[Q,0] and lastreceiver2[1]==receiverarray[Q,1] and lastreceiver2[2]==receiverarray[Q,2]):
                              tempreceiver=HUGE
                  if (tempreceiver < dxreceiver):      
                        dxreceiver=tempreceiver
                        receiverpoint[0]=receiverarray[Q,0]
                        receiverpoint[1]=receiverarray[Q,1]
                        receiverpoint[2]=receiverarray[Q,2]
                  elif (tempreceiver== dxreceiver and tempreceiver != HUGE):
                        receivercheck=tempreceiver
                        print('receivercheck',receivercheck)
                        if(receiverarray[Q,0]==receiverpoint[0] and receiverarray[Q,1]==receiverpoint[1] and receiverarray[Q,2]==receiverpoint[2]):
                              doublehit=0
                        else:
                              receiverpoint2[0]=receiverarray[Q,0]
                              receiverpoint2[1]=receiverarray[Q,1]
                              receiverpoint2[2]=receiverarray[Q,2]
                              doublehit=1
            Q += 1 
#     Check Intersection with ground plane
            GROUNDN=GroundABC
            GROUNDVD=GROUNDn[0]*F[0]+GROUNDN[1]*F[1]+GROUNDN[2]*F[2]
            if (groundhit==1):
                  dxground=huge
            elif (GROUNDVD!=0.0):
                  GROUNDVO=((GROUNDn[0]*veci[0]+GROUNDn[1]*veci[1]+GROUNDn[2]*veci[2])+GROUNDD)
                  dxground1=(-1.000)*GROUNDVO*(1.000)/GROUNDVD
                  dxground=dxground1
                  Vecip1=veci+dxground*F
                  tmp=(GROUNDabc[0]*Vecip1[0]+GROUNDabc[1]*Vecip1[1]+GROUNDabc[2]*Vecip1[2]+GROUNDD)                  
                  if (dxground < 0.0):
                        dxground=HUGE
            else:
                  dxground=HUGE
#     Check intersection with building
            dxbuilding=HUGE
            hit=0
            planehit=0
#     Check intersection with Boxes
            Q=0
            while Q <  boxnumber:
#got to here

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

                  dxnear, dxfar, hit, planehit=BOX(Boxarraynear[Q,0:2], Boxarrayfar[Q,0:2],Veci,F)
                  if (dxnear.lt.dxbuilding)then
                     dxbuilding=dxnear
C                     print*, 'dxbuilding',dxnear
                     Vecip1=veci+dxbuilding*F
                     whichbox=Q
                     Call PLANE(Vecip1, boxarraynear(whichbox,1:3), 
     *                    boxarrayfar(whichbox,1:3), planehit, nbox)
                  endif
C                  print*, 'dxbuilding',dxbuilding
                  Q+=1
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
                  receiverhit=1
                  checkdirection=F
                  
                  if(doublehit.eq.1)then
                     receiverhit=2
                  endif
                  hitcount=hitcount+1
                  print*, 'hit receiver',sum,tmpsum,receiverpoint
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
                        print*, 'this happens outputarray'
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
                  print*, 'assigned to outputarray1'
                  print*, outputarray1(1,2), outputarray1(1,3), 
     *                 outputarray1(1,4)
                  Call receiverHITFUNC(sizefft,outputarray1,
     *                 arraysize,temparray)
                  print*, 'receiverHITFUNC completed'
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
                  print*, 'got to the end of receiver hit'
               endif
C     If the ray hits the ground then bounce off the ground and continue
               if (abs(dx-dxground).lt.10.0**(-13.0)) then
                  Vecip1=veci+dxground*F
                  tmp=(GROUNDabc(1)*Vecip1(1)+GROUNDabc(2)*Vecip1(2)+
     *                 GROUNDabc(3)*Vecip1(3)+GROUNDD)
                  if(tmp.ne.GROUNDD) Vecip1(3)=0.0
                  print*,'hit ground'
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
               endif
               
C     if the ray hits the building then change the direction and continue
               if (dx.eq.dxbuilding) then
                  Vecip1=veci+dx*F
                  veci=Vecip1
                  print*, 'hit building'
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
            endif
         else
C     If there was no interaction with buildings then proceed with one step. 
            tmpsum=tmpsum+h
            Vecip1=veci+(h)*F
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
               if (phaseinitial(W).GT.PI) then
                  phaseinitial(W)=phaseinitial(W)-twopi
               endif
 23         CONTINUE
         endif
      I += 1
 10   CONTINUE
      print*, 'finished ray', ray
      ray += 1
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
               print*, 'finished patch', D, 'of',PatchNo
 55         CONTINUE
            print*, arraysize,PatchNo,sizefft
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
               print*, 'finished receiver', D, 'of', arraysize
 50         CONTINUE
 53      CONTINUE
      endif
C     Reconstruct the time signal
      CALL TIMERECONSTRUCT(sizefft, timearray, arraysize, temparray, 
     *     timetemparray)
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
            print*, 'finished time',timetemparray(1,W,4)
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
