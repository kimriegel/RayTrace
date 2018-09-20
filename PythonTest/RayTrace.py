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
import Functions as fun
import cmath as cm
import math as m
#Formerly RayTrace_KAR 
#Cutting down on files -G
# port and import receiver file
receiverhit=0
groundhit=0
#np.seterr(divide='ignore', invalid ='ignore')
#only here because window gets clogged in errors. Remove when done

# Initialize counters 
PI=3.1415965358979323846
XJ=(0.0,1.0)
radius2 = PF.radius**2
twopi= 2.0*PI
S=1
K=0
raysum=0
#read in input file
with open(PF.INPUTFILE) as IPFile:
      inputsignal=np.loadtxt(IPFile)
K=len(inputsignal)
HUGE=1000000.0
F=np.zeros([1,3])
alphanothing=np.array([0.0])

# Allocate the correct size to the signal and fft arrays

outputsignal=np.zeros(int(K/2+1))
# not working, fix it -Will
inputarray=np.zeros((int(K/2),3))
#take the fft of the input signal with fftw

sizefft=K
#print('sizefft is ',sizefft)
sizeffttwo=int(sizefft/2)
outputsignal=np.fft.fft(inputsignal,sizefft)
timearray=np.zeros(sizefft)
#print('sizeffttwo is ',sizeffttwo)
ampinitial=np.zeros(sizeffttwo)
phaseinitial=np.zeros(sizeffttwo)

#       Create initial signal 

airabsorb = np.zeros(sizeffttwo)
K=0
while (K<sizeffttwo):
    inputarray[K,0]=(K+1)*PF.Fs/2*1/(sizeffttwo)
    inputarray[K,1]=abs(outputsignal[K]/sizefft)
    inputarray[K,2]=np.arctan2(np.imag(outputsignal[K]/sizefft),np.real(outputsignal[K]/sizefft))
    airabsorb[K]=fun.ABSORPTION(PF.ps,inputarray[K,0],PF.hr,PF.Temp)
    K+=1
K=0
while (K<sizefft):
    timearray[K]=(K)*1/PF.Fs
    K+=1

#close output signal and input #somehow

#       Set initial values
Vinitial=[PF.xinitial,PF.yinitial,PF.zinitial]
xiinitial=m.cos(PF.phi)*m.sin(PF.theta)
ninitial=m.sin(PF.phi)*m.sin(PF.theta)
zetainitial=m.cos(PF.theta)
length=m.sqrt(xiinitial*xiinitial+ninitial*ninitial+zetainitial*zetainitial)
Finitial=np.array([xiinitial,ninitial,zetainitial])
# print everything 
tmp=(Finitial[0]*Vinitial[0]+Finitial[1]*Vinitial[1]+Finitial[2]*Vinitial[2])
PLANEABC=[Finitial[0],Finitial[1],Finitial[2],tmp]
#print(Vinitial)
# print(xiinitial)
# print(ninitial)
# print(zetainitial)
# print(length)
# print(tmp)
# print(PLANEABC)
#Our validation check 

#       Create initial boom array

yspace=PF.boomspacing*abs(m.cos(PF.phi))
#print('yspace: ',yspace)
zspace=PF.boomspacing*abs(m.sin(PF.theta))
#print('zspace: ',zspace)
if (PF.xmin == PF.xmax):
   RAYMAX=int((PF.ymax-PF.ymin)/yspace)*int((PF.zmax-PF.zmin)/zspace)
elif(PF.ymin == PF.ymax):
   RAYMAX=int((PF.xmax-PF.xmin)/xspace)*int((PF.zmax-PF.zmin)/zspace)
elif(PF.zmin == PF.zmax):
   RAYMAX=int((PF.ymax-PF.ymin)/yspace)*int((PF.xmax-PF.xmin)/xspace)
boomarray = np.zeros((RAYMAX,2))
#print('boomarray: ', boomarray)
print(RAYMAX , ' is the RAYMAX')
boomarray,sizex,sizey,sizez=fun.InitialGrid(PF.boomspacing,PLANEABC[0],PLANEABC[1],PLANEABC[2],PLANEABC[3],PF.theta,PF.phi,PF.xmin,PF.ymin,PF.zmin,PF.xmax,PF.ymax,PF.zmax,RAYMAX)


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
RPS.receiverarray1=None
if RPS.planenum >=2 :
    RPS.receiverarray2=None
if RPS.planenum >=3 :
    RPS.receiverarray3=None
if RPS.planenum >=4 :
    RPS.receiverarray4=None
if RPS.planenum >=5 :
    RPS.receiverarray5=None
if RPS.planenum >=6 :
    RPS.receiverarray6=None
if RPS.planenum >=7 :
    RPS.receiverarray7=None
#############################################################

#       Initialize normalization factor 
normalization=(PI*radius2)/(PF.boomspacing**2)
temparray=np.zeros((RPS.arraysize,sizeffttwo,6))
timetemparray=np.zeros((RPS.arraysize,sizefft,5))
for D in range(1,RPS.arraysize):
    for W in range(1,int(sizeffttwo)):
        temparray[D,W,0]=RPS.receiverarray[D,0]
        temparray[D,W,1]=RPS.receiverarray[D,1]
        temparray[D,W,2]=RPS.receiverarray[D,2]
        temparray[D,W,3]=inputarray[W,0]
        temparray[D,W,4]=0.0 
        temparray[D,W,5]=0.0 

#       Define ground plane

groundheight=0.000000000
GROUNDABC=[0.000000000,0.000000000,1.00000000]
GROUNDD=-groundheight
nground=np.array([0.0,0.0,1.0])
alphaground=np.zeros(sizeffttwo)
#     Allocate absorption coefficients for each surface for each frequency

for D in range(1,sizeffttwo):
    if inputarray[D,1] >= 0.0 or inputarray[D,1] < 88.0 :
        alphaground[D]=PF.tempalphaground[1]
    elif inputarray[D,1] >= 88.0 or inputarray[D,1] < 177.0 :
        alphaground[D]=PF.tempalphaground[2]
    elif inputarray[D,1] >= 177.0 or inputarray[D,1] < 355.0 :
        alphaground[D]=PF.tempalphaground[3]
    elif inputarray[D,1] >= 355.0 or inputarray[D,1] < 710.0 :
        alphaground[D]=PF.tempalphaground[4]
    elif inputarray[D,1] >= 710.0 or inputarray[D,1] < 1420.0 :
        alphaground[D]=PF.tempalphaground[5]
    elif inputarray[D,1] >= 1420.0 or inputarray[D,1] < 2840.0 :
        alphaground[D]=PF.tempalphaground[6]
    elif inputarray[D,1] >= 2840.0 or inputarray[D,1] < 5680.0 :
        alphaground[D]=PF.tempalphaground[7]
    elif inputarray[D,1] >= 5680.0 or inputarray[D,1] < inputarray[sizeffttwo,1]:
        alphaground[D]=PF.tempalphaground[8]

alphabuilding = np.zeros((PF.absorbplanes,sizeffttwo))

for W in range(1,PF.absorbplanes):
    for D in range(1,PF.absorbplanes):
        if  inputarray[D,1] >= 0.0 or inputarray[D,1] < 88.0:
            alphabuilding[W,D]=PF.tempalphabuilding[W,1]
        elif inputarray[D,1] >= 88.0 or inputarray[D,1] < 177.0:
            alphabuilding[W,D] = PF.tempalphabuilding [W,2]


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
if PF.radiosity == 1 :
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

#     Loop through the intial ray locations
ray=0
while ray <= RAYMAX:
      hitcount=0
      tmpsum=0.0
      doublehit=0
      W=0
      while W < sizeffttwo:

            ampinitial[W]=inputarray[W,0]/normalization
            phaseinitial[W]=inputarray[W,1]
            W+=1
      Vinitial=[boomarray[ray,0],boomarray[ray,1],boomarray[ray,2]]
      if (PF.h < (2*PF.radius)): 
            print('h is less than 2r')
            break
      F=Finitial
      print('F1',F)
      veci=Vinitial
      #print(Vinitial, ' is Vinitial')
      #print('so veci is ', veci)
# Making small steps along the ray path.  For each step we should return, 
# location, phase and amplitude
      I=0
      while I < PF.IMAX:
            dxreceiver=HUGE
# Find the closest sphere and store that as the distance
#           print(veci)
            Q=0
            while (Q < RPS.arraysize):
                  tempreceiver=fun.SPHERECHECK(RPS.receiverarray[Q],radius2,F,veci)
                  if (receiverhit > 1):
                        if(lastreceiver[0]==RPS.receiverarray[Q,0] and lastreceiver[1]==receiverarray[Q,1] and lastreceiver[2] == receiverarray[Q,2]):
                              tempreceiver=HUGE
                        if(F[0] == checkdirection[0] and F[1] == checkdirection[1] and F[2] == checkdirection[2]):
                              OC[0]=RPS.receiverarray[Q,0]-veci[0]
                              OC[1]=RPS.receiverarray[Q,1]-veci[1]
                              OC[2]=RPS.receiverarray[Q,2]-veci[2] 
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
 #                       print('receivercheck',receivercheck)
                        if(receiverarray[Q,0]==receiverpoint[0] and receiverarray[Q,1]==receiverpoint[1] and receiverarray[Q,2]==receiverpoint[2]):
                              doublehit=0
                        else:
                              receiverpoint2[0]=receiverarray[Q,0]
                              receiverpoint2[1]=receiverarray[Q,1]
                              receiverpoint2[2]=receiverarray[Q,2]
                              doublehit=1
                  Q += 1 
#     Check Intersection with ground plane
            GROUNDN=GROUNDABC
            GROUNDVD=GROUNDN[0]*F[0]+GROUNDN[1]*F[1]+GROUNDN[2]*F[2]
            if (groundhit==1):
                  dxground=HUGE
            elif (GROUNDVD!=0.0):
                  GROUNDVO=((GROUNDN[0]*veci[0]+GROUNDN[1]*veci[1]+GROUNDN[2]*veci[2])+GROUNDD)
                  dxground1=(-1.000)*GROUNDVO*(1.000)/GROUNDVD
                  dxground=dxground1
                  Vecip1=veci+dxground*np.array(F)
                  tmp=(GROUNDABC[0]*Vecip1[0]+GROUNDABC[1]*Vecip1[1]+GROUNDABC[2]*Vecip1[2]+GROUNDD)                  
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
            while (Q < BG.Boxnumber):
                  dxnear, dxfar, hit, planehit=fun.BOX(BG.Boxarraynear[Q], BG.Boxarrayfar[Q],veci,F)
                  if (dxnear < dxbuilding):
                        dxbuilding=dxnear
                        Vecip1=(veci+dxbuilding)*np.array(F)
                        #print('debug for vecip1 ',veci,dxbuilding,np.array(F))
                        whichbox=Q
                        #print('Chunk of inputs coming ',Vecip1, BG.Boxarraynear[whichbox],BG.Boxarrayfar[whichbox], planehit)
                        nbox=fun.PLANE(Vecip1, BG.Boxarraynear[whichbox],BG.Boxarrayfar[whichbox], planehit)
                        #print('nbox from Boxnumber ', nbox)
                  Q+=1
#     Check intersection with Triangles
            if(BG.TriangleNumber > 0):
                  Q=0
                  while (Q < BG.TriangleNumber):
                        dxnear, behind = fun.Polygon(veci,F,Q,3,TriangleNumber,PointNumbers,Trianglearray,BuildingPoints,normal,FaceNormalNo,FaceNormals)
                        if (dxnear < dxbuilding):
                              dxbuilding=dxnear
                              nbox=normal
                              #print('nbox from triangles',nbox)
                              whichbox=Q
                        Q+=1
#    Check intersection with Squares
            if(BG.SquareNumber>0):
                  Q=0
                  while (Q < BG.SquareNumber):
                        dxnear, behind=Polygon(veci,F,Q,4,SquareNumber,PointNumbers,SquareArray,BuildingPoints,normal,FaceNormalNo,FaceNormals)
                        if (dxnear < dxbuilding):
                              dxbuilding=dxnear
                              nbox=normal
                              #print('nbox from squares', nbox)
                              whichbox=Q
                        Q+=1
            buildinghit=0
            receiverhit=0
            groundhit=0
#     Check to see if ray hits within step size
            if (dxreceiver < PF.h or dxground < PF.h or dxbuilding < PF.h):
                  dx=min(dxreceiver,dxground,dxbuilding)
                  tmpsum=tmpsum+dx
#     if the ray hits a receiver, store in an array.  If the ray hits twice
#     Create two arrays to store in.
                  if (dx==dxreceiver):
                        sum=sum+1
                        Vecip1=veci+dx*F
                        veci=Vecip1
                        receiverhit=1
                        checkdirection=F
                        if(doublehit==1):
                              receiverhit=2
                        hitcount=hitcount+1
                        #print('hit receiver',sum,tmpsum,receiverpoint)
                        W=0
                        while (W < sizeffttwo):
                              m=airabsorb(W)
                              lamb=soundspeed/inputarray[W,1]
                              phasefinal=phaseinitial[W]-(twopi*dx)/lamb   
                              ampfinal=ampinitial[W]*(1-alphanothing[W])*np.exp(-m*dx)
                              ampinitial[W]=ampfinal
                              phaseinitial[W]=phasefinal%twopi
                              if (phaseinitial[W]>=PI):
                                    phaseinitial[W]=phaseinitial[W]-twopi
                              if(doublehit==1):
                                    if(W==1):
                                          outputarray1=np.zeros(sizeffttwo, 6)
                                    if(W==1):
                                          dhoutputarray1=np.zeros(sizeffttwo,6)  
                                    outputarray1[W,0]=inputarray[W,0]
                                    outputarray1[W,1]=receiverpoint[0]
                                    outputarray1[W,2]=receiverpoint[1]
                                    outputarray1[W,3]=receiverpoint[2]
                                    outputarray1[W,4]=ampinitial[W]/2.0
                                    outputarray1[W,5]=phaseinitial[W]
                                    dhoutputarray1[W,0]=inputarray[W,0]
                                    dhoutputarray1[W,1]=receiverpoint2[0]
                                    dhoutputarray1[W,2]=receiverpoint2[1]
                                    dhoutputarray1[W,3]=receiverpoint2[2]
                                    dhoutputarray1[W,4]=ampinitial[W]/2.0
                                    dhoutputarray1[W,5]=phaseinitial[W]
                                    lastreceiver[0]=receiverpoint[0]
                                    lastreceiver[1]=receiverpoint[1]
                                    lastreceiver[2]=receiverpoint[2]
                                    lastreceiver2[0]=receiverpoint2[0]
                                    lastreceiver2[1]=receiverpoint2[1]
                                    lastrecever2[2]=receiverpoint2[2]
                              else:
 #                                   print('this happens outputarray')
                                    if(W==1):
                                          outputarray1=np.zeros(sizeffttwo,6)
                                    outputarray1[W,0]=inputarray(W,1)
                                    outputarray1[W,1]=receiverpoint(1)
                                    outputarray1[W,2]=receiverpoint(2)
                                    outputarray1[W,3]=receiverpoint(3)
                                    outputarray1[W,4]=ampinitial(W)
                                    outputarray1[W,5]=phaseinitial(W)
                                    lastreceiver[0]=receiverpoint[0]
                                    lastreceiver[1]=receiverpoint[1]
                                    lastreceiver[2]=receiverpoint[2]
 #                       print('assigned to outputarray1')
 #                       print(outputarray1[0,1], outputarray1[1,3], outputarray1[1,4])
                        temparray=fun.receiverHITFUNC(sizefft,outputarray1,arraysize)
 #                       print('receiverHITFUNC completed')
                        if (doublehit==1):
                              temparray=fun.receiverHITFUNC(sizefft,dhoutputarray1,arraysize,temparray)
                              count+=1
                        count+=1
                        outputarray1=None
                        if(doublehit==1):
                              dhoutputarray1=None
 #                       print('got to the end of receiver hit')
#     If the ray hits the ground then bounce off the ground and continue
                  if (abs(dx-dxground)< 10.0**(-13.0)):
                        Vecip1=veci+dxground*F
                        print('F is ', F)
                        tmp=(GROUNDABC[0]*Vecip1[0]+GROUNDABC[1]*Vecip1[1]+GROUNDABC[2]*Vecip1[2]+GROUNDD)
                        if(tmp != GROUNDD): 
                              Vecip1[2]=0.0
  #                      print('hit ground')
                        veci=Vecip1
  #                      print 'nground', nground[0]
                        dot1=(F[0]*nground[0]+F[1]*nground[1]+F[2]*nground[2])
                        n2=(nground[0]*nground[0]+nground[1]*nground[1]+nground[2]*nground[2])
                        #print('nground[0]', nground[0])
                        #print('nground[1]', nground[1])
                        #print('nground[2]', nground[2])                #Good
  #                      print (nground)
                        r=np.around(F-2.0*(dot1/n2)*nground,8)
                        length=np.sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])
                        F=np.around([r[0],r[1],r[2]],8)
                        print('F is now ', F)
                        print('r as follows: ', r[0],r[1],r[2])
                        groundhit=1
                        twopidx=twopi*dxground
#     Loop through all the frequencies
                        W=0
                        while (W<sizeffttwo):
                              m=airabsorb[W]
                              lamb=PF.soundspeed/inputarray[W,0]
                              phasefinal=phaseinitial[W]-(twopidx)/lamb
                              ampfinal=ampinitial[W]*(1.0-alphaground[W])*(1.0-diffusionground)*np.exp(-m*dxground)
                              phaseinitial[W]=phasefinal%twopi
                              if (phaseinitial[W]>PI):
                                    phaseinitial[W]=phaseinitial[W]-twopi
                              if(PF.radiosity==1 and (diffusionground!=0.0)):
                                    Q=0
                                    while (Q<PatchNo):
                                          if (formfactors[0,Q,1]==1):
                                                if(veci[0]<=(patcharray[Q,W,0]+0.5*patcharray[Q,W,3]) and veci[0]>=(patcharray[Q,W,0]-0.5*patcharray[Q,W,3])):
                                                      if(veci[1]<=(patcharray[Q,W,1]+0.5*patcharray[Q,W,4]) and veci[1]>=(patcharray[Q,W,1]-0.5*patcharray[Q,W,4])):
                                                            if(veci[2]<=(patcharray[Q,W,2]+0.5*patcharray[Q,W,5]) and veci[2]>=(patcharray[Q,W,2]-0.5*patcharray[Q,W,5])):
                                                                  temp2=complex(abs(patcharray[Q,W,6])*np.exp(XJ*patcharray[Q,W,7]))
                                                                  temp3=complex(abs(ampinitial[W]*(1.0-alphaground[W])*diffusionground*exp(-m*dxground))*exp(1j*phasefinal))
                                                                  temp4=temp2+temp3
                                                                  patcharray[Q,W,6]=abs(temp4)
                                                                  patcharray[Q,W,7]=np.arctan(temp4.imag,temp4.real)
                                          Q+=1
                              ampinitial[W]=ampfinal   
                              W+=1                        
               
#     if the ray hits the building then change the direction and continue
                  if (dx==dxbuilding):
                        Vecip1=veci+dx*np.array(F)
                        veci=Vecip1
                        #print('hit building')
                        n2=(nbox[0]*nbox[0]+nbox[1]*nbox[1]+nbox[2]*nbox[2])
                        #print('nbox for nbuilding: ',nbox)
                        #print('n2', n2)
                        #print('n2: ', n2)
                        nbuilding=nbox/np.sqrt(n2)
                        dot1=(F[0]*nbuilding[0]+F[1]*nbuilding[1]+F[2]*nbuilding[2])
                        r=F-2.0*(dot1/n2)*nbuilding
                        length=np.sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])
                        F=[r[0],r[1],r[2]]
                        buildinghit=1
                        twopidx=twopi*dx
                        W=0
                        while (W<sizeffttwo):
                              if(PF.complexabsorption==1):
                                    if (PF.absorbplanes==2):
                                          if(veci[2]>0.0 and veci[2]< height1):
                                                alpha=alphabuilding[0,W]
                                          elif(veci[2]>height1 and veci[2]<=height2):
                                                alpha=alphabuilding[1,W]
                                    if(PF.absorbplanes==3):
                                          if(veci[2]>height2 and veci[2] <=height3):
                                                alpha=alphabuilding[2,W]
                                    if(PF.absorbplanes==4):
                                          if(veci[2]>height3):
                                                alpha=alphabuilding(4,W)
                              else:
                                    alpha=alphabuilding[0,W]
                              m=airabsorb[W]
                              lamb=PF.soundspeed/inputarray[W,0]                
                              phasefinal=phaseinitial[W]-(twopidx)/lamb
                              ampfinal=ampinitial[W]*(1.0-alpha)*(1.0-diffusion)*np.exp(-m*dx)
                              phaseinitial[W]=phasefinal%twopi
                              if (phaseinitial[W]>PI):
                                    phaseinitial[W]=phaseinitial[W]-twopi
                              ampinitial[W]=ampfinal
                              W+=1
#COME BACK TO THIS FOR RADIOSITY
#                         if(radiosity.eq.1)then
# C     Loop through all patches if radiosity is turned on.
#                            DO 29 Q=1,PatchNo
#                               if (formfactors(1,Q,2).eq.2.0)then
#                                  if((veci(1).le.(patcharray(Q,W,1)+0.5*
#      *                                patcharray(Q,W,4))).and.(veci(1)
#      *                                .ge.(patcharray(Q,W,1)-0.5*
#      *                                patcharray(Q,W,4))))then
#                                     if((veci(2).le.(patcharray(Q,W,2)+
#      *                                   0.5*patcharray(Q,W,5))).and.
#      *                                   (veci(2).ge.(patcharray(Q,W,2)-
#      *                                   0.5*patcharray(Q,W,5))))then
#                                        if((veci(3).le.(patcharray(Q,W,3)
#      *                                      +0.5*patcharray(Q,W,6)))
#      *                                      .and.(veci(3).ge.(patcharray
#      *                                      (Q,W,3)-0.5*patcharray(Q,W,6
#      *                                      ))))then
#                                           temp2=cmplx(abs(patcharray
#      *                                         (Q,W,7))*exp(XJ*
#      *                                         patcharray(Q,W,8)))
#                                           temp3=abs(ampinitial(W)*(1.0-
#      *                                         alpha)*diffusion*exp
#      *                                         (-m*dx))*exp(XJ*
#      *                                         phaseinitial(W))
#                                           temp4=temp2+temp3
#                                           patcharray(Q,W,7)=abs(temp4)
#                                           patcharray(Q,W,8)=ATAN2(
#      *                                         imagpart(temp4),realpart(
#      *                                         temp4))
#                                           GOTO 27
#                                        endif
#                                     endif
#                                  endif
#  27                              CONTINUE
#                               endif
#  29                     CONTINUE
#                      endif 
                              
            else:
#     If there was no interaction with buildings then proceed with one step. 
                  tmpsum=tmpsum+PF.h
                  Vecip1=veci+(PF.h)*np.array(F)
                  veci=Vecip1
                  twopih=twopi*PF.h
                  W=0
                  while (W<sizeffttwo):
#     Loop through all frequencies. 
                        m=airabsorb[W]
                        lamb=PF.soundspeed/inputarray[W,0]
                        phasefinal=phaseinitial[W]-(twopih)/lamb
                        ampfinal=ampinitial[W]*(1-alphanothing)*np.exp(-m*PF.h)
                        ampinitial[W]=ampfinal                  
                        phaseinitial[W]=phasefinal%twopi
                        if (phaseinitial[W] > PI):
                              phaseinitial[W]=phaseinitial[W]-twopi
                        W+=1
            I += 1
      print('finished ray', ray + 1)
      ray += 1
#     Once all rays are complete.  Deallocate all arrays that are no longer needed
#      boomarray=None
#      receiverarray=None
#      ampinitial=None
#      phaseinitial=None
      # Gk=np.zeros(PatchNo,sizeffttwo)
      # Gkminus1=np.zeros(PatchNo,sizeffttwo)
#COMEBACK FOR RADIOSITY
#       if (radiosity.eq.1)then
# C     If radiosity is turned on then do the energy exchange. 
#          KMAX=3
#          Npatch=10
#          DO 53 K=1,KMAX
#             DO 55 D=1,PatchNo
#                DO 54 W=1,sizeffttwo
#                   if (K.eq.1)then
#                      Gk(D,W)=cmplx(abs(patcharray(D,W,7)
#      *                    )*exp(XJ*patcharray(D,W,8)))
#                   else
#                      Gkminus1(D,W)=Gk(D,W)
#                      DO 56 I=1,PatchNo
#                         if (I.eq.1)then
#                            Gk(D,W)=0.0
#                         else
#                            if (formfactors(D,I,2).eq.1.0)then 
#                               alpha=alphaground(W)
#                            elseif (formfactors(D,I,2).eq.2.0)then
#                               if( complexabsorption.eq.1)then
#                                  if(absorbplanes.eq.2)then
#                                  if(patcharray(I,1,3).gt.0.0.and.
#      *                                patcharray(I,1,3).le.height1)then
#                                     alpha=alphabuilding(1,W)
#                                  elseif(patcharray(I,1,3).gt.height1
#      *                                   .and.patcharray(I,1,3).le.
#      *                                   height2)then
#                                     alpha=alphabuilding(2,W)
#                                  endif
#                                  endif
#                                  if(absorbplanes.eq.3)then
#                                  if(patcharray(I,1,3).gt.height2
#      *                                   .and.patcharray(I,1,3).le.
#      *                                   height3)then
#                                     alpha=alphabuilding(3,W)
#                                  endif
#                                  endif 
#                                  if(absorbplanes.eq.4)then
#                                  if(patcharray(I,1,3).gt.height3)
#      *                                   then
#                                     alpha=alphabuilding(4,W)
#                                  endif
#                                  endif
#                               else
#                                  alpha=alphabuilding(1,W)
#                               endif
#                            endif
#                            m=airabsorb(W)
#                            temp2=(1-alpha)*exp(-m*formfactors(D,I,3))*
#      *                          formfactors(D,I,1)*Gkminus1(D,W)*exp(
#      *                          -XJ*twopi*inputarray(W,1)*formfactors
#      *                          (D,I,3)/soundspeed)
#                            Gk(D,W)=Gk(D,W)+temp2
#                         endif
#  56                  CONTINUE
#                   endif
#  54            CONTINUE
#                print*, 'finished patch', D, 'of',PatchNo
#  55         CONTINUE
#             print*, arraysize,PatchNo,sizefft
# C     Do energy exchange with other receivers
#             DO 50 D=1,arraysize
#                DO 51 Q=1,PatchNo
#                   Rlm=0.0
#                   DO 58 I=1,Npatch
#                      DO 59 J=1,Npatch
#                         DO 60 S=1,Npatch
#                            tmp1=((patcharray(Q,1,1)-.5*patcharray
#      *                          (Q,1,4)+(patcharray(Q,1,4)/Npatch)
#      *                          *(I-.5))-temparray(D,1,1))*((patcharray(
#      *                          Q,1,1)-.5*patcharray(Q,1,4)+(patcharray(
#      *                          Q,1,4)/Npatch)*(I-.5))-temparray(D,1,1))
#                            tmp2=((patcharray(Q,1,2)-.5*patcharray
#      *                          (Q,1,5)+(patcharray(Q,1,5)/Npatch)
#      *                          *(J-.5))-temparray(D,1,2))*((patcharray(
#      *                          Q,1,2)-.5*patcharray(Q,1,5)+(patcharray(
#      *                          Q,1,5)/Npatch)*(J-.5))-temparray(D,1,2))
#                            tmp3=((patcharray(Q,1,3)-.5*patcharray
#      *                          (Q,1,6)+(patcharray(Q,1,6)/Npatch)
#      *                          *(S-.5))-temparray(D,1,3))*((patcharray(
#      *                          Q,1,3)-.5*patcharray(Q,1,6)+(patcharray(
#      *                          Q,1,6)/Npatch)*(S-.5))-temparray(D,1,3))
#                            Rlm=Rlm+1.0/(NPatch*Npatch*Npatch)*sqrt(tmp1+
#      *                          tmp2+tmp3)
#  60                     Continue
#  59                  CONTINUE
#  58               CONTINUE
#                   Patchlength=sqrt(((patcharray(Q,1,1)-temparray(D,1,1))
#      *                 *(patcharray(Q,1,1)-temparray(D,1,1)))+((patchar
#      *                 ray(Q,1,2)-temparray(D,1,2))*(patcharray(Q,1,2)
#      *                 -temparray(D,1,2)))+((patcharray(Q,1,3)-temparray
#      *                (D,1,3))*(patcharray(Q,1,3)-temparray
#      *                (D,1,3))))
#                   Finitial(1)=(patcharray(Q,1,1)-temparray(D,1,1))
#                   Finitial(2)=(patcharray(Q,1,2)-temparray(D,1,2))
#                   Finitial(3)=(patcharray(Q,1,3)-temparray(D,1,3))
#                   F=(/Finitial(1)/Rlm,Finitial(2)/Rlm,
#      *                 Finitial(3)/Rlm/)
#                   dxbuilding=HUGE
# C     Check to see that the receiver is visible by patches. 
#                   DO 57 I=1,boxnumber,1 
#                      call BOX(boxarraynear(I,1:3),
#      *                    Boxarrayfar(I,1:3),temparray(D,1,1:3),
#      *                    F,dxnear,dxfar,hit, planehit)
#                      if (dxnear.lt.dxbuilding)then
#                         dxbuilding=dxnear
#                      endif
#                      if ((temparray(D,1,1).gt.boxarraynear(I,1)
#      *                    .and.temparray(D,1,2).gt.boxarraynear(I,2)
#      *                    .and.temparray(D,1,3).gt.boxarraynear(I,3))
#      *                    .and.(temparray(D,1,1).lt.boxarrayfar(I,1)
#      *                    .and.temparray(D,1,2).lt.boxarrayfar(I,2)
#      *                    .and.temparray(D,1,3).lt.boxarrayfar(I,3)))
#      *                    then
#                         dxbuilding=-1.0
#                      endif
#  57               CONTINUE
#                   if(TriangleNumber.gt.0)then
#                      DO 67 I=1,TriangleNumber,1 
#                         call Polygon(temparray(D,1,1:3),F,I,3,
#      *                       TriangleNumber,PointNumbers,Trianglearray,
#      *                       BuildingPoints,normal,FaceNormalNo,
#      *                       FaceNormals,dxnear,behind)
#                         if (dxnear.lt.dxbuilding)then
#                            dxbuilding=dxnear
#                         endif
#  67                  CONTINUE
#                   endif
#                   if(SquareNumber.gt.0)then
#                      DO 68 I=1,SquareNumber,1 
#                         call Polygon(temparray(D,1,1:3),F,I,4,
#      *                       SquareNumber,PointNumbers,SquareArray,
#      *                       BuildingPoints,normal,FaceNormalNo,
#      *                       FaceNormals,dxnear,behind)
#                         if (dxnear.lt.dxbuilding)then
#                            dxbuilding=dxnear
#                         endif
#  68                  CONTINUE
#                   endif
#                   if(PolyBuilding.gt.0)then
#                      DO 69 I=1,PolyBuilding
# C     Check that the receivers are not inside the building
#                         CALL INSIDECHECK(temparray(D,1,1:3),
#      *                       BuildingPoints,PointNumbers,F,Trianglearray
#      *                       ,SquareArray,TriangleNumber,SquareNumber,
#      *                       TriangleSequence,SquareSequence,Triangles,
#      *                       Squares,PolyBuilding,I,inside,FaceNormalNo,
#      *                       FaceNormals)
#                         if(inside.eq.1)then
#                            dxbuilding=-1.0
#                         endif
#  69                  CONTINUE
#                   endif
#                   vec3=FaceNormals(int(patcharray(Q,1,10)),1:3)
#                   PIRlm2=PI*Rlm*Rlm
#                   DO 52 W=1,sizeffttwo
#                      if (Gk(Q,W).ne.0.0)then   
#                         if(dxbuilding.lt.Patchlength)then
#                            Ek=0.0
#                         elseif(dxbuilding.ge.Patchlength)then
#                            length=sqrt(vec3(1)*vec3(1)+
#      *                          vec3(2)*vec3(2)+vec3(3)*vec3(3))
#                            cosxilm=(vec3(1)*F(1)+vec3(2)*F(2)+vec3(3)*
#      *                          F(3))/(Rlm*length)
#                            m=airabsorb(W)
#                            Ek=exp(-m*Rlm)*((cosxilm*Gk(Q,W)*exp(-XJ*
#      *                          twopi*inputarray(W,1)*Rlm/soundspeed))/
#      *                          (PIRlm2))
#                         endif
#                         temp2=cmplx(abs(temparray(D,W,5))*exp(XJ*
#      *                       temparray(D,W,6)))
#                         temp3=Ek+temp2
#                         temparray(D,W,5)=ABS(temp3)
#                         temparray(D,W,6)=ATAN2(imagpart(temp3)
#      *                       ,realpart(temp3))
#                      endif
#  52               CONTINUE
#  51            CONTINUE
#                print*, 'finished receiver', D, 'of', arraysize
#  50         CONTINUE
#  53      CONTINUE
#       endif

#Reconstruct the time signal
timetemparray=fun.TIMERECONSTRUCT(sizefft, timearray, RPS.arraysize, temparray)
#     Write out time signatures for each receiver. 
OPFile=open(PF.OUTPUTFILE,"w")
true=fun.Header(PF.OUTPUTFILE)
#OPFile=open(PF.OUTPUTFILE,"a")
if(RPS.planenum>=1):
      W=0
      while (W<sizefft):
            true=fun.TimeHeader(OPFile,timetemparray[0,W,3],RPS.sizex1,RPS.sizey1,RPS.sizez1,RPS.planename1)
            D=0
            while (D<RPS.arraysize1):
                  OPFile.write('%d,%d,%d,%d'%(timetemparray[D,W,0],timetemparray[D,W,1],timetemparray[D,W,2],timetemparray[D,W,4]))
                  D+=1
 #                print('finished time',timetemparray[0,W,3])
            W+=1
                  
if(RPS.planenum>=2):
      W=0
      while (W<sizefft):
            true=fun.TimeHeader(OPFile,timetemparray[0,W,3],sizex2,sizey2,sizez2,planename2)
            D=RPS.arraysize1
            while (D<RPS.arraysize1+RPS.arraysize2):
                  OPFile.write(timetemparray[D,W,0:2],timetemparray[D,W,4])
                  D+=1
            W+=1
if(RPS.planenum>=3):
      W=0
      while (W<sizefft):
            true=fun.TimeHeader(OPFile,timetemparray[0,W,3],sizex3,sizey3,sizez3,planename3)
            D=arraysize1+arraysize2
            while (D<(arraysize1+arraysize2+arraysize3)):
                  OPFile.write(timetemparray[D,W,0:2],timetemparray[D,W,4])
                  D+=1
            W+=1
if(RPS.planenum>=4):
      W=0
      while (W<sizefft):
            true=fun.TimeHeader(OPFile,timetemparray[0,W,3],sizex4,sizey4,sizez4,planename4)
            D=arraysize1+arraysize2+arraysize3
            while (D<(arraysize1+arraysize2+arraysize3+arraysize4)):
                  OPFile.write(timetemparray[D,W,0:2],timetemparray[D,W,4])
                  D+=1
            W+=1
if(RPS.planenum>=5):
      W=0
      while (W<sizefft):
            true=fun.TimeHeader(OPFile,timetemparray[0,W,3],sizex5,sizey5,sizez5,planename5)
            D=arraysize1+arraysize2+arraysize3+arraysize4
            while(D<(arraysize1+arraysize2+arraysize3+arraysize4+arraysize5)):
                  OPFile.write(timetemparray[D,W,0:2],timetemparray[D,W,4])
                  D+=1
            W+=1
if(RPS.planenum>=6):
      W=0
      while (W<sizefft):
            true=fun.TimeHeader(OPFile,timetemparray[0,W,3],sizex6,sizey6,sizez6,planename6)
            D=arraysize1+arraysize2+arraysize3+arraysize4+arraysize5
            while (D<arraysize1+arraysize2+arraysize3+arraysize4+arraysize5+arraysize6):
                  OPFile.write(timetemparray[D,W,0:2],timetemparray[D,W,4])
                  D+=1
            W+=1
if(RPS.planenum>=7):
      W=0
      while (W<sizefft):
            true=fun.TimeHeader(OPFile,timetemparray[1,W,4],sizex7,sizey7,sizez7,planename7)
            D=arraysize1+arraysize2+arraysize3+arraysize4+arraysize5+arraysize6
            while (D<arraysize1+arraysize2+arraysize3+arraysize4+arraysize5+arraysize6+arraysize7):
                  OPFile.write(timetemparray[D,W,0:2],timetemparray[D,W,4])
                  D+=1
            W+=1
OPFile.close()