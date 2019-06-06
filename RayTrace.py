# RayTrace
# version 1.0.10

# Kimberly Lefkowitz created this program to propagate sonic booms around
# large structures, and to graduate. It is a ray tracing model that 
# will include specular and diffuse reflections. It will print out the
# sound field at ear height, at relevant microphone locations, and at 
# the building walls. It will read in the fft of a sonic boom signature.

# Dr. Riegel, William Costa, and George Seaton porting program from Fortran to python

# Initialize variables and functions
import numpy as np

import Parameterfile as PF
import BuildingGeometry as BG
import Functions as fun
import ReceiverPointSource as RPS
#import GeometryParser as BG

import time
import memory_profiler as mem 

t = time.time()
      
# What it does not do
"""
      Interacts with geometry parser
      Have a way of reading in complex geometries - Yes, but not yet integrated
      Anything resembling radiosity
"""

def initial_signal(signalLength,fftOutput):
    """
    Making the array for the initial signals.
    Input sizefft and outputsignal
    """
    signalLength2 = int(signalLength//2)    # Making sizeffttwo in this function and setting it as an int again jsut to be sure
    outputFrequency = np.zeros((signalLength2,3))    #Making output array    equivalent to inputarray in old code
    throwarray = np.arange(1,signalLength2 + 1)     #Helps get rid of for-loops in old version

    outputFrequency[:,0] = throwarray * PF.Fs / signalLength #Tried simplifying the math a bit from original
    outputFrequency[:,1] = abs(fftOutput[1:1+signalLength2]/signalLength) #Only go up to sizeffttwo 
    outputFrequency[:,2] = np.arctan2( np.imag(fftOutput[1:1+signalLength2]/signalLength) , np.real(fftOutput[1:1+signalLength2]/signalLength)) 

    return outputFrequency

def updateFreq(dx,alpha,diffusion):
      """
      Update ray phase and amplitude
      """
      global phase,amplitude        # works directly

      twopidx = twopi * dx
      ein  = phase - (twopidx/lamb)
      zwei = ein % twopi
      masque = zwei> PI
      drei = masque * zwei - twopi
 
      phase =np.where(masque,drei,ein)
      amplitude *= ((1.0-alpha) * (1.0-diffusion) * np.exp(-airabsorb*dx))

def vex(y,z):
    """The x coordinate of the ray 
    Used for veci"""
    return (D-Finitial[1]*y-Finitial[2]*z)/Finitial[0]

# port and import receiver file
receiverhit=0
groundhit=0

# Initialize counters 
PI = np.pi
twopi = PI*2
XJ=complex(0.0,1.0)
radius2 = PF.radius**2
raysum=0

# Initiailize receiver variables
lastreceiver = np.zeros(3)
lastreceiver2 = np.zeros(3)
receiverpoint  = np.zeros(3)
receiverpoint2 = np.zeros(3)
#OC = np.empty(3)

# Read in input file
with open(PF.INPUTFILE) as IPFile:
      inputsignal=np.loadtxt(IPFile)
K=len(inputsignal)
masque = inputsignal>0
HUGE=1000000.0

# Allocate the correct size to the signal and fft arrays
sizefft=K
sizeffttwo=sizefft//2
outputsignal=np.fft.rfft(inputsignal,sizefft)

#       Create initial signal 
frecuencias = initial_signal(sizefft,outputsignal)      # Equivalent to inputarray in original
airabsorb=fun.ABSORPTION(PF.ps,frecuencias[:,0],PF.hr,PF.Temp)        #sizeffttwo
lamb = PF.soundspeed/frecuencias[:,0]     # Used for updating frequencies in update function
timearray = np.arange(K) /PF.Fs

#       Set initial values
Vinitial =  np.array([PF.xinitial,PF.yinitial,PF.zinitial])
xiinitial  =np.cos(PF.phi)*np.sin(PF.theta)
ninitial   =np.sin(PF.phi)*np.sin(PF.theta)
zetainitial=np.cos(PF.theta)
length =    np.sqrt(xiinitial*xiinitial+ninitial*ninitial+zetainitial*zetainitial)
Finitial=np.array([xiinitial,ninitial,zetainitial])
#tmp=(Finitial[0]*Vinitial[0]+Finitial[1]*Vinitial[1]+Finitial[2]*Vinitial[2])
D = np.dot(Finitial,Vinitial)   #equivalent to tmp

#PLANEABC=np.array([Finitial[0],Finitial[1],Finitial[2],tmp])

#       Create initial boom array
#yspace=PF.boomspacing*abs(np.cos(PF.phi))
#zspace=PF.boomspacing*abs(np.sin(PF.theta))
#if (PF.xmin == PF.xmax):
#   RAYMAX=int((PF.ymax-PF.ymin)/yspace)*int((PF.zmax-PF.zmin)/zspace)
#elif(PF.ymin == PF.ymax):
#   RAYMAX=int((PF.xmax-PF.xmin)/xspace)*int((PF.zmax-PF.zmin)/zspace)
#elif(PF.zmin == PF.zmax):
#   RAYMAX=int((PF.ymax-PF.ymin)/yspace)*int((PF.xmax-PF.xmin)/xspace)
#print(RAYMAX , ' is the RAYMAX')
#boomarray,sizex,sizey,sizez=fun.InitialGrid(PF.boomspacing,PLANEABC[0],PLANEABC[1],PLANEABC[2],PLANEABC[3],PF.theta,PF.phi,PF.xmin,PF.ymin,PF.zmin,PF.xmax,PF.ymax,PF.zmax,RAYMAX)

#x = (PF.xmin,PF.xmax) 
#y = (PF.ymin,PF.ymax) 
#z = (PF.zmin,PF.zmax) 
#
#yspace=PF.boomspacing*abs(np.cos(PF.phi))
#zspace=PF.boomspacing*abs(np.sin(PF.theta))
#if (PF.xmin == PF.xmax):
#   raymax=int((PF.ymax-PF.ymin)/yspace)*int((PF.zmax-PF.zmin)/zspace)
#print(raymax , ' is the raymax')
#
#j = np.arange(1,1+int((y[1]-y[0])//yspace))
#k = np.arange(1,1+int((z[1]-z[0])//zspace))
#
#rayy = y[0] + j*yspace
#rayz = z[0] + k*zspace

yspace=PF.boomspacing*abs(np.cos(PF.phi))
zspace=PF.boomspacing*abs(np.sin(PF.theta))
if (PF.xmin == PF.xmax):
   raymax=int((PF.ymax-PF.ymin)/yspace)*int((PF.zmax-PF.zmin)/zspace)
print(raymax , ' is the raymax')

j = np.arange(1,1+int((PF.ymax-PF.ymin)//yspace))
k = np.arange(1,1+int((PF.zmax-PF.zmin)//zspace))

rayy = PF.ymin + j*yspace
rayz = PF.zmin + k*zspace

boomcarpet = ((vex(y,z),y,z) for z in rayz for y in rayy )

#for ray in boomcarpet:
#    # Positioning for rays along initial grid
#    veci = ray
    #pass

#     Create a receiver array, include a receiver file. 
alphanothing = np.zeros(sizeffttwo)

# Making specific receiver points using receiver module
RPS.Receiver.initialize(PF.RecInput)
ears = RPS.Receiver.rList           #easier to write
for R in ears:          #hotfix
      R.magnitude = np.zeros(sizeffttwo)
      R.direction = np.zeros(sizeffttwo)

#       Initialize normalization factor 
normalization=(PI*radius2)/(PF.boomspacing**2) 

outputarray1=np.zeros((sizeffttwo,6))
dhoutputarray1=np.zeros((sizeffttwo,6))  

#       Define ground plane
groundheight=0.000000000
GROUNDABC=np.array([0.000000000,0.000000000,1.00000000])
GROUNDD=-groundheight
groundD= -groundheight
nground=np.array([0.0,0.0,1.0])

#     Allocate absorption coefficients for each surface for each frequency
alphaground=np.zeros(sizeffttwo)
for D in range(0,sizeffttwo):       #This loop has a minimal impact on performance
    if   frecuencias[D,0] >= 0.0 or    frecuencias[D,0] < 88.0 :
        alphaground[D]=PF.tempalphaground[0]
    elif frecuencias[D,0] >= 88.0 or   frecuencias[D,0] < 177.0 :
        alphaground[D]=PF.tempalphaground[1]
    elif frecuencias[D,0] >= 177.0 or  frecuencias[D,0] < 355.0 :
        alphaground[D]=PF.tempalphaground[2]
    elif frecuencias[D,0] >= 355.0 or  frecuencias[D,0] < 710.0 :
        alphaground[D]=PF.tempalphaground[3]
    elif frecuencias[D,0] >= 710.0 or  frecuencias[D,0] < 1420.0 :
        alphaground[D]=PF.tempalphaground[4]
    elif frecuencias[D,0] >= 1420.0 or frecuencias[D,0] < 2840.0 :
        alphaground[D]=PF.tempalphaground[5]
    elif frecuencias[D,0] >= 2840.0 or frecuencias[D,0] < 5680.0 :
        alphaground[D]=PF.tempalphaground[6]
    elif frecuencias[D,0] >= 5680.0 or frecuencias[D,0] < frecuencias[sizeffttwo,0]:
        alphaground[D]=PF.tempalphaground[7]

alphabuilding = np.zeros((PF.absorbplanes,sizeffttwo))
for W in range(PF.absorbplanes):        #These also look minimal
    for D in range(sizeffttwo):
        if   frecuencias[D,0] >= 0.0   or  frecuencias[D,0] < 88.0:
            alphabuilding[W,D]    =    PF.tempalphabuilding[W,0]
        elif frecuencias[D,0] >= 88.0  or  frecuencias[D,0] < 177.0:
            alphabuilding[W,D]    =    PF.tempalphabuilding[W,1]
        elif frecuencias[D,0] >= 177.0 or  frecuencias[D,0] < 355.0 :
            alphabuilding[W,D]    =    PF.tempalphabuilding[W,2]
        elif frecuencias[D,0] >= 355.0 or  frecuencias[D,0] < 710.0 :
            alphabuilding[W,D]    =    PF.tempalphabuilding[W,3]
        elif frecuencias[D,0] >= 710.0 or  frecuencias[D,0] < 1420.0 :
            alphabuilding[W,D]    =    PF.tempalphabuilding[W,4]
        elif frecuencias[D,0] >= 1420.0 or frecuencias[D,0] < 2840.0 :
            alphabuilding[W,D]    =    PF.tempalphabuilding[W,5]
        elif frecuencias[D,0] >= 2840.0 or frecuencias[D,0] < 5680.0 :
            alphabuilding[W,D]    =    PF.tempalphabuilding[W,6]
        elif frecuencias[D,0] >= 5680.0 or frecuencias[D,0] < frecuencias[sizeffttwo,0] :
            alphabuilding[W,D]    =    PF.tempalphabuilding[W,7]

#        Mesh the patches for the environment.  Include patching file. 
diffusionground = 0.0
if PF.radiosity:  # If it exists as a non-zero number
      import SingleBuildingGeometry
      diffusion = PF.radiosity
else:
      diffusion = 0.0

# Begin tracing
#     Loop through the intial ray locations

#ray = 606                     # @ PF.boomspacing = 1
#ray = 455174                 # @ PF.boomspacing = 0.06
#if ray:                 #for debugging
#for ray in range(RAYMAX):
      #veci = boomarray[ray,:]                                     # Position
print('began rays')
raycounter = 0
for ray in boomcarpet:
      veci = ray      # initial ray position
      hitcount=0
      doublehit=0
      amplitude = frecuencias[:,1]/normalization
      #print('ini: ', list(amplitude[-5:]))
      phase=frecuencias[:,2]
      if (PF.h < (2*PF.radius)): 
            print('h is less than 2r')
      #      break
      F = np.array(Finitial)                                      # Direction
      for I in range(PF.IMAX):      # Making small steps along the ray path.  For each step we should return, location, phase and amplitude
            dxreceiver=HUGE
            # Find the closest sphere and store that as the distance
            for R in ears:
                  # The way that tempreceiver works now, it's only used here and only should be used here. It's not defined inside the receiver because it's ray dependant.
                  tempreceiver = R.SphereCheck(radius2,F,veci)    # Distrance to receiver
                  if (receiverhit >= 1):  #if you hit a receiver last time, don't hit it again
                        if np.all(R.position ==lastreceiver):
                              tempreceiver=HUGE
                        if np.all(F == checkdirection):
                              OC = R.position - veci
                              OCLength = np.dot(OC,OC)
                              if(OCLength < radius2):
                                    tempreceiver=HUGE
                  if(receiverhit >= 2):
                        if np.all(R.position == lastreceiver):
                              tempreceiver=HUGE
                  if (tempreceiver < dxreceiver):   
                        R.dxreceiver=tempreceiver
                        dxreceiver=tempreceiver
                        receiverpoint= R.position
                  elif (tempreceiver == dxreceiver and tempreceiver != HUGE):
                        receivercheck=tempreceiver          
                        if np.all(R.position==receiverpoint):
                              doublehit=0
                        else:
                              #receiverpoint2 = R.position
                              R2 = R
                              doublehit=1
                              print('double hit')

            #     Check Intersection with ground plane
            GROUNDN=GROUNDABC
            GROUNDVD=GROUNDN[0]*F[0]+GROUNDN[1]*F[1]+GROUNDN[2]*F[2]
            if (groundhit==1):
                  dxground=HUGE
            elif (GROUNDVD!=0.0):
                  GROUNDVO=((GROUNDN[0]*veci[0]+GROUNDN[1]*veci[1]+GROUNDN[2]*veci[2])+GROUNDD)
                  dxground1=(-1.000)*GROUNDVO*(1.000)/GROUNDVD
                  dxground=dxground1
                  Vecip1=veci+dxground*F
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
            for Q in range(0,BG.Boxnumber):
                  dxnear, dxfar, hit, planehit=fun.BOX(BG.Boxarraynear[Q], BG.Boxarrayfar[Q],veci,F)
                  if (dxnear < dxbuilding):
                        dxbuilding=dxnear
                        Vecip1=veci+np.multiply(dxbuilding,F)
                        whichbox=Q
                        nbox=fun.PLANE(Vecip1, BG.Boxarraynear[whichbox],BG.Boxarrayfar[whichbox], planehit)
            #     Check intersection with Triangles
            if(BG.TriangleNumber > 0):
                  for Q in range(0, BG.TriangleNumber):
                        dxnear, behind = fun.Polygon(veci,F,Q,3,TriangleNumber,PointNumbers,Trianglearray,BuildingPoints,normal,FaceNormalNo,FaceNormals)
                        if (dxnear < dxbuilding):
                              dxbuilding=dxnear
                              nbox=normal
                              whichbox=Q
            #     Check intersection with Squares
            if(BG.SquareNumber>0):
                  for Q in range(0,BG.SquareNumber):
                        dxnear, behind=Polygon(veci,F,Q,4,SquareNumber,PointNumbers,SquareArray,BuildingPoints,normal,FaceNormalNo,FaceNormals)
                        if (dxnear < dxbuilding):
                              dxbuilding=dxnear
                              nbox=normal
                              whichbox=Q
            buildinghit=0
            receiverhit=0
            groundhit=0

            #     Check to see if ray hits within step size
            if (dxreceiver < PF.h or dxground < PF.h or dxbuilding < PF.h):
                  dx=min(dxreceiver,dxground,dxbuilding)
                  #     if the ray hits a receiver, store in an array.  If the ray hits two, create two arrays to store in.
                  for R in ears:
                        if dx == R.dxreceiver:
                              #print('Ray ',ray +1,' hit receiver ',R.recNumber,' at step ',I)
                              print('Ray ',ray +1,' hit receiver ',R.recNumber)
                              veci += (dx*F)
                              receiverhit=1
                              checkdirection=F
                              if(doublehit==1):
                                    receiverhit=2
                              hitcount=hitcount+1
                              updateFreq(dx,alphanothing,0)
                              lastreceiver = receiverpoint
                              outputarray1[:,0] = frecuencias[:,0]
                              outputarray1[:,1:4] = receiverpoint[:]
                              outputarray1[:,5] = phase[:]
                              #print(list(ears[1].magnitude))
                              if doublehit == 1 :
                                    #R2 = R      #Supposed to be other R, but just a placeholder for now
                                    R.on_Hit(amplitude/2,phase)
                                    R2.on_Hit(amplitude/2,phase)
                              else:
                                    R.on_Hit(amplitude,phase)

                              #if(doublehit==1):
                              #      outputarray1[:,4]=amplitude[:]/2.0
                              #      dhoutputarray1[:,0]=inputarray[:,0]
                              #      dhoutputarray1[:,1:4]=receiverpoint2[:]
                              #      dhoutputarray1[:,4]=amplitude[:]/2.0
                              #      dhoutputarray1[:,5]=phase[:]
                              #      lastreceiver2 = receiverpoint2
                              #else:
                              #      outputarray1[:,4]=amplitude[:]
                              #temparray=fun.receiverHITFUNC(sizefft,outputarray1,RPS.arraysize,temparray)   # looks like it does the same thing as on_Hit. Here later
                              #R.on_Hit(amplitude,phase)
                              #if (doublehit==1):
                              #      #temparray=fun.receiverHITFUNC(sizefft,dhoutputarray1,RPS.arraysize,temparray) #Using objects may circumvent the need to have this, but it stays for now
                              #      count+=1
                              #count+=1
                  if (abs(dx-dxground)< 10.0**(-13.0)):                  #     If the ray hits the ground then bounce off the ground and continue
                        veci += (dxground*F)
                        tmp = np.dot(GROUNDABC,veci)
                        if(tmp != GROUNDD): 
                              veci[2] = 0
                        print('hit ground at ',I)
                        dot1 = np.dot(F,nground)
                        n2 = np.dot(nground,nground)
                        F -= (2.0*(dot1/n2 *nground))
                        length = np.sqrt(np.dot(F,F))
                        groundhit=1
                        twopidx=twopi*dxground
                        #     Loop through all the frequencies
                        updateFreq(dxground,alphaground,diffusionground)
                        if(PF.radiosity==1 and (diffusionground!=0.0)):
                              for Q in range (0,PatchNo):
                                    if (formfactors[0,Q,1]==1):
                                          if(veci[0]<=(patcharray[Q,W,0]+0.5*patcharray[Q,W,3]) and veci[0]>=(patcharray[Q,W,0]-0.5*patcharray[Q,W,3])):
                                                if(veci[1]<=(patcharray[Q,W,1]+0.5*patcharray[Q,W,4]) and veci[1]>=(patcharray[Q,W,1]-0.5*patcharray[Q,W,4])):
                                                      if(veci[2]<=(patcharray[Q,W,2]+0.5*patcharray[Q,W,5]) and veci[2]>=(patcharray[Q,W,2]-0.5*patcharray[Q,W,5])):
                                                            temp2=complex(abs(patcharray[Q,W,6])*np.exp(XJ*patcharray[Q,W,7]))
                                                            temp3=complex(abs(amplitude[W]*(1.0-alphaground[W])*diffusionground*exp(-m*dxground))*exp(1j*phasefinal))
                                                            temp4=temp2+temp3
                                                            patcharray[Q,W,6]=abs(temp4)
                                                            patcharray[Q,W,7]=np.arctan(temp4.imag,temp4.real)
                  if (dx==dxbuilding):                  #     if the ray hits the building then change the direction and continue
                        veci += (dx*F)
                        print('hit building at step ',I)
                        n2 = np.dot(nbox,nbox)
                        nbuilding=nbox/np.sqrt(n2)
                        dot1= np.dot(F,nbuilding)
                        F -= (2.0*(dot1/n2 *nbuilding))
                        length = np.sqrt(np.dot(F,F))
                        buildinghit=1
                        if PF.complexabsorption:
                              if PF.absorbplanes==2:
                                    if (veci[2]>0.0) and (veci[2]<height1):
                                          alpha = alphabuilding[0,:]
                                    elif(veci[2]>height1 and veci[2]<=height2):
                                          alpha=alphabuilding[1,:]
                              if(PF.absorbplanes==3):
                                    if(veci[2]>height2 and veci[2] <=height3):
                                          alpha=alphabuilding[2,:]
                              if(PF.absorbplanes==4):
                                    if(veci[2]>height3):
                                          alpha=alphabuilding[4,:]
                        else:
                              alpha=alphabuilding[0,:]
                        updateFreq(dx,alpha,diffusion) 
            else:     #     If there was no interaction with buildings then proceed with one step. 
                  veci += (PF.h*F)
                  updateFreq(PF.h,alphanothing,0)
      raycounter +=1
      print('finished ray', raycounter)

# Radiosity removed for readability


#Reconstruct the time signal and print to output file
for R in ears:
      R.timeReconstruct(sizefft)

#print(list(ears[1].magnitude))
#print(list(ears[1].direction))
#print(list(ears[1].signal))
      #TIMERECONSTRUCT(sizefft, timearray, arraysize, temparray, timetemparray)
#print(timearray[:5])   #magnitude, initial pressure,direction, timearray, timesignal

print('Writing to output file')
fileid = PF.outputfile 
with open (fileid,'w') as f:
      fun.Header(fileid)

with open (fileid,'a') as f:
      for w in range(sizefft):
            RPS.Receiver.timeHeader(f,timearray[w],w)
print('time: ',time.time()-t)
#