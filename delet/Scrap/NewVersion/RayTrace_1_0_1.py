# RayTrace
# Python 3.7.0 64-bit

# Kimberly Lefkowitz created this program to propagate sonic booms around
# large structures, and to graduate. It is a ray tracing model that 
# will include specular and diffuse reflections. It will print out the
# sound field at ear height, at relevant microphone locations, and at 
# the building walls. It will read in the fft of a sonic boom signature.

# Dr. Riegel, William Costa, and Kory Seaton porting program from Fortran to python

##### This version only exists to check the wild theories that I think about all day.

##### As I progress through this I will keep one version with my scrap paper
##### The next version will remove that so it looks less jumbled

import time
t = time.time()
# Linking other python files with initialized variables and functions
import Parameterfile_methods as PF
#import BuildingGeometry as BG  #Unused so far. May be replaced with a mesh file anyway
import numpy as np
import BuildingGeometry as BG
import Functions_Testing as fun
import math 
print('Import: ', time.time()-t)
t = time.time()

def initial_signal(signalLength,fftOutput):
    """
    Making the array for the initial signals.
    Input sizefft and outputsignal
    """
    signalLength2 = int(signalLength//2)    # Making sizeffttwo in this function and setting it as an int again jsut to be sure
    outputFrequency = np.zeros((signalLength2,3))    #Making output array    equivalent to inputarray in old code
    throwarray = np.arange(1,signalLength2 + 1)     #Helps get rid of for-loops in old version

    outputFrequency[:,0] = throwarray * PF.Fs / signalLength #Tried simplifying the math a bit from original
    outputFrequency[:,1] = abs(fftOutput[:signalLength2]/signalLength) #Only go up to sizeffttwo 
    outputFrequency[:,2] = np.arctan2( np.imag(fftOutput[:signalLength2]/signalLength) , np.real(fftOutput[:signalLength2]/signalLength)) 

    return outputFrequency

def stretch3(terry):
    """
    Stetches a 1d array by 3 
    """
    terry2 = np.append(terry[:,None],terry[:,None],axis=1) #I didn't know how to integrate both at once
    terry = np.append(terry2,terry[:,None],axis=1)
    #print(terry)
    return terry

def stretch(terry):
    """
    Stretches a 1d array by whatever
    """
    scrap = np.ones(3)
    temp = terry[:,None] * scrap
    #print(temp)
    return temp

#Formerly RayTrace_KAR        cutting down on files -G
# port and import receiver file
receiverhit=0
groundhit=0

        # Initialize counters 
PI=3.1415965358979323846
XJ=complex(0.0,1.0)
radius2 = PF.radius ** 2
twopi = 2.0 * PI
twopih = twopi * PF.h
S=1     #Used for radiosity
K=0
raysum=0

        # Initiailize receiver variables
lastreceiver  =  np.zeros(3)
lastreceiver2 =  np.zeros(3)
OC            =  np.zeros(3)

        # Read in input file
with open(PF.INPUTFILE) as IPFile:
      inputsignal=np.loadtxt(IPFile)
K=len(inputsignal)
HUGE=1000000.0

#F=np.empty([1,3])   #Not used anymore, but helps when I forget to replace F with Fvector

sizefft = int(K)
sizeffttwo = int(K//2)

""" Take the fft of the input signal with fftw """
outputsignal=np.fft.fft(inputsignal,sizefft)

timearray=np.zeros(sizefft)
#ampinitial=np.empty(sizeffttwo)
#phaseinitial=np.empty(sizeffttwo)
 
        # Create initial signal
#t = time.time()

frecuencias = initial_signal(sizefft,outputsignal)      # Equivalent to inputarray in original

airabsorb = fun.ABSORPTION(PF.ps,frecuencias[:,0],PF.hr,PF.Temp)
timearray = np.arange(sizefft) * 1 /PF.Fs

""" From Below """
m=airabsorb[:]      #c   
lamb=PF.soundspeed/ frecuencias[:,0]   #c

#print('Created initial signal: ',time.time()-t)

boomarray = PF.boomarray
sizex = PF.sizex
sizey = PF.sizey
sizez = PF.sizez

        # Create a receiver array, include a receiver file. 

#alphanothing = np.zeros(sizeffttwo)    # Buglog 12/27/18 Again (Do we use this)
alphanothing = 0.0 # Buglog 2/24/19 (Including it in name only but I don't think we use it oterwise)

#import ReceiverPointSource as RPS 
import ReceiverPointMethods as RPS
#receiverarray = np.array(RPS.receiverarray)
receiverarray = np.array(RPS.Receiver.Array)
receiverpoint = np.array([0.,0.,0.])
sum = 0

        # These are seemingly unused
receiverarray1=None
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

#       Initialize normalization factor 
normalization=(PI*radius2)/(PF.boomspacing * PF.boomspacing) 
timetemparray=np.zeros((RPS.Receiver.arraysize,sizefft,5))
#print(frecuencias)


#for D in range(RPS.Receiver.arraysize):
#    RPS.Receiver.rList[D] = frecuencias[:,0]    
    # No longer needs to define position. Those are constant with the object
    # No longer need to define magnitude and direction. They default to 0.0
    # Only need to define first line of inputarray(frecuencias aqui)
        # Since this is the same for all I will use it as its own default.

RPS.Receiver.de_frecuencias(frecuencias[:,0])
        #print(RPS.Receiver.rList[0].initial_frequency) # It checks out

#print(frecuencias)

#print(RPS.Receiver.rList)

#        Create ground plane
groundheight=0.000000000
GROUNDABC=np.array([0.000000000,0.000000000,1.00000000])
GROUNDD=-groundheight
nground=np.array([0.0,0.0,1.0])
alphaground=np.zeros(sizeffttwo)

for D in range(1,sizeffttwo):       #This loop has a minimal impact on performance
    if (frecuencias[D,1] >= 0.0 or frecuencias[D,1] < 88.0 ):
        alphaground[D]=PF.tempalphaground[1]
    elif (frecuencias[D,1] >= 88.0 or frecuencias[D,1] < 177.0 ):
        alphaground[D]=PF.tempalphaground[2]
    elif (frecuencias[D,1] >= 177.0 or frecuencias[D,1] < 355.0 ):
        alphaground[D]=PF.tempalphaground[3]
    elif (frecuencias[D,1] >= 355.0 or frecuencias[D,1] < 710.0 ):
        alphaground[D]=PF.tempalphaground[4]
    elif (frecuencias[D,1] >= 710.0 or frecuencias[D,1] < 1420.0 ):
        alphaground[D]=PF.tempalphaground[5]
    elif (frecuencias[D,1] >= 1420.0 or frecuencias[D,1] < 2840.0 ):
        alphaground[D]=PF.tempalphaground[6]
    elif (frecuencias[D,1] >= 2840.0 or frecuencias[D,1] < 5680.0 ):
        alphaground[D]=PF.tempalphaground[7]
    elif (frecuencias[D,1] >= 5680.0 or frecuencias[D,1] < frecuencias[sizeffttwo,1]):
        alphaground[D]=PF.tempalphaground[8]

alphabuilding = np.zeros((PF.absorbplanes,sizeffttwo))

for W in range(1,PF.absorbplanes):        #These also look minimal
    for D in range(1,PF.absorbplanes):
        if (frecuencias[:,1] >= 0.0 or frecuencias[:,1] < 88.0 ):
            alphabuilding[W,D]=PF.tempalphabuilding[W,1]
        elif (frecuencias[:,1] >= 88.0 or frecuencias[:,1] < 177.0 ):
            alphabuilding[W,D] = PF.tempalphabuilding [W,2]
        elif (frecuencias[:,1] >= 177.0 or frecuencias[:,1] < 355.0 ):
            alphabuilding[W,D] = PF.tempalphabuilding[W,3]
        elif (frecuencias[:,1] >= 355.0 or frecuencias[:,1] < 710.0 ):
            alphabuilding[W,D] = PF.tempalphabuilding[W,4]
        elif (frecuencias[:,1] >= 710.0 or frecuencias[:,1] < 1420.0 ):
            alphabuilding[W,D] = PF.tempalphabuilding[W,5]
        elif (frecuencias[:,1] >= 1420.0 or frecuencias[:,1] < 2840.0 ):
            alphabuilding[W,D] = PF.tempalphabuilding[W,6]
        elif (frecuencias[:,1] >= 2840.0 or frecuencias[:,1] < 5680.0 ):
            alphabuilding[W,D] = PF.tempalphabuilding[W,7]
        elif (frecuencias[:,1] >= 5680.0 or frecuencias[:,1] < frecuencias[sizeffttwo,1]):
            alphabuilding[W,D] = PF.tempalphabuilding[W,8]

#        Mesh the patches for the environment.  Include patching file. 
if PF.radiosity == 1 :
    #import SingleBuildingGeometry
    import BuildingGeometry
    diffusion=PF.percentdiffuse
    diffusionground=0.0
else:
    diffusion=0.0
    diffusionground=0.0

count=0

print('Congrats on initializing everything: ', time.time()-t)

print('began rays')
#ray = 606
ampinitial = np.zeros(sizeffttwo)
phaseinitial = np.zeros(sizeffttwo)
######################################################################################################################3
# Rays start here.
# I'm testing outside of a loop
# Remember to attach each run to the coresponding receiver via methods. 1/8

hitcount=0
tmpsum=0.0
doublehit=0

#t = time.time()
#for W in range(0,sizeffttwo):
#    ampinitial[W]=frecuencias[W,0]/normalization
#    phaseinitial[W]=frecuencias[W,1]
#print('Old ray frequencies: ',time.time()-t)

#t = time.time()
ampinitial = frecuencias[:,0]/normalization * np.ones((PF.RAYMAX,1))
phaseinitial = frecuencias[:,1] * np.ones((PF.RAYMAX,1))
#print('New ray frequencies: 'time.time()-t)

#Vinitial=np.array([boomarray[ray-1,0],boomarray[ray-1,1],boomarray[ray-1,2]])     # Set value equal to last itteration
#Vinitial = np.array( boomarray[1:,:3])  # horizontal from second point to end, vertical has 3 columns 
#print(Vinitial)
#Vinitial = np.array([boomarray[:PF.RAYMAX,0],boomarray[:PF.RAYMAX,1],boomarray[:PF.RAYMAX,2]])  # Take directly from boomarray up to the RayMax
    # Isn't boomarray only up to raymax anyway, you may be asking.
    # Yea, probably but I feel like I'll break it again while making this. So that limit is defined twice
if (PF.h < (2*PF.radius)): 
    print('h is less than 2r')
    #break  #Ends function here (Hopefully)
xiInitial = math.cos(PF.phi)*math.sin(PF.theta)
nInitial  = math.sin(PF.phi)*math.sin(PF.theta)
zetaInitial=math.cos(PF.theta)
Fvector = np.array([xiInitial,nInitial,zetaInitial])    #=Finitial
#veci = np.array([boomarray[:PF.RAYMAX,0],boomarray[:PF.RAYMAX,1],boomarray[:PF.RAYMAX,2]]) # Just cuz 
veci = boomarray[:PF.RAYMAX,:3]
#print(boomarray[:,:3])
#print(boomarray.shape)
# Note: I think Vinitial works best with Fvector if you convert Fvector as Fvector[:,None]  -George 1/1

# Original: Vinitial=np.array([boomarray[ray-1,0],boomarray[ray-1,1],boomarray[ray-1,2]])     #Where code diverges # <- I found what this means. It's just a Gandhi bug
# Mine: Vinitial = np.array([boomarray[:RayMax,0],boomarray[:RayMax,1],boomarray[:RayMax,2]])   Take directly from boomarray up to the RayMax



#for I in range(0,PF.IMAX):
#for I in range(0,3):
dxreceiver=HUGE
# Find the closest sphere and store that as the distance
t = time.time()
for Q in range(0,RPS.Receiver.arraysize):
    tempreceiver = fun.SPHERECHECK(RPS.Receiver.Array[Q],radius2,Fvector,veci)  #finally passed this. 

    #if (receiverhit >= 1):  #if you hit a receiver last time, don't hit it again
    # These are hard to check now, considering all the values are Huge anyway.
    # Anyway, here's a prototype that might do everything
    tempreceiver = np.where(np.all(lastreceiver == RPS.Receiver.Array[Q,:],axis=0),HUGE,tempreceiver)
    #print(tempreceiver)
    #print(RPS.Receiver.Array[Q,:])

    RPS.Receiver.rList[Q].tempreceiver = tempreceiver
    #print(RPS.Receiver.rList[Q].tempreceiver)
    
    #print(np.all(lastreceiver == RPS.Receiver.Array[Q,:],axis=0))   #horizontal axis

        # This one requires us to calculate the Fvector and the checkdirection as similarly sized arrays
        # Otherwise they won't work independently of each other. So I'll come back to it

    #    if(
    #    F[0]==checkdirection[0] and
    #    F[1] == checkdirection[1] and 
    #    F[2] == checkdirection[2]):
    #        OC[0]=RPS.Receiver.Array[Q,0]-veci[0]
    #        OC[1]=RPS.Receiver.Array[Q,1]-veci[1]
    #        OC[2]=RPS.Receiver.Array[Q,2]-veci[2] 
    #        OCLength=OC[0]*OC[0]+OC[1]*OC[1]+OC[2]*OC[2]
    #        if(OCLength < radius2):
    #            tempreceiver=HUGE
    #if(receiverhit >= 2):
        #if(lastreceiver2[0]== RPS.Receiver.Array[Q,0] and lastreceiver2[1]==RPS.Receiver.Array[Q,1] and lastreceiver2[2]==RPS.Receiver.Array[Q,2]):
        #print('recHit > 2 happens')        #This does not happen, good
        #    tempreceiver=HUGE
        
        #tempreceiver = np.where(np.all(receiverpoint == RPS.Receiver.Array[Q,:],axis=0),HUGE,tempreceiver)

""" Create a temporary receiver to hold these truth values, stretch it to the needed size (3) and copy rArray values where true """
#print(tempreceiver < dxreceiver)
#trashreceiver = tempreceiver < dxreceiver
#trashreceiver2 = np.append(trashreceiver[:,None],trashreceiver[:,None],axis=1) #I didn't know how to integrate both at once
#trashreceiver = np.append(trashreceiver2,trashreceiver[:,None],axis=1)
trashreceiver = stretch(tempreceiver < dxreceiver)

receiverpoint = np.where(trashreceiver,RPS.Receiver.Array[Q,:],receiverpoint[:])
#print(trashreceiver)
#print(receiverpoint)

# Since we're doing this in relation to individual receivers I'm skipping everything with double hit -2/7


        # GroundIntersect
""" Check intersection with ground plane """ 
GROUNDN=GROUNDABC
GROUNDVD = np.dot(GROUNDN,Fvector)

# dxground has not yet been defined
#dxground = HUGE
#dxground = np.where((groundhit==1),HUGE,dxground)   # Takes effect after a change in groundhit but besides that looks like it'll cause bugs
#print(dxground)
#dxground = np.where(    (GROUNDVD != 0.0 ),   groundfunc  ,     dxground )

                  #GROUNDVO=((GROUNDN[0]*veci[0]+GROUNDN[1]*veci[1]+GROUNDN[2]*veci[2])+GROUNDD)
                  #dxground1=(-1.000)*GROUNDVO*(1.000)/GROUNDVD
                  #dxground=dxground1
                  #Vecip1=veci+dxground*np.array(F)
                  #tmp=(GROUNDABC[0]*Vecip1[0]+GROUNDABC[1]*Vecip1[1]+GROUNDABC[2]*Vecip1[2]+GROUNDD)                  
                  #if (dxground < 0.0):
                  #      dxground=HUGE

                # Yea this is going to require nesting and just revisiting in general

#"""     Check intersection with building
#        Should probably nest the rest inside this   """
#dxbuilding=HUGE
#hit=0
#planehit=0
#
#"""     Check intersection with Boxes     """

#for Q in range(BG.Boxnumber):
#      dxnear, dxfar, hit, planehit=fun.BOX(BG.Boxarraynear[Q], BG.Boxarrayfar[Q],veci,Fvector)
#      print(veci)
#      if (dxnear < dxbuilding):
#            dxbuilding=dxnear
#            Vecip1=veci+np.multiply(dxbuilding,F)
#            whichbox=Q
#            nbox=fun.PLANE(Vecip1, BG.Boxarraynear[whichbox],BG.Boxarrayfar[whichbox], planehit)

# Has problems running array math and may get outdated in Will's update
# Moving on


#for I in range(PF.RAYMAX):
#    pass
for I in range(3):
    """ Step """
    tmpsum += PF.h
    #Vecip1=veci+(PF.h)*Fvector[:,None]
    #veci=Vecip1
    #print(veci)
    #veci += (PF.h*Fvector[:,None]) #Does that work better?
    veci += (PF.h * Fvector)
    Vecip1 = veci
    # c is constant f is for changing(flux) 
    """
            c can be made in a single array, flux is more complex
    """
    #m=airabsorb[:]      #c   
    #lamb=PF.soundspeed/ frecuencias[:,0]   #c
    phasefinal=phaseinitial[:]-(twopih)/lamb[:]   #f    Not a placeholder anymore Neither is that V
    phaseinitial[:]=phasefinal % twopi          #f    Size needed:(Raymax,sizeffttwo) Doesn't work yet. Just a placeholder
            # Muda Da
    #ampfinal=ampinitial[:]*(1- alphanothing )*np.exp(-m*PF.h)  #f
    #ampinitial[:]=ampfinal[:]                
    ampinitial *= ((1 - alphanothing) * np.exp( PF.h* -m))  #f  
    #print(ampinitial)
    phaseinitial[:] = np.where((phaseinitial > PI),     (phaseinitial - twopi),     phaseinitial)
    ##if (phaseinitial[W] > PI):
    ##      phaseinitial[W]=phaseinitial[W]-twopi
    
#print(veci)
    
    

## 1)                        
##                       for W in range(0,sizeffttwo):
#                              if(PF.complexabsorption==1):
#                                    if (PF.absorbplanes==2):
#                                          if(veci[2]>0.0 and veci[2]< height1):
#                                                alpha=alphabuilding[0,W]
#                                          elif(veci[2]>height1 and veci[2]<=height2):
#                                                alpha=alphabuilding[1,W]
#                                    if(PF.absorbplanes==3):
#                                          if(veci[2]>height2 and veci[2] <=height3):
#                                                alpha=alphabuilding[2,W]
#                                    if(PF.absorbplanes==4):
#                                          if(veci[2]>height3):
#                                                alpha=alphabuilding(4,W)
#                              else:
#                                    alpha=alphabuilding[0,W]
#                              m=airabsorb[W]
#                              lamb=PF.soundspeed/inputarray[W,0]                
#                              phasefinal=phaseinitial[W]-(twopidx)/lamb
#                              ampfinal=ampinitial[W]*(1.0-alpha)*(1.0-diffusion)*np.exp(-m*dx)
#                              ampinitial[W]=ampfinal 
#                              phaseinitial[W]=phasefinal%twopi
#                              if (phaseinitial[W]>PI):
#                                    phaseinitial[W]=phaseinitial[W]-twopi
#
#
## 2)            else:
#                  # If there was no interaction with buildings then proceed with one step. 
#                  tmpsum=tmpsum+PF.h
#                  Vecip1=veci+(PF.h)*np.array(F)
#                  veci=Vecip1
#                  #twopih=twopi*PF.h
#                  # Loop through all frequencies.
#                  for W in range (0,sizeffttwo):
#                        m=airabsorb[W]
#                        lamb=PF.soundspeed/inputarray[W,0]
#                        phasefinal=phaseinitial[W]-(twopih)/lamb
#                        ampfinal=ampinitial[W]*(1-alphanothing)*np.exp(-m*PF.h)
#                        ampinitial[W]=ampfinal[W]                 
#                        phaseinitial[W]=phasefinal%twopi
#                        if (phaseinitial[W] > PI):
#                              phaseinitial[W]=phaseinitial[W]-twopi
