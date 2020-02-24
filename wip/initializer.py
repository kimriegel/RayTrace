# wip

# This file initializes all of the important variables.  
# This file is called from raycaster.py 


# include built in math library

import numpy as np

# initialize all user defined variables.

import Parameterfile as PF

#include all functions needed for ray trace
import Functions as fun

# initialie all receivers. 
import ReceiverPointSource as RPS

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

def boomplane(radius,A,B,C,D,theta,phi,xmin,ymin,zmin,xmax,ymax,zmax,arraysize):
    """    This function creates an equally spaced grid of size "step" apart    """
    receiverarray=np.zeros((arraysize,3))
    yspace=radius*abs(np.cos(phi))
    zspace=radius*abs(np.sin(theta))
    count = 0
    i=0
    j=0
    if xmin == xmax:
        for i in range(int((zmax-zmin)//zspace)):        
            for j in range(int((ymax-ymin)//yspace)):
                receiverarray[count,0]=(D-B*(ymin+(j+1)*yspace)-C*(zmin+(i+1)*zspace))/A
                receiverarray[count,1]=ymin+(j+1)*yspace
                receiverarray[count,2]=zmin+(i+1)*zspace
                count=count+1
        sizex=int((ymax-ymin)//(yspace))
        sizey=int((zmax-zmin)//zspace)
        sizez=1
    if ymin == ymax:
        for i in range(int(xmax-xmin)//(xspace)):
            for j in range(int((zmax-zmin)//(zspace))):
                receiverarray[count,0]=xmin+(i+1)*xspace
                receiverarray[count,1]=(D-A*(xmin+(i+1)*xspace)-C*(zmin+(j+1)*zspace))/B
                receiverarray[count,2]=zmin+j*zspace
                count=count+1
        sizex=int((zmax-zmin)/zspace)
        sizey=int((xmax-xmin)/(xspace))
        sizez=1
    if zmin == zmax:
        for i in range(int((xmax-xmin)/(xspace))):
            for j in range(int((ymax-ymin)/(yspace))):
                receiverarray[count,0]=xmin+(i+1)*xspace  
                receiverarray[count,1]=ymin+(j+1)*yspace
                receiverarray[count,2]=(D-A*(xmin+(i+1)*xspace)-B*(ymin+(j+1)*yspace))/C
                count=count+1
        sizex=int((xmax-xmin)/(xspace))
        sizey=int((ymax-ymin)/yspace)
        sizez=1
    return receiverarray, sizex, sizey, sizez

def absorption(ps,freq,hr,Temp):
    """
    This function computes the air absorption for a given frequency, 
    ambient pressure, relative humidity and temperature.
    """
# Define all variables and reference values
    ps0=1.0
    hr=20.0
    T0=293.15
    T01=273.16
    F=freq/ps

# Compute all relevant parameters
    psat=ps0*10**(-6.8346*(T01/Temp)**1.261+4.6151)
    h=ps0*(hr/ps)*(psat/ps0)
    FrN=1/ps0*(T0/Temp)**(1/2)*(9+280*h*np.exp(-4.17*((T0/Temp)**(1/3)-1)))
    FrO=1/ps0*(24+4.04*10**4*h*((.02+h)/(.391+h)))
    term1=0.01275*(np.exp((-2239.1/Temp))/(FrO+F**2/FrO))
    term2=0.1068*(np.exp(-3352/Temp)/(FrN+F**2/FrN))
    ABSORPTION=ps0*F**2*((1.84*10**(-11.0)*(Temp/T0)**(0.5)*ps0)+(Temp/T0)**(-5.0/2.0)*(term1+term2))

    return ABSORPTION

# port and import receiver file
receiverhit=0
groundhit=0

# Initialize counters 
PI = np.pi
twopi = PI*2
XJ=(0.0,1.0)
radius2 = PF.radius**2
raysum=0

# Initiailize receiver variables
lastreceiver = np.empty(3)
lastreceiver2 = np.empty(3)
OC = np.empty(3)

# Read in input file
with open(PF.INPUTFILE) as IPFile:
      inputsignal=np.loadtxt(IPFile)
K=len(inputsignal)
HUGE=1000000.0

# Allocate the correct size to the signal and fft arrays
sizefft=K
sizeffttwo=sizefft//2
outputsignal=np.fft.fft(inputsignal,sizefft)
#ampinitial=np.empty(sizeffttwo)
#phaseinitial=np.empty(sizeffttwo)

#       Create initial signal 
frecuencias = initial_signal(sizefft,outputsignal)      # Equivalent to inputarray in original
airabsorb=absorption(PF.ps,frecuencias[:,0],PF.hr,PF.Temp)        #sizeffttwo
lamb = PF.soundspeed/frecuencias[:,0]     # Used for updating frequencies in update function
timearray = np.arange(K) /PF.Fs

#       Set initial values
Vinitial =  np.array([PF.xinitial,PF.yinitial,PF.zinitial])
xiinitial  =np.cos(PF.phi)*np.sin(PF.theta)
ninitial   =np.sin(PF.phi)*np.sin(PF.theta)
zetainitial=np.cos(PF.theta)
length =    np.sqrt(xiinitial*xiinitial+ninitial*ninitial+zetainitial*zetainitial)
Finitial=np.array([xiinitial,ninitial,zetainitial])
tmp=(Finitial[0]*Vinitial[0]+Finitial[1]*Vinitial[1]+Finitial[2]*Vinitial[2])
PLANEABC=np.array([Finitial[0],Finitial[1],Finitial[2],tmp])

#       Create initial boom array
yspace=PF.boomspacing*abs(np.cos(PF.phi))
zspace=PF.boomspacing*abs(np.sin(PF.theta))
if (PF.xmin == PF.xmax):
   RAYMAX=int((PF.ymax-PF.ymin)/yspace)*int((PF.zmax-PF.zmin)/zspace)
elif(PF.ymin == PF.ymax):
   RAYMAX=int((PF.xmax-PF.xmin)/xspace)*int((PF.zmax-PF.zmin)/zspace)
elif(PF.zmin == PF.zmax):
   RAYMAX=int((PF.ymax-PF.ymin)/yspace)*int((PF.xmax-PF.xmin)/xspace)
boomarray = np.zeros((RAYMAX,2))
print(RAYMAX , ' is the RAYMAX')
boomarray,sizex,sizey,sizez=boomplane(PF.boomspacing,PLANEABC[0],PLANEABC[1],PLANEABC[2],PLANEABC[3],PF.theta,PF.phi,PF.xmin,PF.ymin,PF.zmin,PF.xmax,PF.ymax,PF.zmax,RAYMAX)

#     Create a receiver array, include a receiver file. 
alphanothing = np.zeros(sizeffttwo)

# Making specific receiver points using receiver module
RPS.Receiver.initialize(PF.RecInput)
ears = RPS.Receiver.rList           #easier to write
for R in ears:          #hotfix
      R.magnitude = np.zeros(sizeffttwo)
      R.direction = np.zeros(sizeffttwo)
RPS.arraysize = RPS.Receiver.arraysize
receiverpoint  = np.zeros(3)
receiverpoint2 = np.zeros(3)

#       Initialize normalization factor 
normalization=(PI*radius2)/(PF.boomspacing**2) 
temparray=np.empty((    RPS.Receiver.arraysize,sizeffttwo,6))
timetemparray=np.zeros((RPS.Receiver.arraysize,sizefft,5))

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
for D in range(1,sizeffttwo):       #This loop has a minimal impact on performance
    if   frecuencias[D,1] >= 0.0 or    frecuencias[D,1] < 88.0 :
        alphaground[D]=PF.tempalphaground[1]
    elif frecuencias[D,1] >= 88.0 or   frecuencias[D,1] < 177.0 :
        alphaground[D]=PF.tempalphaground[2]
    elif frecuencias[D,1] >= 177.0 or  frecuencias[D,1] < 355.0 :
        alphaground[D]=PF.tempalphaground[3]
    elif frecuencias[D,1] >= 355.0 or  frecuencias[D,1] < 710.0 :
        alphaground[D]=PF.tempalphaground[4]
    elif frecuencias[D,1] >= 710.0 or  frecuencias[D,1] < 1420.0 :
        alphaground[D]=PF.tempalphaground[5]
    elif frecuencias[D,1] >= 1420.0 or frecuencias[D,1] < 2840.0 :
        alphaground[D]=PF.tempalphaground[6]
    elif frecuencias[D,1] >= 2840.0 or frecuencias[D,1] < 5680.0 :
        alphaground[D]=PF.tempalphaground[7]
    elif frecuencias[D,1] >= 5680.0 or frecuencias[D,1] < frecuencias[sizeffttwo,1]:
        alphaground[D]=PF.tempalphaground[8]

alphabuilding = np.zeros((PF.absorbplanes,sizeffttwo))
for W in range(1,PF.absorbplanes):        #These also look minimal
    for D in range(1,PF.absorbplanes):
        if   frecuencias[D,1] >= 0.0   or  frecuencias[D,1] < 88.0:
            alphabuilding[W,D]=PF.tempalphabuilding[W,1]
        elif frecuencias[D,1] >= 88.0  or  frecuencias[D,1] < 177.0:
            alphabuilding[W,D] = PF.tempalphabuilding [W,2]
        elif frecuencias[D,1] >= 177.0 or  frecuencias[D,1] < 355.0 :
            alphabuilding[W,D] = tempalphabuilding[W,3]
        elif frecuencias[D,1] >= 355.0 or  frecuencias[D,1] < 710.0 :
            alphabuilding[W,D] = tempalphabuilding[W,4]
        elif frecuencias[D,1] >= 710.0 or  frecuencias[D,1] < 1420.0 :
            alphabuilding[W,D] = tempalphabuilding[W,5]
        elif frecuencias[D,1] >= 1420.0 or frecuencias[D,1] < 2840.0 :
            alphabuilding[W,D] = tempalphabuilding[W,6]
        elif frecuencias[D,1] >= 2840.0 or frecuencias[D,1] < 5680.0 :
            alphabuilding[W,D] = tempalphabuilding[W,7]
        elif frecuencias[D,1] >= 5680.0 or frecuencias[D,1] < frecuencias[sizeffttwo,1] :
            alphabuilding[W,D] = tempalphabuilding[W,8]

#        Mesh the patches for the environment.  Include patching file. 
diffusionground = 0.0
if PF.radiosity:  # If it exists as a non-zero number
      import SingleBuildingGeometry
      diffusion = PF.radiosity
else:
      diffusion = 0.0
