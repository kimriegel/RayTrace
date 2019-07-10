# RayTrace
# version 1.1.0

# Kimberly A. Riegel, PHD created this program to propagate sonic booms 
# around Large structures, and to graduate. It is a ray tracing model 
# that  will include specular and diffuse reflections. It will print
# out the sound field at ear height, at relevant microphone locations,
# and at the building walls. It will read in the fft of a sonic boom 
# signature.

# Dr. Riegel, William Costa, and George Seaton porting program from Fortran to python

# Initialize variables and functions
import numpy as np

import Parameterfile as Pf
import BuildingGeometry as Bg
import Functions as Fun
import ReceiverPointSource as Rps
# import GeometryParser as Bg

import time

t = time.time()
      
# What it does not do
"""
      Interacts with geometry parser
      Have a way of reading in complex geometries - Yes, but not yet integrated
      Anything resembling radiosity
"""

def initial_signal(signal_length, fft_output):
    """
    Making the array for the initial signals.
    Input sizeFFT and output_signal
    """
    signal_length2 = int(signal_length // 2)  # Making sizeFFTTwo and setting it as an int again just to be sure
    output_frequency = np.zeros((signal_length2, 3))    # Making output array equivalent to input_array in old code
    throw_array = np.arange(1, signal_length2 + 1)     # Helps get rid of for-loops in old version

    output_frequency[:, 0] = throw_array * Pf.Fs / signal_length  # Tried simplifying the math a bit from original
    output_frequency[:, 1] = abs(fft_output[1:1 + signal_length2] / signal_length)  # Only go up to sizeFFTTwo
    output_frequency[:, 2] = np.arctan2(np.imag(fft_output[1:1 + signal_length2] / signal_length),
                                        np.real(fft_output[1:1 + signal_length2] / signal_length))

    return output_frequency

def update_freq(dx_update, alpha_update, diffusion_update):
    """
    Update ray phase and amplitude
    """
    global phase, amplitude        # works directly

    two_pi_dx_update = twopi * dx_update
    ein = phase - (two_pi_dx_update / lamb)
    zwei = ein % twopi
    masque = zwei > PI
    drei = masque * zwei - twopi
 
    phase = np.where(masque, drei, ein)
    amplitude *= ((1.0 - alpha_update) * (1.0 - diffusion_update) * np.exp(-airAbsorb * dx_update))


def vex(y, z):
    """The x coordinate of the ray 
    Used for veci"""
    return (D - FInitial[1] * y - FInitial[2] * z) / FInitial[0]


# port and import receiver file
receiverHit = 0
groundHit = 0

# Initialize counters 
PI = np.pi
twopi = PI*2
XJ = complex(0.0, 1.0)
radius2 = Pf.radius**2
raySum = 0

# Initialize receiver variables
lastReceiver = np.zeros(3)
lastReceiver2 = np.zeros(3)
receiverPoint = np.zeros(3)
receiverPoint2 = np.zeros(3)
# OC = np.empty(3)

# Read in input file
inputSignal = np.loadtxt(Pf.INPUTFILE)
K = len(inputSignal)
# masque = inputSignal > 0
HUGE = 1000000.0

# Allocate the correct size to the signal and fft arrays
sizeFFT = K
sizeFFTTwo = sizeFFT // 2
outputSignal = np.fft.rfft(inputSignal, sizeFFT)

# Create initial signal
frecuencias = initial_signal(sizeFFT, outputSignal)      # Equivalent to inputArray in original
airAbsorb = Fun.absorption(Pf.ps, frecuencias[:, 0], Pf.hr, Pf.Temp)   # sizeFFTTwo
lamb = Pf.soundspeed/frecuencias[:, 0]     # Used for updating frequencies in update function
timeArray = np.arange(K) / Pf.Fs

#       Set initial values
vInitial = np.array([Pf.xinitial, Pf.yinitial, Pf.zinitial])
xiInitial = np.cos(Pf.phi) * np.sin(Pf.theta)
nInitial = np.sin(Pf.phi) * np.sin(Pf.theta)
zetaInitial = np.cos(Pf.theta)
length = np.sqrt(xiInitial * xiInitial + nInitial * nInitial + zetaInitial * zetaInitial)
FInitial = np.array([xiInitial, nInitial, zetaInitial])
D = np.dot(FInitial, vInitial)   # equivalent to tmp

#       Create initial boom array
#  Roll this all into a function later
ySpace = Pf.boomspacing * abs(np.cos(Pf.phi))
zSpace = Pf.boomspacing * abs(np.sin(Pf.theta))
if Pf.xmin == Pf.xmax:
    rayMax = int((Pf.ymax - Pf.ymin) / ySpace) * int((Pf.zmax - Pf.zmin) / zSpace)
    print(rayMax, ' is the rayMax')

j = np.arange(1, 1 + int((Pf.ymax-Pf.ymin) // ySpace))
k = np.arange(1, 1 + int((Pf.zmax-Pf.zmin) // zSpace))
rayY = Pf.ymin + j * ySpace
rayZ = Pf.zmin + k * zSpace

boomCarpet = ((vex(y, z), y, z) for z in rayZ for y in rayY)
# Create a receiver array, include a receiver file.
alphaNothing = np.zeros(sizeFFTTwo)

# Making specific receiver points using receiver module
Rps.Receiver.initialize(Pf.RecInput)
ears = Rps.Receiver.rList           # easier to write
for R in ears:          # hotfix
    R.magnitude = np.zeros(sizeFFTTwo)
    R.direction = np.zeros(sizeFFTTwo)

#       Initialize normalization factor 
normalization = (PI*radius2)/(Pf.boomspacing**2)

outputArray1 = np.zeros((sizeFFTTwo, 6))
dHOutputArray1 = np.zeros((sizeFFTTwo, 6))

#       Define ground plane
groundHeight = 0.000000000
GroundABC = np.array([0.000000000, 0.000000000, 1.00000000])
GroundD = -groundHeight
nGround = np.array([0.0, 0.0, 1.0])

#     Allocate absorption coefficients for each surface for each frequency
alphaGround = np.zeros(sizeFFTTwo)
for D in range(0, sizeFFTTwo):       # This loop has a minimal impact on performance
    if frecuencias[D, 0] >= 0.0 or frecuencias[D, 0] < 88.0:
        alphaGround[D] = Pf.tempalphaground[0]
    elif frecuencias[D, 0] >= 88.0 or frecuencias[D, 0] < 177.0:
        alphaGround[D] = Pf.tempalphaground[1]
    elif frecuencias[D, 0] >= 177.0 or frecuencias[D, 0] < 355.0:
        alphaGround[D] = Pf.tempalphaground[2]
    elif frecuencias[D, 0] >= 355.0 or frecuencias[D, 0] < 710.0:
        alphaGround[D] = Pf.tempalphaground[3]
    elif frecuencias[D, 0] >= 710.0 or frecuencias[D, 0] < 1420.0:
        alphaGround[D] = Pf.tempalphaground[4]
    elif frecuencias[D, 0] >= 1420.0 or frecuencias[D, 0] < 2840.0:
        alphaGround[D] = Pf.tempalphaground[5]
    elif frecuencias[D, 0] >= 2840.0 or frecuencias[D, 0] < 5680.0:
        alphaGround[D] = Pf.tempalphaground[6]
    elif frecuencias[D, 0] >= 5680.0 or frecuencias[D, 0] < frecuencias[sizeFFTTwo, 0]:
        alphaGround[D] = Pf.tempalphaground[7]

alphaBuilding = np.zeros((Pf.absorbplanes, sizeFFTTwo))
for W in range(Pf.absorbplanes):        # These also look minimal
    for D in range(sizeFFTTwo):
        if frecuencias[D, 0] >= 0.0 or frecuencias[D, 0] < 88.0:
            alphaBuilding[W, D] = Pf.tempalphabuilding[W, 0]
        elif frecuencias[D, 0] >= 88.0 or frecuencias[D, 0] < 177.0:
            alphaBuilding[W, D] = Pf.tempalphabuilding[W, 1]
        elif frecuencias[D, 0] >= 177.0 or frecuencias[D, 0] < 355.0:
            alphaBuilding[W, D] = Pf.tempalphabuilding[W, 2]
        elif frecuencias[D, 0] >= 355.0 or frecuencias[D, 0] < 710.0:
            alphaBuilding[W, D] = Pf.tempalphabuilding[W, 3]
        elif frecuencias[D, 0] >= 710.0 or frecuencias[D, 0] < 1420.0:
            alphaBuilding[W, D] = Pf.tempalphabuilding[W, 4]
        elif frecuencias[D, 0] >= 1420.0 or frecuencias[D, 0] < 2840.0:
            alphaBuilding[W, D] = Pf.tempalphabuilding[W, 5]
        elif frecuencias[D, 0] >= 2840.0 or frecuencias[D, 0] < 5680.0:
            alphaBuilding[W, D] = Pf.tempalphabuilding[W, 6]
        elif frecuencias[D, 0] >= 5680.0 or frecuencias[D, 0] < frecuencias[sizeFFTTwo, 0]:
            alphaBuilding[W, D] = Pf.tempalphabuilding[W, 7]

D = np.dot(FInitial, vInitial)   # Hotfix  We used this name right above

#        Mesh the patches for the environment.  Include patching file. 
diffusionGround = 0.0
if Pf.radiosity:  # If it exists as a non-zero number
    #    import SingleBuildingGeometry
    diffusion = Pf.radiosity
else:
    diffusion = 0.0

rayCounter = 0

# These are for debugging, Uncomment this block and comment out the for loop below
# ray = 606                     # @ Pf.boomSpacing = 1
# for i in range(606):
#      ray =      next(boomCarpet)
#      rayCounter += 1
#
# if ray:
# Begin tracing
#print('Memory (before) : ' + str(mem.memory_usage()) + 'MB')
checkDirection = [0, 0, 0]
nBox = [0, 0, 0]
veci = np.array([0, 0, 0])
print('began rays')
for ray in boomCarpet:              # Written like this for readability
    veci = ray      # initial ray position
    hitCount = 0
    doubleHit = 0

    amplitude = frecuencias[:, 1]/normalization
    phase = frecuencias[:, 2]
    if Pf.h < (2*Pf.radius):
        print('h is less than 2r')
        break
    F = np.array(FInitial)                                      # Direction
    for I in range(Pf.IMAX):      # Making small steps along the ray path.
        # For each step we should return, location, phase and amplitude
        dxReceiver = HUGE
        # Find the closest sphere and store that as the distance
        for R in ears:
            # The way that tempReceiver works now, it's only used here and only should be used here.
            # It's not defined inside the receiver because it's ray dependant.
            tempReceiver = R.SphereCheck(radius2, F, veci)    # Distance to receiver
            if receiverHit >= 1:  # if you hit a receiver last time, don't hit it again
                if np.all(R.position == lastReceiver):
                    tempReceiver = HUGE
                if np.all(F == checkDirection):
                    OC = R.position - veci
                    OCLength = np.dot(OC, OC)
                    if OCLength < radius2:
                        tempReceiver = HUGE
            if receiverHit >= 2:
                if np.all(R.position == lastReceiver):
                    tempReceiver = HUGE
            if tempReceiver < dxReceiver:
                dxReceiver = tempReceiver
                receiverPoint = R.position
            elif tempReceiver == dxReceiver and tempReceiver != HUGE:
                receiverCheck = tempReceiver

# We need to double check that double hit actually works.  R2 is not really
# a thing, we should make sure it is doing what we want.
                if np.all(R.position == receiverPoint):
                    doubleHit = 0
                else:
                    R2 = R
                    doubleHit = 1
                    print('double hit')

            #     Check Intersection with ground plane
        GroundN = GroundABC
        GroundVD = GroundN[0] * F[0] + GroundN[1] * F[1] + GroundN[2] * F[2]
        if groundHit == 1:
            dxGround = HUGE
        elif GroundVD != 0.0:
            GroundVO = ((GroundN[0] * veci[0] + GroundN[1] * veci[1] + GroundN[2] * veci[2]) + GroundD)
            dxGround1 = (-1.000) * GroundVO * 1.000 / GroundVD
            dxGround = dxGround1
            Vecip1 = veci + dxGround * F
            tmp = (GroundABC[0] * Vecip1[0] + GroundABC[1] * Vecip1[1] + GroundABC[2] * Vecip1[2] + GroundD)
            if dxGround < 0.0:
                dxGround = HUGE
        else:
            dxGround = HUGE

        #     Check intersection with building
        dxBuilding = HUGE
        hit = 0
        planeHit = 0
        #     Check intersection with Boxes
        for Q in range(0, Bg.BoxNumber):
            dxNear, dxFar, hit, planeHit = Fun.box(Bg.BoxArrayNear[Q], Bg.BoxArrayFar[Q], veci, F)
            if dxNear < dxBuilding:
                dxBuilding = dxNear
                Vecip1 = veci + np.multiply(dxBuilding, F)
                whichBox = Q
                nBox = Fun.plane(Vecip1, Bg.BoxArrayNear[whichBox], Bg.BoxArrayFar[whichBox], planeHit)
        # This part doesn't really work well.  We have not incorporated it.
        # Eventually all interactions will be triangles anyway so I'm leaving it here to be updated.

        #   Check intersection with Triangles
#        if Bg.TriangleNumber > 0:
#            for Q in range(0, Bg.TriangleNumber):
#                dxNear, behind = Fun.Polygon(veci, F, Q, 3, Bg.TriangleNumber, Bg.PointNumbers, Bg.TriangleArray,
#                                             Bg.BuildingPoints, normal, FaceNormalNo, FaceNormals)
#                if dxNear < dxBuilding:
#                    dxBuilding = dxNear
#                    nBox = normal
#                    whichBox = Q
        #     Check intersection with Squares
#        if Bg.SquareNumber > 0:
#            for Q in range(0, Bg.SquareNumber):
#                dxNear, behind = Fun.Polygon(veci, F, Q, 4, SquareNumber,
        #                PointNumbers, SquareArray, BuildingPoints,
        #                normal, FaceNormalNo, FaceNormals)
#                if dxNear < dxBuilding:
#                    dxBuilding = dxNear
#                    nBox = normal
#                    whichBox = Q
        buildingHit = 0
        receiverHit = 0
        groundHit = 0

        #     Check to see if ray hits within step size
        if dxReceiver < Pf.h or dxGround < Pf.h or dxBuilding < Pf.h:
            dx = min(dxReceiver, dxGround, dxBuilding)
            #  if the ray hits a receiver, store in an array.  If the ray hits two, create two arrays to store in.
            for R in ears:
                if dx == R.dx:
                    # print('Ray ',ray +1,' hit receiver ',R.recNumber,' at step ',I)
                    print('Ray ', rayCounter, ' hit receiver ', R.recNumber)
                    veci += (dx * F)
                    receiverHit = 1
                    checkDirection = F
                    if doubleHit == 1:
                        receiverHit = 2
                    hitCount = hitCount + 1
                    update_freq(dx, alphaNothing, 0)
                    lastReceiver = receiverPoint
                    outputArray1[:, 0] = frecuencias[:, 0]
                    outputArray1[:, 1:4] = receiverPoint[:]
                    outputArray1[:, 5] = phase[:]
                    # print(list(ears[1].magnitude))
                    if doubleHit == 1:
                        # R2 = R      #Supposed to be other R, but just a placeholder for now
                        R.on_Hit(amplitude/2, phase)
                        R2.on_Hit(amplitude/2, phase)
                    else:
                        R.on_Hit(amplitude, phase)

                    # if(doubleHit==1):
                    #      outputArray1[:,4]=amplitude[:]/2.0
                    #      dHOutputArray1[:,0]=inputArray[:,0]
                    #      dHOutputArray1[:,1:4]=receiverPoint2[:]
                    #      dHOutputArray1[:,4]=amplitude[:]/2.0
                    #      dHOutputArray1[:,5]=phase[:]
                    #      lastReceiver2 = receiverPoint2
                    # else:
                    #      outputArray1[:,4]=amplitude[:]
                    # tempArray=Fun.receiverHITFUNC(sizeFFT,outputArray1,Rps.arraySize,tempArray)
                    # looks like it does the same thing as on_Hit. Here later
                    # R.on_Hit(amplitude,phase)
                    # if (doubleHit==1):
                    #      tempArray=Fun.receiverHITFUNC(sizeFFT,dHOutputArray1,Rps.arraySize,tempArray)
                    #      Using objects may circumvent the need to have this, but it stays for now
                    #      count+=1
                    # count+=1

            if abs(dx - dxGround) < 10.0**(-13.0):  # If the ray hits the ground then bounce and continue
                veci += (dxGround * F)

                tmp = np.dot(GroundABC, veci)
                if tmp != GroundD:
                    veci[2] = 0
                print('hit ground at ', I,veci)
                dot1 = np.dot(F, nGround)
                n2 = np.dot(nGround, nGround)
                F -= (2.0 * (dot1 / n2 * nGround))
                length = np.sqrt(np.dot(F, F))
                groundHit = 1
                twoPiDx = twopi * dxGround
                #     Loop through all the frequencies
                update_freq(dxGround, alphaGround, diffusionGround)
#                if Pf.radiosity == 1 and (diffusionGround != 0.0):
#                    for Q in range(0, PatchNo):
#                        if formFactors[0, Q, 1] == 1:
#                            if (veci[0] <= (patchArray[Q, W, 0] + 0.5 * patchArray[Q, W, 3]) and
            #                            veci[0]>=(patchArray[Q, W, 0] - 0.5 * patchArray[Q, W, 3])):
#                                if veci[1] <= (patchArray[Q, W, 1] + 0.5 * patchArray[Q, W, 4]) and
            #                                veci[1]>=(patchArray[Q, W, 1] - 0.5 * patchArray[Q, W, 4]):
#                                    if veci[2] <= (patchArray[Q, W, 2] + 0.5 * patchArray[Q, W, 5]) and
            #                                    veci[2]>=(patchArray[Q, W, 2] - 0.5 * patchArray[Q, W, 5]):
#                                        temp2 = complex(abs(patchArray[Q, W, 6])*np.exp(XJ*patchArray[Q, W, 7]))
#                                        temp3 = complex(abs(amplitude[W] * (1.0 - alphaGround[W]) * diffusionGround *
            #                                        exp(-m * dxGround)) * exp(1j * phaseFinal))
#                                        temp4 = temp2 + temp3
#                                        patchArray[Q, W, 6] = abs(temp4)
#                                        patchArray[Q, W, 7] = np.arctan(temp4.imag,temp4.real)
            if dx == dxBuilding:   # if the ray hits the building then change the direction and continue
                veci += (dx * F)
                print('hit building at step ', I,veci)
                n2 = np.dot(nBox, nBox)
                nBuilding = nBox / np.sqrt(n2)
                dot1 = np.dot(F, nBuilding)
                F -= (2.0 * (dot1 / n2 * nBuilding))
                length = np.sqrt(np.dot(F, F))
                buildingHit = 1
                # We need to look into complex absorption and see if this is really the best way.
#                if Pf.complexAbsorption:
#                    if Pf.absorbPlanes == 2:
#                        if (veci[2] > 0.0) and (veci[2] < height1):
#                            alpha = alphaBuilding[0, :]
#                        elif veci[2] > height1 and veci[2] <= height2:
#                            alpha = alphaBuilding[1, :]
#                    if Pf.absorbPlanes == 3:
#                        if veci[2] > height2 and veci[2] <= height3:
#                            alpha = alphaBuilding[2, :]
#                    if Pf.absorbPlanes == 4:
#                        if veci[2] > height3:
#                            alpha = alphaBuilding[4, :]
#                else:
                alpha = alphaBuilding[0, :]
                update_freq(dx, alpha, diffusion)
        else:  # If there was no interaction with buildings then proceed with one step.
            veci += (Pf.h * F)
            update_freq(Pf.h, alphaNothing, 0)
    rayCounter += 1
    print('finished ray', rayCounter)

# Radiosity removed for readability

# Reconstruct the time signal and print to output file
for R in ears:
    R.timeReconstruct(sizeFFT)

print('Writing to output file')
fileid = Pf.outputfile
with open(fileid, 'w') as f:
    Fun.header(fileid)

with open (fileid, 'a') as f:
    for w in range(sizeFFT):
        Rps.Receiver.timeHeader(f, timeArray[w], w)
print('time: ', time.time()-t)
#


testfile = 'OutputTest.txt'
with open (testfile,'w') as test:
      rec = Rps.Receiver.rList[1]
      #print(rec.recNumber,file=test)
      for w in range(sizeFFT):
            print(rec.signal[w],file=test)
