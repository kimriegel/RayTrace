import Parameterfile as Pf
import numpy as np


Fs = 24000.0     #sampling frequency
xinitial = 145.0   #boom starts
yinitial = 35.0
boomspacing = 0.6   # .6   # in between all testing points
h = 10.0           # step size in meters

outputfile = "PythonTestSimple" + str(boomspacing) + "_withdecl.txt"

radiosity = 1
complexabsorption = 0

INPUTFILE = "input/inputNASABOOM1.txt"
RecInput = "Env/Receivers/PointReceivers.txt"
ipname = 'Env/SimpleEMBuilding/SingleBuilding.obj' #geometry building 3D #blender ENVIRONMENT (*ENVIRONMENT/geometry file)


tempalphabuilding = np.zeros([Pf.absorbplanes, 8])

if complexabsorption == 1:
    tempalphabuilding[0] = [0.55, 0.55, 0.25, 0.18, 0.12, 0.07, 0.04, 0.04]
else:
    tempalphabuilding = np.zeros([Pf.absorbplanes, 8])

tempalphaground = [0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.03, 0.03]
tempalphabuilding[0, :] = [0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.03, 0.03]

percentdiffuse = 0.0


zinitial = 0.0
radius = 0.15
soundspeed = 348.537
ps = 1.0
Temp = 302.182778
hr = 20.0
theta = 1.6863372
phi = 3.44458181

def create_parameters(Fs, xinitial, yinitial, boomspacing, h):
    zinitial = 0.0
    radius = 0.15
    soundspeed = 348.537
    ps = 1.0
    Temp = 302.182778
    hr = 20.0
    theta = 1.6863372
    phi = 3.44458181

    outputfile = "PythonTestSimple" + str(boomspacing) + "_withdecl.txt"

    radiosity = 1
    complexabsorption = 0

    INPUTFILE = "input/inputNASABOOM1.txt"
    RecInput = "Env/Receivers/PointReceivers.txt"
    ipname = "Env/SimpleEMBuilding/SingleBuilding.obj"

    tempalphabuilding = np.zeros([Pf.absorbplanes, 8])

    if complexabsorption == 1:
        tempalphabuilding[0] = [0.55, 0.55, 0.25, 0.18, 0.12, 0.07, 0.04, 0.04]
    else:
        tempalphabuilding = np.zeros([Pf.absorbplanes, 8])

    tempalphaground = [0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.03, 0.03]
    tempalphabuilding[0, :] = [0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.03, 0.03]

    percentdiffuse = 0.0

    return {
        "Fs": Fs,
        "xinitial": xinitial,
        "yinitial": yinitial,
        "zinitial": zinitial,
        "boomspacing": boomspacing,
        "h": h,
        "outputfile": outputfile,
        "INPUTFILE": INPUTFILE,
        "RecInput": RecInput,
        "ipname": ipname
    }



# Fs = 24000.0     #sampling frequency
# xinitial = 145.0   #boom starts
# yinitial = 35.0    
# zinitial = 0.0
# radius = .15       #reciever 
# soundspeed = 348.537   
# ps = 1.0    #atmospheric pressure 
# Temp = 302.182778

# hr = 20.0   # relative humidity %%%
# theta = 1.6863372   #angle that the boom comes in at 
# phi = 3.44458181  #NESW 
# boomspacing = 0.6   # .6   # in between all testing points 
# # boomspacing= 0.035
# # boomspacing= 0.1
# # boomspacing= 1

# h = 10.0           # step size in meters

# outputfile = "PythonTestSimple" + str(boomspacing) + "_withdecl.txt"       # debugging? pressure signatures


# # Turn Radiosity on or off.  This will include diffuse reflections
# radiosity = 1   #different kind of reflection 
# # Turn on complex absorption
# complexabsorption = 0


# INPUTFILE = "input/inputNASABOOM1.txt" #pressure signature of boom
# RecInput = "Env/Receivers/PointReceivers.txt" #location of recievers (MICROPHONES) xyz coordinates (*RECIEVER FILE)
# ipname = 'Env/SimpleEMBuilding/SingleBuilding.obj' #geometry building 3D #blender ENVIRONMENT (*ENVIRONMENT/geometry file)

# tempalphabuilding = np.zeros([Pf.absorbplanes, 8])
# if complexabsorption == 1:
#     tempalphabuilding[0] = [0.55, 0.55, 0.25, 0.18, 0.12, 0.07, 0.04, 0.04]
# else:
#     tempalphabuilding = np.zeros([Pf.absorbplanes, 8])
# # Enter an array for absorption of alpha ground octave bands between
# # 63 and 8000
# tempalphaground = [0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.03, 0.03]
# # Enter an array for absorption of Alpha Building octave bands between
# # 63 and 8000
# tempalphabuilding[0, :] = [0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.03, 0.03]
# # what percentage of the energy is reflected diffusely between 0,1
# percentdiffuse = 0.0
