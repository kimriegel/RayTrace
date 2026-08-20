import Parameterfile as Pf
import numpy as np

#set to be replaced 
Fs = None
xinitial = None
yinitial = None
zinitial = None
radius = None
soundspeed = None
ps = None
Temp = None
hr = None
theta = None
phi = None
boomspacing = None
h = None

outputfile = None


radiosity = 1
complexabsorption = 0

INPUTFILE = "input/inputNASABOOM1.txt" #pressure signature of boom
RecInput = "Env/Receivers/PointReceivers.txt" #location of recievers (MICROPHONES) xyz coordinates (*RECIEVER FILE)
ipname = 'Env/SimpleEMBuilding/SingleBuilding.obj' #geometry building 3D #blender ENVIRONMENT (*ENVIRONMENT/geometry file)


tempalphabuilding = np.zeros([Pf.absorbplanes, 8])

if complexabsorption == 1:
    tempalphabuilding[0] = [0.55, 0.55, 0.25, 0.18, 0.12, 0.07, 0.04, 0.04]
else:
    tempalphabuilding = np.zeros([Pf.absorbplanes, 8])
# Enter an array for absorption of alpha ground octave bands between
# 63 and 8000
tempalphaground = [0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.03, 0.03]
# Enter an array for absorption of Alpha Building octave bands between
# 63 and 8000
tempalphabuilding[0, :] = [0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.03, 0.03]
# what percentage of the energy is reflected diffusely between 0,1
percentdiffuse = 0.0

# keep everything above and then create this function 
# this is to take new values and replace the parameter values (the none) with them 
def set_parameters(
    new_Fs,
    new_xinitial,
    new_yinitial,
    new_zinitial,
    new_radius,
    new_soundspeed,
    new_ps,
    new_Temp,
    new_hr,
    new_theta,
    new_phi,
    new_boomspacing,
    new_h
):

    global Fs
    global xinitial
    global yinitial
    global zinitial
    global radius
    global soundspeed
    global ps
    global Temp
    global hr
    global theta
    global phi
    global boomspacing
    global h
    global outputfile

    Fs = new_Fs
    xinitial = new_xinitial
    yinitial = new_yinitial
    zinitial = new_zinitial
    radius = new_radius
    soundspeed = new_soundspeed
    ps = new_ps
    Temp = new_Temp
    hr = new_hr
    theta = new_theta
    phi = new_phi
    boomspacing = new_boomspacing
    h = new_h

    outputfile = ("PythonTestSimple" + str(boomspacing) + "_withdecl.txt") # debugging? pressure signatures

    print("New Parameters updated")
    print("Fs:", Fs)
    print("xinitial:", xinitial)
    print("boomspacing:", boomspacing)
    print("h:", h)



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