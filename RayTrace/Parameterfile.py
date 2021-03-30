import numpy as np

"""
Input Boom, Receiver Locations, and Env Files
"""
INPUTFILE = "input/inputNASABOOM1.txt"
RecInput = "Env/Receivers/PointReceivers.txt"
ipname = 'Env/SimpleEMBuilding/SingleBuildingTest.obj'
"""
Signal Constants
"""
Fs = 24000.0
#time = .01 #Check if unused.

"""
Boom Positioning Variables
"""
xinitial = 145.0
yinitial = 35.0
zinitial = 0.0
theta = 1.6863372
phi = 3.44458181

xmin = -1
ymin = 30.0
zmin = 0.0
xmax = -1
ymax = 100.0
zmax = 25.0
IMAX = 75
h = 10.0

"""
Boom Spacing 
Default: 0.035
"""
# boomspacing = 0.6
# boomspacing= 0.035
# boomspacing= 0.1
boomspacing= 1

"""
Receiver Variables
"""
radius = .15

"""
Atmospheric Variables
"""
soundspeed = 348.537
ps = 1.0
Temp = 302.182778
hr = 20.0

"""
Output Variables 
"""
outputfile = "PythonTestSimple" + str(boomspacing) + ".txt"
graphName = "TestGraph"

"""
Radiosity - Include diffuse reflections
Binary 1-On, 0-Off
"""
radiosity = 1

"""
Complex Absorption
Binary 1-On, 0 - Off
"""
complexabsorption = 0

absorbplanes = 1

tempalphabuilding = np.zeros([absorbplanes, 8])
if complexabsorption == 1:
    
    # tempalphabuilding[1] = [0.55,0.55,0.25,0.18,0.12,0.07,0.04,0.04]
    tempalphabuilding[:, 1] = [0.55, 0.55, 0.25, 0.18, 0.12, 0.07, 0.04, 0.04]
else:
    tempalphabuilding = np.zeros([1, 8])
# Enter an array for absorption of alpha ground octave bands between
# 63 and 8000
tempalphaground = [0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.03, 0.03]
# Enter an array for absorption of Alpha Building octave bands between
# 63 and 8000
tempalphabuilding[0, :] = [0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.03, 0.03]
# what percentage of the energy is reflected diffusely between 0,1
percentdiffuse = 0.0
