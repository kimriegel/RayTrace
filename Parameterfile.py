#     BigBertha
INPUTFILE="inputNASABOOM1.txt"
Fs=24000.0
xinitial=-14.0
yinitial=5.0
zinitial=0.0
radius=1.5
soundspeed=343.0
ps=1.0
Temp=293.15
time=.01
hr=20.0
theta=1.92
phi=0.79
boomspacing=.6
xmin=-1
ymin=-24.0
zmin=12.0
xmax=-1
ymax=45.0
zmax=57.0
IMAX=75
h=10.0
absorbplanes=1
#allocate(tempalphabuilding(absorbplanes,8))
#Find way to rephrase
OUTPUTFILE='FortranTest3.dat'
#Turn Radiosity on or off.  This will include diffuse reflections
radiosity=0
#Turn on complex absorption
complexabsorption=0


#if(complexabsorption == 1):
#    tempalphabuilding(2,1:8)=[0.55,0.55,0.25,0.18,0.12,0.07,0.04,0.04]
##Enter an array for absorption of alpha ground octave bands between
##63 and 8000
#tempalphaground=[0.01,0.01,0.01,0.02,0.02,0.02,0.03,0.03]
##Enter an array for absorption of Alpha Building octave bands between
##63 and 8000
#tempalphabuilding(1,1:8)=[0.01,0.01,0.01,0.02,0.02,0.02,0.03,0.03]
##what percentage of the energy is reflected diffusely between 0,1
#percentdiffuse=0.5

#Broken all down to:
complexabsorption = 1
if complexabsorption == 1:
    tempalphaground = [(0.55,0.55,0.25,0.18,0.12,0.07,0.04,0.04),
    (0.01,0.01,0.01,0.02,0.02,0.02,0.03,0.03),(0.01,0.01,0.01,0.02,0.02,0.02,0.03,0.03)]
print(tempalphaground)
