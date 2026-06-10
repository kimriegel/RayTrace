#     BigBertha
#   Same as NASABOOM1EMParameterFile
import numpy as np




time = .01   # time between samples -- no change 

#everything above but time and add h 

xmin = -1
ymin = 30.0
zmin = 0.0
xmax = -1                   #  } size of the portion of the boom being used
ymax = 100.0
zmax = 25.0
IMAX = 75        #number of iterations 

absorbplanes = 1  
# allocate(tempalphabuilding(absorbplanes,8))
# Find way to rephrase
# outputfile = 'PythonTest1.txt'
graphName = "TestGraph"                                     # No not use full file extension here
# Will's
# outputfile = "PythonTestEnv" + str(boomspacing) + ".txt"       # debugging






# Broken all down to:
# complexabsorption = 1
# if complexabsorption == 1:
#     tempalphaground=np.array([[0.55,0.55,0.25,0.18,0.12,0.07,0.04,0.04],
#     [0.01,0.01,0.01,0.02,0.02,0.02,0.03,0.03],[0.01,0.01,0.01,0.02,0.02,0.02,0.03,0.03]])
# print(tempalphaground)

# if __name__ == "__main__":      #being lazy. You can run from here now
#    import RayTrace
# pass
