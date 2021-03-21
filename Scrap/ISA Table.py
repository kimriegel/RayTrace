# ISA Table
# Making the graph just for fun
import numpy as np


#aFeet = [*range(40,-2, -1)]     # altitude in ft
#aFeet = np.arange(40,-2, -1)     # altitude in ft
#altitude = aFeet/ 3.28083
#print(altitude)

altitude = np.arange(40,-2, -1) / 3.28083   #altitude in meters
print(altitude)

#The international Standard Atmosphere parameters (temperature, pressure, density) 
# can be provided as a function of the altitude under a tabulated form

#rho = p/ (R*T)
#delta = p/p0