"""
make 
list of temperatures    (c)
list of heights         (z)

one function that does math slightly differently if fz is <0 or >0
"""

import numpy as np
#plane data
#Normal = np.array((0,0,1))
#v0 = np.array((0,0,1))
#v1 = np.array((1,0,1))
#v2 = np.array((0,1,1))
#v1 = v0
#v2 = v1

#N = (v1-v0).cross(v2-v0)
#N = np.cross((v1-v0),(v2-v0))
#print(N)

#reimplement paralel check

class Strata:
    zmax = 10
    zmin = 0
    dz = -1
    z = [i for i in range(zmax,zmin,dz)]
    # find new technique to handle changing dz
    N = np.array((0,0,-1))  #normal used for calculating hit in boundaries   

    @classmethod
    def boundaryPass(cls,ray):
        if ray.direction[2] <= 0:    #if going down
            cls.N = np.array((0,0,1))   
        else:
            cls.N = np.array((0,0,-1))   #check if faster to make a generic one or to do whole thing every time
        # add hit math here in same function

#ray data
rayposition = np.array((0,0,2))
raydirection = np.array((0,0,-1))
orig = rayposition

# plane data
    #check if faster to make a generic one with replacements or to do whole thing every time

z0 = 1 # placeholder
#print(z)
v0 = np.array((0,0,z0))

if raydirection[2] <= 0:    #if going down
    N = np.array((0,0,1))   
else:
    N = np.array((0,0,-1))   #check if faster to make a generic one or to do whole thing every time


d = N.dot(v0)
t = N.dot(orig) + d # positive if ray faces plane
print(t)

P = t * raydirection + orig
print(P)

