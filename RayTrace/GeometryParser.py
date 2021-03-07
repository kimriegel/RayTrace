# this function detects when it crosses a plane
# much of it is hardcoded, but that can all be fixes
import numpy as np
import pywavefront as pwf
from Parameterfile import ipname
import Functions as f
from Parameterfile import h as step_size     # temporary, just used to make sure we do not overstep

# ipname = 'Env/duckscaled.obj'
ipfile = pwf.Wavefront(ipname)    # Read in geometry file
env = pwf.ObjParser(ipfile, ipname, strict=False, encoding="utf-8",
                    create_materials=True, collect_faces=True, parse=True, cache=False)
vertices = env.wavefront.vertices                                           # useful
faces = env.mesh.faces       # list of keys to vertices
# Boxnumber = 1     # supposed to import from s, come back to this later
# Is this similar to Will's bands?
# Boxarraynear=np.array([10,10,0])
# Boxarrayfar= np.array([64.4322,46.9316,8.2423])

# mesh = [np.array((vertices[f[0]],vertices[f[1]],vertices[f[2]])) for f in env.mesh.faces]

mesh = [np.array((
        (vertices[f[0]][0], vertices[f[0]][1], vertices[f[0]][2]),
        (vertices[f[1]][0], vertices[f[1]][1], vertices[f[1]][2]),
        (vertices[f[2]][0], vertices[f[2]][1], vertices[f[2]][2])))
        for f in env.mesh.faces]    # Brute force technique just to get source to x,y,z

# for face in mesh:
#    #print(face)
#    foo = collisionCheck(face,veci,F)
# print(foo)
# if foo == True:
#    print(foo,hit,'dxbuilding: ',dxBuilding)
#
#
# # start here   (debugging)
# myFaces = []
#    # trying to make more usable faces
# #print(mesh)
# for f in env.mesh.faces:
#    myFaces.append((vertices[f[0]],vertices[f[1]],vertices[f[2]]))
# #print(myFaces)
#
# collisionCheck(([[0,1,5],[1,0,5],[0,0,5]]),np.array([0,0,10]),np.array([0,-1,-1]))
