import numpy as np

class Receiver:
    initial_frequency = None    # Gives this name to all receivers
    rList = []  # See append_list
    arraysize = 1     # Number of receivers. Supposed to be 5 for this test    # I use it only for printing receiver label on graphs


    def __init__(self, position):
        """
        Create and defines position of receiver
        """
        self.position = np.array(position)
        self.recNumber = int(Receiver.arraysize)         # planned for debugging but we don't seem to use it
        Receiver.rList.append(self)  # See append_list
        # Initial values
        self.pressure = 0
        self.magnitude = 0
        self.direction = 0

        Receiver.arraysize += 1

ipfile = "PointReceivers.txt"

with open(ipfile) as vertex:        
    #Read in 2d-array of coordinates from file
    rho = np.genfromtxt(vertex)
#    print(rho)

import pywavefront as pwf

ipname = 'Env/SingleBuilding.obj'
ipfile = pwf.Wavefront(ipname)    # Read in geometry file
env = pwf.ObjParser(ipfile,ipname, strict=False, encoding="utf-8", 
        create_materials=True, collect_faces=True, parse=True, cache=False)
vertices = env.wavefront.vertices                                           # useful
faces = env.mesh.faces                                                      # list of keys to vertices

#mesh = [np.array((vertices[f[0]],vertices[f[1]],vertices[f[2]])) for f in env.mesh.faces]
#mesh = [(vertices[f[0]],vertices[f[1]],vertices[f[2]]) for f in env.mesh.faces]
#print(list(mesh))
#print(mesh[0][0])
"""
Make a list of receivers from vertex list
Run receiver list thru geometry code
"""
#for v in vertices:
#    print(v)
for v in vertices:
    Receiver(v)
print(len(Receiver.rList))

#print(vertices)

#mesh = [(vertices[f[0]],vertices[f[1]],vertices[f[2]]) for f in env.mesh.faces]
#rMesh = [(vertices[f[0]],vertices[f[1]],vertices[f[2]]) for f in env.mesh.faces]
rMesh = [(Receiver.rList[f[0]],Receiver.rList[f[1]],Receiver.rList[f[2]]) for f in env.mesh.faces]
#print(rMesh)
for q in rMesh:
    print(q)


#