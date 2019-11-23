import numpy as np

def faceNormal(face):
    a = np.array(face[0])
    b = np.array(face[1])
    c = np.array(face[2])
    d = np.cross((b-a),(c-a))   # [D]irection
    e = d//np.sqrt(d.dot(d))
    return e

def edgeTest(triangle,P,N):
    """
    Checks if some point P is inside a triangle, uses a given Normal
    """

    edge = ((triangle[1] - triangle[0]), 
            (triangle[2] - triangle[1]), 
            (triangle[0] - triangle[2]))

    chi =  ((P - triangle[0]), 
            (P - triangle[1]), 
            (P - triangle[2]))

    sha = ( N.dot(np.cross(edge[0],chi[0])) > 0, 
            N.dot(np.cross(edge[1],chi[1])) > 0, 
            N.dot(np.cross(edge[2],chi[2])) > 0)

    return np.all(sha)

def collisionCheck(FACE,VECI,F):
    """
    find if a ray hits the face for our mesh function

    the caps lock just reinforces that the variables are only used inside this function set
    """
    HUGE = 1000000.0
    F = np.array(F)         # hotfix
    N = faceNormal(FACE)    # compute plane normal
            # Finding intersection [P]oint
    # parallel check
    NF = N.dot(F)        # rayDir in notes, plane normal dot F
    isParallel = (abs(NF) < epsilon)    # bool, vD in old code
    #print('NF ',NF)
    if isParallel:
        return HUGE
        #print('parallel','\n',FACE)
#    d = np.dot(N,FACE[2])   # is tri[0] and v0 in notes
    w = VECI-FACE[2]
    si= -N.dot(w)/NF
    # find distance between origin and intersect
#    t = -(np.dot(N,VECI) + d) / NF          # dx, distance that ray travels

    if (si < 0):         # ray starts behind the face, break
        #print('ray behind face','\n',FACE)
        return HUGE,N
    elif (si > stepSize):
        #print('too far away, ignoring','\n',FACE)
        return HUGE,N    # does not hit inside step, ignore it
    else:               # if and only if it hits within the step then
        p = VECI + (si * F)
        a = np.cross(FACE[1]-FACE[0],p-FACE[0])
        b = np.cross(FACE[2]-FACE[1],p-FACE[1])
        c = np.cross(FACE[0]-FACE[2],p-FACE[2])
        if (a.dot(N) < 0):
            return HUGE,N
        elif (b.dot(N) < 0):
            return HUGE,N
        elif (c.dot(N) < 0):
            return HUGE, N
        else:
            # ideally have the hit calculation run here, but come back later
            # hit calculation assumes that distances between verts are very very small
            #for Rec in FACE:
            #    Rec.on_Hit
            return si, N

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

    def on_Hit(self, amplitude, phase):
        """ 
        My version of old receiver hit function. 
        Modifies direction and magnitude of rays with respect to each receiver
        """
        XJ = complex(0,1)
        # print('initiating hit function')

        temp1 = abs(self.magnitude) * np.exp(XJ*self.direction)
        temp2 = abs(amplitude[:])   * np.exp(XJ*phase[:])
        # print(temp2.shape)
        # print(list(temp2[-20:]))
        # print(list(temp2[:]))
        # print(list(phase))
        temp3 = temp1 + temp2 

        self.magnitude =  abs(temp3)                                 
        self.direction =  np.arctan2(np.imag(temp3) , np.real(temp3))
        # See bug log 3/13 for what happened with positions checks

ipfile = "PointReceivers.txt"

with open(ipfile) as vertex:        
    #Read in 2d-array of coordinates from file
    rho = np.genfromtxt(vertex)

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