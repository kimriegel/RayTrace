# this function detects when it crosses a plane
# much of it is hardcoded, but that can all be fixes
import numpy as np
import pywavefront as pwf
#from Parameterfile import h as stepSize     #temporary, just used to make sure we do not overstep
stepSize = 100
#FaceNormals = [(-1,0,0),(0,1,0),(1,0,0),(0,-1,0),(0,0,1)]  Desired
epsilon = 1e-6  # how small angle between ray and plane has to be to count as parallel

def faceNormal(face):
    a = np.array(face[0])
    b = np.array(face[1])
    c = np.array(face[2])
    d = np.cross((b-a),(c-a))   # [D]irection
    #normal = d/np.sqrt(d.dot(d))   # doesn't seem like we need to normalize
    # Where n is [n]ot [a] [n]umber, return 0, else return n
    #return normal
    return d

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

#def LinePlaneCollision(planeNormal, planePoint, rayDirection, rayPoint, epsilon=1e-6):
# 
#	ndotu = planeNormal.dot(rayDirection)
#	if abs(ndotu) < epsilon:
#		raise RuntimeError("no intersection or line is within plane")
# 
#	w = rayPoint - planePoint
#	si = -planeNormal.dot(w) / ndotu
#	Psi = w + si * rayDirection + planePoint
#	return Psi

#def foo(FACE,VECI,F):
def collisionCheck(FACE,VECI,F):
    """
    find if a ray hits the face for our mesh function

    the caps lock just reinforces that the variables are only used inside this function set
    """
    F = np.array(F)         # hotfix
    N = faceNormal(FACE)    # compute plane normal
            # Finding intersection [P]oint
    # parallel check
    NF = np.dot(N,F)        # rayDir in notes, plane normal dot F
    isParallel = (abs(NF) < epsilon)    # bool, vD in old code
    #print('NF ',NF)
    if isParallel:
        #print('parallel',FACE)
        return False        #ray does not hit, find an output to express that

    d = np.dot(N,FACE[0])   # is tri[0] and v0 in notes

    # find distance between origin and intersect
    t = -(np.dot(N,VECI) + d) / NF          # dx, distance that ray travels
    #print('t is ', t)
    if (t < 0):         # ray starts behind the face, break
        #print('ray behind face',FACE)
        return False
    elif (t > stepSize):
        #print('too far away, ignoring',FACE)
        return False    # does not hit inside step, ignore it

    else:               # if and only if it hits within the step then
        p = VECI + (t * F)
        #print('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
        #print(p,FACE)
        isHit = edgeTest(FACE,p,N)
        print(isHit)
        if isHit == True:
            print('hits at ', p,FACE)
        return isHit
        #return p        # should return p as it is the dx, just a placeholder for now


 
#if __name__ == "__main__":
#
#    #ipname = 'Env/SingleBuilding.obj'
#    ipname = 'SingleBuilding.obj'
#    #ipname = 'TwoWalls.obj'
#    ipfile = pwf.Wavefront(ipname)    # Read in geometry file
#    env = pwf.ObjParser(ipfile,ipname, strict=False, encoding="utf-8", 
#            create_materials=True, collect_faces=True, parse=True, cache=False)
#    vertices = env.wavefront.vertices
#    faces = env.mesh.faces
#    #normals = env.normals               # Delete normal[5] needed
#    #print(faces)
#    #print(normals)
#    faceNormals = env.normals
#    print(len(vertices))
#    print(len(faces))
#    print(len(faceNormals))
#    print(faceNormals[0])
#    faceNormalNo = len(faceNormals)
#    Boxnumber = 1     # supposed to import from s, come back to this later
#        # Is this similar to Will's bands?
#    Boxarraynear=np.array([10,10,0])
#    Boxarrayfar= np.array([64.4322,46.9316,8.2423])
#    
#    #for v in vertices:
#    #    print(v)
#        #pass
#    #print(faces)
#    #print(vertices[4])
#    #for f in faces:
#    #    print(f)
#    #for n in normals:
#    #    print(n)

#ipname = 'Env/SingleBuilding.obj'
#ipname = 'Env/SingleBuildingTest.obj'
ipname = 'Env/SingleBuilding.obj'
#ipname = 'TwoWalls.obj'
ipfile = pwf.Wavefront(ipname)    # Read in geometry file
env = pwf.ObjParser(ipfile,ipname, strict=False, encoding="utf-8", 
        create_materials=True, collect_faces=True, parse=True, cache=False)
vertices = env.wavefront.vertices                                           # useful
faces = env.mesh.faces                                                      # list of keys to vertices


Boxnumber = 1     # supposed to import from s, come back to this later
    # Is this similar to Will's bands?
Boxarraynear=np.array([10,10,0])
Boxarrayfar= np.array([64.4322,46.9316,8.2423])

TriangleNumber=0
SquareNumber=0
PolyBuilding=0

mesh = [np.array((vertices[f[0]],vertices[f[1]],vertices[f[2]])) for f in env.mesh.faces]

veci = np.array([65.64138185, 26.74269388,  8.06987805])
F =    np.array([-0.94808515, -0.29638514,  0.11528397])


#for face in mesh:
##    print(face)
#    foo = collisionCheck(face,veci,F)
    #if foo == True:
    #    print(foo,hit,'dxbuilding: ',dxBuilding)





# start here
myFaces = []
    # trying to make more usable faces
#print(mesh)
for f in env.mesh.faces:
    myFaces.append((vertices[f[0]],vertices[f[1]],vertices[f[2]]))
#print(myFaces)

print(env.mesh.faces[0])
print(env.wavefront.vertices[0])
print(env.wavefront.vertices[0][0])
AAAA = [np.array((vertices[f[0]],vertices[f[1]],vertices[f[2]])) for f in env.mesh.faces]
#print(AAAA[0])
#for i in mesh:
#    print(i)
#    print(i[0])
#    print(i[:,0])
#BruteForce = 
#face = myFaces
#print(  vertices[f[0]][0],
#        vertices[f[1]][0],
#        vertices[f[2]][0])
###print(  vertices[f[0]][0],      #YES
###        vertices[f[0]][1],
###        vertices[f[0]][2])

print(  vertices[f[0]][0],      #YES
        vertices[f[0]][2],
        vertices[f[0]][1])

print(
        (vertices[f[0]][0], vertices[f[0]][2], vertices[f[0]][1]),
        (vertices[f[1]][0], vertices[f[1]][2], vertices[f[1]][1]),
        (vertices[f[2]][0], vertices[f[2]][2], vertices[f[2]][1])
)

print(
    
        np.array((vertices[f[0]][0], vertices[f[0]][2], vertices[f[0]][1])),
        np.array((vertices[f[1]][0], vertices[f[1]][2], vertices[f[1]][1])),
        np.array((vertices[f[2]][0], vertices[f[2]][2], vertices[f[2]][1])) 

)

print(
type((        
        np.array((vertices[f[0]][0], vertices[f[0]][2], vertices[f[0]][1])),
        np.array((vertices[f[1]][0], vertices[f[1]][2], vertices[f[1]][1])),
        np.array((vertices[f[2]][0], vertices[f[2]][2], vertices[f[2]][1])) 
))
)

print(
            [
        np.array((vertices[f[0]][0], vertices[f[0]][2], vertices[f[0]][1])),
        np.array((vertices[f[1]][0], vertices[f[1]][2], vertices[f[1]][1])),
        np.array((vertices[f[2]][0], vertices[f[2]][2], vertices[f[2]][1])) 
        ]
)
force = [
        np.array((vertices[f[0]][0], vertices[f[0]][2], vertices[f[0]][1])),
        np.array((vertices[f[1]][0], vertices[f[1]][2], vertices[f[1]][1])),
        np.array((vertices[f[2]][0], vertices[f[2]][2], vertices[f[2]][1])) 
        ]
print(force[0])

#print([
#        [
#        (vertices[f[0]][0], vertices[f[0]][2], vertices[f[0]][1]),
#        (vertices[f[1]][0], vertices[f[1]][2], vertices[f[1]][1]),
#        (vertices[f[2]][0], vertices[f[2]][2], vertices[f[2]][1])]
#    for f in env.mesh.faces
#    ])

force = ([
        [
        (vertices[f[0]][0], vertices[f[0]][2], vertices[f[0]][1]),
        (vertices[f[1]][0], vertices[f[1]][2], vertices[f[1]][1]),
        (vertices[f[2]][0], vertices[f[2]][2], vertices[f[2]][1])]
    for f in env.mesh.faces
    ])

print(force[0])
print(force[0][1])
print(force)
BRUTEFORCE = 

([
        [
        ((vertices[f[0]][0], vertices[f[0]][2], vertices[f[0]][1])),
        ((vertices[f[1]][0], vertices[f[1]][2], vertices[f[1]][1])),
        ((vertices[f[2]][0], vertices[f[2]][2], vertices[f[2]][1]))]
    for f in env.mesh.faces
    ])


mesh = [np.array((vertices[f[0]],vertices[f[1]],vertices[f[2]])) for f in env.mesh.faces]
#mesh = [np.array((vertices[f[:,0]],vertices[f[:,1]],vertices[f[:,2]])) for f in env.mesh.faces]        #nope


#print(myFaces)

#####"""
#####Checks if and where a ray hits a plane
#####"""
######epsilon = 0.0
#####planePoint = np.array(mesh[0])
#####normalPlane = np.array(faceNormal(mesh[0]))
######veci = np.array((0, 0, 0))
######F = np.array((1, -1, 1))
######F = np.array((1, 1, 1))
#####veci = np.array((10, 40, 0))
#####F = np.array((0, 0, 1))
#####
#####print(planePoint)
#####collisionCheck(planePoint,veci,F)
#####
####################print(np.array(faceNormals[0]).dot((1,1,1)))    # n dot u   or vD in old
###################vD = np.array(faceNormal(mesh[0])).dot(F)
###################if vD <= epsilon:
###################    print("Nyet, no hit")
###################    #return
#################### a donde
###################ein = veci - planePoint                         # w
###################zwei = -np.dot(normalPlane, ein)/ vD            # si
###################drei = ein + zwei * F + planePoint              # psi
####################return drei
#################### el fin

#if __name__=="__main__":
#	#Define plane
#	planeNormal = np.array([0, 0, 1])
#	planePoint = np.array([0, 0, 5]) #Any point on the plane
# 
#	#Define ray
#	rayDirection = np.array([0, -1, -1])
#	rayPoint = np.array([0, 0, 10]) #Any point along the ray
# 
#	Psi = LinePlaneCollision(planeNormal, planePoint, rayDirection, rayPoint)
#	print ("intersection at", Psi)



#     Check intersection with Boxes
#for Q in range(0, Bg.BoxNumber):
#    dxNear, dxFar, hit, planeHit = Fun.box(Bg.BoxArrayNear[Q], Bg.BoxArrayFar[Q], veci, F)
#    if dxNear < dxBuilding:
#        dxBuilding = dxNear
#        Vecip1 = veci + np.multiply(dxBuilding, F)
#        whichBox = Q
#        nBox = Fun.plane(Vecip1, Bg.BoxArrayNear[whichBox], Bg.BoxArrayFar[whichBox], planeHit)
## This part doesn't really work well.  We have not incorporated it.
## Eventually all interactions will be triangles anyway so I'm leaving it here to be updated.

#   Check intersection with Triangles
#if Bg.TriangleNumber > 0:
#    for Q in range(0, Bg.TriangleNumber):
#        dxNear, behind = Fun.Polygon(veci, F, Q, 3, Bg.TriangleNumber, Bg.PointNumbers, Bg.TriangleArray,
#                                     Bg.BuildingPoints, normal, FaceNormalNo, FaceNormals)
#        if dxNear < dxBuilding:
#            dxBuilding = dxNear
#            nBox = normal
#            whichBox = Q
#     Check intersection with Squares

###def tri(veci, F, Q, Number, PointNumbers, PolyArray, v, normal,FaceNormalNo,vn,dxbuilding,behind):
###   """
###   Lab notebook 5/16
###   This is an attempt to merge box and polygon into one function since
###   we are working entirely in triangular meshes now
###   ********************************Untested***********************************
###   A 1:1 translation was made from Fortran. This is the closest match to
###   what we are trying to do with triangle geometry. However there is no
###   readily available geometry file to test this.
###   [No Description given in Fortran]
###   """
###   size = 3
###   G = np.zeros((size,2))  # (3,2)
###   # inits
###   NC=0
###   behind=0
###   normal = vn[PolyArray[Q,1],:]
###   #normal[0]=FaceNormals[int(PolyArray[Q,1]),0]
###   #normal[1]=FaceNormals[int(PolyArray[Q,1]),1]
###   #normal[2]=FaceNormals[int(PolyArray[Q,1]),2]
###   d=-np.dot(normal,ValueError(PolyArray[Q,2]))
###   Vd=np.dot(normal,F)
###   if Vd >= 0.0:
###       dxbuilding = HUGE
###   V0= -(np.dot(normal,veci)+d)
###   t=V0/Vd
###   if(t < 0.0):
###       dxbuilding=HUGE
###       behind = 1
###       #Stage 1
###   intersection = veci + F*t
###   maximum = max(abs(normal))
###       # G: What if two normal values are the same? Anyway:
###   if(maximum == abs(normal[0])):
###       for P in range(size):
###           G[P,:] = (intersection[1]-v[int(PolyArray[Q,1+P]),1]
###                     ,intersection[2]-v[int(PolyArray[Q,1+P]),2])
###   elif (maximum == normal[1]):
###       for P in range(size):
###           G[P,:] = (intersection[0]-v[int(PolyArray[Q,1+P]),0]
###                     ,intersection[2]-v[int(PolyArray[Q,1+P]),2])
###   elif (maximum == normal[2]):
###       for P in range(size):
###           G[P,:] = (intersection[0]-v[int(PolyArray[Q,1+P]),0]
###                     ,intersection[1]-v[int(PolyArray[Q,1+P]),1])
###   #Stage 2
###   for P in range(size):
###       if P == size:
###           if G[P,1] < 0.0:
###               SH = -1
###           else:
###               SH = 1
###           if G[0,1] < 0.0:
###               NSH = -1
###           else:
###               NSH = 1
###       else:
###           if G[P,1] < 0.0:
###               SH = -1
###           else:
###               SH = 1
###           if G[P+1,2] < 0.0:
###               NSH = -1
###           else:
###               NSH = 1
###       if SH != NSH:
###           if (P == size):
###               if (G[P,0] > 0.0) and (G[0,0]>0.0):
###                   NC += 1
###               elif (G[P,0]> 0.0) or (G[0,0] > 0.0):
###                   if (G[P,0]-(G[P,1]*(G[P+1,0]-G[P,0])/(G[P+1,1]-G[P,1]))) > 0.0:
###                       NC += 1
###           else:
###               if (G[P,0] > 0.0) and (G[P+1,0] > 0.0):
###                   NC += 1
###               elif (G[P,0] > 0.0) or (G[P+1,1] > 0.0):
###                   if (G[P,0]-(G[P,1]*(G[P+1,0]-G[P,0])/(G[P+1,1]-G[P,1]))) > 0.0:
###                       NC += 1
###       odd = NC % 2    #get remainder to find if odd or not
###       # This was this way in original fortran
###       if odd:
###           dxbuilding = t
###       else:
###           dxbuilding = HUGE
###       return dxbuilding,behind