# this function detects when it crosses a plane
# much of it is hardcoded, but that can all be fixes
import numpy as np
import pywavefront as pwf
#FaceNormals = [(-1,0,0),(0,1,0),(1,0,0),(0,-1,0),(0,0,1)]  Desired

def isHit(near,far,veci,F):
    """
    Returns a bool if the surface is hit or not
    Standin for the old box function for now
    Goal is to look for dxnear 
        aka dxbuilding
    
    This is very important to register that something is hit during a 
    step
    """

def LinePlaneCollision(planeNormal, planePoint, rayDirection, rayPoint, epsilon=1e-6):
 
	ndotu = planeNormal.dot(rayDirection)
	if abs(ndotu) < epsilon:
		raise RuntimeError("no intersection or line is within plane")
 
	w = rayPoint - planePoint
	si = -planeNormal.dot(w) / ndotu
	Psi = w + si * rayDirection + planePoint
	return Psi
 
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
ipname = 'SingleBuilding.obj'
#ipname = 'TwoWalls.obj'
ipfile = pwf.Wavefront(ipname)    # Read in geometry file
env = pwf.ObjParser(ipfile,ipname, strict=False, encoding="utf-8", 
        create_materials=True, collect_faces=True, parse=True, cache=False)
vertices = env.wavefront.vertices                                           # useful
faces = env.mesh.faces                                                      # list of keys to vertices
#normals = env.normals               # Delete normal[5] needed
#print(faces)
#print(normals)
faceNormals = env.normals                                                   # ???
##print(len(vertices))
##print(len(faces))
##print(faces)
##print(len(faceNormals))
##print(vertices[0])
##print(faces[0])
##print(faceNormals[0])
faceNormalNo = len(faceNormals)
Boxnumber = 1     # supposed to import from s, come back to this later
    # Is this similar to Will's bands?
Boxarraynear=np.array([10,10,0])
Boxarrayfar= np.array([64.4322,46.9316,8.2423])

TriangleNumber=0
SquareNumber=0
PolyBuilding=0

# start here
myFaces = []
    # trying to make more usable faces
for f in env.mesh.faces:
    myFaces.append((vertices[f[0]],vertices[f[1]],vertices[f[2]]))
    #print(vertices[f[0]],vertices[f[1]],vertices[f[2]])
    
    #print(f[0])
    #print(type(f[0]))
    #x = f[0]
    #print(vertices[f[0]])

    #print(int(f[0]))
    #print(vertices(int(f[0])))
    #pass

#print(myFaces[0])
#print(myFaces)
face = myFaces

## Testing if keys point towards the same vertex
#print(faces)
#print(faces[0:2])
#print(faces[0],faces[1])
#print(faces[0][2],faces[1][2])
#print(faces[0][2] is(faces[1][2]))
## Yes

## Testing if actual objects point towards same vertex
#print(myFaces)
#print(myFaces[0:2])
#print(myFaces[0],myFaces[1])
#print(myFaces[0][2],myFaces[1][2])
#print(myFaces[0][2] is(myFaces[1][2]))
#print(id(myFaces[0][2]),id(myFaces[1][2]))
#print(id(myFaces[0][2]) == id(myFaces[1][2]))
##print(myFaces[0][2][0],myFaces[1][2][0])
## Yes, it identifies the coordinate as the same object (existing in memory)

# Trying to make normals (help?)
    # Dir = (B-A) cross (C-A)
    # Normal = Dir/len(dir)
print(myFaces[0][2]) 
A = np.array(myFaces[0][0])
B = np.array(myFaces[0][1])
C = np.array(myFaces[0][2])
Direction = np.cross((B-A),(C-A))
print(B-A)
print(C-A)
print(Direction) 
Normal = Direction/abs(Direction)
print(np.isnan(Normal))
N = np.where(np.isnan(Normal),    0,  Normal )      # gettings rid of nans
print(Normal)
print(N)

def faceNormal(face):
    a = np.array(face[0])
    b = np.array(face[1])
    c = np.array(face[2])
    d = np.cross((b-a),(c-a))   # [D]irection
    n = d/abs(d)
    # Where n is [n]ot [a] [n]umber, return 0, else return n
    normal = np.where(np.isnan(n),0,n)  
    return tuple(normal)

for i in range(len(env.normals)):
    print('Viejo: ', env.normals[i])
    print('Nueva: ', faceNormal(myFaces[i]))

"""
Checks if and where a ray hits a plane
"""
# si
epsilon = 0.0
planePoint = np.array(face[0])
normalPlane = np.array(faceNormals[0])
F = np.array((1, -1, 1))
veci = np.array((0, 0, 0))
#print(np.array(faceNormals[0]).dot((1,1,1)))    # n dot u   or vD in old
vD = np.array(faceNormals[0]).dot(F)
if vD <= epsilon:
    print("Nyet, no hit")
    #return
# a donde
ein = veci - planePoint                         # w
zwei = -np.dot(normalPlane, ein)/ vD            # si
drei = ein + zwei * F + planePoint              # psi
#return drei
# el fin

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