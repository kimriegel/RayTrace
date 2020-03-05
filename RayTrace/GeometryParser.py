# this function detects when it crosses a plane
# much of it is hardcoded, but that can all be fixes
import numpy as np
import pywavefront as pwf
#from Parameterfile import h as stepSize     #temporary, just used to make sure we do not overstep
stepSize = 10
#FaceNormals = [(-1,0,0),(0,1,0),(1,0,0),(0,-1,0),(0,0,1)]  Desired
epsilon = 1e-6  # how small angle between ray and plane has to be to count as parallel
HUGE = 1000000.0



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
        #isHit = edgeTest(FACE,p,N)
        #print(isHit)
#        if isHit == True:
#            print('hits at ', p, si)
            return si, N
        #return p        # should return p as it is the dx, just a placeholder for now

ipname = 'Env/SingleBuilding.obj'
ipfile = pwf.Wavefront(ipname)    # Read in geometry file
env = pwf.ObjParser(ipfile,ipname, strict=False, encoding="utf-8", 
        create_materials=True, collect_faces=True, parse=True, cache=False)
#vertices = env.wavefront.vertices                                           # useful but twice

uselessV = env.wavefront.vertices                                           # useful
vertices = [uselessV[v] for v in range(len(uselessV)//2)]


faces = env.mesh.faces                                                      # list of keys to vertices

#Boxnumber = 1     # supposed to import from s, come back to this later
    # Is this similar to Will's bands?
#Boxarraynear=np.array([10,10,0])
#Boxarrayfar= np.array([64.4322,46.9316,8.2423])

#mesh = [np.array((vertices[f[0]],vertices[f[1]],vertices[f[2]])) for f in env.mesh.faces]

mesh = [np.array((
        (vertices[f[0]][0], vertices[f[0]][2], vertices[f[0]][1]),
        (vertices[f[1]][0], vertices[f[1]][2], vertices[f[1]][1]),
        (vertices[f[2]][0], vertices[f[2]][2], vertices[f[2]][1])))
    for f in env.mesh.faces]    # Brute force technique just to get source to x,y,z

#for face in mesh:
#    #print(face)
#    foo = collisionCheck(face,veci,F)
    #print(foo)
    #if foo == True:
    #    print(foo,hit,'dxbuilding: ',dxBuilding)
#
#
## start here   (debugging)

#######################################################################################################

###print(vertices)
##with open("foo.py",'w') as f:
##    #print(('vertices = ',vertices),file=f)
##    #print('vertices = ')
##    #for v in vertices:
##    f.write('vertices = [\n')
##    for v in range(len(vertices)//2):
##        #f.write(str(v))
##        #f.write(str(('\t',vertices[v])))
##        print(( vertices[v]), file=f)
##        #f.write(str(('\t',vertices[v])))
##        #f.write('\n,')
##    f.write(']')
#mlem = [x for x in range(10//2)]
#print(mlem)
#mlem = [v for v in range(len(vertices)//2)]
#mlem

#print(vertices)

#          print('\t%f\t%f\t%f\t%f' %(R.position[0],R.position[1],R.position[2],R.signal[w]),file=f)    #time signal

#with open("foo.py",'w') as f:
#    f.write('#\n')
#    f.write('')
