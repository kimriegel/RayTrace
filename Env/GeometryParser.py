# this function detects when it crosses a plane
# much of it is hardcoded, but that can all be fixes
import numpy as np
import pywavefront as pwf
from Parameterfile import ipname
from Parameterfile import h as step_size     # temporary, just used to make sure we do not overstep
# FaceNormals = [(-1,0,0),(0,1,0),(1,0,0),(0,-1,0),(0,0,1)]  Desired
epsilon = 1e-6  # how small angle between ray and plane has to be to count as parallel


def face_normal(face):
    a = np.array(face[0])
    b = np.array(face[1])
    c = np.array(face[2])
    d = np.cross((b-a), (c-a))   # [D]irection
#    e = d//np.sqrt(d.dot(d))
#    print('d1',d)
    return d


def edge_test(triangle, p, n):
    """
    Checks if some point P is inside a triangle, uses a given Normal
    """

    edge = ((triangle[1] - triangle[0]),
            (triangle[2] - triangle[1]),
            (triangle[0] - triangle[2]))

    chi = ((p - triangle[0]),
           (p - triangle[1]),
           (p - triangle[2]))

    sha = (n.dot(np.cross(edge[0], chi[0])) > 0,
           n.dot(np.cross(edge[1], chi[1])) > 0,
           n.dot(np.cross(edge[2], chi[2])) > 0)

    return np.all(sha)


def collision_check(face, veci, f):

    # find if a ray hits the face for our mesh function

    # the caps lock just reinforces that the variables are only used inside this function set

    huge = 1000000.0
#    F = np.array(F)         # hotfix
    n = face_normal(face)    # compute plane normal
    # Finding intersection [P]oint
    # parallel check
    n_f = n.dot(f)        # rayDir in notes, plane normal dot F
    is_parallel = (abs(n_f) < epsilon)    # bool, vD in old code
    # print('NF ',NF)
    if is_parallel:
        return huge        # ray does not hit, find an output to express that
    w = veci-face[2]
    si = -n.dot(w)/n_f
    # find distance between origin and intersect
    if si < 0:         # ray starts behind the face, break
        return huge, n
    elif si > step_size:
        return huge, n    # does not hit inside step, ignore it
    else:               # if and only if it hits within the step then
        p = veci + (si * f)
        a = np.cross(face[1]-face[0], p-face[0])
        b = np.cross(face[2]-face[1], p-face[1])
        c = np.cross(face[0]-face[2], p-face[2])
        if a.dot(n) < 0:
            return huge, n
        elif b.dot(n) < 0:
            return huge, n
        elif c.dot(n) < 0:
            return huge, n
        else:
            return si, n


def face_normal_array(face):
    a = np.array(face)[:, 0]
    b = np.array(face)[:, 1]
    c = np.array(face)[:, 2]
    d = np.cross((b-a), (c-a))  # [D]irection
    d_n=np.sqrt(np.einsum("ij,ij->i",d,d))
    e=np.zeros([len(d),3])
    for i in np.arange(len(d)):
        e[i] = np.divide(d[i],d_n[i] )
#    #print('d', len(d), len(np.sqrt(np.einsum("ij,ij->i", d, d))))
    return e


def collision_check2(face, veci, f):

    # find if a ray hits the face for our mesh function

    # the caps lock just reinforces that the variables are only used inside this function set

    si = np.array([])
    tmp = np.zeros([3, len(face), 3])
    huge = 1000000.0
    n = face_normal_array(face)    # compute plane normal
#    print('n',n)
    # Finding intersection [P]oint
    # parallel check
    nf = np.dot(n, f)        # rayDir in notes, plane normal dot F
#    si[(np.where((abs(NF) < epsilon)))] = HUGE
#    isParallel = (abs(NF) < epsilon)    # bool, vD in old code

#    if isParallel:
#        return HUGE        #ray does not hit, find an output to express that
    w = veci-np.array(face)[:, 2]
    si = -np.einsum('ij,ij->i', n, w)/nf
#    print('si',si)
    p = veci + si[:, np.newaxis]*f
#    print('p',p)
    tmp[0] = np.subtract(p, np.array(face)[:, 0])
    tmp[1] = np.subtract(p, np.array(face)[:, 1])
    tmp[2] = np.subtract(p, np.array(face)[:, 2])
    a = np.cross((np.array(face)[:, 1]-np.array(face)[:, 0]), tmp[0, :])
    b = np.cross((np.array(face)[:, 2]-np.array(face)[:, 1]), tmp[1, :])
    c = np.cross((np.array(face)[:, 0]-np.array(face)[:, 2]), tmp[2, :])
    cond = (np.einsum('ij,ij->i', a, n) < 0) | (np.einsum('ij,ij->i', b, n) < 0) | (np.einsum('ij,ij->i', c, n) < 0)
    si[np.where(cond | (abs(nf) < epsilon) | (si < 0.0))] = huge
    index = np.argmin(si)
#    print('si',si)
    return si[index], n[index]


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
