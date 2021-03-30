# this function detects when it crosses a plane
# much of it is hardcoded, but that can all be fixes
import numpy as np
import pywavefront as pwf
from Parameterfile import ipname
#import Functions as f
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

# Move everything below here to Terrain

epsilon = 1e-6  # how small angle between ray and plane has to be to count as parallel

def collision_check2(face, veci, f):
    """find if a ray hits the face for our mesh function
    Only exists here because of weird python things. Will eventually move to Terrain
    """

    si = np.array([])
    tmp = np.zeros([3, len(face), 3])
    huge = 1000000.0
    n = face_normal_array(face)    # compute plane normal
    # Finding intersection [P]oint
    # parallel check  WE DON'T ACTUALLY CHECK FOR PARALLEL HERE!!!  WE SHOULD DO THAT
    nf = np.dot(n, f)        # rayDir in notes, plane normal dot F
    w = veci-np.array(face)[:, 2]
    si = -np.einsum('ij,ij->i', n, w)/nf         # data compression thing
    p = veci + si[:, np.newaxis]*f
    tmp[0] = np.subtract(p, np.array(face)[:, 0])
    tmp[1] = np.subtract(p, np.array(face)[:, 1])
    tmp[2] = np.subtract(p, np.array(face)[:, 2])
    a = np.cross((np.array(face)[:, 1]-np.array(face)[:, 0]), tmp[0, :])
    b = np.cross((np.array(face)[:, 2]-np.array(face)[:, 1]), tmp[1, :])
    c = np.cross((np.array(face)[:, 0]-np.array(face)[:, 2]), tmp[2, :])
    cond = (np.einsum('ij,ij->i', a, n) < 0) | (np.einsum('ij,ij->i', b, n) < 0) | (np.einsum('ij,ij->i', c, n) < 0)
    si[np.where(cond | (abs(nf) < epsilon) | (si < 0.0) | (si > step_size))] = huge
    # print('si',si)
    index = np.argmin(si)
    return si[index], n[index]

def face_normal_array(face):
    a = np.array(face)[:, 0]
    b = np.array(face)[:, 1]
    c = np.array(face)[:, 2]
    d = np.cross((b-a), (c-a))   # [D]irection
    e=np.sqrt(np.einsum('...i,...i', d, d))
    return d/e[:,None]