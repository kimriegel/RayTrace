import numpy as np
import pywavefront as pwf
from Parameterfile import ipname
from Parameterfile import h as step_size

epsilon = 1e-6  # How small the angle between ray and plane has to be to count as parallel

def face_normal(face):
    a = np.array(face[0])
    b = np.array(face[1])
    c = np.array(face[2])
    return np.cross((b - a), (c - a))  # [D]irection

def collision_check(face, veci, f):
    huge = 1000000.0
    n = face_normal(face)
    n_f = np.dot(n, f)
    if abs(n_f) < epsilon:
        return huge, n
    w = veci - face[2]
    si = -np.dot(n, w) / n_f
    if si < 0 or si > step_size:
        return huge, n
    p = veci + (si * f)
    a = np.cross(face[1] - face[0], p - face[0])
    b = np.cross(face[2] - face[1], p - face[1])
    c = np.cross(face[0] - face[2], p - face[2])
    if a.dot(n) < 0 or b.dot(n) < 0 or c.dot(n) < 0:
        return huge, n
    return si, n

def collision_check2Test(mesh, veci, f):
    closest_distance = float('inf')
    closest_face_index = -1
    closest_normal = None
    #print("veci collision",veci)
    for i, face in enumerate(mesh):
        distance, normal = collision_check(face, veci, f)
        if distance < closest_distance:
            closest_distance = distance
            closest_face_index = i
            closest_normal = normal
    
    #print(f"closest_distance: {closest_distance}, closest_face_index: {closest_face_index}, closest_normal: {closest_normal}")
    return closest_distance, closest_face_index, closest_normal




# Initialize the mesh from the obj file
ipfile = pwf.Wavefront(ipname)
env = pwf.ObjParser(ipfile, ipname, strict=False, encoding="utf-8",
                    create_materials=True, collect_faces=True, parse=True, cache=False)
vertices = env.wavefront.vertices
faces = env.mesh.faces

mesh = [np.array((
    (vertices[f[0]][0], vertices[f[0]][1], vertices[f[0]][2]),
    (vertices[f[1]][0], vertices[f[1]][1], vertices[f[1]][2]),
    (vertices[f[2]][0], vertices[f[2]][1], vertices[f[2]][2])))
    for f in env.mesh.faces]
# print("mesh Geometry Parser",mesh)
