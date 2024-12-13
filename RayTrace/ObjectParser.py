import numpy as np


def parse_obj_file(obj_file_path):
    vertices = []
    faces = []
    current_material = None
    face_material_mapping = {}  

    with open(obj_file_path, 'r') as file:
        for line in file:
            if line.startswith('usemtl'):
                current_material = line.strip().split()[1]
            elif line.startswith('v '):
                vertices.append(list(map(float, line.strip().split()[1:])))
            elif line.startswith('f '):
                face = []
                for vertex in line.strip().split()[1:]:
                    face.append(int(vertex.split('/')[0]) - 1)
                faces.append(face)
                face_material_mapping[len(faces) - 1] = current_material  # Associate the current face with its material

    mesh = [np.array([vertices[face[0]], vertices[face[1]], vertices[face[2]]]) for face in faces]
    # print("mesh Object Parser",mesh)
    return mesh, face_material_mapping

