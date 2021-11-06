# William Costa
# 03/19/19: Creating an environment class to work in conjunction with RayTrace
# This program will take .obj files as input and retrieve the necessary data
# to conduct ray-environment interactions. All environments should be triangulated


# .obj files and .mtl files are needed to run the code to completion
import numpy as np
import pywavefront as pwf
from pywavefront import ObjParser
import itertools
# import Parameterfile_methods


class Environment:

    # Uses .obj files as input and defines the environment as a series of triangulated faces

    def __init__(self, file_name):
        self.wavefront = pwf.Wavefront(file_name)
        environment = ObjParser(self.wavefront, file_name, strict=False, encoding="utf-8", create_materials=False,
                            collect_faces=True, parse=True, cache=False)
        # print(environment.parse_f)     # supposed to add .mtl file if it doesn't exist

        self.normals = environment.normals
        numVert = len(environment.wavefront.vertices)
        self.vertices = environment.wavefront.vertices[0:numVert//2]
        self.faces = environment.mesh.faces

    def IndexVertices(self):
        length = len(self.vertices)
        self.IndexedVertices = {}
        for i in range(1,length+1):
            self.IndexedVertices[i] = self.vertices[i-1]

    def IndexedFaces(self):
        length=len(self.faces)
        self.IndexedFaces={}
        for i in range(1,length+1):
            self.IndexedFaces[i] = self.faces[i-1]

    def SortVert(self, vertices, axis):

        # Sorts the list self.vertices into the list self.sortvert. List sorted by axis X-0, Y-1, Z-2

        self.sortedvertices = {k: v for k, v in sorted(self.IndexedVertices.items(), key=lambda item: item[1][axis])}
    def bands(self, bandNumber):
        zMin = self.sortedvertices[0][2]
        zMax = self.sortedvertices[-1][2]
        bandWidth = (zMax-zMin)/bandNumber
        bandDictionary = {}
        for i in range(1,bandNumber+1):
            for val,group in itertools.groupby(self.sortedvertices, lambda x: x[2] >= zMin+(bandWidth*(i-1)) and x[2] <= (zMin+(bandWidth*i))):
                if val:
                    bandDictionary[i]=list(group)

if __name__ == "__main__":          # What I'm writing now
    """
    when running file from here it will do this, else nothing
    """
    environ = Environment('Env/monkey.obj')

    environ.IndexVertices()
    environ.SortVert(environ.vertices,2)
    environ.IndexedFaces()
    print(environ.IndexedFaces)
    #environ.bands(10)
    #{10: (0.5, 0.09375, 0.6875)}

    #vertices = {1: (x1,y1,z1), 2: (x2,y2,z2)}
    #faces = {1: [vertices[1], vertices[2], vertices[2]]}
