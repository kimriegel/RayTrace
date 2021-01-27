# William Costa
# 03/19/19: Creating an environment class to work in conjunction with RayTrace
# This program will take .obj files as input and retrieve the necessary data
# to conduct ray-environment interactions. All environments should be triangulated


# .obj files and .mtl files are needed to run the code to completion
import numpy as np
import pywavefront as pwf
from pywavefront import ObjParser
# import Parameterfile_methods


class Environment:

    # Uses .obj files as input and defines the environment as a series of triangulated faces

    def __init__(self, file_name):
        self.wavefront = pwf.Wavefront(file_name)
        environment = ObjParser(self.wavefront, file_name, strict=False, encoding="utf-8", create_materials=False,
                            collect_faces=True, parse=True, cache=False)
        # print(environment.parse_f)     # supposed to add .mtl file if it doesn't exist

        self.normals = environment.normals
        self.vertices = environment.wavefront.vertices[0:len(environ.wavefront.vertices)//2]
        self.faces = environment.mesh.faces
        # self.vertlength = [0,len(self.vertices)]
        # self.sortvert=dict([len(self.vertices),self.vertices])
    # def intersection(self,ray,faces):

    def sortvert(self, vertices, axis):

        # Sorts the list self.vertices into the list self.sortvert. List sorted by axis X-0, Y-1, Z-2

        self.sortvert = []

        def axissort(elem):
            return elem[0][axis]
        for index in range(0, len(vertices)):
            self.sortvert.append([vertices[index], index])
            self.sortvert.sort(axissort)
        self.axismin = self.sortvert[0][0][axis]
        self.axismax = self.sortvert[len(self.sortvert)-1][0][axis]
        self.axisheight = self.axismax-self.axismin
        self.bandwidth = self.axisheight/256

    def rayinteraction(self, ray, axis, divisions):

        subvert = []
        subfaces = []
        count = 0
        if ray[axis] > self.axismax or ray[axis] < self.axismin:
            pass
        else:
            subvert = self.sortvert
            # print(len(subvert))
            bandwidth = self.axisheight
            for divide in range(0, divisions):
                if bandwidth <= 2*self.bandwidth:
                    pass
                elif ray[axis] < subvert[len(subvert)//2][0][axis]:
                    subvert = subvert[0:len(subvert)//2]
                    # print(len(subvert))
                else:
                    # Rain drop     subvert chop chop
                    subvert = subvert[len(subvert)//2:len(subvert)]
                    # print(len(subvert))
                axismin = subvert[0][0][axis]
                axismax = subvert[len(subvert)-1][0][axis]
                axisheight = axismax-axismin
                bandwidth = axisheight
        for vertex in range(0, len(subvert)):
            vertindex = subvert[vertex][1]
            for x in range(0, len(self.faces)):
                if vertindex in self.faces[x]:
                    subfaces.append(self.faces[x])
        for face in range(0, len(subfaces)):
            a = subfaces[face][0]
            b = subfaces[face][1]
            c = subfaces[face][2]
            # print(self.sortvert)

        print(subvert)
        # print(len(subvert))
        print(subfaces)
        # print(len(subfaces))
        return


if __name__ != "__main__":      # Old code

    environ = Environment('/Users/lovelace/Will Costa Version/monkey.obj')
    Environment.sortvert(environ.vertices, 2)
    Environment.rayinteraction([10, 20, 0], 2, 100)

if __name__ == "__main__":          # What I'm writing now
    """
    when running file from here it will do this, else nothing
    """
    # env = environment('EnvTest\SingleBuildingGeometry.obj')
    env = Environment('EnvTest/SingleBuilding.obj')
    # print(env.vertices)
    # trying to target specific vertex from face
    # print(env.faces)
    # print(env.faces[1][2])
    # print(env.vertices[env.faces[1][2]])        # works
    # testVert = env.vertices[env.faces[1][2]]
    # print(testVert)                             # output is value of specified vertex in list
    print(env.normals)
