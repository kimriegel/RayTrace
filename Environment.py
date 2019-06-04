# William Costa
# 03/19/19: Creating an environment class to work in conjunction with RayTrace
# This program will take .obj files as input and retrieve the necessary data
# to conduct ray-environment interactions. All environments should be triangulated


# .obj files and .mtl files are needed to run the code to completion
import numpy as np
import pywavefront as pwf
from pywavefront import ObjParser
import Parameterfile_methods as pfm

# I think my environment class should be more malleable. There should be a set number of methods that all 
# environment objects should have, buildings and receivers alike, and then specific methods for special cases

class environment():
    """
    Uses .obj files as input and defines the environment as a series of triangulated faces

    bandwidth will be defined by two times the step length 'h' as defined in the parameter file,
    along whichever axis is sorted.

    fail conditions:
        -if the bandwidth is much greater than the maximum distance between the highest and lowest points, all
        lists generated will be of maximum size and speed will be compromised
        -the bandwidth becomes 0.0, in which case the vertices/faces will not be called properly  for ray interaction
        -duplicate vertices will be indexed differently, causing potential errors when calling the faces.
    current efforts:
        - determine appropriate general case parameters to prevent running at maximum array size
        - incorporate general ray-plane interaction
        - determine how best to incorporate receivers
    """

    def __init__(self,file_name):
        self.wavefront=pwf.Wavefront(file_name)
        environment=ObjParser(self.wavefront,file_name, strict=False, encoding="utf-8", create_materials=False, collect_faces=True, parse=True, cache=False)
        environment.parse_f
        self.normals=environment.normals
        self.vertices=environment.wavefront.vertices[0:len(environment.wavefront.vertices)//2]
        self.faces=environment.mesh.faces
        #self.vertlength=[0,len(self.vertices)]
        #self.sortvert=dict([len(self.vertices),self.vertices])
    #def intersection(self,ray,faces):
    def sortvert(self,vertices,axis):
        '''
        Sorts the list self.vertices into the list self.sortvert. List sorted by axis X-0, Y-1, Z-2
        '''
        self.sortvert=[]
        def axissort(elem):
            return elem[0][axis]
        for index in range(0,len(vertices)):
            self.sortvert.append([vertices[index],index])
            self.sortvert.sort(axissort)
        self.axismin=self.sortvert[0][0][axis]
        self.axismax=self.sortvert[len(self.sortvert)-1][0][axis]
        self.axisheight=self.axismax-self.axismin
        self.bandwidth=2*pfm.h
        
    def rayinteraction(self,ray,axis,divisions):
        '''
        creates 
        '''
        subvert=[]
        subfaces=[]
        if ray[axis]>self.axismax or ray[axis]<self.axismin: # if the ray is above or below the max/min, no interaction
            pass
        else:
            subvert=self.sortvert # creates a sorted subset of the vertices
            bandwidth=self.axisheight 
            for divide in range(0,divisions):
                if bandwidth <= self.bandwidth:
                    pass
                elif ray[axis]<subvert[len(subvert)//2][0][axis]:
                    subvert=subvert[0:len(subvert)//2]
                    #print(len(subvert))
                else:
                    # Rain drop     subvert chop chop
                    subvert=subvert[len(subvert)//2:len(subvert)]
                    #print(len(subvert))
                axismin=subvert[0][0][axis]
                axismax=subvert[len(subvert)-1][0][axis]
                axisheight=axismax-axismin
                bandwidth=axisheight
        for vertex in range(0,len(subvert)):
            vertindex=subvert[vertex][1]
            count=0
            for x in range(0,len(self.faces)):
                if self.faces[x] in subfaces: # This sanitation check is slow, alternative methods?
                    pass
                elif vertindex in self.faces[x]:
                    subfaces.append(self.faces[x])
        for face in range(0,len(subfaces)):
            A=subfaces[face][0]
            B=subfaces[face][1]
            C=subfaces[face][2]
            #print(self.sortvert)
        #print(bandwidth)
        #print(self.bandwidth)
        #print(axisheight)
        #print(self.axisheight)
        #print(subvert)
        #print(len(subvert))
        #print(subfaces)
        print(len(subfaces))
        return

environment=environment('/Users/lovelace/RayTrace/monkey.obj')
environment.sortvert(environment.vertices,2)
#print(len(environment.vertices))
#print(environment.vertices)
#print(environment.sortvert)
print(len(environment.faces))
environment.rayinteraction([10,20,0],2,100)
#print(environment.bandwidth)
#print(environment.sortvert)