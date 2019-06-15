# William Costa
# 03/19/19: Creating an environment class to work in conjunction with RayTrace
# This program will take .obj files as input and retrieve the necessary data
# to conduct ray-environment interactions. All environments should be triangulated


# .obj files and .mtl files are needed to run the code to completion
import numpy as np
import pywavefront as pwf
from pywavefront import ObjParser
import Parameterfile_methods as pfm

class environment():
    """
    Uses .obj files as input and defines the environment as a series of triangulated faces
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
            self.sortvert.sort(key=axissort)
        self.axismin=self.sortvert[0][0][axis]
        self.axismax=self.sortvert[len(self.sortvert)-1][0][axis]
        self.axisheight=self.axismax-self.axismin
        self.bandwidth=pfm.h*2 #sets bandwidth to 2x the step length
    def rayinteraction(self,veci,F,axis):
        '''
        creates 
        '''
        subvert=[]
        subfaces=[]
        count=0
        if veci[axis]>self.axismax or veci[axis]<self.axismin:
            pass
        else:
            subvert=self.sortvert
            bandwidth=self.axisheight
            while bandwidth >self.bandwidth:
                if veci[axis]<subvert[len(subvert)//2][0][axis]:
                    subvert=subvert[0:len(subvert)//2]
                else:
                    subvert=subvert[len(subvert)//2:len(subvert)]
                axismin=subvert[0][0][axis]
                axismax=subvert[len(subvert)-1][0][axis]
                axisheight=axismax-axismin
                bandwidth=axisheight
        for vertex in range(0,len(subvert)):
            vertindex=subvert[vertex][1]
            for x in range(0,len(self.faces)):
                if vertindex in self.faces[x]:
                    subfaces.append(self.faces[x])
        for face in range(0,len(subfaces)): # Using ray-plane algorithm from Haines chapter 3
            A=subfaces[face][0]
            B=subfaces[face][1]
            C=subfaces[face][2]
            V1=np.array(self.vertices[A]) #These create arrays of the vertices for the face
            V2=np.array(self.vertices[B])
            V3=np.array(self.vertices[C])
            L1=V2-V1 # calculates the two vectors using V1 as the reference vertex
            L2=V3-V1
            normal=np.cross(L1,L2) # calculates the normal vector to the plane
            D=np.dot(normal,V1) # calculates plane equation D: Ax+By+Cz+D=0
            vd=np.dot(normal,F) # dot product between normal and ray direction
            if vd==0: # ray is parallel to plane and no intersection occurs. ## special case??
                pass
            else:
                v0=-(np.dot(normal,veci)+D) 
                t=v0/vd # distance from ray origin to plane intersection
                if t<0: # ray intersection behind ray origin
                    pass
                else:
                    ri=veci+(F*t) # calculates ray intersection
                    print('ray-plane intersection =' , ri)
                    if vd<0: # Adjusts normal such that it points back towards ray-origin.
                        rn=normal
                        print('surface normal =' ,rn)
                    else:
                        rn=-normal
                        print('surface normal =' ,rn)

        return

#environment=environment('monkey.obj')
#environment.sortvert(environment.vertices,2)
#F=np.array([1,0,1])
#environment.rayinteraction([10,20,0],F,2)
