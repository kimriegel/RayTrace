# William Costa
# 03/19/19: Creating an environment class to work in conjunction with RayTrace
# This program will take .obj files as input and retrieve the necessary data
# to conduct ray-environment interactions. All environments should be triangulated


# .obj files and .mtl files are needed to run the code to completion
import numpy as np
import pywavefront as pwf
from pywavefront import ObjParser
import Parameterfile as pf


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
        self.vertices=environment.wavefront.vertices[0:len(environment.wavefront.vertices)//2]
        self.faces=environment.mesh.faces
        self.t=100
    def sortvert(self,vertices,axis):
        '''
        Sorts the list self.vertices into the list self.sortvert. List sorted by axis X-0, Z-1, Y-2
        '''
        self.axis=axis
        self.sortvert=[]
        def axissort(elem):
            return elem[0][axis]
        for index in range(0,len(vertices)):
            self.sortvert.append([vertices[index],index])
            self.sortvert.sort(key=axissort)
        self.axismin=self.sortvert[0][0][axis]
        self.axismax=self.sortvert[len(self.sortvert)-1][0][axis]
        self.axisheight=self.axismax-self.axismin
        self.bandwidth=pf.h*2 #sets bandwidth to 2x the step length
    # 07/10/19: Considering dividing "rayinteraction" into two functions:
    # 1: RayIntersection: calculates the new Veci to set as dxbuilding
    # 2: RayHit: when the ray hits the plane, change the direction.
    # For the moment, this may simplify incorporating Environment into RayTrace.
    
    def RayIntersection(self, veci,F):
        '''
        veci is the ray position as defined in RayTrace.py
        F is the ray direction as defined in RayTrace.py
        '''
        subvert=[]
        subfaces=[]
        distances=np.array([100000])
        rayaxis=0 # index used for (x,y,z) ordered ray coordinate
        if self.axis == 1:
            rayaxis=2
        elif self.axis == 2:
            rayaxis=1
        else:
            pass
        if veci[rayaxis]>self.axismax or veci[rayaxis]<self.axismin: # if the ray is above or below the max/min, no interaction
            pass
        else:
            subvert=self.sortvert # creates a sorted subset of the vertices
            bandwidth=self.axisheight #establishes a bandwidth to be compared to self.bandwidth
            while bandwidth > self.bandwidth:
                if veci[rayaxis]<subvert[len(subvert)//2][0][self.axis]:
                    subvert=subvert[0:len(subvert)//2]
                else:
                    subvert=subvert[len(subvert)//2:len(subvert)]
                axismin=subvert[0][0][self.axis]
                axismax=subvert[len(subvert)-1][0][self.axis]
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
            self.V1=np.array(self.vertices[A]) #These create arrays of the vertices for the face
            self.V2=np.array(self.vertices[B])
            self.V3=np.array(self.vertices[C])
            L1=self.V2-self.V1 # calculates the two vectors using V1 as the reference vertex
            L2=self.V3-self.V1
            normal=np.cross(L1,L2)
            self.unitnormal=normal/np.sqrt(np.dot(normal,normal)) # calculates the normal vector to the plane
            D=np.dot(self.unitnormal,self.V1) # calculates plane equation D: Ax+By+Cz+D=0
            self.vd=np.dot(self.unitnormal,F) # dot product between normal and ray direction
            if self.vd==0: # ray is parallel to plane and no intersection occurs. ## special case??
                self.t=1000000 #HOTFIX
                pass
            else:
                v0=-(np.dot(self.unitnormal,veci)+D)
                self.t=v0/self.vd # distance from ray origin to plane intersection
                if self.t<=pf.h and self.t>0:
                    distances=np.append(distances,self.t)
        return min(distances)
    
    def RayHit(self,veci,F,distance):
        ''' "RayHit is the one with intersection, you should put a 3-quote note" -G.K. Seaton '''

        if distance<0: # ray intersection behind ray origin
            print('Ray Intersection behind Ray Origin')
            pass
        else:
            adjustment=F*distance
            ri=veci+(F*distance) # calculates ray intersection
            if self.vd<0: # Adjusts normal such that it points back towards ray-origin.
                rn=self.unitnormal
            else:
                rn=-self.unitnormal
            dominant=np.argmax(self.unitnormal) # Haines 3.2, coordinate w/ greatest magnitude
            uv1=np.delete(self.V1,dominant) # translation to UV coordinates
            uv2=np.delete(self.V2,dominant)
            uv3=np.delete(self.V3,dominant)
            riuv=np.delete(ri,dominant) # ray intersection UV coordinates
            uv1p=uv1-riuv #uv1prime, etc. adjusted ray intersection to coordinate system origin
            uv2p=uv2-riuv
            uv3p=uv3-riuv
            nc=0 #number of crossings
            sh=0 # sign holder
            nsh=0 # next sign holder
            # first edge test
            if uv1p[1]<0:
                sh=-1
            else:
                sh=1
            if uv2p[1]<0:
                nsh=-1
            else:
                nsh=1
            if sh!=nsh:
                if uv1p[0]>0 and uv2p[0]:
                    nc=nc+1
                elif uv1p[0]>0 or  uv2p[0]>0:
                    if uv1p[0]-uv1p[1]*(uv2p[0]-uv1p[0])/(uv2p[1]-uv1p[1])>0:
                        nc=nc+1
                sh=nsh
            #second edge test
            if uv2p[1]<0:
                sh=-1
            else:
                sh=1
            if uv3p[1]<0:
                nsh=-1
            else:
                nsh=1
            if sh!=nsh:
                if uv2p[0]>0 and uv3p[0]:
                    nc=nc+1
                elif uv2p[0]>0 or  uv3p[0]>0:
                    if uv2p[0]-uv2p[1]*(uv3p[0]-uv2p[0])/(uv3p[1]-uv2p[1])>0:
                        nc=nc+1
                sh=nsh
            #third edge test
            if uv3p[1]<0:
                sh=-1
            else:
                sh=1
            if uv1p[1]<0:
                nsh=-1
            else:
                nsh=1
            if sh!=nsh:
                if uv3p[0]>0 and uv1p[0]:
                    nc=nc+1
                elif uv3p[0]>0 or  uv1p[0]>0:
                    if uv1p[0]-uv3p[1]*(uv1p[0]-uv3p[0])/(uv1p[1]-uv3p[1])>0:
                        nc=nc+1
                sh=nsh
            if nc%2==0:
                pass
            if nc%2!=0:
                rn2=np.dot(rn,rn)
                nbuilding=rn/np.sqrt(rn2)
                dot1=np.dot(F,nbuilding)
                F=F-(2.0*(dot1/rn2*nbuilding))
                length=np.sqrt(np.dot(F,F))  
                print('veci', veci, 'ri', ri)    
        return veci, F