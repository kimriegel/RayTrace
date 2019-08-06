# William Costa
# 03/19/19: Creating an environment class to work in conjunction with RayTrace
# This program will take .obj files as input and retrieve the necessary data
# to conduct ray-environment interactions. All environments should be triangulated


# .obj files and .mtl files are needed to run the code to completion
import numpy as np
import pywavefront as pwf
from pywavefront import ObjParser
import Parameterfile as pf

HUGE = 1000000.0


def ground(height, point, normal):
    global gHeight, gABC,gNormal
    gHeight = -height
    gABC = point
    gNormal = normal
    return 

#def edgetest(aleph,bet):
#    nc = 0
#    if aleph[1]<0:
#        sh=-1
#    else:
#        sh=1
#    if bet[1]<0:
#        nsh=-1
#    else:
#        nsh=1
#    if sh!=nsh:
#        if aleph[0]>0 and bet[0]:
#            nc=nc+1
#        elif aleph[0]>0 or  bet[0]>0:
#            if aleph[0]-aleph[1]*(bet[0]-aleph[0])/(bet[1]-aleph[1])>0:
#                nc=nc+1
#        sh=nsh




# Order of Operations:
# 1. The ray location and direction are determined
# 2. The distance between the ray and an object is determined
# 3. If the object is within one step of the ray, the point of intersection between the ray and object is determined
# 4. THe point of intersection will become the new ray location, and the ray direction will be adjusted accordingly

class environment():
    """
    Uses .obj files as input and defines the environment as a series of triangulated faces

    bandwidth will be defined by two times the step length 'h' as defined in the parameter file,
    along whichever axis is sorted.
    """

    def __init__(self,file_name):
        self.wavefront=pwf.Wavefront(file_name)
        environment=ObjParser(self.wavefront,file_name, strict=False, encoding="utf-8", create_materials=False, collect_faces=True, parse=True, cache=False)
        environment.parse_f
        #self.vertices=environment.wavefront.vertices[0:len(environment.wavefront.vertices)//2]
        self.vertices=environment.wavefront.vertices[0:len(environment.wavefront.vertices)]
        #seeing if twice the faces prevents it from ignoring surfaces -G 7/24
        self.faces=environment.mesh.faces
        self.t=100
    def SortVertices(self,vertices,axis):
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
    def Boundaries(self):
        '''creates an array that defines the minimum and maximum axis values of the imported environment(building).
            in the format boundary=([xmin,xmax],[ymin,ymax],[zmin,zmax]) '''
        bounds=np.zeros((3,2))
        for vertex in self.sortvert:
            # x-minimum
            if vertex[0][0]<bounds[0][0]:
                bounds[0][0]=vertex[0][0]
            # x-maximum
            elif vertex[0][0]>bounds[0][1]:
                bounds[0][1]=vertex[0][0]
            else:
                pass
            #y-minimum
            if vertex[0][1]<bounds[1][0]:
                bounds[1][0]=vertex[0][1]
            #y-maximum
            elif vertex[0][1]>bounds[1][1]:
                bounds[1][1]=vertex[0][1]
            else:
                pass
            #z-minimum
            if vertex[0][2]<bounds[2][0]:
                bounds[2][0]=vertex[0][2]
            #z-maximum
            elif vertex[0][2]>bounds[2][1]:
                bounds[2][1]=vertex[0][2]
            else:
                pass
        return bounds

    def RayIntersection(self, veci,F, bounds):
        """
        veci is the ray position as defined in RayTrace.py
        F is the ray direction as defined in RayTrace.py
        """
        print('Calculating Ray Intersection')
        #print('veci', veci, 'F', F)
       # adjustment=[[-2*pf.h,+2*pf.h]]
       # shell = bounds+adjustment
       # if veci[0]<shell[0][0] or veci[0]>shell[0][1] or veci[1]<shell[1][0] or veci[1]>shell[1][1] or veci[2]<shell[2][0] or veci[2]>shell[2][1]:
       #     self.t=HUGE
       #     return self.t
        #stepped = veci+(F*pf.h)
        #print('stepped', stepped)
        zeros=np.zeros(3)
        subvert=[]
        subfaces=[]
        self.planes=[]
        distances=np.array([HUGE])
        rayaxis=0 # index used for (x,y,z) ordered ray coordinate
        if self.axis == 1:
            rayaxis=2
        elif self.axis == 2:
            rayaxis=1
        else:
            pass
        if veci[rayaxis]>self.axismax or veci[rayaxis]<self.axismin: # if the ray is above or below the max/min, no interaction
            print('no interaction')
            t=HUGE
            return t
        else:
            subvert=self.sortvert # creates a sorted subset of the vertices
            bandwidth=self.axisheight#establishes a bandwidth to be compared to self.bandwidth
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
                if vertindex in self.faces[x] and self.faces[x] not in subfaces:
                    subfaces.append(self.faces[x])
        for face in subfaces: # Using ray-plane algorithm from Haines chapter 3
            #print('Face', face)
            A=face[0]
            B=face[1]
            C=face[2]
            V1=np.array([self.vertices[A][0],self.vertices[A][2],self.vertices[A][1]]) #These create arrays of the vertices for the face, reordered to match, x,y,z format
            V2=np.array([self.vertices[B][0],self.vertices[B][2],self.vertices[B][1]])
            V3=np.array([self.vertices[C][0],self.vertices[C][2],self.vertices[C][1]])

            # Distance between veci and a point on the plane
            
            L1=V2-V1 # calculates the two vectors using V1 as the reference vertex
            L2=V3-V1
            normal=np.cross(L1,L2) +zeros # adding zero gets rid of negative zero error.
            unitnormal=normal/np.sqrt(np.dot(normal,normal)) # calculates the normal vector to the plane
            D=-1*np.dot(unitnormal,V2) # calculates plane equation D: Ax+By+Cz+D=0

            vd=np.dot(unitnormal,F) # dot product between normal and ray direction
            
            line1=V1-veci
            line2=V2-veci
            line3=V3-veci
            magnitude1=np.sqrt(np.dot(line1,line1))
            magnitude2=np.sqrt(np.dot(line2,line2))
            magnitude3=np.sqrt(np.dot(line3,line3))
            lines=[line1, line2, line3]
            magnitudes = [magnitude1, magnitude2, magnitude3]
            stepped=veci+(pf.h*F)
            print('stepped', stepped)
            step1=V1-stepped
            step2=V2-stepped
            step3=V3-stepped
            steps=[step1, step2, step3]
            print('steps +lines', step1, line1)

            if vd==0: # ray is parallel to plane and no intersection occurs. ## special case??
                t=HUGE #HOTFIX
                print('face', face, 'parallel, does not hit')
                pass
            else:

                v0=-(np.dot(unitnormal,veci)+D)
                #print('D',D,'self.vd', self.vd, 'v0', v0)
                t=v0/vd # distance from ray origin to plane intersection
                #rint('distance', t, 'plane', face)
                for magnitude in magnitudes:
                    if magnitude<=pf.h and magnitude1>0:
                        distances=np.append(distances,magnitude)
                        self.planes.append([face, vd, unitnormal, V1, V2, V3])
                #if t<=pf.h and t>0:
                #    distances=np.append(distances,t)
                #    self.planes.append([face, vd, unitnormal, V1, V2, V3])
                    #print('distances', distances)
                    #print('dxbuilding', min(distances))                                                                                                                                                            
        return min(distances)
    
        ###Because of the building shape and triangulation, there are 2 planes that are possibly hit
        # they share the same plane. RayHit should properly determine which of the two was hit. 
        
    
    def RayHit(self,veci,F,distance):
        ''' "RayHit is the one with intersection, you should put a 3-quote note" -G.K. Seaton '''
        #print('RayHit Began')
        count=0
        for plane in self.planes:
            count+=1
            face=plane[0]
            vd=plane[1]
            unitnormal=plane[2]
            v1=plane[3]
            v2=plane[4]
            v3=plane[5]
            if distance<0: # ray intersection behind ray origin
                print('Ray Intersection behind Ray Origin')
                pass
            else:
                adjustment=F*distance
                ri=veci+(F*distance) # calculates ray intersection
                if vd<0: # Adjusts normal such that it points back towards ray-origin.
                    rn=unitnormal
                else:
                    rn=-unitnormal
                dominant=np.argmax(abs(unitnormal))# Haines 3.2, coordinate w/ greatest magnitude
                uv1=np.delete(v1,dominant) # translation to UV coordinate
                uv2=np.delete(v2,dominant)
                uv3=np.delete(v3,dominant)
                riuv=np.delete(ri,dominant) # ray intersection UV coordinates
                uv1p=uv1-riuv #uv1prime, etc. adjusted ray intersection to coordinate system origin
                uv2p=uv2-riuv
                uv3p=uv3-riuv
                nc=0 #number of crossings
                sh=0 # sign holder
                nsh=0 # next sign holder
                # first edge test
                #print('Plane', count)
                #print('veci', veci)
                #print('first edge test')
                #print('unitnormal', unitnormal, 'absolute unit normal', abs(unitnormal))
                #print('dominant', dominant)
                #print('ri', ri, 'riuv', riuv)
                #print('face', face, 'triangle vertices', v1, v2, v3)
                #print('uv coordinates', uv1p, uv2p, uv3p)
                #print('uv1p', uv1p, 'uv2p', uv2p)
                if uv1p[1]<0:
                    sh=-1
                    #print('sh=-1')
                else:
                    sh=1
                    #print('sh=1')
                if uv2p[1]<0:
                    nsh=-1
                    #print('nsh=-1')
                else:
                    nsh=1
                    #print('nsh=1')
                if sh!=nsh:
                    if uv1p[0]>0 and uv2p[0]>0:
                        nc=nc+1
                        #print('nc=',nc)
                    elif uv1p[0]>0 or  uv2p[0]>0:
                        if uv1p[0]-uv1p[1]*(uv2p[0]-uv1p[0])/(uv2p[1]-uv1p[1])>0:
                            nc=nc+1
                            #print('nc=', nc)
                    sh=nsh
                    #print('sh=nsh=',sh, nsh, 'end first edge test')
                #second edge test
                #print('second edge test')
                #print('riuv', riuv)
                #print('uv2p', uv2p, 'uv3p', uv3p)
                if uv2p[1]<0:
                    sh=-1
                    #print('sh=-1')
                else:
                    sh=1
                    #print('sh=1')
                if uv3p[1]<0:
                    nsh=-1
                    #print('nsh=-1')
                else:
                    nsh=1
                    #print('nsh=1')
                if sh!=nsh:
                    if uv2p[0]>0 and uv3p[0]>0:
                        nc=nc+1
                        #print('nc=',nc)
                    elif uv2p[0]>0 or  uv3p[0]>0:
                        if uv2p[0]-uv2p[1]*(uv3p[0]-uv2p[0])/(uv3p[1]-uv2p[1])>0:
                            nc=nc+1
                            #print('nc=',nc)
                    sh=nsh
                    #print('sh=nsh=',sh,nsh,'end second edge test')
                #third edge test
                #print('third edge test')
                #print('riuv', riuv)
                #print('uv3p', uv3p, 'uv1p', uv1p)
                if uv3p[1]<0:
                    sh=-1
                    #print('sh=-1')
                else:
                    sh=1
                    #print('sh=1')
                if uv1p[1]<0:
                    nsh=-1
                    #print('nsh=-1')
                else:
                    nsh=1
                    #print('nsh=1')
                if sh!=nsh:
                    #print('BIGTEST')
                    #print('uv3p', uv3p, 'uv1p', uv1p)
                    if uv3p[0]>0 and uv1p[0]>0:
                        nc=nc+1
                        #print('nc=',nc)
                    elif uv3p[0]>0 or  uv1p[0]>0:
                        if uv1p[0]-uv3p[1]*(uv1p[0]-uv3p[0])/(uv1p[1]-uv3p[1])>0:
                            nc=nc+1
                           #print('nc=',nc)
                    sh=nsh
                    #print('sh=nsh=',sh,nsh,'end third edge test')
                if nc%2==0:
                    #print('number of crossings =', nc,'even, not in face', face)
                    #print('finished plane', count)
                    pass
                if nc%2!=0:
                    #print('ODD!!!!!!!!!!!!!!! number of crossings=', nc, 'odd, in face', face)

                    rn2=np.dot(rn,rn)
                    nbuilding=rn/np.sqrt(rn2)
                    dot1=np.dot(F,nbuilding)
                    F=F-(2.0*(dot1/rn2*nbuilding))
                    length=np.sqrt(np.dot(F,F))  
                    #print('veci', veci, 'ri', ri)    
                    print('finished plane', count)
            #print('RayHit Complete')
        return veci, F
#######################################################################333        
#triFace = range(3)
#for edge in range(triFace):
#    # bounce number
#    boz += edgetest(uvp[edge],uvp[(edge + 1) % 3])        # a side and its next side Iloops around
#
#                    #if __uv1p[0]__-uv3p[1]  #was bugged and may be covered by edge function now
############################################################
if __name__ == "__main__":
    # You can run main trace from here now
    import RayTrace_costa       #I'm kinda lazy -G
