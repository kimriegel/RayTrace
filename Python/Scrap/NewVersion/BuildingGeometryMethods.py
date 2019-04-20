# William Costa

#Class/Methods Practice With Geometry Files

#FaceNormalNo=5
#FaceNormals = [(-1,0,0),(0,1,0),(1,0,0),(0,-1,0),(0,0,1)]
##^Will's Code
#import numpy as np
#Boxnumber=1
#Boxarraynear=np.zeros([Boxnumber,3])
#Boxarrayfar=np.zeros([Boxnumber,3])
#Boxarraynear[0]=[10,10,0]
#Boxarrayfar[0]=[64.4322,46.9316,8.2423]
#TriangleNumber=0
#SquareNumber=0
#PolyBuilding=0
import numpy as np
class box():
    """
    Geometry Class with attributes origin, boxnumber,normals, diagonal

    Initial box is 1x1x1.
    """
    def __init__(self,boxnumber,origin,endpoint):
        self.boxnumber=boxnumber
        self.origin=np.array(origin)
        self.endpoint=np.array(endpoint)
        #origin=np.array((0,0,0))
        #oxnumber=1
        self.normals=np.array(((-1,0,0),(0,1,0),(1,0,0),(0,-1,0),(0,0,1),(0,0,-1)))
        #endpoint=np.array((1,1,1))
        self.diagonal= self.endpoint-self.origin
    def edit_origin(self,x,y,z):
        '''
        Adjust origin coordinates by adding x,y,z values to current origin
        Adjusts diagonal attribute.
        '''

        self.origin=self.origin + np.array((x,y,z))
        self.diagonal=self.endpoint-self.origin
    def new_origin(self,x,y,z):
        '''
        Replace current origin with new origin (x,y,z)
        Adjusts diagonal attribute.
        '''
        self.origin=np.array((x,y,z))
        self.diagonal=self.endpoint-self.origin
    def edit_endpoint(self,x,y,z):
        '''
        Adjust the location of the diagonal endpoint.
        Adds (x,y,z) to current diagonal endpoint coordinates
        Adjusts diagonal attribute.
        '''
        self.endpoint=self.endpoint+np.array((x,y,z))
        self.diagonal=self.endpoint-self.origin
    def new_endpoint(self,x,y,z):
        '''
        Replace current diagonal endpoint with (x,y,z) coordinates.

        Adjusts diagonal attribute.
        '''
        self.endpoint=np.array((x,y,z))
        self.diagonal=self.endpoint-self.origin
    def inplane(self,x,y,z):
        '''
        Determines whether point(x,y,z) is in any of the planes of
        the box.
        '''
        naughts=np.array((self.origin,self.endpoint,self.endpoint,self.origin,self.endpoint,self.origin))
        point=np.array((x,y,z))
        for a in range(0,6):
            if np.dot(self.normals[a],(point-naughts[a]))==0 and ((any(point>self.endpoint)==False and any(point>self.origin)==True) or (any(point>self.endpoint)==True and any(point>self.origin==False))):
                print(point,"is in the plane",a+1)
                print(self.normals[a])
                planehit=True
            elif np.dot(self.normals[a],(point-naughts[a]))==0 and (point==self.endpoint or point==self.origin):
                print(point,"is in the plane",a+1)
                print(self.normals[a])
                planehit=True
            else:
                print(point, "is not on plane",a+1, "of this box")
                print(self.normals[a])
                planehit=False

building=box(1,[10,10,0],[64.4322,46.9316,8.2423])
#building.new_origin(10,10,0)
#building.new_endpoint(64.4322,46.9316,8.2423)
#print(building.origin[0],building.endpoint[1])
#boxx=box(1,[0,0,0],[1,1,1])
#print(boxx.origin,boxx.diagonal,boxx.endpoint)
