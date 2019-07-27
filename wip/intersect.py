# just pseudocode
# don't run

#for face in range(0,len(subfaces)): # Using ray-plane algorithm from Haines chapter 3
for face in subfaces: # Using ray-plane algorithm from Haines chapter 3
    #A=subfaces[face][0]
    #B=subfaces[face][1]
    #C=subfaces[face][2]
    A = face[0]
    B = face[1]
    C = face[2]
    
    # Need to see a sample to work on this
    #self.V1=np.array(self.vertices[A]) #These create arrays of the vertices for the face
    #self.V2=np.array(self.vertices[B])
    #self.V3=np.array(self.vertices[C])

    L1=self.V2-self.V1 # calculates the two vectors using V1 as the reference vertex
    L2=self.V3-self.V1
    normal=np.cross(L1,L2)

    self.unitnormal=normal/np.sqrt(np.dot(normal,normal)) # calculates the normal vector to the plane

    D=np.dot(self.unitnormal,self.V1) # calculates plane equation D: Ax+By+Cz+D=0

    self.vd=np.dot(self.unitnormal,F) # dot product between normal and ray direction

    if self.vd==0: # ray is parallel to plane and no intersection occurs. ## special case??
        self.t=Huge #HOTFIX
        #pass
    else:
        v0=-(np.dot(self.unitnormal,veci)+D)
        self.t=v0/self.vd # distance from ray origin to plane intersection
        if self.t<=pf.h and self.t>0:
            distances=np.append(distances,self.t)