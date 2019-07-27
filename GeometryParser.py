
import numpy as np
import pywavefront as pwf
#FaceNormals = [(-1,0,0),(0,1,0),(1,0,0),(0,-1,0),(0,0,1)]  Desired

def isHit(near,far,veci,F):
    """
    Returns a bool if the surface is hit or not
    Standin for the old box function for now
    Goal is to look for dxnear 
        aka dxbuilding
    
    This is very important to register that something is hit during a 
    step
    """


if __name__ == "__main__":

    ipname = 'SingleBuilding.obj'
    #ipname = 'TwoWalls.obj'
    ipfile = pwf.Wavefront(ipname)    # Read in geometry file
    env = pwf.ObjParser(ipfile,ipname, strict=False, encoding="utf-8", 
            create_materials=True, collect_faces=True, parse=True, cache=False)
    vertices = env.wavefront.vertices
    faces = env.mesh.faces
    #normals = env.normals               # Delete normal[5] needed
    #print(faces)
    #print(normals)
    FaceNormals = env.normals
    FaceNormalNo = len(FaceNormals)
    Boxnumber = 1     # supposed to import from s, come back to this later
        # Is this similar to Will's bands?
    Boxarraynear=np.array([10,10,0])
    Boxarrayfar= np.array([64.4322,46.9316,8.2423])
    
    #for v in vertices:
    #    print(v)
        #pass
    #print(faces)
    #print(vertices[4])
    #for f in faces:
    #    print(f)
    #for n in normals:
    #    print(n)

TriangleNumber=0
SquareNumber=0
PolyBuilding=0
FaceNormalNo=5
