
import numpy as np
import pywavefront as pwf
FaceNormalNo=5
#FaceNormals = [(-1,0,0),(0,1,0),(1,0,0),(0,-1,0),(0,0,1)]  Desired
Boxarraynear=[10,10,0]
Boxarrayfar =[64.4322,46.9316,8.2423]

if __name__ == "__main__":

    ipname = 'SingleBuilding.obj'
    ipfile = pwf.Wavefront(ipname)    # Read in geometry file
    env = pwf.ObjParser(ipfile,ipname, strict=False, encoding="utf-8", 
            create_materials=True, collect_faces=True, parse=True, cache=False)
    #print(env)
    vertices = env.wavefront.vertices
    faces = env.mesh.faces
    normals = env.normals               # Delete normal[5] needed
    
    #for v in vertices:
    #    print(v)
        #pass
    #print(faces)
    #print(vertices[4])
    #for f in faces:
    #    print(f)
    #for n in normals:
    #    print(n)


Boxnumber=1
TriangleNumber=0
SquareNumber=0
PolyBuilding=0
