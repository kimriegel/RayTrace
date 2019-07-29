import pywavefront as pwf
from pywavefront import ObjParser

scene=pwf.Wavefront('/Users/lovelace/Running_Version/SingleBuildingGeometry.obj')

#visualization.draw(scene)

building =ObjParser(scene,'/Users/lovelace/Running_Version/SingleBuildingGeometry.obj', strict=False, encoding="utf-8", create_materials=True, collect_faces=True, parse=True, cache=False)
print(building.normals)
vertices=[]
for vertex in range(0,len(building.wavefront.vertices)//2):
    vertices.append(building.wavefront.vertices[vertex])
print(vertices)
#print(building.wavefront.vertices)
building.parse_f
print(building.mesh.faces)
