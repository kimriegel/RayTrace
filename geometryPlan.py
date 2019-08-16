"""
__Phase 1__

read in from file
translate to list of objects
    List of vertices
    list of faces to be used in geometry checking
    Face contains vertex coordinates and normal
        normal can be calculated
First check plane intersection
Then check if within ranges
Finally check normals
"""

"""
__Phase 2__

Divide geometries into bands
Geometry outputs to intermidiate graph to allow code to work w complex geometry
"""


#Objects contain vertices, faces and normals
#Divide list based on position of faces
#save to output file
#"""
# saving to output file allows checks where different numbers of rays are used
# saving also avoids bloat from procedurally generating turbulence

# does the parser /need/ to cut list in half?
"""""""""""
 Done
"""""""""""

"""
After casting then save data to file
    Make sure file can be easily read in python
Create graphs based on output file
    Has trouble doing multiple graphs at once, look into it
"""

