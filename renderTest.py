# credit to Victoria Zhong for helping me find this algorithm 

#A,B,C are coordinates
# D is distance from (0,0,0) to plane, parallel to plane normal
# N  is the plane Normal
# v0 is the triangle's vertex
# 
# P = O + (t*R)
#
# D = N.dot(v0)     #########
#  
# P is point of intersection
# O is ray origin               veci in old code
# t is distance from O to P     dx in old code
# R is ray direction            F in old code
# 
# We are given O and R of the ray, also (a,b,c) and D. 
# P can be found with the known data
# We then find t using P,O, and R
# 
# t = - (N.dot(O) + D) / N.dot(R)   ###########
#
# We perform a check to see that t is less than the step length
# 
# ***Before t***
#   Compute N dot R 
#     if 0 then ray and plane are parallel
#     you would also be dividing by 0 
#
# if t < 0:
#   ray origin is behind the tri
#   no intersection
# 
#  



import numpy as np



# Triangle
#def triTest(v):
def triTest(tri,P,N):
    """
    triangle is a triple of the vertices that will be checked
    v is a triple of the coordinates inside the tri
    P is the intersection point with the ray

    returns where the intersection is inside the triangle or not
    """
    
    #edge = ((tri[1] - tri[0]), (tri[2] - tri[1]), (tri[0] - tri[2]))
    edge0 = tri[1] - tri[0]
    edge1 = tri[2] - tri[1]
    edge2 = tri[0] - tri[2]

    # tri[0] is v0 *****
    #C = ((P - tri[0]), (P - tri[1]), (P - tri[2]))
    c0 = P - tri[0]
    c1 = P - tri[1]
    c2 = P - tri[2]

    (N.dot( edge0.cross(C0) ), 
    N.dot( edge1.cross(C1) ), 
    N.dot( edge2.cross(C2) ))




edge0 = v1 - v0
edge1 = v2 - v1
edge2 = v0 - v2

c0 = P - v0
c1 = P - v1
c2 = P - v2

# just a placeholder for now
# eventually turn into an array to use linear math to speed this up
N.dot( edge0.cross(C0) )
N.dot( edge1.cross(C1) )
N.dot( edge2.cross(C2) )
