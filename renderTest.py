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

def faceNormal(face):
    a = np.array(face[0])
    b = np.array(face[1])
    c = np.array(face[2])
    d = np.cross((b-a),(c-a))   # [D]irection
    normal = d/np.sqrt(d.dot(d))
    # Where n is [n]ot [a] [n]umber, return 0, else return n
    return tuple(normal)

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

    # check if products are positive
    
    sha = ( N.dot(edge0.cross(c0)) > 0, 
            N.dot(edge1.cross(c1)) > 0, 
            N.dot(edge2.cross(c2)) > 0)

    # variables were not previously defined as arrays. 
    # This is a backup if they need to be, otherwise useless
    #sha =   ( np.dot(N,np.cross(edge0,c0) > 0,
    #        ( np.dot(N,np.cross(edge1,c1) > 0,
    #        ( np.dot(N,np.cross(edge2,c2) > 0))
    
    return np.all(sha) # return whether all conditions were met or not

def edgeTest(triangle,P,N):
    """     Checks if some point P is inside a triangle, uses a given Normal    """
    edge = np.array((triangle[1] - triangle[0]), 
                    (triangle[2] - triangle[1]), 
                    (triangle[0] - triangle[2]))

    chi = np.array( (P - triangle[0]), 
                    (P - triangle[1]), 
                    (P - triangle[2]))
    
    sha = ( N.dot(edge[0].cross(chi[0])) > 0, 
            N.dot(edge[1].cross(chi[1])) > 0, 
            N.dot(edge[2].cross(chi[2])) > 0)


# compute normal
def foo(FACE,VECI,F):
    """
    find if a ray hits the face for our mesh function
    the caps lock just reinforces that the variables are only used inside this function set
    """
    N = faceNormal(FACE)    # compute normal
    # find intersection [P]oint
        # parallel check
    eps = 0.01              # how small it has to be not to count, exists inside function for debugging
    NF = np.dot(N,F)        # rayDir in notes, plane normal dor F
    isParallel = (abs(NF) < eps)    # bool
    if isParallel:
        print('parallel')
        return               #ray does not hit, find an output to express that

    d = np.dot(N,FACE[0])   # is tri[0] and v0 in notes

    # find distance between origin and intersect
    t = (np.dot(N,VECI) + d) / NF
    if (t < 0):         # ray starts behind the face, break
        return

    p = VECI + (t * F)
    if (p > stepSize):
        return          # does not hit inside step, ignore it

    else:
        isHit = edgeTest()
    


# find intersection [P]oint
    #check if parallel
eps = 0.01              # how small it has to be not to count
rayDir = np.dot(N,F)    #plane normal dot F
# ^ this is n dot r   break here if 0
isParallel = (abs(rayDir) < eps )
# return if is Parallel or not

# find distance from origin to the face
d = np.dot(N,v0)

# find t
# the letter 'O' is a terrible variable name, why did someone write that?
t = (np.dot(N,origin) + d ) / rayDir
isBehind = (t < 0)      # check if ray is starting behind face
# break? if behind or?

# find intersection
p = origin + (t * F) 
# if (p > stepSize):
#    ignore this then

# if and only if p is within our step should we continue

#Now we find if p is inside face








stepSize = 10 # Pf.h


#edge0 = v1 - v0
#edge1 = v2 - v1
#edge2 = v0 - v2
#
#c0 = P - v0
#c1 = P - v1
#c2 = P - v2
#
## just a placeholder for now
## eventually turn into an array to use linear math to speed this up
#N.dot( edge0.cross(C0) )
#N.dot( edge1.cross(C1) )
#N.dot( edge2.cross(C2) )

testFace = ((0,0,0),(1,0,0),(1,1,0))
testRay = (0,0,1)
rayF = (0,0,-1)