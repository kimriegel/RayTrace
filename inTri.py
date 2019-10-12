def inTri(triangle,p):

    a = triangle[0]
    b = triangle[1]
    c = triangle[2]

    s1 = c[1] - a[1]
    s2 = c[0] - a[0]
    s3 = b[1] - a[1]

    s4 = p[1] - a[1]

    w1 = ((a[0] * s1 + s4 * s2 - p[0] * s1) / 
            (s3 * s2 - (b[0] - a[0]) * s1))
    if s1 == 0:
        s1 = 1
    w2 = (s4 - w1 * s3) / s1

    return w1 >= 0 and w2 > 0 and (w1 + w2) <= 1

tri = ((-1,-1),(1,1),(1,-1))
#tri = ((-1,0),(0,1),(1,0))
p1 = (0.5,0.5)
p2 = (1,0.5)
p3 = (0.5,1)
p4 = (0,1) 
p5 = (-1,-1)

print(inTri(tri,p1))
print(inTri(tri,p2))
print(inTri(tri,p3))
print(inTri(tri,p4))
print(inTri(tri,p5))