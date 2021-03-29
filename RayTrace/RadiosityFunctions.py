# $Revision: 1.0 $ : $Date:$ : $Author: kim
# ***------------------------------------------------------------***
# ***-------------------------- PATCHES -------------------------***
# ***------------------------------------------------------------***

#     This function creates patch lengths using a geometric sum.

import numpy as np


def patches_short(patch_min, patch_max, n, q, dd):
    int(n)
    k = ((patch_max-patch_min)/2)*(1-q)/(1-q**(n/2))
    for m in range(n):
        if m <= (n/2):
            dd[m] = k*q**(m-1)
        elif (m <= n) and (m > (n/2)):
            dd[m] = k*q**(n-m)


def create_patch_array(dd_m, dd_l, n_m, n_l, origin, termius, patch_array, slope, b, slope1, b1, normal):
    """
    origin is x1,y1,z1
    termius is x2,y2,z2
    """
    int(n_m)
    int(n_l)
    # This function creates an array of patches for each plane
    d = -np.dot(origin, normal)

    # ###################### Z Plane ############################
    if origin[2] == termius[2]:
        count = 0
        print('ZPLANE', dd_l, n_m, origin[0], origin[1], origin[2])
        for idxL, L in enumerate(dd_l):
            for idxm, m in enumerate(dd_m):
                x = origin[0] - m/2 + sum(dd_m[:idxm])
                y = origin[1] - L/2 + sum(dd_l[:idxL])
                z = (-d - normal[0]*x - normal[1]*y)/normal[2]
                if idxm == 0:
                    zcenter = z
                if idxL == 0:
                    ddz = z-origin[2]
                else:
                    ddz = 0.5*(z-zcenter)
                if (y < slope*x+b) and (y > slope1*x+b1):
                    patch_array[count] = [x, y, z, dd_m[idxm], dd_l[idxL], ddz]
                    count += 1
    # ###################### Y Plane ############################
    elif origin[1] == termius[1]:
        count = 0
        print('YPLANE', dd_l, n_m, origin[0], origin[1], origin[2])
        for idxL, L in enumerate(dd_l):
            for idxm, m in enumerate(dd_m):
                x = origin[0] - m/2 + sum(dd_m[:idxm])
                z = origin[2] - L/2 + sum(dd_l[:idxL])
                y = (-d - normal[0]*x-normal[2]*z)/normal[1]
                if idxm == 0:
                    ycenter = y
                if idxL == 0:
                    ddy = y-origin[1]
                else:
                    ddy = 0.5*(y-ycenter)
                if(z < slope*x+b) and (z > slope1*x+b1):
                    patch_array[count] = [x, y, z, dd_m[idxm], ddy, dd_l[idxL]]
                    count += 1

    # ###################### X Plane ############################
    elif origin[0] == termius[0]:
        count = 0
        print('XPLANE', dd_l, n_m, origin[0], origin[1], origin[2])
        for idxL, L in enumerate(dd_l):
            for idxm, m in enumerate(dd_m):
                y = origin[1] - m/2 + sum(dd_m[:idxm])
                z = origin[2] - L/2 + sum(dd_l[:idxL])
                x = (-d - normal[2]*z-normal[1]*y)/normal[0]
                if idxm == 0:
                    xcenter = x
                if idxL == 0:
                    ddx = x-origin[0]
                else:
                    ddx = 0.5*(x-xcenter)
                if(z < slope*y+b) and (z > slope1*y+b1):
                    patch_array[count] = [x, y, z, ddx, dd_m[idxm], dd_l[idxL]]
                    count += 1


def perp_form_factor(patch_array, patch_no, sizeffttwo, formfactors, q, w, pi, face_normals, face_normal_no):
    s_1 = np.zeros(3)
    s_2 = np.zeros(3)

    dlnlm = np.sqrt((patch_array[q, 0, 0]-patch_array[w, 0, 0])**2+(patch_array[q, 0, 1]-patch_array[q, 0, 1])**2 +
                    (patch_array[q, 0, 2]-patch_array[q, 0, 2])**2)
    formfactors = np.zeros((patch_no, patch_no, 3))
    formfactors[q, w, 2] = dlnlm

    vec1 = face_normals[int(patch_array[q, 0, 9])]
    vec2 = face_normals[int(patch_array[w, 0, 9])]

    length1 = np.sqrt(vec1[0]**2+vec1[1]**2+vec1[2]**2)
    length2 = np.sqrt(vec2[0]**2+vec2[1]**2+vec2[2]**2)

    s_1[0] = patch_array[w, 0, 0]-patch_array[q, 0, 0]
    s_1[1] = patch_array[w, 0, 1]-patch_array[q, 0, 1]
    s_1[2] = patch_array[w, 0, 2]-patch_array[q, 0, 2]
    s1_length = np.sqrt(s_1[0]**2+s_1[1]**2+s_1[2]**2)

    s_2[0] = patch_array[q, 0, 0]-patch_array[w, 0, 0]
    s_2[1] = patch_array[q, 0, 1]-patch_array[w, 0, 1]
    s_2[2] = patch_array[q, 0, 2]-patch_array[w, 0, 2]
    s2_length = np.sqrt(s_2[0]**2+s_2[1]**2+s_2[2]**2)
 
    costheta1 = (vec1[0]*s_1[0]+vec1[1]*s_1[1]+vec1[2]*s_1[2])/(length1*s1_length)
    costheta2 = (vec2[0]*s_2[0]+vec2[1]*s_2[1]+vec2[2]*s_2[2])/(length2*s2_length)

    minimum = min(patch_array[w, 0, 3], patch_array[w, 0, 4], patch_array[w, 0, 5])

    if minimum == patch_array[w, 0, 3]:
        area = patch_array[w, 0, 4]*patch_array[w, 0, 5]
    elif minimum == patch_array[w, 0, 4]:
        area = patch_array[w, 0, 3]*patch_array[w, 0, 5]
    elif minimum == patch_array[w, 0, 5]:
        area = patch_array[w, 0, 4]*patch_array[w, 0, 3]
    formfactors[q, w, 0] = (costheta1*costheta2)*area/(np.pi*s1_length**2)
