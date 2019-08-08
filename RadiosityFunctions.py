# $Revision: 1.0 $ : $Date:$ : $Author: kim

#***------------------------------------------------------------***
#***-------------------------- PATCHES -------------------------***
#***------------------------------------------------------------***

#     This function creates patch lengths using a geometric sum.

import numpy as np

def PATCHESSHORT(min,max,N,q,dd):
    int(N)
    k = ((max-min)/2)*(1-q)/(1-q**(N/2))
    #print(k)
    for m in range (N):
        if m <= (N/2):
            #print(k*q**(m-1))
            dd[m] = k*q**(m-1)
        elif (m <= N) and (m > (N/2)):
            dd[m] = k*q**(N-m)

def CREATEPATCHARRAY(ddm,ddL,Nm,Nl,origin,termius,patcharray,slope,b,slope1,b1,normal):
    """
    origin is x1,y1,z1
    termius is x2,y2,z2
    """
    int(Nm)
    int(Nl)
    #This function creates an array of patches for each plane
    d = -np.dot(origin,normal)

    ####################### Z Plane ############################
    if origin[2] == termius[2]:
        count = 0
        p =0
        print('ZPLANE', ddL, Nm, origin[0], origin[1], origin[2])
        for idxL, L in enumerate(ddL):
            for idxm, m in enumerate(ddm):
                x = origin[0] - m/2 + sum(ddm[:idxm])
                y = origin[1] - L/2 + sum(ddL[:idxL])
                z = (-d -normal[0]*x -normal[1]*y)/normal[2]
                if idxm == 0:
                    zcenter=z
                if idxL == 0:
                    ddz=z-origin[2]
                else:
                    ddz = 0.5*(z-zcenter)
                if (y < slope*x+b) and (y > slope1*x+b1):
                    patcharray[count] = [x,y,z,ddm[idxm],ddL[idxL],ddz]
                    count+=1
    ####################### Y Plane ############################
    elif origin[1] == termius[1]:
        count = 0
        print('YPLANE',ddL, Nm, origin[0], origin[1], origin[2])
        for idxL, L in enumerate(ddL):
            for idxm, m in enumerate(ddm):
                x = origin[0] - m/2 + sum(ddm[:idxm])
                z = origin[2] - L/2 + sum(ddL[:idxL])
                y = (-d -normal[0]*x-normal[2]*z)/normal[1]
                if idxm == 0:
                    ycenter = y
                if idxL == 0:
                    ddy = y-origin[1]
                else:
                    ddy = 0.5*(y-ycenter)
                if(z < slope*x+b) and (z > slope1*x+b1):
                    patcharray[count] = [x,y,z,ddm[idxm],ddy,ddL[idxL]]
                    count+=1

    ####################### X Plane ############################
    elif origin[0] == termius[0]:
        count=0
        print('XPLANE', ddL, Nm, origin[0], origin[1], origin[2])
        for idxL, L in enumerate(ddL):
            for idxm, m in enumerate(ddm):
                y = origin[1] - m/2 + sum(ddm[:idxm])
                z = origin[2] - L/2 + sum(ddL[:idxL])
                x = (-d -normal[2]*z-normal[1]*y)/normal[0]
                if idxm == 0:
                    xcenter=x
                if idxL == 0:
                    ddx= x-origin[0]
                else:
                    ddx= 0.5*(x-xcenter)
                if(z < slope*y+b) and (z > slope1*y+b1):
                    patcharray[count] = [x,y,z,ddx,ddm[idxm],ddL[idxL]]
                    count+=1

def PERPFORMFACTOR(patcharray,PatchNo,sizeffttwo,formfactors,Q,W,PI,FaceNormals,FaceNormalNo):
    S1 = np.zeros(3)
    S2 = np.zeros(3)

    dlnlm =np.sqrt((patcharray[Q,0,0]-patcharray[W,0,0])**2+(patcharray[Q,0,1]-patcharray[W,0,1])**2+(patcharray[Q,0,2]-patcharray[W,0,2])**2)
    formfactors = np.zeros((PatchNo,PatchNo,3))
    formfactors[Q,W,2]=dlnlm

    vec1=FaceNormals[int(patcharray[Q,0,9])]
    vec2=FaceNormals[int(patcharray[W,0,9])]

    length1=np.sqrt(vec1[0]**2+vec1[1]**2+vec1[2]**2)  
    length2=np.sqrt(vec2[0]**2+vec2[1]**2+vec2[2]**2)

    S1[0]=patcharray[W,0,0]-patcharray[Q,0,0]
    S1[1]=patcharray[W,0,1]-patcharray[Q,0,1]
    S1[2]=patcharray[W,0,2]-patcharray[Q,0,2]
    S1length=np.sqrt(S1[0]**2+S1[1]**2+S1[2]**2)

    S2[0]=patcharray[Q,0,0]-patcharray[W,0,0]
    S2[1]=patcharray[Q,0,1]-patcharray[W,0,1]
    S2[2]=patcharray[Q,0,2]-patcharray[W,0,2]
    S2length=np.sqrt(S2[0]**2+S2[1]**2+S2[2]**2)
 
    costheta1=(vec1[0]*S1[0]+vec1[1]*S1[1]+vec1[2]*S1[2])/(length1*S1length)
    costheta2=(vec2[0]*S2[0]+vec2[1]*S2[1]+vec2[2]*S2[2])/(length2*S2length)

    minimum=min(patcharray[W,0,3],patcharray[W,0,4],patcharray[W,0,5])

    if minimum == patcharray[W,0,3]:
        area = patcharray[W,0,4]*patcharray[W,0,5]       
    elif minimum == patcharray[W,0,4]:
        area = patcharray[W,0,3]*patcharray[W,0,5]
    elif minimum == patcharray[W,0,5]:
        area = patcharray[W,0,4]*patcharray[W,0,3]
    formfactors[Q,W,0]=(costheta1*costheta2)*area/(np.pi*S1length**2)