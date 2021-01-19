####################################
#Needs Validation along with RadiosityFunction
#Needs Validation along with RadiosityFunction
#Needs Validation along with RadiosityFunction
#Needs Validation along with RadiosityFunction
import numpy as np
np.set_printoptions(threshold=np.inf)
import RadiosityFunctions as F 
import Parameterfile as Pf

FaceNormalNo = 5
FaceNormals = [(-1, 0, 0), (0, 1, 0), (1, 0, 0), (0, -1, 0), (0, 0, 1)]
# ^Will's Code

BoxNumber = 1
BoxArrayNear = np.zeros([BoxNumber, 3])
BoxArrayFar = np.zeros([BoxNumber, 3])
BoxArrayNear[0] = [10, 10, 0]
BoxArrayFar[0] = [64.4322, 46.9316, 8.2423]
TriangleNumber = 0
SquareNumber = 0
PolyBuilding = 0
sizeFFTTwo = 18000

PatchNox = 4
xlimit = np.zeros(PatchNox)
Nx = np.zeros(PatchNox-1)

xlimit[0]=0.0
xlimit[1]=10.0
xlimit[2]=64.4322
xlimit[3]=100.0


Nx[0]=10
Nx[1]=10
Nx[2]=10

qx=3.0
      
PatchNoy=4
ylimit = np.zeros(PatchNoy)
Ny = np.zeros(PatchNoy -1)

ylimit[0]=0.0
ylimit[1]=10.0
ylimit[2]=46.9316
ylimit[3]=60.0

Ny[0]=10
Ny[1]=10
Ny[2]=10

qy=3.0

PatchNoz=2
zlimit= np.zeros(PatchNoz)
Nz = np.zeros(PatchNoz-1)
zlimit[0]=0.0
zlimit[1]=8.2423

Nz[0]=10

qz=3.0

PatchNo = int(Nx[0]*Ny[0]+Nx[1]*Ny[0]+Nx[2]*Ny[0]+Nx[0]*Ny[1]+Nx[1]*Ny[1]+Nx[2]*Ny[1]+Nx[0]*Ny[2]+Nx[1]*Ny[2]+Nx[2]*Ny[2]+2*Nx[1]*Nz[0]+2*Ny[1]*Nz[0])
#PatchNo must be integer to put inside np.zeros
patcharray = np.zeros((PatchNo,sizeFFTTwo,10))
#patcharray has index of PatchNo with sizeFFTTwo amount of array in each with 10 zeros

#print(patcharray)
formfactors = np.zeros((PatchNo,PatchNo,3))
#patcharray = 0.0

#Let's say Dr.Riegel's diagram is the shape creating which I believe it is.
#Let's say when Planes go from 1 to 3, thats x-axis
#Let's say whne Planes go from 1 to 6, thats y-axis

#If these coords are correct, then the following code is 

W=0 #What is it? Placeholder
slope=0
slope1=0
a = 1
Q=0

#################################################### Plane 1 ##############################################################
ddx1 = np.zeros((int(Nx[0])),np.float16)
ddy1 = np.zeros((int(Ny[0])),np.float16)
increment=1
tempsize=Nx[0]*Ny[0]


patcharraytemp = np.zeros((int(Nx[0])*int(Ny[0]),6)) 
b=ylimit[1]
b1=ylimit[0]
origin=[xlimit[0],ylimit[0],zlimit[0]] 
terminal=[xlimit[1],ylimit[1],zlimit[0]]

F.PATCHESSHORT(xlimit[0],xlimit[1],int(Nx[0]),qx,ddx1)
F.PATCHESSHORT(ylimit[0],ylimit[1],int(Ny[0]),qy,ddy1)
F.CREATEPATCHARRAY(ddx1,ddy1,int(Nx[0]),int(Ny[0]),origin,terminal,patcharraytemp,slope,b,slope1,b1,FaceNormals[4])

while Q + increment < (increment + tempsize*a - 1):
        while W < sizeFFTTwo:
                patcharray[Q,W,0]=patcharraytemp[Q-increment+1,0]
                patcharray[Q,W,1]=patcharraytemp[Q-increment+1,1]
                patcharray[Q,W,2]=patcharraytemp[Q-increment+1,2]
                patcharray[Q,W,3]=patcharraytemp[Q-increment+1,3]
                patcharray[Q,W,4]=patcharraytemp[Q-increment+1,4]
                patcharray[Q,W,5]=patcharraytemp[Q-increment+1,5]
                patcharray[Q,W,6]=0.0
                patcharray[Q,W,7]=0.0
                patcharray[Q,W,8]=1.0
                patcharray[Q,W,9]=4.0
                W +=1
        W = 0
        Q += 1

a+=1

#################################################### Plane 2 ##############################################################
ddx1 = np.zeros((int(Nx[1])),np.float16)
ddy1 = np.zeros((int(Ny[0])),np.float16)
increment = Q
tempsize = Nx[1]*Ny[0]

patcharraytemp= np.zeros((int(Nx[1])*int(Ny[0]),6))
b = ylimit[1]
b1 = ylimit[0]
origin= [xlimit[1],ylimit[0],zlimit[0]]                      
terminal=[xlimit[2],ylimit[1],zlimit[0]]                     

F.PATCHESSHORT(xlimit[1],xlimit[2],int(Nx[1]),qx,ddx1)
F.PATCHESSHORT(ylimit[0],ylimit[1],int(Ny[0]),qy,ddy1)
F.CREATEPATCHARRAY(ddx1,ddy1,Nx[1],Ny[0],origin,terminal,patcharraytemp,slope,b,slope1,b1,FaceNormals[4])

while Q + increment < (increment + tempsize*a - 1):
        while W < sizeFFTTwo:
                patcharray[Q,W,0]=patcharraytemp[Q-increment,0]
                patcharray[Q,W,1]=patcharraytemp[Q-increment,1]
                patcharray[Q,W,2]=patcharraytemp[Q-increment,2]
                patcharray[Q,W,3]=patcharraytemp[Q-increment,3]
                patcharray[Q,W,4]=patcharraytemp[Q-increment,4]
                patcharray[Q,W,5]=patcharraytemp[Q-increment,5]
                patcharray[Q,W,6]=0.0
                patcharray[Q,W,7]=0.0
                patcharray[Q,W,8]=2.0
                patcharray[Q,W,9]=4.0
                W +=1
        W=0
        Q += 1
a+=1


################################################### Plane 3 ##############################################################
ddx1 = np.zeros((int(Nx[2])),np.float16)
ddy1 = np.zeros((int(Ny[0])),np.float16)
increment = Q
tempsize = Nx[2]*Ny[0]

patcharraytemp= np.zeros((int(Nx[2])*int(Ny[0]),6))
b = ylimit[1]
b1 = ylimit[0]
origin=[xlimit[2],ylimit[0],zlimit[0]]
terminal=[xlimit[3],ylimit[1],zlimit[0]]

F.PATCHESSHORT(xlimit[2],xlimit[3],int(Nx[2]),qx,ddx1)
F.PATCHESSHORT(ylimit[0],ylimit[1],int(Ny[0]),qy,ddy1)
F.CREATEPATCHARRAY(ddx1,ddy1,int(Nx[2]),int(Ny[0]),origin,terminal,patcharraytemp,slope,b,slope1,b1,FaceNormals[4])

while Q + increment < (increment + tempsize*a - 1):
        while W < sizeFFTTwo:
                patcharray[Q,W,0]=patcharraytemp[Q-increment,0]
                patcharray[Q,W,1]=patcharraytemp[Q-increment,1]
                patcharray[Q,W,2]=patcharraytemp[Q-increment,2]
                patcharray[Q,W,3]=patcharraytemp[Q-increment,3]
                patcharray[Q,W,4]=patcharraytemp[Q-increment,4]
                patcharray[Q,W,5]=patcharraytemp[Q-increment,5]
                patcharray[Q,W,6]=0.0
                patcharray[Q,W,7]=0.0
                patcharray[Q,W,8]=3.0
                patcharray[Q,W,9]=4.0
                W +=1
        W=0
        Q += 1
a+=1
#################################################### Plane 4 ##############################################################
ddx1=np.zeros((int(Nx[0])),np.float16)
ddy1=np.zeros((int(Ny[1])),np.float16)
increment = Q
tempsize = Nx[0] *Ny[1]

patcharraytemp= np.zeros((int(Nx[0])*int(Ny[1]),6))
b = ylimit[2]
b1 = ylimit[1]
origin=[xlimit[0],ylimit[1],zlimit[0]]
terminal=[xlimit[1],ylimit[2],zlimit[0]]

F.PATCHESSHORT(xlimit[0],xlimit[1],int(Nx[0]),qx,ddx1)
F.PATCHESSHORT(ylimit[1],ylimit[2],int(Ny[1]),qy,ddy1)
F.CREATEPATCHARRAY(ddx1,ddy1,int(Nx[0]),int(Ny[1]),origin,terminal,patcharraytemp,slope,b,slope1,b1,FaceNormals[4])
while Q + increment < (increment + tempsize*a - 1):
        while W < sizeFFTTwo:
                patcharray[Q,W,0]=patcharraytemp[Q-increment,0]
                patcharray[Q,W,1]=patcharraytemp[Q-increment,1]
                patcharray[Q,W,2]=patcharraytemp[Q-increment,2]
                patcharray[Q,W,3]=patcharraytemp[Q-increment,3]
                patcharray[Q,W,4]=patcharraytemp[Q-increment,4]
                patcharray[Q,W,5]=patcharraytemp[Q-increment,5]
                patcharray[Q,W,6]=0.0
                patcharray[Q,W,7]=0.0
                patcharray[Q,W,8]=4.0
                patcharray[Q,W,9]=4.0
                W +=1
        W=0
        Q += 1
a+=1
#################################################### Plane 5 ##############################################################
ddx1=np.zeros((int(Nx[2])),np.float16)
ddy1=np.zeros((int(Ny[1])),np.float16)
increment = Q
tempsize = Nx[2] * Ny[1]

patcharraytemp= np.zeros((int(Nx[2])*int(Ny[1]),6))
b = ylimit[2]
b1 = ylimit[1]
origin=[xlimit[2],ylimit[1],zlimit[0]]
terminal=[xlimit[3],ylimit[2],zlimit[0]]

F.PATCHESSHORT(xlimit[2],xlimit[3],int(Nx[2]),qx,ddx1)
F.PATCHESSHORT(ylimit[1],ylimit[2],int(Ny[1]),qy,ddy1)
F.CREATEPATCHARRAY(ddx1,ddy1,Nx[2],Ny[1],origin,terminal,patcharraytemp,slope,b,slope1,b1,FaceNormals[4])


while Q + increment < (increment + tempsize*a - 1):
        while W < sizeFFTTwo:
                patcharray[Q,W,0]=patcharraytemp[Q-increment,0]
                patcharray[Q,W,1]=patcharraytemp[Q-increment,1]
                patcharray[Q,W,2]=patcharraytemp[Q-increment,2]
                patcharray[Q,W,3]=patcharraytemp[Q-increment,3]
                patcharray[Q,W,4]=patcharraytemp[Q-increment,4]
                patcharray[Q,W,5]=patcharraytemp[Q-increment,5]
                patcharray[Q,W,6]=0.0
                patcharray[Q,W,7]=0.0
                patcharray[Q,W,8]=5.0
                patcharray[Q,W,9]=4.0
                W +=1
        W=0
        Q += 1

#################################################### Plane 6 ##############################################################
ddx1=np.zeros((int(Nx[0])),np.float16)
ddy1=np.zeros((int(Ny[2])),np.float16)
increment = Q
tempsize = Nx[0] * Ny[2]

patcharraytemp= np.zeros((int(Nx[0])*int(Ny[2]),6))
b = ylimit[3]
b1 = ylimit[2]
origin=[xlimit[0],ylimit[2],zlimit[0]]
terminal=[xlimit[1],ylimit[3],zlimit[0]]

F.PATCHESSHORT(xlimit[0],xlimit[1],int(Nx[0]),qx,ddx1)
F.PATCHESSHORT(ylimit[2],ylimit[3],int(Ny[2]),qy,ddy1)
F.CREATEPATCHARRAY(ddx1,ddy1,Nx[0],Ny[2],origin,terminal,patcharraytemp,slope,b,slope1,b1,FaceNormals[4])

while Q + increment < (increment + tempsize*a - 1):
        while W < sizeFFTTwo:
                patcharray[Q,W,0]=patcharraytemp[Q-increment,0]
                patcharray[Q,W,1]=patcharraytemp[Q-increment,1]
                patcharray[Q,W,2]=patcharraytemp[Q-increment,2]
                patcharray[Q,W,3]=patcharraytemp[Q-increment,3]
                patcharray[Q,W,4]=patcharraytemp[Q-increment,4]
                patcharray[Q,W,5]=patcharraytemp[Q-increment,5]
                patcharray[Q,W,6]=0.0
                patcharray[Q,W,7]=0.0
                patcharray[Q,W,8]=6.0
                patcharray[Q,W,9]=4.0
                W +=1
        W=0
        Q += 1
a+=1
#################################################### Plane 7 ##############################################################
ddx1=np.zeros((int(Nx[1])),np.float16)
ddy1=np.zeros((int(Ny[2])),np.float16)
increment = Q
tempsize = Nx[1] * Ny[2]

patcharraytemp= np.zeros((int(Nx[1])*int(Ny[2]),6))
b = ylimit[3]
b1 = ylimit[2]
origin=[xlimit[1],xlimit[2],zlimit[0]]
terminal=[xlimit[2],ylimit[3],zlimit[0]]

F.PATCHESSHORT(xlimit[1],xlimit[2],int(Nx[1]),qx,ddx1)
F.PATCHESSHORT(ylimit[2],ylimit[3],int(Ny[2]),qy,ddy1)
F.CREATEPATCHARRAY(ddx1,ddy1,Nx[1],Ny[2],origin,terminal,patcharraytemp,slope,b,slope1,b1,FaceNormals[4])


while Q + increment < (increment + tempsize*a - 1):
        while W < sizeFFTTwo:
                patcharray[Q,W,0]=patcharraytemp[Q-increment,0]
                patcharray[Q,W,1]=patcharraytemp[Q-increment,1]
                patcharray[Q,W,2]=patcharraytemp[Q-increment,2]
                patcharray[Q,W,3]=patcharraytemp[Q-increment,3]
                patcharray[Q,W,4]=patcharraytemp[Q-increment,4]
                patcharray[Q,W,5]=patcharraytemp[Q-increment,5]
                patcharray[Q,W,6]=0.0
                patcharray[Q,W,7]=0.0
                patcharray[Q,W,8]=7.0
                patcharray[Q,W,9]=4.0
                W +=1
        W=0
        Q += 1
a+=1
#################################################### Plane 8 ##############################################################
ddx1=np.zeros((int(Nx[2])),np.float16)
ddy1=np.zeros((int(Ny[2])),np.float16)
increment = Q
tempsize = Nx[2] * Ny[2]

patcharraytemp= np.zeros((int(Nx[2])*int(Ny[2]),6))
b = ylimit[3]
b1 = ylimit[2]
origin=[xlimit[2],ylimit[2],zlimit[0]]
terminal=[xlimit[3],ylimit[3],zlimit[0]]

F.PATCHESSHORT(xlimit[2],xlimit[3],int(Nx[2]),qx,ddx1)
F.PATCHESSHORT(ylimit[2],ylimit[3],int(Ny[2]),qy,ddy1)
F.CREATEPATCHARRAY(ddx1,ddy1,Nx[2],Ny[2],origin,terminal,patcharraytemp,slope,b,slope1,b1,FaceNormals[4])


while Q + increment < (increment + tempsize*a - 1):
        while W < sizeFFTTwo:
                patcharray[Q,W,0]=patcharraytemp[Q-increment,0]
                patcharray[Q,W,1]=patcharraytemp[Q-increment,1]
                patcharray[Q,W,2]=patcharraytemp[Q-increment,2]
                patcharray[Q,W,3]=patcharraytemp[Q-increment,3]
                patcharray[Q,W,4]=patcharraytemp[Q-increment,4]
                patcharray[Q,W,5]=patcharraytemp[Q-increment,5]
                patcharray[Q,W,6]=0.0
                patcharray[Q,W,7]=0.0
                patcharray[Q,W,8]=8.0
                patcharray[Q,W,9]=4.0
                W +=1
        W=0
        Q += 1
a+=1
#################################################### Plane 9 ##############################################################
ddx1=np.zeros((int(Nx[1])),np.float16)
ddz1=np.zeros((int(Nz[0])),np.float16)
increment = Q
tempsize = Nx[1] * Nz[0]

patcharraytemp= np.zeros((int(Nx[1])*int(Nz[0]),6))
b = zlimit[1]
b1 = zlimit[0]
origin=[xlimit[1],ylimit[1],zlimit[0]]
terminal=[xlimit[2],ylimit[1],zlimit[1]]

F.PATCHESSHORT(xlimit[1],xlimit[2],int(Nx[1]),qx,ddx1)
F.PATCHESSHORT(zlimit[0],zlimit[1],int(Nz[0]),qz,ddz1)
F.CREATEPATCHARRAY(ddx1,ddz1,Nx[1],Nz[0],origin,terminal,patcharraytemp,slope,b,slope1,b1,FaceNormals[3])


while Q + increment < (increment + tempsize*a - 1):
        while W < sizeFFTTwo:
                patcharray[Q,W,0]=patcharraytemp[Q-increment,0]
                patcharray[Q,W,1]=patcharraytemp[Q-increment,1]
                patcharray[Q,W,2]=patcharraytemp[Q-increment,2]
                patcharray[Q,W,3]=patcharraytemp[Q-increment,3]
                patcharray[Q,W,4]=patcharraytemp[Q-increment,4]
                patcharray[Q,W,5]=patcharraytemp[Q-increment,5]
                patcharray[Q,W,6]=0.0
                patcharray[Q,W,7]=0.0
                patcharray[Q,W,8]=9.0
                patcharray[Q,W,9]=3.0
                W +=1
        W=0
        Q += 1
a+=1
#################################################### Plane 10 #############################################################
ddy1=np.zeros((int(Ny[1])),np.float16)
ddz1=np.zeros((int(Nz[0])),np.float16)
increment = Q
tempsize=Ny[1]*Nz[0]

patcharraytemp= np.zeros((int(Ny[1])*int(Nz[0]),6))
b = zlimit[1]
b1 = zlimit[0]
origin=[xlimit[2],ylimit[1],zlimit[0]]
terminal=[xlimit[2],ylimit[2],zlimit[1]]

F.PATCHESSHORT(ylimit[1],ylimit[2],int(Ny[1]),qy,ddy1)
F.PATCHESSHORT(zlimit[0],zlimit[1],int(Nz[0]),qz,ddz1)
F.CREATEPATCHARRAY(ddy1,ddz1,Ny[1],Nz[0],origin,terminal,patcharraytemp,slope,b,slope1,b1,FaceNormals[2])

while Q + increment < (increment + tempsize*a - 1):
        while W < sizeFFTTwo:
                patcharray[Q,W,0]=patcharraytemp[Q-increment,0]
                patcharray[Q,W,1]=patcharraytemp[Q-increment,1]
                patcharray[Q,W,2]=patcharraytemp[Q-increment,2]
                patcharray[Q,W,3]=patcharraytemp[Q-increment,3]
                patcharray[Q,W,4]=patcharraytemp[Q-increment,4]
                patcharray[Q,W,5]=patcharraytemp[Q-increment,5]
                patcharray[Q,W,6]=0.0
                patcharray[Q,W,7]=0.0
                patcharray[Q,W,8]=10.0
                patcharray[Q,W,9]=2.0
                W +=1
        W=0
        Q += 1
a+=1
#################################################### Plane 11 #############################################################
ddy1=np.zeros((int(Ny[1])),np.float16)
ddz1=np.zeros((int(Nz[0])),np.float16)
increment=Q
tempsize=Ny[1]*Nz[0]

patcharraytemp= np.zeros((int(Ny[1])*int(Nz[0]),6))
b = zlimit[1]
b1 = zlimit[0]
origin=[xlimit[1],ylimit[1],zlimit[0]]
terminal=[xlimit[1],ylimit[2],zlimit[1]]

F.PATCHESSHORT(ylimit[1],ylimit[2],int(Ny[1]),qy,ddx1)
F.PATCHESSHORT(zlimit[0],zlimit[1],int(Nz[0]),qz,ddz1)
F.CREATEPATCHARRAY(ddy1,ddz1,Ny[1],Nz[0],origin,terminal,patcharraytemp,slope,b,slope1,b1,FaceNormals[0])

while Q + increment < (increment + tempsize*a - 1):
        while W < sizeFFTTwo:
                patcharray[Q,W,0]=patcharraytemp[Q-increment,0]
                patcharray[Q,W,1]=patcharraytemp[Q-increment,1]
                patcharray[Q,W,2]=patcharraytemp[Q-increment,2]
                patcharray[Q,W,3]=patcharraytemp[Q-increment,3]
                patcharray[Q,W,4]=patcharraytemp[Q-increment,4]
                patcharray[Q,W,5]=patcharraytemp[Q-increment,5]
                patcharray[Q,W,6]=0.0
                patcharray[Q,W,7]=0.0
                patcharray[Q,W,8]=11.0
                patcharray[Q,W,9]=0.0
                W +=1
        W=0
        Q += 1
a+=1
#################################################### Plane 12 #############################################################
ddx1=np.zeros((int(Nx[2])),np.float16)
ddz1=np.zeros((int(Ny[0])),np.float16)
increment=Q
tempsize=Nx[1] * Nz[0]

patcharraytemp= np.zeros((int(Nx[1])*int(Ny[0]),6))
b = zlimit[1]
b1 = zlimit[0]
origin=[xlimit[1],ylimit[2],zlimit[0]]
terminal=[xlimit[2],ylimit[2],zlimit[1]]

F.PATCHESSHORT(xlimit[1],xlimit[2],int(Nx[1]),qx,ddx1)
F.PATCHESSHORT(zlimit[0],zlimit[1],int(Nz[0]),qz,ddz1)
F.CREATEPATCHARRAY(ddx1,ddz1,Nx[1],Nz[0],origin,terminal,patcharraytemp,slope,b,slope1,b1,FaceNormals[1])

while Q + increment < (increment + tempsize*a - 1):
        while W < sizeFFTTwo:
                patcharray[Q,W,0]=patcharraytemp[Q-increment,0]
                patcharray[Q,W,1]=patcharraytemp[Q-increment,1]
                patcharray[Q,W,2]=patcharraytemp[Q-increment,2]
                patcharray[Q,W,3]=patcharraytemp[Q-increment,3]
                patcharray[Q,W,4]=patcharraytemp[Q-increment,4]
                patcharray[Q,W,5]=patcharraytemp[Q-increment,5]
                patcharray[Q,W,6]=0.0
                patcharray[Q,W,7]=0.0
                patcharray[Q,W,8]=12.0
                patcharray[Q,W,9]=1.0
                W +=1
        W=0
        Q += 1

a+=1
#################################################### Plane 13 #############################################################

ddy1=np.zeros((int(Ny[1])),np.float16)
ddx1=np.zeros((int(Nx[1])),np.float16)
increment = Q
tempsize=Ny[1]*Nx[1]

patcharraytemp= np.zeros((int(Ny[1])*int(Nx[1]),6))
b = ylimit[2]
b1 = ylimit[1]
origin=[xlimit[1],ylimit[1],zlimit[1]]
terminal=[xlimit[2],ylimit[2],zlimit[1]]

F.PATCHESSHORT(ylimit[1],ylimit[2],int(Ny[1]),qy,ddy1)
F.PATCHESSHORT(xlimit[1],xlimit[2],int(Nx[1]),qz,ddx1)
F.CREATEPATCHARRAY(ddx1,ddy1,Nx[1],Ny[1],origin,terminal,patcharraytemp,slope,b,slope1,b1,FaceNormals[4])

while Q + increment < (increment + tempsize*a - 1):
        while W < sizeFFTTwo:
                patcharray[Q,W,0]=patcharraytemp[Q-increment,0]
                patcharray[Q,W,1]=patcharraytemp[Q-increment,1]
                patcharray[Q,W,2]=patcharraytemp[Q-increment,2]
                patcharray[Q,W,3]=patcharraytemp[Q-increment,3]
                patcharray[Q,W,4]=patcharraytemp[Q-increment,4]
                patcharray[Q,W,5]=patcharraytemp[Q-increment,5]
                patcharray[Q,W,6]=0.0
                patcharray[Q,W,7]=0.0
                patcharray[Q,W,8]=13.0
                patcharray[Q,W,9]=4.0
                W +=1
        W=0
        Q += 1
a+=1

# print('end',Q)
# a = 300
# b = a + 200
# while a < b:
#         print(a)
#         print(patcharray[a,0])
#         print(patcharray[a,17999])
#         a+=1

Q = 0

while Q < PatchNo:
        while W < PatchNo:
                if patcharray[W,0,8] <= 7.0:
                        formfactors[Q,W,1] = 1.0
                if patcharray[W,0,8] >= 9.0:
                        formfactors[Q,W,1] = 2.0
                if patcharray[Q,0,8] == .10 and (patcharray[W,0,8] ==  8.0 or patcharray[W,0,8] == 9.0 or patcharray[W,0,8] == 13.0 or patcharray[W,0,8] == 14.0 or patcharray[W,0,8] == 15.0):
                        F.PERPFORMFACTOR(patcharray,PatchNo,sizeFFTTwo,formfactors,Q,W,np.pi,FaceNormals,FaceNormalNo)
                elif(patcharray[Q,0,8] == 2.0 and patcharray[W,0,8] == 12.0):
                        F.PERPFORMFACTOR(patcharray,PatchNo,sizeFFTTwo,formfactors,Q,W,np.pi,FaceNormals,FaceNormalNo)
                elif(patcharray[Q,0,8] == 3.0 and (patcharray[W,0,8] == 13.0 or patcharray[W,0,8] == 14.0)):
                        F.PERPFORMFACTOR(patcharray,PatchNo,sizeFFTTwo,formfactors,Q,W,np.pi,FaceNormals,FaceNormalNo)
                elif(patcharray[Q,0,8] == 4.0 and patcharray[W,0,8] == 15.0):
                        F.PERPFORMFACTOR(patcharray,PatchNo,sizeFFTTwo,formfactors,Q,W,np.pi,FaceNormals,FaceNormalNo)
                elif(patcharray[Q,0,8] == 5.0 and (patcharray[W,0,8] == 10.0 or patcharray[W,0,8] == 11.0 or patcharray[W,0,8] == 12.0 or patcharray[W,0,8]== 13.0 or patcharray[W,0,8]==14.0 or patcharray[W,0,8]== 15.0)):
                        F.PERPFORMFACTOR(patcharray,PatchNo,sizeFFTTwo,formfactors,Q,W,np.pi,FaceNormals,FaceNormalNo)
                elif(patcharray[Q,0,8]== 8.0 and (patcharray[W,0,8]==1.0)):
                        F.PERPFORMFACTOR(patcharray,PatchNo,sizeFFTTwo,formfactors,Q,W,np.pi,FaceNormals,FaceNormalNo)
                elif(patcharray[Q,0,8]==9.0 and patcharray[W,0,8]==1.0):
                        F.PERPFORMFACTOR(patcharray,PatchNo,sizeFFTTwo,formfactors,Q,W,np.pi,FaceNormals,FaceNormalNo)
                elif(patcharray[Q,0,8]==10.0 and patcharray[W,0,9] == 5.0 ):
                        F.PERPFORMFACTOR(patcharray,PatchNo,sizeFFTTwo,formfactors,Q,W,np.pi,FaceNormals,FaceNormalNo)
                elif(patcharray[Q,0,8]==11.0 and patcharray[W,0,8]==5.0):
                        F.PERPFORMFACTOR(patcharray,PatchNo,sizeFFTTwo,formfactors,Q,W,np.pi,FaceNormals,FaceNormalNo)
                elif(patcharray[Q,0,8]==12.0 and (patcharray[W,0,8]== 1.0 or patcharray[W,0,8]==2.0 or patcharray[W,0,8]==5.0)):
                        F.PERPFORMFACTOR(patcharray,PatchNo,sizeFFTTwo,formfactors,Q,W,np.pi,FaceNormals,FaceNormalNo)
                elif(patcharray[Q,0,8] == 13.0 and (patcharray[W,0,8]==1.0 or patcharray[W,0,8]==3.0 or patcharray[W,0,8]== 5.0)):
                        F.PERPFORMFACTOR(patcharray,PatchNo,sizeFFTTwo,formfactors,Q,W,np.pi,FaceNormals,FaceNormalNo)
                elif(patcharray[Q,0,8]==14.0 and (patcharray[W,0,8]==1.0 or patcharray[W,0,8]== 3.0 or patcharray[W,0,8]== 5.0 or patcharray[W,0,8]== 13.0)):
                        F.PERPFORMFACTOR(patcharray,PatchNo,sizeFFTTwo,formfactors,Q,W,np.pi,FaceNormals,FaceNormalNo)
                elif(patcharray[Q,0,9]== 15.0 and (patcharray[W,0,8]==1.0 or patcharray[W,0,8]== 4.0 or patcharray[W,0,8]==5.0)):
                        F.PERPFORMFACTOR(patcharray,PatchNo,sizeFFTTwo,formfactors,Q,W,np.pi,FaceNormals,FaceNormalNo)
                else:
                        formfactors[Q,W,0] = 0.0
                        formfactors[Q,W,2] = 0.0
                W+= 1
        W=0
        Q+=1
                        
                