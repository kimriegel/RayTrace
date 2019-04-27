####################################
#unfinished
import numpy as np
import Functions as F 
PatchNox=4
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
ylimit = np.zero(PatchNoy)
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

PatchNo = Nx[0]*Ny[0]+Nx[1]*Ny[0]+Nx[2]*Ny[0]+Nx[0]*Ny[1]+Nx[1]*Ny[1]+Nx[2]*Ny[1]+Nx[0]*Ny[2]+Nx[1]*Ny[2]+Nx[2]*Ny[2]+2*Nx[1]*Nz[0]+2*Ny[1]*Nz[0]
patcharray = np.zeros((PatchNo,sizeffttwo,10))
formfactors = ((PatchNo,PatchNo,3))
patcharray = 0.0

#     This creates a patches for the Building. 
increment=1
slope=0
slope1=0
tempsize=Nx[0]*Ny[0]
ddx1 = np.zeros(Nx[0])
ddy1 = np.zeros(Ny[0])
patcharraytemp = np.zeros((Nx[0]*Ny[0],5)) 
b=ylimit[1]
b1=ylimit[0]
F.PATCHESSHORT(xlimit[0],xlimit[1],Nx[0],qx,ddx1)
F.PATCHESSHORT(ylimit[0],ylimit[1],Ny[0],qy,ddy1)
F.CREATEPATCHARRAY(ddx1,ddy1,Nx[0],Ny[0],xlimit[0],xlimit[1],ylimit[0],ylimit[1],zlimit[0],zlimit[0],patcharraytemp,slope,b,slope1,b1,count,FaceNormals[4,0:2])




while Q + increment < (increment + tempsize - 1):
        while W < sizeffttwo:
                patcharray[Q,W,0]=patcharraytemp[Q-increment,0]
                patcharray[Q,W,1]=patcharraytemp[Q-increment,1]
                patcharray[Q,W,2]=patcharraytemp[Q-increment,2]
                patcharray[Q,W,3]=patcharraytemp[Q-increment,3]
                patcharray[Q,W,4]=patcharraytemp[Q-increment,4]
                patcharray[Q,W,5]=patcharraytemp[Q-increment,5]
                patcharray[Q,W,6]=0.0
                patcharray[Q,W,7]=0.0
                patcharray[Q,W,8]=1.0
                patcharray[Q,W,9]=5.0
                W +=1
        Q += 1
np.clear(ddx1)
np.clear(ddy1)
increment = Q
tempsize = Nx[2] *Ny[0]
np.clear(patcharraytemp)
ddx1=np.zeros(Nx[2])
ddy1=np.zeros(Ny[0])
patcharraytemp= np.zeros((Nx[2]*Ny[0],5))
b = ylimit[1]
b1 = ylimit[0]

F.PATCHESSHORT(xlimit[2],xlimit[3],Nx[2],qx,ddx1)
F.PATCHESSHORT(ylimit[0],ylimit[1],Ny[0],qy,ddy1)
F.CREATEPATCHARRAY(ddx1,ddy1,Nx[2],Ny[0],xlimit[2],xlimit[3],ylimit[0],ylimit[1],zlimit[0],zlimit[0],patcharraytemp,slope,b,slope1,b1,count,FaceNormals[4,0:2])

 
while Q + increment < (increment + tempsize - 1):
        while W < sizeffttwo:
                patcharray[Q,W,0]=patcharraytemp[Q-increment,0]
                patcharray[Q,W,1]=patcharraytemp[Q-increment,1]
                patcharray[Q,W,2]=patcharraytemp[Q-increment,2]
                patcharray[Q,W,3]=patcharraytemp[Q-increment,3]
                patcharray[Q,W,4]=patcharraytemp[Q-increment,4]
                patcharray[Q,W,5]=patcharraytemp[Q-increment,5]
                patcharray[Q,W,6]=0.0
                patcharray[Q,W,7]=0.0
                patcharray[Q,W,8]=3.0
                patcharray[Q,W,9]=5.0
                W +=1
        Q += 1

np.clear(ddx1)
np.clear(ddy1)
increment = Q
tempsize = Nx[0] *Ny[1]
np.clear(patcharraytemp)
ddx1=np.zeros(Nx[0])
ddy1=np.zeros(Ny[1])
patcharraytemp= np.zeros((Nx[0]*Ny[1],5))
b = ylimit[2]
b1 = ylimit[1]

F.PATCHESSHORT(xlimit[0],xlimit[1],Nx[0],qx,ddx1)
F.PATCHESSHORT(ylimit[1],ylimit[2],Ny[1],qy,ddy1)
F.CREATEPATCHARRAY(ddx1,ddy1,Nx[0],Ny[1],xlimit[0],xlimit[1],ylimit[1],ylimit[2],zlimit[0],zlimit[0],patcharraytemp,slope,b,slope1,b1,count,FaceNormals[4,0:2])

while Q + increment < (increment + tempsize - 1):
        while W < sizeffttwo:
                patcharray[Q,W,0]=patcharraytemp[Q-increment,0]
                patcharray[Q,W,1]=patcharraytemp[Q-increment,1]
                patcharray[Q,W,2]=patcharraytemp[Q-increment,2]
                patcharray[Q,W,3]=patcharraytemp[Q-increment,3]
                patcharray[Q,W,4]=patcharraytemp[Q-increment,4]
                patcharray[Q,W,5]=patcharraytemp[Q-increment,5]
                patcharray[Q,W,6]=0.0
                patcharray[Q,W,7]=0.0
                patcharray[Q,W,8]=4.0
                patcharray[Q,W,9]=5.0
                W +=1
        Q += 1

np.clear(ddx1)
np.clear(ddy1)
increment = Q
tempsize = Nx[2] * Ny[1]
np.clear(patcharraytemp)
ddx1=np.zeros(Nx[2])
ddy1=np.zeros(Ny[1])
patcharraytemp= np.zeros((Nx[2]*Ny[1],5))
b = ylimit[2]
b1 = ylimit[1]

F.PATCHESSHORT(xlimit[2],xlimit[3],Nx[2],qx,ddx1)
F.PATCHESSHORT(ylimit[1],ylimit[2],Ny[1],qy,ddy1)
F.CREATEPATCHARRAY(ddx1,ddy1,Nx[2],Ny[1],xlimit[2],xlimit[3],ylimit[1],ylimit[2],zlimit[0],patcharraytemp,slope,b,slope1,b1,count,FaceNormals[4,0:2])

while Q + increment < (increment + tempsize - 1):
        while W < sizeffttwo:
                patcharray[Q,W,0]=patcharraytemp[Q-increment,0]
                patcharray[Q,W,1]=patcharraytemp[Q-increment,1]
                patcharray[Q,W,2]=patcharraytemp[Q-increment,2]
                patcharray[Q,W,3]=patcharraytemp[Q-increment,3]
                patcharray[Q,W,4]=patcharraytemp[Q-increment,4]
                patcharray[Q,W,5]=patcharraytemp[Q-increment,5]
                patcharray[Q,W,6]=0.0
                patcharray[Q,W,7]=0.0
                patcharray[Q,W,8]=5.0
                patcharray[Q,W,9]=5.0
                W +=1
        Q += 1



np.clear(ddx1)
np.clear(ddy1)
increment = Q
tempsize = Nx[0] * Ny[2]
np.clear(patcharraytemp)
ddx1=np.zeros(Nx[0])
ddy1=np.zeros(Ny[2])
patcharraytemp= np.zeros((Nx[0]*Ny[2],5))
b = ylimit[3]
b1 = ylimit[2]

F.PATCHESSHORT(xlimit[0],xlimit[1],Nx[0],qx,ddx1)
F.PATCHESSHORT(ylimit[2],ylimit[3],Ny[2],qy,ddy1)
F.CREATEPATCHARRAY(ddx1,ddy1,Nx[0],Ny[2],xlimit[0],xlimit[1],ylimit[2],ylimit[3],zlimit[0],zlimit[0],patcharraytemp,slope,b,slope1,b1,count,FaceNormals[4,0:2])

while Q + increment < (increment + tempsize - 1):
        while W < sizeffttwo:
                patcharray[Q,W,0]=patcharraytemp[Q-increment,0]
                patcharray[Q,W,1]=patcharraytemp[Q-increment,1]
                patcharray[Q,W,2]=patcharraytemp[Q-increment,2]
                patcharray[Q,W,3]=patcharraytemp[Q-increment,3]
                patcharray[Q,W,4]=patcharraytemp[Q-increment,4]
                patcharray[Q,W,5]=patcharraytemp[Q-increment,5]
                patcharray[Q,W,6]=0.0
                patcharray[Q,W,7]=0.0
                patcharray[Q,W,8]=6.0
                patcharray[Q,W,9]=5.0
                W +=1
        Q += 1

np.clear(ddx1)
np.clear(ddy1)
increment = Q
tempsize = Nx[1] * Ny[2]
np.clear(patcharraytemp)
ddx1=np.zeros(Nx[1])
ddy1=np.zeros(Ny[2])
patcharraytemp= np.zeros((Nx[1]*Ny[2],5))
b = ylimit[3]
b1 = ylimit[2]

F.PATCHESSHORT(xlimit[1],xlimit[2],Nx[1],qx,ddx1)
F.PATCHESSHORT(ylimit[2],ylimit[3],Ny[2],qy,ddy1)
F.CREATEPATCHARRAY(ddx1,ddy1,Nx[1],Ny[2],xlimit[1],xlimit[2],ylimit[2],ylimit[3],zlimit[0],zlimit[0],patcharraytemp,slope,b,slope1,b1,count,FaceNormals[4,0:2])

while Q + increment < (increment + tempsize - 1):
        while W < sizeffttwo:
                patcharray[Q,W,0]=patcharraytemp[Q-increment,0]
                patcharray[Q,W,1]=patcharraytemp[Q-increment,1]
                patcharray[Q,W,2]=patcharraytemp[Q-increment,2]
                patcharray[Q,W,3]=patcharraytemp[Q-increment,3]
                patcharray[Q,W,4]=patcharraytemp[Q-increment,4]
                patcharray[Q,W,5]=patcharraytemp[Q-increment,5]
                patcharray[Q,W,6]=0.0
                patcharray[Q,W,7]=0.0
                patcharray[Q,W,8]=7.0
                patcharray[Q,W,9]=5.0
                W +=1
        Q += 1

np.clear(ddx1)
np.clear(ddy1)
increment = Q
tempsize = Nx[2] * Ny[2]
np.clear(patcharraytemp)
ddx1=np.zeros(Nx[2])
ddy1=np.zeros(Ny[2])
patcharraytemp= np.zeros((Nx[2]*Ny[2],5))
b = ylimit[3]
b1 = ylimit[2]

F.PATCHESSHORT(xlimit[2],xlimit[3],Nx[2],qx,ddx1)
F.PATCHESSHORT(ylimit[2],ylimit[3],Ny[2],qy,ddy1)
F.CREATEPATCHARRAY(ddx1,ddy1,Nx[2],Ny[2],xlimit[2],xlimit[3],ylimit[2],ylimit[3],zlimit[0],zlimit[0],patcharraytemp,slope,b,slope1,b1,count,FaceNormals[4,0:2])

while Q + increment < (increment + tempsize - 1):
        while W < sizeffttwo:
                patcharray[Q,W,0]=patcharraytemp[Q-increment,0]
                patcharray[Q,W,1]=patcharraytemp[Q-increment,1]
                patcharray[Q,W,2]=patcharraytemp[Q-increment,2]
                patcharray[Q,W,3]=patcharraytemp[Q-increment,3]
                patcharray[Q,W,4]=patcharraytemp[Q-increment,4]
                patcharray[Q,W,5]=patcharraytemp[Q-increment,5]
                patcharray[Q,W,6]=0.0
                patcharray[Q,W,7]=0.0
                patcharray[Q,W,8]=8.0
                patcharray[Q,W,9]=5.0
                W +=1
        Q += 1

np.clear(ddx1)
np.clear(ddy1)
increment = Q
tempsize = Nx[1] * Ny[0]
np.clear(patcharraytemp)
ddx1=np.zeros(Nx[1])
ddy1=np.zeros(Ny[0])
patcharraytemp= np.zeros((Nx[1]*Ny[0],5))
b = ylimit[1]
b1 = ylimit[0]

F.PATCHESSHORT(xlimit[1],xlimit[2],Nx[1],qx,ddx1)
F.PATCHESSHORT(ylimit[0],ylimit[1],Ny[0],qy,ddy1)
F.CREATEPATCHARRAY(ddx1,ddy1,Nx[1],Nz[0],xlimit[1],xlimit[2],ylimit[1],ylimit[2],zlimit[0],zlimit[1],patcharraytemp,slope,b,slope1,b1,count,FaceNormals[3,0:2])

while Q + increment < (increment + tempsize - 1):
        while W < sizeffttwo:
                patcharray[Q,W,0]=patcharraytemp[Q-increment,0]
                patcharray[Q,W,1]=patcharraytemp[Q-increment,1]
                patcharray[Q,W,2]=patcharraytemp[Q-increment,2]
                patcharray[Q,W,3]=patcharraytemp[Q-increment,3]
                patcharray[Q,W,4]=patcharraytemp[Q-increment,4]
                patcharray[Q,W,5]=patcharraytemp[Q-increment,5]
                patcharray[Q,W,6]=0.0
                patcharray[Q,W,7]=0.0
                patcharray[Q,W,8]=9.0
                patcharray[Q,W,9]=4.0
                W +=1
        Q += 1

np.clear(ddx1)
np.clear(ddy1)
increment = Q
tempsize = Nx[1] * Ny[0]
np.clear(patcharraytemp)
ddx1=np.zeros(Nx[1])
ddy1=np.zeros(Ny[0])
patcharraytemp= np.zeros((Nx[1]*Ny[0],5))
b = ylimit[1]
b1 = ylimit[0]

F.PATCHESSHORT(xlimit[1],xlimit[2],Nx[1],qx,ddx1)
F.PATCHESSHORT(ylimit[1],ylimit[1],Ny[0],qy,ddy1)
F.CREATEPATCHARRAY(ddx1,ddy1,Nx[1],Nz[0],xlimit[1],xlimit[2],ylimit[1],ylimit[1],zlimit[0],zlimit[1],patcharraytemp,slope,b,slope1,b1,count,FaceNormals[3,0:2])

while Q + increment < (increment + tempsize - 1):
        while W < sizeffttwo:
                patcharray[Q,W,0]=patcharraytemp[Q-increment,0]
                patcharray[Q,W,1]=patcharraytemp[Q-increment,1]
                patcharray[Q,W,2]=patcharraytemp[Q-increment,2]
                patcharray[Q,W,3]=patcharraytemp[Q-increment,3]
                patcharray[Q,W,4]=patcharraytemp[Q-increment,4]
                patcharray[Q,W,5]=patcharraytemp[Q-increment,5]
                patcharray[Q,W,6]=0.0
                patcharray[Q,W,7]=0.0
                patcharray[Q,W,8]=9.0
                patcharray[Q,W,9]=4 .0
                W +=1
        Q += 1

#left off here