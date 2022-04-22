import numpy as np
import matplotlib.pyplot as plt
import random as ran
import math

ran.seed()
# define all inputs
#x_i=np.array([2.0,3.0,5.0])
x_i=np.array([2.0,3.0,5.0]) #initial position
theta=2.36 #Angle of ray
phi=0.0
#print(x_i)
d=np.array([np.cos(phi)*np.sin(theta),np.sin(theta)*np.sin(phi),np.cos(theta)])  #direction vector
nground=np.array([0,0,1])
magnitude = np.linalg.norm(d)
h=1.0  #step size
x_f=np.zeros((10,3))  #array of coordinates
#print(x_f[0])
#compute Ray Step
k=0
#diffusion coefficent
Diffusion=0.5
while (k<10):
    x_f[k]= x_i+h*d #First element of array
    #print('x_f',x_f,k)
    if (x_f[k,2] <= 0): #Ray has intersected with plane
        l = -x_i[2]/ d[2] #How far to hit the plane
        x_f[k] = x_i + l*d #Do not go past plane
        RanDiffuse = ran.random()
        print(d)
        if RanDiffuse <= Diffusion:
            z1=ran.random()
            z2=ran.random()
            theta=math.acos(math.sqrt(z1))
            phi=2*math.pi*z2
            d=np.array([np.cos(phi)*np.sin(theta),np.sin(theta)*np.sin(phi),np.cos(theta)])  #direction vector
            print("Diffuse")
        else:
            dot1 = np.dot(d, nground) #normal of the ground
            n2 = np.dot(nground, nground) #magnitude of the ground
            d -= (2.0 * (dot1 / n2 * nground))
        print(d)
       
        
    #else: l <= 0
    #continue  #Ray does not intersect with plane
    x_i = x_f[k]
    k+=1
    #print(x_f)

    print('k',k)


plt.plot(x_f[:,0],x_f[:,2])

plt.show()