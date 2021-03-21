import numpy as np
import time
import warnings
from matplotlib import pyplot as plt 
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning) 
begin = time.time() #begins time measurements

start = np.array([6.0,7.0, 500.0]) #start position
a_vec = start
h = 1.00 #step
to = 20 #starting temperature
alpha = -.0039 #Temperature coefficient 
co = 331.4 + 0.6*alpha* #co goes to cix, and ci
#v initial =  np.array([331+.6*t,]) 
# ci = np.array([co,co,co + alpha * start[2]])
phi = 0 #Angle phi position
theta = 4.36 # Angle theta position 4.36 originally
d = np.array([np.cos(phi)*np.sin(theta),np.sin(phi)*np.sin(theta), np.cos(theta)]) #Directional vector
print("The start location", start) # prints loop
n = np.array([0.0,0.0,1.0]) 
k = 0 #Sets our step counter
while k < 7000:
    final = start + h*d
    # Z_1 = start[2]
    # X_1 = start[0]
    ciz = c0[2] + (alpha *position[2])
    ci = np.array([co,co,co + alpha * final[2]])
    Zeta_f = d[2]/ci[2] - ((h*alpha)/(ci[2]**2))
    d = np.array(([d[0],d[1], ci[2]*Zeta_f ]) )
    #start issue
    if final[2] <= 0:
        d_building = -start[2]/d[2]
        final = start+(d_building*d)
        print("We have reached the floor", final) # Ray reaches floor    #Issue: Keeps returning back to previous loop
        d = d - 2*((np.dot(d,n))/((np.linalg.norm(n))**2))*n
        print("Dot product", np.dot(d,n))
        print("This is d",d)
        print("xx")
        print("Post bounce", final)
    print(final)
    k = k+1
    start = final
    plt.plot(final[0],final[2],'bo')
    # Z_2= final[2]
    # X_2 = final[0]
    # print("The slope is:", (Z_2-Z_1)/(X_2-X_1))
plt.xlabel('longitude')
plt.ylabel('Altitude')
plt.show()
time.sleep(1)
end = time.time()
print(f"Total runtime of the program is {end - begin}") 

#creating "on & off" feature. Use if statements, maybe .input?
#c = -.0039z + (.6To + 331.4) where paranthesis is intercept