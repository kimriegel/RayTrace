import numpy as np
import time
import warnings
from matplotlib import pyplot as plt 

warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning) 
begin = time.time() #begins time measurements

import matplotlib.font_manager as fm
# Font
stdfont = fm.FontProperties()
stdfont.set_family('serif')
stdfont.set_name('Times New Roman')
stdfont.set_size(20)

def rayBend(position, direction, n, alpha, zModThing):
    """
    Ray bending as originally designed by Christian Gomez
    Modified by George Seaton on 3/21/21 to get a feel for bending
    Input: Positon, direction, normal, alpha, and []
    Output: Position and direction
    """

    if (position[2]+ (alpha * zModThing)) <= 0:
        d_building = -position[2]/direction[2]
        #final = position+(d_building*d)
        position += (d_building* direction)
        direction -=  (2*((np.dot(direction,n))/((np.linalg.norm(n))**2))*n)
        #print("We have reached the floor", final) # Ray reaches floor    #Issue: Keeps returning back to previous loop
        #d = d - 2*((np.dot(d,n))/((np.linalg.norm(n))**2))*n
        #print("Dot product", np.dot(d,n))
        #print("This is d",d)
        #print("xx")
        #print("Post bounce", final)

    else:
        position[2] += (alpha * zModThing)
        direction[2] = ((direction[2]/position[2])-(h*alpha/position[2]**2) ) * position[2]
        #Zeta_f = d[2]/ci[2] - ((h*alpha)/(ci[2]**2))
        #d = np.array(([d[0],d[1], ci[2]*Zeta_f ]) )

    #yield position, direction      #generator
    return position, direction      #function


Imax= 7000

position = np.array([6.0,7.0, 500.0]) #start position
h = 1.00 #step
t0 = 20 #starting temperature
alpha = -.0039 #Temperature coefficient 
#alpha = -600 #Temperature coefficient 
co = 331.4 + 0.6*alpha #co goes to cix, and ci
#v initial =  np.array([331+.6*t,]) 
# ci = np.array([co,co,co + alpha * start[2]])
phi = 0 #Angle phi position
theta = 4.36 # Angle theta position 4.36 originally
d = np.array([np.cos(phi)*np.sin(theta),
              np.sin(phi)*np.sin(theta),
              np.cos(theta)]) #Directional vector
print("The start location", position) # prints loop
n = np.array([0.0,0.0,1.0]) 
xAxis = np.zeros(Imax)
yAxis = np.zeros(Imax)

for i in range(Imax):
    ciz = co + (alpha *position[2])                 #since only z of ci was used
    position, d = rayBend(position,d,n,alpha,ciz)
    xAxis[i] = position[0]                          #runs plt easier on my computer this way
    xAxis[i] = position[2]
    #plt.plot(position[0],position[2],'bo')         # kept just in case


                                                    #kept older code for comparison
#for k in range(700):
#    final = start + h*d     # position += dx*direction
#    # Z_1 = start[2]
#    # X_1 = start[0]
#    ci = np.array([co,co,co + alpha * final[2]])
#    Zeta_f = d[2]/ci[2] - ((h*alpha)/(ci[2]**2))
#    d = np.array(([d[0],d[1], ci[2]*Zeta_f ]) )
#    #start issue
#    if final[2] <= 0:
#        d_building = -start[2]/d[2]
#        final = start+(d_building*d)
#        print("We have reached the floor", final) # Ray reaches floor    #Issue: Keeps returning back to previous loop
#        d = d - 2*((np.dot(d,n))/((np.linalg.norm(n))**2))*n
#        print("Dot product", np.dot(d,n))
#        print("This is d",d)
#        print("xx")
#        print("Post bounce", final)
#    print(final)
#    start = final
#    plt.plot(final[0],final[2],'bo')
#    # Z_2= final[2]
#    # X_2 = final[0]
#    # print("The slope is:", (Z_2-Z_1)/(X_2-X_1))

end = time.time()
plt.figure(num=i, dpi=120, facecolor='#e6e6fa', edgecolor='r')  # lavender
plt.plot(xAxis, yAxis, '#780303')
plt.grid(True)
plt.xlabel('Longitude []',fontproperties=stdfont)
plt.ylabel('Altitude []' ,fontproperties=stdfont)
plt.title('Longitude vs Altitude of ray being tested ',
          fontproperties=stdfont,
          fontsize=16,
          fontweight='bold')
plt.show()
#time.sleep(1)
print(f"Total runtime of the program is {end - begin}") 

#creating "on & off" feature. Use if statements, maybe .input?
#c = -.0039z + (.6To + 331.4) where paranthesis is intercept


# Just here because we mentioned derivatives and integrals


def dydx_fwd(x,y):
	"""
	Forward differentiation
	"""
	
	y2 = y[:-1]
	y1 = y[1:]
	dy = y2-y1
	
	x2 = x[:-1]
	x1 = x[1:]
	dx = x2-x1

	dydx = dy/dx
	return dydx	

def dydx_ctr(x,y):
	"""
	Central differentiation
	"""
	
	y2 = y[:-1]
	y1 = y[1:]
	dy = y2-y1
	
	x2 = x[:-1]
	x1 = x[1:]
	dx = 2 * (x2-x1)

	dydx = dy/dx
	return dydx


def trapezoidal(a,b,fx):
    """
    a and b are the ranges to integrate
    n is the number of integrals to integrate along
    fx is the array of numbers to be integrated
    """
    n = len(fx) # calculate n using length of array instead of have it hardcoded
    dx = (b-a) / n # change along x axis
    fx_inner = np.sum(fx[1:-1]) * 2     # The inner range of numbers, which is the function multiplied by two
 
    fx_outer = fx[0] + fx[-1]   # The outer two numbers, which are the the function
 
    integral = (dx/2) * (fx_outer + fx_inner)
     
    return integral
 
def simpson(a,b,fx):
    """
    a and b are the ranges to integrate
    n is the number of integrals to integrate along
    """
    n = len(fx)             # calculate n using length of array instead of have it hardcoded
    dx = (b-a)/n            # change along x axis
    fx_outer = fx[0] + fx[-1]
 
    fx_4s = 0.0             #Sum of numbers multiplied by 4
    fx_2s = 0.0             #Sum of numbers multiplied by 2
 
    #Add all the numbers that multiply by 4
    for _ in range(1,n,2):    # from second place to penultimate
        fx_4s += (4*fx[_])
 
    #Add all the numbers that multiply by 2
    for _ in range(2,n-1,2):    
        fx_2s += (2*fx[_])
    integral = (dx/3) * (fx_outer + fx_2s +fx_4s )
 
    return integral