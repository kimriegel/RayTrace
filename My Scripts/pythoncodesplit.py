import numpy as np 
import matplotlib as mpl
import matplotlib.pyplot as plt


time_arrayP = np.array([])
pressure1_arrayP = np.array([])
pressure2_arrayP = np.array([])
pressure3_arrayP = np.array([])
pressure4_arrayP = np.array([])
pressure5_arrayP = np.array([])

#These first two "with opens" open the first file selected, and then pulls data to create both an array 
#for both time and pressure which can then be plotted

with open('Constant_Comp_Simple0.06.txt', 'r' , encoding = 'utf8') as f:
    for line in f:
        if 'ZONE T=" Single Point "' in line:
            time = f.readline()
            time = time.split(' ')
            time = float(time[2])
            time_arrayP = np.append(time_arrayP, time)
f.close()           

with open('Constant_Comp_Simple0.06.txt', 'r', encoding = 'utf8' ) as f:
    for line in f:
        if 'DT=(SINGLE SINGLE SINGLE SINGLE )' in line:
            #Receiver1
            pressure = f.readline()
            pressure1 = pressure.split('\t')
            pressure1 = float(pressure1[4])
            pressure1_arrayP = np.append(pressure1_arrayP, pressure1)
            #Receiver2 
            pressure = f.readline()
            pressure2 = pressure.split('\t')
            pressure2 = float(pressure2[4])
            pressure2_arrayP = np.append(pressure2_arrayP, pressure2)
            #Receiver3
            pressure = f.readline()
            pressure3 = pressure.split('\t')
            pressure3 = float(pressure3[4])
            pressure3_arrayP = np.append(pressure3_arrayP, pressure3)
            #Receiver4
            pressure = f.readline()
            pressure4 = pressure.split('\t')
            pressure4 = float(pressure4[4])
            pressure4_arrayP = np.append(pressure4_arrayP, pressure4)
            #Receiver5
            pressure = f.readline()
            pressure5 = pressure.split('\t')
            pressure5 = float(pressure5[4])
            pressure5_arrayP = np.append(pressure5_arrayP, pressure5)
f.close()


#Plotting Receiver1
plt.figure(layout = 'constrained')
plt.plot(time_arrayP, pressure1_arrayP, color = 'red', alpha = 0.5)
plt.title("Pressure vs Time of Receiver 1")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#Plotting Receiver 2 
plt.figure(layout = 'constrained')
plt.plot(time_arrayP, pressure2_arrayP, color = 'red', alpha = 0.5)
plt.title("Pressure vs Time of Receiver 2")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#Plotting Receiver 3 
plt.figure(layout = 'constrained')
plt.plot(time_arrayP, pressure3_arrayP, color = 'red', alpha = 0.5)
plt.title("Pressure vs Time of Receiver 3")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#Plotting Receiver 4 
plt.figure(layout = 'constrained')
plt.plot(time_arrayP, pressure4_arrayP, color = 'red', alpha = 0.5)
plt.title("Pressure vs Time of Receiver 4")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#Plotting Receiver 5
plt.figure(layout = 'constrained')
plt.plot(time_arrayP, pressure5_arrayP, color = 'red', alpha = 0.5)
plt.title("Pressure vs Time of Receiver 5")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()