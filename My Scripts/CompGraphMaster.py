import numpy as np 
import matplotlib.pyplot as plt


time_array = np.array([])
pressure1_array = np.array([])
pressure2_array = np.array([])
pressure3_array = np.array([])
pressure4_array = np.array([])
pressure5_array = np.array([])

#These first two "with opens" open the first file selected, and then pulls data to create both an array 
#for both time and pressure which can then be plotted

with open('Constant_Comp_Simple0.06.txt', 'r' , encoding = 'utf8') as f:
    for line in f:
        if 'ZONE T=" Single Point "' in line:
            time = f.readline()
            time = time.split(' ')
            time = float(time[2])
            time_array = np.append(time_array, time)
f.close()           


with open('Constant_Comp_Simple0.06.txt', 'r', encoding = 'utf8' ) as f:
    for line in f:
        if 'DT=(SINGLE SINGLE SINGLE SINGLE )' in line:
            pressure = f.readline()
            pressure1 = pressure[30:39]
            pressure1 = float(pressure[30:39])
            pressure1_array = np.append(pressure1_array, pressure1)
            pressure = f.readline()
            pressure2 = float(pressure[30:39])
            pressure2_array = np.append(pressure2_array, pressure2)
            pressure = f.readline()
            pressure3 = float(pressure[30:39])
            pressure3_array = np.append(pressure3_array, pressure3)
            pressure = f.readline()
            pressure4 = float(pressure[30:39])
            pressure4_array = np.append(pressure4_array, pressure4)
            pressure = f.readline()
            pressure5 = float(pressure[30:39])
            pressure5_array = np.append(pressure5_array, pressure5)
f.close()

#The next two open statements pull data from the second file to do the exact same thing as 
#the first two loops, but this time for the FORTRAN code

time_arrayF = np.array([])
pressure1_arrayF = np.array([])
pressure2_arrayF = np.array([])
pressure3_arrayF = np.array([])
pressure4_arrayF = np.array([])
pressure5_arrayF = np.array([])

with open('Fortran_simpleBuildingComp.dat', 'r' , encoding = 'utf8') as f:
    for line in f:
        if 'ZONE T="Single Point        "' in line:
            time = f.readline()
            time = float(time[26:45])
            time_arrayF = np.append(time_arrayF, time)
f.close()           

with open('Fortran_simpleBuildingComp.dat', 'r', encoding = 'utf8' ) as f:
    for line in f:
        if 'DT=(SINGLE SINGLE SINGLE SINGLE )' in line:
            pressure = f.readline()
            pressure1 = float(pressure[51:68])
            pressure1_arrayF = np.append(pressure1_arrayF, pressure1)
            pressure = f.readline()
            pressure2 = float(pressure[51:68])
            pressure2_arrayF = np.append(pressure2_arrayF, pressure2)
            pressure = f.readline()
            pressure3 = float(pressure[51:68])
            pressure3_arrayF = np.append(pressure3_arrayF, pressure3)
            pressure = f.readline()
            pressure4 = float(pressure[51:68])
            pressure4_arrayF = np.append(pressure4_arrayF, pressure4)
            pressure = f.readline()
            pressure5 = float(pressure[51:68])
            pressure5_arrayF = np.append(pressure5_arrayF, pressure5)
f.close()





#Plotting Receiver1

plt.figure(layout = 'constrained')
plt.plot(time_array, pressure1_array, color = 'red', alpha = 0.5)
plt.plot(time_arrayF, pressure1_arrayF, color = 'blue', alpha = 0.5)
plt.title("Pressure vs Time of Receiver 1")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#Plotting Receiver 2 
plt.figure(layout = 'constrained')
plt.plot(time_array, pressure2_array, color = 'red', alpha = 0.5)
plt.plot(time_arrayF, pressure2_arrayF, color = 'blue', alpha = 0.5)
plt.title("Pressure vs Time of Receiver 2")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#Plotting Receiver 3 
plt.figure(layout = 'constrained')
plt.plot(time_array, pressure3_array, color = 'red', alpha = 0.5)
plt.plot(time_arrayF, pressure3_arrayF, color = 'blue', alpha = 0.5)
plt.title("Pressure vs Time of Receiver 3")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#Plotting Receiver 4 
plt.figure(layout = 'constrained')
plt.plot(time_array, pressure4_array, color = 'red', alpha = 0.5)
plt.plot(time_arrayF, pressure4_arrayF, color = 'blue', alpha = 0.5)
plt.title("Pressure vs Time of Receiver 4")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#Plotting Receiver 5
plt.figure(layout = 'constrained')
plt.plot(time_array, pressure5_array, color = 'red', alpha = 0.5)
plt.plot(time_arrayF, pressure5_arrayF, color = 'blue', alpha = 0.5)
plt.title("Pressure vs Time of Receiver 5")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()