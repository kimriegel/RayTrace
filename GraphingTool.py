# This program will take data from two different files
# one from Python and one from Fortran and graph them against each other for comparison.
# The Python file should be read from the top and aligned with the "P" arrays.
# The Fortran file should be read second and aligned with the "F" arrays.

import numpy as np
import matplotlib.pyplot as plt

# Python Reader 1 Code 


# Creating the neccesary blank arrays for the program to work, one array for the time values and one for each of the
# five pressure receivers
time_arrayP = np.array([])
pressure1_arrayP = np.array([])
pressure2_arrayP = np.array([])
pressure3_arrayP = np.array([])
pressure4_arrayP = np.array([])
pressure5_arrayP = np.array([])

# These first two "with opens" open the first file selected, and then pulls data to create both an array
# for both time and pressure which can then be plotted

# This loops finds the line before the neccesary values, readline then moves to the line with the values,
# uses the split function to skip the spacing between the string values read back, time[2] is the time value
# floats that value, then appends it to the time array created before.

with open('StratPython_comparison_0.06.txt', 'r', encoding='utf8') as f:
    for line in f:
        if 'ZONE T=" Single Point "' in line:
            time = f.readline()
            time = time.split(' ')
            time = float(time[2])
            time_arrayP = np.append(time_arrayP, time)
f.close()

# This loop does the same thing as the one above, but uses the loop to find the pressure values for
# each receiver, and then appends them into separate arrays for each.

with open('StratPython_comparison_0.06.txt', 'r', encoding='utf8') as f:
    for line in f:
        if 'DT=(SINGLE SINGLE SINGLE SINGLE )' in line:
            # Receiver1
            pressure = f.readline()
            pressure1 = pressure.split('\t')
            pressure1 = float(pressure1[4])
            pressure1_arrayP = np.append(pressure1_arrayP, pressure1)
            print(pressure1_arrayP)
            # Receiver2
            pressure = f.readline()
            pressure2 = pressure.split('\t')
            pressure2 = float(pressure2[4])
            pressure2_arrayP = np.append(pressure2_arrayP, pressure2)
            # Receiver3
            pressure = f.readline()
            pressure3 = pressure.split('\t')
            pressure3 = float(pressure3[4])
            pressure3_arrayP = np.append(pressure3_arrayP, pressure3)
            # Receiver4
            pressure = f.readline()
            pressure4 = pressure.split('\t')
            pressure4 = float(pressure4[4])
            pressure4_arrayP = np.append(pressure4_arrayP, pressure4)
            # Receiver5
            pressure = f.readline()
            pressure5 = pressure.split('\t')
            pressure5 = float(pressure5[4])
            pressure5_arrayP = np.append(pressure5_arrayP, pressure5)
f.close()
# End of Python Reader 1

# Python Reader 2

time_arrayP2 = np.array([])
pressure1_arrayP2 = np.array([])
pressure2_arrayP2 = np.array([])
pressure3_arrayP2 = np.array([])
pressure4_arrayP2 = np.array([])
pressure5_arrayP2 = np.array([])

with open('StratPython_comparison_0.1.txt', 'r', encoding='utf8') as f:
    for line in f:
        if 'ZONE T=" Single Point "' in line:
            time = f.readline()
            time = time.split(' ')
            time = float(time[2])
            time_arrayP2 = np.append(time_arrayP2, time)
f.close()

with open('StratPython_comparison_0.1.txt', 'r', encoding='utf8') as f:
    for line in f:
        if 'DT=(SINGLE SINGLE SINGLE SINGLE )' in line:
            # Receiver1
            pressure = f.readline()
            pressure1 = pressure.split('\t')
            pressure1 = float(pressure1[4])
            pressure1_arrayP2 = np.append(pressure1_arrayP2, pressure1)
            # Receiver2
            pressure = f.readline()
            pressure2 = pressure.split('\t')
            pressure2 = float(pressure2[4])
            pressure2_arrayP2 = np.append(pressure2_arrayP2, pressure2)
            # Receiver3
            pressure = f.readline()
            pressure3 = pressure.split('\t')
            pressure3 = float(pressure3[4])
            pressure3_arrayP2 = np.append(pressure3_arrayP2, pressure3)
            # Receiver4
            pressure = f.readline()
            pressure4 = pressure.split('\t')
            pressure4 = float(pressure4[4])
            pressure4_arrayP2 = np.append(pressure4_arrayP2, pressure4)
            # Receiver5
            pressure = f.readline()
            pressure5 = pressure.split('\t')
            pressure5 = float(pressure5[4])
            pressure5_arrayP2 = np.append(pressure5_arrayP2, pressure5)
f.close()
# End Python Reader 2

# Python Reader 3 

time_arrayP3 = np.array([])
pressure1_arrayP3 = np.array([])
pressure2_arrayP3 = np.array([])
pressure3_arrayP3 = np.array([])
pressure4_arrayP3 = np.array([])
pressure5_arrayP3 = np.array([])

with open('StratPython_comparison_0.3.txt', 'r', encoding='utf8') as f:
    for line in f:
        if 'ZONE T=" Single Point "' in line:
            time = f.readline()
            time = time.split(' ')
            time = float(time[2])
            time_arrayP3 = np.append(time_arrayP3, time)
f.close()

with open('StratPython_comparison_0.3.txt', 'r', encoding='utf8') as f:
    for line in f:
        if 'DT=(SINGLE SINGLE SINGLE SINGLE )' in line:
            # Receiver1
            pressure = f.readline()
            pressure1 = pressure.split('\t')
            pressure1 = float(pressure1[4])
            pressure1_arrayP3 = np.append(pressure1_arrayP3, pressure1)
            # Receiver2
            pressure = f.readline()
            pressure2 = pressure.split('\t')
            pressure2 = float(pressure2[4])
            pressure2_arrayP3 = np.append(pressure2_arrayP3, pressure2)
            # Receiver3
            pressure = f.readline()
            pressure3 = pressure.split('\t')
            pressure3 = float(pressure3[4])
            pressure3_arrayP3 = np.append(pressure3_arrayP3, pressure3)
            # Receiver4
            pressure = f.readline()
            pressure4 = pressure.split('\t')
            pressure4 = float(pressure4[4])
            pressure4_arrayP3 = np.append(pressure4_arrayP3, pressure4)
            # Receiver5
            pressure = f.readline()
            pressure5 = pressure.split('\t')
            pressure5 = float(pressure5[4])
            pressure5_arrayP3 = np.append(pressure5_arrayP3, pressure5)
f.close()
#End Python Reader 3 

#Python Reader 4

time_arrayP4 = np.array([])
pressure1_arrayP4 = np.array([])
pressure2_arrayP4 = np.array([])
pressure3_arrayP4 = np.array([])
pressure4_arrayP4 = np.array([])
pressure5_arrayP4 = np.array([])

with open('StratPython_comparison_0.6.txt', 'r', encoding='utf8') as f:
    for line in f:
        if 'ZONE T=" Single Point "' in line:
            time = f.readline()
            time = time.split(' ')
            time = float(time[2])
            time_arrayP4 = np.append(time_arrayP4, time)
f.close()

with open('StratPython_comparison_0.6.txt', 'r', encoding='utf8') as f:
    for line in f:
        if 'DT=(SINGLE SINGLE SINGLE SINGLE )' in line:
            # Receiver1
            pressure = f.readline()
            pressure1 = pressure.split('\t')
            pressure1 = float(pressure1[4])
            pressure1_arrayP4 = np.append(pressure1_arrayP4, pressure1)
            # Receiver2
            pressure = f.readline()
            pressure2 = pressure.split('\t')
            pressure2 = float(pressure2[4])
            pressure2_arrayP4 = np.append(pressure2_arrayP4, pressure2)
            # Receiver3
            pressure = f.readline()
            pressure3 = pressure.split('\t')
            pressure3 = float(pressure3[4])
            pressure3_arrayP4 = np.append(pressure3_arrayP4, pressure3)
            # Receiver4
            pressure = f.readline()
            pressure4 = pressure.split('\t')
            pressure4 = float(pressure4[4])
            pressure4_arrayP4 = np.append(pressure4_arrayP4, pressure4)
            # Receiver5
            pressure = f.readline()
            pressure5 = pressure.split('\t')
            pressure5 = float(pressure5[4])
            pressure5_arrayP4 = np.append(pressure5_arrayP4, pressure5)
f.close()
#End Python reader 4 

#Python Plots 

# Plotting Receiver1
plt.figure(layout='constrained')
plt.plot(time_arrayP, pressure1_arrayP, color='red', alpha=0.5)
plt.plot(time_arrayP2, pressure1_arrayP2, color='blue', alpha=0.5)
plt.plot(time_arrayP3, pressure1_arrayP3, color = 'green', alpha = 0.5)
plt.plot(time_arrayP4, pressure1_arrayP4, color = 'black', alpha = 0.5)
plt.title("Pressure vs Time of Receiver 1 Pythons")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#Plotting Receiver 2 

plt.figure(layout='constrained')
plt.plot(time_arrayP, pressure2_arrayP, color='red', alpha=0.5)
plt.plot(time_arrayP2, pressure2_arrayP2, color='blue', alpha=0.5)
plt.plot(time_arrayP3, pressure2_arrayP3, color = 'green', alpha = 0.5)
plt.plot(time_arrayP4, pressure2_arrayP4, color = 'black', alpha = 0.5)
plt.title("Pressure vs Time of Receiver 2 Pythons")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#Plotting R3
plt.figure(layout='constrained')
plt.plot(time_arrayP, pressure3_arrayP, color='red', alpha=0.5)
plt.plot(time_arrayP2, pressure3_arrayP2, color='blue', alpha=0.5)
plt.plot(time_arrayP3, pressure3_arrayP3, color = 'green', alpha = 0.5)
plt.plot(time_arrayP4, pressure3_arrayP4, color = 'black', alpha = 0.5)
plt.title("Pressure vs Time of Receiver 3 Pythons")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#Plotting R4
plt.figure(layout='constrained')
plt.plot(time_arrayP, pressure4_arrayP, color='red', alpha=0.5)
plt.plot(time_arrayP2, pressure4_arrayP2, color='blue', alpha=0.5)
plt.plot(time_arrayP3, pressure4_arrayP3, color = 'green', alpha = 0.5)
plt.plot(time_arrayP4, pressure4_arrayP4, color = 'black', alpha = 0.5)
plt.title("Pressure vs Time of Receiver 4 Pythons")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#Plotting R5
plt.figure(layout='constrained')
plt.plot(time_arrayP, pressure5_arrayP, color='red', alpha=0.5)
plt.plot(time_arrayP2, pressure5_arrayP2, color='blue', alpha=0.5)
plt.plot(time_arrayP3, pressure5_arrayP3, color = 'green', alpha = 0.5)
plt.plot(time_arrayP4, pressure5_arrayP4, color = 'black', alpha = 0.5)
plt.title("Pressure vs Time of Receiver 5 Pythons")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()



# Fortran Reader Code

# Creates the neccesary blank arrays for the program.

time_arrayF = np.array([])
pressure1_arrayF = np.array([])
pressure2_arrayF = np.array([])
pressure3_arrayF = np.array([])
pressure4_arrayF = np.array([])
pressure5_arrayF = np.array([])

# Mirrors the code that appends values to the time array, but is tailored for the different spacing and readout of the
# Fortran text file as opposed to the Python text file.

with open('OrigFortran_comparison_006.dat', 'r', encoding='utf8') as f:
    for line in f:
        if 'ZONE T="Single Point        "' in line:
            time = f.readline()
            time = time.split('  ')
            time = float(time[1])
            time_arrayF = np.append(time_arrayF, time)
f.close()

# Same thing as above, but for the pressure receivers.

with open('OrigFortran_comparison_006.dat', 'r', encoding='utf8') as f:
    for line in f:
        if 'DT=(SINGLE SINGLE SINGLE SINGLE )' in line:
            # Receiver 1
            pressure = f.readline()
            pressure1 = pressure.split('     ')
            pressure1 = float(pressure1[3])
            pressure1_arrayF = np.append(pressure1_arrayF, pressure1)
            # Receiver2
            pressure = f.readline()
            pressure2 = pressure.split('     ')
            pressure2 = float(pressure2[3])
            pressure2_arrayF = np.append(pressure2_arrayF, pressure2)
            # Receiver3
            pressure = f.readline()
            pressure3 = pressure.split('     ')
            pressure3 = float(pressure3[3])
            pressure3_arrayF = np.append(pressure3_arrayF, pressure3)
            # Receiver4
            pressure = f.readline()
            pressure4 = pressure.split('     ')
            pressure4 = float(pressure4[3])
            pressure4_arrayF = np.append(pressure4_arrayF, pressure4)
            # Receiver5
            pressure = f.readline()
            pressure5 = pressure.split('     ')
            pressure5 = float(pressure5[3])
            pressure5_arrayF = np.append(pressure5_arrayF, pressure5)
f.close()

# End of Fortran Reader Code 1

#Fortran Reader 2

time_arrayF2 = np.array([])
pressure1_arrayF2 = np.array([])
pressure2_arrayF2 = np.array([])
pressure3_arrayF2 = np.array([])
pressure4_arrayF2 = np.array([])
pressure5_arrayF2 = np.array([])

# Mirrors the code that appends values to the time array, but is tailored for the different spacing and readout of the
# Fortran text file as opposed to the Python text file.

with open('OrigFortran_comparison_03.dat', 'r', encoding='utf8') as f:
    for line in f:
        if 'ZONE T="Single Point        "' in line:
            time = f.readline()
            time = time.split('  ')
            time = float(time[1])
            time_arrayF2 = np.append(time_arrayF2, time)
f.close()

# Same thing as above, but for the pressure receivers.

with open('OrigFortran_comparison_03.dat', 'r', encoding='utf8') as f:
    for line in f:
        if 'DT=(SINGLE SINGLE SINGLE SINGLE )' in line:
            # Receiver 1
            pressure = f.readline()
            pressure1 = pressure.split('     ')
            pressure1 = float(pressure1[3])
            pressure1_arrayF2 = np.append(pressure1_arrayF2, pressure1)
            # Receiver2
            pressure = f.readline()
            pressure2 = pressure.split('     ')
            pressure2 = float(pressure2[3])
            pressure2_arrayF2 = np.append(pressure2_arrayF2, pressure2)
            # Receiver3
            pressure = f.readline()
            pressure3 = pressure.split('     ')
            pressure3 = float(pressure3[3])
            pressure3_arrayF2 = np.append(pressure3_arrayF2, pressure3)
            # Receiver4
            pressure = f.readline()
            pressure4 = pressure.split('     ')
            pressure4 = float(pressure4[3])
            pressure4_arrayF2 = np.append(pressure4_arrayF2, pressure4)
            # Receiver5
            pressure = f.readline()
            pressure5 = pressure.split('     ')
            pressure5 = float(pressure5[3])
            pressure5_arrayF2 = np.append(pressure5_arrayF2, pressure5)
f.close()
#End Fortran Reader 2

#Fortran Reader 3 
time_arrayF3 = np.array([])
pressure1_arrayF3 = np.array([])
pressure2_arrayF3 = np.array([])
pressure3_arrayF3 = np.array([])
pressure4_arrayF3 = np.array([])
pressure5_arrayF3 = np.array([])

# Mirrors the code that appends values to the time array, but is tailored for the different spacing and readout of the
# Fortran text file as opposed to the Python text file.

with open('OrigFortran_comparison_06.dat', 'r', encoding='utf8') as f:
    for line in f:
        if 'ZONE T="Single Point        "' in line:
            time = f.readline()
            time = time.split('  ')
            time = float(time[1])
            time_arrayF3 = np.append(time_arrayF3, time)
f.close()

# Same thing as above, but for the pressure receivers.

with open('OrigFortran_comparison_06.dat', 'r', encoding='utf8') as f:
    for line in f:
        if 'DT=(SINGLE SINGLE SINGLE SINGLE )' in line:
            # Receiver 1
            pressure = f.readline()
            pressure1 = pressure.split('     ')
            pressure1 = float(pressure1[3])
            pressure1_arrayF3 = np.append(pressure1_arrayF3, pressure1)
            # Receiver2
            pressure = f.readline()
            pressure2 = pressure.split('     ')
            pressure2 = float(pressure2[3])
            pressure2_arrayF3 = np.append(pressure2_arrayF3, pressure2)
            # Receiver3
            pressure = f.readline()
            pressure3 = pressure.split('     ')
            pressure3 = float(pressure3[3])
            pressure3_arrayF3 = np.append(pressure3_arrayF3, pressure3)
            # Receiver4
            pressure = f.readline()
            pressure4 = pressure.split('     ')
            pressure4 = float(pressure4[3])
            pressure4_arrayF3 = np.append(pressure4_arrayF3, pressure4)
            # Receiver5
            pressure = f.readline()
            pressure5 = pressure.split('     ')
            pressure5 = float(pressure5[3])
            pressure5_arrayF3 = np.append(pressure5_arrayF3, pressure5)
f.close()
#End Reader 3 

#Reader 4 Fortran
time_arrayF4 = np.array([])
pressure1_arrayF4 = np.array([])
pressure2_arrayF4 = np.array([])
pressure3_arrayF4 = np.array([])
pressure4_arrayF4 = np.array([])
pressure5_arrayF4 = np.array([])

# Mirrors the code that appends values to the time array, but is tailored for the different spacing and readout of the
# Fortran text file as opposed to the Python text file.

with open('OrigFortran_comparison_01.dat', 'r', encoding='utf8') as f:
    for line in f:
        if 'ZONE T="Single Point        "' in line:
            time = f.readline()
            time = time.split('  ')
            time = float(time[1])
            time_arrayF4 = np.append(time_arrayF4, time)
f.close()

# Same thing as above, but for the pressure receivers.

with open('OrigFortran_comparison_01.dat', 'r', encoding='utf8') as f:
    for line in f:
        if 'DT=(SINGLE SINGLE SINGLE SINGLE )' in line:
            # Receiver 1
            pressure = f.readline()
            pressure1 = pressure.split('     ')
            pressure1 = float(pressure1[3])
            pressure1_arrayF4 = np.append(pressure1_arrayF4, pressure1)
            # Receiver2
            pressure = f.readline()
            pressure2 = pressure.split('     ')
            pressure2 = float(pressure2[3])
            pressure2_arrayF4 = np.append(pressure2_arrayF4, pressure2)
            # Receiver3
            pressure = f.readline()
            pressure3 = pressure.split('     ')
            pressure3 = float(pressure3[3])
            pressure3_arrayF4 = np.append(pressure3_arrayF4, pressure3)
            # Receiver4
            pressure = f.readline()
            pressure4 = pressure.split('     ')
            pressure4 = float(pressure4[3])
            pressure4_arrayF4 = np.append(pressure4_arrayF4, pressure4)
            # Receiver5
            pressure = f.readline()
            pressure5 = pressure.split('     ')
            pressure5 = float(pressure5[3])
            pressure5_arrayF4 = np.append(pressure5_arrayF4, pressure5)
f.close()
# #End Reader 4 

#Plots for Fortran

# Plotting Receiver1
plt.figure(layout='constrained')
plt.plot(time_arrayF, pressure1_arrayF, color='red', alpha=0.5)
plt.plot(time_arrayF2, pressure1_arrayF2, color='blue', alpha=0.5)
plt.plot(time_arrayF3, pressure1_arrayF3, color = 'green', alpha = 0.5)
plt.plot(time_arrayF4, pressure1_arrayF4, color = 'black', alpha = 0.5)
plt.title("Pressure vs Time of Receiver 1 Fortran")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#Plotting Receiver 2 

plt.figure(layout='constrained')
plt.plot(time_arrayF, pressure2_arrayF, color='red', alpha=0.5)
plt.plot(time_arrayF2, pressure2_arrayF2, color='blue', alpha=0.5)
plt.plot(time_arrayF3, pressure2_arrayF3, color = 'green', alpha = 0.5)
plt.plot(time_arrayF4, pressure2_arrayF4, color = 'black', alpha = 0.5)
plt.title("Pressure vs Time of Receiver 2 Fortran")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#Plotting R3
plt.figure(layout='constrained')
plt.plot(time_arrayF, pressure3_arrayF, color='red', alpha=0.5)
plt.plot(time_arrayF2, pressure3_arrayF2, color='blue', alpha=0.5)
plt.plot(time_arrayF3, pressure3_arrayF3, color = 'green', alpha = 0.5)
plt.plot(time_arrayF4, pressure3_arrayF4, color = 'black', alpha = 0.5)
plt.title("Pressure vs Time of Receiver 3 Fortran")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#Plotting R4
plt.figure(layout='constrained')
plt.plot(time_arrayF, pressure4_arrayF, color='red', alpha=0.5)
plt.plot(time_arrayF2, pressure4_arrayF2, color='blue', alpha=0.5)
plt.plot(time_arrayF3, pressure4_arrayF3, color = 'green', alpha = 0.5)
plt.plot(time_arrayF4, pressure4_arrayF4, color = 'black', alpha = 0.5)
plt.title("Pressure vs Time of Receiver 4 Fortran")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#Plotting R5
plt.figure(layout='constrained')
plt.plot(time_arrayF, pressure5_arrayF, color='red', alpha=0.5)
plt.plot(time_arrayF2, pressure5_arrayF2, color='blue', alpha=0.5)
plt.plot(time_arrayF3, pressure5_arrayF3, color = 'green', alpha = 0.5)
plt.plot(time_arrayF4, pressure5_arrayF4, color = 'black', alpha = 0.5)
plt.title("Pressure vs Time of Receiver 5 Fortran")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()
#End of Fortran Plots


#Fortran vs Python Plots

# Plotting Fortran vs Python R1's

#0.06
plt.figure(layout='constrained')
plt.plot(time_arrayF, pressure1_arrayF, color='blue', alpha=0.5)
plt.plot(time_arrayP, pressure1_arrayP, color='red', alpha=0.5)
plt.title("Pressure vs Time of Receiver 1 Fortran v Python 0.03")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#0.1
plt.figure(layout='constrained')
plt.plot(time_arrayF4, pressure1_arrayF4, color='blue', alpha=0.5)
plt.plot(time_arrayP2, pressure1_arrayP2, color='red', alpha=0.5)
plt.title("Pressure vs Time of Receiver 1 Fortran v Python 0.1")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#0.3
plt.figure(layout='constrained')
plt.plot(time_arrayF2, pressure1_arrayF2, color='blue', alpha=0.5)
plt.plot(time_arrayP3, pressure1_arrayP3, color = 'red', alpha = 0.5)
plt.title("Pressure vs Time of Receiver 1 Fortran vs Python 0.3")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#0.6
plt.figure(layout='constrained')
plt.plot(time_arrayF3, pressure1_arrayF3, color='blue', alpha=0.5)
plt.plot(time_arrayP4, pressure1_arrayP4, color='red', alpha=0.5)
plt.title("Pressure vs Time of Receiver 1 Fortran vs Python 0.6")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()


#Plotting Receiver 2 Python vs Fortran

#0.06
plt.figure(layout='constrained')
plt.plot(time_arrayF, pressure2_arrayF, color='blue', alpha=0.5)
plt.plot(time_arrayP, pressure2_arrayP, color='red', alpha=0.5)
plt.title("Pressure vs Time of Receiver 2 Fortran vs Python 0.03")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#0.1
plt.figure(layout='constrained')
plt.plot(time_arrayF4, pressure2_arrayF4, color='blue', alpha=0.5)
plt.plot(time_arrayP2, pressure2_arrayP2, color='red', alpha=0.5)
plt.title("Pressure vs Time of Receiver 2 Fortran v Python 0.1")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#0.3
plt.figure(layout='constrained')
plt.plot(time_arrayF2, pressure2_arrayF2, color='blue', alpha=0.5)
plt.plot(time_arrayP3, pressure2_arrayP3, color = 'red', alpha = 0.5)
plt.title("Pressure vs Time of Receiver 2 Fortran v Python 0.3")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#0.6
plt.figure(layout='constrained')
plt.plot(time_arrayF3, pressure2_arrayF3, color='blue', alpha=0.5)
plt.plot(time_arrayP4, pressure2_arrayP4, color = 'red', alpha = 0.5)
plt.title("Pressure vs Time of Receiver 2 Fortran v Python 0.6")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#Plotting Receiver 3 F v P 

# 0.06
plt.figure(layout='constrained')
plt.plot(time_arrayF, pressure3_arrayF, color='blue', alpha=0.5)
plt.plot(time_arrayP, pressure3_arrayP, color='red', alpha=0.5)
plt.title("Pressure vs Time of Receiver 3 Fortran v Python 0.03")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#0.1
plt.figure(layout='constrained')
plt.plot(time_arrayF4, pressure3_arrayF4, color='blue', alpha=0.5)
plt.plot(time_arrayP2, pressure3_arrayP2, color='red', alpha=0.5)
plt.title("Pressure vs Time of Receiver 3 Fortran v Python 0.1")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#0.3
plt.figure(layout='constrained')
plt.plot(time_arrayF2, pressure3_arrayF2, color='blue', alpha=0.5)
plt.plot(time_arrayP3, pressure3_arrayP3, color='red', alpha=0.5)
plt.title("Pressure vs Time of Receiver 3 Fortran v Python 0.3")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#0.6
plt.figure(layout='constrained')
plt.plot(time_arrayF3, pressure3_arrayF3, color='blue', alpha=0.5)
plt.plot(time_arrayP4, pressure3_arrayP4, color='red', alpha=0.5)
plt.title("Pressure vs Time of Receiver 3 Fortran v Python 0.6")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#Plotting Receiver 5 Python v Fortran

#Plotting R5

#0.06
plt.figure(layout='constrained')
plt.plot(time_arrayF, pressure5_arrayF, color='blue', alpha=0.5)
plt.plot(time_arrayP, pressure5_arrayP, color='red', alpha=0.5)
plt.title("Pressure vs Time of Receiver 5 Fortran v Python 0.03")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#0.1
plt.figure(layout='constrained')
plt.plot(time_arrayF4, pressure5_arrayF4, color='blue', alpha=0.5)
plt.plot(time_arrayP2, pressure5_arrayP2, color='red', alpha=0.5)
plt.title("Pressure vs Time of Receiver 5 Fortran v Python 0.1")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#0.3
plt.figure(layout='constrained')
plt.plot(time_arrayF2, pressure5_arrayF2, color='blue', alpha=0.5)
plt.plot(time_arrayP3, pressure5_arrayP3, color='red', alpha=0.5)
plt.title("Pressure vs Time of Receiver 5 Fortran v Python 0.3")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()

#0.6
plt.figure(layout='constrained')
plt.plot(time_arrayF3, pressure5_arrayF3, color='blue', alpha=0.5)
plt.plot(time_arrayP4, pressure5_arrayP4, color='red', alpha=0.5)
plt.title("Pressure vs Time of Receiver 5 Fortran v Python 0.6")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.xlim(0.0, 0.8)
plt.grid()
plt.show()







# # Using MatPlotlib to graph both sets of data against each other

# # Plotting Receiver1
# plt.figure(layout='constrained')
# plt.plot(time_arrayP, pressure1_arrayP, color='red', alpha=0.5)
# plt.plot(time_arrayF, pressure1_arrayF, color='blue', alpha=0.5)
# #plt.plot(time_arrayP2, pressure1_arrayP2, color = 'green', alpha = 0.5)
# plt.title("Pressure vs Time of Receiver 1")
# plt.xlabel("Time [s]")
# plt.ylabel("Pressure [Pa]")
# plt.xlim(0.0, 0.8)
# plt.grid()
# plt.show()

# # Plotting Receiver2
# plt.figure(layout='constrained')
# plt.plot(time_arrayP, pressure2_arrayP, color='red', alpha=0.5)
# plt.plot(time_arrayF, pressure2_arrayF, color='blue', alpha=0.5)
# #plt.plot(time_arrayP2, pressure2_arrayP2, color = 'green', alpha = 0.5)
# plt.title("Pressure vs Time of Receiver 2")
# plt.xlabel("Time [s]")
# plt.ylabel("Pressure [Pa]")
# plt.xlim(0.0, 0.8)
# plt.grid()
# plt.show()

# # Plotting Receiver 3
# plt.figure(layout='constrained')
# plt.plot(time_arrayP, pressure3_arrayP, color='red', alpha=0.5)
# plt.plot(time_arrayF, pressure3_arrayF, color='blue', alpha=0.5)
# #plt.plot(time_arrayP2, pressure3_arrayP2, color = 'green', alpha = 0.25)
# plt.title("Pressure vs Time of Receiver 3")
# plt.xlabel("Time [s]")
# plt.ylabel("Pressure [Pa]")
# plt.xlim(0.0, 0.8)
# plt.grid()
# plt.show()

# # Plotting Receiver 4
# plt.figure(layout='constrained')
# plt.plot(time_arrayP, pressure4_arrayP, color='red', alpha=0.5)
# plt.plot(time_arrayF, pressure4_arrayF, color='blue', alpha=0.5)
# #plt.plot(time_arrayP2, pressure4_arrayP2, color = 'green', alpha = 0.25)
# plt.title("Pressure vs Time of Receiver 4")
# plt.xlabel("Time [s]")
# plt.ylabel("Pressure [Pa]")
# plt.xlim(0.0, 0.8)
# plt.grid()
# plt.show()

# # Plotting Receiver 5
# plt.figure(layout='constrained')
# plt.plot(time_arrayP, pressure5_arrayP, color='red', alpha=0.5)
# plt.plot(time_arrayF, pressure5_arrayF, color='blue', alpha=0.5)
# #plt.plot(time_arrayP2, pressure5_arrayP2, color = 'green', alpha = 0.25)
# plt.title("Pressure vs Time of Receiver 5")
# plt.xlabel("Time [s]")
# plt.ylabel("Pressure [Pa]")
# plt.xlim(0.0, 0.8)
# plt.grid()
# plt.show()

