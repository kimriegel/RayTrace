import matplotlib.pyplot as plt  # for graphing
import numpy as np
import math as m


# For N wave
time_array=np.linspace(1,1000)*m.pi/180.
pressure = np.sin(time_array)
pressure2= np.cos(time_array)

plt.figure(num=1, figsize=(5.20, 5.80), dpi=120, facecolor='#e6e6fa', edgecolor='g')  # lavender
plt.grid(True)
plt.plot(time_array, pressure, 'm',time_array,pressure2,'b')
    # Labeling axes
plt.xlabel('Time [s]')
plt.ylabel('Pressure [Pa]')
plt.title('Pressure vs Time of Receiver ',
            fontsize=10,
            fontweight='bold')

    # Saving
    # plt.savefig(Pf.graphName + str(i) + '.png', facecolor='#eeeeee')    # grey
    # plt.savefig(Pf.graphName + str(i) + '.png', facecolor='#e0dae6')    # muted lilac
plt.show()
print('Graph time: ', time.time() - t)