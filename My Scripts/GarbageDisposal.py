import numpy as np

f = open('Constant_Comp_Simple0.06.txt')
count = 0

garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()

time_array = np.array([])

print(garbage[26:34])

garbage2 = float(garbage[26:34])

time_array = np.append(time_array, garbage2)

print(time_array)
print(garbage2)

garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()

print(garbage)

garbage2 = float(garbage[26:34])

time_array = np.append(time_array, garbage2)

print(time_array)
print(garbage2)

garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()
garbage = f.readline()

print(garbage)

garbage2 = float(garbage[26:34])

time_array = np.append(time_array, garbage2)

print(time_array)
print(garbage2)