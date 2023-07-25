import numpy as np  # matrices and arrays
import matplotlib.pyplot as plt  # for graphing


with open("Putty.txt") as f:
    contents = f.read()
    count = contents.count("xcoor")

print(count)

#putty = open("Putty.txt",'r')


#with open('Putty.txt', 'r', encoding='utf8') as f:
    #for line in f:
        
