import numpy as np

# Goal: make a series of stratified planes along the z axis

normal= (0,0,1)
temperature = np.arange(20,-1,-1)    #dummy values
c= 348                  #also dummy
h=1
height = np.arange(20,-1,-h)    # Guess

#print(temperature)
yeg = height[0]        # brute force

# This algorithm hopefully works well without bugs
"""
We are using a triangle with three infinite points.
Inf,Inf,Z   Inf,-Inf,Z  -Inf,-Inf,Z
I need to check if -Inf,Inf works
"""

