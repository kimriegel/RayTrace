# Smarter than the average ray

import numpy as np

class Ray:
    """
    To increase functionality and available perfomance increases
    """
    frequency = None
    direction = None
    phase = None
    amplitude = None
    imax = 75

    rayList = [] # for debugging, to make sure rays are deleted after use

    # frequency, direction, phase, amplitude are independent to each ray, but start the same
    # position is diffent at all points

    def __init__(self, position ):
        """
        Create and defines position of ray        
        """

        Ray.rayList.append(self)  # See append_list

        self.position = np.array(position)
        self.frequency = np.array(Ray.frequency )
        self.direction = np.array(Ray.direction )
        self.phase     = np.array(Ray.phase     )
        self.amplitude = np.array(Ray.amplitude )

        self.step = (foo(i) for i in range(Ray.imax))

    def reset(self, position ):
        """
        I have not mastered classes yet
        The goal of this function is to create a new ray at the end of 
         the previous ray's life while not making multiple objects

        Delete when this is figured out

        """
        Ray.rayList.append(self)  # See append_list

        self.position = np.array(position)
        self.frequency = np.array(Ray.frequency )
        self.direction = np.array(Ray.direction )
        self.phase     = np.array(Ray.phase     )
        self.amplitude = np.array(Ray.amplitude )

        self.step = (foo(i) for i in range(Ray.imax))


    def foo(self,data):
        """
        Checks the nearest surface and updates ray data based on results
        """
        pass

    @classmethod
    def initialize(cls, frequency, direction, phase, amplitude):
        """
        Add values of the values that start all rays as constants
        """

        cls.frequency = frequency
        cls.direction = direction
        cls.phase     = phase
        cls.amplitude = amplitude
       


# Initializing Rays

#r_1 = Ray((93.4428213, 28.8397178, 0.151))
#r_2 = Ray((64.5832, 28.5998, 0.151))

#I=3
#for k in range(3):
#    step= (i*i for i in range(I))   #aka rayPath
#    #print(next(step))
#    for j in step:
#        print(j)

print(Ray.rayList)
Jeff = Ray((1,1,1))
print(Ray.rayList)
Jeff = Ray((1,1,1))
print(Ray.rayList)
Jeff.reset((1,1,1))
print(Ray.rayList)