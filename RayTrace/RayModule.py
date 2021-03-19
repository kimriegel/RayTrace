# Smarter than the average ray

import numpy as np
import config 

class RayModule:
    """
    To increase functionality and available perfomance increases
    """
    #Default values
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
        Ray Module
        Create and defines position of ray        
        """

        RayModule.rayList.append(self)  #debug

        self.position = np.array(position)
        self.frequency = np.array(RayModule.frequency )
        self.direction = np.array(RayModule.direction )
        self.phase     = np.array(RayModule.phase     )
        self.amplitude = np.array(RayModule.amplitude )

        self.step = (foo(i) for i in range(RayModule.imax))

    def reset(self, position ):
        """
        I have not mastered classes yet
        The goal of this function is to create a new ray at the end of 
         the previous ray's life while not making multiple objects

        Delete when this is figured out

        """

        # Can use list to make sure there is still one object but multiple names
        RayModule.rayList.append(self)  #debug

        self.position =  np.array(position)
        self.frequency = np.array(RayModule.frequency )
        self.direction = np.array(RayModule.direction )
        self.phase     = np.array(RayModule.phase     )
        self.amplitude = np.array(RayModule.amplitude )

        #self.step = (foo(i) for i in range(Ray.imax))
        
    def update_freq(self,dx_update, alpha_update, diffusion_update, lamb, air_absorb):
        """
        Update ray phase and amplitude
        """
        two_pi_dx_update = config.twopi * dx_update
        ein = self.phase - (two_pi_dx_update / lamb)
        zwei = ein % config.twopi
        masque = zwei > np.pi
        drei = masque * zwei - config.twopi 

        self.phase = np.where(masque, drei, ein)
        self.amplitude *= ((1.0 - alpha_update) * (1.0 - diffusion_update) * np.exp(-air_absorb * dx_update))


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

#r_1 = RayModule((93.4428213, 28.8397178, 0.151))
#r_2 = RayModule((64.5832, 28.5998, 0.151))

#I=3
#for k in range(3):
#    step= (i*i for i in range(I))   #aka rayPath
#    #print(next(step))
#    for j in step:
#        print(j)

## Sucessfully uses only one ray
#print(RayModule.rayList)
#Jeff = RayModule((1,1,1))
#print(RayModule.rayList)
#Jeff = RayModule((1,1,1))
#print(RayModule.rayList)
#Jeff.reset((1,1,1))
#print(RayModule.rayList)
#
#x=Jeff.position
#print(x is Jeff.position)
#print(Jeff.position)
#Jeff.position += 1
#print(Jeff.position)
#print(x)
#
##