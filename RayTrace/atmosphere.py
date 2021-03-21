"""
This will contain any fututre atmospheric conditions such as:
Temperature
Turbulence

By default this will:
Stratify the environment to decrease comp time
Destroy steps as they currently exist
Look Pretty
"""

class Receiver:
    """
    Receiver class with attributes of position (x,y,z)
    Functionality for pressure not yet added.
    """
    rList = []  # See append_list

    def __init__(self, position):
        """
        Create and defines position of receiver
        
        Works automatically when class is called
        """
        self.position = np.array(position)
        self.recNumber = int(Receiver.arraysize)         # planned for debugging but we don't seem to use it
        Receiver.rList.append(self)  # See append_list
        # Initial values
        self.pressure = 0
        self.magnitude = 0
        self.direction = 0



import numpy as np


strataMod = np.arange(1,10 +1,0.5)          # The modfier to be multiplied by temperature
temp = 297                              # In Kelvin
print(strataMod)
print((strataMod * temp))