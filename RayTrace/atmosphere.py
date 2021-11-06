import numpy as np
"""
This will contain any fututre atmospheric conditions such as:
Temperature
Turbulence

By default this will:
Stratify the environment to decrease comp time
Destroy steps as they currently exist
Look Pretty
"""

class Strata:
    rList = []  # See append_list


    def __init__(self, max=25, min=0, size=5):
        """
        Create and defines position of receiver
        
        Works automatically when class is called
        """
        #Receiver.rList.append(self)  # See append_list
        
        self.size = size
        self.altitudes = np.arange([min,max,size])


strataMod = np.arange(1,10 +1,0.5)          # The modfier to be multiplied by temperature
temp = 297                              # In Kelvin
print(strataMod)
print((strataMod * temp))