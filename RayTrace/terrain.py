"""
Terrain class
Includes buildings, ground, any physical env tat isn't atmosphere

Reflection
Replace ray hit (geometry) and box/plane funcs

Going to apply to entire obj file included
Triangles

Find triangles on same plane

1. dxbuilding is going to be this
2.default ground plane is not 0, there were bugs making it zero
    2a. Do not make this zero
"""
import config 
import numpy as np


"""
Action plan

Steps
1. Create the redirect method
2. Identify and initialize faces as terrain in collision check
3. Store terrain in globally accessible storage
4. Check stored terrain for ray interaction before collision check 2
5. Take F of ray and see if line intersects with face defined in terrain 
    if yes, skip check and return dxbuilding
    Will not eliminate need for dxbuilding but should increase speed of czechs

    
"""

class Terrain:

    def __init__(self, position):
        pass

    def dampen(self, absorption):
        """reduces ampliude on hit"""
        pass

    def redirect(self,direction):
        """bounce on hit"""
        pass

