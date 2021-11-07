# This file defines atmosphere class with several options:
#
#            Constant
#            Standard Atmosphere
#            File input
import Parameterfile as Pf
import numpy as np


class Atmosphere:

    def __init__(self, ground_temp, strat_height, type):
        if type == 1:
            self.strata = np.linspace(0, Pf.zmax, Pf.zmax/strat_height+1)
            self.sound_speed = np.ones(len(self.strata))*331.3 + 0.606*(ground_temp-273)
        elif type == 2:
            # This is a place holder
            self.strata = np.linspace(0, Pf.zmax, Pf.zmax / strat_height)
        elif type == 3:
            # This is a place holder
            self.strata = np.linspace(0, Pf.zmax, Pf.zmax / strat_height)

