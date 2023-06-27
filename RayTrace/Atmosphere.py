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
            # this is a constant temperature atmosphere
            self.strata = np.linspace(0.0, Pf.zmax, int(Pf.zmax/strat_height+1))
            self.sound_speed = np.ones(len(self.strata))*331.3 + 0.606*(ground_temp-273.15)
        elif type == 2:
            # This is based on the ISO standard in meters and Celcius. This is only valid up to the tropopause @ 20000m
            self.strata = np.linspace(0, Pf.zmax, int(Pf.zmax / strat_height + 1))
            temp = np.array([])
            for strat in self.strata:
                if strat <= 11000:
                    temp = np.append(temp,ground_temp-6.5*strat/1000)
                else:
                    temp=np.append(temp,ground_temp - 6.5 * 11000 / 1000)
            self.sound_speed = 331.3 + 0.606 * (temp - 273.15)
        elif type == 3:
            # This is a place holder
            self.strata = np.linspace(0, Pf.zmax, Pf.zmax / strat_height)

