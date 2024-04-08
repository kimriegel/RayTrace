# This file defines atmosphere class with several options:
#
#            Constant
#            Standard Atmosphere
#            File input
import Parameterfile as Pf
import matplotlib.pyplot as plt
import numpy as np
import RayTrace as rt
from netCDF4 import Dataset

class Atmosphere:

    def __init__(self, ground_temp, strat_height, type):
        if type == 1:
            print("This is a debug line for type 1")
            # this is a constant temperature atmosphere
            self.strata = np.linspace(0.0, Pf.zmax, int(Pf.zmax/strat_height+1))
            self.sound_speed = np.ones(len(self.strata))*331.3 + 0.606*(ground_temp-273.15)
        elif type == 2:
            # This is based on the ISO standard in meters and Celcius. This is only valid up to the tropopause @ 20000m
            self.strata = np.linspace(0, Pf.zmax, int(Pf.zmax / strat_height + 1))
            temp = np.array([])
            print("This is a debug line for type 2")
            for strat in self.strata:
                if strat <= 11000:
                    temp = np.append(temp, ground_temp-6.5*strat/1000)
                else:
                    temp = np.append(temp, ground_temp - 6.5 * 11000 / 1000)
            self.sound_speed = 331.3 + 0.606 * (temp - 273.15)
        elif type == 3:
            # temp with geopotential heights, first strata a little above building, estimate temperature to ground
            # figure linear point difference between points y=mx+b etc,
            # generalize linear equation between files.

            # code for user input for specific file requests.
            data_month = input("Enter the month in ALL Capital letters: ")
            data_year = input("Enter the last 2 digits of the year: ")
            data_time = input("Enter the first two digits of the Zulu time (00Z, 06Z, 12Z, 18Z): ")
            FILE_NAME = data_month + data_year + data_time
            print("This is a debug line for type 3")
            # input file directory containing all the nCDF4 files here
            # NOTE: Variable names should probably be changed to class match class inputs.
            # important variables right now are: Geopotential, temp, and wind.
            # proabbly future important variables would append vspeed into account

            # takes user input and parses through input directory to look for proper .nc files

            data_dir = r'RayTrace/AtmosphereProfiles/'
            file_name_stuff=r''+data_dir+FILE_NAME+'.grb2.nc'
            print('file name', file_name_stuff)
            f = Dataset(data_dir+FILE_NAME+'.grb2.nc')
            temp = f.variables['TMP_L100_Avg']
            long = f.variables['lon']
            level0 = f.variables['level0']
            wind = f.variables['U_GRD_L100_Avg']
            vspeed = f.variables['V_GRD_L100_Avg']
            geopotential = f.variables['HGT_L100_Avg']

            # opens up the variable data set and indexes into arrays

            temperature_array = temp[0, :, 0, 0]
            wind_array = wind[0, :, 0, 0]
            vspeed_array = vspeed[0, :, 0, 0]
            geopotential_array = geopotential[0, :, 0, 0]

            # assigns each column from the proper index into its own array.
            # wind and vspeed are for future iterations

            tempIntegersArray = np.array(temperature_array)
            windIntegersArray = np.array(wind_array)
            vspeedIntegersArray = np.array(vspeed_array)
            geopotentialIntegersArray = np.array(geopotential_array)

            print("Temp integers array before insertion: ", tempIntegersArray)

            print("temp at potential height size: ", len(tempIntegersArray), tempIntegersArray)
            print("geo potential height size: ", len(geopotentialIntegersArray), geopotentialIntegersArray)

            tempGradientArray = []
            tempGradientArrayPost = []
            windGradientArray = []
            vspeedGradientArray = []
            geopotentialGradientArray = []
            geopotentialGradientArrayPost = []

            # calculates change in temp array

            print("Temp integers Array Range: ", tempIntegersArray)

            # add null value into temp at first ele

            for i in range(1, len(tempIntegersArray) - 1):
                diffTemp = tempIntegersArray[i + 1] - tempIntegersArray[i]
                tempGradientArray.append(diffTemp)
            print("Size and values Temp Gradient array in : ", len(tempGradientArray), tempGradientArray)

            tempSlope = []

            # Calculates change in geopotential array

            for j in range(1, len(geopotentialIntegersArray) - 1):
                diffGeoPot = geopotentialIntegersArray[j] - geopotentialIntegersArray[j-1]
                geopotentialGradientArray.append(diffGeoPot)
            print("Length of geopotential gradient array, and geopotarray in : ", len(geopotentialGradientArray), geopotentialGradientArray)

            # calculates slope divided by 1000 to give lapse rate per 1000 meters

            for i in range(1, len(tempGradientArray)):
                slope = (tempGradientArray[i]/geopotentialGradientArray[i])
                tempSlope.append(slope)

            print("Rate of change in Temp : ", tempSlope)
            print("Slope from ground to first geopotential height: ", tempSlope[0])

            # calculates ground_temp by assuming slope between first 2 points of both temp and geopot height continues to the ground
            # assigns it to ground_temp

            print("tempint", tempIntegersArray[0], "geoInt", geopotentialIntegersArray[0], "tempSlope", tempSlope[0])
            ground_temp = tempIntegersArray[0] - (geopotentialIntegersArray[0] * tempSlope[0])

            print("Ground temp is approximately: ", ground_temp) #due to python rounding
            np.insert(tempIntegersArray, 0, ground_temp)
            print(tempIntegersArray)
            self.strata = np.linspace(0, geopotentialIntegersArray[-1], int(Pf.zmax/ strat_height + 1))
            print("strat height: ", int(geopotentialIntegersArray[-1] / strat_height) + 1)
            print("Max height: ", geopotentialIntegersArray[-1], "strata: ", self.strata)

            # calculates sound speed and appends to array

            print(tempIntegersArray[0])

            tempCelsius = np.array([])

            for i in range(1, len(tempIntegersArray) + 1):
                tempCelsius = np.append(tempCelsius, tempIntegersArray[i-1] - 273.15)

            print("Temperature in Celsius: ", tempCelsius, "Celsius temp length: ", len(tempCelsius))

            self.sound_speed = 331.3 + 0.606 * tempCelsius

            plt.plot(tempIntegersArray, geopotentialIntegersArray)
            plt.xlabel('Temperature (K)')
            plt.ylabel('Geopotential Height (m)')
            plt.title(data_month + " " + data_year + " " + data_time + " Temperature as a function of Geopotential Height")
            plt.grid(True)
            plt.show()

            plt.plot(tempCelsius, self.sound_speed)
            plt.xlabel('Temperature (C)')
            plt.ylabel('Speed of Sound (m/s)')
            plt.title(data_month + " " + data_year + " " + data_time + " Speed of Sound as a function of Temperature")
            plt.grid(True)
            plt.show()

            plt.plot(self.sound_speed, geopotentialIntegersArray)
            plt.xlabel('Speed of Sound (m/s)')
            plt.ylabel('Geopotential Height (m)')
            plt.title(data_month + " " + data_year + " " + data_time + " Speed of Sound as a function of Geopotential Height")
            plt.grid(True)
            plt.show()


