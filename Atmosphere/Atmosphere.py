# This file defines atmosphere class with several options:
#
#            Constant
#            Standard Atmosphere
#            File input
import Parameterfile as Pf
import numpy as np
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
            # This is a place holder
            # temp with geopotential heights, first strata a little above building, estimate temperature to ground
            # figure linear point difference between points y=mx+b etc,
            # if lowest height, extrapolate what 0 is.
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
            geopotential_array = geopotential[0, :, 0,0]

            # assigns each column from the proper index into its own array.

            tempIntegersArray = np.array(temperature_array)
            windIntegersArray = np.array(wind_array)
            vspeedIntegersArray = np.array(vspeed_array)
            geopotentialIntegersArray = np.array(geopotential_array)

            print("Temp integers array before insertion: ", tempIntegersArray)

            #inserts the user input ground temperature in order to gauge a slope and compare against standard.

            geopotentialIntegersArray = np.insert(geopotentialIntegersArray, 0, 0)

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
                # geopotentialGradientArray.append(diffGeoPot)
            print("Size and values Temp Gradient array in : ", len(tempGradientArray), tempGradientArray)
            # print("Size and Values geopotential gradient array: ", len(geopotentialGradientArray))

            tempSlope = []


            # Calculates change in geopotential array up to 10,207.7ft

            for j in range(1, len(geopotentialIntegersArray) - 1):
                diffGeoPot = geopotentialIntegersArray[j] - geopotentialIntegersArray[j-1]
                geopotentialGradientArray.append(diffGeoPot)
            print("Length of geopotential gradient array, and geopotarray in : ", len(geopotentialGradientArray), geopotentialGradientArray)
                # slope = geopotentialGradientArray[j] / tempGradientArray[j]

            # calculates slope divided by 1000 meters. should have 0.006

            for i in range(1, len(tempGradientArray)):
                slope = (geopotentialGradientArray[i] / tempGradientArray[i]) / 1000
                tempSlope.append(slope)

            print("Rate of change in : ", tempSlope)
            print("Slope from ground to first geopotential height: ", tempSlope[0])
                # tempGeoSlope = geopotentialGradientArray[i] / tempGradientArray[i]

            # calculates ground_temp and assigns it to ground_temp

            ground_temp = tempIntegersArray[0] + (geopotentialIntegersArray[1] * tempSlope[0])/1000

            print("Ground temp is approximately: ", ground_temp)

            #
            # print("tempIntegersArray: ", len(tempIntegersArray))
            # print("tempgradientarray: ", len(tempGradientArray))
            # print("geopotentialgradarray length: ", len(geopotentialGradientArray))
            # print("temp height slope: ", len(tempSlope))

            # slope = np.mean(tempGeoSlope)

            self.strata = np.linspace(0, geopotentialIntegersArray[-1],int((geopotentialIntegersArray[-1] / strat_height) + 1))
            self.sound_speed = np.ones(len(self.strata)) * 331.3 + 0.606*(ground_temp-273.15)
            # self.sound_speed =

           # print("Temperature slope: ", slope)

            # need to figure out a way to take the above, and fit them for parameter file in order to work the array
            # behavior over the building domains.

            # self.strata = np.linspace(0, Pf.zmax, Pf.zmax / strat_height)

            # Find out what strat height in paramter file means, and what Temp in parameter file means.
            # is it a height interval for it to iterate the ray over or?
            # is it the temperature at the lowest ground height etc?

# strat_height is
# temp is