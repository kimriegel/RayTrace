# RayTrace
# version 1.1.0

# Kimberly A. Riegel, PHD created this program to propagate sonic booms around
# large structures, and to graduate. It is a ray tracing model that 
# will include specular and diffuse reflections. It will print out the
# sound field at ear height, at relevant microphone locations, and at 
# the building walls. It will read in the fft of a sonic boom signature.

# Dr. Riegel, William Costa, and George Seaton porting program from Fortran to python

# Initialize variables and functions
import numpy as np  # matrices and arrays
import matplotlib.pyplot as plt  # for graphing

import Parameterfile as Pf
import Functions as Fun
import ReceiverPointSource as Rps  # For receivers
import GeometryParser as Gp


# import GeometryParser as Bg

import time  # Time checks
t = time.time()
phase = 0
amplitude = 0
print(Pf.Fs)


# What it does not do
"""
      Interacts with geometry parser
      Have a way of reading in complex geometries - Yes, but not yet integrated
      Anything resembling radiosity
"""


def initial_signal(signal_length, fft_output):
    """
    Making the array for the initial signals.
    Input size_fft_two and output_signal
    """
    signal_length2 = int(signal_length // 2)  # Making size_fft_two and setting it as an int again just to be sure
    output_frequency = np.zeros((signal_length2, 3))  # Making output array equivalent to input_array in old code
    throw_array = np.arange(1, signal_length2 + 1)  # Helps get rid of for-loops in old version

    output_frequency[:, 0] = throw_array * Pf.Fs / signal_length  # Tried simplifying the math a bit from original
    output_frequency[:, 1] = abs(fft_output[1:1 + signal_length2] / signal_length)  # Only go up to size_ftt_two
    output_frequency[:, 2] = np.arctan2(np.imag(fft_output[1:1 + signal_length2] / signal_length),
                                        np.real(fft_output[1:1 + signal_length2] / signal_length))

    return output_frequency


def update_freq(dx_update, alpha_update, diffusion_update, lamb, air_absorb):
    """
    Update ray phase and amplitude
    """
    global phase, amplitude  # works directly
    two_pi_dx_update = twopi * dx_update
    ein = phase - (two_pi_dx_update / lamb)
    zwei = ein % twopi
    masque = zwei > np.pi
    drei = masque * zwei - twopi
 
    phase = np.where(masque, drei, ein)
    amplitude *= ((1.0 - alpha_update) * (1.0 - diffusion_update) * np.exp(-air_absorb * dx_update))


def vex(d, f_initial, y, z):
    """The x coordinate of the ray 
    Used for veci"""
    return (d - f_initial[1] * y - f_initial[2] * z) / f_initial[0]


def main():

    global phase
    global amplitude
    global twopi
    twopi = np.pi * 2
    t = time.time()

    # port and import receiver file
    receiver_hit = 0
    ground_hit = 0
    building_hit = 0

    # Initialize counters
    xj = complex(0.0, 1.0)
    radius2 = Pf.radius**2
    ray_sum = 0

    # Initialize receiver variables
    last_receiver = np.zeros(3)
    last_receiver2 = np.zeros(3)
    receiver_point = np.zeros(3)
    receiver_point2 = np.zeros(3)

    # Read in input file
    input_signal = np.loadtxt(Pf.INPUTFILE)
    k = len(input_signal)
    # masque = input_signal > 0
    huge = 1000000.0

    # Allocate the correct size to the signal and fft arrays
    size_fft = k
    size_fft_two = size_fft // 2
    output_signal = np.fft.rfft(input_signal, size_fft)

    # Create initial signal
    frecuencias = initial_signal(size_fft, output_signal)      # Equivalent to inputArray in original
    air_absorb = Fun.absorption(Pf.ps, frecuencias[:, 0], Pf.hr, Pf.Temp)   # size_fft_two
    lamb = Pf.soundspeed/frecuencias[:, 0]     # Used for updating frequencies in update function
    time_array = np.arange(k) / Pf.Fs

    #       Set initial values
    v_initial = np.array([Pf.xinitial, Pf.yinitial, Pf.zinitial])
    xi_initial = np.cos(Pf.phi) * np.sin(Pf.theta)
    n_initial = np.sin(Pf.phi) * np.sin(Pf.theta)
    zeta_initial = np.cos(Pf.theta)
    length = np.sqrt(xi_initial * xi_initial + n_initial * n_initial + zeta_initial * zeta_initial)
    f_initial = np.array([xi_initial, n_initial, zeta_initial])
    d4 = np.dot(f_initial, v_initial)   # equivalent to tmp
    #       Create initial boom array
    #  Roll this all into a function later
    y_space = Pf.boomspacing * abs(np.cos(Pf.phi))
    z_space = Pf.boomspacing * abs(np.sin(Pf.theta))
    if Pf.xmin == Pf.xmax:
        ray_max = int((Pf.ymax - Pf.ymin) / y_space) * int((Pf.zmax - Pf.zmin) / z_space)
        print(ray_max, ' is the ray_max')

    j = np.arange(1, 1 + int((Pf.ymax-Pf.ymin) // y_space))
    k_2 = np.arange(1, 1 + int((Pf.zmax-Pf.zmin) // z_space))
    ray_y = Pf.ymin + j * y_space
    ray_z = Pf.zmin + k_2 * z_space

    boom_carpet = ((vex(d4, f_initial, y, z), y, z) for z in ray_z for y in ray_y)
    # Create a receiver array, include a receiver file.
    alpha_nothing = np.zeros(size_fft_two)

    # Making specific receiver points using receiver module
    Rps.Receiver.initialize(Pf.RecInput)
    ears = Rps.Receiver.rList           # easier to write
    for R in ears:          # hotfix
        R.magnitude = np.zeros(size_fft_two)
        R.direction = np.zeros(size_fft_two)
    temp_receiver = np.array(np.zeros(len(ears)))
    #       Initialize normalization factor
    normalization = (np.pi*radius2)/(Pf.boomspacing**2)

    output_array1 = np.zeros((size_fft_two, 6))
    dh_output_array1 = np.zeros((size_fft_two, 6))

    #       Define ground plane
    ground_height = 0.000000000
    ground_n = np.array([0.000000000, 0.000000000, 1.00000000])
    ground_d = -ground_height

    #     Allocate absorption coefficients for each surface for each frequency
    alpha_ground = np.zeros(size_fft_two)
    for D1 in range(0, size_fft_two):       # This loop has a minimal impact on performance
        if frecuencias[D1, 0] >= 0.0 or frecuencias[D1, 0] < 88.0:
            alpha_ground[D1] = Pf.tempalphaground[0]
        elif frecuencias[D1, 0] >= 88.0 or frecuencias[D1, 0] < 177.0:
            alpha_ground[D1] = Pf.tempalphaground[1]
        elif frecuencias[D1, 0] >= 177.0 or frecuencias[D1, 0] < 355.0:
            alpha_ground[D1] = Pf.tempalphaground[2]
        elif frecuencias[D1, 0] >= 355.0 or frecuencias[D1, 0] < 710.0:
            alpha_ground[D1] = Pf.tempalphaground[3]
        elif frecuencias[D1, 0] >= 710.0 or frecuencias[D1, 0] < 1420.0:
            alpha_ground[D1] = Pf.tempalphaground[4]
        elif frecuencias[D1, 0] >= 1420.0 or frecuencias[D1, 0] < 2840.0:
            alpha_ground[D1] = Pf.tempalphaground[5]
        elif frecuencias[D1, 0] >= 2840.0 or frecuencias[D1, 0] < 5680.0:
            alpha_ground[D1] = Pf.tempalphaground[6]
        elif frecuencias[D1, 0] >= 5680.0 or frecuencias[D1, 0] < frecuencias[size_fft_two, 0]:
            alpha_ground[D1] = Pf.tempalphaground[7]

    alpha_building = np.zeros((Pf.absorbplanes, size_fft_two))
    for W in range(Pf.absorbplanes):        # These also look minimal
        for D2 in range(size_fft_two):
            if frecuencias[D2, 0] >= 0.0 or frecuencias[D2, 0] < 88.0:
                alpha_building[W, D2] = Pf.tempalphabuilding[W, 0]
            elif frecuencias[D2, 0] >= 88.0 or frecuencias[D2, 0] < 177.0:
                alpha_building[W, D2] = Pf.tempalphabuilding[W, 1]
            elif frecuencias[D2, 0] >= 177.0 or frecuencias[D2, 0] < 355.0:
                alpha_building[W, D2] = Pf.tempalphabuilding[W, 2]
            elif frecuencias[D2, 0] >= 355.0 or frecuencias[D2, 0] < 710.0:
                alpha_building[W, D2] = Pf.tempalphabuilding[W, 3]
            elif frecuencias[D2, 0] >= 710.0 or frecuencias[D2, 0] < 1420.0:
                alpha_building[W, D2] = Pf.tempalphabuilding[W, 4]
            elif frecuencias[D2, 0] >= 1420.0 or frecuencias[D2, 0] < 2840.0:
                alpha_building[W, D2] = Pf.tempalphabuilding[W, 5]
            elif frecuencias[D2, 0] >= 2840.0 or frecuencias[D2, 0] < 5680.0:
                alpha_building[W, D2] = Pf.tempalphabuilding[W, 6]
            elif frecuencias[D2, 0] >= 5680.0 or frecuencias[D2, 0] < frecuencias[size_fft_two, 0]:
                alpha_building[W, D2] = Pf.tempalphabuilding[W, 7]

    # This does not appear to be used, so I commented it out -- r0ml
    # D = np.dot(f_initial, v_initial)   # Hotfix  We used this name right above

    #        Mesh the patches for the environment.  Include patching file.
    diffusion_ground = 0.0
    if Pf.radiosity:  # If it exists as a non-zero number
        #    import SingleBuildingGeometry
        diffusion = Pf.radiosity
    else:
        diffusion = 0.0

    ray_counter = 0

    if Pf.h < (2 * Pf.radius):
        print('h is less than 2r')
        raise SystemExit

    # These are for debugging, Uncomment this block and comment out the for loop below
    # ray = 606                     # @ Pf.boomSpacing = 1
    # for i in range(606):
    #      ray =      next(boom_carpet)
    #      ray_counter += 1
    #
    # if ray:
    # Begin tracing
    check_direction = [0, 0, 0]
    n_box = [0, 0, 0]
    veci = np.array([0, 0, 0])
    print('began rays')
    for ray in boom_carpet:              # Written like this for readability
        veci = ray      # initial ray position
        hit_count = 0
        double_hit = 0

        amplitude = frecuencias[:, 1]/normalization
        phase = frecuencias[:, 2]

        f = np.array(f_initial)                                      # Direction
        for I in range(Pf.IMAX):      # Making small steps along the ray path.
            # For each step we should return, location, phase and amplitude
            dx_receiver = huge
            # Find the closest sphere and store that as the distance
            i = 0
            for R in ears:
                # The way that tempReceiver works now, it's only used here and only should be used here.
                # It's not defined inside the receiver because it's ray dependant.
                temp_receiver[i] = R.sphere_check(radius2, f, veci)    # Distance to receiver
                i += 1

                # if receiver_hit >= 1:  # if you hit a receiver last time, don't hit it again
                #     if np.all(R.position == last_receiver):
                #         tempReceiver = huge
                #     if np.all(f == check_direction):
                #         OC = R.position - veci
                #         OCLength = np.dot(OC, OC)
                #         if OCLength < radius2:
                #             tempReceiver = huge
                # if receiver_hit >= 2:
                #     if np.all(R.position == last_receiver):
                #         tempReceiver = huge
                # if tempReceiver < dx_receiver:
                #     dx_receiver = tempReceiver
                #     receiver_point = R.position
                # elif tempReceiver == dx_receiver and tempReceiver != huge:
                #     receiverCheck = tempReceiver

    # We need to double check that double hit actually works.  R2 is not really
    # a thing, we should make sure it is doing what we want.
    #                 if np.all(R.position == receiver_point):
    #                     double_hit = 0
    #                 else:
    #                     R2 = R
    #                     double_hit = 1
    #                     print('double hit')
            temp_receiver[np.where((temp_receiver < (10.0**(-13.0))))] = huge
            tmp = np.argmin(temp_receiver)
            dx_receiver = temp_receiver[tmp]
            if dx_receiver != huge:
                receiver_point = ears[tmp].position

                #     Check Intersection with ground plane
            ground_vd = np.dot(ground_n, f)
#            ground_vd = ground_n[0] * f[0] + ground_n[1] * f[1] + ground_n[2] * f[2]
            if ground_hit == 1:
                dx_ground = huge
            elif ground_vd != 0.0:
                ground_vo = ((np.dot(ground_n, veci)) + ground_d)
                dx_ground = -ground_vo / ground_vd
                if dx_ground < 0.0:
                    dx_ground = huge
            else:
                dx_ground = huge

            #     Check intersection with building
            # dx_building = huge
#            hit=0
#            planeHit = 0
            #     Check intersection with Boxes
            #      for Q in range(0, Bg.BoxNumber):
            #          dxNear, dxFar, hit, planeHit = Fun.box(Bg.BoxArrayNear[Q], Bg.BoxArrayFar[Q], veci, f)
            #          if dxNear < dx_building:
            #              dx_building = dxNear
            #              Vecip1 = veci + np.multiply(dx_building, f)
            #              whichBox = Q
            #              n_box = Fun.plane(Vecip1, Bg.BoxArrayNear[whichBox], Bg.BoxArrayFar[whichBox], planeHit)
            #   Implement Geometry parser
            if building_hit == 1:
                dx_building = huge
            else:
                dx_building, n_box = Gp.collision_check2(Gp.mesh, veci, f)
                # for face in Gp.mesh:
                #     dxnear, nTemp = Gp.collisionCheck(face, veci, f)
                #     if dxnear < dx_building:
                #         dx_building1 = dxnear
                #         n_box1 = nTemp
                # if (ray_counter == 606):
                #     print('original',dx_building,n_box)

            # This part doesn't really work well.  We have not incorporated it.
            # Eventually all interactions will be triangles anyway so I'm leaving it here to be updated.

            #   Check intersection with Triangles
            #        if Bg.TriangleNumber > 0:
            #            for Q in range(0, Bg.TriangleNumber):
            #                dxNear, behind = Fun.Polygon(veci, f, Q, 3, Bg.TriangleNumber, Bg.PointNumbers,
            #                Bg.TriangleArray,
            #                                             Bg.BuildingPoints, normal, FaceNormalNo, FaceNormals)
            #                if dxNear < dx_building:
            #                    dx_building = dxNear
            #                    n_box = normal
            #                    whichBox = Q
            #     Check intersection with Squares
            #        if Bg.SquareNumber > 0:
            #            for Q in range(0, Bg.SquareNumber):
            #                dxNear, behind = Fun.Polygon(veci, f, Q, 4, SquareNumber,
            #                PointNumbers, SquareArray, BuildingPoints,
            #                normal, FaceNormalNo, FaceNormals)
            #                if dxNear < dx_building:
            #                    dx_building = dxNear
            #                    n_box = normal
            #                    whichBox = Q
            building_hit = 0
            receiver_hit = 0
            ground_hit = 0

            #     Check to see if ray hits within step size
            if dx_receiver < Pf.h or dx_ground < Pf.h or dx_building < Pf.h:
                dx = min(dx_receiver, dx_ground, dx_building)
                #  if the ray hits a receiver, store in an array.  If the ray hits two, create two arrays to store in.
        #        for R in ears:
                if dx == dx_receiver:
                    print('Ray ', ray_counter, ' hit receiver ', R.recNumber)
                    veci += (dx * f)
                    # receiver_hit = 1
                    # checkDirection = f
                    # if double_hit == 1:
                    #    receiver_hit = 2
                    hit_count = hit_count + 1
                    update_freq(dx, alpha_nothing, 0, lamb, air_absorb)
                    # last_receiver = receiver_point
                    output_array1[:, 0] = frecuencias[:, 0]
                    output_array1[:, 1:4] = receiver_point[:]
                    output_array1[:, 5] = phase[:]
                    # if double_hit == 1:
                    #    # R2 = R      #Supposed to be other R, but just a placeholder for now
                    #    R.on_hit(amplitude/2, phase)
                    #    R2.on_hit(amplitude/2, phase)
                    # else:
                    ears[tmp].on_hit(amplitude, phase)

                    # if(double_hit==1):
                    #      output_array1[:,4]=amplitude[:]/2.0
                    #      dh_output_array1[:,0]=inputArray[:,0]
                    #      dh_output_array1[:,1:4]=receiver_point2[:]
                    #      dh_output_array1[:,4]=amplitude[:]/2.0
                    #      dh_output_array1[:,5]=phase[:]
                    #      last_receiver2 = receiver_point2
                    # else:
                    #      output_array1[:,4]=amplitude[:]
                    # tempArray=Fun.receiverHITFUNC(size_fft,output_array1,Rps.arraySize,tempArray)
                    # looks like it does the same thing as on_hit. Here later
                    # R.on_hit(amplitude,phase)
                    # if (double_hit==1):
                    #      tempArray=Fun.receiverHITFUNC(size_fft,dh_output_array1,Rps.arraySize,tempArray)
                    #      Using objects may circumvent the need to have this, but it stays for now
                    #      count+=1
                    # count+=1

                if abs(dx - dx_ground) < 10.0**(-13.0):  # If the ray hits the ground then bounce and continue
                    veci += (dx_ground * f)
                    tmp = np.dot(ground_n, veci)
                    if tmp != ground_d:
                        veci[2] = 0
                    print('hit ground at ', I)
                    dot1 = np.dot(f, ground_n)
                    n2 = np.dot(ground_n, ground_n)
                    f -= (2.0 * (dot1 / n2 * ground_n))
#                    length = np.sqrt(np.dot(f, f))
                    ground_hit = 1
#                    twoPiDx = np.pi * 2 * dx_ground
                    #     Loop through all the frequencies
                    update_freq(dx_ground, alpha_ground, diffusion_ground, lamb, air_absorb)
    #                if Pf.radiosity == 1 and (diffusion_ground != 0.0):
    #                    for Q in range(0, PatchNo):
    #                        if formFactors[0, Q, 1] == 1:
    #                            if (veci[0] <= (patchArray[Q, W, 0] + 0.5 * patchArray[Q, W, 3]) and
                #                            veci[0]>=(patchArray[Q, W, 0] - 0.5 * patchArray[Q, W, 3])):
    #                                if veci[1] <= (patchArray[Q, W, 1] + 0.5 * patchArray[Q, W, 4]) and
                    #                                veci[1]>=(patchArray[Q, W, 1] - 0.5 * patchArray[Q, W, 4]):
    #                                    if veci[2] <= (patchArray[Q, W, 2] + 0.5 * patchArray[Q, W, 5]) and
                    #                                    veci[2]>=(patchArray[Q, W, 2] - 0.5 * patchArray[Q, W, 5]):
    #                                        temp2 = complex(abs(patchArray[Q, W, 6])*np.exp(xj*patchArray[Q, W, 7]))
    #                                        temp3 = complex(abs(amplitude[W] * (1.0 - alphaGround[W]) *
                #                                        diffusion_ground *
                #                                        exp(-m * dx_ground)) * exp(1j * phaseFinal))
    #                                        temp4 = temp2 + temp3
    #                                        patchArray[Q, W, 6] = abs(temp4)
    #                                        patchArray[Q, W, 7] = np.arctan(temp4.imag,temp4.real)
                if dx == dx_building:   # if the ray hits the building then change the direction and continue
                    veci += (dx * f)
                    print('hit building at step ', I)
                    n2 = np.dot(n_box, n_box)
                    n_building = n_box / np.sqrt(n2)
                    n3 = np.dot(n_building, n_building)
                    dot1 = np.dot(f, n_building)
                    f -= (2.0 * (dot1 / n3 * n_building))

#                    length = np.sqrt(np.dot(f, f))
                    building_hit = 1
                    # We need to look into complex absorption and see if this is really the best way.
    #                if Pf.complexAbsorption:
    #                    if Pf.absorbPlanes == 2:
    #                        if (veci[2] > 0.0) and (veci[2] < height1):
    #                            alpha = alpha_building[0, :]
    #                        elif veci[2] > height1 and veci[2] <= height2:
    #                            alpha = alpha_building[1, :]
    #                    if Pf.absorbPlanes == 3:
    #                        if veci[2] > height2 and veci[2] <= height3:
    #                            alpha = alpha_building[2, :]
    #                    if Pf.absorbPlanes == 4:
    #                        if veci[2] > height3:
    #                            alpha = alpha_building[4, :]
    #                else:
                    alpha = alpha_building[0, :]
                    update_freq(dx, alpha, diffusion, lamb, air_absorb)
            else:  # If there was no interaction with buildings then proceed with one step.
                veci += (Pf.h * f)
                update_freq(Pf.h, alpha_nothing, 0, lamb, air_absorb)
        ray_counter += 1
        print('finished ray', ray_counter)

    # Radiosity removed for readability

    # Reconstruct the time signal and print to output file
    for R in ears:
        R.time_reconstruct(size_fft)

    print('Writing to output file')
    fileid = Pf.outputfile
    with open(fileid, 'w') as file:
        Fun.header(fileid)

    with open(fileid, 'a') as file:
        for w in range(size_fft):
            Rps.Receiver.timeHeader(file, time_array[w], w)
    print('time: ', time.time()-t)

    # Outputting graphs
    t = time.time()

    # ######################################################################
    # Will eventually be moved to a receiver function,
    # here now for ease of access of others reading this
    # ######################################################################
    import matplotlib.font_manager as fm
    # Font
    stdfont = fm.FontProperties()
    stdfont.set_family('serif')
    stdfont.set_name('Times New Roman')
    stdfont.set_size(20)

    for R in ears:
        # For N wave
        pressure = R.signal
        i = R.recNumber
        # plt.figure(i)
        # plt.figure(num = i, figsize=(19.20, 10.80), dpi=120, facecolor='#eeeeee', edgecolor='r')   # grey
        # plt.figure(num = i, figsize=(19.20, 10.80), dpi=120, facecolor='#e0dae6', edgecolor='r')   # muted lilac
        plt.figure(num=i, figsize=(19.20, 10.80), dpi=120, facecolor='#e6e6fa', edgecolor='r')  # lavender
        # plt.plot(time_array,pressure,'r--')
        plt.grid(True)
        plt.plot(time_array, pressure, '#780303')
        # Labeling axes
        plt.xlabel('Time [s]', fontproperties=stdfont)
        plt.ylabel('Pressure [Pa]', fontproperties=stdfont)
        plt.title('Pressure vs Time of Receiver ' + str(i),
                  fontproperties=stdfont,
                  fontsize=26,
                  fontweight='bold')

        # Saving
        # plt.savefig(Pf.graphName + str(i) + '.png', facecolor='#eeeeee')    # grey
        # plt.savefig(Pf.graphName + str(i) + '.png', facecolor='#e0dae6')    # muted lilac
        plt.savefig(Pf.graphName + str(i) + '.png', facecolor='#e6e6fa')  # lavender
        print('Saved receiver', i)
    print('Graph time: ', time.time() - t)
