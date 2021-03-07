# Table of Contents:
# 1.   absorption
# 2.   receiver_hit_func
# 3.   header
# 4.   sphere_check
# 5.   cross
# 6.   tri
# 7.   plane
# 8.   box
# 9.   rotation
# 10. initial_signal
# 11. update_freq
# 12. vec
import numpy as np
import math as m
import Parameterfile as Pf
import config

HUGE = 1000000.0
XJ = complex(0, 1)


def absorption(ps, freq, hr, temp):
    '''
    This function computes the air absorption for a given frequency,
    ambient pressure, relative humidity, and temperature.

    ps, hr, and temp are defined in Parameterfile.py
    '''
    if temp == 0 or ps == 0:
        raise ValueError('Cannot divide by Zero!')
    if temp < 0:
        raise ValueError('Temperature is in Kelvin. Non-negative values only.')

    # Define all variables and reference values
    ps0 = 1.0
    t0 = 293.15
    t01 = 273.16
    f = freq/ps

    # Compute all relevant parameters
    p_sat = ps0*10**(-6.8346 * (t01 / temp) ** 1.261 + 4.6151)
    h = ps0*(hr/ps)*(p_sat/ps0)
    fr_n = 1 / ps0 * (t0 / temp) ** (1 / 2) * (9 + 280 * h * m.exp(-4.17 * ((t0 / temp) ** (1 / 3) - 1)))
    fr_o = 1/ps0*(24+4.04*10**4*h*((.02+h)/(.391+h)))
    term1 = 0.01275*(m.exp((-2239.1 / temp)) / (fr_o + f ** 2 / fr_o))
    term2 = 0.1068*(m.exp(-3352 / temp) / (fr_n + f ** 2 / fr_n))
    absorb = ps0*f**2*((1.84 * 10 ** (-11.0) * (temp / t0) ** 0.5 * ps0) + (temp / t0) **
                       (-5.0 / 2.0) * (term1 + term2))
    return absorb


def receiver_hit_func(size_fft, output_array, array_size, temp_array):
    '''
    This Function adds the pressures from a ray when it hits the receiver.
    '''

    # Define arrays with numpy zeros function
    xj = complex(0, 1)
    
    # Add new pressures to existing pressures in temp_array
    # First Look for the correct location.
    for D in range(0, array_size):
        # If the location is the same, loop through the frequency and add current values with new values.
        for W in range(0, (size_fft // 2)):
            temp1 = abs(temp_array[D, W, 4]) * m.e ** (xj * temp_array[D, W, 5])
            temp2 = abs(output_array[W, 4]) * m.e ** (xj * output_array[W, 5])
            temp3 = temp1+temp2
            temp_array[D, W, 4] = abs(temp3)  # magnitude
            temp_array[D, W, 5] = np.arctan2(np.imag(temp3), np.real(temp3))  # direction
    # imag_part and real_part in original code fortran functions
    # using numpy .real and .imag function to achieve same result (Probably)

    return temp_array


def header(output_file):
    '''
    This function prints the header for the tecplot data.
    '''

    f = open(output_file, "w")

    f.write('TITLE = "Pressure at earlevel"\n')
    f.write('VARIABLES = "X[m]" "Y[m]" "Z[m]" "P[Pa]"\n')
    f.write('TEXT\n')
    f.write('CS=FRAME\n')
    f.write('X=71.9660948264,Y=82.9866270431\n')
    f.write('C=BLACK\n')
    f.write('S=LOCAL\n')
    f.write('HU=POINT\n')
    f.write('LS=1 AN=MIDCENTER\n')
    f.write('BX=Filled BXM=60 LT=0.1 BXO=BLACK BXF=WHITE\n')
    f.write('F=HELV\n')
    f.write('H=20 A=0\n')
    f.write('MFC=""\n')
    f.write('CLIPPING=CLIPTOVIEWPORT\n')
    f.write('T="Time = &(SOLUTIONTIME%4f)" \n')
    printer = 0
    return printer


def sphere_check(sc, sr2, f, veci):
    '''
    This function performs a check wwhether a ray hits a sphere. If
    it does hit the function returns the distance to the sphere
    '''

    huge_check = 1000000.0
    oc = np.zeros(3)

    oc[0] = sc[0] - veci[0]
    oc[1] = sc[1] - veci[1]
    oc[2] = sc[2] - veci[2]
    l2_oc = np.dot(oc, oc)
    tca = np.dot(oc, f)
    t2hc = sr2 - l2_oc + tca ** 2
    if l2_oc == sr2:
        dx = huge_check
    elif tca < 0.0:
        dx = huge_check
    elif t2hc < 0.0:
        dx = huge_check
    else:
        dx = tca - (t2hc**(1/2))
    return dx


def cross(a, b):
    '''
    This function calculates a cross product of A and B and returns normal
    '''

    if len(a) != 3 or len(b) != 3:
        raise ValueError('Input must be 3D vector.')

    normal = np.zeros(3)
    normal[0] = a[1] * b[2] - a[2] * b[1]
    normal[1] = a[2] * b[0] - a[0] * b[2]
    normal[2] = a[0] * b[1] - a[1] * b[0]
    length = (normal[0]**2.0 + normal[1]**2 + normal[2] ** 2)**(1/2)
    if length != 0.0:
        normal = normal/length
    return normal


def tri(veci, f, q, number, point_numbers, poly_array, v, normal, face_normal_no, vn, dxbuilding, behind):
    """
    Lab notebook 5/16
    This is an attempt to merge box and polygon into one function since
    we are working entirely in triangular meshes now
    ********************************Untested***********************************
    A 1:1 translation was made from Fortran. This is the closest match to
    what we are trying to do with triangle geometry. However there is no
    readily available geometry file to test this.
    [No Description given in Fortran]
    """
    size = 3
    g = np.zeros((size, 2))  
    # inits
    nc = 0
    behind = 0
    normal = vn[poly_array[q, 1], :]
    d = -np.dot(normal, ValueError(poly_array[q, 2]))
    v_d = np.dot(normal, f)
    if v_d >= 0.0:
        dxbuilding = HUGE
    v_0 = -(np.dot(normal, veci)+d)
    t = v_0/v_d
    if t < 0.0:
        dxbuilding = HUGE
        behind = 1
        # Stage 1
    intersection = veci + f*t
    maximum = max(abs(normal))

    if maximum == abs(normal[0]):
        for P in range(size):
            g[P, :] = (intersection[1]-v[int(poly_array[q, 1+P]), 1], intersection[2]-v[int(poly_array[q, 1+P]), 2])
    elif maximum == normal[1]:
        for P in range(size):
            g[P, :] = (intersection[0]-v[int(poly_array[q, 1+P]), 0], intersection[2]-v[int(poly_array[q, 1+P]), 2])
    elif maximum == normal[2]:
        for P in range(size):
            g[P, :] = (intersection[0]-v[int(poly_array[q, 1+P]), 0], intersection[1]-v[int(poly_array[q, 1+P]), 1])
    # Stage 2
    for P in range(size):
        if P == size:
            if g[P, 1] < 0.0:
                sh = -1
            else:
                sh = 1
        elif (g[P, 0] > 0.0) or (g[P+1, 1] > 0.0):
            if (g[P, 0]-(g[P, 1]*(g[P+1, 0]-g[P, 0])/(g[P+1, 1]-g[P, 1]))) > 0.0:
                nc += 1
        odd = nc % 2    # get remainder to find if odd or not
        if odd:
            dxbuilding = t
        else:
            dxbuilding = HUGE
        return dxbuilding, behind


def plane(vecip1, b1, b2, plane_hit):
    '''
    This function calculates the normal at the hitpoint of a box.
    '''

    n_box = [0, 0, 0]
    if plane_hit == 1:
        if vecip1[0] == b1[0]:
            point2 = [b1[0], b1[1], b2[2]]
            point3 = [b1[0], b2[1], b1[2]]
            n_box = cross(np.subtract(point2, b1), np.subtract(point3, b1))

        elif vecip1[0] == b2[0]:
            point1 = (b2[0], b1[1], b1[2])
            point2 = (b2[0], b1[1], b2[2])
            point3 = (b2[0], b2[1], b1[2])

            n_box = cross(np.subtract(point3, point1), np.subtract(point2, point1))

    if plane_hit == 2:

        if vecip1[1] == b1[1]:
            point2 = (b2[0], b1[1], b1[2])
            point3 = (b1[0], b1[1], b2[2])
            n_box = cross(np.subtract(point2, b1), np.subtract(point3, b1))
        elif vecip1[1] == b2[1]:
            point1 = (b1[0], b2[1], b1[0])
            point2 = (b1[0], b2[1], b2[2])
            point3 = (b2[0], b2[1], b1[2])
            n_box = cross(np.subtract(point2, point1), np.subtract(point3, point1))
    if plane_hit == 3:
        if vecip1[2] == b1[2]:
            point2 = (b2[0], b1[1], b1[2])
            point3 = (b1[0], b2[1], b1[2])
            n_box = cross(np.subtract(point3, b1), np.subtract(point2, b1))
        elif vecip1[2] == b2[2]:
            point2 = (b1[0], b2[1], b2[2])
            point3 = (b2[0], b1[1], b2[2])
            n_box = cross(np.subtract(point2, b2), np.subtract(point3, b2))

    return n_box


def box(b1, b2, vecip1, f):
    '''
    Defunct with meshes

    This function checks to see if the ray hits a box. It determines which
    plane the ray hits.
    t1x is the distance to the close side
    t2x is the distance to the far side
    '''

    hit = 5
    huge_box = 1000000.0
    dx_near = -huge_box
    dx_far = huge_box
    temp_f = f
    plane_hit = 0
    if (f[0] == 0.0) or (f[1] == 0.0) or (f[2] == 0.0):
        if f[0] == 0.0:
            if (vecip1[0] < b1[0]) or (vecip1[0] > b2[0]):
                hit = 0
                dx_near = huge_box
                return dx_near, dx_far, hit, plane_hit
        if f[1] == 0.0:
            if (vecip1[1] < b1[1]) or (vecip1[1] > b2[1]):
                hit = 0
                dx_near = huge_box
                return dx_near, dx_far, hit, plane_hit
        if f[2] == 0.0:
            if (vecip1[2] < b1[2]) or (vecip1[2] > b2[2]):
                hit = 0
                dx_near = huge_box
                return dx_near, dx_far, hit, plane_hit
    if hit != 0.0:

        if f[0] == 0.0:
            temp_f[0] = 1.0
        if f[1] == 0.0:
            temp_f[1] = 1.0
        if f[2] == 0.0:
            temp_f[2] = 1.0
        if f[0] != 0.0:
            t1x = (b1[0] - vecip1[0]) / temp_f[0]
            t2_x = (b2[0] - vecip1[0]) / temp_f[0]
            if t1x > t2_x:
                tmp = t1x
                t1x = t2_x
                t2_x = tmp
            if t1x > dx_near:
                dx_near = t1x
            if t2_x < dx_far:
                dx_far = t2_x
            if dx_near > dx_far:
                hit = 0
                dx_near = huge_box
                return dx_near, dx_far, hit, plane_hit
            elif dx_far < 0.0:
                hit = 0
                dx_near = huge_box
                return dx_near, dx_far, hit, plane_hit
        if f[1] != 0.0:
            t1_y = (b1[1] - vecip1[1]) / temp_f[1]
            t2_y = (b2[1] - vecip1[1]) / temp_f[1]
            if t1_y > t2_y:
                tmp = t1_y
                t1_y = t2_y
                t2_y = tmp
            if t1_y > dx_near:
                dx_near = t1_y
            if t2_y < dx_far:
                dx_far = t2_y
            if dx_near > dx_far:
                hit = 0
                dx_near = huge_box
                return dx_near, dx_far, hit, plane_hit
            elif dx_far < 0.0:
                hit = 0
                dx_near = huge_box
                return dx_near, dx_far, hit, plane_hit


        if f[2] != 0.0:
            t1_z = (b1[2] - vecip1[2]) / temp_f[2]
            t2_z = (b2[2] - vecip1[2]) / temp_f[2]
            if t1_z > t2_z:
                tmp = t1_z
                t1_z = t2_z
                t2_z = tmp
            if t1_z > dx_near:
                dx_near = t1_z
            if t2_z < dx_far:
                dx_far = t2_z
            if dx_near > dx_far:
                hit = 0
                dx_near = huge_box
                return dx_near, dx_far, hit, plane_hit
            elif dx_far < 0:
                hit = 0
                dx_near = huge_box
                return dx_near, dx_far, hit, plane_hit
            elif dx_near < 0:
                hit = 0
                dx_near = huge_box
                return dx_near, dx_far, hit, plane_hit
    if hit != 0:
        if dx_near < dx_far:
            hit = 1
            if dx_near == t1x:
                plane_hit = 1
            if dx_near == t1_y:
                plane_hit = 2
            if dx_near == t1_z:
                plane_hit = 3
    return dx_near, dx_far, hit, plane_hit


def rotation(axis, angle, rotation_matrix):
    '''
    What DOES this do?
    Do we call this?
    '''

    rotation_matrix[1, 1] = axis[1] ** 2 + (1 - axis[1] ** 2) * m.cos(angle)
    rotation_matrix[1, 2] = axis[1] * axis[2] * (1 - m.cos(angle)) + axis[3] * m.sin(angle)
    rotation_matrix[1, 3] = axis[1] * axis[3] * (1 - m.cos(angle)) - axis[2] * m.sin(angle)
    rotation_matrix[2, 1] = axis[1] * axis[2] * (1 - m.cos(angle)) - axis[3] * m.sin(angle)
    rotation_matrix[2, 2] = axis[2] ** 2 + (1 - axis[2] ** 2) * m.cos(angle)
    rotation_matrix[2, 3] = axis[2] * axis[3] * (1 - m.cos(angle)) + axis[1] * m.sin(angle)
    rotation_matrix[3, 1] = axis[1] * axis[3] * (1 - m.cos(angle)) + axis[2] * m.sin(angle)
    rotation_matrix[3, 2] = axis[2] * axis[3] * (1 - m.cos(angle)) - axis[1] * m.sin(angle)
    rotation_matrix[3, 3] = axis[3] ** 2 + (1 - axis[3] ** 2) * m.cos(angle)
    return rotation_matrix


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
    two_pi_dx_update = config.twopi * dx_update
    ein = phase - (two_pi_dx_update / lamb)
    zwei = ein % config.twopi
    masque = zwei > np.pi
    drei = masque * zwei - config.twopi 

    phase = np.where(masque, drei, ein)
    amplitude *= ((1.0 - alpha_update) * (1.0 - diffusion_update) * np.exp(-air_absorb * dx_update))


def vex(d, f_initial, y, z):
    """The x coordinate of the ray 
    Used for veci"""
    return (d - f_initial[1] * y - f_initial[2] * z) / f_initial[0]
