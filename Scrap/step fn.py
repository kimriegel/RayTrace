#for I in range(Pf.IMAX):      # Making small steps along the ray path.
# add position, direction, phase, amplitude


I=3
for k in range(3):
    step= (i*i for i in range(I))   #aka rayPath
    #print(next(step))
    for j in step:
        print(j)

#boom_carpet = ((vex(d4, f_initial, y, z), y, z) for z in ray_z for y in ray_y)
def step():
    """
    Each step along the path a single ray takes.
    Finds distance to all surface types
    Finds closest surface
    Updates ray data and direction based on answer
    """

    # Using building as sample as it will be changed least in full release
    #    #   Implement Geometry parser
    #if building_hit == 1:   #Avoid hitting building twice
    #    dx_building = huge
    #else:
    dx_building, n_box = Gp.collision_check2(Gp.mesh, veci, f)

    # Find distance to building and plane normal

    #     Check to see if ray hits within step size
    if dx_receiver < Pf.h or dx_ground < Pf.h or dx_building < Pf.h:
        dx = min(dx_receiver, dx_ground, dx_building)
        
        if dx == dx_building:   # if the ray hits the building then change the direction and continue
            veci += (dx * f)
            #print('hit building at step ', I)
            n2 = np.dot(n_box, n_box)
            n_building = n_box / np.sqrt(n2)
            n3 = np.dot(n_building, n_building)
            dot1 = np.dot(f, n_building)
            f -= (2.0 * (dot1 / n3 * n_building))
            building_hit = 1
            alpha = alpha_building[0, :]
            update_freq(dx, alpha, diffusion, lamb, air_absorb)
    else:  # If there was no interaction with buildings then proceed with one step.
        veci += (Pf.h * f)
        update_freq(Pf.h, alpha_nothing, 0, lamb, air_absorb)



rayPath = (step(i) for i in range(Pf.IMAX) )  
for I in rayPath:
#for I in rayPath:
#rayPath = (x for step in range(Pf.IMAX) )  
#x is output, not sure there even is one
# step is everything that happens within I loop now


#