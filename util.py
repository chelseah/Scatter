#!/usr/bin/env python
def check_for_bad_dt(sim):
    bad_dt = 0

    par = sim.particles
    p0 = sim.particles[0]                                      #star
    p = sim.particles[1]                            #planet
    dx = p.x - p0.x
    dy = p.y - p0.y
    dz = p.z - p0.z
    mh = (p.m/(3*p0.m))**(1./3.)
    rhill_p = mh*(dx*dx + dy*dy + dz*dz)**(0.5)               #hill radius planet
    rh_sum = rhill_p
    vmax = 0                                                  #max relative velocity squared
    for i in xrange(2,sim.N):
        dx_i = par[i].x - p0.x
        dy_i = par[i].y - p0.y
        dz_i = par[i].z - p0.z
        mh_i = (par[i].m/(3*p0.m))**(1./3.)
        rhill_p_i = mh_i*(dx_i*dx_i + dy_i*dy_i + dz_i*dz_i)**(0.5)
        rh_sum+=rhill_p_i
        dvx = par[i].vx - p.vx
        dvy = par[i].vy - p.vy
        dvz = par[i].vz - p.vz
        v = (dvx*dvx + dvy*dvy + dvz*dvz)**(0.5)
        print 'v',v,dvx,dvy,dvz
        if v > vmax:
            vmax = v
    HSR = sim.ri_hybarid.switch_radius
    min_dt = HSR*rh_sum / vmax
    print 'vmax',vmax,HSR*rh_sum,min_dt,HSR
    if(min_dt < 4*sim.dt):
        bad_dt = 1                                           #factor of 4 for extra wiggle room
    return bad_dt

def set_hill(mass_pl):
    #set up the giant planets so that they will be instable at one point
    i_hill=0
    a_pl=np.zeros(3)
    while (a_pl[0]<1):
        while (i_hill<1):
            spacing=np.random.rand(3)
            print 'spacing,',spacing
            a_inner=a_min+(a_max-a_min)*np.random.rand()
            a_max_i=a_max
            print 'amax_i,',a_max_i
            a_aux=a_inner+(a_max_i-a_inner)*spacing
            print 'a_aux,',a_aux
            a_pl=np.sort(a_aux)

            R_hill_12=((mass_pl[0]+mass_pl[1])/3)**(1./3.)*(a_pl[0]+a_pl[1])/2.
            R_hill_23=((mass_pl[1]+mass_pl[2])/3)**(1./3.)*(a_pl[1]+a_pl[2])/2.
            diff_a_12=abs(a_pl[1]-a_pl[0])
            diff_a_23=abs(a_pl[2]-a_pl[1])

            if(diff_a_12>k_Hill*R_hill_12 and diff_a_23>k_Hill*R_hill_23):
                break
    return [a_inner,a_pl]

