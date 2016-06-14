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


