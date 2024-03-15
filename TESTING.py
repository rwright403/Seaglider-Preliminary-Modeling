import numpy as np
import matplotlib.pyplot as plt

def intometer(x):
    return 0.0254*x

####CONSTANTS:
RHO_PVC = 1380 #kg/m^3
GRAVITY = 9.81 #m/s^2
TIMESTEP = 0.5 #s
BE_EXTENSION_STEP = intometer(0.625)
#TODO: ISSUE WITH IT BEING 1.75


#input all dimensions as imperial!
class BuoyancyEngine:
    def __init__(self, id, od, cont_len, travel_len, laspeed, midpoint):
        self.id = intometer(id)
        self.od = intometer(od)
        self.cont_len = intometer(cont_len)
        self.travel_len = intometer(travel_len)
        self.midpoint = midpoint #passed in in metric!
        self.laspeed = laspeed #m/s

        #displaced volume
        self.V_cont = 0.25*(np.pi*self.od**2)*self.cont_len
        self.V_mid = 0.25*np.pi*(self.id**2)*self.midpoint + self.V_cont

        #self.V_ext = 0.25*np.pi*(self.id**2)*(midpoint+self.cont_len) #BUG: use inner diameter for whole length, likely the issue here

        self.allowable_ext = midpoint

        self.mass = RHO_PVC*0.25*np.pi*(self.od**2-self.id**2)*self.cont_len

class PressureHull:
    def __init__(self,id,od,len,percent_stability):
        self.id = id
        self.od = od
        self.l_hull = len
        self.percent_stability = percent_stability
        
        self.stability = percent_stability*self.od
        self.V_hull = 0.25*np.pi*self.l_hull*self.od**2

        self.mass = RHO_PVC*(np.pi/4*(self.od**2-self.id**2)*self.l_hull)

def seaglider_trajectory(rho_water, be, hull, midpoint, m_glider, hfoil_coeff):
    max_speed = 0
    glide_period = 0
    output_arr = []

    ###iteratively solve BE extension distances w sim in first while loop
    max_allowable_depth = -20 #m
    sim_depth = 0 #m

    ###remember to save the final valid numbers, probably use a break statement
    while(sim_depth > max_allowable_depth):
        ####START Trajectory Sim in second while loop 
        i = 2

        x_accel_arr = [0, 0]
        y_accel_arr = [0, 0]

        x_vel_arr = [0, 0]
        y_vel_arr = [0, 0]
        
        x_disp_arr = [0, 0]
        y_disp_arr = [0, 0]

        time_arr = [0, 0,TIMESTEP]

        max_speed = 0

        s_x_change_down_arr = []
        s_x_change_up_arr = []        
        be_pos_arr = [0, 0]

        #setup current_len
        current_len = midpoint
        #print(current_len)
        diving = True
        contunuing_criteria = True

        while (contunuing_criteria == True):
            #calc be volume
            V_be = be.V_cont + (0.25*np.pi*be.id**2)*current_len

            #solve glide angle
            theta = np.arctan( ( hull.l_hull * (rho_water*GRAVITY*(0.25*np.pi*be.id**2)*current_len) ) / (m_glider*GRAVITY*hull.stability) )

            #solve lift
            L = 0.5*rho_water*hfoil_coeff*(x_vel_arr[-1]**2 + y_vel_arr[-1]**2)

            #solve Fy, Fx
            if(y_disp_arr[i-1] <= x_disp_arr[i-2]):
                F_y = rho_water*GRAVITY*(hull.V_hull + V_be ) - m_glider*GRAVITY + L*np.sin(theta)
            else:
                F_y = rho_water*GRAVITY*(hull.V_hull + V_be ) - m_glider*GRAVITY - L*np.sin(theta)
            
            F_x = L*np.cos(theta) 

            #solve acceleration
            y_accel_arr.append(F_y/m_glider)
            x_accel_arr.append(F_x/m_glider)

            #solve velocity with trapezoidal method to approx integral
            y_vel_arr.append(y_vel_arr[i-1] + np.trapz(y=y_accel_arr[:i+1],dx=TIMESTEP))
            x_vel_arr.append(x_vel_arr[i-1] + np.trapz(y=x_accel_arr[:i+1],dx=TIMESTEP))

            #solve displacement with trapezoidal method to approx integral
            y_disp_arr.append(y_disp_arr[i-1] + np.trapz(y=y_vel_arr[:i+1],dx=TIMESTEP))
            x_disp_arr.append(x_disp_arr[i-1] + np.trapz(y=x_vel_arr[:i+1],dx=TIMESTEP))

            #print(y_accel_arr[-1], 39.3701*current_len, y_vel_arr[-1])

            #print(y_vel_arr[-1],x_vel_arr[-1])
            #print(180/np.pi *theta, F_x, F_y, y_accel_arr[-1], x_vel_arr[i], y_vel_arr[i])
            print("F_y: ", F_y, "lift: ", L, "y accel: ", y_accel_arr[-1], "y_vel: ", y_vel_arr[-1], "y_disp: ", y_disp_arr[-1], 39.3701*current_len)
            #print(rho_water*GRAVITY*(hull.V_hull + V_be ), m_glider*GRAVITY, L*np.sin(theta))

            if( current_len <= (midpoint - be.allowable_ext) and diving == True):
                diving = False
                s_x_change_up_arr.append(x_disp_arr[i])
                #print("going up", s_x)

            if( current_len >= (midpoint + be.allowable_ext) and diving == False):
                diving = True
                contunuing_criteria = False
                s_x_change_down_arr.append(x_disp_arr[i])
                #print("going down", s_x)

            if diving == True:
                current_len -= be.laspeed*TIMESTEP
            else:
                current_len += be.laspeed*TIMESTEP

            #print(theta, F_y, v_y, v_x, s_x, current_len)
            
            time_arr.append(time_arr[i] + TIMESTEP)
            be_pos_arr.append(39.3701*current_len)
            i+= 1

              
        
        sim_depth = np.min(y_disp_arr)
        max_speed = np.max(x_vel_arr)
        glide_period = x_disp_arr[-1]

        #after sim check max allowable depth. If we havent hit save, if we have break out of while loop
        #if(sim_depth < max_allowable_depth):
        #    be.allowable_ext = max_be_extension

        ###iterate!!!!!
        be.allowable_ext+= BE_EXTENSION_STEP

        #print(sim_depth)



    print("allowable extension: ",be.allowable_ext*39.3701, " (in)")
    print("max fwd speed: ",max_speed, " (m/s)")
    print("glide period: ",glide_period, '\n')

    print("be id: ",39.3701*be.id, " (in)")
    print("hull id: ",39.3701*hull.id, " (in)")
    print("total glider mass: ", m_glider, " (kg)")


    output_arr.append(be.allowable_ext)
    output_arr.append(max_speed)
    output_arr.append(glide_period)

    output_arr.append(be)
    output_arr.append(hull)


    plt.subplot(1,2,1)  
    plt.plot(x_disp_arr, y_disp_arr, label='Seaglider Position', color='blue')
    for i in s_x_change_down_arr:
        plt.axvline(x=i, color='r', linestyle='--', label='change up')
    for i in s_x_change_up_arr:
        plt.axvline(x=i, color='g', linestyle='--', label='change up')

    plt.xlabel('X-Pos (m)')
    plt.ylabel('Y-Pos (m)')

    plt.subplot(1,2,2)
    plt.plot(x_disp_arr, be_pos_arr, label='BE Position Over X-Pos', color='blue')
    plt.xlabel('X-Pos (m)')
    plt.ylabel('B.E. Pos (in)')
    plt.show()

    return output_arr




####inputs:
hfoil_coeff = 0.008 #Area of wing * Coeff Lift
percent_stability = 0.15 #%
rho_water =997 #kg/m^3
midpoint = intometer(4)
internal_mass = 11 #kg

#geometry
buoyeng = BuoyancyEngine(4.5,5.515,10.3, 8, 0.08*0.0254,midpoint)

hull_id = intometer(5.0)
hull_od = intometer(5.5)

hull_len = (internal_mass + buoyeng.mass - rho_water*buoyeng.V_mid) / ( np.pi *0.25* (rho_water*(hull_od**2) - RHO_PVC*(hull_od**2 - hull_id**2)) )
print(buoyeng.mass)

print("hull_len: ", 39.3701*hull_len, " (in)")

preshull = PressureHull( hull_id, hull_od, hull_len, percent_stability)
m_glider = preshull.mass + buoyeng.mass +internal_mass #1.75 #rho_water*(preshull.V_hull+buoyeng.V_ext)+1.75

seaglider_trajectory(rho_water, buoyeng, preshull, midpoint, m_glider, hfoil_coeff)
#print(39.3701*buoyeng.allowable_ext)

#print(preshull.V_hull+  (0.25*np.pi*be.id**2)*current_len)