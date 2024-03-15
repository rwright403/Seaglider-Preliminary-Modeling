import numpy as np
import matplotlib.pyplot as plt

#Constants
RHO_PVC = 1380 #kg/m^3
RHO_WATER =997 #kg/m^3
GRAVITY = 9.81 #m/s^2
TIMESTEP = 0.01 #s

def intometer(x):
    return 0.0254*x

class BuoyancyEngine:
    def __init__(self, id, od, cont_len, travel_len, laspeed, midpoint):
        self.id = intometer(id)
        self.od = intometer(od)
        self.cont_len = intometer(cont_len)
        self.travel_len = intometer(travel_len)
        self.midpoint = midpoint #passed in in metric!
        self.laspeed = laspeed #m/s

        print("init be: midpoint, travel len", self.midpoint, self.travel_len)

        #displaced volume
        self.V_cont = 0.25*(np.pi*self.od**2)*self.cont_len
        self.V_mid = 0.25*np.pi*(self.id**2)*self.midpoint + self.V_cont

        #self.V_ext = 0.25*np.pi*(self.id**2)*(midpoint+self.cont_len) #BUG: use inner diameter for whole length, likely the issue here

        self.allowable_ext = midpoint

        self.mass = RHO_PVC*0.25*np.pi*(self.od**2-self.id**2)*self.cont_len
        #print("be mass: ", self.mass)


class PressureHull:
    def __init__(self,id,od,len,percent_stability):
        self.id = id
        self.od = od
        self.l_hull = len
        self.percent_stability = percent_stability

        #print(self.od)
        
        self.stability = percent_stability*self.od
        self.V_hull = 0.25*np.pi*self.l_hull*self.od**2

        self.mass = RHO_PVC*(np.pi/4*(self.od**2-self.id**2)*self.l_hull)
        #print("hull mass: ", self.mass)

class SeagliderFSM:
    def __init__(self, be, hull, m_glider):
        self.be = be
        self.hull = hull
        self.m_glider = m_glider

        self.current_len = midpoint

        self.i = 1

        self.x_accel_arr = [0, 0]
        self.y_accel_arr = [0, 0]

        self.x_vel_arr = [0, 0]
        self.y_vel_arr = [0, 0]
        
        self.x_disp_arr = [0, 0]
        self.y_disp_arr = [0, 0]

        self.time_arr = [0, 0]

        self.max_speed = 0

        self.s_x_change_down_arr = []
        self.s_x_change_up_arr = []        
        self.be_pos_arr = [4,4]

        #setup current_len
        self.current_len = midpoint
        self.V_be = 0
        self.theta = 0

        self.F_y = 0
        self.F_x = 0

        #cg location
        self.X_cg = ( (hull.V_hull*0.5*hull.l_hull) + (be.V_cont*(hull.l_hull+0.5*be.cont_len)) + ((be.V_mid-be.V_cont)*(hull.l_hull+be.cont_len+0.5*be.midpoint)) ) / (hull.V_hull + be.V_mid)
        #print ("X_cg: ", 39.3701*self.X_cg, " (in)")
        #print("hull: ", hull.V_hull, " l_hull: ", hull.l_hull, " be.V_cont  ", be.V_cont, " hull.l_hull: ", hull.l_hull," be.cont_len: ", be.cont_len, " be.V_mid: ",be.V_mid, " be.V_cont: ", be.V_cont)
        self.X_cb = self.X_cg

        # Define the states
        self.states = ['hold_be_pos_up', 'hold_be_pos_down', 'move_down_deccel', 'move_up_accel', 'move_down_accel', 'move_up_deccel', 'end']
        # Define the initial state
        self.current_state = 'move_down_accel'

    def calc_theta(self, current_len, V_be):
        self.X_cb = ( (hull.V_hull*0.5*hull.l_hull) + (be.V_cont*(hull.l_hull+0.5*be.cont_len)) + ((be.V_mid-be.V_cont)*(hull.l_hull+be.cont_len+0.5*current_len)) ) / (hull.V_hull + V_be)       
        delta_X_cb = self.X_cb-self.X_cg
        self.theta = np.arccos(hull.stability/np.sqrt(delta_X_cb**2 +hull.stability**2))
        #print("theta: ", theta)


    def kinematics(self):
        #solve acceleration
        self.y_accel_arr.append(self.F_y/m_glider)
        self.x_accel_arr.append(self.F_x/m_glider)

        #print(self.F_x/m_glider, self.y_accel_arr[-1])

        #solve velocity with trapezoidal method to approx integral
        self.y_vel_arr.append(self.y_vel_arr[self.i-1] + np.trapz(y=self.y_accel_arr[self.i-1:self.i+1],dx=TIMESTEP))
        self.x_vel_arr.append(self.x_vel_arr[self.i-1] + np.trapz(y=self.x_accel_arr[self.i-1:self.i+1],dx=TIMESTEP))

        #print(self.y_vel_arr[-1])

        #solve displacement with trapezoidal method to approx integral
        self.y_disp_arr.append(self.y_disp_arr[self.i-1] + np.trapz(y=self.y_vel_arr[self.i-1:self.i+1],dx=TIMESTEP))
        self.x_disp_arr.append(self.x_disp_arr[self.i-1] + np.trapz(y=self.x_vel_arr[self.i-1:self.i+1],dx=TIMESTEP))

        #"y accel: ", self.y_accel_arr[-1],
        print("x accel: ", self.x_accel_arr[-1], "y accel: ", self.y_accel_arr[-1], "x vel: ", self.x_vel_arr[-1],"y disp: ", self.y_disp_arr[-1], "x disp: ", self.x_disp_arr[-1])
        self.i+=1
        #print(self.y_disp_arr[-1])


    def trajectory(self, dist_to_retract, time_hold_down, dist_to_extend, time_hold_up):
        self.X_cb = ( (hull.V_hull*0.5*hull.l_hull) + (be.V_cont*(hull.l_hull+0.5*be.cont_len)) + ((be.V_mid-be.V_cont)*(hull.l_hull+be.cont_len+0.5*be.midpoint)) ) / (hull.V_hull + be.V_mid)
        self.X_cg = ( (hull.V_hull*0.5*hull.l_hull) + (be.V_cont*(hull.l_hull+0.5*be.cont_len)) + ((be.V_mid-be.V_cont)*(hull.l_hull+be.cont_len+0.5*be.midpoint)) ) / (hull.V_hull + be.V_mid)

        self.dist_to_retract = dist_to_retract
        self.time_hold_down = time_hold_down
        self.dist_to_extend = dist_to_extend
        self.time_hold_up = time_hold_up

        print("start trajectory", self.current_state)


        if self.current_state == 'hold_be_pos_up':
            print("state: hold be pos", "F_y: ", self.F_y, "current_len: ", self.current_len)
            counter = 0
            print("holding time (s): ", self.time_hold_up)

            while counter <= self.time_hold_up:
                
                self.L = 0.5*RHO_WATER*hfoil_coeff*(self.x_vel_arr[-1]**2 + self.y_vel_arr[-1]**2)
                self.D = 0.5*drag_coeff*RHO_WATER*(self.x_vel_arr[-1]**2 + self.y_vel_arr[-1]**2)*(0.25*np.pi*hull.od**2)


                self.F_y = RHO_WATER*GRAVITY*(hull.V_hull + self.V_be ) - m_glider*GRAVITY - self.L*np.cos(self.theta) + self.D*np.sin(self.theta)
                self.F_x = self.L*np.sin(self.theta) - self.D*np.cos(self.theta)
                #print("F_y (N): ", self.F_y, "F_x (N): ", self.F_x, "current len (in): ", self.current_len*39.3701)

                self.kinematics()
                #iterate counter
                counter += TIMESTEP
                #print("theta! ", (180/np.pi)*self.theta, "F_y: ", self.F_y)
                self.be_pos_arr.append(39.3701*self.current_len)
                self.time_arr.append(self.time_arr[-1]+TIMESTEP)
            
            print("state entering: move up deccel", self.x_accel_arr[-1], self.F_y)
            self.current_state = 'move_up_deccel'
            print(self.current_state)

        if self.current_state == 'hold_be_pos_down':
            print("state: hold be pos", "F_y: ", self.F_y, "current_len: ", self.current_len)
            counter = 0
            print("holding time (s): ", self.time_hold_down)

            while counter <= self.time_hold_down:
                self.L = 0.5*RHO_WATER*hfoil_coeff*(self.x_vel_arr[-1]**2 + self.y_vel_arr[-1]**2)
                self.D = 0.5*drag_coeff*RHO_WATER*(self.x_vel_arr[-1]**2 + self.y_vel_arr[-1]**2)*(0.25*np.pi*hull.od**2)
                
                self.F_y = RHO_WATER*GRAVITY*(hull.V_hull + self.V_be ) - m_glider*GRAVITY + self.L*np.cos(self.theta) + self.D*np.sin(self.theta)
                self.F_x = self.L*np.sin(self.theta) - self.D*np.cos(self.theta)
                #print("F_y (N): ", self.F_y, "F_x (N): ", self.F_x, "current len (in): ", self.current_len*39.3701)

                self.kinematics()
                #iterate counter
                counter += TIMESTEP
                self.be_pos_arr.append(39.3701*self.current_len)
                self.time_arr.append(self.time_arr[-1]+TIMESTEP)
            
            print("state entering: move down deccel", self.x_accel_arr[-1], self.F_y)
            self.current_state = 'move_down_deccel'

            

        elif self.current_state == 'move_down_deccel':
            print("state: move down deccel")

            while self.current_len < midpoint:

                self.V_be = be.V_cont + (0.25*np.pi*be.id**2)*self.current_len
                self.calc_theta(self.current_len, self.V_be)

                self.L = 0.5*RHO_WATER*hfoil_coeff*(self.x_vel_arr[-1]**2 + self.y_vel_arr[-1]**2)
                self.D = 0.5*drag_coeff*RHO_WATER*(self.x_vel_arr[-1]**2 + self.y_vel_arr[-1]**2)*(0.25*np.pi*hull.od**2)


                self.F_y = RHO_WATER*GRAVITY*(hull.V_hull + self.V_be ) - m_glider*GRAVITY + self.L*np.cos(self.theta) + self.D*np.sin(self.theta)
                self.F_x = self.L*np.sin(self.theta) - self.D*np.cos(self.theta)
                #print("F_y (N): ", self.F_y, "F_x (N): ", self.F_x, "current len (in): ", self.current_len*39.3701)
                #print("theta! ", (180/np.pi)*theta, "F_y: ", self.F_y)
                #print("F_y (N): ", self.F_y, "y lift: ", self.L*np.cos(self.theta), "y drag: ", self.D*np.sin(self.theta), "Fg: ", - m_glider*GRAVITY, "Fb: ", RHO_WATER*GRAVITY*(hull.V_hull + self.V_be ), self.current_len)
                
                self.kinematics()

                #iterate current_len
                self.current_len += be.laspeed*TIMESTEP
                self.be_pos_arr.append(39.3701*self.current_len)
                self.time_arr.append(self.time_arr[-1]+TIMESTEP)
                #print(self.F_y, 39.3701*self.current_len)

            print("state entering: move up accel", self.x_disp_arr[-1])
            self.current_state = 'move_up_accel'
    

        elif self.current_state == 'move_up_accel':
            print("state: move up accel")
            while self.current_len < dist_to_extend:
                
                self.V_be = be.V_cont + (0.25*np.pi*be.id**2)*self.current_len
                self.calc_theta(self.current_len, self.V_be)
            
                self.L = 0.5*RHO_WATER*hfoil_coeff*(self.x_vel_arr[-1]**2 + self.y_vel_arr[-1]**2)
                self.D = 0.5*drag_coeff*RHO_WATER*(self.x_vel_arr[-1]**2 + self.y_vel_arr[-1]**2)*(0.25*np.pi*hull.od**2)


                self.F_y = RHO_WATER*GRAVITY*(hull.V_hull + self.V_be ) - m_glider*GRAVITY - self.L*np.cos(self.theta) - self.D*np.sin(self.theta)
                self.F_x = self.L*np.sin(self.theta) - self.D*np.cos(self.theta)

                #print("F_y (N): ", self.F_y, "F_x (N): ", self.F_x, "current len (in): ", self.current_len*39.3701)
                #print("theta! ", (180/np.pi)*theta, "F_y: ", self.F_y)

                #print("F_y (N): ", self.F_y, "y lift: ", -self.L*np.cos(self.theta), "y drag: ", -self.D*np.sin(self.theta), "Fg: ", - m_glider*GRAVITY, "Fb: ", RHO_WATER*GRAVITY*(hull.V_hull + self.V_be ), self.current_len)
                
                
                self.kinematics()

                #iterate self.current_len
                self.current_len += be.laspeed*TIMESTEP
                self.be_pos_arr.append(39.3701*self.current_len)
                self.time_arr.append(self.time_arr[-1]+TIMESTEP)

            print("state entering: hold_be_pos_up, ", self.x_disp_arr[-1])
            self.current_state = 'hold_be_pos_up'
            

    
        elif self.current_state == 'move_down_accel':
            print("state: move down accel")
            while self.current_len > dist_to_retract:
                self.V_be = be.V_cont + (0.25*np.pi*be.id**2)*self.current_len
                self.calc_theta(self.current_len, self.V_be)
                #theta = np.arctan( ( hull.l_hull * (RHO_WATER*GRAVITY*(0.25*np.pi*be.id**2)*self.current_len) ) / (m_glider*GRAVITY*hull.stability) )
                self.L = 0.5*RHO_WATER*hfoil_coeff*(self.x_vel_arr[-1]**2 + self.y_vel_arr[-1]**2)
                self.D = 0.5*drag_coeff*RHO_WATER*(self.x_vel_arr[-1]**2 + self.y_vel_arr[-1]**2)*(0.25*np.pi*hull.od**2)
                
                self.F_y = RHO_WATER*GRAVITY*(hull.V_hull + self.V_be ) - m_glider*GRAVITY + self.L*np.cos(self.theta) + self.D*np.sin(self.theta)
                self.F_x = self.L*np.sin(self.theta) - self.D*np.cos(self.theta)

        
                #print("F_y (N): ", self.F_y, "F_x (N): ", self.F_x, "current len (in): ", self.current_len*39.3701)
                #print("theta! ", (180/np.pi)*theta, "F_y: ", self.F_y)
                
                self.kinematics()

                #iterate self.current_len
                self.current_len -= be.laspeed*TIMESTEP
                self.be_pos_arr.append(39.3701*self.current_len)
                self.time_arr.append(self.time_arr[-1]+TIMESTEP)
                #print("LIFT: ", self.L, self.x_vel_arr[-1], self.y_vel_arr[-1]**2)

            print("state entering: hold_be_pos_down", self.x_disp_arr[-1])
            self.current_state = 'hold_be_pos_down'
            


        elif self.current_state == 'move_up_deccel':
            print("state: move up deccel ", self.current_len, self.x_disp_arr[-1])
            while self.current_len > midpoint:
                self.V_be = be.V_cont + (0.25*np.pi*be.id**2)*self.current_len
                self.calc_theta(self.current_len, self.V_be)
                #theta = np.arctan( ( hull.l_hull * (RHO_WATER*GRAVITY*(0.25*np.pi*be.id**2)*self.current_len) ) / (m_glider*GRAVITY*hull.stability) )
                self.L = 0.5*RHO_WATER*hfoil_coeff*(self.x_vel_arr[-1]**2 + self.y_vel_arr[-1]**2)
                self.D = 0.5*drag_coeff*RHO_WATER*(self.x_vel_arr[-1]**2 + self.y_vel_arr[-1]**2)*(0.25*np.pi*hull.od**2)
                
                self.F_y = RHO_WATER*GRAVITY*(hull.V_hull + self.V_be ) - m_glider*GRAVITY - self.L*np.cos(self.theta) - self.D*np.sin(self.theta)
                self.F_x = self.L*np.sin(self.theta) - self.D*np.cos(self.theta)

                #print("F_y (N): ", self.F_y, "F_x (N): ", self.F_x, "current len (in): ", self.current_len*39.3701)
                #print("theta! ", (180/np.pi)*theta, "F_y: ", self.F_y)
                
                self.kinematics()

                #iterate self.current_len
                self.current_len -= be.laspeed*TIMESTEP
                self.be_pos_arr.append(39.3701*self.current_len)
                self.time_arr.append(self.time_arr[-1]+TIMESTEP)

            print("state entering: move down accel", self.x_disp_arr[-1])
            self.current_state = 'end'

        elif self.current_state == 'end':
            return 0


    def get_state(self):
        return self.current_state
    
    def set_state(self, new_state):
        self.current_state = new_state
    

####inputs:
hfoil_coeff = 0.0007 #Area of wing * Coeff Lift
drag_coeff = 0 #estimated from hull outer geometry ss
percent_stability = 0.15 #% diam
midpoint = intometer(4)
internal_mass = 11 #kg

#geometry
#def __init__(self, id, od, cont_len, travel_len, laspeed, midpoint):
be = BuoyancyEngine(3,3.5,10.3, 8, 0.5*0.0254,midpoint)

hull_id = intometer(5.0)
hull_od = intometer(5.563)

hull_len = (internal_mass + be.mass - RHO_WATER*be.V_mid) / ( np.pi *0.25* (RHO_WATER*(hull_od**2) - RHO_PVC*(hull_od**2 - hull_id**2)) )

print("hull_len: ", 39.3701*hull_len, " (in)")

#def __init__(self,id,od,len,percent_stability):
hull = PressureHull( hull_id, hull_od, hull_len, percent_stability)

m_glider = hull.mass + be.mass +internal_mass


max_allowable_depth = -50 #m
sim_depth = 0 #m
state = 'str'


#setup trajectory parameters
start_dist_to_retract = intometer(0)
start_time_hold_down = 2
start_dist_to_extend = intometer(8)
start_time_hold_up = 6.75

nom_dist_to_retract = intometer(0)
nom_time_hold_down = 2
nom_dist_to_extend = intometer(8)
nom_time_hold_up = 2

###remember to save the final valid numbers
while(sim_depth > max_allowable_depth):
    #def __init__(self, be, hull, m_glider):
    seapup = SeagliderFSM(be, hull, m_glider)
    print("\n")
    print("AAAA STARTING SIM!!!")
    print("\n")
    ###loop through fsm

    #starting from surface with no speed
    while(seapup.get_state() != 'end'):
        #def trajectory(self, dist_to_retract, time_hold_down, dist_to_extend, time_hold_up)
        seapup.trajectory(start_dist_to_retract, start_time_hold_down, start_dist_to_extend, start_time_hold_up)

    print("\n")
    print("split")
    print("\n")

    
    #normal yo
    for i in range(6):
        seapup.set_state('move_down_accel')
        while(seapup.get_state() != 'end'):
            #def trajectory(self, dist_to_retract, time_hold_down, dist_to_extend, time_hold_up)
            seapup.trajectory(nom_dist_to_retract, nom_time_hold_down, nom_dist_to_extend, nom_time_hold_up)
    


    #note the following breaks it, am i not sensitive enough with the inputs??? its weird because it seems nom at first then some error randomly appears
    """
    print(' ')
    print("split", seapup.x_disp_arr[-1])
    print(' ')

    seapup.set_state('move_down_accel')
    while(seapup.get_state() != 'end'):
        #def trajectory(self, dist_to_retract, time_hold_down, dist_to_extend, time_hold_up)
        seapup.trajectory(nom_dist_to_retract, nom_time_hold_down, nom_dist_to_extend, nom_time_hold_up)
    """

    nom_time_hold_up += 2
    nom_time_hold_down += 2

    state = 'str'

    sim_depth = np.min(seapup.y_disp_arr)
    max_speed = np.max(seapup.x_vel_arr)
    glide_period = seapup.x_disp_arr[-1]

    break



print("allowable extension: ",be.allowable_ext*39.3701, " (in)")
print("max fwd speed: ",max_speed, " (m/s)")
print("glide period: ",glide_period, '\n')

print("be id: ",39.3701*be.id, " (in)")
print("hull id: ",39.3701*hull.id, " (in)")
print("total glider mass: ", m_glider, " (kg)")



#data
"""
plt.subplot(1,3,1)  
plt.plot(seapup.x_disp_arr, seapup.y_disp_arr, label='Seaglider Position', color='blue')
plt.xlabel('X-Pos (m)')
plt.ylabel('Y-Pos (m)')

"""
for i in seapup.s_x_change_down_arr:
    plt.axvline(x=i, color='r', linestyle='-', label='change up')
for i in seapup.s_x_change_up_arr:
    plt.axvline(x=i, color='g', linestyle='-', label='change up')
"""

plt.subplot(1,3,2)
plt.plot(seapup.x_disp_arr, seapup.be_pos_arr, label='BE Position Over X-Pos', color='blue')
plt.xlabel('X-Pos (m)')
plt.ylabel('B.E. Pos (in)')



plt.subplot(1,3,3)
plt.plot(seapup.time_arr, seapup.be_pos_arr, label='BE Position Over X-Pos', color='blue')
plt.xlabel('Time (s)')
plt.ylabel('B.E. Pos (in)')


plt.show()
"""

#uncomment for big data

plt.subplot(2,6,1)  
plt.plot(seapup.x_disp_arr, seapup.y_disp_arr, label='Seaglider Position', color='blue')
plt.xlabel('X-Pos (m)')
plt.ylabel('Y-Pos (m)')

"""
for i in seapup.s_x_change_down_arr:
    plt.axvline(x=i, color='r', linestyle='-', label='change up')
for i in seapup.s_x_change_up_arr:
    plt.axvline(x=i, color='g', linestyle='-', label='change up')
"""

plt.subplot(2,6,2)
plt.plot(seapup.time_arr, seapup.x_disp_arr, label='Time Vs X-Pos', color='blue')
plt.xlabel('Time (s)')
plt.ylabel('X-Pos (m)')

plt.subplot(2,6,3)
plt.plot(seapup.time_arr, seapup.y_disp_arr, label='Time Vs Y-Pos', color='blue')
plt.xlabel('Time (s)')
plt.ylabel('Y-Pos (m)')

plt.subplot(2,6,4)
plt.plot(seapup.time_arr, seapup.x_vel_arr, label='Time Vs X-Vel', color='blue')
plt.xlabel('Time (s)')
plt.ylabel('X-Vel (m/s)')

plt.subplot(2,6,5)
plt.plot(seapup.time_arr, seapup.y_vel_arr, label='Time Vs Y-Vel', color='blue')
plt.xlabel('Time (s)')
plt.ylabel('Y-Vel (m/s)')

plt.subplot(2,6,6)
plt.plot(seapup.time_arr, seapup.x_accel_arr, label='Time Vs X-Accel', color='blue')
plt.xlabel('Time (s)')
plt.ylabel('X-Accel (m/s^2)')

plt.subplot(2,6,7)
plt.plot(seapup.time_arr, seapup.y_accel_arr, label='Time Vs Y-Accel', color='blue')
plt.xlabel('Time (s)')
plt.ylabel('Y-Accel (m/s^2)')

plt.subplot(2,6,8)
plt.plot(seapup.be_pos_arr, seapup.x_disp_arr, label='BE Position Vs X-Pos', color='blue')
plt.xlabel('B.E. Pos (in)')
plt.ylabel('X-Pos (m)')

plt.subplot(2,6,9)
plt.plot(seapup.be_pos_arr, seapup.x_accel_arr, label='BE Position Vs x-Accel', color='blue')
plt.xlabel('B.E. Pos (in)')
plt.ylabel('X-accel (m/s^2)')

plt.subplot(2,6,10)
plt.plot(seapup.be_pos_arr, seapup.y_vel_arr, label='BE Position Vs Y-Vel', color='blue')
plt.xlabel('B.E. Pos (in)')
plt.ylabel('Y-Vel (m/s)')

plt.subplot(2,6,11)
plt.plot(seapup.be_pos_arr, seapup.y_accel_arr, label='BE Position Vs Y-accel', color='blue')
plt.xlabel('B.E. Pos (in)')
plt.ylabel('Y-Accel (m/s^2)')

plt.subplot(2,6,12)
plt.plot(seapup.time_arr, seapup.be_pos_arr, label='Time Vs BE Position', color='blue')
plt.xlabel('Time (S)')
plt.ylabel('B.E. Pos (in)')


plt.show()
