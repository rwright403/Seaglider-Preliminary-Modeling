import numpy as np
import matplotlib.pyplot as plt

#Constants
RHO_PVC = 1380 #kg/m^3
RHO_WATER =997 #kg/m^3
GRAVITY = 9.81 #m/s^2
TIMESTEP = 0.5 #s

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
    def __init__(self):
        
        self.current_len = midpoint

        self.s_y_prev = 0 #m
        self.s_y = 0 #m
        self.s_x_prev = 0 #m
        self.s_x = 0 #m

        self.v = 1 #m/s
        self.v_x_arr = []
        self.max_speed = 0

        self.s_x_arr = []
        self.s_y_arr = []
        self.s_x_change_down_arr = []
        self.s_x_change_up_arr = []

        self.time_arr = []
        self.be_pos_arr = []
        #--> going down (L+)
        self.L = 0.5*RHO_WATER*hfoil_coeff*self.v**2

        self.F_y = 0
        self.F_x = 0

        # Define the states
        self.states = ['hold_be_pos', 'move_down_deccel', 'move_up_accel', 'move_down_accel', 'move_up_deccel']
        # Define the initial state
        self.current_state = 'move_down_accel'







    def transition(self):
        if self.current_state == 'hold_be_pos':
            counter = 0
            while counter <= holding_time:
                #calculate forces
                self.v_y = (self.F_y/m_glider)*TIMESTEP
                self.v_x = (self.F_x/m_glider)*TIMESTEP

                self.v_x_arr.append(self.v_x)

                self.s_y = self.s_y_prev + self.v_y
                self.s_x = self.s_x_prev + self.v_x

                self.s_x_arr.append(self.s_x)
                self.s_y_arr.append(self.s_y)

                self.s_y_prev = self.s_y
                self.s_x_prev = self.s_x
                #iterate counter
                counter += TIMESTEP
            
            if(self.F_y < 0):
                self.current_state = 'move_down_deccel'
            else:
                self.current_state = 'move_up_deccel'

    
        elif self.current_state == 'move_down_deccel':
            while self.current_len < midpoint:
                V_be = be.V_cont + (0.25*np.pi*be.id**2)*self.current_len
                theta = np.arctan( ( hull.l_hull * (RHO_WATER*GRAVITY*(0.25*np.pi*be.id**2)*self.current_len) ) / (m_glider*GRAVITY*hull.stability) )
                self.F_y = RHO_WATER*GRAVITY*(hull.V_hull + V_be ) - m_glider*GRAVITY + self.L*np.sin(theta)
                self.F_x = self.L*np.cos(theta)
                
                self.v_y = (self.F_y/m_glider)*TIMESTEP
                self.v_x = (self.F_x/m_glider)*TIMESTEP

                self.v_x_arr.append(self.v_x)

                self.s_y = self.s_y_prev + self.v_y
                self.s_x = self.s_x_prev + self.v_x

                self.s_x_arr.append(self.s_x)
                self.s_y_arr.append(self.s_y)

                self.s_y_prev = self.s_y
                self.s_x_prev = self.s_x

                #iterate current_len
                self.current_len += be.laspeed*TIMESTEP

            self.current_state = 'move_up_accel'
    

    
        elif self.current_state == 'move_up_accel':
            while self.current_len < be.travel_len:
                V_be = be.V_cont + (0.25*np.pi*be.id**2)*self.current_len
                theta = np.arctan( ( hull.l_hull * (RHO_WATER*GRAVITY*(0.25*np.pi*be.id**2)*self.current_len) ) / (m_glider*GRAVITY*hull.stability) )
                self.F_y = RHO_WATER*GRAVITY*(hull.V_hull + V_be ) - m_glider*GRAVITY - self.L*np.sin(theta)
                self.F_x = self.L*np.cos(theta)
                
                self.v_y = (self.F_y/m_glider)*TIMESTEP
                self.v_x = (self.F_x/m_glider)*TIMESTEP

                self.v_x_arr.append(self.v_x)

                self.s_y = self.s_y_prev + self.v_y
                self.s_x = self.s_x_prev + self.v_x

                self.s_x_arr.append(self.s_x)
                self.s_y_arr.append(self.s_y)

                self.s_y_prev = self.s_y
                self.s_x_prev = self.s_x

                #iterate self.current_len
                self.current_len += be.laspeed*TIMESTEP

            self.current_state = 'hold_be_pos'

    


        elif self.current_state == 'move_down_accel':
            while self.current_len > 0:
                V_be = be.V_cont + (0.25*np.pi*be.id**2)*self.current_len
                theta = np.arctan( ( hull.l_hull * (RHO_WATER*GRAVITY*(0.25*np.pi*be.id**2)*self.current_len) ) / (m_glider*GRAVITY*hull.stability) )
                self.F_y = RHO_WATER*GRAVITY*(hull.V_hull + V_be ) - m_glider*GRAVITY + self.L*np.sin(theta)
                self.F_x = self.L*np.cos(theta)
                
                self.v_y = (self.F_y/m_glider)*TIMESTEP
                self.v_x = (self.F_x/m_glider)*TIMESTEP

                self.v_x_arr.append(self.v_x)

                self.s_y = self.s_y_prev + self.v_y
                self.s_x = self.s_x_prev + self.v_x

                self.s_x_arr.append(self.s_x)
                self.s_y_arr.append(self.s_y)

                self.s_y_prev = self.s_y
                self.s_x_prev = self.s_x

                #iterate self.current_len
                self.current_len -= be.laspeed*TIMESTEP

            self.current_state = 'hold_be_pos'
            


        elif self.current_state == 'move_up_deccel':
            while self.current_len < midpoint:
                V_be = be.V_cont + (0.25*np.pi*be.id**2)*self.current_len
                theta = np.arctan( ( hull.l_hull * (RHO_WATER*GRAVITY*(0.25*np.pi*be.id**2)*self.current_len) ) / (m_glider*GRAVITY*hull.stability) )
                self.F_y = RHO_WATER*GRAVITY*(hull.V_hull + V_be ) - m_glider*GRAVITY - self.L*np.sin(theta)
                self.F_x = self.L*np.cos(theta)
                
                self.v_y = (self.F_y/m_glider)*TIMESTEP
                self.v_x = (self.F_x/m_glider)*TIMESTEP

                self.v_x_arr.append(self.v_x)

                self.s_y = self.s_y_prev + self.v_y
                self.s_x = self.s_x_prev + self.v_x

                self.s_x_arr.append(self.s_x)
                self.s_y_arr.append(self.s_y)

                self.s_y_prev = self.s_y
                self.s_x_prev = self.s_x

                #iterate self.current_len
                self.current_len -= be.laspeed*TIMESTEP

            self.current_state = 'move_down_deccel'


    def get_state(self):
        return self.current_state
    

####inputs:
hfoil_coeff = 0.008 #Area of wing * Coeff Lift
percent_stability = 0.15 #%
midpoint = intometer(4)
internal_mass = 11 #kg

#geometry
be = BuoyancyEngine(4.5,5.515,10.3, 8, 0.08*0.0254,midpoint)

hull_id = intometer(5.0)
hull_od = intometer(5.5)

hull_len = (internal_mass + be.mass - RHO_WATER*be.V_mid) / ( np.pi *0.25* (RHO_WATER*(hull_od**2) - RHO_PVC*(hull_od**2 - hull_id**2)) )

print("hull_len: ", 39.3701*hull_len, " (in)")

hull = PressureHull( hull_id, hull_od, hull_len, percent_stability)

m_glider = hull.mass + be.mass +internal_mass


# Create an instance of the TrafficLightFSM
seapup = SeagliderFSM()


max_allowable_depth = -20 #m
sim_depth = 0 #m
holding_time = 0
###remember to save the final valid numbers
while(sim_depth > max_allowable_depth):
    
    ###loop through fsm
    holding_time += 5
    seapup.transition()
    sim_depth = min(seapup.s_y_arr)





plt.subplot(1,2,1)  
plt.plot(seapup.s_x_arr, seapup.s_y_arr, label='Seaglider Position', color='blue')
for i in seapup.s_x_change_down_arr:
    plt.axvline(x=i, color='r', linestyle='--', label='change up')
for i in seapup.s_x_change_up_arr:
    plt.axvline(x=i, color='g', linestyle='--', label='change up')

plt.xlabel('X-Pos (m)')
plt.ylabel('Y-Pos (m)')

plt.subplot(1,2,2)
plt.plot(seapup.s_x_arr, seapup.s_y_arr, label='BE Position Over X-Pos', color='blue')
plt.xlabel('X-Pos (m)')
plt.ylabel('B.E. Pos (in)')
plt.show()