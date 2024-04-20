#Ryan Wright, UVEEC Seaglider Mech Lead 2023-2024
import numpy as np
import matplotlib.pyplot as plt

#Constants
RHO_PVC = 1380 #kg/m^3
RHO_WATER =997 #kg/m^3
GRAVITY = 9.81 #m/s^2
TIMESTEP = 0.1 #s

AOA_C_LIFT_C_DRAG_COEFFS = [(0,0,0),(0.5,0.01,0.01),(15,0.9,0.15),(25,1.75,0.4),(35,1.4,1.0),(50,1.1,1.5),(65,0.75,1.6),(75,0.45,1.75),(90,0.05,1.9)]

AOA_LOOKUP, C_LIFT_LOOKUP, C_DRAG_LOOKUP = zip(*AOA_C_LIFT_C_DRAG_COEFFS)



"""
fucntion: intometer: converts a dimension from inches to meters
input: a dimension in inches
output: the dimension converted to meters
"""
def intometer(x):
    return 0.0254*x



"""
function: linear_interpolation: 
input: 
output: the interpolated value for hull coeff of drag
"""
def linear_interpolation(x,lookup_arr):
    return (np.pi/180)*np.interp(np.abs((180/np.pi)*x), AOA_LOOKUP, lookup_arr)





"""
class: Hydrofoil
Purpose: this class stores all the information about the hydrofoil geometry for force balance calculations
inputs: the aspect ratio of the wing and the span (long dimension perpendicular to gliders main axis)
atributes: 
Aspect ratio, span (m) and wing area (m^2)
"""
class Hydrofoil:
    def __init__(self,AR,span):
        self.AR = AR
        self.span = intometer(span) #m
        
        self.width = self.span/self.AR #m

        self.Area = self.span*self.width #m^2
        print("AREA: ", self.Area)



"""
class: BuoyancyEngine
Purpose: this class stores all the information about the buoyancy engine geometry and mass for force balance calculations

inputs: inner B.E. diameter, outer B.E. diameter, fully contracted length, length of travel, speed of actuator, midpoint distance from contracted point
(inputs passed in in imperial)

attributes: 
all inputs 
volume at the contracted length and volume when B.E. is extended to midpoint
mass of expansion fitting (just pvc not considering linear actuator mass in this num)
"""
class BuoyancyEngine:
    def __init__(self, id, od, cont_len, travel_len, laspeed, midpoint):
        self.id = intometer(id) #m
        self.od = intometer(od) #m
        self.cont_len = intometer(cont_len) #m
        self.travel_len = intometer(travel_len) #m
        self.midpoint = midpoint #already passed in in metric!
        self.laspeed = laspeed #m/s

        print("init be: midpoint, travel len", self.midpoint, self.travel_len)

        #displaced volumes
        self.V_cont = 0.25*(np.pi*self.od**2)*self.cont_len #m^3
        self.V_mid = 0.25*np.pi*(self.id**2)*self.midpoint + self.V_cont #m^3

        self.mass = RHO_PVC*0.25*np.pi*(self.od**2-self.id**2)*self.cont_len #kg

        """debug print statements"""
        #print("be mass: ", self.mass)



"""
class: PressureHull
Purpose: this class stores all the information about the hull geometry and mass for force balance calculations

inputs: P.H. id, od, length. Also percent stability, which is used to calculate the distance between the cb and cg based on hull diameter
(all dimensions passed in in metric because metric used to calculate the hull length to get a neutrally buoyant glider for a given mass)

attributes: 
all inputs, 
calculated glider stability using product of % stability and hull diam
geometry like cross sectional area and volume
mass of pressure hull (pvc pipe mass)
"""
class PressureHull:
    def __init__(self,id,od,len,percent_stability):
        self.id = id #m
        self.od = od #m
        self.l_hull = len #m
        self.percent_stability = percent_stability #%
        
        self.stability = percent_stability*self.od #m
        self.area = 0.25*np.pi*self.od**2 #m
        self.surface_area = (2*self.area + self.od*np.pi*len)*0.3333 #TODO: ADD A PROPER SURFACE AREA CALC
        self.V_hull = self.area*self.l_hull #m

        self.mass = RHO_PVC*(np.pi/4*(self.od**2-self.id**2)*self.l_hull) #kg

        """debug print statements"""
        #print(self.od)
        #print("hull mass: ", self.mass)



"""
class: Seaglider FSM
Purpose: This is the class that calculates the seaglider trajectory. There are six different states:
['hold_be_pos_up', 'hold_be_pos_down', 'move_down_deccel', 'move_up_accel', 'move_down_accel', 'move_up_deccel', 'end']

inputs: buoyancy engine object, pressure hull object, glider mass, hydrofoil object

attributes:
the buoyancy engine, hydrofoil, pressure hull objects as well as the glider's internal mass
arrays to store accel, velo, disp

methods:
__init__
calc_theta
calc_coeff_lift
calc_hydro_force
kinematics
trajectory
get_state
set_state
"""
class SeagliderFSM:

    """
    __init__
    Purpose: 
    creates attributes ^
    setup storage arrays used for numpy plotting
    calculate the location of the center of buoyancy at b.e. midpoint (also cg)
    set the current state to 'move_down_deccel'
    """
    def __init__(self, be, hull, m_glider, hydrofoil):
        self.be = be
        self.hull = hull
        self.m_glider = m_glider
        self.hydrofoil = hydrofoil

        self.current_len = midpoint

        self.i = 0

        self.x_accel_arr = [0]
        self.y_accel_arr = [0]

        self.x_vel_arr = [0.0001]
        self.y_vel_arr = [-0.001]
        
        self.x_disp_arr = [0]
        self.y_disp_arr = [0]

        self.time_arr = [0]

        self.max_speed = 0

        self.s_x_change_down_arr = []
        self.s_x_change_up_arr = []        
        self.be_pos_arr = [4]

        #setup current_len
        self.current_len = midpoint
        self.V_be = 0

        self.theta = 0
        self.phi = -np.arctan(self.y_vel_arr[-1]/self.x_vel_arr[-1])
        self.aoa = 0

        self.F_y = 0
        self.F_x = 0


        self.aoa_arr = [0]
        self.phi_arr = [0]
        self.theta_arr = [0]
        self.delta_X_cb_arr = [0]

        #solve cg location: we are assuming cg is placed such that glider is horizontal at midpoint b.e.
        self.X_cb = ( (hull.V_hull*0.5*hull.l_hull) + (be.V_cont*(hull.l_hull+0.5*be.cont_len)) + ((be.V_mid-be.V_cont)*(hull.l_hull+be.cont_len+0.5*be.midpoint)) ) / (hull.V_hull + be.V_mid)
        self.X_cg = self.X_cb 

        #states:['hold_be_pos_up', 'hold_be_pos_down', 'move_down_deccel', 'move_up_accel', 'move_down_accel', 'move_up_deccel', 'end']
        # Define the initial state
        self.current_state = 'move_down_accel'

    """
    method: calc theta
    Purpose: calculate the glide angle of the glider. This is determined by the position of cb relative to cg.(new cb (changes as b.e. volume changes))
    Inputs: buoyancy engine length (measured abs length from contraction point), the current volume of the buoyancy engine
    
    Returns: nothing
    Updates: self.theta --> unsigned (positive) value of the glide angle
    """
    #TODO: GET RID OF VBE AS AN INPUT SINCE IT IS STORED
    def calc_theta(self, current_len, V_be):
        self.X_cb = ( (hull.V_hull*0.5*hull.l_hull) + (be.V_cont*(hull.l_hull+0.5*be.cont_len)) + ((self.V_be-be.V_cont)*(hull.l_hull+be.cont_len+0.5*current_len)) ) / (hull.V_hull + self.V_be)       
        delta_X_cb = self.X_cb - self.X_cg

        self.delta_X_cb_arr.append(delta_X_cb)

        self.theta = np.arccos(hull.stability/np.sqrt(delta_X_cb**2 +hull.stability**2))

        if(delta_X_cb >=0):
            self.theta = np.arccos(hull.stability/np.sqrt(delta_X_cb**2 +hull.stability**2))
        else:
            self.theta = -np.arccos(hull.stability/np.sqrt(delta_X_cb**2 +hull.stability**2))

        self.theta_arr.append((180/np.pi)*self.theta)

        """debug print statements"""
        #print("theta: ", self.theta)
        #print("theta: ",(180/np.pi)*self.theta," phi: ",(180/np.pi)*self.phi)

    """
    method: calc coeff lift
    Purpose: uses the empirical korvin theory for flat plate airfoils to calculate the wing lift coefficent (see hydrofoil project charter)
    Inputs: hydrofoil object ()

    Returns: nothing
    Updates: self.lift_coeff to access aspect ratio
    
    def calc_coeff_lift(self, hydrofoil):
        self.lift_coeff = 0.73 * (np.pi*np.sin(self.aoa))/(1+2/hydrofoil.AR)+ (2*np.pi*np.sin(self.aoa)**2)/(np.pi+4)

        #debug print statements
        #print(self.lift_coeff)


    method: calc hydro force
    Purpose: this calculates the lift and drag on the glider
    calls calc theta to update the glider's glide angle, then uses glide angle in lift coeff eqn
    
    Returns: nothing
    Updates: lift --> self.L    drag --> self.D ***NOTE: these are absolute valued magnitudes, and signed and made into components in the respective states
    (calls calc_theta and calc_coeff_lift which then update self.X_cb, self.theta and self.lift_coeff respectfully)
    """
    def calc_hydro_force(self, current_len, V_be, hydrofoil):
        self.calc_theta(current_len,V_be)

        velocity_magnitude = (self.x_vel_arr[-1]**2 + self.y_vel_arr[-1]**2)**0.5

        self.aoa = self.phi-self.theta

        #TODO: PROJECTED AREA

        self.L = 0.5*RHO_WATER*linear_interpolation(self.aoa,C_LIFT_LOOKUP)*hydrofoil.Area*velocity_magnitude**2
        self.D = 0.5*RHO_WATER*linear_interpolation(self.aoa,C_DRAG_LOOKUP)*hull.surface_area*velocity_magnitude**2

        """debug print statements"""
        #print("theta: ", self.theta, "lift coeff: ", self.lift_coeff, "lift: ", self.L, "velo: ", (self.x_vel_arr[-1]**2 + self.y_vel_arr[-1]**2)**0.5, "net y force: ", self.F_y, "be pos: ", 39.3701*self.current_len)
        #print("aoa: ",(180/np.pi)*self.aoa, "theta: ",(180/np.pi)*self.theta, "phi: ",(180/np.pi)*self.phi, "lift coeff: ", linear_interpolation(self.aoa,C_LIFT_LOOKUP), "drag coeff: ",linear_interpolation(self.aoa,C_DRAG_LOOKUP))


    """
    method: kinematics
    Purpose: uses the net force to solve acceleration, then the trapezoidal method to approx integrals for velo and displ.
    
    Returns: nothing
    Updates: x,y accel, velo, displ couter i, which is used for appending
    Appends: ^ values to plotting arrays.
    """
    def kinematics(self):
        
        self.F_y = RHO_WATER*GRAVITY*(hull.V_hull + self.V_be ) - m_glider*GRAVITY + self.L*np.cos(np.abs(self.phi)) + self.D*np.sin(np.abs(self.phi))
        self.F_x = self.L*np.sin(np.abs(self.phi)) - self.D*np.cos(np.abs(self.phi))
        #print("Fx: ", self.F_x,"x lift: ", self.L*np.sin(np.abs(self.phi)), "x drag: ", self.D*np.cos(np.abs(self.phi)))

        if(self.aoa < 0):
            self.F_y = RHO_WATER*GRAVITY*(hull.V_hull + self.V_be ) - m_glider*GRAVITY + self.L*np.cos(np.abs(self.phi)) + self.D*np.sin(np.abs(self.phi))
            #print("Fy: ", f"{self.F_y:.5f}", "Fb: ", f"{RHO_WATER*GRAVITY*(hull.V_hull + self.V_be):.5f}", "Fg: ",f"{-m_glider*GRAVITY:.5f}", "y lift: ", f"{self.L*np.cos(np.abs(self.phi)):.5f}", "y drag", f"{self.D*np.sin(np.abs(self.phi)):.5f}")
        else:
            self.F_y = RHO_WATER*GRAVITY*(hull.V_hull + self.V_be ) - m_glider*GRAVITY - self.L*np.cos(np.abs(self.phi)) - self.D*np.sin(np.abs(self.phi))
            #print("Fy: ", f"{self.F_y:.5f}", "Fb: ", f"{RHO_WATER*GRAVITY*(hull.V_hull + self.V_be):.5f}", "Fg: ",f"{-m_glider*GRAVITY:.5f}", "y lift: ", f"{-self.L*np.cos(np.abs(self.phi)):.5f}", "y drag", f"{-self.D*np.sin(np.abs(self.phi)):.5f}")
        
        #solve acceleration
        self.y_accel_arr.append(self.F_y/m_glider)
        self.x_accel_arr.append(self.F_x/m_glider)

        #solve velocity with trapezoidal method to approx integral
        self.y_vel_arr.append(self.y_vel_arr[-1] + np.trapz(y=self.y_accel_arr[-2:],dx=TIMESTEP))
        self.x_vel_arr.append(self.x_vel_arr[-1] + np.trapz(y=self.x_accel_arr[-2:],dx=TIMESTEP))

        #solve displacement with trapezoidal method to approx integral
        self.y_disp_arr.append(self.y_disp_arr[-1] + np.trapz(y=self.y_vel_arr[-2:],dx=TIMESTEP))
        self.x_disp_arr.append(self.x_disp_arr[-1] + np.trapz(y=self.x_vel_arr[-2:],dx=TIMESTEP))

        #iterate
        self.i+=1 


        self.phi = np.arctan(self.y_vel_arr[-1]/self.x_vel_arr[-1])

        self.phi_arr.append((180/np.pi)*self.phi)
        #now we have solved new velocity, so we can solve new phi

        

        #BUG: figure out transient

        

        self.aoa_arr.append((180/np.pi)*self.aoa)
        


        """debug print statements"""
        #print(self.F_x/m_glider, self.y_accel_arr[-1])
        #print(self.y_vel_arr[-1])
        #print(self.y_disp_arr[-1])
        #print(self.theta)
        #"y accel: ", self.y_accel_arr[-1],
        #print(self.x_disp_arr[-1],self.F_y)
        #print( (self.x_vel_arr[-1]**2+self.y_vel_arr[-1]**2)**0.5)
        #print("x accel: ", self.x_accel_arr[-1], "y accel: ", self.y_accel_arr[-1], "x vel: ", self.x_vel_arr[-1],"y vel: ", self.y_vel_arr[-1], "x disp: ", self.x_disp_arr[-1])
        #print(self.F_y, " | ",self.y_accel_arr[-2],self.y_accel_arr[-1])
        #print("prev: ", self.y_vel_arr[-1], "trap: ", np.trapz(y=self.y_accel_arr[-2:],dx=TIMESTEP))
        #print(self.F_y, " | ",self.y_accel_arr[-2],self.y_accel_arr[-1])
        #print("velo: ", (self.x_vel_arr[-1]**2 + self.y_vel_arr[-1]**2)**0.5, "x velo: ", self.x_vel_arr[-1], "y velo: ",self.y_vel_arr[-1], "THETA: ", (180/np.pi)*self.theta )
        #print("\n")
        #print("THETA: ", f"{(180/np.pi)*self.theta:.5f}", "PHI: ", f"{(180/np.pi)*self.phi:.5f}", "AOA: ", f"{(180/np.pi)*self.aoa:.5f}")
        #print( "v_x: ", self.x_vel_arr[-1], "v_y: ", self.y_vel_arr[-1])
        #print("i: ", self.i)
        #print("\n")



    """
    method: trajectory
    Purpose: calculate the 
    
    Returns: nothing
    Updates:
    """
    def trajectory(self, dist_to_retract, time_hold_down, dist_to_extend, time_hold_up):
        
        self.dist_to_retract = dist_to_retract
        self.time_hold_down = time_hold_down
        self.dist_to_extend = dist_to_extend
        self.time_hold_up = time_hold_up
        print(time_hold_up)

        print("start trajectory", self.current_state)





        """
        state: move_down_accel
        previous state: move_up_deccel   next state: hold_be_pos_down

        Returns: nothing
        Updates: net forces F_y, F_x, buoyancy engine position (current len)
        Appends: buoyancy engine position and time
        """
        if self.current_state == 'move_down_accel':
            print("state: move down accel")
            while self.current_len > dist_to_retract:

                self.V_be = be.V_cont + (0.25*np.pi*be.id**2)*self.current_len
                self.calc_hydro_force(self.current_len, self.V_be, self.hydrofoil)
                """
                self.F_y = RHO_WATER*GRAVITY*(hull.V_hull + self.V_be ) - m_glider*GRAVITY + self.L*np.cos(np.abs(self.phi)) + self.D*np.sin(np.abs(self.phi))
                self.F_x = self.L*np.sin(np.abs(self.phi)) - self.D*np.cos(np.abs(self.phi))
                """
                #print("Fy: ", f"{self.F_y:.5f}", "Fb: ", f"{RHO_WATER*GRAVITY*(hull.V_hull + self.V_be):.5f}", "Fg: ",f"{-m_glider*GRAVITY:.5f}", "y lift: ", f"{self.L*np.cos(np.abs(self.phi)):.5f}", "y drag", f"{self.D*np.sin(np.abs(self.phi)):.5f}")
                #print("Fx: ", self.F_x,"x lift: ", self.L*np.sin(np.abs(self.phi)), "x drag: ", self.D*np.cos(np.abs(self.phi)))

                self.kinematics()

                #iterate self.current_len
                self.current_len -= be.laspeed*TIMESTEP
                self.be_pos_arr.append(39.3701*self.current_len)
                self.time_arr.append(self.time_arr[-1]+TIMESTEP)

                """debug print statements"""
                
                #print("x lift: ", self.L*np.sin(self.theta), "x drag: ", - self.D*np.cos(self.theta))
                #print("F_y (N): ", self.F_y, "F_x (N): ", self.F_x, "current len (in): ", self.current_len*39.3701)
                #print("theta! ", (180/np.pi)*theta, "F_y: ", self.F_y)
                #print("LIFT: ", self.L, self.x_vel_arr[-1], self.y_vel_arr[-1]**2)
                #print("theta! ", (180/np.pi)*self.theta, "F_y: ", self.F_y)
            print("state entering: hold_be_pos_down", self.x_disp_arr[-1])
            self.current_state = 'hold_be_pos_down'
            

            """
            state: hold_be_pos_up
            previous state: move_up_accel    next state: move_up_deccel
            Purpose: calculate the force balance for when the buoyancy engine is expanding past the midpoint
        
            Returns: nothing
            Updates: net forces F_y, F_x
            Appends: buoyancy engine position and time
            """
        elif self.current_state == 'hold_be_pos_up':

            print("state: hold be pos", "F_y: ", self.F_y, "current_len: ", self.current_len)
            print("holding time UP!!! (s): ", self.time_hold_up)
            
            current_time = 0
            print("holding time up (s): ", self.time_hold_up)

            while current_time < self.time_hold_up:

                self.calc_hydro_force(self.current_len, self.V_be, self.hydrofoil)
                self.kinematics()

                #iterate/move forward current_time
                current_time += TIMESTEP

                self.be_pos_arr.append(39.3701*self.current_len)
                self.time_arr.append(self.time_arr[-1]+TIMESTEP)

                """debug print statements"""
                #print("Fy: ", self.F_y, "Fb: ",RHO_WATER*GRAVITY*(hull.V_hull + self.V_be), "Fg: ",-m_glider*GRAVITY, "y lift: ", -self.L*np.cos(self.theta), "y drag", -self.D*np.sin(self.theta))
                #print("theta! ", (180/np.pi)*self.theta, "F_y: ", self.F_y)
                #print("x lift: ", self.L*np.sin(self.theta), "x drag: ", - self.D*np.cos(self.theta), "velo: ", (self.x_vel_arr[-1]**2 + self.y_vel_arr[-1]**2)**0.5, "theta: ", (180/np.pi)*self.theta)
                #print("F_y (N): ", self.F_y, "F_x (N): ", self.F_x, "current len (in): ", self.current_len*39.3701)
                #print("theta! ", (180/np.pi)*self.theta, "F_y: ", self.F_y)
        
            print("state entering: move up deccel", self.x_accel_arr[-1], self.F_y)
            self.current_state = 'move_up_deccel'



            """
            state: hold_be_pos_down
            previous state: move_down_accel    next state: move_down_deccel
        
            Returns: nothing
            Updates: net forces F_y, F_x
            Appends: buoyancy engine position and time
            """
        elif self.current_state == 'hold_be_pos_down':
            print("state: hold be pos", "F_y: ", self.F_y, "current_len: ", self.current_len)
            current_time = 0
            print("holding time down (s): ", self.time_hold_down)

            while current_time < self.time_hold_down:

                self.calc_hydro_force(self.current_len, self.V_be, self.hydrofoil)
                self.kinematics()

                #iterate/move forward current_time
                current_time += TIMESTEP
                self.be_pos_arr.append(39.3701*self.current_len)
                self.time_arr.append(self.time_arr[-1]+TIMESTEP)

                """debug print statements"""
                #print("F_y (N): ", self.F_y, "F_x (N): ", self.F_x, "current len (in): ", self.current_len*39.3701)
                #print("Fy: ", self.F_y, "Fb: ",RHO_WATER*GRAVITY*(hull.V_hull + self.V_be), "Fg: ",-m_glider*GRAVITY, "y lift: ", +self.L*np.cos(self.theta), "y drag", self.D*np.sin(self.theta))
                #print("x lift: ", self.L*np.sin(self.theta), "x drag: ", - self.D*np.cos(self.theta))
                #print("theta! ", (180/np.pi)*self.theta, "F_y: ", self.F_y)

            print("state entering: move down deccel", self.x_accel_arr[-1], self.F_y)
            self.current_state = 'move_down_deccel'



            
            """
            state: move_down_deccel
            previous state: hold_be_pos_up     next state: move_up_accel

            Returns: nothing
            Updates: net forces F_y, F_x, buoyancy engine position (current len)
            Appends: buoyancy engine position and time
            """
        elif self.current_state == 'move_down_deccel':
            print("state: move down deccel")
            while self.current_len < midpoint:

                self.V_be = be.V_cont + (0.25*np.pi*be.id**2)*self.current_len
                self.calc_hydro_force(self.current_len, self.V_be, self.hydrofoil)
                """
                self.F_y = RHO_WATER*GRAVITY*(hull.V_hull + self.V_be ) - m_glider*GRAVITY + self.L*np.cos(np.abs(self.phi)) + self.D*np.sin(np.abs(self.phi))
                self.F_x = self.L*np.sin(np.abs(self.phi)) - self.D*np.cos(np.abs(self.phi))            
                """
                self.kinematics()

                #iterate current_len
                self.current_len += be.laspeed*TIMESTEP
                self.be_pos_arr.append(39.3701*self.current_len)
                self.time_arr.append(self.time_arr[-1]+TIMESTEP)

                """debug print statements"""
                #print("F_y (N): ", self.F_y, "F_x (N): ", self.F_x, "current len (in): ", self.current_len*39.3701)
                #print("theta! ", (180/np.pi)*self.theta, "F_y: ", self.F_y, "y lift: ", +self.L*np.cos(self.theta), "y drag", +self.D*np.sin(self.theta), "velo: ", (self.x_vel_arr[-1]**2 + self.y_vel_arr[-1]**2)**0.5)
                #print("Fy: ", self.F_y, "Fb: ",RHO_WATER*GRAVITY*(hull.V_hull + self.V_be), "Fg: ",-m_glider*GRAVITY, "y lift: ", +self.L*np.cos(self.theta), "y drag", +self.D*np.sin(self.theta))
                #print("x lift: ", self.L*np.sin(self.theta), "x drag: ", - self.D*np.cos(self.theta))
                #print(self.F_y, 39.3701*self.current_len)

            print("state entering: move up accel", self.x_disp_arr[-1])
            self.current_state = 'move_up_accel'
    



            """
            state: move_up_accel
            previous state: move_down_deccel    next state: hold_be_pos_up

            Returns: nothing
            Updates: net forces F_y, F_x, buoyancy engine position (current len)
            Appends: buoyancy engine position and time
            """
        elif self.current_state == 'move_up_accel':
            print("state: move up accel")
            while self.current_len < dist_to_extend:
                
                self.V_be = be.V_cont + (0.25*np.pi*be.id**2)*self.current_len
                
                """
                self.F_y = RHO_WATER*GRAVITY*(hull.V_hull + self.V_be ) - m_glider*GRAVITY - self.L*np.cos(self.theta) - self.D*np.sin(self.theta)
                self.F_x = self.L*np.sin(self.theta) - self.D*np.cos(self.theta)
                """
                self.kinematics()

                self.calc_hydro_force(self.current_len, self.V_be, self.hydrofoil)

                #iterate self.current_len
                self.current_len += be.laspeed*TIMESTEP
                self.be_pos_arr.append(39.3701*self.current_len)
                self.time_arr.append(self.time_arr[-1]+TIMESTEP)

                """debug print statements"""
                #print("Fy: ", self.F_y, "Fb: ",RHO_WATER*GRAVITY*(hull.V_hull + self.V_be), "Fg: ",-m_glider*GRAVITY, "y lift: ", -self.L*np.cos(self.theta), "y drag", -self.D*np.sin(self.theta))
                #print("x lift: ", self.L*np.sin(self.theta), "x drag: ", - self.D*np.cos(self.theta))
                #print("F_y (N): ", self.F_y, "F_x (N): ", self.F_x, "current len (in): ", self.current_len*39.3701)
                #print("theta! ", (180/np.pi)*theta, "F_y: ", self.F_y)
                #print("F_y (N): ", self.F_y, "y lift: ", -self.L*np.cos(self.theta), "y drag: ", -self.D*np.sin(self.theta), "Fg: ", - m_glider*GRAVITY, "Fb: ", RHO_WATER*GRAVITY*(hull.V_hull + self.V_be ), "x velo: ", self.x_vel_arr[-1])
                #print("theta! ", (180/np.pi)*self.theta, "F_y: ", self.F_y)

            print("state entering: hold_be_pos_up, ", self.x_disp_arr[-1])
            self.current_state = 'hold_be_pos_up'
            print(self.current_state)
            


            



            """
            state: move_up_deccel
            previous state: hold_be_pos_up   next state: end

            Returns: nothing
            Updates: net forces F_y, F_x, buoyancy engine position (current len)
            Appends: buoyancy engine position and time
            """
        elif self.current_state == 'move_up_deccel':
            print("state: move up deccel ", self.current_len, self.x_disp_arr[-1])
            while self.current_len > midpoint:

                self.V_be = be.V_cont + (0.25*np.pi*be.id**2)*self.current_len
                self.calc_hydro_force(self.current_len, self.V_be, self.hydrofoil)
                """
                self.F_y = RHO_WATER*GRAVITY*(hull.V_hull + self.V_be ) - m_glider*GRAVITY - self.L*np.cos(self.theta) - self.D*np.sin(self.theta)
                self.F_x = self.L*np.sin(self.theta) - self.D*np.cos(self.theta)
                """
                self.kinematics()

                #iterate self.current_len
                self.current_len -= be.laspeed*TIMESTEP
                self.be_pos_arr.append(39.3701*self.current_len)
                self.time_arr.append(self.time_arr[-1]+TIMESTEP)

                """debug print statements"""
                #print("Fy: ", self.F_y, "Fb: ",RHO_WATER*GRAVITY*(hull.V_hull + self.V_be), "Fg: ",-m_glider*GRAVITY, "y lift: ", -self.L*np.cos(self.theta), "y drag", -self.D*np.sin(self.theta))
                #print("x lift: ", self.L*np.sin(self.theta), "x drag: ", - self.D*np.cos(self.theta), "theta: ", self.theta)
                #print("F_y (N): ", self.F_y, "F_x (N): ", self.F_x, "current len (in): ", self.current_len*39.3701)
                #print("theta! ", (180/np.pi)*theta, "F_y: ", self.F_y)
                #print("theta! ", (180/np.pi)*self.theta, "F_y: ", self.F_y)

            print("state entering: end", self.x_disp_arr[-1])
            self.current_state = 'end'




            """
            state: end
            previous state: move_up_deccel   no next state

            Returns: 0
            Serves to end the yo
            """
        elif self.current_state == 'end':
            return 0


    """
    method: get_state
    Purpose: get the current state of the trajectory method
    
    Returns: the current state
    Updates: nothing
    """
    def get_state(self):
        return self.current_state
    

    """
    method: set_state
    Purpose: set the current state of the trajectory method 
    
    Returns: nothing
    Updates: the current state
    """
    def set_state(self, new_state):
        self.current_state = new_state
    





"""
THIS IS THE MAIN BODY OF THE CODE
this script works under the assumption that the linear actuator in the b.e. is so fast that it can reach the endpoint before it needs to reverse to avoid crush depth
the main unknown are the speeds to hold the b.e. We solve this ITERATIVELY. Every iteration increases the time we are holding the b.e. pos until we hold it so long the glider reaches operational depth.

the first loop checks if the trajectory has reached max allowable depth
the second loop(s) run the trajectory 
"""

####inputs/setup:
hydrofoil = Hydrofoil(3, 12) # aspect ratio and span (in inches)

drag_coeff = 0 #estimated from hull outer geometry ss
percent_stability = 0.15 #% diam
midpoint = intometer(4)
internal_mass = 11 #kg

#glider part objects:
#def __init__(self, id, od, cont_len, travel_len, laspeed, midpoint):
be = BuoyancyEngine(3,3.5,10.3, 8, 0.5*0.0254,midpoint)

#setup hull
hull_id = intometer(5.0)
hull_od = intometer(5.563)
hull_len = (internal_mass + be.mass - RHO_WATER*be.V_mid) / ( np.pi *0.25* (RHO_WATER*(hull_od**2) - RHO_PVC*(hull_od**2 - hull_id**2)) ) #hull length is solved for to create a neutrally buoyant glider
print("hull_len: ", 39.3701*hull_len, " (in)")

#def __init__(self,id,od,len,percent_stability): #NOTE: passing dimensions in in meters
hull = PressureHull( hull_id, hull_od, hull_len, percent_stability)

#this is the total glider mass
m_glider = hull.mass + be.mass +internal_mass

#setting the maximum depth to simulate to
max_allowable_depth = -50 #m
sim_depth = 0 #m

#setup trajectory parameters --> either distance to go to or times to hold down
start_dist_to_retract = intometer(0) #m
start_time_hold_down = 4 #s
start_dist_to_extend = intometer(8) #m
start_time_hold_up = 4 #s


nom_dist_to_retract = intometer(0) #m
nom_time_hold_down = 4 #s
nom_dist_to_extend = intometer(8) #m
nom_time_hold_up = 1000000000 #s



"""
This loop is where the trajectory fsm is run
we 
"""
while(sim_depth > max_allowable_depth):
    
    #def __init__(self, be, hull, m_glider):
    seapup = SeagliderFSM(be, hull, m_glider, hydrofoil)

    #TODO: DELETE WHEN SCRIPT WORKS AND NO MORE DEBUGGING
    print("\n")
    print("STARTING SIM!!!")
    print("\n")

    ###loop through fsm, #starting from surface with no speed for first run, second run we already have speed and are gliding normally
    while(seapup.get_state() != 'end'):
        #def trajectory(self, dist_to_retract, time_hold_down, dist_to_extend, time_hold_up)
        seapup.trajectory(start_dist_to_retract, start_time_hold_down, start_dist_to_extend, start_time_hold_up)



    print("\n")
    print("entering second")
    print("\n")

    
    #normal yo
    
    for i in range(1):
        seapup.set_state('move_down_accel')
        while(seapup.get_state() != 'hold_be_pos_up'):
            #def trajectory(self, dist_to_retract, time_hold_down, dist_to_extend, time_hold_up)
            seapup.trajectory(nom_dist_to_retract, nom_time_hold_down, nom_dist_to_extend, nom_time_hold_up)

        #TODO: REMOVE AFTER DEBUGGING
    


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

    #iterate and reset variables
    nom_time_hold_up += 2
    nom_time_hold_down += 2

    state = 'start'

    #this is for testing purposes otherwise we get stuck in loop BUG: get stuck in loop
    break

sim_depth = np.min(seapup.y_disp_arr)
max_speed = np.max(seapup.x_vel_arr)
glide_period = seapup.x_disp_arr[-1]

#print a summary of noteworthy params
print("max fwd speed: ",max_speed, " (m/s)")
print("glide period: ",glide_period, '\n')

print("be id: ",39.3701*be.id, " (in)")
print("hull id: ",39.3701*hull.id, " (in)")
print("total glider mass: ", m_glider, " (kg)")













#big graphing:
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

plt.subplot(2,8,1)  
plt.plot(seapup.x_disp_arr, seapup.y_disp_arr, label='Seaglider Position', color='blue')
plt.xlabel('X-Pos (m)')
plt.ylabel('Y-Pos (m)')

"""
for i in seapup.s_x_change_down_arr:
    plt.axvline(x=i, color='r', linestyle='-', label='change up')
for i in seapup.s_x_change_up_arr:
    plt.axvline(x=i, color='g', linestyle='-', label='change up')
"""

plt.subplot(2,8,2)
plt.plot(seapup.x_disp_arr, seapup.time_arr, label='X-Pos Vs Time', color='blue')
plt.ylabel('Time (s)')
plt.xlabel('X-Pos (m)')

plt.subplot(2,8,3)
plt.plot(seapup.time_arr, seapup.y_disp_arr, label='Time Vs Y-Pos', color='blue')
plt.xlabel('Time (s)')
plt.ylabel('Y-Pos (m)')

plt.subplot(2,8,4)
plt.plot(seapup.time_arr, seapup.x_vel_arr, label='Time Vs X-Vel', color='blue')
plt.xlabel('Time (s)')
plt.ylabel('X-Vel (m/s)')

plt.subplot(2,8,5)
plt.plot(seapup.time_arr, seapup.y_vel_arr, label='Time Vs Y-Vel', color='blue')
plt.xlabel('Time (s)')
plt.ylabel('Y-Vel (m/s)')

plt.subplot(2,8,6)
plt.plot(seapup.time_arr, seapup.x_accel_arr, label='Time Vs X-Accel', color='blue')
plt.xlabel('Time (s)')
plt.ylabel('X-Accel (m/s^2)')

plt.subplot(2,8,7)
plt.plot(seapup.time_arr, seapup.y_accel_arr, label='Time Vs Y-Accel', color='blue')
plt.xlabel('Time (s)')
plt.ylabel('Y-Accel (m/s^2)')

plt.subplot(2,8,8)
plt.plot(seapup.x_disp_arr, seapup.be_pos_arr, label='X-Pos Vs BE Position ', color='blue')
plt.ylabel('B.E. Pos (in)')
plt.xlabel('X-Pos (m)')

plt.subplot(2,8,9)
plt.plot(seapup.be_pos_arr, seapup.x_accel_arr, label='BE Position Vs x-Accel', color='blue')
plt.xlabel('B.E. Pos (in)')
plt.ylabel('X-accel (m/s^2)')

plt.subplot(2,8,10)
plt.plot(seapup.be_pos_arr, seapup.y_vel_arr, label='BE Position Vs Y-Vel', color='blue')
plt.xlabel('B.E. Pos (in)')
plt.ylabel('Y-Vel (m/s)')

plt.subplot(2,8,11)
plt.plot(seapup.be_pos_arr, seapup.y_accel_arr, label='BE Position Vs Y-accel', color='blue')
plt.xlabel('B.E. Pos (in)')
plt.ylabel('Y-Accel (m/s^2)')

plt.subplot(2,8,12)
plt.plot(seapup.time_arr, seapup.be_pos_arr, label='Time Vs BE Position', color='blue')
plt.xlabel('Time (S)')
plt.ylabel('B.E. Pos (in)')

plt.subplot(2,8,13)
plt.plot(seapup.time_arr, seapup.aoa_arr, label='Time Vs aoa (deg)', color='blue')
plt.xlabel('Time (S)')
plt.ylabel('aoa (deg)')

plt.subplot(2,8,14)
plt.plot(seapup.time_arr, seapup.phi_arr, label='Time Vs phi (deg)', color='blue')
plt.xlabel('Time (S)')
plt.ylabel('phi (deg)')

plt.subplot(2,8,15)
plt.plot(seapup.time_arr, seapup.theta_arr, label='Time Vs theta (deg)', color='blue')
plt.xlabel('Time (S)')
plt.ylabel('theta (deg)')

plt.subplot(2,8,16)
plt.plot(seapup.time_arr, seapup.delta_X_cb_arr, label='Time Vs delta cb (m)', color='blue')
plt.xlabel('Time (S)')
plt.ylabel('delta X_cb (m)')



plt.show()

#