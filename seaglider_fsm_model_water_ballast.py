#Ryan Wright, UVEEC Seaglider Mech Lead 2023-2024
import numpy as np
import matplotlib.pyplot as plt

#Constants
RHO_PVC = 1380 #kg/m^3
RHO_WATER = 997 #kg/m^3
RHO_ALU = 2710 #kg/m^3
GRAVITY = 9.81 #m/s^2
TIMESTEP = 0.05 #s


#NOTE: 5.5" coeff
#AOA_C_LIFT_C_DRAG_COEFFS = [(0,0,0),(10,0.60129,0.16156),(20,1.10264,0.41936),(30,1.68006,0.86424),(40,2.34058,1.62498),(50,2.71196,2.56013),(60,2.61962,3.44442),(70,2.03519,4.06102),(80,1.15359,4.41248),(90,0.04466,4.43128)]

#NOTE: 6.5" coeff w 18" wetted area
#AOA_C_LIFT_C_DRAG_COEFFS = [(0,0,0.13), (10, 0.8744, 0.63502), (20, 1.31262, 1.68018), (30, 1.59276, 2.91821), (40, 1.86567, 4.92178), (50, 2.01134, 7.52057), (60, 1.80005, 9.36923), (70, 1.38445, 10.70453), (80, 0.86694, 11.76384), (90, 0.02080, 12.05325)]

#NOTE: 6.5" coeff w 24" wetted area
AOA_C_LIFT_C_DRAG_COEFFS = [(0,0,0.13), (10, 0.8744, 0.63502), (20, 1.31262, 1.68018), (30, 1.59276, 2.91821), (40, 1.86567, 4.92178), (50, 2.01134, 7.52057), (60, 1.80005, 9.36923), (70, 1.38445, 10.70453), (80, 0.86694, 11.76384), (90, 0.02080, 12.05325)]

AOA_LOOKUP, C_LIFT_LOOKUP, C_DRAG_LOOKUP = zip(*AOA_C_LIFT_C_DRAG_COEFFS)



"""
function: intometer: converts a dimension from inches to meters
input: a dimension in inches
output: the dimension converted to meters
"""
def intometer(x):
    return 0.0254*x



"""
function: linear_interpolation: 
input: 
output: the interpolated value for glider coeff of drag and lift in degrees
"""
def linear_interpolation(x,lookup_arr):
    y = ((180/np.pi)*x)
    z =  np.interp(np.abs(y), AOA_LOOKUP, lookup_arr)
    #print(z)
    return z





"""
class: Hydrofoil
Purpose: this class stores all the information about the hydrofoil geometry for force balance calculations
inputs: reference area (m^2)

attributes:
all inputs
"""
class Hydrofoil:
    def __init__(self,ref_area):
        self.ref_area = ref_area #m^2
        print("AREA: ", self.ref_area)





"""
class: BallastTank
Purpose: this class stores all the information about the pressure hull and mass for force balance calculation

inputs: tank length, od tank, id tank, mass of half tank water, distance of tank center from tail of glider
(input dimensions passed in as imperial, masses in kg)

attributes:
all inputs
tank volume
mass of emtpy tank assuming aluminum construction (calculated from volume and density)
mass of the tank for the glider to have neutral buoyancy
full mass of the tank
current mass of the tank NOTE: this is updated throughout the trajectory

method:
update_mass --> This updates the mass of the water tank and returns the new total glider mass.
"""
class BallastTank:
    def __init__(self, l_tank, od_tank, id_tank, m_half_water, x_tank):
        self.l_tank = intometer(l_tank) #m
        self.od_tank = intometer(od_tank) #m
        self.id_tank = intometer(id_tank) #m
        self.m_half_water = m_half_water #kg
        self.x_tank = intometer(x_tank) #m 

        self.tank_volume = self.l_tank*0.25*np.pi*(self.id_tank**2)  #m^3
        self.m_empty_tank = RHO_ALU*(0.25*np.pi*( (self.od_tank**2) - (self.id_tank**2) ))*self.l_tank #kg
        
        self.m_neutral_buoy = self.m_empty_tank + self.m_half_water #kg
        self.m_full = self.m_neutral_buoy + self.m_half_water #kg

        self.m_current = self.m_neutral_buoy #kg
    
    def update_mass(self, delta_water_mass, m_glider):
        self.m_current = self.m_current + delta_water_mass
        return m_glider + delta_water_mass



"""
class: TCSmicropump
Purpose: class stores information about every pump we are considering
inputs: name of pump
(input dimensions passed in as imperial, masses in kg)

attributes:
for each respective pump:
    flowrate
"""
class TCSmicropump:
    def __init__(self, name):

        if(name == "MGD3000F"):
            self.q = 0.00005833333333
        elif(name == "MGD3000S"):
            self.q = 0.00004166666667
        elif(name == "MGD2000F"):
            self.q = 0.00003833333333
        elif(name == "MGD1000F"):
            self.q = 0.00001916666667


"""
class: PressureHull
Purpose: this class stores all the information about the hull geometry and mass for force balance calculations

inputs: P.H. id, od, length, (percent stability, which is used to calculate the distance between the cb and cg based on hull diameter) and cfd ref area
(all dimensions passed in in metric because metric used to calculate the hull length to get a neutrally buoyant glider for a given mass)

attributes: 
all inputs, 
calculated glider stability using product of % stability and hull diam
geometry like cross sectional area and volume
reference area from cfd for drag coeffs
mass of pressure hull (pvc pipe mass)
"""
class PressureHull:
    def __init__(self,id,od,len,stability,ref_area):
        self.id = intometer(id) #m
        self.od = intometer(od) #m
        self.l_hull = intometer(len) #m
        self.stability = intometer(stability) #m

        self.area = 0.25*np.pi*self.od**2 #m^2
        self.ref_area = ref_area #m^2 

        self.x_hull = self.l_hull/2
        self.V_displacement = self.area*self.l_hull #m
        self.mass = RHO_PVC*(np.pi/4*(self.od**2-self.id**2)*self.l_hull) #kg

        """debug print statements"""
        #print(self.od)
        #print("hull mass: ", self.mass)

"""
class: MovingMass
Purpose: stores info about remaining mass in glider (not wate tank), can be developed into a true moving mass later (might help w control system???)

inputs: mass, distance of center of mass to glider tail
(input dimensions passed in METRIC!, masses in kg)

attributes:
all inputs
"""
class MovingMass:
    def __init__(self, m_int, x_int):
        self.m_int = m_int
        self.x_int = x_int


"""
class: Seaglider FSM
Purpose: This is the class that calculates the seaglider trajectory. There are six different states:
['hold_be_pos_up', 'hold_be_pos_down', 'move_down_deccel', 'move_up_accel', 'move_down_accel', 'move_up_deccel', 'end']

inputs: buoyancy engine object, pressure hull object, glider mass, hydrofoil object

attributes:
the water tamk,pump. moving mass, hydrofoil, pressure hull objects as well as the glider's internal mass
arrays to store accel, velo, disp (and graph)

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
    def __init__(self, tank, pump, movingmass, hull, m_glider, hydrofoil):
        self.tank = tank
        self.pump = pump
        self.movingmass = movingmass
        self.hull = hull
        self.m_glider = m_glider
        self.hydrofoil = hydrofoil

        self.delta_water_mass = TIMESTEP*pump.q*RHO_WATER

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
        self.tank_mass_arr = [self.tank.m_current]



        self.theta = 0
        self.phi = -np.arctan(self.y_vel_arr[-1]/self.x_vel_arr[-1])
        self.aoa = 0

        self.F_y = 0
        self.F_x = 0


        self.aoa_arr = [0]
        self.phi_arr = [0]
        self.theta_arr = [0]
        self.delta_X_cg_arr = [0]

        #solve cb location: we are assuming cb is in middle of hull (neglecting non pressure hull displacement)
        self.X_cb = self.hull.x_hull
        self.X_cg = self.X_cb
        
        self.sim_depth = 0

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
    def calc_theta(self):
        self.X_cg = ( self.hull.x_hull*self.hull.mass + self.movingmass.x_int*self.movingmass.m_int + self.tank.x_tank*self.tank.m_current)/self.m_glider 
        
        delta_X_cg = self.X_cb - self.X_cg

        self.delta_X_cg_arr.append(delta_X_cg)

        self.theta = np.arctan(delta_X_cg/self.hull.stability) #TODO: DOUBLE CHECK WHEN NOT TIRED BUT LIKELY MAKES SENSE

        self.theta_arr.append((180/np.pi)*self.theta)

        """debug print statements"""
        #print("theta: ", (180/np.pi)*self.theta, "delta_X_cg: ", delta_X_cg)
        #print("theta: ",(180/np.pi)*self.theta," phi: ",(180/np.pi)*self.phi)
        #print(self.y_accel_arr[-1], delta_X_cg, self.theta)

    """
    method: calc hydro force
    Purpose: this calculates the lift and drag on the glider
    calls calc theta to update the glider's glide angle, then uses glide angle in lift coeff eqn
    
    Returns: nothing
    Updates: lift --> self.L    drag --> self.D ***NOTE: these are absolute valued magnitudes, and signed and made into components in the respective states
    (calls calc_theta and calc_coeff_lift which then update self.X_cb, self.theta and self.lift_coeff respectfully)
    """
    def calc_hydro_force(self):
        self.calc_theta()

        velocity_magnitude = (self.x_vel_arr[-1]**2 + self.y_vel_arr[-1]**2)**0.5

        self.aoa = self.phi-self.theta #TODO: check this sign

        self.L = 0.5*RHO_WATER*linear_interpolation(self.aoa,C_LIFT_LOOKUP)*self.hydrofoil.ref_area*velocity_magnitude**2
        self.D = 0.5*RHO_WATER*linear_interpolation(self.aoa,C_DRAG_LOOKUP)*self.hull.ref_area*velocity_magnitude**2


        """debug print statements"""
        #print("theta: ", self.theta, "lift: ", self.L, "velo: ", (self.x_vel_arr[-1]**2 + self.y_vel_arr[-1]**2)**0.5, "net y force: ", self.F_y, "be pos: ", self.m_glider)
        #print("aoa: ",f"{(180/np.pi)*self.aoa:.2f}", "theta: ",f"{(180/np.pi)*self.theta:.2f}", "phi: ",f"{(180/np.pi)*self.phi:.2f}", "lift coeff: ", linear_interpolation(self.aoa,C_LIFT_LOOKUP), "drag coeff: ",linear_interpolation(self.aoa,C_DRAG_LOOKUP))


    """
    method: kinematics
    Purpose: uses the net force to solve acceleration, then the trapezoidal method to approx integrals for velo and displ.
    
    Returns: nothing
    Updates: x,y accel, velo, displ couter i, which is used for appending
    Appends: ^ values to plotting arrays.
    """
    def kinematics(self):


        self.F_x = self.L*np.sin(np.abs(self.phi)) - self.D*np.cos(np.abs(self.phi))
        #print("Fx: ", self.F_x,"x lift: ", self.L*np.sin(np.abs(self.phi)), "x drag: ", -self.D*np.cos(np.abs(self.phi)))

        if(self.aoa < 0):
            self.F_y = RHO_WATER*GRAVITY*(self.hull.V_displacement) - self.m_glider*GRAVITY + self.L*np.cos(np.abs(self.phi)) + self.D*np.sin(np.abs(self.phi))
            #print("Fy: ", f"{self.F_y:.5f}", "Fb: ", f"{RHO_WATER*GRAVITY*(hull.V_displacement):.5f}", "Fg: ",f"{-self.m_glider*GRAVITY:.5f}", "y lift: ", f"{self.L*np.cos(np.abs(self.phi)):.5f}", "y drag", f"{self.D*np.sin(np.abs(self.phi)):.5f}")
        else:
            self.F_y = RHO_WATER*GRAVITY*(self.hull.V_displacement) - self.m_glider*GRAVITY - self.L*np.cos(np.abs(self.phi)) - self.D*np.sin(np.abs(self.phi))
            #print("Fy: ", f"{self.F_y:.5f}", "Fb: ", f"{RHO_WATER*GRAVITY*(hull.V_displacement):.5f}", "Fg: ",f"{-self.m_glider*GRAVITY:.5f}", "y lift: ", f"{-self.L*np.cos(np.abs(self.phi)):.5f}", "y drag", f"{-self.D*np.sin(np.abs(self.phi)):.5f}")
        
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
        self.aoa_arr.append((180/np.pi)*self.aoa)
        


        """debug print statements"""
        #print(self.m_glider, self.tank.m_current)
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
    Purpose: calculate the glider trajectory
    
    Returns: nothing
    Updates:
    current tank mass
    accel, velo displacement, time and tank mass arrays
    """
    def trajectory(self,time_hold_down, time_hold_up):
        
        self.time_hold_down = time_hold_down
        self.time_hold_up = time_hold_up

        """
        state: move_down_accel
        previous state: move_up_deccel   next state: hold_be_pos_down

        Returns: nothing
        Updates: net forces F_y, F_x, tank mass
        Appends: tank mass and time
        """
        if self.current_state == 'move_down_accel':
            print("state: move down accel")
            while self.tank.m_current < self.tank.m_full:

                self.m_glider = self.tank.update_mass(self.delta_water_mass, self.m_glider)

                self.calc_hydro_force()
                self.kinematics()

                self.tank_mass_arr.append(self.tank.m_current)
                self.time_arr.append(self.time_arr[-1]+TIMESTEP)

                """debug print statements"""
                
                #print("x lift: ", self.L*np.sin(self.theta), "x drag: ", - self.D*np.cos(self.theta))
                #print("F_y (N): ", self.F_y, "F_x (N): ", self.F_x, "current len (in): ", self.current_len*39.3701)
                #print("theta! ", (180/np.pi)*theta, "F_y: ", self.F_y)
                #print("LIFT: ", self.L, self.x_vel_arr[-1], self.y_vel_arr[-1]**2)
                #print("theta! ", (180/np.pi)*self.theta, "F_y: ", self.F_y)
                #print("Fy: ", f"{self.F_y:.5f}", "Fb: ", f"{RHO_WATER*GRAVITY*(hull.V_hull + self.V_be):.5f}", "Fg: ",f"{-m_glider*GRAVITY:.5f}", "y lift: ", f"{self.L*np.cos(np.abs(self.phi)):.5f}", "y drag", f"{self.D*np.sin(np.abs(self.phi)):.5f}")
                #print("Fx: ", self.F_x,"x lift: ", self.L*np.sin(np.abs(self.phi)), "x drag: ", self.D*np.cos(np.abs(self.phi)))
            
            print("state entering: hold_be_pos_down", self.time_arr[-1])
            self.current_state = 'hold_be_pos_down'
            

            """
            state: hold_be_pos_up
            previous state: move_up_accel    next state: move_up_deccel
            Purpose: calculate the force balance for when the buoyancy engine is expanding past the midpoint
        
            Returns: nothing
            Updates: net forces F_y, F_x
            Appends: tank mass and time
            """
        elif self.current_state == 'hold_be_pos_up':
            print("holding time up (s): ", self.time_hold_up)

            current_time = 0

            while current_time < self.time_hold_up:

                self.calc_hydro_force()
                self.kinematics()

                #iterate/move forward current_time
                current_time += TIMESTEP
                self.tank_mass_arr.append(self.tank.m_current)
                self.time_arr.append(self.time_arr[-1]+TIMESTEP)

                """debug print statements"""
                #print("Fy: ", self.F_y, "Fb: ",RHO_WATER*GRAVITY*(hull.V_hull + self.V_be), "Fg: ",-m_glider*GRAVITY, "y lift: ", -self.L*np.cos(self.theta), "y drag", -self.D*np.sin(self.theta))
                #print("theta! ", (180/np.pi)*self.theta, "F_y: ", self.F_y)
                #print("x lift: ", self.L*np.sin(self.theta), "x drag: ", - self.D*np.cos(self.theta), "velo: ", (self.x_vel_arr[-1]**2 + self.y_vel_arr[-1]**2)**0.5, "theta: ", (180/np.pi)*self.theta)
                #print("F_y (N): ", self.F_y, "F_x (N): ", self.F_x, "current len (in): ", self.current_len*39.3701)
                #print("theta! ", (180/np.pi)*self.theta, "F_y: ", self.F_y)
        
            print("state entering: move up deccel", self.time_arr[-1])
            self.current_state = 'move_up_deccel'



            """
            state: hold_be_pos_down
            previous state: move_down_accel    next state: move_down_deccel
        
            Returns: nothing
            Updates: net forces F_y, F_x,
            Appends: tank mass and time
            """
        elif self.current_state == 'hold_be_pos_down':
            print("holding time down (s): ", self.time_hold_down)

            current_time = 0

            while current_time < self.time_hold_down:

                self.calc_hydro_force()
                self.kinematics()

                #iterate/move forward current_time
                current_time += TIMESTEP
                self.tank_mass_arr.append(self.tank.m_current)
                self.time_arr.append(self.time_arr[-1]+TIMESTEP)

                """debug print statements"""
                #print("F_y (N): ", self.F_y, "F_x (N): ", self.F_x, "current len (in): ", self.current_len*39.3701)
                #print("Fy: ", self.F_y, "Fb: ",RHO_WATER*GRAVITY*(hull.V_hull + self.V_be), "Fg: ",-m_glider*GRAVITY, "y lift: ", +self.L*np.cos(self.theta), "y drag", self.D*np.sin(self.theta))
                #print("x lift: ", self.L*np.sin(self.theta), "x drag: ", - self.D*np.cos(self.theta))
                #print("theta! ", (180/np.pi)*self.theta, "F_y: ", self.F_y)

            print("state entering: move down deccel", self.time_arr[-1])
            self.current_state = 'move_down_deccel'



            
            """
            state: move_down_deccel
            previous state: hold_be_pos_up     next state: move_up_accel

            Returns: nothing
            Updates: net forces F_y, F_x, tank mass
            Appends: tank mass and time
            """
        elif self.current_state == 'move_down_deccel':
            print("state: move down deccel")
            while self.tank.m_current > tank.m_neutral_buoy:

                self.m_glider = self.tank.update_mass(-self.delta_water_mass, self.m_glider)
                
                self.calc_hydro_force()
                self.kinematics()

                self.time_arr.append(self.time_arr[-1]+TIMESTEP)
                self.tank_mass_arr.append(self.tank.m_current)

                """debug print statements"""
                #print("F_y (N): ", self.F_y, "F_x (N): ", self.F_x, "current len (in): ", self.current_len*39.3701)
                #print("theta! ", (180/np.pi)*self.theta, "F_y: ", self.F_y, "y lift: ", +self.L*np.cos(self.theta), "y drag", +self.D*np.sin(self.theta), "velo: ", (self.x_vel_arr[-1]**2 + self.y_vel_arr[-1]**2)**0.5)
                #print("Fy: ", self.F_y, "Fb: ",RHO_WATER*GRAVITY*(hull.V_hull + self.V_be), "Fg: ",-m_glider*GRAVITY, "y lift: ", +self.L*np.cos(self.theta), "y drag", +self.D*np.sin(self.theta))
                #print("x lift: ", self.L*np.sin(self.theta), "x drag: ", - self.D*np.cos(self.theta))
                #print(self.F_y, 39.3701*self.current_len)

            print("state entering: move up accel", self.time_arr[-1])
            self.sim_depth = self.y_disp_arr[-1]
            #print("y disp arr", self.sim_depth, self.y_disp_arr[-1])
            self.current_state = 'move_up_accel'
    



            """
            state: move_up_accel
            previous state: move_down_deccel    next state: hold_be_pos_up

            Returns: nothing
            Updates: net forces F_y, F_x, tank mass
            Appends: tank mass and time
            """
        elif self.current_state == 'move_up_accel':
            print("state: move up accel")
            while self.tank.m_current > self.tank.m_empty_tank:
                
                self.m_glider = self.tank.update_mass(-self.delta_water_mass, self.m_glider)

                self.calc_hydro_force()
                self.kinematics()

                self.time_arr.append(self.time_arr[-1]+TIMESTEP)
                self.tank_mass_arr.append(self.tank.m_current)

                """debug print statements"""
                #print("Fy: ", self.F_y, "Fb: ",RHO_WATER*GRAVITY*(hull.V_hull + self.V_be), "Fg: ",-m_glider*GRAVITY, "y lift: ", -self.L*np.cos(self.theta), "y drag", -self.D*np.sin(self.theta))
                #print("x lift: ", self.L*np.sin(self.theta), "x drag: ", - self.D*np.cos(self.theta))
                #print("F_y (N): ", self.F_y, "F_x (N): ", self.F_x, "current len (in): ", self.current_len*39.3701)
                #print("theta! ", (180/np.pi)*theta, "F_y: ", self.F_y)
                #print("F_y (N): ", self.F_y, "y lift: ", -self.L*np.cos(self.theta), "y drag: ", -self.D*np.sin(self.theta), "Fg: ", - m_glider*GRAVITY, "Fb: ", RHO_WATER*GRAVITY*(hull.V_hull + self.V_be ), "x velo: ", self.x_vel_arr[-1])
                #print("theta! ", (180/np.pi)*self.theta, "F_y: ", self.F_y)

            print("state entering: hold_be_pos_up, ", self.time_arr[-1])
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
            print("state: move up deccel")
            while self.tank.m_current < self.tank.m_neutral_buoy:

                self.m_glider = self.tank.update_mass(self.delta_water_mass, self.m_glider)

                self.calc_hydro_force()
                self.kinematics()

                self.time_arr.append(self.time_arr[-1]+TIMESTEP)
                self.tank_mass_arr.append(self.tank.m_current)

                """debug print statements"""
                #print("Fy: ", self.F_y, "Fb: ",RHO_WATER*GRAVITY*(hull.V_hull + self.V_be), "Fg: ",-m_glider*GRAVITY, "y lift: ", -self.L*np.cos(self.theta), "y drag", -self.D*np.sin(self.theta))
                #print("x lift: ", self.L*np.sin(self.theta), "x drag: ", - self.D*np.cos(self.theta), "theta: ", self.theta)
                #print("F_y (N): ", self.F_y, "F_x (N): ", self.F_x, "current len (in): ", self.current_len*39.3701)
                #print("theta! ", (180/np.pi)*theta, "F_y: ", self.F_y)
                #print("theta! ", (180/np.pi)*self.theta, "F_y: ", self.F_y)

            print("state entering: end", self.time_arr[-1])
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
the main unknown in the trajectory is the time to hold water in the tank so we sink enough to reach our desired depth.
We solve this ITERATIVELY. Every iteration increases the time we are holding water in the tank until we hold it so long the glider reaches operational depth.

the first loop checks if the trajectory has reached max allowable depth
the second loop(s) run the trajectory 
"""


###INPUT DIRECTLY FROM TCS MICROPUMP SIZING SPREADSHEET AND GORDON'S CFD SPREADSHEET
tank = BallastTank(25.89929537, 4.5,  4, 0.75, 28.346472)
pump = TCSmicropump("MGD2000F")
hull = PressureHull( 5.709, 6.625, 36.4375, 1, 0.02227)
hydrofoil = Hydrofoil(0.06968) # reference area in m^2

m_glider = RHO_WATER*hull.V_displacement
print(m_glider)

m_int = m_glider-( tank.m_neutral_buoy + hull.mass )
x_int = ( hull.x_hull*(m_glider-hull.mass)-tank.x_tank*(tank.m_neutral_buoy) )/ m_int

movingmass = MovingMass(m_int, x_int)

#setting the maximum depth to simulate to
max_allowable_depth = -50 #m
sim_depth = 0 #m

#setup trajectory parameters --> either distance to go to or times to hold down
start_time_hold_down = 0 #s
start_time_hold_up = 0 #s



nom_time_hold_down = 0 #s
nom_time_hold_up = 0 #s




while(sim_depth > max_allowable_depth):
    
    sim_depth = 0
    seapup = SeagliderFSM(tank, pump, movingmass, hull, m_glider, hydrofoil)

    #TODO: DELETE WHEN SCRIPT WORKS AND NO MORE DEBUGGING
    print("\n")
    print("STARTING SIM!!!")
    print("\n")

    ###loop through fsm, #starting from surface with no speed for first run, second run we already have speed and are gliding normally
    while(seapup.get_state() != 'end'):
        seapup.trajectory(start_time_hold_down, start_time_hold_up)
        sim_depth = seapup.sim_depth


    print("\n")
    print(seapup.sim_depth)
    print("\n")

    
    #normal yo
    """
    for i in range(2):
        seapup.set_state('move_down_accel')
        while(seapup.get_state() != 'end'):
            seapup.trajectory(nom_time_hold_down, nom_time_hold_up)

        #TODO: REMOVE AFTER DEBUGGING
    """

    """
    print(' ')
    print("split", seapup.x_disp_arr[-1])
    print(' ')

    seapup.set_state('move_down_accel')
    while(seapup.get_state() != 'end'):
        #def trajectory(self, time_hold_down, time_hold_up)
        seapup.trajectory(nom_time_hold_down, nom_time_hold_up)
    """

    #iterate and reset variables
    start_time_hold_up += 12
    start_time_hold_down += 12

    state = 'start'

    #this is for testing purposes otherwise we get stuck in loop BUG: get stuck in loop
    #NOTE: DONT FORGET THIS WHEN DEBUGGING
    break


sim_depth = np.min(seapup.y_disp_arr)
max_speed = np.max(seapup.x_vel_arr)
glide_period = seapup.x_disp_arr[-1]

#print a summary of noteworthy params
print("max fwd speed: ",max_speed, " (m/s)")
print("glide period: ",glide_period, '(m)')












"""
#big graphing:

plt.subplot(1,1,1)  
plt.plot(seapup.x_disp_arr, seapup.y_disp_arr, label='Seaglider Position', color='blue')
plt.xlabel('X-Pos (m)')
plt.ylabel('Y-Pos (m)')


for i in seapup.s_x_change_down_arr:
    plt.axvline(x=i, color='r', linestyle='-', label='change up')
for i in seapup.s_x_change_up_arr:
    plt.axvline(x=i, color='g', linestyle='-', label='change up')

plt.show()


plt.subplot(1,3,2)
plt.plot(seapup.x_disp_arr, seapup.tank_mass_arr, label='BE Position Over X-Pos', color='blue')
plt.xlabel('X-Pos (m)')
plt.ylabel('B.E. Pos (in)')



plt.subplot(1,3,3)
plt.plot(seapup.time_arr, seapup.tank_mass_arr, label='BE Position Over X-Pos', color='blue')
plt.xlabel('Time (s)')
plt.ylabel('B.E. Pos (in)')


plt.show()
"""

#uncomment for big data

plt.subplot(2,8,1)  
plt.plot(seapup.x_disp_arr, seapup.y_disp_arr, label='Seaglider Position', color='blue')
plt.xlabel('X-Pos (m)')
plt.ylabel('Y-Pos (m)')


for i in seapup.s_x_change_down_arr:
    plt.axvline(x=i, color='r', linestyle='-', label='change up')
for i in seapup.s_x_change_up_arr:
    plt.axvline(x=i, color='g', linestyle='-', label='change up')


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
plt.plot(seapup.x_disp_arr, seapup.tank_mass_arr, label='X-Pos Vs BE Position ', color='blue')
plt.ylabel('Tank Mas (kg)')
plt.xlabel('X-Pos (m)')

plt.subplot(2,8,9)
plt.plot(seapup.tank_mass_arr, seapup.x_accel_arr, label='BE Position Vs x-Accel', color='blue')
plt.xlabel('Tank Mas (kg)')
plt.ylabel('X-accel (m/s^2)')

plt.subplot(2,8,10)
plt.plot(seapup.tank_mass_arr, seapup.y_vel_arr, label='BE Position Vs Y-Vel', color='blue')
plt.xlabel('Tank Mas (kg)')
plt.ylabel('Y-Vel (m/s)')

plt.subplot(2,8,11)
plt.plot(seapup.tank_mass_arr, seapup.y_accel_arr, label='BE Position Vs Y-accel', color='blue')
plt.xlabel('Tank Mas (kg)')
plt.ylabel('Y-Accel (m/s^2)')

plt.subplot(2,8,12)
plt.plot(seapup.time_arr, seapup.tank_mass_arr, label='Time Vs BE Position', color='blue')
plt.xlabel('Time (S)')
plt.ylabel('Tank Mas (kg)')

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
plt.plot(seapup.time_arr, seapup.delta_X_cg_arr, label='Time Vs delta cb (m)', color='blue')
plt.xlabel('Time (S)')
plt.ylabel('delta X_cg (m)')

plt.show()