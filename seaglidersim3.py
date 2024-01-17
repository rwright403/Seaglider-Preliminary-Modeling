import numpy as np
from pymoo.factory import get_problem, get_algorithm
from pymoo.optimize import minimize
import matplotlib.pyplot as plt

def intometer(x):
    return 0.0254*x

####CONSTANTS:
RHO_PVC = 1380 #kg/m^3
GRAVITY = 9.81 #m/s^2
TIMESTEP = 0.5 #s
BE_EXTENSION_STEP = intometer(0.625) #TODO: FILL IN WITH REAL VALUES LATER!!!!!


#input all dimensions as imperial!
class BuoyancyEngine:
    def __init__(self, id, od, cont_len, travel_len, laspeed, midpoint):
        self.id = intometer(id)
        self.od = intometer(od)
        self.cont_len = intometer(cont_len)
        self.travel_len = intometer(travel_len)
        self.laspeed = laspeed #m/s

        #displaced volume
        self.V_cont = 0.25*(np.pi*self.od**2)*self.cont_len
        self.V_ext = 0.25*np.pi*(self.id**2)*(intometer(midpoint)+self.cont_len)

        self.allowable_ext = intometer(midpoint)

        self.mass = RHO_PVC*0.25*np.pi*(self.od**2-self.id**2)*self.cont_len
        #print("be mass: ", self.mass)

class TubeSize:
    def __init__(self,id,od):
        self.dimensions = [id, od]

class PressureHull:
    def __init__(self,tubesize,len,percent_stability):
        self.id = intometer(tubesize.dimensions[0])
        self.od = intometer(tubesize.dimensions[1])
        self.l_hull = len
        self.percent_stability = percent_stability

        #print(self.od)
        
        self.stability = percent_stability*self.od
        self.V_hull = 0.25*np.pi*self.l_hull*self.od**2

        self.mass = RHO_PVC*(np.pi/4*(self.od**2-self.id**2)*self.l_hull)
        #print("hull mass: ", self.mass)


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

        ###SETUP CONSTANTS
        s_y_prev = 0 #m
        s_y = 0 #m
        s_x_prev = 0 #m
        s_x = 0 #m

        v = 1 #m/s
        #v_y_prev = 0.00001 #m/s
        #v_x_prev = v #m/s --> starting from top so all velocity in x

        v_x_arr = []
        max_speed = 0

        s_x_arr = []
        s_y_arr = []
        s_x_change_down_arr = []
        s_x_change_up_arr = []

        time_arr = []
        be_pos_arr = []
        #--> going down (L+)
        L = 0.5*rho_water*hfoil_coeff*v**2


        #setup current_len
        current_len = midpoint
        time = 0

        diving = True
        contunuing_criteria = True
        while (contunuing_criteria == True):
            #calc be volume
            V_be = be.V_cont + (0.25*np.pi*be.id**2)*current_len

            #solve glide angle
            theta = np.arctan( ( hull.l_hull * (rho_water*GRAVITY*(0.25*np.pi*be.id**2)*current_len) ) / (m_glider*GRAVITY*hull.stability) )

            #solve Fy, Fx
            if(s_y <= s_y_prev):
                F_y = rho_water*GRAVITY*(hull.V_hull + V_be ) - m_glider*GRAVITY + L*np.sin(theta)
            else:
                F_y = rho_water*GRAVITY*(hull.V_hull + V_be ) - m_glider*GRAVITY - L*np.sin(theta)
            
            F_x = L*np.cos(theta)

            #print(rho_water*GRAVITY*(hull.V_hull + V_be ), m_glider*GRAVITY, L*np.sin(theta), F_y)
            #print(s_x, F_y)

            #"integrate" w timestep
            v_y = (F_y/m_glider)*TIMESTEP
            v_x = (F_x/m_glider)*TIMESTEP

            v_x_arr.append(v_x)
            #print(v_x)

            #print(v_y,v_x)

            s_y = s_y_prev + v_y
            s_x = s_x_prev + v_x

            s_x_arr.append(s_x)
            s_y_arr.append(s_y)

            s_y_prev = s_y
            s_x_prev = s_x


            if( current_len <= (midpoint - be.allowable_ext) and diving == True):
                diving = False
                s_x_change_up_arr.append(s_x)
                #print("going up", s_x)

            if( current_len >= (midpoint + be.allowable_ext) and diving == False):
                diving = True
                contunuing_criteria = False
                s_x_change_down_arr.append(s_x)
                #print("going down", s_x)

            if diving == True:
                current_len -= be.laspeed*TIMESTEP
            else:
                current_len += be.laspeed*TIMESTEP

            #print(theta, F_y, v_y, v_x, s_x, current_len)
            time = time + TIMESTEP
            time_arr.append(time)
            be_pos_arr.append(39.3701*current_len)

        #for debugging
        """
        plt.subplot(1,2,1)  
        plt.plot(s_x_arr, s_y_arr, label='Seaglider Position', color='blue')
        for i in s_x_change_down_arr:
            plt.axvline(x=i, color='r', linestyle='--', label='change up')
        for i in s_x_change_up_arr:
            plt.axvline(x=i, color='g', linestyle='--', label='change up')

        plt.xlabel('X-Pos (m)')
        plt.ylabel('Y-Pos (m)')

        plt.subplot(1,2,2)
        plt.plot(s_x_arr, be_pos_arr, label='BE Position Over X-Pos', color='blue')
        plt.xlabel('X-Pos (m)')
        plt.ylabel('B.E. Pos (m)')
        plt.show()
        """

        """
        sim_depth = np.min(s_y_arr)
        max_speed = np.max(v_x_arr)
        glide_period = (4/3)*s_x_arr[-1]
        """

        ###iterate!!!!!
        be.allowable_ext+= BE_EXTENSION_STEP
        if(be.allowable_ext>be.travel_len):
            break
        #print(sim_depth)
###TODO: STILL BUGGY CHECK INPUTS BUG

    print("allowable extension: ",be.allowable_ext)
    print("max fwd speed: ",max_speed)
    print("glide period: ",glide_period)

    output_arr = []
    #[be.allowable_ext, -max_speed, -glide_period]
    output_arr.append(be.allowable_ext)
    output_arr.append(-max_speed)
    output_arr.append(-glide_period)

    output_arr.append(be)
    output_arr.append(hull)
    return output_arr




problem = get_problem("zdt1")

####inputs:
hfoil_coeff = 0.008 #Area of wing * Coeff Lift
percent_stability = 0.2 #%
rho_water =997 #kg/m^3
midpoint = intometer(4)
internal_mass = 8.05 #kg

#setup midpoint!!!!

#list of buoyancy engines
BuoyancyEngineList = [
    #BuoyancyEngine(  id, od, cont_len, travel_len, laspeed, midpoint):
    BuoyancyEngine( 0.5, 1.15, 5.7, 4, 0.08*0.0254, midpoint),
    BuoyancyEngine( 0.75, 1.37, 5.7, 4, 0.08*0.0254, midpoint),
    BuoyancyEngine( 1, 1.68, 5.7, 4, 0.08*0.0254, midpoint),
    BuoyancyEngine( 1.25, 2.04, 5.7, 4, 0.08*0.0254, midpoint),
    BuoyancyEngine( 1.5, 2.29, 5.7, 4, 0.08*0.0254, midpoint),
    BuoyancyEngine( 2, 2.82, 10.3, 8, 0.08*0.0254, midpoint),
    BuoyancyEngine( 2.5, 3.42, 10.3, 8, 0.08*0.0254, midpoint),
    BuoyancyEngine( 3, 4, 10.3, 8, 0.08*0.0254, midpoint),
    BuoyancyEngine( 3.5, 4.63, 10.3, 8, 0.08*0.0254, midpoint),
    BuoyancyEngine( 4, 5.15, 10.3, 8, 0.08*0.0254, midpoint),
    BuoyancyEngine( 5, 6.26, 10.3, 8, 0.08*0.0254, midpoint),
]

#list of nom tube dimensions
TubeSizeList = [
    TubeSize(4,5),
    TubeSize(5,5.563),
    TubeSize(6,6.625),
]


"""
#NOTE: we want a fast seaglider to overcome currents, 
but we dont want it too fast that the BE speed is too sensitive
--> put some bounds on the actuation range?
---> probably just want to target a minimum mass? dont really want to optimize 
"""

# Define your problem as a pymoo Problem
class GliderOptimizationProblem(problem):
    def __init__(self, BuoyancyEngineList, TubeSizeList):
        super().__init__(
            n_var=4,  # Number of design variables
            n_obj=3,  # Number of objectives
            n_constr=1,  # Number of constraints
            #mass, midpoint, buoyancy engine list, tube list
            xl=np.array([9, 0, BuoyancyEngineList[0],TubeSizeList[0]]),  # Lower bounds for design variables
            xu=np.array([14, 4, BuoyancyEngineList[-1],TubeSizeList[-1]]),  # Upper bounds for design variables
            elementwise_evaluation=True
        )

    def _evaluate(self, x, out, *args, **kwargs):
        # Extract design variables
        #m_glider, buoyancy engine list, tube list, midpoint = x

        # Simulation logic using the design variables
        # ... (Replace the constants with the corresponding design variables)
        #seaglider_trajectory(rho_water, buoyeng, preshull, midpoint, m_glider, hfoil_coeff)
        #objective functions!

        #[be.allowable_ext, -max_speed, -glide_period]
        out["F"] = seaglider_trajectory(rho_water, buoyeng, preshull, midpoint, m_glider, hfoil_coeff)
        #want to maximize glider speed? --> see notebook

        # Constraints:
        #out["G"] = 


# Create an instance of your optimization problem
problem = GliderOptimizationProblem(BuoyancyEngineList,TubeSizeList)

# Choose genetic algorithm as the optimization algorithm
algorithm = get_algorithm("nsga2",
    pop_size=100,
    crossover_prob = 0.9,
    mutation_prob=1.0/problem.n_var,
    eliminate_duplicates=True
)

# Optimize the problem using the selected algorithm
result = minimize(problem,
                  algorithm,
                  termination=('n_gen', 100),
                  seed=1,
                  save_history=True,
                  verbose=True)

# Extract and print the optimal design variables
optimal_design_variables = result.X
print("Optimal Design Variables:", optimal_design_variables)

# Update the simulation logic with the optimal design variables


# Simulation logic with optimized design variables
# ... (Replace the constants with the optimal design variables)

# Example: Re-run the simulation with the optimized design variables
# s_x_arr_opt, s_y_arr_opt, time_arr_opt, be_pos_arr_opt = simulate_glider_trajectory(l_hull_opt, m_glider_opt, hfoil_coeff_opt, l_travel_opt)
# ...

# Print additional information or outputs as needed
# ...
"""
# Plot the results for both the original and optimized scenarios
plt.subplot(1, 2, 1)
plt.plot(s_x_arr, s_y_arr, label='Original Seaglider Position', color='blue')
plt.axvline(x=s_x_change_down, color='r', linestyle='--', label='change up')
plt.xlabel('X-Pos (m)')
plt.ylabel('Y-Pos (m)')
plt.title('Original Seaglider Trajectory')

# Plot the optimized trajectory
plt.plot(s_x_arr_opt, s_y_arr_opt, label='Optimized Seaglider Position', color='green')

plt.legend()

plt.subplot(1, 2, 2)
plt.plot(time_arr, be_pos_arr, label='Original BE Position Over Time', color='blue')
plt.xlabel('Time (s)')
plt.ylabel('B.E. Pos (m)')
plt.title('Original BE Position Over Time')

# Plot the optimized BE position over time
plt.plot(time_arr_opt, be_pos_arr_opt, label='Optimized BE Position Over Time', color='green')

plt.legend()

plt.show()
"""