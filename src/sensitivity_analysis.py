import numpy as np
import matplotlib.pyplot as plt

from src import constants

from src.SeagliderComponents.BuoyancyEngine import BuoyancyEngine
from src.SeagliderComponents.PressureHull import PressureHull

from src.seaglider_trajectory_incremental import seaglider_trajectory

def intometer(x):
    return 0.0254*x

class DataClass:
    def __init__(self, s_x_arr, s_y_arr, be_pos_arr):
        self.s_x_arr = s_x_arr
        self.s_y_arr = s_y_arr
        self.be_pos_arr = be_pos_arr




def sensitivityAnalysis(hfoil_coeff,m_glider,rho_water,be,ph):

    #run sim
    s =seaglider_trajectory(rho_water, be, ph, constants.midpoint, m_glider, hfoil_coeff)
    print(constants.GRAVITY*m_glider, constants.GRAVITY*constants.rho_water*(ph.V_hull+be.V_mid))
    s.incremental_sim()

    #return arrays for plotting
    return DataClass(s.s_x_arr, s.s_y_arr, s.be_pos_arr)


def produce_time_graphs(big_data,i_arr):
    j=1
    while(j<3):

        if(j==1):
            k=0
            for a in big_data:
                plt.subplot(1,2,j)
                plt.plot(a.s_x_arr,a.s_y_arr, label= f'{i_arr[k]}')
                k+=1
            plt.legend()
            plt.xlabel('X-Pos (m)')
            plt.ylabel('Y-Pos (m)')
            plt.title('Seaglider Position')
            plt.grid(True)

        if(j==2):
            k=0
            for a in big_data:
                plt.subplot(1,2,j)
                plt.plot(a.s_x_arr,a.be_pos_arr, label= f'{i_arr[k]}')
                k+=1
            plt.legend()
            plt.xlabel('X-Pos (m)')
            plt.ylabel('B.E. Pos (m)')
            plt.title('BE Position Over X-Pos')
            plt.grid(True)

        j += 1

    plt.show()



def update_i(i):
    i_arr.append(i)
    return i + (constants.max_bound-constants.min_bound)/(constants.num_iterations -1)


i = constants.min_bound
big_data = []
i_arr = []
#malding if statements, you hate to see it for each input, then call function with sensitivity analysis inside

if constants.sensitivity_var=="hfoil_coeff":
    while(i<=constants.max_bound):
        big_data.append( sensitivityAnalysis(constants.sensitivity_var,constants.m_glider,constants.rho_water,constants.buoyeng,constants.preshull) )
        i = update_i(i)
    produce_time_graphs(big_data,i_arr)

if constants.sensitivity_var=="internal_mass":
    while(i<=constants.max_bound):

        #need to recalculate mass since hull sized based on payload mass for neutral buoyancy at be midpoint
        sensitivity_updated_hull_len = (i + constants.buoyeng.mass - constants.rho_water*constants.buoyeng.V_mid) / ( np.pi *0.25* (constants.rho_water*(constants.hull_od**2) - constants.RHO_PVC*(constants.hull_od**2 - constants.hull_id**2)) )

        sensitivity_updated_preshull = PressureHull( constants.hull_id, constants.hull_od, sensitivity_updated_hull_len, constants.percent_stability)
        m_glider = sensitivity_updated_preshull.mass + constants.buoyeng.mass + i 

        big_data.append( sensitivityAnalysis(constants.hfoil_coeff,m_glider,constants.rho_water,constants.buoyeng,constants.preshull) )
        i = update_i(i)
    produce_time_graphs(big_data,i_arr)

if constants.sensitivity_var=="rho_water":
    while(i<=constants.max_bound):
        big_data.append( sensitivityAnalysis(constants.hfoil_coeff,constants.m_glider,i,constants.buoyeng,constants.preshull) )
        i = update_i(i)
    produce_time_graphs(big_data,i_arr)

if constants.sensitivity_var=="be":

    BuoyancyEngineList = [
        #BuoyancyEngine(  id, od, cont_len, travel_len, laspeed, midpoint):
        BuoyancyEngine( intometer(0.5), intometer(1.15), intometer(5.7), intometer(4), 0.08*0.0254, constants.midpoint),
        BuoyancyEngine( intometer(0.75), intometer(1.37), intometer(5.7), intometer(4), 0.08*0.0254, constants.midpoint),
        BuoyancyEngine( intometer(1), intometer(1.68), intometer(5.7), intometer(4), 0.08*0.0254, constants.midpoint),
        BuoyancyEngine( intometer(1.25), intometer(2.04), intometer(5.7), intometer(4), 0.08*0.0254, constants.midpoint),
        BuoyancyEngine( intometer(1.5), intometer(2.29), intometer(5.7), intometer(4), 0.08*0.0254, constants.midpoint),
        BuoyancyEngine( intometer(2), intometer(2.82), intometer(10.3), intometer(8), 0.08*0.0254, constants.midpoint),
        BuoyancyEngine( intometer(2.5), intometer(3.42), intometer(10.3), intometer(8), 0.08*0.0254, constants.midpoint),
        BuoyancyEngine( intometer(3), intometer(4), intometer(10.3), intometer(8), 0.08*0.0254, constants.midpoint),
        BuoyancyEngine( intometer(3.5), intometer(4.63), intometer(10.3), intometer(8), 0.08*0.0254, constants.midpoint),
        BuoyancyEngine( intometer(4), intometer(5.15), intometer(10.3), intometer(8), 0.08*0.0254, constants.midpoint),
        BuoyancyEngine( intometer(5), intometer(6.26), intometer(10.3), intometer(8), 0.08*0.0254, constants.midpoint),
    ]

    for b in BuoyancyEngineList:

        #need to adjust hull dim

        hull_len = (constants.internal_mass + b.mass - constants.rho_water*b.V_mid) / ( np.pi *0.25* (constants.rho_water*(constants.hull_od**2) - constants.RHO_PVC*(constants.hull_od**2 - constants.hull_id**2)) )

        ph = PressureHull( constants.hull_id, constants.hull_od, hull_len, constants.percent_stability)

        m_glider = ph.mass + b.mass + constants.internal_mass

        print(constants.GRAVITY*m_glider, constants.GRAVITY*constants.rho_water*(ph.V_hull+b.V_mid))
        print("VOLUME!!!",ph.V_hull,b.V_mid)

    
        big_data.append( sensitivityAnalysis(constants.hfoil_coeff,m_glider,constants.rho_water,b,ph) )
        i = update_i(i)
    produce_time_graphs(big_data,i_arr)

