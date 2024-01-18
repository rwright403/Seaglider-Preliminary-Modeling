import numpy as np

from src.SeagliderComponents.BuoyancyEngine import BuoyancyEngine
from src.SeagliderComponents.PressureHull import PressureHull

def intometer(x):
    return 0.0254*x


####CONSTANTS:

GRAVITY = 9.81 #m/s^2
TIMESTEP = 0.5 #s
BE_EXTENSION_STEP = intometer(0.0375)
max_allowable_depth = -15 #m

RHO_PVC = 1380 #kg/m^3
GRAVITY = 9.81 #m/s^2
TIMESTEP = 0.5 #s


####inputs:
hfoil_coeff = 0.008 #Area of wing * Coeff Lift
percent_stability = 0.2 #%
rho_water =997 #kg/m^3
midpoint = intometer(4)
internal_mass = 10 #kg

#geometry
buoyeng = BuoyancyEngine(intometer(4),intometer(5.15),intometer(10.3), intometer(8), 0.08*0.0254, midpoint)

hull_id = intometer(5.0)
hull_od = intometer(5.5)

hull_len = (internal_mass + buoyeng.mass - rho_water*buoyeng.V_mid) / ( np.pi *0.25* (rho_water*(hull_od**2) - RHO_PVC*(hull_od**2 - hull_id**2)) )

preshull = PressureHull( hull_id, hull_od, hull_len, percent_stability)

m_glider = preshull.mass + buoyeng.mass +internal_mass 



#Buoyancy Engine Combinations:
"""
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
"""

#list of nom tube dimensions
"""
TubeSizeList = [
    TubeSize(4,5),
    TubeSize(5,5.563),
    TubeSize(6,6.625),
]
"""

###SENSITIVITY ANALYSIS INFORMATION!!!!
"""
Variables to analyize:
hfoil_coeff
internal_mass
rho_water

be
"""

sensitivity_var = "be"
min_bound = 6
max_bound = 16
num_iterations=5