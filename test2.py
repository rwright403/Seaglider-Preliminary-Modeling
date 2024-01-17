import numpy as np
import matplotlib.pyplot as plt

def intometer(x):
    return 0.0254*x

####CONSTANTS:
RHO_PVC = 1380 #kg/m^3
GRAVITY = 9.81 #m/s^2
TIMESTEP = 0.5 #s
BE_EXTENSION_STEP = 0.006 #TODO: FILL IN WITH REAL VALUES LATER!!!!!

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


class PressureHull:
    def __init__(self,id,od,len,percent_stability):
        self.id = id
        self.od = od
        self.l_hull = len
        self.percent_stability = percent_stability
        
        self.stability = percent_stability*self.od
        self.V_hull = 0.25*np.pi*self.l_hull*self.od**2

        self.mass = RHO_PVC*(np.pi/4*(self.od**2-self.id**2)*self.l_hull)
        print("hull mass: ", self.mass)


####inputs:
hfoil_coeff = 0.008 #Area of wing * Coeff Lift
percent_stability = 0.2 #%
rho_water =997 #kg/m^3
midpoint = 4
m_internal = 8.05 #kg

#geometry
buoyeng = BuoyancyEngine(3.0,4,10.3, 8, 0.08*0.0254,midpoint)

hull_id = intometer(4.0)
hull_od = intometer(4.5)

hull_len = (m_internal + buoyeng.mass - rho_water*buoyeng.V_ext) / ( 3.14 *0.25* (rho_water*(hull_od**2) - RHO_PVC*(hull_od**2 - hull_id**2)) )

preshull = PressureHull( hull_id, hull_od, hull_len, percent_stability)

m_glider = m_internal+buoyeng.mass+preshull.mass

print(m_glider)