from src import constants
import numpy as np



class PressureHull:
    def __init__(self,id,od,len,percent_stability):
        self.id = id
        self.od = od
        self.l_hull = len
        self.percent_stability = percent_stability

        #print(self.od)
        
        self.stability = percent_stability*self.od
        self.V_hull = 0.25*np.pi*self.l_hull*self.od**2

        self.mass = constants.RHO_PVC*(np.pi/4*(self.od**2-self.id**2)*self.l_hull)
        #print("hull mass: ", self.mass)