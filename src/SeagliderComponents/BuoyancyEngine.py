from src import constants
import numpy as np



class BuoyancyEngine:
    def __init__(self, id, od, cont_len, travel_len, laspeed, midpoint):
        self.id = id
        self.od = od
        self.cont_len = cont_len
        self.travel_len = travel_len
        self.laspeed = laspeed #m/s
        self.midpoint = midpoint #passed in in metric!

        #displaced volume
        self.V_cont = 0.25*(np.pi*self.od**2)*self.cont_len
        self.V_mid = 0.25*np.pi*(self.id**2)*self.midpoint + self.V_cont

        #self.V_ext = 0.25*np.pi*(self.id**2)*(midpoint+self.cont_len) #BUG: use inner diameter for whole length, likely the issue here

        self.allowable_ext = midpoint

        self.mass = constants.RHO_PVC*0.25*np.pi*(self.od**2-self.id**2)*self.cont_len
        #print("be mass: ", self.mass)


#print("inside buoyancye engin1!!")