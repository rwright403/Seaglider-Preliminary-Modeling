import numpy as np
import matplotlib.pyplot as plt

from src import constants

from src.SeagliderComponents.BuoyancyEngine import BuoyancyEngine
from src.SeagliderComponents.PressureHull import PressureHull

def intometer(x):
    return 0.0254*x

class seaglider_trajectory():
    def __init__(self,rho_water, be, hull, midpoint, m_glider, hfoil_coeff):
        
        self.rho_water = rho_water
        self.be = be
        self.hull = hull
        self.midpoint = midpoint
        self.m_glider = m_glider
        self.hfoil_coeff = hfoil_coeff

        self.s_x_arr = []
        self.s_y_arr = []
        self.time_arr = []
        self.be_pos_arr = []

        


    def incremental_sim(self):
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

            self.s_x_arr.clear()
            self.s_y_arr.clear()
            self.time_arr.clear()
            self.be_pos_arr.clear()
            
            s_x_change_down_arr = []
            s_x_change_up_arr = []

            #--> going down (L+)
            L = 0.5*self.rho_water*constants.hfoil_coeff*v**2


            #setup current_len
            current_len = constants.midpoint
            #print("current length: ", current_len)
            time = 0

            diving = True
            contunuing_criteria = True
            while (contunuing_criteria == True):

                V_be = self.be.V_cont + (0.25*np.pi*self.be.id**2)*current_len

                #solve glide angle
                theta = np.arctan( ( self.hull.l_hull * (self.rho_water*constants.GRAVITY*(0.25*np.pi*self.be.id**2)*current_len) ) / (constants.m_glider*constants.GRAVITY*self.hull.stability) )

                #solve Fy, Fx
                if(s_y <= s_y_prev):
                    F_y = self.rho_water*constants.GRAVITY*(self.hull.V_hull + V_be ) - constants.m_glider*constants.GRAVITY + L*np.sin(theta)
                else:
                    F_y = self.rho_water*constants.GRAVITY*(self.hull.V_hull + V_be ) - constants.m_glider*constants.GRAVITY - L*np.sin(theta)
                
                F_x = L*np.cos(theta) #NOTE: Tried adding 0.7 to sim drag this worked, separate term seemed to nuke simnot sure????

                #"integrate" w timestep
                v_y = (F_y/constants.m_glider)*constants.TIMESTEP
                v_x = (F_x/constants.m_glider)*constants.TIMESTEP

                v_x_arr.append(v_x)
                #print(v_x)

                #print("vertical speed: ", v_y,"horizontal speed: ", v_x)

                s_y = s_y_prev + v_y
                s_x = s_x_prev + v_x

                self.s_x_arr.append(s_x)
                self.s_y_arr.append(s_y)

                s_y_prev = s_y
                s_x_prev = s_x


                if( current_len <= (constants.midpoint - self.be.allowable_ext) and diving == True):
                    diving = False
                    s_x_change_up_arr.append(s_x)
                    #print("going up", s_x)

                if( current_len >= (constants.midpoint + self.be.allowable_ext) and diving == False):
                    diving = True
                    contunuing_criteria = False
                    s_x_change_down_arr.append(s_x)
                    #print("going down", s_x)

                if diving == True:
                    current_len -= self.be.laspeed*constants.TIMESTEP
                else:
                    current_len += self.be.laspeed*constants.TIMESTEP

                #print(theta, F_y, v_y, v_x, s_x, current_len)
                time = time + constants.TIMESTEP
                self.time_arr.append(time)
                self.be_pos_arr.append(39.3701*current_len)
           
            
            sim_depth = np.min(self.s_y_arr)
            max_speed = np.max(v_x_arr)
            glide_period = (4/3)*self.s_x_arr[-1]

            #after sim check max allowable depth. If we havent hit save, if we have break out of while loop
            #if(sim_depth < max_allowable_depth):
            #    be.allowable_ext = max_be_extension
            """
            plt.subplot(1,2,1)  
            plt.plot(self.s_x_arr, self.s_y_arr, label='Seaglider Position', color='blue')
            for i in s_x_change_down_arr:
                plt.axvline(x=i, color='r', linestyle='--', label='change up')
            for i in s_x_change_up_arr:
                plt.axvline(x=i, color='g', linestyle='--', label='change up')

            plt.xlabel('X-Pos (m)')
            plt.ylabel('Y-Pos (m)')

            plt.subplot(1,2,2)
            plt.plot(self.s_x_arr, self.be_pos_arr, label='BE Position Over X-Pos', color='blue')
            plt.xlabel('X-Pos (m)')
            plt.ylabel('B.E. Pos (m)')
            plt.show()
            ###iterate!!!!!
            self.be.allowable_ext+= constants.BE_EXTENSION_STEP
            """
            #print(sim_depth)

        plt.subplot(1,2,1)  
        plt.plot(self.s_x_arr, self.s_y_arr, label='Seaglider Position', color='blue')
        for i in s_x_change_down_arr:
            plt.axvline(x=i, color='r', linestyle='--', label='change up')
        for i in s_x_change_up_arr:
            plt.axvline(x=i, color='g', linestyle='--', label='change up')

        plt.xlabel('X-Pos (m)')
        plt.ylabel('Y-Pos (m)')

        plt.subplot(1,2,2)
        plt.plot(self.s_x_arr, self.be_pos_arr, label='BE Position Over X-Pos', color='blue')
        plt.xlabel('X-Pos (m)')
        plt.ylabel('B.E. Pos (m)')
        plt.show()


        print("allowable extension: ",self.be.allowable_ext)
        print("max fwd speed: ",max_speed)
        print("glide period: ",glide_period,'\n')


        #output an array of key statistics
        output_arr.append(self.be.allowable_ext)
        output_arr.append(max_speed)
        output_arr.append(glide_period)

        output_arr.append(self.be)
        output_arr.append(self.hull)
        return output_arr