import numpy as np
import matplotlib.pyplot as plt

###INPUTS

#Constants
rho_w = 997 #kg/m^3
rho_pvc = 1380 #kg/m^3
g = 9.81 #m/s^2

#Pressure Hull
diam_hull = 4.5 #in
l_hull = 1.23698 #m #NOTE: SENSITIVE INPUT

diam_hull = 0.0254*diam_hull #m

V_hull = np.pi*l_hull*(diam_hull/2)**2 #m^3

#Buoyancy Engine
diam_be = 4.07 #in
l_be_cont = 10.3 #in
l_travel = 8 #in
l_be_mid = l_be_cont+0.5*l_travel #in

diam_be = 0.0254*diam_be #m
l_be_cont = 0.0254*l_be_cont #m
l_travel = 0.0254*l_travel #m
l_be_mid = 0.0254*l_be_mid #m

V_be_cont = np.pi*l_be_cont*(diam_be/2)**2 #m^3
V_be_mid = np.pi*(l_be_cont+l_travel/2)*(diam_be/2)**2 #m^3
#print("START", V_be_mid, V_hull)
#V_be_exp = np.pi*(l_be_cont+l_travel)*(diam_be/2)**2 #m^3

#Hydrofoil setup
hfoil_coeff = 0.015 #Area of wing * Coeff Lift
v = 2 #m/s

#Other
percent_stability = 0.1 #%
stability = diam_hull * percent_stability #m
m_glider = 15.6 #kg #TODO: WANT TO CALCULATE THIS AS A FUNCTION OF LENGTH INSTEAD OF HARDCODING also sensitive input
dt = 0.25 #s
delta_be = 0.01 #in

#delta_be = 0.0254*delta_be #m
i = 0

###SETUP
s_y_prev = 0 #m
s_x_prev = 0 #m

v_y_prev = 0.00001 #m/s
v_x_prev = v #m/s --> starting from top so all velocity in x

L = 0.5*rho_w*hfoil_coeff*v**2
delta_l = 0

s_x_arr = []
s_y_arr = []

###TODO: FORGOT L_BE CONTRACT IN BUOYANCY

###START Trajectory Sim --> going down (L+)
while(i < 2*np.pi):
    #Update glider
    V_be = V_be_mid + np.pi*delta_l*(diam_be/2)**2 #m^3
    F_b = rho_w*g*(V_hull + V_be)

    #TODO: THETA CALCULATION IS WRONG!!!!! ALWAYS AT 90 DEG
    #L IS WRONG!!!!!!
    theta = np.arctan( ( ( (l_hull -0.6906110606)+(l_be_mid+delta_l)/2 ) * (rho_w*g*(V_be-V_be_mid) ) ) / (m_glider*g*stability) )
    print("theta", (180/np.pi)*theta, "delta_l", delta_l, "V_be", V_be, "y lift", L*np.sin(theta))
    F_y = rho_w*g*(V_hull + V_be ) - m_glider*g + L*np.sin(theta)
    #print("theta", (180/np.pi)*theta, "force", F_y, "y lift", L*np.sin(theta), "buoyancy", rho_w*g*(V_hull+ V_be ), "gravity", -m_glider*g)

    #Find resultant forces 
    F_x = L*np.cos(theta)

    #"integrate" w timestep
    v_y = (F_y/m_glider)*dt
    v_x = (F_x/m_glider)*dt

    v = np.sqrt(v_x**2+v_y**2)
    #print (v_x, v_y)
    
    
    v_prev = np.sqrt(v_x_prev**+v_y**2)
    v_vector = np.array([v_x, v_y])
    v_prev_vector = np.array([v_x_prev,v_y_prev]) #TODO: FIX, DOESNT SEEM TO LIKE INPUTS LESS THAN ZERO

    delta_theta = v_vector.dot(v_prev_vector)
    #print(delta_theta)
    

    s_y = s_y_prev + v_y*dt
    s_x = s_x_prev + v_x*dt

    #graph array and setup for next loop
    s_x_arr.append(s_x)
    s_y_arr.append(s_y)

    s_y_prev = s_y
    s_x_prev = s_x

    v_y_prev = v_y
    v_x_prev = v_x

    i = i +0.01
    delta_l = -(0.5*l_travel)*np.sin(i)




###CONTINUE Trajectory Sim --> going up (L-)
"""
i=0
while(i < np.pi):
    #Update glider
    V_be = V_be_mid + delta_l*diam_be #m^3

    F_b = rho_w*g*(V_hull+V_be)
    print("Fb", F_b, "l", delta_l, "V_be", V_be)

    theta = -np.tanh( ((l_hull+delta_l/2)*rho_w*g*V_be) / (m_glider*g*stability) )

    #Find resultant forces 
    F_y = rho_w*g*(V_hull+V_be) - m_glider*g - L*np.sin(theta)
    #print("Fy", F_y)

    F_x = L*np.cos(theta)

    #"integrate" w timestep
    v_y = v_y_prev + (F_y/m_glider)*dt
    v_x = v_x_prev + (F_x/m_glider)*dt

    s_y = s_y_prev + v_y*dt
    s_x = s_x_prev + v_x*dt

    #graph array and setup for next loop
    s_x_arr.append(s_x)
    s_y_arr.append(s_y)

    s_y_prev = s_y
    s_x_prev = s_x

    i = i + 0.05
    delta_l = l_travel*np.sin(i)
    #print("delta_l", delta_l)
"""


###GRAPH
plt.plot(s_x_arr, s_y_arr, label='Seaglider Position', color='blue')
plt.xlabel('X-Pos (m)')
plt.ylabel('Y-Pos (m)')
plt.show()