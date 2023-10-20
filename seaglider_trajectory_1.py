import numpy as np
import matplotlib.pyplot as plt

###INPUTS

#Constants
rho_w = 997 #kg/m^3
g = 9.81 #m/s^2

#Pressure Hull
diam_hull = 4.5 #in
l_hull = 0.70358 #in

diam_hull = 0.0254*diam_hull #m
l_hull = 0.0254*l_hull #m

V_hull = np.pi*l_hull*(diam_hull/2)**2 #m^3

#Buoyancy Engine
diam_be = 4.07 #in
l_be_cont = 10.3 #in
l_travel = 4 #in

diam_be = 0.0254*diam_be #m
l_be_cont = 0.0254*l_be_cont #m
l_travel = 0.0254*l_travel #m

V_be_cont = np.pi*l_be_cont*(diam_be/2)**2 #m^3

#Hydrofoil setup
hfoil_coeff = 0.0025 #Area of wing * Coeff Lift
v = 2 #m/s

#Other
percent_stability = 0.1 #%
stability = diam_hull * percent_stability #m
m_glider = 13.6 #kg
dt = 0.1 #s
delta_be = 0.01 #in

#delta_be = 0.0254*delta_be #m
x = 0

###SETUP
s_y_prev = 0 #m
s_x_prev = 0 #m

v_y_prev = 0 #m/s
v_x_prev = v #m/s --> starting from top so all velocity in x

L = 0.5*rho_w*hfoil_coeff*v**2
delta_l = l_travel/2

s_x_arr = []
s_y_arr = []

###TODO: MIXING UP BUOYANCY ENGINE POSITIONS AND LIFT SIGN IN LOOP

###START Trajectory Sim --> going down (L+)
while(x < 2*np.pi):
    #Update glider
    V_be = V_be_cont + delta_l*diam_be #m^3

    F_b = rho_w*g*(V_hull+V_be)

    theta = np.tanh( ((l_hull+delta_l/2)*rho_w*g*V_be) / (m_glider*g*stability) )

    #Find resultant forces 
    F_y = rho_w*g*(V_hull+V_be) - m_glider*g + L*np.sin(theta)
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

    x = x +0.05
    delta_l = l_travel*np.sin(x+np.pi/2)

###CONTINUE Trajectory Sim --> going up (L-)
print(x)
delta_l = l_travel/2
x = 0
while(x < 2*np.pi):
    #Update glider
    V_be = V_be_cont + delta_l*diam_be #m^3

    F_b = rho_w*g*(V_hull+V_be)

    theta = np.tanh( ((l_hull+delta_l/2)*rho_w*g*V_be) / (m_glider*g*stability) )

    #Find resultant forces 
    F_y = rho_w*g*(V_hull+V_be) - m_glider*g - L*np.sin(theta)
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

    x = x + 0.05
    delta_l = l_travel*np.cos(x-np.pi)
    print(delta_l)

###GRAPH
plt.plot(s_x_arr, s_y_arr, label='Seaglider Position', color='blue')
plt.xlabel('X-Pos (m)')
plt.ylabel('Y-Pos (m)')
plt.show()