""" SEAGLIDER 2024!!! """
from src import constants

def run():
    user_input = 0
    while(user_input ==0):
        print("\n")
        print("1 --> Trajectory Sim")
        print("2 --> Sensitivity Analysis")
        print("3 --> Optimization")

        user_input = input("Enter number to select analysis: ")

    if user_input =='1':
        print("Running Trajectory Sim from inputs in constants.py")
        from src.seaglider_trajectory_incremental import seaglider_trajectory
        s =seaglider_trajectory(constants.rho_water, constants.buoyeng, constants.preshull, constants.midpoint, constants.m_glider, constants.hfoil_coeff)
        s.incremental_sim()

    if user_input =='2':
        print("Running Sensitivity Analysis")
        from src.sensitivity_analysis import sensitivityAnalysis


    if user_input =='3':
        print("here")