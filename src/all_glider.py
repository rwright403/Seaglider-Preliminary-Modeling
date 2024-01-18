import numpy as np
import matplotlib.pyplot as plt

from src import constants

from src.SeagliderComponents.BuoyancyEngine import BuoyancyEngine
from src.SeagliderComponents.PressureHull import PressureHull


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


class TubeSize:
    def __init__(self,id,od):
        self.dimensions = [id, od]

TubeSizeList = [
    TubeSize(4,5),
    TubeSize(5,5.563),
    TubeSize(6,6.625),
]



for b in BuoyancyEngineList:

    for t in TubeSizeList:










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


