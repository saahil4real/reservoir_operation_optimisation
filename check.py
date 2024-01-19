from platypus import *
# from platypus.algorithms import NSGAII
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

df = pd.read_excel('ravishankar.ods', engine='odf')

# inflow = np.asarray(df["Inflow \n(MCM)"])
# evap = np.asarray(df["Evaporation\n(MCM)"])

# print(inflow)
print(df.columns)
# print(evap)