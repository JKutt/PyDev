import numpy as np
data_directory = "../testDataIP/"
import matplotlib.pyplot as plt

from simpegEMIP.StretchedExponential import SEInvProblem, SESurvey
from SimPEG import *

time_gates = "30:50,50:70,70:100,100:130,130:170,170:210,210:260,260:320,320:390,390:470,470:560,560:660,660:770,770:890"
data = "0.780       0.626       0.522       0.456       0.401       0.347       0.309       0.273       0.236       0.208       0.180       0.155       0.135       0.124"
Vp = 24.388
time_gates = np.array([np.array(gate.split(":"), dtype=float) for gate in time_gates.split(',')])
data = np.array(data.split(), dtype=float)
times = np.sqrt(time_gates[:,0] * time_gates[:,1]) * 1e-3
print(time_gates[:,0] * time_gates[:,1])