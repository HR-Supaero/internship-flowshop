# Runs no wait flowshop optimization
import numpy as np

import random

from nowait_flowshop import MIP_model, CP_model

import time

##### Data #####
nb_machines = 5
nb_jobs = 10

processing_times = np.array([[random.randint(1, 10) for j in range(nb_jobs)] for i in range(nb_machines)])



##### Optimization #####

# Measure and print solving time for MIP model
print('-----------MIP-----------')
start_time = time.time()
mip_solution = MIP_model(processing_times, plotGantt=0)
mip_time = time.time() - start_time


# Measure and print solving time for CP model
print('-----------CP-----------')
start_time = time.time()
cp_solution = CP_model(processing_times, plot_gantt=0)
cp_time = time.time() - start_time

print(f"MIP Solving Time: {mip_time} seconds")
print(f"CP Solving Time: {cp_time} seconds")







