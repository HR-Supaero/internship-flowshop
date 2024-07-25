# This code will run an optimization on a no wait flowshop with setup times

import numpy as np
import random
from nowait_setup_flowshop import MIP_model, CP_model
import time

##### Data #####
nb_machines = 5
nb_jobs = 10

processing_times = np.array([[random.randint(1, 10) for j in range(nb_jobs)] for i in range(nb_machines)])
setup_times = {i: np.array([[random.randint(1,5) for j in range(nb_jobs+1)] for k in range(nb_jobs+1)]) for i in range(nb_machines)}
for i in range(nb_machines):
    for j in range(nb_jobs+1):
        setup_times[i][j][j] = 0
        setup_times[i][j][0] = 0

##### Optimization #####

# Measure and print solving time for MIP model
print('-----------MIP-----------')
start_time = time.time()
MIP_model(processing_times, setup_times, plotGantt=0)
mip_time = time.time() - start_time

# Measure and print solving time for CP model
print('-----------CP-----------')
start_time = time.time()
CP_model(processing_times, setup_times, plot_gantt=0)
cp_time = time.time() - start_time

print(f"MIP Solving Time: {mip_time} seconds")
print(f"CP Solving Time: {cp_time} seconds")






