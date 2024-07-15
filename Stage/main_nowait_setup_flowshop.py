# This code will run an optimization on a no wait flowshop with setup times

import numpy as np
import random
import nowait_setup_flowshop as nowait_setup_flowshop

##### Data #####

processing_times = np.array([[3, 5, 3], [5, 4, 1]])
nb_jobs = len(processing_times[0])
nb_machines = len(processing_times)

setup_times = {i: np.array([[random.randint(1,5) for j in range(nb_jobs+1)] for k in range(nb_jobs+1)]) for i in range(nb_machines)}
#setup_times = {0: np.array([[0, 2, 2, 3],[0, 0, 2, 2],[0, 3, 0, 5],[0, 1, 4, 0]]), 1: np.array([[0, 1, 3, 4],[0, 0, 2, 2],[0, 5, 0, 5],[0, 3, 4, 0]])}
#setup_times = {0: np.array([[1, 3, 3, 2], [1, 4, 3, 4], [5, 1, 1, 5], [5, 2, 3, 1]]), 1: np.array([[2, 2, 3, 2], [5, 4, 3, 2], [5, 3, 2, 2], [1, 3, 1, 1]])}


##### Optimization #####
model, processing_times, idle_times, setup_times = nowait_setup_flowshop.MIP_model(processing_times, setup_times, 1)
print(nowait_setup_flowshop.cmax(processing_times, setup_times))

#nowait_setup_flowshop.CP_model(processing_times, setup_times, 1)





