import itertools
import robust_nwt_rsetup_flowshop as robust_nwt_rsetup_flowshop
import numpy as np
import random

##### Data #####
nb_machines = 2
nb_jobs = 3

nom_processing_times = np.array([[random.randint(1,2) for j in range(nb_jobs)] for i in range(nb_machines)])


max_dev = np.zeros(np.shape(nom_processing_times))

setup_times = {i: np.array([[random.randint(1,1) for j in range(nb_jobs+1)] for k in range(nb_jobs+1)]) for i in range(nb_machines)}
for i in range(nb_machines):
    for j in range(nb_jobs+1):
        setup_times[i][j][j] = 0
        setup_times[i][j][0] = 0

setup_dev = setup_times

#setup_times = {0: np.array([[0, 5, 5, 4], [0, 0, 5, 3], [0, 1, 0, 4], [0, 4, 5, 0]]), 1: np.array([[0, 5, 2, 2], [0, 0, 1, 1], [0, 4, 0, 4], [0, 1, 5, 0]])}
process_budget = 1
setup_budget = 2
cmax, job_sequence = robust_nwt_rsetup_flowshop.MIP_model_extended(nom_processing_times, max_dev, setup_times, setup_dev, process_budget, setup_budget)
print(robust_nwt_rsetup_flowshop.checker(cmax, job_sequence, nom_processing_times, max_dev, setup_times, setup_dev, process_budget, setup_budget))
