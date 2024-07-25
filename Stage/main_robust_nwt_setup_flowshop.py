# Optimization of the robust flowshop problem
# Follows an extended model
import numpy as np
import random
import time
import math
from robust_nwt_setup_flowshop import (
    MIP_model_extended,
    CP_model_extended,
    MIP_MIP_optimization,
    CP_CP_optimization,
    MIP_CP_optimization,
    CP_MIP_optimization,
    setup_checker,
    generate_some_scenariis
)

##### Data #####
nb_machines = 5
nb_jobs = 10

nom_processing_times = np.array([[random.randint(1, 10) for j in range(nb_jobs)] for i in range(nb_machines)])
max_dev = np.array([[random.randint(1, 5) for j in range(nb_jobs)] for i in range(nb_machines)])
budget = int(np.ceil(nb_machines * nb_jobs * 0.05))
setup_times = {i: np.array([[random.randint(1, 5) for j in range(nb_jobs + 1)] for k in range(nb_jobs + 1)]) for i in range(nb_machines)}

# Ensure no self setup times
for i in range(nb_machines):
    for j in range(nb_jobs + 1):
        setup_times[i][j][j] = 0
        setup_times[i][j][0] = 0

scenario_pool = generate_some_scenariis(1, budget, nb_machines, nb_jobs)

##### Functions #####

# Measure solving time for each method and print results
def measure_and_print(method_name, method_function, *args):
    print(f'--------{method_name}--------')
    start_time = time.time()
    cmax, job_sequence = method_function(*args)
    solving_time = time.time() - start_time
    robustness = setup_checker(job_sequence, nom_processing_times, max_dev, setup_times, budget, cmax)
    print(f'Solving Time: {solving_time} seconds')
    print(f'Robustness: {robustness}')
    print(f'Sequence: {job_sequence}')
    print(f'Makespan: {cmax}')
    print()
    return solving_time

##### Main Function #####

def main():
    solving_times = []
    
    #solving_times.append(measure_and_print('MIP', MIP_model_extended, nom_processing_times, max_dev, setup_times, budget))
    #solving_times.append(measure_and_print('CP', CP_model_extended, nom_processing_times, max_dev, setup_times, budget))
    #solving_times.append(measure_and_print('MIPMIP', MIP_MIP_optimization, scenario_pool, nom_processing_times, max_dev, setup_times, budget))
    #solving_times.append(measure_and_print('MIPCP', CP_CP_optimization, scenario_pool, nom_processing_times, max_dev, setup_times, budget))
    #solving_times.append(measure_and_print('CPMIP', MIP_CP_optimization, scenario_pool, nom_processing_times, max_dev, setup_times, budget))
    solving_times.append(measure_and_print('CPCP', CP_MIP_optimization, scenario_pool, nom_processing_times, max_dev, setup_times, budget))
    

main()

