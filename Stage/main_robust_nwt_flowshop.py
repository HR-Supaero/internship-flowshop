# Optimization of the robust flowshop problem
# Follows an extended model
import numpy as np
from robust_nwt_flowshop import(
    MIP_extended,
    CP_model_extended,
    MIP_MIP_optimization,
    CP_CP_optimization,
    MIP_CP_optimization,
    CP_MIP_optimization,
    check_robustness,
    generate_some_scenariis
)
import random
import time

##### Data #####
nb_machines = 5
nb_jobs = 10
nom_processing_times = np.array([[random.randint(1,10) for j in range(5)] for i in range(3)])

max_dev = np.array([[random.randint(1,5) for j in range(5)] for i in range(3)])
budget = int(np.ceil(nb_machines*nb_jobs*0.2))

nb_machines = len(nom_processing_times)
nb_jobs = len(nom_processing_times[0])

scenario_pool = generate_some_scenariis(1, budget, nb_machines, nb_jobs)
##### Optimization #####

# Measure solving time for each method and print results
def measure_and_print(method_name, method_function, *args):
    print(f'--------{method_name}--------')
    start_time = time.time()
    try:
        cmax, sequence = method_function(*args)
        solving_time = time.time() - start_time
        robustness = check_robustness(nom_processing_times, max_dev, sequence, budget, cmax)
        print(f'Solving Time: {solving_time} seconds')
        print(f'Robustness: {robustness}')
        print(f'Sequence: {sequence}')
        print(f'Makespan: {cmax}')
        print()
        return solving_time
    except Exception as e:
        solving_time = None
        print(f'An error occurred: {e}')
        print()
        return solving_time

# Main function to run and compare all optimization methods
def main():
    solving_times = []
    
    solving_times.append(measure_and_print('MIP', MIP_extended, nom_processing_times, max_dev, budget))
    #solving_times.append(measure_and_print('CP', CP_model_extended, nom_processing_times, max_dev, budget))
    solving_times.append(measure_and_print('MIPMIP', MIP_MIP_optimization, scenario_pool, nom_processing_times, max_dev, budget))
    solving_times.append(measure_and_print('CPMIP', CP_MIP_optimization, scenario_pool, nom_processing_times, max_dev, budget))
    #solving_times.append(measure_and_print('CPCP', CP_CP_optimization, scenario_pool, nom_processing_times, max_dev, budget))
    #solving_times.append(measure_and_print('MIPCP', MIP_CP_optimization, scenario_pool, nom_processing_times, max_dev, budget))

main()



