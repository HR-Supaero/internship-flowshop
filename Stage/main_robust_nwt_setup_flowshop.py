import numpy as np
import random
import robust_nwt_setup_flowshop as robust_nwt_setup_flowshop

##### Data #####

nom_processing_times = np.array([[random.randint(1,10) for j in range(5)] for i in range(3)])
nb_machines = len(nom_processing_times)
nb_jobs = len(nom_processing_times[0])

max_dev = np.ones(np.shape(nom_processing_times))

setup_times = {i: np.array([[random.randint(1,5) for j in range(nb_jobs+1)] for k in range(nb_jobs+1)]) for i in range(nb_machines)}


#setup_times = {0: np.array([[0, 5, 5, 4], [0, 0, 5, 3], [0, 1, 0, 4], [0, 4, 5, 0]]), 1: np.array([[0, 5, 2, 2], [0, 0, 1, 1], [0, 4, 0, 4], [0, 1, 5, 0]])}
budget = 6

scenario_pool = robust_nwt_setup_flowshop.generate_some_scenariis(10, budget, nb_machines, nb_jobs)
#print(robust_nwt_setup_flowshop.generate_all_scenariis(budget, nb_machines, nb_jobs))


##### Optimizations #####


#scenario_pool = robust_nwt_setup_flowshop.generate_some_scenariis(5, budget, nb_machines, nb_jobs)

#job_sequence, cmax = robust_nwt_setup_flowshop.MIP_model_extended(nom_processing_times, max_dev, setup_times, budget)


#cmax, job_sequence = robust_nwt_setup_flowshop.CP_model_extended(nom_processing_times, max_dev, setup_times, budget)
#master_cmax, job_sequence = robust_nwt_setup_flowshop.MIP_master_problem(scenario_pool, nom_processing_times, max_dev, setup_times)
#print(f'Master sequence makespan {master_cmax} under {job_sequence}')

#adverse_cmax, adverse_scenario = robust_nwt_setup_flowshop.MIP_adverse_problem(job_sequence, nom_processing_times, max_dev, setup_times, budget)
#print(f'Adverse makespan {adverse_cmax}')


#robust_nwt_setup_flowshop.MIP_adverse_problem(job_sequence, nom_processing_times, max_dev, setup_times, budget)

#robust_nwt_setup_flowshop.CP_adverse_problem(job_sequence, nom_processing_times, max_dev, setup_times, budget)

cmax, job_sequence = robust_nwt_setup_flowshop.MIP_MIP_optimization(scenario_pool, nom_processing_times, max_dev, setup_times, budget)

#cmax, job_sequence = robust_nwt_setup_flowshop.CP_CP_optimization(scenario_pool, nom_processing_times, max_dev, setup_times, budget)

cmax, job_sequence = robust_nwt_setup_flowshop.MIP_CP_optimization(scenario_pool, nom_processing_times, max_dev, setup_times, budget)

#cmax, job_sequence = robust_nwt_setup_flowshop.CP_MIP_optimization(scenario_pool, nom_processing_times, max_dev, setup_times, budget)

#print(robust_nwt_setup_flowshop.generate_all_scenariis(budget, nb_machines, nb_jobs))
#print(robust_nwt_setup_flowshop.setup_checker(job_sequence, nom_processing_times, max_dev, setup_times, budget, cmax))

