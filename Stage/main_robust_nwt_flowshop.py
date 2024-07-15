# Optimization of the robust flowshop problem
# Follows an extended model
import numpy as np
import robust_nwt_flowshop as robust_nwt_flowshop
import random

##### Data #####

nom_processing_times = np.array([[random.randint(1,10) for j in range(5)] for i in range(3)])

max_dev = nom_processing_times
budget = 2


nb_machines = len(nom_processing_times)
nb_jobs = len(nom_processing_times[0])

scenario_pool = robust_nwt_flowshop.generate_some_scenariis(4, budget, nb_machines, nb_jobs)
##### Optimization #####

print('--------MIP--------')
model, sequence, cmax = robust_nwt_flowshop.MIP_extended(nom_processing_times, max_dev, budget)
print(robust_nwt_flowshop.check_robustness(nom_processing_times, max_dev, sequence, budget, cmax))

print('--------CP--------')
#cmax, sequence = robust_nwt_flowshop.CP_model_extended(nom_processing_times, max_dev, budget)
print(robust_nwt_flowshop.check_robustness(nom_processing_times, max_dev, sequence, budget, cmax))

print('--------MIPMIP--------')
cmax, sequence = robust_nwt_flowshop.MIP_MIP_optimization(scenario_pool, nom_processing_times, max_dev, budget)
print(robust_nwt_flowshop.check_robustness(nom_processing_times, max_dev, sequence, budget, cmax))

print('--------CPCP--------')
#cmax, sequence = robust_nwt_flowshop.CP_CP_optimization(scenario_pool, nom_processing_times, max_dev, budget)
print(robust_nwt_flowshop.check_robustness(nom_processing_times, max_dev, sequence, budget, cmax))

print('--------MIPCP--------')
cmax, sequence = robust_nwt_flowshop.MIP_CP_optimization(scenario_pool, nom_processing_times, max_dev, budget)
print(robust_nwt_flowshop.check_robustness(nom_processing_times, max_dev, sequence, budget, cmax))

print('--------CPMIP--------')
cmax, sequence = robust_nwt_flowshop.CP_MIP_optimization(scenario_pool, nom_processing_times, max_dev, budget)
print(robust_nwt_flowshop.check_robustness(nom_processing_times, max_dev, sequence, budget, cmax))


