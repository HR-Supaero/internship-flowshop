# Runs no wait flowshop optimization
import numpy as np

import random

from nowait_flowshop import MIP_model, CP_model

##### Data #####

# Small instance
#processing_times = np.array([[random.randint(1, 10) for j in range(3)] for i in range(2)])

# Medium instance
#processing_times = np.array([[random.randint(1, 10) for j in range(6)] for i in range(4)])

# Large instance
processing_times = np.array([[random.randint(1, 10) for j in range(10)] for i in range(5)])

##### Optimization #####

print('-----------MIP-----------')
MIP_model(processing_times, plotGantt=1)

print('-----------CP-----------')
CP_model(processing_times, plot_gantt=1)







