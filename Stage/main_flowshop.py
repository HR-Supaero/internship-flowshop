# Runs classic flowshop optimization

import random

import numpy as np

from flowshop import MIP_model

##### Data #####

# Small instance
#processing_times = np.array([[random.randint(1, 10) for j in range(3)] for i in range(2)])

# Medium instance
#processing_times = np.array([[random.randint(1, 10) for j in range(6)] for i in range(4)])

# Large instance
processing_times = np.array([[random.randint(1, 10) for j in range(10)] for i in range(5)])


##### Run optimizations #####

MIP_model(processing_times, 1)



