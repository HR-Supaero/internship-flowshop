import random
import gurobipy as gp
import matplotlib
import numpy as np
import flowshop as flowshop
from itertools import combinations
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import this
#Processing time for each job on each machine




#Setup time for each machine and for each job and for the job preceding it
#The first index corresponds to the preceeding job, the second represents the current job, the last, the machine
#The index 0 for i and l correspond to the dummy job

d = {
    'spam': 'ham',
    'knights': 'lumberjack',
    'bacon': 'sausage'
}



pt = np.array([[3., 5., 5.],
 [2., 4., 4.],
 [4., 4., 4.],
 [3., 3., 6.]])
#Number of jobs (w/o the dummy job)

#flowshop.standardFlowshop(nom_pt, 1)

A = np.array([[1, 2, 3],
              [1, 2, 3]])

B = np.array([[3, 1, 2],
              [3, 1, 2]])



print(np.array([1, 2]) + np.array([1, 2]))


