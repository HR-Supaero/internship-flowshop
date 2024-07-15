#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
from docplex.cp.model import *
from docplex.cp.model import CpoModel

import numpy as np

from SP_fct import *

# In[2]:



#context.solver.local.execfile="/Users/lhoussin/CPLEX_Studio2211/cpoptimizer/bin/x86-64_osx/cpoptimizer"


# In[4]:


# data
nb_jobs=4
nb_machines=3

sequence=[1,0,2,3]

budget=4

#P[i,j] processing time of job i on machine j

P= np.matrix([[3, 2,2 ], [2, 1,3], [2, 3,3],[1, 4,1]])
SetUp= np.matrix([[2,3,1 ], [3,3,4], [3,1,2],[1, 1,1]])
devP= np.matrix([[3, 2,2 ], [2, 1,3], [2, 3,3],[1, 4,1]])

print(P)
print(SetUp)
print(P[0,0])


[a,b]=SP(sequence,2,nb_jobs, nb_machines, P,SetUp,devP,1)
    
print(a)
print(b)
    
    



