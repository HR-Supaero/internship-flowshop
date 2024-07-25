#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
from docplex.cp.model import *
from docplex.cp.model import CpoModel

import numpy as np



#from docplex.mp.check_list import run_docplex_check_list
#run_docplex_check_list()

#context.solver.local.execfile="/Users/lhoussin/CPLEX_Studio2211/cpoptimizer/bin/x86-64_osx/cpoptimizer"



def SP(sequence,budget,nb_jobs,nb_machines,P,SetUp,devP,visu):

    mdl = CpoModel()

    Cmax = mdl.integer_var(name = "Cmax") # makespan
    tasks = [[mdl.interval_var(name='tasks_{}{}'.format(i,j), size=int(P[i,j])) for j in range(nb_machines)] for i in range(nb_jobs)]
    SetUptasks = [[mdl.interval_var(name='SetUptasks_{}{}'.format(i,j), size=int(SetUp[i,j])) for j in range(nb_machines)] for i in range(nb_jobs)]
    DeviationTask = [[mdl.interval_var(name='deviations_{}{}'.format(i,j), size=int(devP[i,j]), optional=True) for j in range(nb_machines)] for i in range(nb_jobs)]

       
    ev = mdl.integer_var_list(nb_jobs,name = "ev") # slack



    # Every SetUptasks of first job is scheduled at 0
    # $$ StartOf(SetUpTasks_{0,m})=0 \quad \forall m\in \mathcal{M}$$

    # In[19]:


    # Every SetUptasks of first job is scheduled at 0
    for m in range( nb_machines):
        mdl.add(mdl.start_of(SetUptasks[sequence[0]][m])==0)


    # A SetUptasks is always scheduled before the corresponding task
    # $$ EndBeforeStart(SetUpTasks_{j,m}, Tasks_{j,m})  \quad \forall m\in \mathcal{M}$$

    # In[7]:


    #  precedence constraints SetUp -> task
    for j in range(nb_jobs):
        for m in range( nb_machines):
            mdl.add(mdl.end_before_start(SetUptasks[sequence[j]][m ], tasks[sequence[j]][m]))


    # The deviation of task occurrs at the end of the task
    # $$ StartAtEnd(DeviationTask_{j,m},Tasks_{j,m}) \quad \forall m\in \mathcal{M},\forall j\in \mathcal{J}$$

    # In[8]:


    #  Must stick constraints task -> Devtask
    for j in range(nb_jobs):
        for m in range( nb_machines):
            mdl.add(mdl.start_at_end(DeviationTask[sequence[j]][m], tasks[sequence[j]][m]))


    # A SetUptask(j,m) starts at the end the task(j-1,m) and  task(j-1,m) can be extended (deviation)
    # $$ StartOf(SetUpTasks_{j,m})=\max(EndOf(Tasks_{j-1,m}),EndOf(DeviationTask{j-1,m})) \quad \forall m\in \mathcal{M},\forall j\in \mathcal{J}$$

    # In[20]:


    # precedence  constraints  setUp(i).start = max( task(i-1).end, DeviationTask(i-1).end)
    for j in range(1,nb_jobs):
        for m in range( nb_machines):
            mdl.add(mdl.start_of(SetUptasks[sequence[j]][m]) == mdl.max(mdl.end_of(tasks[sequence[j-1]][m]),mdl.end_of(DeviationTask[sequence[j-1]][m])))
                


    # No wait precedence constraint
    # $$ StartOf(Tasks_{j,m})= \max(EndOf(Tasks_{j,m-1}),EndOf(DeviationTask{j,m-1})) \quad \forall m\in \mathcal{M},\forall j\in \mathcal{J}$$

    # In[21]:


    # NoWait precedence constraints
    for j in range(nb_jobs):
        for m in range(1, nb_machines):
            mdl.add(mdl.start_of(tasks[sequence[j]][m ]) == mdl.max(mdl.end_of(tasks[sequence[j]][m-1]) , mdl.end_of(DeviationTask[sequence[j]][m-1]) ))


    # Linking jobs with no slack at least one time
    # $$ EV_{j}= \min_{m\in\mathcal{M}}(StartOf(Tasks_{j,m})-EndOf(SetUpTasks{j,m})) \quad
    # \forall j\in \mathcal{J}$$
    # $$\max_{j\in \mathcal{J}} EV_{j}=0$$

    # In[11]:


    # Alternative way to link jobs with  slack=0 at least one time
    for i in range(nb_jobs):
        mdl.add(ev[i]==mdl.min([mdl.start_of(tasks[sequence[i]][m])-mdl.end_of(SetUptasks[sequence[i]][m]) for m in range(nb_machines)]))

    mdl.add(mdl.max([ev[i] for i in range(nb_jobs)])==0)


    # Budget constraint
    # $$ \sum_{m\in\mathcal{M}}\sum_{j\in\mathcal{J}} PresenceOf(DeviationTask_{j,m})\leq budget$$

    # In[12]:


    # Budget constraints method 3

    mdl.add(mdl.sum([mdl.presence_of(DeviationTask[sequence[j]][m])  for j in range(nb_jobs) for m in range(nb_machines)]) <= budget)
           # D = mdl.sum([s.level * mdl.presence_of(wtasks[(h, s)]) for s in SKILLS for h in HOUSES])


    # Makespan expression (M= last machine)
    #
    # $$ C_{max} = \max(\max_{j\in\mathcal{J}} (EndOf(Tasks_{j,M})), \max_{j\in\mathcal{J}} (EndOf(DeviationTasks_{j,M})) $$

    # In[22]:


    # Cmax contraint           TODO: OK
    mdl.add(Cmax == mdl.max([mdl.end_of(tasks[i][nb_machines-1]) for i in range(nb_jobs)]+ [mdl.end_of(DeviationTask[i][nb_machines-1]) for i in range(nb_jobs)] ))


    # In[14]:


    # creation de l'objectif
    mdl.maximize(Cmax)


    # In[15]:


    # Solve the model
    print("\nSolving model....")
    msol = mdl.solve(TimeLimit=10,LogVerbosity = 'Quiet')
    print("done")


    # In[24]:


    objective_SP = msol.get_value(Cmax) # on recupere la valeur optimale du makespan

    #sol_tasks = [[msol.get_value(tasks[i][j]) for j in range(nb_machines)] for i in range(nb_jobs)]
    var_sol = [[msol.get_var_solution(tasks[i][j]) for j in range(nb_machines)] for i in range(nb_jobs)]
    var_sol_SetUp = [[msol.get_var_solution(SetUptasks[i][j]) for j in range(nb_machines)] for i in range(nb_jobs)]
    var_sol_Dev = [[msol.get_var_solution(DeviationTask[i][j]) for j in range(nb_machines)] for i in range(nb_jobs)]
    #var_sol_deviation_nb= msol.get_value(cumul_deviation)
    ev_sol= [msol.get_value(ev[i]) for i in range(nb_jobs-1)]

    print(ev_sol)
    print("Cmax=",objective_SP)
    #print(var_sol_deviation_nb)
    print("----")

    #for i in range(nb_jobs):
    #    for j in range(nb_machines):
    #        print("i,j", i , j)
            #print(var_sol_SetUp[i][j].get_start())
            #print(var_sol[i][j].get_start())
            #print(var_sol_Dev[i][j].get_start())
            
    #print(ev_sol)

    xi=np.zeros((nb_jobs,nb_machines))

    for i in range(nb_jobs):
        for j in range(nb_machines):
            if var_sol_Dev[i][j].get_start():
                xi[i,j]=1
    print("xi=",xi)

    # In[26]:

    if visu==1:
        ########################## GANTT VISU ###########################################
        #https://www.geeksforgeeks.org/python-basic-gantt-chart-using-matplotlib/

        # Importing the matplotlib.pyplot
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors

        #print(mcolors.TABLEAU_COLORS)
        color_names = list(mcolors.TABLEAU_COLORS)

        plt.rcParams["figure.figsize"] = (14,3)
        plt.rcParams["figure.dpi"] = 200
         
        # Declaring a figure "gnt"
        fig, gnt = plt.subplots()

        # Setting Y-axis limits
        gnt.set_ylim(0, 5)
         
        # Setting X-axis limits
        gnt.set_xlim(0, objective_SP+5)
         
        # Setting labels for x-axis and y-axis
        gnt.set_xlabel('Time')
        gnt.set_ylabel('Machine')
        #gnt.legend(title='alpha='+str(alpha))
         
        # Setting ticks on y-axis
        gnt.set_yticks([1, 2, 3, 4])
        # Labelling tickes of y-axis
        gnt.set_yticklabels(['M4', 'M3', 'M2', 'M1'])

        grid_x_ticks = np.arange(0, objective_SP+5, 1)
        gnt.set_xticks(grid_x_ticks)
        # Setting graph attribute
        gnt.grid(axis='x',alpha=0.3)
        #plt.grid(axis='x', color='0.95')





        for i in range(nb_jobs):

            for j in range(nb_machines):


                    gnt.broken_barh([(var_sol[i][j].get_start(), var_sol[i][j].get_size())], (3-j+1-0.2, 0.4),facecolors =(color_names[i]), edgecolor=('black'))
                    print(f'{var_sol[i][j].get_start()} woooooola')
                    gnt.text( var_sol[i][j].get_start() +  P[i,j]/2 , 4-j, str(i), ha='center', va='center')
                    #print(sequence[i])
                    #print(var_sol[i][j].get_start() , P[i,j])
                    
                    gnt.broken_barh([(var_sol_SetUp[i][j].get_start(), var_sol_SetUp[i][j].get_size())], (3-j+1-0.2, 0.4),facecolors =(color_names[i]), edgecolor=('red'),linewidth=1, hatch='///',)
                    if var_sol_Dev[i][j].get_start():
                        gnt.broken_barh([(var_sol_Dev[i][j].get_start(), var_sol_Dev[i][j].get_size())], (3-j+1-0.2, 0.4),facecolors =(color_names[i]), edgecolor=('cyan'),linewidth=1, hatch='ooo',)
                    

        plt.savefig("gantt1.png")
        #plt.savefig("gantt1.pdf", format='pdf', dpi=1200)

        #plt.show()

    return(xi,objective_SP)
    



