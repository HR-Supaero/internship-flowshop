# This code implements MIP model for a no wait flowshop with
# sequence dependent setup times

import gurobipy as gp
from gurobipy import GRB
import numpy as np
import numpy.linalg as alg
from docplex.cp.model import *
import flowshop_gantt

##### Helper functions ######

def reorder_sequence(
        job_sequence: np.ndarray, 
        mat: np.ndarray
) -> np.ndarray:
    """
     This function reorders the columns of a matrix along a job sequence
     
     job_sequence: the order of jobs in the sequence

     mat: the matrix that we want to reorder
    
     returns: the reorganised matrix

    """
    nb_jobs = len(mat[0])
    permu_mat = np.zeros((nb_jobs, nb_jobs))
    for j in range(nb_jobs):
        permu_mat[j, job_sequence[j]] = 1
    permu_mat = alg.inv(permu_mat)

    return np.dot(mat, permu_mat)

def cmax(
        processing_times,
        setup_times
):
    """
    This function calculates the makespan of the no wait flowhshop
    with sequence dependent setups
    WARNING: the setup matrix is the setup times corresponding to the given
    sequence of jobs !
    """
    nb_jobs = len(processing_times[0])
    nb_machines = len(processing_times)
    def completion_difference(
            processing_times, 
            setup_times, 
            sequence_index
    ):
        max_list = []
        if sequence_index == 0:
            for i in range(nb_machines):
                max_list.append(sum(processing_times[l][sequence_index] for l in range(i, nb_machines)) + setup_times[i][sequence_index])
            return max(max_list)
        for i in range(nb_machines):
            max_list.append(sum(processing_times[l][sequence_index] - processing_times[l][sequence_index-1] for l in range(i, nb_machines)) + processing_times[i][sequence_index-1] + setup_times[i][sequence_index])
        return max(max_list)
    
    makespan = 0
    for k in range(nb_jobs):
        makespan += completion_difference(processing_times, setup_times, k)
    
    return makespan

##### Models #####

def MIP_model(
        processing_times: np.ndarray, 
        setup_times: dict[tuple, int], 
        plotGantt: int
) -> tuple[gp.Model, np.ndarray, np.ndarray, np.ndarray]:
    """
     This function takes a matrix representing production times
     and a dictionnary representing the setup times
     and returns an optimized model of the no wait setup flowshop
    
     processing_times: the processing times time matrix, 
    
     setup_times: the setup times dictionnary
           
     plot_gantt: Define this value to 1 if wishing to plot the Gantt chart
    
     returns: an optimized gurobipy model 
              the reordered processing times along the optimal sequence
              the idle times matrix
              the setup times matrix for the optimal sequence
    
    """

    if type(processing_times) != np.ndarray:
        raise TypeError("The processing times format is incorrect")
    
    if type(setup_times) != dict:
        raise TypeError("The setup times format is incorrect")
    

    ##### Data #####

    nb_jobs = len(processing_times[0])      # dummy job is not counted
    nb_machines = len(processing_times)
    A = 100000                              # arbitrary large number

    
            
    ##### Initializing the model #####

    m = gp.Model('Flowshop_wsetup')
    m.Params.IntegralityFocus = 1
    m.Params.OutputFlag = 0

    ##### Variables #####

    # idle[i, j, k] idle time on machine i 
    # between the end of job j and the start of job k 
    # assuming k follows j in the sequence
    # j = 0 corresponds to the dummy job
    idle = m.addVars(nb_machines, nb_jobs + 1, nb_jobs + 1, lb = 0.0, ub = GRB.INFINITY, name = "idle")

    # difference[j, k] time difference between 
    # the completion times of job j and job k on the last machine
    # assuming k follows j in the sequence
    # j = 0 corresponds to the dummy job
    difference = m.addVars(nb_jobs + 1, nb_jobs + 1, lb = 0.0, ub = GRB.INFINITY, name = 'difference')

    # u[j] corresponds to the position of job j in the sequence
    # j = 0 corresponds to the dummy job
    u = m.addVars(nb_jobs + 1, vtype = GRB.INTEGER, name = 'u')

    #z[j,k] equals 1 if job k follows job j in the sequence, 0 otherwise
    z = m.addVars(nb_jobs + 1, nb_jobs + 1, vtype = GRB.BINARY, name = 'z')


    ##### Constraints #####

    # Edge cases
    # The idle time, completion time difference, following decision variables
    # cant be defined if the following job is the dummy job
    # or if the current and following job are the same

    m.addConstrs(idle[i, j, 0] == 0 
                 for i in range(nb_machines) 
                 for j in range(nb_jobs+1))
    
    m.addConstrs(idle[i, j, j] == 0 
                 for i in range(nb_machines) 
                 for j in range(nb_jobs+1))

    m.addConstrs(difference[j,j] == 0 
                 for j in range(nb_jobs+1))

    m.addConstrs(difference[j, 0] == 0 
                 for j in range(nb_jobs+1))

    m.addConstrs(z[j,j] == 0 
                 for j in range(nb_jobs+1))

    m.addConstrs(z[j,0] == 0 
                 for j in range(nb_jobs+1))


    # The idle times between jobs has to be bigger than the setup times
    m.addConstrs(idle[i,j,k] >= setup_times[i][j][k] 
                 for i in range(nb_machines) 
                 for j in range(nb_jobs+1) 
                 for k in range(nb_jobs+1))

    # JAML constraints
    m.addConstrs(processing_times[i][j-1] 
                + idle[i,j,k] 
                == 
                 processing_times[i-1][k-1] 
                + idle[i-1,j,k] 
                 for i in range(1, nb_machines) 
                 for j in range(1, nb_jobs + 1) 
                 for k in range(1, nb_jobs+1) 
                 if j!=k)

    # JAML constraints for the first job in the sequence
    m.addConstrs(idle[i,0,k] == 
                 idle[i-1,0,k] 
                 + processing_times[i-1][k-1] 
                 for i in range(1, nb_machines) 
                 for k in range(1, nb_jobs+1))

    # completion time difference > 0 
    # if and only if j is followed by k, otherwise its null
    m.addConstrs(difference[j,k] >= idle[nb_machines-1,j,k] 
                                    + processing_times[nb_machines-1][k-1] 
                                    - A*(1-z[j,k]) 
                                    for j in range(nb_jobs + 1) 
                                    for k in range(1, nb_jobs + 1))

    # Each job can only have one job after it
    m.addConstrs(gp.quicksum(z[j,k] for k in range(nb_jobs + 1)) <= 1 
                 for j in range(nb_jobs+1))
    # Each job has one job preceding it apart from the dummy job
    m.addConstrs(gp.quicksum(z[j,k] for j in range(nb_jobs + 1)) >= 1 
                 for k in range(1,nb_jobs+1))

    # Defines the job position variable
    m.addConstrs(u[j] - u[k] + nb_jobs*z[j,k] + (nb_jobs-2)*z[k,j] <= nb_jobs-1 
                 for j in range(1, nb_jobs+1) 
                 for k in range(1, nb_jobs+1))
    m.addConstr(u[0] == 0)
    m.addConstrs(u[j] >= 1 for j in range(1,nb_jobs+1))
    m.addConstrs(u[j] <= nb_jobs for j in range(1,nb_jobs+1))


    ##### Objective #####

    m.setObjective(
        gp.quicksum(difference[j,k] 
                    for j in range(nb_jobs+1) 
                    for k in range(nb_jobs+1)), 
                    GRB.MINIMIZE )

    ##### Optimization #####
    m.optimize()

    ##### Retrieve data #####
    for v in m.getVars():
        if v.X != 0:
            print(v)

    # Retrieve the optimal sequence and the makespan
    opt_sequence = []
    for v in m.getVars():
        if v.VarName[0] == 'u':
            opt_sequence.append((v.VarName, int(v.X)))
    opt_sequence = sorted(opt_sequence, key=lambda elem: elem[1])
    def retrieve_sequencePos(varName):
        return int(varName[0].split('[')[1].split(']')[0])
    for k in range(len(opt_sequence)):
        opt_sequence[k] = retrieve_sequencePos(opt_sequence[k])
    print(f'one optimal sequence is {opt_sequence}')

    job_sequence = opt_sequence.copy()         # Will contain the optimal sequence w/o the dummy job
    for k in range(len(job_sequence)):
        job_sequence[k] = job_sequence[k] - 1  # job indexes in opt_sequence are offseted by 1 due to the dummy job in index 0
    job_sequence.remove(-1)
    
    makespan = m.getObjective().getValue()
    if not makespan.is_integer():
        raise ValueError("The makespan is not of integer value")
    print(f'The makespan is {makespan}')

    # Reordering the processing times along the given optimal sequence
    print(f'given processing times:\n', processing_times)
    processing_times = reorder_sequence(job_sequence, processing_times)
    print(f'reordered processing times\n', processing_times)

    # Retrieving the idle times
    idle_mat= np.zeros(np.shape(processing_times))
    for i in range(nb_machines):      # The values in opt_sequence are offseted by 1 
        for j in range(nb_jobs):       
            idle_mat[i,j] = idle[i, opt_sequence[j], opt_sequence[j+1]].X
    print(f'Idle times:\n', idle_mat)

    # Retrieving the setup times
    setup_mat = np.zeros(np.shape(processing_times))
    for i in range(nb_machines):     
        for j in range(nb_jobs):       
            setup_mat[i,j] = setup_times[i][opt_sequence[j]][opt_sequence[j+1]]
    print(f'setup_times times:\n', setup_mat)

    if plotGantt:
        flowshop_gantt.gantt_setup(processing_times, setup_mat, idle_mat, makespan)

    return m, processing_times, idle_mat, setup_mat

def CP_model(processing_times, setup_times, plot_gantt):


    ##### Data #####
    nb_machines = len(processing_times)
    nb_jobs = len(processing_times[0])

    dummy_processing_times = np.zeros((nb_machines, 1))
    processing_times = np.hstack((dummy_processing_times, processing_times))

    for i in range(nb_machines):
        for j in range(nb_jobs):
             setup_times[i][j][j] = 0
             setup_times[i][j][0] = 0

    ##### Initialize #####
    mdl = CpoModel()

    ##### Variables #####

    tasks = [[interval_var(size=int(processing_times[i][j]), name='task_{}_{}'.format(i,j)) for j in range(nb_jobs+1)] for i in range(nb_machines)]

    sequence = [sequence_var(tasks[i], name='sequence_{}'.format(i)) for i in range(nb_machines)]

    setup_mat = [transition_matrix(setup_times[i]) for i in range(nb_machines)]

    Cmax = integer_var(name='Cmax')

    ##### Constraints #####

    mdl.add(start_of(tasks[i][0]) == 0 for i in range(nb_machines))

    mdl.add(greater_or_equal(Cmax, end_of(tasks[nb_machines-1][j])) for j in range(nb_jobs+1))

    mdl.add(end_at_start(tasks[i][j], tasks[i+1][j]) for i in range(nb_machines-1) for j in range(nb_jobs+1))

    mdl.add(no_overlap(sequence[i], setup_mat[i]) for i in range(nb_machines))

    mdl.add(same_sequence(sequence[0], sequence[i]) for i in range(1, nb_machines))

    ##### Objective #####

    mdl.minimize(Cmax)

    ##### Solve ! #####

    msol = mdl.solve()


    ##### Retrieve the solutions #####
    makespan = msol.get_value(Cmax)

    seqc_var = msol.get_value(sequence[1])
    def retrieve_sequence_position(var_name):
        return int(var_name.split('_')[2])
    opt_sequence = []
    for s in seqc_var:
        opt_sequence.append(retrieve_sequence_position(s.get_name()))
        
    job_sequence = opt_sequence.copy()         # Will contain the optimal sequence w/o the dummy job
    for k in range(len(job_sequence)):
        job_sequence[k] = job_sequence[k] - 1  # job indexes in opt_sequence are offseted by 1 due to the dummy job in index 0
    job_sequence.remove(-1)
    print(f'one optimal sequence is: {job_sequence}')

    opt_sequence.remove(0)
    opt_sequence.insert(0, 0)
    print(f'opt sequence w: dummy: {opt_sequence}')
    
    slvd_tasks = [[msol.get_var_solution(tasks[i][j]) for j in range(nb_jobs+1)] for i in range(nb_machines)]
    for i in range(nb_machines):
        for j in range(nb_jobs+1):
            print(f'before {slvd_tasks[i][j].get_name()} start: {slvd_tasks[i][j].get_start()}, size = {slvd_tasks[i][j].get_size()}')
    for i in range(nb_machines):
        del slvd_tasks[i][0]
        
    slvd_setup_mat = np.zeros((nb_machines, nb_jobs))
    for i in range(nb_machines):
        for j in range(nb_jobs):
            slvd_setup_mat[i][j] = setup_times[i][opt_sequence[j]][opt_sequence[j+1]]
    print(f'The setup times are: \n {slvd_setup_mat}')

    if plot_gantt:
        flowshop_gantt.ganttCP_setup(slvd_tasks, slvd_setup_mat, makespan)



