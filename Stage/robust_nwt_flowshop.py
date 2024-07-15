
import gurobipy as gp
from gurobipy import GRB

import flowshop_gantt
import nowait_flowshop as nowait_flowshop

from itertools import combinations, permutations

import numpy as np
import numpy.linalg as alg

import random
import math

from docplex.cp.model import *
from docplex.cp.model import CpoModel

############################
##### Helper functions #####
############################

def reorder_sequence(
        job_sequence: np.ndarray, 
        mat: np.ndarray
) -> np.ndarray:
    """Reorders the columnns of a matrix along a given sequence"""
    return mat[:, job_sequence]

def generate_some_scenariis(
        nb_scenariis: int,
        budget: int, 
        nb_machines: int, 
        nb_jobs: int
) -> list[np.ndarray]:
    """Generates a list of the desired number of scenarios"""

    list_scen = []
    
    for gamma in range(budget+1):
        possible_ones = list(combinations(range(nb_machines*nb_jobs), gamma))
        for indx_comb in possible_ones:
            scenario = np.zeros((nb_machines, nb_jobs))
            for index in indx_comb:
                row = index // nb_jobs
                col = index % nb_jobs
                scenario[row, col] = 1
            list_scen.append(scenario.copy())

    sampled_scenariis = random.sample(list_scen, nb_scenariis)
    return sampled_scenariis

def generate_all_scenariis(
        budget: int, 
        nb_machines: int, 
        nb_jobs: int
) -> list[np.ndarray]:
    """This code returns an array of all possible sceanriis"""

    list_scen = []
    
    for gamma in range(budget+1):
        possible_ones = list(combinations(range(nb_machines*nb_jobs), gamma))
        for indx_comb in possible_ones:
            scenario = np.zeros((nb_machines, nb_jobs))
            for index in indx_comb:
                row = index // nb_jobs
                col = index % nb_jobs
                scenario[row, col] = 1
            list_scen.append(scenario.copy())

    return list_scen

def calculate_processing_times(
        scenario: np.ndarray, 
        nom_production_times: np.ndarray, 
        max_deviations: np.ndarray
) -> np.ndarray:
    """returns the production time matrix in a certain scenario"""
    return nom_production_times + np.multiply(scenario, max_deviations)

def generate_all_processing_times(
        nom_processing_times, 
        max_dev, 
        scenario_pool
) -> list[np.ndarray]:
    """Generates every possible processing times under some scenario pool"""
    processing_times = []       
    for scen in scenario_pool:
        pt_s = calculate_processing_times(scen, nom_processing_times, max_dev)
        processing_times.append(pt_s.copy())
    return processing_times

def check_robustness(
        nom_processing_times: np.ndarray,
        max_dev: np.ndarray, 
        job_sequence: np.ndarray, 
        budget: int, 
        worst_makespan: int
) -> bool:
    
    """
    Checks if under a given sequence, 
    all scenarios lead to a makespan lower than worst case makespan.
    Then checks if for other sequences, there exists a longer makespan
    """
    
    ##### Data #####
        

    nb_machines = len(nom_processing_times)
    nb_jobs = len(nom_processing_times[0])
    scenariis = generate_all_scenariis(budget, nb_machines, nb_jobs)
    
    def evaluate_sceanriis(scenariis, nom_processing_times, max_dev, worst_makespan, job_sequence):
        """
        For the given sequence, check for all given scenarios 
        if makespans are lower than the worst case makespan
        """

        nom_processing_times = reorder_sequence(job_sequence, nom_processing_times)
        max_dev = reorder_sequence(job_sequence, max_dev)

        for scen in scenariis:
            scen = reorder_sequence(job_sequence, scen)
            processing_times = calculate_processing_times(scen, nom_processing_times, max_dev)
            makespan = nowait_flowshop.cmax(processing_times)
            #print(f'The makespan is {Cmax} for scenario \n {scen}')
            if makespan > worst_makespan:
                print(f'The scenario \n {scen} \n violates the worst case makespan with a makespan of {makespan}')
                return False    
    
    def find_worst_makespan(opt_sequence, nom_processing_times, max_dev, worst_makespan):
        """ For other sequences, find a greater makespan """
        other_sequences = list(permutations(range(nb_jobs)))
        other_sequences.remove(tuple(opt_sequence))

        for seq in other_sequences:

            nom_processing_times = reorder_sequence(seq, nom_processing_times)
            max_deviations = reorder_sequence(seq, max_deviations)
            longer_makespan_found = False
            worst_makespan_for_seq = 0

            for scen in scenariis:
                scen = reorder_sequence(seq, scen)
                processing_times = calculate_processing_times(scen, nom_processing_times, max_deviations)
                makespan = nowait_flowshop.cmax(processing_times)
                if makespan >= worst_makespan:
                    longer_makespan_found = True
                    worst_makespan_for_seq = makespan
                    break
            if not longer_makespan_found:
                print(f'All makespans for sequence {seq} are shorter to the worst makespan.')
            return False
        
        
        if not evaluate_sceanriis(scenariis, nom_processing_times, max_dev, worst_makespan, job_sequence):
            return False
        
        if not find_worst_makespan(job_sequence, nom_processing_times, max_dev, worst_makespan):
            return False
        
    return True
    
def is_scenario_in_pool(
        scenario: np.ndarray, 
        scenario_pool: list[np.ndarray]
    ) -> bool:
    """Check if a scenario is present in a scenario pool"""
    return any(np.array_equal(scenario, scen) for scen in scenario_pool)

def retrieve_opt_seq_MIP(
        slvd_model: gp.Model
    ) -> list:
    """retrieve the optimal sequence of the model"""

    var_list = []
    for v in slvd_model.getVars():
        if v.X != 0:
            if (v.VarName[0] == "x"):
                var_list.append(v.VarName)
    nb_jobs = len(var_list)

    def retrive_sequence_position(var_name):
        return int(var_name.split(',')[1].split(']')[0])
    
    def retrieve_job_number(var_name):
       return int(var_name.split(',')[0].split('[')[1])
    
    opt_sequence = sorted(var_list, key = retrive_sequence_position)
    for k in range(nb_jobs):
        opt_sequence[k] = retrieve_job_number(opt_sequence[k])

    return opt_sequence

def retrieve opt_seq_CP()
##################
##### Models #####
##################

def MIP_extended(
        nom_processing_times: np.ndarray, 
        max_dev: np.ndarray, 
        budget: int
) -> tuple[gp.Model, np.ndarray, int]:
    """
    This function optimizes a robust flowshop model following an extended
    MIP formulation.
    Meaning all possible scenarios under the budget are generated    

    Parameters:

    nom_processing_times: the nominal processing times time matrix, 
    
    max_deviations: the maximum deviation values matrix
    
    budget: the uncertainty budget
    
    returns: an optimized gurobipy model\n 
             one optimal sequence\n
             the worst makespan under the optimal sequence
    
    """

    ##### Data #####

    nb_machines = len(nom_processing_times)
    nb_jobs = len(nom_processing_times[0])
    
    scenariis = generate_all_scenariis(budget, nb_machines, nb_jobs)
    nb_scenariis = len(scenariis)
 
    processing_times = generate_all_processing_times(nom_processing_times, max_dev, scenariis)

    ###### Model #####
    m = gp.Model('Extended_Robust_Flowshop')
    m.Params.IntegralityFocus = 1
    m.Params.OutputFlag = 0

    ###### Variables #####

    # Job position variable
    x = m.addVars(nb_jobs, nb_jobs, vtype = GRB.BINARY, name = 'x')

    # Auxilary variable to reorder processing times along the sequence
    ord_processing_times = m.addVars(nb_machines, nb_jobs, nb_scenariis, name = 'ord_processing_times')
    for s in range(nb_scenariis):
        m.addConstrs(
            ord_processing_times[i, k, s] == gp.quicksum(processing_times[s][i][j]*x[j,k] for j in range(nb_jobs)) 
            for k in range(nb_jobs) 
            for i in range(nb_machines))

    # Idle time on machine i BEFORE the processing of job j
    idle = m.addVars(nb_machines, nb_jobs, nb_scenariis, lb =0.0, name = 'idle')

    # Makespan
    makespan = m.addVar(lb =0.0, vtype = GRB.CONTINUOUS, name = 'Cmax')

    ##### Constraints #####

    # Each position in the sequence is occupied by only one job
    m.addConstrs(gp.quicksum(x[j, k] for j in range(nb_jobs)) == 1 for k in range(nb_jobs))
    # A job appears only once in the sequence
    m.addConstrs(gp.quicksum(x[j, k] for k in range(nb_jobs)) == 1 for j in range(nb_jobs))

    for s in range(nb_scenariis):

        # JAML constraints
        m.addConstrs(
            ord_processing_times[i, k+1, s] + idle[i, k+1, s] == 
            ord_processing_times[i+1, k, s] + idle[i+1, k+1, s] 
            for i in range(nb_machines-1) 
            for k in range(nb_jobs-1))
        
        # Idle time for the first job in the sequence        
        m.addConstrs(
            idle[i, 0, s] ==  gp.quicksum(ord_processing_times[l, 0, s] for l in range(i))
            for i in range(nb_machines))

        # Makespan definition
        m.addConstr(makespan >= gp.quicksum(
            ord_processing_times[nb_machines-1, k, s] + idle[nb_machines-1, k, s] 
            for k in range(nb_jobs)))


    ##### Objective function #####
    m.setObjective(makespan, GRB.MINIMIZE)

    ##### Optimization #####
    m.optimize()

    # Retrieve the optimal sequence

    opt_sequence = retrieve_opt_seq_MIP(m)
    print(f'one optimal sequence is {opt_sequence}')

    # Retrieve the worst case makespan
    worst_makespan = m.getObjective().getValue()
    print(f'For a budget of {budget} the worst case makespan is {makespan}')

    return worst_makespan, opt_sequence

def CP_model_extended(
        nom_processing_times: np.ndarray, 
        max_dev: np.ndarray, 
        budget: int
    ) -> tuple[int, np.ndarray]:
    """
    This function optimizes a robust flowshop model following an extended
    CP formulation.
    Meaning all possible scenarios under the budget are generated    

    Parameters:

    nom_processing_times: the nominal processing times time matrix, 
    
    max_deviations: the maximum deviation values matrix
    
    budget: the uncertainty budget
    
    returns: one optimal sequence\n
             the worst makespan under the optimal sequence
    
    """

    ##### Data #####

    nb_machines = len(nom_processing_times)
    nb_jobs = len(nom_processing_times[0])
    scenariis = generate_all_scenariis(budget, nb_machines, nb_jobs)
    nb_scenariis = len(scenariis)

    ##### Initialization #####

    mdl = CpoModel()

    ##### Variables ######

    tasks = integer_var_list(nb_scenariis)
    sequences = []
    for s in range(nb_scenariis):
        processing_times = nom_processing_times + np.multiply(scenariis[s],max_dev) 
        tasks[s] = [
            [interval_var(name='task_{}_{}_{}'.format(i,j,s), size=int(processing_times[i,j])) 
             for j in range(nb_jobs)] 
             for i in range(nb_machines)]
        
        sequence = [sequence_var(tasks[s][i], name='sequence_{}_{}'.format(i,s)) for i in range(nb_machines)]
        sequences.append(sequence)

    Cmax = integer_var(name='Cmax')

    ##### Constraints #####

    for s in range(nb_scenariis):

        mdl.add(greater_or_equal(Cmax, end_of(tasks[s][nb_machines-1][j])) for j in range(nb_jobs))

        mdl.add(end_at_start(tasks[s][i][j], tasks[s][i+1][j]) for i in range(nb_machines-1) for j in range(nb_jobs))

        mdl.add(no_overlap(sequences[s][i]) for  i in range(nb_machines))

        mdl.add(same_sequence(sequences[s][0], sequences[s][i]) for i in range(nb_machines))

        mdl.add(same_sequence(sequences[0][i], sequences[s][i]) for i in range(nb_machines))

    ##### Objective #####

    mdl.minimize(Cmax)

    ##### Solving #####
    msol = mdl.solve()

    ##### Retrieving the data #####

    worst_makespan = msol.get_value(Cmax)

    seqc_var = msol.get_value(sequences[0][0])
    def retrieve_sequence_position(var_name):
        return int(var_name.split('_')[2])
    opt_sequence = []
    for s in seqc_var:
        opt_sequence.append(retrieve_sequence_position(s.get_name()))
    print(opt_sequence)

    return worst_makespan, opt_sequence

def MIP_master_problem(scenario_pool, nom_processing_times: np.ndarray, 
        max_dev: np.ndarray, 
) -> tuple[gp.Model, np.ndarray, int]:
    """
    This function implements the master of the no wait flowshop
    following a MIP formulation 

    Parameters:
    nom_processing_times: the nominal processing times time matrix, 
    
    max_dev: the maximum deviation values matrix
    
    
    returns: one optimal sequence \n  
             the worst makespan under the optimal sequence
    
    """

    ##### Data #####

    nb_machines = len(nom_processing_times)
    nb_jobs = len(nom_processing_times[0])
    nb_scenariis = len(scenario_pool)
 
    processing_times = generate_all_processing_times(nom_processing_times, max_dev, scenario_pool)

    ###### Model #####
    m = gp.Model()
    m.Params.IntegralityFocus = 1
    m.Params.OutputFlag = 0

    ###### Variables #####

    # Job position variable
    x = m.addVars(nb_jobs, nb_jobs, vtype = GRB.BINARY, name = 'x')

    # Auxilary variable to reorder processing times along the sequence
    ord_processing_times = m.addVars(nb_machines, nb_jobs, nb_scenariis, name = 'ord_processing_times')
    for s in range(nb_scenariis):
        m.addConstrs(
            ord_processing_times[i, k, s] == gp.quicksum(processing_times[s][i][j]*x[j,k] for j in range(nb_jobs)) 
            for k in range(nb_jobs) 
            for i in range(nb_machines) )

    # Idle time on machine i BEFORE the processing of job j
    idle = m.addVars(nb_machines, nb_jobs, nb_scenariis, lb =0.0, name = 'idle')

    # Makespan
    makespan = m.addVar(lb =0.0, vtype = GRB.CONTINUOUS, name = 'Cmax')

    ##### Constraints #####

    # Each position in the sequence is occupied by only one job
    m.addConstrs(gp.quicksum(x[j, k] for j in range(nb_jobs)) == 1 for k in range(nb_jobs))
    # A job appears only once in the sequence
    m.addConstrs(gp.quicksum(x[j, k] for k in range(nb_jobs)) == 1 for j in range(nb_jobs))

    for s in range(nb_scenariis):

        # JAML constraints
        m.addConstrs(
            ord_processing_times[i, k+1, s] + idle[i, k+1, s] == 
            ord_processing_times[i+1, k, s] + idle[i+1, k+1, s] 
            for i in range(nb_machines-1) 
            for k in range(nb_jobs-1))
        
        # Idle time for the first job in the sequence        
        m.addConstrs(
            idle[i, 0, s] ==  gp.quicksum(ord_processing_times[l, 0, s] for l in range(i)) 
            for i in range(nb_machines))

        # Makespan definition
        m.addConstr(
            makespan >= gp.quicksum(
                ord_processing_times[nb_machines-1, k, s] + idle[nb_machines-1, k, s] 
                for k in range(nb_jobs)))


    ##### Objective function #####
    m.setObjective(makespan, GRB.MINIMIZE)

    ##### Optimization #####
    m.optimize()

    # Retrieve the optimal sequence 
    opt_sequence = retrieve_opt_seq_MIP(m)
    print(f'one optimal sequence is {opt_sequence}')

    # Retrieve the worst case makespan
    worst_makespan = m.getObjective().getValue()
    #print(f'For a budget of {budget} the worst case makespan is {makespan}')

    return worst_makespan, opt_sequence

def MIP_adverse_problem(
        job_sequence: np.ndarray, 
        nom_processing_times: np.ndarray,
        max_dev: np.ndarray, 
        budget: int,
        plotgantt: int
) -> tuple[int, np.ndarray]:
    """
    Implements the adverse problem in a MIP formulation
    """

    ##### Data #####
    nom_processing_times = reorder_sequence(job_sequence, nom_processing_times)
    max_dev = reorder_sequence(job_sequence, max_dev)

    nb_machines = len(nom_processing_times)
    nb_jobs = len(nom_processing_times[0])
    M = 1000             # An arbitrary large number

    ##### Initialize the model #####
    model = gp.Model('AP_robust_flowshop')
    model.Params.IntegralityFocus = 1
    model.Params.OutputFlag = 0
    model.Params.InfUnbdInfo = 1

    ##### Variables #####

    # WARNING: the scenario will be calculated for the given sequence
    # khi[i,k] = 1 if the procesing time for the kth job in the sequence
    # on the ith machine deviates from its nominal value, 0 otherwise
    khi = model.addVars(nb_machines, nb_jobs, vtype = GRB.BINARY, name = 'khi')

    # Idle times before the treatment of jobs
    idle = model.addVars(nb_machines, nb_jobs-1, lb = 0.0, ub = GRB.INFINITY, name = 'idle')

    # Auxilary variable for definig the idle times
    z = model.addVars(nb_machines, nb_jobs-1, vtype=GRB.BINARY, name='z')

    ##### Constraints #####

    # Ensures that for one job, at least one idle time is null 
    for j in range(nb_jobs-1):
        model.addConstrs((idle[i,j] <= M *(1- z[i,j]) for i in range(nb_machines)))
        model.addConstr(gp.quicksum(z[i,j] for i in range(nb_machines)) >= 1)


    # Budget constraint
    model.addConstr(gp.quicksum(khi[i,j] for i in range(nb_machines) for j in range(nb_jobs)) <= budget)

    # JAML constraints
    model.addConstrs(

        nom_processing_times[i, k+1] 
        + khi[i, k+1]*max_dev[i, k+1] 
        + idle[i,k] 
        == 
        nom_processing_times[i+1, k] 
        + khi[i+1, k]*max_dev[i+1, k]  
        + idle[i+1,k] 

        for i in range(nb_machines-1) 
        for k in range(nb_jobs-1))
    
    

    ##### Objective #####

    # Maximize the makespan
    model.setObjective(
        gp.quicksum(nom_processing_times[i, 0] + khi[i, 0]*max_dev[i, 0] for i in range(nb_machines)) 
        + gp.quicksum(idle[nb_machines-1, j] for j in range(nb_jobs-1))
        + gp.quicksum(nom_processing_times[nb_machines-1, j] + khi[nb_machines-1, j]*max_dev[nb_machines-1, j] for j in range(1, nb_jobs)),
        GRB.MAXIMIZE)

    ##### Optimization time ! #####
    model.optimize()

    ##### Retrieve data #####

    # Makespan
    worst_makespan = model.getObjective().getValue()

    # Retrieving the idle times in a matrix
    idle_mat = np.zeros(np.shape(nom_processing_times))
    for i in range(1, nb_machines):
        idle_mat[i,0] = sum(nom_processing_times[k,0] + khi[k, 0].X*max_dev[k, 0] for k in range(i))
    for i in range(nb_machines):
        for j in range(1, nb_jobs):
            idle_mat[i, j] = idle[i, j-1].X
    #print(f'idle times\n', idle_mat)

    # Retrieving the processing times
    pt = np.zeros((nb_machines, nb_jobs))
    for i in range(nb_machines):
        for j in range(nb_jobs):
            pt[i, j] = nom_processing_times[i, j] + khi[i,j].X*max_dev[i,j]
    
    # Retrieving the critical scenario in the correct sequence
    scen = np.zeros(np.shape(nom_processing_times))
    scen_sequence = np.zeros(np.shape(job_sequence))
    for j in range(nb_jobs):
        scen_sequence[job_sequence[j]] = j
    scen_sequence = scen_sequence.astype(int)
    
    for i in range(nb_machines):
        for j in range(nb_jobs):
            scen[i,j] = khi[i,j].X
    print(f'One critcal scenario is:\n{reorder_sequence(scen_sequence, scen)}')
    scen = reorder_sequence(scen_sequence, scen)


    if plotgantt:
        flowshop_gantt.ganttRobust(nom_processing_times, max_dev, scen, idle_mat, worst_makespan)

    return worst_makespan, scen

def CP_master_problem(
        scenario_pool: list[np.ndarray], 
        nom_processing_times: np.ndarray, 
        max_dev: np.ndarray
    ) -> tuple[int, np.ndarray]:

    """Implements the master problem in a CP model"""
    ##### Data #####

    nb_machines = len(nom_processing_times)
    nb_jobs = len(nom_processing_times[0])
    nb_scenariis = len(scenario_pool)

    ##### Initialization #####

    mdl = CpoModel()

    ##### Variables ######

    tasks = integer_var_list(nb_scenariis)
    sequences = []
    for s in range(nb_scenariis):
        processing_times = nom_processing_times + np.multiply(scenario_pool[s],max_dev) 
        tasks[s] = [
            [interval_var(name='task_{}_{}_{}'.format(i,j,s), size=int(processing_times[i,j])) 
             for j in range(nb_jobs)] 
             for i in range(nb_machines)]
        sequence = [sequence_var(tasks[s][i], name='sequence_{}_{}'.format(i,s)) for i in range(nb_machines)]
        sequences.append(sequence)

    Cmax = integer_var(name='Cmax')

    ##### Constraints #####

    for s in range(nb_scenariis):

        mdl.add(greater_or_equal(Cmax, end_of(tasks[s][nb_machines-1][j])) for j in range(nb_jobs))

        mdl.add(end_at_start(tasks[s][i][j], tasks[s][i+1][j]) for i in range(nb_machines-1) for j in range(nb_jobs))

        mdl.add(no_overlap(sequences[s][i]) for i in range(nb_machines))

        mdl.add(same_sequence(sequences[s][0], sequences[s][i]) for i in range(nb_machines))

        mdl.add(same_sequence(sequences[0][i], sequences[s][i]) for i in range(nb_machines))

    ##### Objective #####

    mdl.minimize(Cmax)

    ##### Solving #####
    msol = mdl.solve()

    ##### Retrieving the data #####

    worst_makespan = msol.get_value(Cmax)

    seqc_var = msol.get_value(sequences[0][0])
    def retrieve_sequence_position(var_name):
        return int(var_name.split('_')[2])
    opt_sequence = []
    for s in seqc_var:
        opt_sequence.append(retrieve_sequence_position(s.get_name()))
    print(opt_sequence)

    return worst_makespan, opt_sequence

def CP_adverse_problem(
        job_sequence: np.ndarray, 
        nom_processing_times: np.ndarray, 
        max_dev: np.ndarray, 
        budget: int, 
        plot_gantt: int
    ) -> tuple[int, np.ndarray]:
    """This function implements the adverse problem via a CP model"""
    ##### Data #####

    nb_machines = len(nom_processing_times)
    nb_jobs = len(nom_processing_times[0])
    nom_processing_times = reorder_sequence(job_sequence, nom_processing_times)
    max_dev = reorder_sequence(job_sequence, max_dev)

    ##### Initialization #####

    mdl = CpoModel()

    ##### Variables #####

    cmax = integer_var(name = "Cmax")

    tasks = [[interval_var(name='tasks_{}_{}'.format(i,j), size=int(nom_processing_times[i][j])) for j in range(nb_jobs)] for i in range(nb_machines)]

    dev_tasks = [[interval_var(name='dev_{}_{}'.format(i,j), size=int(max_dev[i][j]), optional=True) for j in range(nb_jobs)] for i in range(nb_machines)]


    slack = mdl.integer_var_list(nb_jobs,name = "slack") 

    ##### Constraints #####

    mdl.add(start_of(tasks[0][0]) == 0)

    mdl.add(sum([presence_of(dev_tasks[i][j]) for j in range(nb_jobs) for i in range(nb_machines) ]) <= budget)

    mdl.add(start_at_end(dev_tasks[i][j], tasks[i][j]) for j in range(nb_jobs) for i in range(nb_machines))

    mdl.add(start_of(tasks[i+1][j]) == max(end_of(tasks[i][j]), end_of(dev_tasks[i][j])) for i in range(nb_machines-1) for j in range(nb_jobs))

    mdl.add(cmax == max([end_of(tasks[nb_machines-1][j]) for j in range(nb_jobs)] + [end_of(dev_tasks[nb_machines-1][j]) for j in range(nb_jobs)]))

    mdl.add(slack[j] == min([start_of(tasks[i][j]) - max(end_of(dev_tasks[i][j-1]), end_of(tasks[i][j-1])) for i in range(nb_machines)]) for j in range(1, nb_jobs))

    mdl.add(mdl.max([slack[j] for j in range(nb_jobs)])==0)

    ##### Objective ######

    mdl.maximize(cmax)

    ##### Solving #####

    msol = mdl.solve(LogVerbosity = 'Quiet')

    slvd_tasks = [[msol.get_var_solution(tasks[i][j]) for j in range(nb_jobs)] for i in range(nb_machines)]
    slvd_devs = [[msol.get_var_solution(dev_tasks[i][j]) for j in range(nb_jobs)] for i in range(nb_machines)]

    makespan = msol.get_value(cmax)

    # Retrieving the critical scenario in the correct sequence
    scen = np.zeros(np.shape(nom_processing_times))
    scen_sequence = np.zeros(np.shape(job_sequence))
    for j in range(nb_jobs):
        scen_sequence[job_sequence[j]] = j
    scen_sequence = scen_sequence.astype(int)

    for i in range(nb_machines):
        for j in range(nb_jobs):
            if slvd_devs[i][j].get_start():
                    scen[i][j] = 1
    scen = reorder_sequence(scen_sequence, scen)

    if plot_gantt:
        flowshop_gantt.ganttCP_robust(slvd_tasks, slvd_devs, makespan)

    return makespan, scen


def MIP_MIP_optimization(scenario_pool, nom_processing_times, max_dev, budget):
    
    ##### Data #####
    nb_machines = len(nom_processing_times)
    nb_jobs = len(nom_processing_times[0])
    nb_scenariis = np.sum([2 ** math.comb(nb_machines * nb_jobs, gamma) for gamma in range(1, budget+1)])

    ##### Run first optimizations #####
    master_cmax, master_sequence = MIP_master_problem(scenario_pool, nom_processing_times, max_dev)
    adverse_cmax, adverse_scenario = MIP_adverse_problem(master_sequence, nom_processing_times, max_dev, budget, 0)

    print(f'initial LB:{master_cmax} intial UB: {adverse_cmax}')
    while adverse_cmax > master_cmax:
        scenario_pool.append(adverse_scenario)
        print(f'Percentage of scenarios explored {(len(scenario_pool)/nb_scenariis)*100}')
        master_cmax, master_sequence = MIP_master_problem(scenario_pool, nom_processing_times, max_dev)
        print(f'LB:{master_cmax} UB:{adverse_cmax}')
        adverse_cmax, adverse_scenario = MIP_adverse_problem(master_sequence, nom_processing_times, max_dev, budget, 0)
    
    
    print(f"Optimal makespan: {master_cmax}")
    print(f"Optimal sequence: {master_sequence}")
    print(f'Final percentage of scenarios explored {(len(scenario_pool)/nb_scenariis)*100}')
    print(f"Number of unexplored scenarios: {nb_scenariis - len(scenario_pool)}")
    
    return master_cmax, master_sequence

def CP_CP_optimization(scenario_pool, nom_processing_times, max_dev, budget):
    
    ##### Data #####
    nb_machines = len(nom_processing_times)
    nb_jobs = len(nom_processing_times[0])
    nb_scenariis = np.sum([2 ** math.comb(nb_machines * nb_jobs, gamma) for gamma in range(1, budget+1)])

    ##### Run first optimizations #####
    master_cmax, master_sequence = CP_master_problem(scenario_pool, nom_processing_times, max_dev)
    adverse_cmax, adverse_scenario = CP_adverse_problem(master_sequence, nom_processing_times, max_dev, budget, 0)

    print(f'initial LB:{master_cmax} intial UB: {adverse_cmax}')
    print(f'scenario pool \n {scenario_pool}')
    while adverse_cmax > master_cmax:
        scenario_pool.append(adverse_scenario)
        print(f'Percentage of scenarios explored {(len(scenario_pool)/nb_scenariis)*100}')
        master_cmax, master_sequence = CP_master_problem(scenario_pool, nom_processing_times, max_dev)
        print(f'LB:{master_cmax} UB:{adverse_cmax}')
        adverse_cmax, adverse_scenario = CP_adverse_problem(master_sequence, nom_processing_times, max_dev, budget, 0)
        
    
    print(f"Optimal makespan: {master_cmax}")
    print(f"Optimal sequence: {master_sequence}")
    print(f"Total scenarios explored: {len(scenario_pool)} out of {nb_scenariis} possible scenarios")
    print(f"Number of unexplored scenarios: {nb_scenariis - len(scenario_pool)}")
    
    return master_cmax, master_sequence

def MIP_CP_optimization(scenario_pool, nom_processing_times, max_dev, budget):
    
    ##### Data #####
    nb_machines = len(nom_processing_times)
    nb_jobs = len(nom_processing_times[0])
    nb_scenariis = np.sum([2 ** math.comb(nb_machines * nb_jobs, gamma) for gamma in range(1, budget+1)])

    ##### Run first optimizations #####
    master_cmax, master_sequence = MIP_master_problem(scenario_pool, nom_processing_times, max_dev)
    adverse_cmax, adverse_scenario = CP_adverse_problem(master_sequence, nom_processing_times, max_dev, budget, 0)

    print(f'initial LB:{master_cmax} intial UB: {adverse_cmax}')
    while adverse_cmax > master_cmax:
        scenario_pool.append(adverse_scenario)
        print(f'Percentage of scenarios explored {(len(scenario_pool)/nb_scenariis)*100}')
        master_cmax, master_sequence = MIP_master_problem(scenario_pool, nom_processing_times, max_dev)
        print(f'LB:{master_cmax} UB:{adverse_cmax}')
        adverse_cmax, adverse_scenario = CP_adverse_problem(master_sequence, nom_processing_times, max_dev, budget, 0)
        
    
    print(f"Optimal makespan: {master_cmax}")
    print(f"Optimal sequence: {master_sequence}")
    print(f'Final percentage of scenarios explored {(len(scenario_pool)/nb_scenariis)*100}')
    print(f"Number of unexplored scenarios: {nb_scenariis - len(scenario_pool)}")
    
    return master_cmax, master_sequence

def CP_MIP_optimization(scenario_pool, nom_processing_times, max_dev, budget):
    
    ##### Data #####
    nb_machines = len(nom_processing_times)
    nb_jobs = len(nom_processing_times[0])
    nb_scenariis = np.sum([2 ** math.comb(nb_machines * nb_jobs, gamma) for gamma in range(1, budget+1)])

    ##### Run first optimizations #####
    master_cmax, master_sequence = CP_master_problem(scenario_pool, nom_processing_times, max_dev)
    adverse_cmax = 10000
    print(f'initial LB:{master_cmax} intial UB: {adverse_cmax}')
    while adverse_cmax > master_cmax:
        adverse_cmax, adverse_scenario = MIP_adverse_problem(master_sequence, nom_processing_times, max_dev, budget, 0)
        scenario_pool.append(adverse_scenario)
        print(f'Percentage of scenarios explored {(len(scenario_pool)/nb_scenariis)*100}')
        master_cmax, master_sequence = CP_master_problem(scenario_pool, nom_processing_times, max_dev)
        print(f'LB:{master_cmax} UB:{adverse_cmax}')

    
    print(f"Optimal makespan: {master_cmax}")
    print(f"Optimal sequence: {master_sequence}")
    print(f'Final percentage of scenarios explored {(len(scenario_pool)/nb_scenariis)*100}')
    print(f"Number of unexplored scenarios: {nb_scenariis - len(scenario_pool)}")
    
    return master_cmax, master_sequence
