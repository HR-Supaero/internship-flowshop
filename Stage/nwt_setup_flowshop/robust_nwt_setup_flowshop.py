from itertools import combinations, permutations
import gurobipy as gp
import random
import numpy as np
import numpy.linalg as alg
from gurobipy import GRB
import flowshop as flowshop
from docplex.cp.model import *
import nowait_setup_flowshop as nowait_setup_flowshop
import math

##### Helper functions #####

def is_scenario_in_list(scenario, scenario_pool):
    """ True if the scenario is in the pool, False otherwise."""
    return any(np.array_equal(scenario, scen) for scen in scenario_pool)


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

    job_sequence = np.array(job_sequence, dtype=int)
    return mat[:, job_sequence]

def calculate_processing_times(
        scenario: np.ndarray, 
        nom_production_times: np.ndarray, 
        max_deviations: np.ndarray
) -> np.ndarray:
    """
    This function returns the production time matrix in a certain scenario
    """
    return nom_production_times + np.multiply(scenario, max_deviations)

def generate_all_scenariis(
        budget: int, 
        nb_machines: int, 
        nb_jobs: int
) -> np.ndarray:
    """
    This code returns an array of all possible sceanriis fitting a certain budget
    """
    list_scen = []
    
    for gamma in range(budget+1):
        possible_ones = list(combinations(range(nb_machines*nb_jobs), gamma))
        for comb in possible_ones:
            scenario = np.zeros((nb_machines, nb_jobs))
            for index in comb:
                row = index // nb_jobs
                col = index % nb_jobs
                scenario[row, col] = 1
            list_scen.append(scenario.copy())
    return list_scen

def generate_some_scenariis(
        nb_scenariis,
        budget, 
        nb_machines, 
        nb_jobs
):
    """Generates a random sample of scenariis"""
    
    list_scen = []
    if nb_scenariis > 2**(nb_machines*nb_jobs):
        raise ValueError("N is greater than the total number of possible scenarios.")
    
    for gamma in range(budget+1):
        possible_ones = list(combinations(range(nb_machines*nb_jobs), gamma))
        for comb in possible_ones:
            scenario = np.zeros((nb_machines, nb_jobs))
            for index in comb:
                row = index // nb_jobs
                col = index % nb_jobs
                scenario[row, col] = 1
            list_scen.append(scenario.copy())
    sampled_scenariis = random.sample(list_scen, nb_scenariis)
    return sampled_scenariis


def setup_checker(job_sequence, nom_processing_times, max_dev, setup_times, budget, worst_makespan):
    """
    Checks the robustness of a given solution
    """

    # Adjust the sequence for setup times
    opt_sequence = job_sequence.copy()
    for k in range(len(job_sequence)):
        opt_sequence[k] += 1
    opt_sequence = np.insert(opt_sequence, 0, 0)

    ##### Data #####
    nb_machines = len(nom_processing_times)
    nb_jobs = len(nom_processing_times[0])
    
    # Reorder matrices based on job sequence
    reordered_nom_times = reorder_sequence(job_sequence, nom_processing_times)
    reordered_max_dev = reorder_sequence(job_sequence, max_dev)
    
    # Retrieve the setup times
    setup_mat = np.zeros(np.shape(reordered_nom_times))
    for i in range(nb_machines):
        for j in range(nb_jobs):
            setup_mat[i, j] = setup_times[i][int(opt_sequence[j])][int(opt_sequence[j + 1])]

    # Check all scenarios for the given sequence
    scenariis = generate_all_scenariis(budget, nb_machines, nb_jobs)
    for scen in scenariis:
        #print(f'------------Scenario: \n {scen} ------------')
        reordered_scen = reorder_sequence(job_sequence, scen)
        processing_times = calculate_processing_times(reordered_scen, reordered_nom_times, reordered_max_dev)
        makespan = nowait_setup_flowshop.cmax(processing_times, setup_mat)
        #print(f'Makespan of {makespan}')
        if makespan > worst_makespan:
            print(f'The scenario \n{scen} violates the worst-case makespan with a makespan of {makespan}')
            return False, scen
        
     # Check other sequences for at least one makespan longer than the worst makespan
    other_sequences = list(permutations(range(nb_jobs)))
    other_sequences.remove(tuple(job_sequence))
    
    for seq in other_sequences:
        #print(f'------------New sequence :{seq}------------')
        offset_seq = list(seq).copy()
        for k in range(len(seq)):
            offset_seq[k] += 1
        offset_seq = np.insert(offset_seq, 0, 0)
        
        reordered_nom_times = reorder_sequence(list(seq), nom_processing_times)
        reordered_max_dev = reorder_sequence(list(seq), max_dev)
        
        longer_makespan_found = False
        setup_mat = np.zeros(np.shape(reordered_nom_times))
        for i in range(nb_machines):
            for j in range(nb_jobs):
                setup_mat[i, j] = setup_times[i][offset_seq[j]][offset_seq[j + 1]]

        for scen in scenariis:
            reordered_scen = reorder_sequence(seq, scen)
            processing_times = calculate_processing_times(reordered_scen, reordered_nom_times, reordered_max_dev)
            makespan = nowait_setup_flowshop.cmax(processing_times, setup_mat)
            if makespan >= worst_makespan:
                longer_makespan_found = True
                #print(f'longer makespan of {makespan} found for scenario:\n {scen}')
                break

        if not longer_makespan_found:
            print(f'All makespans for sequence {seq} are shorter to the worst makespan.')
            return False, seq
    
    return True

##### Models #####

def MIP_model_extended(nom_processing_times, max_dev, setup_times, budget):

    ##### Data #####

    nb_machines = len(nom_processing_times)
    nb_jobs = len(nom_processing_times[0])
    A = 100000
    
    scenarios = generate_all_scenariis(budget, nb_machines, nb_jobs)
    nb_scenariis = len(scenarios)

    processing_times = []       # Array containing every possible processing times
    for scen in scenarios:
        pt_s = calculate_processing_times(scen, nom_processing_times, max_dev)
        processing_times.append(pt_s.copy())

    # Edge cases for the setup times
    for i in range(nb_machines):
        for j in range(nb_jobs+1):
            setup_times[i][j][j] = 0    
            setup_times[i][j][0] = 0 

    ##### Initializing the model #####

    model = gp.Model('robust_Flowshop_wsetup')
    model.Params.IntegralityFocus = 1
    model.Params.OutputFlag = 1

    ##### Variables #####

    # idle[i, j, k] idle time on machine i 
    # between the end of job j and the start of job k 
    # assuming k follows j in the sequence
    # j = 0 corresponds to the dummy job
    idle = model.addVars(nb_machines, nb_jobs + 1, nb_jobs + 1, nb_scenariis, lb = 0.0, ub = GRB.INFINITY, name = "idle")

    # difference[j, k] time difference between 
    # the completion times of job j and job k on the last machine
    # assuming k follows j in the sequence
    # j = 0 corresponds to the dummy job
    difference = model.addVars(nb_jobs + 1, nb_jobs + 1, nb_scenariis, lb = 0.0, ub = GRB.INFINITY, name = 'difference')

    # u[j] corresponds to the position of job j in the sequence
    # j = 0 corresponds to the dummy job
    u = model.addVars(nb_jobs + 1, vtype = GRB.INTEGER, name = 'u')

    #z[j,k] equals 1 if job k follows job j in the sequence, 0 otherwise
    z = model.addVars(nb_jobs + 1, nb_jobs + 1, vtype = GRB.BINARY, name = 'z')

    Cmax = model.addVar(name='Cmax')

    ##### Constraints #####

    # Edge cases
    # The idle time, completion time difference, following decision variables
    # cant be defined if the following job is the dummy job
    # or if the current and following job are the same
    for s in range(nb_scenariis):
        model.addConstrs(idle[i, j, 0, s] == 0 
                    for i in range(nb_machines) 
                    for j in range(nb_jobs+1))
        
        model.addConstrs(idle[i, j, j, s] == 0 
                    for i in range(nb_machines) 
                    for j in range(nb_jobs+1))

        model.addConstrs(difference[j,j, s] == 0 
                    for j in range(nb_jobs+1))

        model.addConstrs(difference[j, 0, s] == 0 
                    for j in range(nb_jobs+1))

    model.addConstrs(z[j,j] == 0 
                 for j in range(nb_jobs+1))

    model.addConstrs(z[j,0] == 0 
                 for j in range(nb_jobs+1))

    for s in range(nb_scenariis):
        # The idle times between jobs has to be bigger than the setup times
        model.addConstrs(idle[i,j,k, s] >= setup_times[i][j][k] 
                    for i in range(nb_machines) 
                    for j in range(nb_jobs+1) 
                    for k in range(nb_jobs+1))

        # JAML constraints
        model.addConstrs(processing_times[s][i][j-1] 
                    + idle[i,j,k, s] 
                    == 
                    processing_times[s][i-1][k-1] 
                    + idle[i-1,j,k,s] 
                    for i in range(1, nb_machines) 
                    for j in range(1, nb_jobs + 1) 
                    for k in range(1, nb_jobs+1) 
                    if j!=k)

        # JAML constraints for the first job in the sequence
        model.addConstrs(idle[i,0,k,s] == 
                    idle[i-1,0,k,s] 
                    + processing_times[s][i-1][k-1] 
                    for i in range(1, nb_machines) 
                    for k in range(1, nb_jobs+1))

        # completion time difference > 0 
        # if and only if j is followed by k, otherwise its null
        model.addConstrs(difference[j,k,s] >= idle[nb_machines-1,j,k, s] 
                                        + processing_times[s][nb_machines-1][k-1] 
                                        - A*(1-z[j,k]) 
                                        for j in range(nb_jobs + 1) 
                                        for k in range(nb_jobs + 1))

        model.addConstr(gp.quicksum(difference[j,k,s] for j in range(nb_jobs+1) for k in range(nb_jobs+1)) <= Cmax)
    # Each job can only have one job after it
    model.addConstrs(gp.quicksum(z[j,k] for k in range(nb_jobs + 1)) <= 1 
                 for j in range(nb_jobs+1))
    # Each job has one job preceding it apart from the dummy job
    model.addConstrs(gp.quicksum(z[j,k] for j in range(nb_jobs + 1)) >= 1 
                 for k in range(1,nb_jobs+1))

    # Defines the job position variable
    model.addConstrs(u[j] - u[k] + nb_jobs*z[j,k] + (nb_jobs-2)*z[k,j] <= nb_jobs-1 
                 for j in range(1, nb_jobs+1) 
                 for k in range(1, nb_jobs+1))
    model.addConstr(u[0] == 0)
    model.addConstrs(u[j] >= 1 for j in range(1,nb_jobs+1))
    model.addConstrs(u[j] <= nb_jobs for j in range(1,nb_jobs+1))


    ##### Objective #####

    model.setObjective(Cmax, GRB.MINIMIZE )

    ##### Optimization #####
    model.optimize()

    ##### Retrieve data #####
    
    # Retrieve the optimal sequence and the makespan
    opt_sequence = []
    for v in model.getVars():
        if v.VarName[0] == 'u':
            opt_sequence.append((v.VarName, int(v.X)))
    opt_sequence = sorted(opt_sequence, key=lambda elem: elem[1])
    def retrieve_sequencePos(varName):
        return int(varName[0].split('[')[1].split(']')[0])
    for k in range(len(opt_sequence)):
        opt_sequence[k] = retrieve_sequencePos(opt_sequence[k])
    
    job_sequence = opt_sequence.copy()         # Will contain the optimal sequence w/o the dummy job
    for k in range(len(job_sequence)):
        job_sequence[k] = job_sequence[k] - 1  # job indexes in opt_sequence are offseted by 1 due to the dummy job in index 0
    job_sequence.remove(-1)
    #print(f'One optimal sequence is {job_sequence}')
    
    makespan = model.getObjective().getValue()
    if not makespan.is_integer():
        raise ValueError("The makespan is not of integer value")
    #print(f'The makespan is {makespan}')
    return makespan, job_sequence

def CP_model_extended(nom_processing_times, max_dev, setup_times, budget):
    ##### Data #####
    nb_machines = len(nom_processing_times)
    nb_jobs = len(nom_processing_times[0])

    scenarios = generate_all_scenariis(budget, nb_machines, nb_jobs)
    nb_scenariis = len(scenarios)

    for i in range(nb_machines):
        for j in range(nb_jobs+1):
             setup_times[i][j][j] = 0
             setup_times[i][j][0] = 0
    print(setup_times)

    ##### Initialize #####
    mdl = CpoModel()

    ##### Variables #####

    tasks = integer_var_list(nb_scenariis)
    sequences = []
    for s in range(nb_scenariis):
        processing_times = nom_processing_times + np.multiply(scenarios[s],max_dev) 
        processing_times = np.hstack((np.zeros((nb_machines,1)), processing_times))
        tasks[s] = [[
            interval_var(
                name='task_{}_{}_{}'.format(i,j,s), 
                size=int(processing_times[i,j])) 
                    for j in range(nb_jobs+1)] 
                    for i in range(nb_machines)]
        sequence = [sequence_var(tasks[s][i], name='sequence_{}_{}'.format(i,s)) for i in range(nb_machines)]
        sequences.append(sequence)

    Cmax = integer_var(name='Cmax')

    setup_mat = [transition_matrix(setup_times[i]) for i in range(nb_machines)]

    ##### Constraints #####

    for s in range(nb_scenariis):
        
        mdl.add(start_of(tasks[s][i][0]) == 0 for i in range(nb_machines))
        
        mdl.add(greater_or_equal(Cmax, end_of(tasks[s][nb_machines-1][j])) for j in range(nb_jobs+1))

        mdl.add(end_at_start(tasks[s][i][j], tasks[s][i+1][j]) for i in range(nb_machines-1) for j in range(nb_jobs+1))

        mdl.add(no_overlap(sequences[s][i], setup_mat[i]) for i in range(nb_machines))

        mdl.add(same_sequence(sequences[s][0], sequences[s][i]) for i in range(1, nb_machines))

        mdl.add(same_sequence(sequences[0][i], sequences[s][i]) for i in range(nb_machines))
        
    ##### Objective #####

    mdl.minimize(Cmax)

    ##### Solve ! #####

    msol = mdl.solve()


    ##### Retrieve the solutions #####
    makespan = msol.get_value(Cmax)

    seqc_var = msol.get_value(sequences[0][0])

    print(f'The setup are {setup_times}')
    def retrieve_sequence_position(var_name):
        return int(var_name.split('_')[2])
    opt_sequence = []
    for s in seqc_var:
        opt_sequence.append(retrieve_sequence_position(s.get_name()))
    print(opt_sequence)

    job_sequence = opt_sequence.copy()         # Will contain the optimal sequence w/o the dummy job
    for k in range(len(job_sequence)):
        job_sequence[k] = job_sequence[k] - 1  # job indexes in opt_sequence are offseted by 1 due to the dummy job in index 0
    job_sequence.remove(-1)
    print(f'one optimal sequence is: {job_sequence}')

    opt_sequence.remove(0)
    opt_sequence.insert(0, 0)
    print(f'opt sequence w: dummy: {opt_sequence}')
     
    slvd_setup_mat = np.zeros((nb_machines, nb_jobs))
    for i in range(nb_machines):
        for j in range(nb_jobs):
            slvd_setup_mat[i][j] = setup_times[i][opt_sequence[j]][opt_sequence[j+1]]
    print(slvd_setup_mat)

    return makespan, job_sequence

def MIP_master_problem(scenario_pool, nom_processing_times, max_dev, setup_times):
##### Data #####

    nb_machines = len(nom_processing_times)
    nb_jobs = len(nom_processing_times[0])
    A = 100000
    nb_scenariis = len(scenario_pool)

    processing_times = []       # Array containing every possible processing times
    for scen in scenario_pool:
        pt_s = calculate_processing_times(scen, nom_processing_times, max_dev)
        processing_times.append(pt_s.copy())

    # Edge cases for the setup times
    for i in range(nb_machines):
        for j in range(nb_jobs+1):
            setup_times[i][j][j] = 0    
            setup_times[i][j][0] = 0 

    ##### Initializing the model #####

    model = gp.Model('robust_Flowshop_wsetup')
    model.Params.IntegralityFocus = 1
    model.Params.OutputFlag = 0

    ##### Variables #####

    # idle[i, j, k] idle time on machine i 
    # between the end of job j and the start of job k 
    # assuming k follows j in the sequence
    # j = 0 corresponds to the dummy job
    idle = model.addVars(nb_machines, nb_jobs + 1, nb_jobs + 1, nb_scenariis, lb = 0.0, ub = GRB.INFINITY, name = "idle")

    # difference[j, k] time difference between 
    # the completion times of job j and job k on the last machine
    # assuming k follows j in the sequence
    # j = 0 corresponds to the dummy job
    difference = model.addVars(nb_jobs + 1, nb_jobs + 1, nb_scenariis, lb = 0.0, ub = GRB.INFINITY, name = 'difference')

    # u[j] corresponds to the position of job j in the sequence
    # j = 0 corresponds to the dummy job
    u = model.addVars(nb_jobs + 1, vtype = GRB.INTEGER, name = 'u')

    #z[j,k] equals 1 if job k follows job j in the sequence, 0 otherwise
    z = model.addVars(nb_jobs + 1, nb_jobs + 1, vtype = GRB.BINARY, name = 'z')

    Cmax = model.addVar(name='Cmax')


    ##### Constraints #####

    # Edge cases
    # The idle time, completion time difference, following decision variables
    # cant be defined if the following job is the dummy job
    # or if the current and following job are the same
    for s in range(nb_scenariis):
        model.addConstrs(idle[i, j, 0, s] == 0 
                    for i in range(nb_machines) 
                    for j in range(nb_jobs+1))
        
        model.addConstrs(idle[i, j, j, s] == 0 
                    for i in range(nb_machines) 
                    for j in range(nb_jobs+1))

        model.addConstrs(difference[j,j, s] == 0 
                    for j in range(nb_jobs+1))

        model.addConstrs(difference[j, 0, s] == 0 
                    for j in range(nb_jobs+1))

    model.addConstrs(z[j,j] == 0 
                 for j in range(nb_jobs+1))

    model.addConstrs(z[j,0] == 0 
                 for j in range(nb_jobs+1))

    for s in range(nb_scenariis):
        # The idle times between jobs has to be bigger than the setup times
        model.addConstrs(idle[i,j,k, s] >= setup_times[i][j][k] 
                    for i in range(nb_machines) 
                    for j in range(nb_jobs+1) 
                    for k in range(nb_jobs+1))

        # JAML constraints
        model.addConstrs(processing_times[s][i][j-1] 
                    + idle[i,j,k, s] 
                    == 
                    processing_times[s][i-1][k-1] 
                    + idle[i-1,j,k,s] 
                    for i in range(1, nb_machines) 
                    for j in range(1, nb_jobs + 1) 
                    for k in range(1, nb_jobs+1) 
                    if j!=k)

        # JAML constraints for the first job in the sequence
        model.addConstrs(idle[i,0,k,s] == 
                    idle[i-1,0,k,s] 
                    + processing_times[s][i-1][k-1] 
                    for i in range(1, nb_machines) 
                    for k in range(1, nb_jobs+1))

        # completion time difference > 0 
        # if and only if j is followed by k, otherwise its null
        model.addConstrs(difference[j,k,s] >= idle[nb_machines-1,j,k, s] 
                                        + processing_times[s][nb_machines-1][k-1] 
                                        - A*(1-z[j,k]) 
                                        for j in range(nb_jobs + 1) 
                                        for k in range(nb_jobs + 1))

        model.addConstr(gp.quicksum(difference[j,k,s] for j in range(nb_jobs+1) for k in range(nb_jobs+1)) <= Cmax)
    # Each job can only have one job after it
    model.addConstrs(gp.quicksum(z[j,k] for k in range(nb_jobs + 1)) <= 1 
                 for j in range(nb_jobs+1))
    # Each job has one job preceding it apart from the dummy job
    model.addConstrs(gp.quicksum(z[j,k] for j in range(nb_jobs + 1)) >= 1 
                 for k in range(1,nb_jobs+1))

    # Defines the job position variable
    model.addConstrs(u[j] - u[k] + nb_jobs*z[j,k] + (nb_jobs-2)*z[k,j] <= nb_jobs-1 
                 for j in range(1, nb_jobs+1) 
                 for k in range(1, nb_jobs+1))
    model.addConstr(u[0] == 0)
    model.addConstrs(u[j] >= 1 for j in range(1,nb_jobs+1))
    model.addConstrs(u[j] <= nb_jobs for j in range(1,nb_jobs+1))


    ##### Objective #####

    model.setObjective(Cmax, GRB.MINIMIZE )

    ##### Optimization #####
    model.optimize()

    ##### Retrieve data #####

    # Retrieve the optimal sequence and the makespan
    opt_sequence = []
    for v in model.getVars():
        if v.VarName[0] == 'u':
            opt_sequence.append((v.VarName, int(v.X)))
    opt_sequence = sorted(opt_sequence, key=lambda elem: elem[1])
    def retrieve_sequencePos(varName):
        return int(varName[0].split('[')[1].split(']')[0])
    for k in range(len(opt_sequence)):
        opt_sequence[k] = retrieve_sequencePos(opt_sequence[k])
    


    job_sequence = opt_sequence.copy()         # Will contain the optimal sequence w/o the dummy job
    for k in range(len(job_sequence)):
        job_sequence[k] = job_sequence[k] - 1  # job indexes in opt_sequence are offseted by 1 due to the dummy job in index 0
    job_sequence.remove(-1)
    #print(f'One optimal sequence is {job_sequence}')
    
    makespan = model.getObjective().getValue()
    if not makespan.is_integer():
        raise ValueError("The makespan is not of integer value")
    
    return makespan, job_sequence    

def MIP_adverse_problem(job_sequence, nom_processing_times, max_dev, setup_times, budget):
    
    ##### Data #####
    nom_processing_times = reorder_sequence(job_sequence, nom_processing_times)
    max_dev = reorder_sequence(job_sequence, max_dev)
    nb_machines = len(nom_processing_times)
    nb_jobs = len(nom_processing_times[0])
    M = 1000             # An arbitrary large number

    sequence_wdummy = job_sequence.copy()
    for k in range(nb_jobs):
        sequence_wdummy[k] += 1
    sequence_wdummy.insert(0, 0)

    setup_mat = np.zeros((nb_machines, nb_jobs))
    for i in range(nb_machines):
        for j in range(nb_jobs):
            setup_mat[i][j] = setup_times[i][sequence_wdummy[j]][sequence_wdummy[j+1]]

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
    idle = model.addVars(nb_machines, nb_jobs, lb = 0.0, ub = GRB.INFINITY, name = 'idle')

    # Auxilary variable for definig the idle times
    z = model.addVars(nb_machines, nb_jobs, vtype=GRB.BINARY, name='z')

    ##### Constraints #####

    # Ensures that for one job, at least one idle time is null 

    model.addConstrs(idle[i,j] >= setup_mat[i][j] for i in range(nb_machines) for j in range(nb_jobs))

    for j in range(nb_jobs):
        model.addConstrs((idle[i,j] - setup_mat[i][j] <= M *(1- z[i,j]) 
                          for i in range(nb_machines)))
        model.addConstr(gp.quicksum(z[i,j] for i in range(nb_machines)) >= 1)

    # Budget constraint
    model.addConstr(gp.quicksum(khi[i,j] for i in range(nb_machines) for j in range(nb_jobs)) <= budget)

    # JAML constraints
    model.addConstrs(

        nom_processing_times[i, j] 
        + khi[i,j]*max_dev[i, j] 
        + idle[i,j] 
            == 
        nom_processing_times[i+1, j-1] 
        + khi[i+1, j-1]*max_dev[i+1, j-1]  
        + idle[i+1,j] 

            for i in range(nb_machines-1) 
            for j in range(1, nb_jobs))

    model.addConstrs(nom_processing_times[i,0] + khi[i,0]*max_dev[i,0] + idle[i,0] == idle[i+1,0] for i in range(nb_machines-1))

    ##### Objective #####

    # Maximize the makespan
    model.setObjective(
        gp.quicksum(
            idle[nb_machines-1, j] 
            for j in range(nb_jobs))
        + gp.quicksum(
            nom_processing_times[nb_machines-1, j] + khi[nb_machines-1, j]*max_dev[nb_machines-1, j] 
            for j in range(nb_jobs)),
            GRB.MAXIMIZE)

    ##### Optimization time ! #####
    model.optimize()

    ##### Retrieve data #####
    
    # Makespan
    worst_makespan = model.getObjective().getValue()

    # Retrieving the idle times in a matrix
    idle_mat = np.zeros(np.shape(nom_processing_times))
    for i in range(nb_machines):
        for j in range(nb_jobs):
            idle_mat[i, j] = idle[i, j].X
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
    #print(f'One critcal scenario is:\n{reorder_sequence(scen_sequence, scen)}')
    scen = reorder_sequence(scen_sequence, scen)

    #print(f'the setup times are \n {setup_mat}')
    #print(f'the makespan is {worst_makespan}')
    #print(f'Processing times :\n {pt}')

    #flowshopGantt.gantt_setup(pt, setup_mat, idle_mat, worst_makespan)

    return worst_makespan, scen

def CP_master_problem(scenario_pool, nom_processing_times, max_dev, setup_times):
    
    ##### Data #####
    nb_machines = len(nom_processing_times)
    nb_jobs = len(nom_processing_times[0])
    nb_scenariis = len(scenario_pool)

    for i in range(nb_machines):
        for j in range(nb_jobs+1):
             setup_times[i][j][j] = 0
             setup_times[i][j][0] = 0
    print(setup_times)

    ##### Initialize #####
    mdl = CpoModel()

    ##### Variables #####

    tasks = integer_var_list(nb_scenariis)
    sequences = []
    for s in range(nb_scenariis):
        processing_times = nom_processing_times + np.multiply(scenario_pool[s],max_dev) 
        processing_times = np.hstack((np.zeros((nb_machines,1)), processing_times))
        tasks[s] = [[
            interval_var(
                name='task_{}_{}_{}'.format(i,j,s), 
                size=int(processing_times[i,j])) 
                    for j in range(nb_jobs+1)] 
                    for i in range(nb_machines)]
        sequence = [sequence_var(tasks[s][i], name='sequence_{}_{}'.format(i,s)) for i in range(nb_machines)]
        sequences.append(sequence)

    Cmax = integer_var(name='Cmax')

    setup_mat = [transition_matrix(setup_times[i]) for i in range(nb_machines)]

    ##### Constraints #####

    for s in range(nb_scenariis):
        
        mdl.add(start_of(tasks[s][i][0]) == 0 for i in range(nb_machines))
        
        mdl.add(greater_or_equal(Cmax, end_of(tasks[s][nb_machines-1][j])) for j in range(nb_jobs+1))

        mdl.add(end_at_start(tasks[s][i][j], tasks[s][i+1][j]) for i in range(nb_machines-1) for j in range(nb_jobs+1))

        mdl.add(no_overlap(sequences[s][i], setup_mat[i]) for i in range(nb_machines))

        mdl.add(same_sequence(sequences[s][0], sequences[s][i]) for i in range(1, nb_machines))

        mdl.add(same_sequence(sequences[0][i], sequences[s][i]) for i in range(nb_machines))
        
    ##### Objective #####

    mdl.minimize(Cmax)

    ##### Solve ! #####

    msol = mdl.solve(LogVerbosity = 'Quiet')


    ##### Retrieve the solutions #####
    makespan = msol.get_value(Cmax)

    seqc_var = msol.get_value(sequences[0][0])

    def retrieve_sequence_position(var_name):
        return int(var_name.split('_')[2])
    opt_sequence = []
    for s in seqc_var:
        opt_sequence.append(retrieve_sequence_position(s.get_name()))
    print(opt_sequence)

    job_sequence = opt_sequence.copy()         # Will contain the optimal sequence w/o the dummy job
    for k in range(len(job_sequence)):
        job_sequence[k] = job_sequence[k] - 1  # job indexes in opt_sequence are offseted by 1 due to the dummy job in index 0
    job_sequence.remove(-1)
    print(f'one optimal sequence is: {job_sequence}')

    opt_sequence.remove(0)
    opt_sequence.insert(0, 0)
     
    slvd_setup_mat = np.zeros((nb_machines, nb_jobs))
    for i in range(nb_machines):
        for j in range(nb_jobs):
            slvd_setup_mat[i][j] = setup_times[i][opt_sequence[j]][opt_sequence[j+1]]
    print(slvd_setup_mat)

    return makespan, job_sequence
 
def CP_adverse_problem(job_sequence, nom_processing_times, max_dev, setup_times, budget):
    
    ##### Data #####
    nom_processing_times = reorder_sequence(job_sequence, nom_processing_times)
    max_dev = reorder_sequence(job_sequence, max_dev)
    nb_machines = len(nom_processing_times)
    nb_jobs = len(nom_processing_times[0])

    sequence_wdummy = job_sequence.copy()
    for k in range(nb_jobs):
        sequence_wdummy[k] += 1
    sequence_wdummy.insert(0, 0)

    setup_mat = np.zeros((nb_machines, nb_jobs))
    for i in range(nb_machines):
        for j in range(nb_jobs):
            setup_mat[i][j] = setup_times[i][sequence_wdummy[j]][sequence_wdummy[j+1]]


    mdl = CpoModel()

    Cmax = integer_var(name = "Cmax") # makespan

    tasks = [[mdl.interval_var(name='tasks_{}_{}'.format(i,j), size=int(nom_processing_times[i,j])) for j in range(nb_jobs)] for i in range(nb_machines)]

    setup_tasks = [[mdl.interval_var(name='SetUptasks_{}_{}'.format(i,j), size=int(setup_mat[i,j])) for j in range(nb_jobs)] for i in range(nb_machines)]

    dev_task = [[mdl.interval_var(name='deviations_{}_{}'.format(i,j), size=int(max_dev[i,j]), optional=True) for j in range(nb_jobs)] for i in range(nb_machines)]

    slack = mdl.integer_var_list(nb_jobs, name = "slack") # slack


    # Every SetUptasks of first job is scheduled at 0
    for i in range(nb_machines):
        mdl.add(start_of(setup_tasks[i][0]) == 0)


    #  precedence constraints SetUp -> task
    for j in range(nb_jobs):
        for i in range(nb_machines):
            mdl.add(end_before_start(setup_tasks[i][j], tasks[i][j]))


    #  Must stick constraints task -> Devtask
    for j in range(nb_jobs):
        for i in range(nb_machines):
            mdl.add(start_at_end(dev_task[i][j], tasks[i][j]))


    # precedence  constraints  setUp(i).start = max( task(i-1).end, DeviationTask(i-1).end)
    for j in range(1,nb_jobs):
        for i in range( nb_machines):
            mdl.add(start_of(setup_tasks[i][j]) == max(end_of(tasks[i][j-1]), end_of(dev_task[i][j-1])))
                

    # NoWait precedence constraints
    for j in range(nb_jobs):
        for i in range(1, nb_machines):
            mdl.add(start_of(tasks[i][j]) == max(end_of(tasks[i-1][j]) , end_of(dev_task[i-1][j]) ))



    # Alternative way to link jobs with  slack=0 at least one time
    for j in range(nb_jobs):
        mdl.add(slack[j] == min([start_of(tasks[i][j]) - end_of(setup_tasks[i][j]) for i in range(nb_machines)]))

    mdl.add(max([slack[j] for j in range(nb_jobs)]) == 0)


    # Budget constraint

    mdl.add(sum([presence_of(dev_task[i][j]) for j in range(nb_jobs) for i in range(nb_machines)]) <= budget)



    # Cmax contraint           
    mdl.add(Cmax == max([end_of(tasks[nb_machines-1][j]) for j in range(nb_jobs)] + [end_of(dev_task[nb_machines-1][j]) for j in range(nb_jobs)] ))

    mdl.maximize(Cmax)

    # Solve the model
    msol = mdl.solve(TimeLimit=10,LogVerbosity = 'Quiet')


    objective_SP = msol.get_value(Cmax) # on recupere la valeur optimale du makespan

    #sol_tasks = [[msol.get_value(tasks[i][j]) for j in range(nb_machines)] for i in range(nb_jobs)]
    var_sol = [[msol.get_var_solution(tasks[i][j]) for j in range(nb_jobs)] for i in range(nb_machines)]
    var_sol_SetUp = [[msol.get_var_solution(setup_tasks[i][j]) for j in range(nb_jobs)] for i in range(nb_machines)]
    var_sol_Dev = [[msol.get_var_solution(dev_task[i][j]) for j in range(nb_jobs)] for i in range(nb_machines)]
    #var_sol_deviation_nb= msol.get_value(cumul_deviation)
    ev_sol= [msol.get_value(slack[j]) for j in range(nb_jobs-1)]


    critical_scen = np.zeros((nb_machines,nb_jobs))

    for j in range(nb_jobs):
        for i in range(nb_machines):
            if var_sol_Dev[i][j].get_start():
                critical_scen[i,j]=1



    return objective_SP, critical_scen

def MIP_MIP_optimization(scenario_pool, nom_processing_times, max_dev, setup_times, budget):
    
    ##### Data #####
    nb_machines = len(nom_processing_times)
    nb_jobs = len(nom_processing_times[0])
    nb_scenariis = np.sum([2 ** math.comb(nb_machines * nb_jobs, gamma) for gamma in range(1, budget+1)])

    ##### Run first optimizations #####
    master_cmax, master_sequence = MIP_master_problem(scenario_pool, nom_processing_times, max_dev, setup_times)
    adverse_cmax, adverse_scenario = MIP_adverse_problem(master_sequence, nom_processing_times, max_dev, setup_times, budget)

    print(f'initial LB:{master_cmax} intial UB: {adverse_cmax}')
    while adverse_cmax > master_cmax:
        scenario_pool.append(adverse_scenario)
        print(f'Percentage of scenarios explored {(len(scenario_pool)/nb_scenariis)*100}')
        master_cmax, master_sequence = MIP_master_problem(scenario_pool, nom_processing_times, max_dev, setup_times)
        print(f'LB:{master_cmax} UB:{adverse_cmax}')
        adverse_cmax, adverse_scenario = MIP_adverse_problem(master_sequence, nom_processing_times, max_dev, setup_times, budget)
        
    
    print(f"Optimal makespan: {master_cmax}")
    print(f"Optimal sequence: {master_sequence}")

    return master_cmax, master_sequence

def CP_CP_optimization(scenario_pool, nom_processing_times, max_dev, setup_times, budget):
    
    ##### Data #####
    nb_machines = len(nom_processing_times)
    nb_jobs = len(nom_processing_times[0])
    nb_scenariis = np.sum([2 ** math.comb(nb_machines * nb_jobs, gamma) for gamma in range(1, budget+1)])

    ##### Run first optimizations #####
    master_cmax, master_sequence = CP_master_problem(scenario_pool, nom_processing_times, max_dev, setup_times)
    adverse_cmax, adverse_scenario = CP_adverse_problem(master_sequence, nom_processing_times, max_dev, setup_times, budget)

    print(f'initial LB:{master_cmax} intial UB: {adverse_cmax}')
    print(f'scenario pool \n {scenario_pool}')
    while adverse_cmax > master_cmax:
        scenario_pool.append(adverse_scenario)
        print(f'Percentage of scenarios explored {(len(scenario_pool)/nb_scenariis)*100}')
        master_cmax, master_sequence = CP_master_problem(scenario_pool, nom_processing_times, max_dev, setup_times)
        print(f'LB:{master_cmax} UB:{adverse_cmax}')
        adverse_cmax, adverse_scenario = CP_adverse_problem(master_sequence, nom_processing_times, max_dev, setup_times, budget)
        
    
    print(f"Optimal makespan: {master_cmax}")
    print(f"Optimal sequence: {master_sequence}")

    
    return master_cmax, master_sequence

def MIP_CP_optimization(scenario_pool, nom_processing_times, max_dev, setup_times, budget):
    
    ##### Data #####
    nb_machines = len(nom_processing_times)
    nb_jobs = len(nom_processing_times[0])
    nb_scenariis = np.sum([2 ** math.comb(nb_machines * nb_jobs, gamma) for gamma in range(1, budget+1)])

    ##### Run first optimizations #####
    master_cmax, master_sequence = MIP_master_problem(scenario_pool, nom_processing_times, max_dev, setup_times)
    adverse_cmax, adverse_scenario = CP_adverse_problem(master_sequence, nom_processing_times, max_dev, setup_times, budget)

    print(f'initial LB:{master_cmax} intial UB: {adverse_cmax}')
    print(f'scenario pool \n {scenario_pool}')
    while adverse_cmax > master_cmax:
        scenario_pool.append(adverse_scenario)
        print(f'Percentage of scenarios explored {(len(scenario_pool)/nb_scenariis)*100}')
        master_cmax, master_sequence = MIP_master_problem(scenario_pool, nom_processing_times, max_dev, setup_times)
        print(f'LB:{master_cmax} UB:{adverse_cmax}')
        adverse_cmax, adverse_scenario = CP_adverse_problem(master_sequence, nom_processing_times, max_dev, setup_times, budget)
        
    
    print(f"Optimal makespan: {master_cmax}")
    print(f"Optimal sequence: {master_sequence}")

    
    return master_cmax, master_sequence

def CP_MIP_optimization(scenario_pool, nom_processing_times, max_dev, setup_times, budget):
    
    ##### Data #####
    nb_machines = len(nom_processing_times)
    nb_jobs = len(nom_processing_times[0])
    nb_scenariis = np.sum([2 ** math.comb(nb_machines * nb_jobs, gamma) for gamma in range(1, budget+1)])

    ##### Run first optimizations #####
    master_cmax, master_sequence = CP_master_problem(scenario_pool, nom_processing_times, max_dev, setup_times)
    adverse_cmax = 10000
    print(f'initial LB:{master_cmax} intial UB: {adverse_cmax}')
    print(f'scenario pool \n {scenario_pool}')
    while adverse_cmax > master_cmax:
        adverse_cmax, adverse_scenario = MIP_adverse_problem(master_sequence, nom_processing_times, max_dev, setup_times, budget)
        scenario_pool.append(adverse_scenario)
        print(f'Percentage of scenarios explored {(len(scenario_pool)/nb_scenariis)*100}')
        master_cmax, master_sequence = CP_master_problem(scenario_pool, nom_processing_times, max_dev, setup_times)
        print(f'LB:{master_cmax} UB:{adverse_cmax}')

    
    print(f"Optimal makespan: {master_cmax}")
    print(f"Optimal sequence: {master_sequence}")

    
    return master_cmax, master_sequence

