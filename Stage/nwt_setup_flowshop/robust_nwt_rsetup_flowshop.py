from itertools import combinations, permutations
import gurobipy as gp
import random
import numpy as np
import numpy.linalg as alg
from gurobipy import GRB
import flowshop as flowshop
from docplex.cp.model import *
import nowait_setup_flowshop

import math

def is_scenario_in_list(scenario, scenario_pool):
    """
    Check if an array (list) is present in a list of arrays (lists).

    Parameters:
    array: The array (list) to check for in the list of arrays (lists).
    list_of_arrays: The list of arrays (lists) in which to check for the array.

    Returns:
    bool: True if the array is in the list of arrays, False otherwise.
    """
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

def calculate_setup_times(scenario, nom_setup_times, max_dev):
    nb_machines = len(nom_setup_times)
    nb_jobs = len(nom_setup_times[0])
    setup_times = {i: np.zeros((nb_jobs+1, nb_jobs+1)) for i in range(nb_machines)}
    for i in range(nb_machines):
        setup_times[i] = nom_setup_times[i] + np.multiply(scenario[i], max_dev[i])
    return setup_times

def setup_dict_to_mtrx(sequence, setup_dict):
    """Morphs under a sequence, a setup dict into a matrix"""


    nb_jobs = len(sequence)
    nb_machines = len(setup_dict)
    sequence = [x + 1 for x in sequence]
    sequence = [0] + sequence

    setup_mtrx = np.zeros((nb_machines, nb_jobs))
    for i in range(nb_machines):
        for j in range(nb_jobs):
            setup_mtrx[i][j] = setup_dict[i][sequence[j]][sequence[j+1]]
    
    return setup_mtrx

def generate_all_process_scenariis(
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

def setup_scenariis_for_seq(sequence, budget, nb_machines, nb_jobs) -> dict[np.ndarray]:
    """Generates all setup scenariis for a given sequence"""
    list_scen = []

    # Prepend 0 and add 1 to every element to include the dummy job
    sequence = [x + 1 for x in sequence]
    sequence = [0] + sequence

    for gamma in range(budget + 1):
            possible_ones = list(combinations(range(nb_machines * nb_jobs), gamma))
            for comb in possible_ones:
                deviations = {i: np.zeros((nb_jobs + 1, nb_jobs + 1)) for i in range(nb_machines)}
                ones_mtrx = np.zeros((nb_machines, nb_jobs))

                for index in comb:
                    machine_ind = index // nb_jobs
                    seq_ind = index % nb_jobs
                    ones_mtrx[machine_ind, seq_ind] = 1

                for i in range(nb_machines):
                    for k in range(nb_jobs):
                        deviations[i][sequence[k]][sequence[k + 1]] = ones_mtrx[i][k]
                list_scen.append(deviations)
    
    return list_scen


def generate_all_setup_scenariis(budget, nb_machines, nb_jobs):
    list_scen = []
    
    # Generate permutations of the sequence excluding the first element (0)
    sequences = list(itertools.permutations(range(nb_jobs)))

    for seq in sequences:
        scen_for_seq = setup_scenariis_for_seq(seq, budget, nb_machines, nb_jobs)
        list_scen = list_scen + scen_for_seq

    return list_scen
                    

def generate_some_scenariis(
        nb_scenariis,
        budget, 
        nb_machines, 
        nb_jobs
):
    
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


def checker(worst_cmax, opt_sequence, nom_process_times, process_dev, nom_setup_times, setup_dev, process_budget, setup_budget):
    """
    Checks for a given solution:\n
    if other makespan under the optimal sequence are lower than the worst cmax\n
    if for all other sequences, there exists a longer makespan

    Parameters:

    worst_cmax: the worst case makespan

    opt_sequence: one optimal job sequence

    nom_process_times: the nominal process times

    process_dev: the deviation values for the processes

    nom_setup_times: the nominal setup times

    setup_dev: the deviation values for the setups

    process_budget: the budget for the process

    setup_budget: the budget for the setups

    Returns:

    True if the given cmax and sequence form a robust solution.
    False otherwise
    """

    ##### Data #####

    nb_machines = len(nom_process_times)
    nb_jobs = len(nom_process_times[0])
    process_scenariis = generate_all_process_scenariis(process_budget, nb_machines, nb_jobs)
    setup_scenariis = setup_scenariis_for_seq(opt_sequence, setup_budget, nb_machines, nb_jobs)
    scenariis = list(itertools.product(process_scenariis, setup_scenariis))

    def check_all_scenariis(worst_cmax, opt_sequence, nom_process_times, process_dev, nom_setup_times, setup_dev, process_budget, setup_budget):
        """Checks if all scenariis under the given sequence lead to a lower cmax"""
    
        
        nom_process_times = reorder_sequence(opt_sequence, nom_process_times)
        process_dev = reorder_sequence(opt_sequence, process_dev)
        
        nom_setup_times = setup_dict_to_mtrx(opt_sequence, nom_setup_times)
        setup_dev = setup_dict_to_mtrx(opt_sequence, setup_dev)
    
        for (process_scen, setup_scen) in scenariis:
            process_scen = reorder_sequence(opt_sequence, process_scen)
            setup_scen = setup_dict_to_mtrx(opt_sequence, setup_scen)
            process_time = calculate_processing_times(process_scen, nom_process_times, process_dev)
            setup_time = calculate_processing_times(setup_scen, nom_setup_times, setup_dev)
            cmax = nowait_setup_flowshop.cmax(process_time, setup_time)
            if cmax > worst_cmax:
                print(f'The scenario \n {process_scenariis, setup_scenariis} \n violates the worst case makespan with a makespan of {cmax}')
                return False
        return True

    
    
    def find_worst_seq(worst_cmax, opt_sequence, nom_process_times, process_dev, nom_setup_times, setup_dev, process_budget, setup_budget):
        
        """Find among every other sequences a worst makespan"""
        other_sequences = list(permutations(range(nb_jobs)))
        other_sequences.remove(tuple(opt_sequence))

        longer_cmax_found = False
        worst_cmax_for_seq = 0

        for seq in other_sequences:
            
            nom_process_mtrx = reorder_sequence(seq, nom_process_times)
            process_dev_mtrx = reorder_sequence(seq, process_dev)
            
            nom_setup_mtrx = setup_dict_to_mtrx(seq, nom_setup_times)
            setup_dev_mtrx = setup_dict_to_mtrx(seq, setup_dev)

            for (process_scen, setup_scen) in scenariis:
                process_scen = reorder_sequence(seq, process_scen)
                setup_scen = setup_dict_to_mtrx(seq, setup_scen)
                process_times = calculate_processing_times(process_scen, nom_process_mtrx, process_dev_mtrx)
                setup_times = calculate_processing_times(setup_scen, nom_setup_mtrx, setup_dev_mtrx)
                cmax = nowait_setup_flowshop.cmax(process_times, setup_times)
                #print(f'The makespan for {seq} under scenario \n {(process_scen, setup_scen)} \n is {cmax}')   
                if cmax >= worst_cmax:
                    longer_cmax_found = True
                    worst_cmax_for_seq = cmax
                    break
            if not longer_cmax_found:
                #print(f'All makespans for sequence {seq} are shorter than {worst_cmax}.')
                return False
            
        return True

    if not check_all_scenariis(worst_cmax, opt_sequence, nom_process_times, process_dev, nom_setup_times, setup_dev, process_budget, setup_budget):
        return False

    if not find_worst_seq(worst_cmax, opt_sequence, nom_process_times, process_dev, nom_setup_times, setup_dev, process_budget, setup_budget):
        print('One sequence leads to a lower worts case makespan')
        return False 

    return True



    
    


def MIP_model_extended(nom_processing_times, max_dev, nom_setup_times, setup_dev, process_budget, setup_budget):

    ##### Data #####

    nb_machines = len(nom_processing_times)
    nb_jobs = len(nom_processing_times[0])
    A = 100000


    process_scenarios = generate_all_process_scenariis(process_budget, nb_machines, nb_jobs)
    nb_process_scen = len(process_scenarios)

    setup_scenarios = generate_all_setup_scenariis(setup_budget, nb_machines, nb_jobs)
    nb_setup_scen = len(setup_scenarios)

    scenariis =list(itertools.product(process_scenarios, setup_scenarios))
    nb_scenariis = len(scenariis)

    processing_times = []       # Array containing every possible processing times
    for scen in scenariis:
        pt_s = calculate_processing_times(scen[0], nom_processing_times, max_dev)
        processing_times.append(pt_s.copy())

    setup_times = []
    for scen in scenariis:
        st_s = calculate_setup_times(scen[1], nom_setup_times, setup_dev)
        setup_times.append(st_s.copy())




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
        model.addConstrs(idle[i,j,k, s] >= setup_times[s][i][j][k] 
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
    print(nb_scenariis, nb_process_scen, nb_setup_scen)
    
    for v in model.getVars():
        if v.VarName[0] == 'z':
            if v.X != 0:
                print(v)
        if v.VarName[0] == 'd':
            if v.X != 0:
                print(v)    
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
    print(f'One optimal sequence is {job_sequence}')
    
    makespan = model.getObjective().getValue()
    if not makespan.is_integer():
        raise ValueError("The makespan is not of integer value")
    print(f'The makespan is {makespan}')

    return makespan, job_sequence

