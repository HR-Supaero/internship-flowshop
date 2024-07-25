# This code implements MIP and CP models of the no wait flowshop

import gurobipy as gp
from gurobipy import GRB

import numpy as np
import numpy.linalg as alg

from docplex.cp.model import *
from docplex.cp.model import CpoModel

import flowshop_gantt

import time

##### Helper functions #####

def reorder_sequence(
        job_sequence: np.ndarray, 
        mat: np.ndarray
) -> np.ndarray:
    """Reorders the columnns of a matrix along a given sequence"""
    return mat[:, job_sequence]

def cmax(
        processing_times: np.ndarray
) -> int:
    """Calculates the makespan of a nowait flowshop"""
    def nowait_completion_time(
            processing_times: np.ndarray,
            machine_index: int,
            job_index: int,
    ) -> int:
        """Calculates the completion time at a given job and machine indexes"""
        nb_machines = len(processing_times)

        # Base case
        if job_index == 0:
            return np.sum(processing_times[l, 0] 
                          for l in range(machine_index + 1))
        
        # General case
        
        overlaps = []

        overlaps.append(
            np.sum(processing_times[l, job_index-1] 
                   for l in range(1, nb_machines))
                   )

        for i in range(1, nb_machines-1):
            overlaps.append(
                np.sum(processing_times[l, job_index] for l in range(i)) 
                + np.sum(processing_times[l, job_index-1] for l in range(i+1, nb_machines))
                )
        
        overlaps.append(
            np.sum(processing_times[l, job_index]
                   for l in range(nb_machines-1))
                )
        
        min_overlap = np.min(overlaps)

        return (nowait_completion_time(processing_times, nb_machines-1, job_index-1) 
                + np.sum(processing_times[l, job_index]for l in range(machine_index+1)) 
                - min_overlap)
    
    nb_machines = len(processing_times)
    nb_jobs = len(processing_times[0])

    return nowait_completion_time(processing_times, nb_machines-1, nb_jobs-1)

##### Models #####

def MIP_model(
        processing_times: np.ndarray, 
        plotGantt: int
) -> None:
    """
     This function takes a matrix representing production times
     and returns an optimized model of the no wait flowshop

    Parameters:
     processing_times: the processing times time matrix
    
     plot_gantt: Define this value to 1 if wishing to plot the Gantt chart
    """
    
    if type(processing_times) != np.ndarray:
        raise TypeError("The processing time format is incorrect")

    ##### Data #####

    nb_jobs = len(processing_times[0])
    nb_machines = len(processing_times)

    ##### Model initialization #####

    m = gp.Model('No_wait_Flowshop')
    m.Params.IntegralityFocus = 1
    m.Params.OutputFlag = 0

    ##### Variables #####

    # x[j,k] = 1 if the j-th job is the k-th of the sequence
    #          0 otherwise
    x = m.addVars(nb_jobs, nb_jobs, vtype = GRB.BINARY, name = "x")

    # Auxilary variables, reorders the processing times along the sequence
    ord_processing_times = m.addVars(nb_machines, nb_jobs, lb = 0.0, name = "ord_processing_times")
    m.addConstrs(ord_processing_times[i,k] == gp.quicksum(x[j,k]*processing_times[i,j] for j in range(nb_jobs)) 
                 for i in range(nb_machines) 
                 for k in range(nb_jobs))

    # idle times after the processing of a job
    idle = m.addVars(nb_machines, nb_jobs, lb = 0.0, ub = GRB.INFINITY, name = "idle")

    ##### Constraints #####

    # Each position in the sequence is occupied by only one job
    m.addConstrs(gp.quicksum(x[j,k] for j in range(nb_jobs)) == 1 for k in range(nb_jobs))
    # A job only appears once in a sequence
    m.addConstrs(gp.quicksum(x[j,k] for k in range(nb_jobs)) == 1 for j in range(nb_jobs))

    # job-adjency-machine-linkage (JAML) constraints
    m.addConstrs(
        ord_processing_times[i, k] + idle[i, k] ==
        ord_processing_times[i+1, k-1] + idle[i+1, k]
        for i in range(nb_machines-1) 
        for k in range(1, nb_jobs))

    m.addConstrs(
        ord_processing_times[i, 0] + idle[i, 0] ==
        idle[i+1, 0]
        for i in range(nb_machines-1))

    ##### Objective function #####

    m.setObjective(
        gp.quicksum(idle[nb_machines-1, j] for j in range(nb_jobs)) 
        + gp.quicksum(ord_processing_times[nb_machines-1, j] for j in range(nb_jobs)),
        GRB.MINIMIZE
        )


    ##### Optimization #####
    m.optimize()

    ##### Retrieving optimized data #####

    # Retrieve the optimal sequence
    var_list = []
    for v in m.getVars():
        if v.X != 0:
            if (v.VarName[0] == "x"):
                var_list.append(v.VarName)
    def retrive_sequence_position(var_name):
        return int(var_name.split(',')[1].split(']')[0])
    def retrieve_job_number(var_name):
       return int(var_name.split(',')[0].split('[')[1])
    
    opt_sequence = sorted(var_list, key = retrive_sequence_position)
    for k in range(nb_jobs):
        opt_sequence[k] = retrieve_job_number(opt_sequence[k])

    # Retrieve makespan
    makespan = m.getObjective().getValue()
    if not makespan.is_integer():
        raise ValueError("The makespan is not of integer value")

    print(f'One optimal sequence is: ', opt_sequence)
    print("The makespan is ",makespan)

    # Reordering the processing times along the given optimal sequence

    print(f'given processing times:\n', processing_times)    
    processing_times = reorder_sequence(opt_sequence, processing_times)
    print(f'reordered processing times\n', processing_times)

    #Retrieving the idle times
    idle_times = np.zeros(np.shape(processing_times))
    for i in range(nb_machines):
        for j in range(nb_jobs):
            idle_times[i, j] = idle[i, j].X
    print(f'idle times\n', idle_times)

    #Plotting the Gantt chart
    if plotGantt:
        flowshop_gantt.gantt(processing_times, idle_times, makespan)


def CP_model(
        processing_times: np.ndarray,
        plot_gantt: int
) -> None:
    
    ##### Data #####

    nb_machines = len(processing_times)
    nb_jobs = len(processing_times[0])

    ##### Initialize the model #####
    mdl = CpoModel()

    ##### Variables #####
    
    # Interval variable for the processing times
    tasks = [[interval_var(name = 'tasks_{}_{}'.format(i,j), size = int(processing_times[i,j])) for j in range(nb_jobs)] for i in range(nb_machines)]
    
    # Sequence variable on each machine
    sequence = [sequence_var(tasks[i], name = 'sequence_{}'.format(i)) for i in range(nb_machines)]

    # Makespan
    Cmax = integer_var(name = 'Cmax')

    ##### Constraints #####

    # Makespan greater than completion of every job on last machine
    mdl.add(greater_or_equal(Cmax, end_of(tasks[nb_machines-1][j])) for j in range(nb_jobs))

    # No wait linkage constraint
    mdl.add(end_at_start(tasks[i][j], tasks[i+1][j]) for i in range(nb_machines-1) for j in range(nb_jobs))

    # Tasks do not overlap on the same machine
    mdl.add(no_overlap(sequence[i]) for i in range(nb_machines))

    # Sequence of jobs identical on each machine
    mdl.add(same_sequence(sequence[0], sequence[i]) for i in range(1, nb_machines)) 

    ##### Objective #####
    mdl.minimize(Cmax)

    ##### Solving #####
  
    msol = mdl.solve()

    ##### Retrieve the solutions #####
    makespan = msol.get_value(Cmax)

    # Retrieve the optimal sequence
    seqc_var = msol.get_value(sequence[0])
    def retrieve_sequence_position(var_name):
        return int(var_name.split('_')[2])
    opt_sequence = []
    for s in seqc_var:
        opt_sequence.append(retrieve_sequence_position(s.get_name()))
    print(opt_sequence)
    
    # Retrieve the solved processing times
    slvd_tasks = [[msol.get_var_solution(tasks[i][j]) for j in range(nb_jobs)] for i in range(nb_machines)]
    
    if plot_gantt:
        flowshop_gantt.ganttCP(slvd_tasks, makespan)


    

