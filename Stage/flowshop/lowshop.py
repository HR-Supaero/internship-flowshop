
# This file contains implementations for the standard flowshop problem

import gurobipy as gp
from gurobipy import GRB

import numpy as np

import flowshop_gantt

##### Helper functions #####

def reorder_sequence(
        job_sequence: np.ndarray, 
        mat: np.ndarray
) -> np.ndarray:
    """Reorders the columnns of a matrix along a given sequence"""
    return mat[:, job_sequence]

def retrieve_permutation(
        matA: np.ndarray, 
        matB: np.ndarray
) -> np.ndarray:
    """Returns the permutation of columns of matA to matB"""
    permut = []
    nb_jobs = len(matA[0])
    nb_machines = len(matA)
    for j in range(nb_jobs):
        for k in range(nb_jobs):
                if k not in permut:
                    i = 0               
                    while i < nb_machines and matA[i, j] == matB[i, k]:
                        i+=1
                    if i == nb_machines:
                        permut.append(k)
                        break
    return permut

def cmax(
        processing_times: np.ndarray
) -> int:
    """Calculates the makespan of a given flowshop"""
    def completion_time(
            machine_indx: int, 
            job_indx: int, 
            processing_times: np.ndarray
    ) -> int:
        """Calculates the completion time at a given job and machine indexes"""

        # Base case for the first job in the sequence
        if job_indx == 0:
            return sum(processing_times[i, 0] for i in range(machine_indx+1))
        
        # Base case for the first machine in the flowshop
        if machine_indx == 0:
            return sum(processing_times[0, j] for j in range(job_indx+1))
        
        return max(
            completion_time(machine_indx-1, job_indx, processing_times), 
            completion_time(machine_indx, job_indx - 1, processing_times)) + processing_times[machine_indx, job_indx]
    
    nb_machines = len(processing_times)
    nb_jobs = len(processing_times[0])
    
    return completion_time(nb_machines-1, nb_jobs-1, processing_times)

##### Model implementations #####

def MIP_model(
        processing_times: np.ndarray, 
        plot_gantt: int
) -> tuple[gp.Model, np.ndarray, np.ndarray]:
    """
     This function takes a matrix representing production times
     and returns an optimized model of the flowshop
    
     Parameters:

     processing_times: the processing times matrix, 
     plot_gantt: Define this value to 1 if wishing to plot the Gantt chart
    
    """

    #print(f'given processing times:\n', processing_times)
    if type(processing_times) != np.ndarray:
        raise TypeError("The processing times format is incorrect")


    ##### Data #####

    nb_machines = len(processing_times)
    nb_jobs = len(processing_times[0])


    ##### Initializing the model #####

    m = gp.Model('standard_flowshop')
    m.Params.IntegralityFocus = 1
    m.Params.OutputFlag = 1


    ##### Variables #####

    # x[j,k] = 1 if the j-th job is the k-th of the sequence
    #          0 otherwise
    x = m.addVars(nb_jobs, nb_jobs, vtype = GRB.BINARY, name = "x")

    # Auxilary variables, reorders the processing times along the sequence
    ord_processing_times = m.addVars(nb_machines, nb_jobs, lb = 0.0, name = "ord_processing_times")
    m.addConstrs(ord_processing_times[i,k] == gp.quicksum(x[j,k]*processing_times[i,j] for j in range(nb_jobs)) 
                 for i in range(nb_machines) 
                 for k in range(nb_jobs))


    # waiting times after the processing of a job
    wait = m.addVars(nb_machines-1, nb_jobs, lb = 0.0, ub = GRB.INFINITY, name = "wait")

    # idle times after the processing of a job
    idle = m.addVars(nb_machines, nb_jobs-1, lb = 0.0, ub = GRB.INFINITY, name = "idle")

    ##### Constraints #####

    # Each position in the sequence is occupied by only one job
    m.addConstrs(gp.quicksum(x[j,k] for j in range(nb_jobs)) == 1 for k in range(nb_jobs))
    # A job only appears once in a sequence
    m.addConstrs(gp.quicksum(x[j,k] for k in range(nb_jobs)) == 1 for j in range(nb_jobs))

    # The waiting time for the first job on each machine is null
    m.addConstrs(wait[i, 0] == 0 for i in range(nb_machines-1))

    # The idle times on the first machine for each job is null
    m.addConstrs(idle[0, k] == 0 for k in range(nb_jobs-1))

    # job-adjency-machine-linkage (JAML) constraints
    m.addConstrs(
        ord_processing_times[i, k+1] + wait[i,k+1] + idle[i,k] == 
        ord_processing_times[i+1, k] + wait[i,k] + idle[i+1,k] 
        for i in range(nb_machines-1) 
        for k in range(nb_jobs-1))


    ##### Objective function #####

    # Calculates the makespan of the flowshop
    m.setObjective(
        gp.quicksum(ord_processing_times[i,0] for i in range(nb_machines-1))
        + gp.quicksum(idle[nb_machines-1, j] for j in range(nb_jobs-1)) 
        + gp.quicksum(ord_processing_times[nb_machines-1, j] for j in range(nb_jobs)),
        GRB.MINIMIZE
        )
    

    ##### Optimization #####

    m.optimize()

    ##### Retrieve optimized data #####

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
        
    #print(f'One optimal sequence is: ', opt_sequence)

    # Retrieve the makespan
    makespan = m.getObjective().getValue()    
    if not makespan.is_integer():
        raise ValueError("The makespan is not of integer value")
    #print("The makespan is ", makespan)

    # Reorders the processing times along the given optimal sequence
    processing_times = reorder_sequence(opt_sequence, processing_times)
    #print(f'reordered processing times\n', processing_times)


    # Retrieving the idle times in a matrix
    idle_mat = np.zeros(np.shape(processing_times))
    for i in range(1, nb_machines):
        idle_mat[i,0] = sum(processing_times[k,0] for k in range(i))
    for i in range(nb_machines):
        for j in range(1, nb_jobs):
            idle_mat[i, j] = idle[i, j-1].X
    #print(f'idle times\n', idle_mat)


    # Plotting the gantt chart
    if plot_gantt:
        flowshop_gantt.gantt(processing_times, idle_mat, makespan)

    return m, processing_times, idle_mat

































        



                
    
 



