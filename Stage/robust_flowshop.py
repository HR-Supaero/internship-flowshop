# helllo it is time for robust flwshop !
from itertools import combinations, permutations
import gurobipy as gp
import random
import numpy as np
import numpy.linalg as alg
from gurobipy import GRB
import flowshop as flowshop
import gantt.flowshop_gantt as flowshop_gantt

##### Helper functions #####

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

def check_criticalness(
        nom_processing_times,
        max_deviations, 
        job_sequence, 
        budget, 
        worst_makespan):
    # This function, given a sequence of jobs, checks if all the scenario fitting this sequence have a makespan lower than the worst case makespan
    nb_machines = len(nom_processing_times)
    nb_jobs = len(nom_processing_times[0])

    # Reorder the production times and deviation values along the job sequence
    opt_nominal_pt = reorder_sequence(job_sequence, nom_processing_times)
    opt_max_dev = reorder_sequence(job_sequence, max_deviations)

    # For the given sequence, check all possible scenarios if they fit the worstcase makespan
    scenariis = generate_all_scenariis(budget, nb_machines, nb_jobs)
    for scen in scenariis:
        scen = reorder_sequence(job_sequence, scen)
        optpt = calculate_processing_times(scen, opt_nominal_pt, opt_max_dev)
        Cmax = flowshop.nowait_cmax(optpt)
        #print(f'The makespan is {Cmax} for scenario \n {scen}')
        if Cmax > worst_makespan:
            print(f'The scenario \n {scen} \n violates the worst case makespan with a makespan of {Cmax}')
            return 0 
        else:
            print('yay !')
            return 1 

##### Models #####

def MIP_extended(
        nom_processing_times: np.ndarray, 
        max_deviations: np.ndarray, 
        budget: int
) -> tuple[gp.Model, np.ndarray, int]:
    """
    This function optimizes a robust flowshop model following an extended
    formulation.
    Meaning all possible scenarios under the budget are generated    

    nom_processing_times: the nominal processing times time matrix, 
    
    max_deviations: the maximum deviation values matrix
    
    budget: the uncertainty budget
    
    returns: an optimized gurobipy model 
             one optimal sequence
             the worst makespan under the optimal sequence
    
    """

    ##### Data #####

    nb_machines = len(nom_processing_times)
    nb_jobs = len(nom_processing_times[0])
    scenariis = generate_all_scenariis(budget, nb_machines, nb_jobs)
    nb_scenariis = len(scenariis)
 
    processing_times = []       # Array containing every possible processing times
    for scen in scenariis:
        pt_s = calculate_processing_times(scen, nom_processing_times, max_deviations)
        processing_times.append(pt_s.copy())


    ###### Model #####
    m = gp.Model('Extended_Robust_Flowshop')
    m.Params.IntegralityFocus = 1
    m.Params.OutputFlag = 1


    ###### Variables #####

    # Job position variable
    x = m.addVars(nb_jobs, nb_jobs, vtype = GRB.BINARY, name = 'x')

    # Auxilary variable to reorder processing times along the sequence
    ord_processing_times = m.addVars(nb_machines, nb_jobs, nb_scenariis, name = 'ord_processing_times')
    for s in range(nb_scenariis):
        m.addConstrs(ord_processing_times[i, k, s] 
                     ==
                      gp.quicksum(
                          processing_times[s][i][j]*x[j,k] 
                          for j in range(nb_jobs)) 
                     for k in range(nb_jobs)
                     for i in range(nb_machines) )

    # Idle time on machine i BEFORE the processing of job j
    idle = m.addVars(nb_machines, nb_jobs, nb_scenariis, lb =0.0, name = 'idle')

    # Waiting time of the jth job after its process on machine i
    wait = m.addVars(nb_machines-1, nb_jobs, nb_scenariis, lb = 0.0, name = 'wait')

    # Makespan
    makespan = m.addVar(lb =0.0, vtype = GRB.CONTINUOUS, name = 'Cmax')

    ##### Constraints #####

    # Each position in the sequence is occupied by only one job
    m.addConstrs(gp.quicksum(x[j, k] for j in range(nb_jobs)) == 1 
                 for k in range(nb_jobs))
    # A job appears only once in the sequence
    m.addConstrs(gp.quicksum(x[j, k] for k in range(nb_jobs)) == 1 
                 for j in range(nb_jobs))

    for s in range(nb_scenariis):

        # JAML constraints
        m.addConstrs(ord_processing_times[i, k+1, s] 
                     + wait[i, k+1, s] 
                     + idle[i, k+1, s] 
                     == 
                     ord_processing_times[i+1, k, s] 
                     + wait[i, k, s] 
                     + idle[i+1, k+1, s] 
                     for i in range(nb_machines-1) 
                     for k in range(nb_jobs-1))
        
        # Idle time for the first job in the sequence        
        m.addConstrs(idle[i, 0, s] ==  
                     gp.quicksum(ord_processing_times[l, 0, s] for l in range(i))
                     for i in range(nb_machines))

        # Makespan definition
        m.addConstr(makespan >= 
                    gp.quicksum(
                        ord_processing_times[nb_machines-1, k, s] 
                        + idle[nb_machines-1, k, s] 
                        for k in range(nb_jobs)))

        # Idle times on first machine are null
        m.addConstrs(idle[0, k, s] == 0 
                     for k in range(nb_jobs))

        # Waiting times for the first job are null
        m.addConstrs(wait[i, 0, s] == 0 
                     for i in range(nb_machines-1))


    ##### Objective function #####
    m.setObjective(makespan, GRB.MINIMIZE)

    ##### Optimization #####
    m.optimize()

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
    #print(f'one optimal sequence is {opt_sequence}')

    # Retrieve the worst case makespan
    worst_makespan = m.getObjective().getValue()
    #print(f'For a budget of {budget} the worst case makespan is {makespan}')

    #Retrieving for each scenario, the process times, and the idle times
    idle_list = []
    permu_mat = np.zeros((nb_jobs, nb_jobs))
    for j in range(nb_jobs):
        permu_mat[j, opt_sequence[j]] = 1
    permu_mat = alg.inv(permu_mat)

    for s in range(nb_scenariis):

        # Reordering the processing times along the given optimal sequence
        for i in range(nb_machines):
            for j in range(nb_jobs):
                processing_times[s][i][j] = ord_processing_times[i, j ,s].X


        # Retrieving the idle times
        idle_list.append(np.zeros(np.shape(processing_times[s])))
        for i in range(nb_machines):
            for j in range(nb_jobs):
                idle_list[s][i][j] = idle[i, j,s].X

    #print(f'The given processing times are \n {nom_processing_times} \n the reordered processing times are \n {np.dot(nom_processing_times, permu_mat)}')
    #print(f'The given deviation values are \n {max_deviations} \n the reordered deviation values are \n {np.dot(max_deviations, permu_mat)}')

    return m, opt_sequence, worst_makespan

def master_callback(model, where):
    """
    Callback function for the master problem of the flowshop formulated in MIP.
    'model' and 'where' arguments are passed by gurobipy when the callback
    is invoked (enfin je crois j'ai pas vérifié).
    The callback is set to stop the optimization process when one
    feasable solution has been found during the MIP solving.
    """
    if where == GRB.Callback.MIPSOL:
        print('Found a new incumbent solution')
        MIPSOLsolcnt = model.cbGet(GRB.Callback.MIPSOL_SOLCNT)
        MIPSOLobj = model.cbGet(GRB.Callback.MIPSOL_OBJ)
        MIPSOLobjbst = model.cbGet(GRB.Callback.MIPSOL_OBJBST)
        MIPSOLobjbnd = model.cbGet(GRB.Callback.MIPSOL_OBJBND)
        MIPSOLnodecnt = model.cbGet(GRB.Callback.MIPSOL_NODCNT)
        print(f'There are {MIPSOLsolcnt} solutions')
        if MIPSOLsolcnt and MIPSOLnodecnt >= 1:
            print('aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa')
            vars = model.getVars()
            try:
                print(f'{vars[0]} {vars[0].Xn}')
                print(model.Status == GRB.INTERRUPTED)
            except AttributeError:
                print(f'Couldnt access the variables values')
            except Exception as e:
                print('Cheh', e)

        print(model.Status == GRB.INTERRUPTED)
        print(model.Status == GRB.OPTIMAL)
        #for v in model.getVars():
            #print('yo')
                    
            #print(f'{v.VarName}{type(v)}')
            #print(f'{v}: {v.Xn}')

def master_problem(
    processing_times, 
    plot_gantt):
    """"
     This function takes a matrix representing production times
     and returns an optimized model of the flowshop
    
     processing_times: the processing times time matrix, 
                       the rows represent the machines
                       the columns the jobs
    
     plot_gantt: Define this value to 1 if wishing to plot the Gantt chart
    
     returns: an optimized gurobipy model 
             the reordered processing times along the optimal sequence
              the idle times matrix
    """

    print(f'given processing times:\n', processing_times)

    if type(processing_times) != np.ndarray:
        raise TypeError("The processing times format is incorrect")

    nb_machines = len(processing_times)
    nb_jobs = len(processing_times[0])

    # Initializing the model
    m = gp.Model('robust_flowshop')
    m.Params.IntegralityFocus = 1
    m.Params.OutputFlag = 0
    m.Params.MIPFocus = 1       # Focus on finding feasible solutions
    #m.Params.LogFile = "standard_log"


    # Variables

    # Decision variable for the sequence 
    # x[j, k] = 1 if the j-th job is the k-th of the sequence
    #          0 otherwise
    x = m.addVars(nb_jobs, nb_jobs, vtype = GRB.BINARY, name = "x")

    # Auxilary variable to reorder the processing times along the sequence
    ordered_processing_times = m.addVars(nb_machines, nb_jobs, lb = 0.0, name = "ordered_processing_times")
    m.addConstrs(ordered_processing_times[i,k] == gp.quicksum(x[j,k]*processing_times[i,j] for j in range(nb_jobs)) 
                 for i in range(nb_machines) 
                 for k in range(nb_jobs))


    # wait[i, k] = time difference between the end of job k on machine i and its start on machine i+1
    wait = m.addVars(nb_machines-1, nb_jobs, lb = 0.0, ub = GRB.INFINITY, name = "wait")

    # idle[i][k] = time on machine i between the end of job k and the start of job k+1
    idle = m.addVars(nb_machines, nb_jobs-1, lb = 0.0, ub = GRB.INFINITY, name = "idle")

    # Constraints

    # Each position in the sequence is occupied by only one job
    m.addConstrs(gp.quicksum(x[j,k] for j in range(nb_jobs)) == 1 
                 for k in range(nb_jobs))
    
    # A job only appears once in a sequence
    m.addConstrs(gp.quicksum(x[j,k] for k in range(nb_jobs)) == 1 
                 for j in range(nb_jobs))

    # The waiting time for the first job on each machine is null
    m.addConstrs(wait[i, 0] <= 0 
                 for i in range(nb_machines-1))

    # The idle times on the first machine for each job is null
    m.addConstrs(idle[0, k] <= 0
                 for k in range(nb_jobs-1))

    # job-adjency-machine-linkage (JAML) constraints
    m.addConstrs(ordered_processing_times[i, k+1] + wait[i,k+1] + idle[i,k] == 
                 ordered_processing_times[i+1, k] + wait[i,k] + idle[i+1,k] 
                 for i in range(nb_machines-1) 
                 for k in range(nb_jobs-1))

    # Objective function
    # Calculates the makespan of the flowshop
    m.setObjective(gp.quicksum(ordered_processing_times[i,0] for i in range(nb_machines-1))
                 + gp.quicksum(idle[nb_machines-1, j] for j in range(nb_jobs-1)) 
                 + gp.quicksum(ordered_processing_times[nb_machines-1, j] 
                                for j in range(nb_jobs)),
                                GRB.MINIMIZE)
    
    #m.optimize()
    m.optimize(master_callback)

    # Prints the optimal sequence and the makespan
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

    print(f'One optimal sequence is: ', opt_sequence)
    makespan = m.getObjective().getValue()
    
    if not makespan.is_integer():
        raise ValueError("The makespan is not of integer value")
    
    #print("The makespan is ", makespan)

    # Reorders the processing times along the given optimal sequence
    permu_mat = np.zeros((nb_jobs, nb_jobs))
    for j in range(nb_jobs):
        permu_mat[j, opt_sequence[j]] = 1
    permu_mat = alg.inv(permu_mat)
    processing_times = np.dot(processing_times, permu_mat)
    print(f'reordered processing times\n', processing_times)

    # Retrieving the idle times in a matrix
    idle_mat = np.zeros(np.shape(processing_times))
    for i in range(1, nb_machines):
        idle_mat[i,0] = sum(processing_times[k,0] for k in range(i))
    for i in range(nb_machines):
        for j in range(1, nb_jobs):
            idle_mat[i, j] = idle[i, j-1].X
    print(f'idle times\n', idle_mat)


    # Plotting the gantt chart
    if plot_gantt:
        flowshop_gantt.gantt(processing_times, idle_mat, makespan)
    

    return m, processing_times, idle_mat
