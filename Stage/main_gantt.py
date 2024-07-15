
#This code will plot a simple Gantt diagramm representing a flowshop.
import numpy as np
import flowshop_gantt as flowshop_gantt
#Data



nom_processing_times = np.array([
    [5, 3, 6, 4, 9], 
    [4, 4, 5, 5, 1], 
    [7, 8, 10, 6, 5]
    ])

#Idle times, each row represent a machine, each column represent the idle time between the job j-1 and j 

idle = np.array([[1, 2, 4, 3],
                 [4, 0, 1, 2],
                 [8, 0, 0, 0]])

setup = np.array([[1, 1, 1, 1],
                  [1, 1, 1, 4],
                  [1, 1, 1, 1]])

max_dev = np.ones(np.shape(nom_processing_times))

scenario = np.zeros(np.shape(nom_processing_times))

makespan = 20

nb_jobs = len(nom_processing_times[0])
nb_machines = len(nom_processing_times)

#flowshopGantt.gantt(pt, idle, makespan, nb_machines, nb_jobs)

#flowshopGantt.ganttSetup(pt, setup, idle, makespan, nb_machines, nb_jobs)

flowshop_gantt.ganttRobust_setup_nowait(nom_processing_times, max_dev, scenario, setup, idle, makespan)


