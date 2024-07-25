#Definig the function that plots the Gantt Diagramm associated with a flowshop

import numpy as np
from docplex import *
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors



def gantt(
        production_times: np.ndarray, 
        idle_times: np.ndarray, 
        makespan: int
) -> None:
    """
    This function plots the Gantt chart associated to a flowshop

     production_times: The processing times

     idle_times: the idle times matrix

     makespan: the makespan of the flowshop
    """

    nb_machines = len(production_times) 
    nb_jobs = len(production_times[0])

    # Defining the canavas
    color_names = list(mcolors.TABLEAU_COLORS)

    plt.rcParams["figure.figsize"] = (14, 3)
    plt.rcParams["figure.dpi"] = 200
            
    # Declaring a figure "gnt"
    fig, gnt = plt.subplots()

    # Setting Y-axis limits
    gnt.set_ylim(0, nb_jobs + 1)
            
    # Setting X-axis limits
    gnt.set_xlim(0, makespan + 5)
            
    # Setting labels for x-axis and y-axis
    gnt.set_xlabel('Time')
    gnt.set_ylabel('Machine')
    #gnt.legend(title='alpha='+str(alpha))
            
    # Setting ticks on y-axis
    gnt.set_yticks([k for k in range(1, nb_machines+1)])
    # Labelling tickes of y-axis
    gnt.set_yticklabels(['M'+ str(nb_machines-k+1) for k in range(1, nb_machines+1)])

    grid_x_ticks = np.arange(0, makespan+5, 1)
    gnt.set_xticks(grid_x_ticks)
    # Setting graph attribute
    gnt.grid(axis='x', alpha=0.3)
    #plt.grid(axis='x', color='0.95')

    for i in range(nb_machines):
        detlaT = 0
        for j in range(nb_jobs):
            gnt.broken_barh([(detlaT + idle_times[i,j], production_times[i,j])], (nb_machines-i-0.2, 0.4), facecolors =(color_names[j]), edgecolor=('black'))
            #gnt.text(detlaT + idle_times[i,j] + production_times[i,j]/2 ,2.25-i+1-0.2, str(j), ha='center', va='center')
            detlaT += idle_times[i,j] + production_times[i,j]

    plt.savefig('gantt.png')
    plt.show()

def gantt_setup(
        production_times: np.ndarray, 
        setup_times: np.ndarray, 
        idle_times: np.ndarray, 
        makespan: int
) -> None:
    """
    This function plots the gantt chart for a flowshop with setup times

     production_times: the processing times

     setup_times: the setup times

     idle_times: the idle times

     makespan: the makespan
    """
    
    nb_machines = len(production_times)
    nb_jobs = len(production_times[0])

    # Defining the canavas
    color_names = list(mcolors.CSS4_COLORS)

    plt.rcParams["figure.figsize"] = (14, 3)
    plt.rcParams["figure.dpi"] = 200
            
    # Declaring a figure "gnt"
    fig, gnt = plt.subplots()

    # Setting Y-axis limits
    gnt.set_ylim(0, nb_jobs+1)
            
    # Setting X-axis limits
    gnt.set_xlim(0, makespan+5)
            
    # Setting labels for x-axis and y-axis
    gnt.set_xlabel('Time')
    gnt.set_ylabel('Machine')
    #gnt.legend(title='alpha='+str(alpha))
            
    # Setting ticks on y-axis
    gnt.set_yticks([k for k in range(1, nb_machines+1)])
    # Labelling tickes of y-axis
    gnt.set_yticklabels(['M' + str(nb_machines-k+1) for k in range(1, nb_machines+1)])

    grid_x_ticks = np.arange(0, makespan+5, 1)
    gnt.set_xticks(grid_x_ticks)
    # Setting graph attribute
    gnt.grid(axis='x', alpha=0.3)
    plt.grid(axis='x', color='0.95')

    for i in range(nb_machines):
        detlaT = 0
        for j in range(nb_jobs):
            gnt.broken_barh([(detlaT, setup_times[i,j])], (nb_machines-i-0.2, 0.4), facecolors = (color_names[j]), edgecolor = ('blue'), hatch = '///')
            detlaT += setup_times[i,j]
            gnt.broken_barh([(detlaT + idle_times[i,j] - setup_times[i,j], production_times[i,j])], (nb_machines-i-0.2, 0.4), facecolors =(color_names[j]), edgecolor=('black'))
            #gnt.text(detlaT + idle_times[i,j] ,2-i+1-0.2, str(i), ha='center', va='center')
            detlaT += idle_times[i,j] - setup_times[i,j] + production_times[i,j]
    
    #plt.savefig('nowait.png')
    plt.show()

def ganttRobust(nom_production_times, max_deviations, scenario, idle_times, worst_makespan):
    ###
    #This function plots the gantt chart for a robust flowshop 
    #nom_production_times: the processing times
    #max_deviations: the maximum deviation values
    #scenario: the matrix representing which values deviate
    #idle_times: the idle_times times
    #makespan: the makespan
    
    nb_machines = len(nom_production_times)
    nb_jobs = len(nom_production_times[0])

    #Defining the canavas
    # Importing the matplotlib.pyplot
    #print(mcolors.TABLEAU_COLORS)
    
    color_names = list(mcolors.CSS4_COLORS)

    plt.rcParams["figure.figsize"] = (14,3)
    plt.rcParams["figure.dpi"] = 200
            
    # Declaring a figure "gnt"
    fig, gnt = plt.subplots()

    # Setting Y-axis limits
    gnt.set_ylim(0, nb_jobs+1)
            
    # Setting X-axis limits
    gnt.set_xlim(0, worst_makespan+10)
            
    # Setting labels for x-axis and y-axis
    gnt.set_xlabel('Time')
    gnt.set_ylabel('Machine')
    #gnt.legend(title='alpha='+str(alpha))
            
    # Setting ticks on y-axis
    gnt.set_yticks([k for k in range(1, nb_machines+1)])
    # Labelling tickes of y-axis
    gnt.set_yticklabels(['M'+ str(nb_machines-k+1) for k in range(1, nb_machines+1)])

    grid_x_ticks = np.arange(0, worst_makespan+10, 1)
    gnt.set_xticks(grid_x_ticks)
    # Setting graph attribute
    gnt.grid(axis='x', alpha=0.3)
    plt.grid(axis='x', color='0.95')

    for i in range(nb_machines):
        detlaT = 0
        for j in range(nb_jobs):
            if scenario[i][j]:
                gnt.broken_barh([(detlaT + idle_times[i,j], nom_production_times[i,j])], (nb_machines-i-0.2, 0.4), facecolors =(color_names[j]), edgecolor=('black'))
                #gnt.text(detlaT + idle_times[i,j] + production_times[i,j]/2 ,2.25-i+1-0.2, str(j), ha='center', va='center')
                detlaT += idle_times[i,j] + nom_production_times[i,j]
                gnt.broken_barh([(detlaT, max_deviations[i, j])], (nb_machines-i-0.2, 0.4), facecolors=(color_names[j]), edgecolor=('black'), hatch='XXX')
                detlaT += max_deviations[i, j]
            else:
                gnt.broken_barh([(detlaT + idle_times[i,j], nom_production_times[i,j])], (nb_machines-i-0.2, 0.4), facecolors =(color_names[j]), edgecolor=('black'))
                #gnt.text(detlaT + idle_times[i,j] + production_times[i,j]/2 ,2.25-i+1-0.2, str(j), ha='center', va='center')
                detlaT += idle_times[i,j] + nom_production_times[i,j]
    
    plt.axvline(worst_makespan, color='red', linestyle='-')

    plt.show()
    return(0)

def ganttRobust_setup_nowait(nom_production_times, max_deviations, scenario, 
                             setup_times, 
                             idle_times, 
                             worst_makespan):

    nb_machines = len(nom_production_times)
    nb_jobs = len(nom_production_times[0])

    # Defining the canavas
    color_names = list(mcolors.CSS4_COLORS)

    plt.rcParams["figure.figsize"] = (14, 3)
    plt.rcParams["figure.dpi"] = 200
            
    # Declaring a figure "gnt"
    fig, gnt = plt.subplots()

    # Setting Y-axis limits
    gnt.set_ylim(0, nb_jobs+1)
            
    # Setting X-axis limits
    gnt.set_xlim(0, worst_makespan+5)
            
    # Setting labels for x-axis and y-axis
    gnt.set_xlabel('Time')
    gnt.set_ylabel('Machine')
    #gnt.legend(title='alpha='+str(alpha))
            
    # Setting ticks on y-axis
    gnt.set_yticks([k for k in range(1, nb_machines+1)])
    # Labelling tickes of y-axis
    gnt.set_yticklabels(['M' + str(nb_machines-k+1) for k in range(1, nb_machines+1)])

    grid_x_ticks = np.arange(0, worst_makespan+5, 1)
    gnt.set_xticks(grid_x_ticks)
    # Setting graph attribute
    gnt.grid(axis='x', alpha=0.3)
    plt.grid(axis='x', color='0.95')    

    for i in range(nb_machines):
        detlaT = 0
        for j in range(nb_jobs):
            gnt.broken_barh([(detlaT, setup_times[i,j])], (nb_machines-i-0.2, 0.4), facecolors = (color_names[j]), edgecolor = ('blue'), hatch = '///')
            detlaT += setup_times[i,j]
            if scenario[i][j]:
                gnt.broken_barh([(detlaT + idle_times[i,j] - setup_times[i,j], nom_production_times[i,j])], (nb_machines-i-0.2, 0.4), facecolors =(color_names[j]), edgecolor=('black'))
                #gnt.text(detlaT + idle_times[i,j] + production_times[i,j]/2 ,2.25-i+1-0.2, str(j), ha='center', va='center')
                detlaT += idle_times[i,j] + nom_production_times[i,j] - setup_times[i,j]
                gnt.broken_barh([(detlaT, max_deviations[i, j])], (nb_machines-i-0.2, 0.4), facecolors=(color_names[j]), edgecolor=('black'), hatch='XXX')
                detlaT += max_deviations[i, j]
            else:
                gnt.broken_barh([(detlaT + idle_times[i,j] - setup_times[i,j], nom_production_times[i,j])], (nb_machines-i-0.2, 0.4), facecolors =(color_names[j]), edgecolor=('black'))
                #gnt.text(detlaT + idle_times[i,j] + production_times[i,j]/2 ,2.25-i+1-0.2, str(j), ha='center', va='center')
                detlaT += idle_times[i,j] + nom_production_times[i,j] - setup_times[i,j] 

    plt.show()               

def ganttCP(
        tasks: np.ndarray,
        makespan: int
):
    
    nb_jobs = len(tasks[0])
    nb_machines = len(tasks)
    color_names = list(mcolors.TABLEAU_COLORS)

    plt.rcParams["figure.figsize"] = (14,3)
    plt.rcParams["figure.dpi"] = 200
        
    # Declaring a figure "gnt"
    fig, gnt = plt.subplots()

    # Setting Y-axis limits
    gnt.set_ylim(0, nb_jobs + 1)
            
    # Setting X-axis limits
    gnt.set_xlim(0, makespan + 5)
            
    # Setting labels for x-axis and y-axis
    gnt.set_xlabel('Time')
    gnt.set_ylabel('Machine')
    #gnt.legend(title='alpha='+str(alpha))
            
    # Setting ticks on y-axis
    gnt.set_yticks([k for k in range(1, nb_machines+1)])
    # Labelling tickes of y-axis
    gnt.set_yticklabels(['M'+ str(nb_machines-k+1) for k in range(1, nb_machines+1)])

    grid_x_ticks = np.arange(0, makespan+5, 1)
    gnt.set_xticks(grid_x_ticks)
    # Setting graph attribute
    gnt.grid(axis='x', alpha=0.3)
    #plt.grid(axis='x', color='0.95')


    for i in range(nb_machines):
        for j in range(nb_jobs):
            gnt.broken_barh([(tasks[i][j].get_start(), tasks[i][j].get_size())], (nb_machines-i-0.2, 0.4), facecolors =(color_names[j]), edgecolor=('black'))
    #plt.savefig("gantt1.png")
    #plt.savefig("gantt1.pdf", format='pdf', dpi=1200)

    plt.show()

def ganttCP_setup(tasks, setup_times, makespan):
    nb_jobs = len(tasks[0])
    nb_machines = len(tasks)
    color_names = list(mcolors.TABLEAU_COLORS)

    plt.rcParams["figure.figsize"] = (14,3)
    plt.rcParams["figure.dpi"] = 200
        
    # Declaring a figure "gnt"
    fig, gnt = plt.subplots()

    # Setting Y-axis limits
    gnt.set_ylim(0, nb_jobs + 1)
            
    # Setting X-axis limits
    gnt.set_xlim(0, makespan + 5)
            
    # Setting labels for x-axis and y-axis
    gnt.set_xlabel('Time')
    gnt.set_ylabel('Machine')
    #gnt.legend(title='alpha='+str(alpha))
            
    # Setting ticks on y-axis
    gnt.set_yticks([k for k in range(1, nb_machines+1)])
    # Labelling tickes of y-axis
    gnt.set_yticklabels(['M'+ str(nb_machines-k+1) for k in range(1, nb_machines+1)])

    grid_x_ticks = np.arange(0, makespan+5, 1)
    gnt.set_xticks(grid_x_ticks)
    # Setting graph attribute
    gnt.grid(axis='x', alpha=0.3)
    #plt.grid(axis='x', color='0.95')

    end_times = [[tasks[i][j].get_end() for j in range(nb_jobs)] for i in range(nb_machines)]
    end_times = np.hstack((np.zeros((nb_machines,1)), end_times))
    for i in range(nb_machines):
        end_times[i].sort()
    setup_start = np.zeros((nb_machines, nb_jobs))
    for i in range(nb_machines):
        max_index = np.argmax(end_times[i])
        setup_start[i] = np.delete(end_times[i], max_index)

    for i in range(nb_machines):
        for j in range(nb_jobs):
            gnt.broken_barh([(tasks[i][j].get_start(), tasks[i][j].get_size())], (nb_machines-i-0.2, 0.4), facecolors =(color_names[j]), edgecolor=('black'))
            gnt.broken_barh([(setup_start[i][j], setup_times[i][j])], (nb_machines-i-0.2, 0.4), facecolors =(color_names[j]), edgecolor=('blue'), hatch=('///')) 
    plt.show()

def ganttCP_robust(tasks, dev_tasks, makespan):
    
    nb_jobs = len(tasks[0])
    nb_machines = len(tasks)
    color_names = list(mcolors.TABLEAU_COLORS)

    plt.rcParams["figure.figsize"] = (14,3)
    plt.rcParams["figure.dpi"] = 200
        
    # Declaring a figure "gnt"
    fig, gnt = plt.subplots()

    # Setting Y-axis limits
    gnt.set_ylim(0, nb_jobs + 1)
            
    # Setting X-axis limits
    gnt.set_xlim(0, makespan + 5)
            
    # Setting labels for x-axis and y-axis
    gnt.set_xlabel('Time')
    gnt.set_ylabel('Machine')
    #gnt.legend(title='alpha='+str(alpha))
            
    # Setting ticks on y-axis
    gnt.set_yticks([k for k in range(1, nb_machines+1)])
    # Labelling tickes of y-axis
    gnt.set_yticklabels(['M'+ str(nb_machines-k+1) for k in range(1, nb_machines+1)])

    grid_x_ticks = np.arange(0, makespan+5, 1)
    gnt.set_xticks(grid_x_ticks)
    # Setting graph attribute
    gnt.grid(axis='x', alpha=0.3)
    #plt.grid(axis='x', color='0.95')

    for i in range(nb_machines):
        for j in range(nb_jobs):
            gnt.broken_barh([(tasks[i][j].get_start(), tasks[i][j].get_size())], (nb_machines-i-0.2, 0.4), facecolors =(color_names[j]), edgecolor=('black'))
            if dev_tasks[i][j].get_start():
                gnt.broken_barh([(dev_tasks[i][j].get_start(), dev_tasks[i][j].get_size())], (nb_machines-i-0.2, 0.4), facecolors =(color_names[j]), edgecolor=('red'))
    #plt.savefig("gantt1.png")
    #plt.savefig("gantt1.pdf", format='pdf', dpi=1200)

    plt.show()


