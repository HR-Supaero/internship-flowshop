# Flowshop research
This is my codes for my gap year's first project

## Problem presentation
  The problem I studied was a variant of the Flowshop problem.
  The standard formulation of the problem goes as follows:
  
  Say we have a fixed sequence of $m$ machines $(M_{1}, \ldots, M_{m})$, and $n$ jobs that have to be processed on the sequence of machines in some order.

  All flowshop variations follows at least those three rules:
  - Two jobs cannot be processed on the same time on the same machine
  - The order in which the jobs are processed is the same on each machine
  - Before being processed on a machine, a job processing has to be completed on the preceeding machine in the sequence.

  The optimization consists of **finding the job sequence $(j_{1}, \ldots, j_{k}, \ldots, j_{n})$ that minimizes the completion time of the last job on the last machine**, otherwise known as the **makespan**.

  There were three main variations of flowshop that I worked on:
  1. The Standard Flowshop
  2. The No-wait Flowshop, a job is being processed without waiting time between each machine
  3. The No-wait with sequence dependent setup time Flowshop, where before the processing of a job, a setup time that depends on the preceeding job is considered

  All variations also feature a robust counterpart, where processing times can randomly deviate from their nominal value. One a way in which the processing times can deviate is called a **scenario**, the goal is now to find the job sequence that minimizes the makespan in the worst case scenario.

## Repository organisation



## Code execution
