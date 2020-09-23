# MPI_OpenMP_Middleware
Processing of car accidents data. <br/>
Perform of three specific queries on a dataset written into a CSV file. <br/>
The dataset is split based on the number of processes/threads present into the system.<br/> 
Each worker will receive almost the same number of bytes to analyse.<br/> 
All the data will be collected by the root and merged.<br/>
<br/>
Both MPI and OpenMP libraries are used in order to obtain the maximum speed-up possible

# Compile

Clone -> Create a directory named Build -> Run the $make command to compile all the project<br/>
If you want to deactivate OpenMP remove the "-fopenmp" flag from the make file<br/>
Use the env var OMP_NUM_THREADS to set the desired number of omp threads<br/>
```
$ export OMP_NUM_THREADS=3
```

# Run

If you want to use MPI run the following command<br/>
```
$ mpirun -np <#ofProcesses> traffic
```
<br/>
It is possible to pass parameters by command line which represent the file's path and the set of queries to be executed: <br/>
*  Only 1 parameter = file's path.<br/>
*  3 parmeters = Only the set of queries to be executed, where each number represents the execution of the relative query (1 = is executed, 0 = is not executed).<br/>
*  Combination of the above.<br/>

Examples of the three possible execution:<br/>

```
$ mpirun -np 6 traffic ../NY_Accidents_Data
$ mpirun -np 2 traffic 1 0 1
$ mpirun -np 4 traffic ../NY_Accidents_Data 0 1 1
```

<br/>
Authors: <br/> 
- Leonardo Romano <br/> 
- Roberto Rocco

