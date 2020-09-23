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

# Run

If you want to use MPI run the following command<br/>
$mpirun -np <#ofProcesses> traffic
<br/>
It is possible to pass parameters by command line to pass the file and the set of queries to be executed:
1. Only 1 parameter = file name, eg: $mpirun -np 2 traffic ../NY_Accidents_Data
2. Only the set of queries to be executed, eg: $mpirun -np 2 traffic 1 0 1, where each number represents the execution of the relative query (1 = is executed, 0 = is not executed)
3. Combination of the above, eg: $mpirun -np 2 traffic ../NY_Accidents_Data 0 1 1

<br/>
Author: Leonardo Romano<br/>

