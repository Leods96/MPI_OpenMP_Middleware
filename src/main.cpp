#include <cstring>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <set>
#include <chrono>

#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "constant.h"
#include "file_management.h"
#include "date.h"
#include "borough.h"
#include "factor.h"

using namespace std;


/*
contains the information on the contributing factor for each line
*/

#ifdef _OPENMP
	void mergeWeekLethalArray(int array[], vector<int> &vec){
		int chunkSize = WEEK_ARRAY_DIM/omp_get_num_threads();
		#pragma omp parallel for schedule(dynamic, chunkSize) shared(array)
		for(int i = 0; i < WEEK_ARRAY_DIM; i++) {
			array[i] += vec[i];
		}
	}
#endif

template <typename T>
void getArrayFromMap(unordered_map<string, T *> &map, T array[]) {
	int i = 0;
	for(auto iter = map.begin(); iter != map.end(); ++iter) {
		memcpy(&array[i], iter -> second, sizeof(T));
		i++;
	}
} 

void printFirstQuery(int result[]) {
	cout << "----------------------------First query----------------------------" << endl << endl;
	for(int i = 0; i < WEEK_ARRAY_DIM; i++) {
		if(i % YEAR_DIM == 0)
			cout << 2012 + (i / YEAR_DIM) << " :" << endl;
		if(result[i] != 0)
			cout << "\t - week " << (i % YEAR_DIM) + 1 << " : " << result[i] << " death" << endl;
	}
}

void printSecondQuery(unordered_map<string, factor_struct *> &factor_map) {
	cout << "----------------------------Second query----------------------------" << endl << endl;
	for(auto iter = factor_map.begin(); iter != factor_map.end(); ++iter) {
		auto p = iter -> second;
		cout << endl << iter -> first << " - #Accidents = " << p -> accidentsNumber << " - Average = " 
		<< (double)p -> lethalAccidentsNumber / (double)p -> accidentsNumber * 100 << endl;
	}
}

void printThirdQuery(unordered_map<string, borough_struct *> &borough_map) {
	cout << "----------------------------Third query----------------------------" << endl << endl;
	for(auto iter = borough_map.begin(); iter != borough_map.end(); ++iter) {
		auto p = iter -> second;
		string name = iter -> first;
		cout << "\nBourgh: " << name;
		for(int i = 0; i < WEEK_ARRAY_DIM; i++) {
			if(i % YEAR_DIM == 0)
				cout << endl << 2012 + (i / YEAR_DIM) << " :" << endl;
			if(p -> weekAccidentsCounter[i] != 0 && p -> weekLethal[i] != 0) {
				cout << "\t- week " << (i % YEAR_DIM) + 1 << " - #Accidents: " << p -> weekAccidentsCounter[i] << " - AvgLethal: ";
				if(p -> weekAccidentsCounter[i] != 0)
					cout << (double)p -> weekLethal[i]/(double)p -> weekAccidentsCounter[i] * 100 << endl;
				else
					cout << "0" << endl;
			}
		}
	}
}

int main(int argc, char * argv[]) {
	/*	
	each year is composed at most by 54 week, each week start from sun to sat, 
	if one year start from saturday the first week is composed by only one day and the first week is burned 
	*/ 	
	int rank, size, i = 0, filesize, starting_offset, limit;
    streampos begin, end;
	string line;

	/*
	Setup for datatype creation, used during the gatherv to pass complex data
	*/
	MPI_Datatype boroughtype, factortype;
	MPI_Datatype typef[4] = {MPI_CHAR, MPI_INT, MPI_INT, MPI_INT}, typeb[3] = {MPI_CHAR, MPI_INT, MPI_INT};
	int blocklenf[4] = {NAME_DIM, 1, 1, 1}, blocklenb[3] = {NAME_DIM, WEEK_ARRAY_DIM, WEEK_ARRAY_DIM};
	MPI_Aint dispf[4], dispb[3];
	MPI_Status stat;

    factor_struct * f;
    borough_struct * b;

    dispf[0] = reinterpret_cast<std::uintptr_t>(f -> name) - reinterpret_cast<std::uintptr_t>(f);
	dispf[1] = reinterpret_cast<std::uintptr_t>(&f -> accidentsNumber) - reinterpret_cast<std::uintptr_t>(f);
    dispf[2] = reinterpret_cast<std::uintptr_t>(&f -> lethalAccidentsNumber) - reinterpret_cast<std::uintptr_t>(f);
    dispf[3] = reinterpret_cast<std::uintptr_t>(&f -> deathsNumber) - reinterpret_cast<std::uintptr_t>(f);

    dispb[0] = reinterpret_cast<std::uintptr_t>(b -> name) - reinterpret_cast<std::uintptr_t>(b);
    dispb[1] = reinterpret_cast<std::uintptr_t>(&b -> weekAccidentsCounter) - reinterpret_cast<std::uintptr_t>(b);
    dispb[2] = reinterpret_cast<std::uintptr_t>(&b -> weekLethal) - reinterpret_cast<std::uintptr_t>(b);

#ifdef _OPENMP
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
#else
    MPI_Init(&argc, &argv); //MPI Initialization
#endif
    /*
    new datatype creation and commit 
    */
    MPI_Type_create_struct(4, blocklenf, dispf, typef, &factortype);
    MPI_Type_commit(&factortype);
    MPI_Type_create_struct(3, blocklenb, dispb, typeb, &boroughtype);
    MPI_Type_commit(&boroughtype);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int root = 0; 

   	ifstream file; 
   	file.open("../../NYPD_Motor_Vehicle_Collisions.csv");
    if (!file.is_open()) {
        cout << "Can not open file" << endl;
        return -1;
    }
	auto start = chrono::high_resolution_clock::now();

	/*
	Calculation of the file size in order to spli the data between the processes
	*/
    begin = file.tellg();
    file.seekg(0, ios::end);
    end = file.tellg();
    filesize = end-begin;
    starting_offset = filesize/size * rank;
    limit = starting_offset + filesize/size;
#ifdef _OPENMP   
	file.close();
#else
	file.seekg(starting_offset);
#endif
	/*
	Parse the file until it is finished or the portion's limit is reached
    */
#ifdef _OPENMP
	vector< unordered_map<string, factor_struct *> > factor_map_vector;
    vector< unordered_map<string, borough_struct *> > borough_map_vector;
    vector< vector<int> > week_lethal_vector;
#else
	unordered_map<string, factor_struct *> factor_map;
	unordered_map<string, borough_struct *> borough_map;
	vector<int> weekLethalCounter;
#endif

   	#pragma omp parallel shared (limit, starting_offset, borough_map_vector, factor_map_vector, week_lethal_vector) private (file, line)
   	{	
		unordered_map<string, factor_struct *> factor_map_local; //Second query structure
    	unordered_map<string, borough_struct *> borough_map_local; //Third query structure
		vector<int> weekLethalCounter_local(WEEK_ARRAY_DIM,0);
   		int private_limit = 0;
#ifdef _OPENMP
   	   	file.open("../../NYPD_Motor_Vehicle_Collisions.csv");
   		int chunkSize = (limit - starting_offset) / omp_get_num_threads();
   		int private_offset = chunkSize  * omp_get_thread_num();
   		private_limit = private_offset + chunkSize;
   		file.seekg(starting_offset + private_offset);
   		#pragma omp critical (b)
   		{
   			borough_map_vector.push_back(borough_map_local);
   			factor_map_vector.push_back(factor_map_local);
   			week_lethal_vector.push_back(weekLethalCounter_local);
   		}
   		#pragma omp barrier
#endif
    	getline(file,line); //Positioning on the beginning of the next line, 0 skip the header line
	    while(getline(file,line)) {
#ifdef _OPENMP
			parseLine(line, borough_map_vector[omp_get_thread_num()], factor_map_vector[omp_get_thread_num()], &(week_lethal_vector[omp_get_thread_num()][0]));
#else 
	    	parseLine(line, borough_map_local, factor_map_local, &(weekLethalCounter_local[0]));
#endif
	    	if(file.tellg() > limit + private_limit) //if we reach the limit -> exit
	    		break;
	    }
	    file.close();
#ifdef _OPENMP
#else
	factor_map = factor_map_local;
	borough_map = borough_map_local;
	weekLethalCounter = weekLethalCounter_local;
#endif
	}
#ifdef _OPENMP
	for(int distance = 1; distance < omp_get_num_threads(); distance *= 2) {
		#pragma omp parallel for 
		for(int i = 0; i < omp_get_num_threads() - distance; i += 2 * distance) {
			mergeWeekLethalArray(&(week_lethal_vector[i][0]), week_lethal_vector[i+distance]);
			mergeFactor(factor_map_vector[i+distance], factor_map_vector[i]);
			mergeBorough(borough_map_vector[i+distance], borough_map_vector[i]);
		}
	}
	unordered_map<string, factor_struct *> factor_map = factor_map_vector[0];
	unordered_map<string, borough_struct *> borough_map = borough_map_vector[0];
	vector<int> weekLethalCounter = week_lethal_vector[0];

#endif
	
	/*
	Gather for first query
    */
    int first_query_res_buf[WEEK_ARRAY_DIM];
    MPI_Reduce(&(weekLethalCounter[0]), first_query_res_buf, WEEK_ARRAY_DIM, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
	/*
	Gatherv for second query
   	*/
   	int * dimension_array = (rank == root) ? (int *) malloc(sizeof(int) * size) : NULL; 
	int gather_size = factor_map.size(), displ_array[size], recv_buf_dim;
    factor_struct factors_for_gatherv[gather_size];
    getArrayFromMap(factor_map, factors_for_gatherv);

    MPI_Gather(&gather_size, 1, MPI_INT, dimension_array, 1, MPI_INT, root, MPI_COMM_WORLD);
    recv_buf_dim = gather_size;
    if(rank == root) {
    	displ_array[0] = 0;
    	for (int i = 1; i < size; i ++) {
    		displ_array[i] = dimension_array[i - 1] + displ_array[i - 1];
    		recv_buf_dim += dimension_array[i];
    	}
    }
    
    factor_struct * factor_recv_buf = (rank == root) ? (factor_struct *) malloc (sizeof(factor_struct) * recv_buf_dim) : NULL;
    MPI_Gatherv(factors_for_gatherv, gather_size, factortype, factor_recv_buf, dimension_array, displ_array, factortype, root, MPI_COMM_WORLD);
   
    if(rank == root)
#ifdef _OPENMP
		mergeFactorRecursive(factor_recv_buf, recv_buf_dim, factor_map);
#else
    	mergeFactor(factor_recv_buf, recv_buf_dim, factor_map);
#endif
	/*
	Gatherv for third query
    */
    gather_size = borough_map.size();
    borough_struct boroughs_for_gatherv[gather_size];
    getArrayFromMap(borough_map, boroughs_for_gatherv);

    MPI_Gather(&gather_size, 1, MPI_INT, dimension_array, 1, MPI_INT, root, MPI_COMM_WORLD);
	recv_buf_dim = gather_size;
	if(rank == root){
		displ_array[0] = 0;
		for(int i = 1; i < size; i++) {
    		displ_array[i] = dimension_array[i - 1] + displ_array[i - 1];
    		recv_buf_dim += dimension_array[i];
		}
	}

	borough_struct * borough_recv_buf = NULL;
	if(rank == root)
		borough_recv_buf = (borough_struct *) malloc (sizeof(borough_struct) * recv_buf_dim);
    MPI_Gatherv(boroughs_for_gatherv, gather_size, boroughtype, borough_recv_buf, dimension_array, displ_array, boroughtype, root, MPI_COMM_WORLD);

	if(rank == root)
#ifdef _OPENMP
    	mergeBoroughRecursive(borough_recv_buf, recv_buf_dim, borough_map);
#else
		mergeBorough(borough_recv_buf, recv_buf_dim, borough_map);
#endif

	/*
	Print of the results
	*/
	if(rank == 0) {
	    auto stop = chrono::high_resolution_clock::now();

	    //printFirstQuery(first_query_res_buf);
	    //printSecondQuery(factor_map);
	    //printThirdQuery(borough_map);

	    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	    cout << "Execution time = " << duration.count() << " milliseconds" << endl;
	}
	MPI_Finalize();
	return 0;
}
