#include "constant.h"
#include <string>
#include <unordered_map>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include "borough.h"
#include <iostream>
#include <set>

#ifdef _OPENMP
#include <omp.h>
#endif

extern int threadnum;

using namespace std;

/*
contains the information on the borough for each line
*/

borough_struct * getBorough(string token, unordered_map<string, borough_struct *> &borough_map) {
	if(borough_map.find(token) != borough_map.end())
		return borough_map.find(token) -> second;
	borough_struct * temp_borough = new borough_struct();
	borough_map.insert({token, temp_borough});
	strcpy(temp_borough -> name, token.c_str());
	return temp_borough;
}

void mergeBorough(borough_struct * b, int dim, unordered_map<string, borough_struct *> &borough_map) {
	//string token = "";
	borough_map.clear();
	for(int i = 0; i < dim; i++) {
		string token(b[i].name);
		//token = b[i].name;
		if(borough_map.find(token) == borough_map.end())
			borough_map.insert({token, &b[i]});
		else {
			borough_struct * temp_borough = borough_map.find(token) -> second;
			#pragma omp parallel for shared(temp_borough, b)
			for(int k = 0; k < WEEK_ARRAY_DIM; k++) {
				temp_borough -> weekAccidentsCounter[k] += b[i].weekAccidentsCounter[k];
				temp_borough -> weekLethal[k] += b[i].weekLethal[k];
			}
		}	
		//token = "";	
	}
}

#ifdef _OPENMP
void mergeBoroughRecursive(borough_struct * f, int dim, unordered_map<string, borough_struct *> &borough_map) {
    vector<unordered_map<string, borough_struct *>> borough_map_array;
    #pragma omp parallel
    {
        int split_dimension = dim / threadnum;
        int displacement = split_dimension * omp_get_thread_num();
        if(omp_get_thread_num() == threadnum - 1)
            split_dimension += dim - (split_dimension * threadnum);
        borough_struct * starting_point = &f[displacement];   
        unordered_map<string, borough_struct *> borough_map_local;
        mergeBorough(starting_point, split_dimension, borough_map_local);

        #pragma omp critical (push)
        {
            borough_map_array.push_back(borough_map_local);
        }

    }    
    for(int distance = 1; distance < threadnum; distance *= 2) {
		#pragma omp parallel for 
		for(int i = 0; i < threadnum - distance; i += 2 * distance) {
			mergeBorough(borough_map_array[i+distance], borough_map_array[i]);
		}
	}
	borough_map = borough_map_array[0];
}
#endif

#ifdef _OPENMP
	void mergeBorough(unordered_map<string, borough_struct *> &map_to_be_merged, unordered_map<string, borough_struct *> &map) {
		for(auto iter = map_to_be_merged.begin(); iter != map_to_be_merged.end(); ++iter) {
			if(map.find(iter -> first) == map.end())
				map.insert({iter -> first, iter -> second});
			else {
				borough_struct * temp = map.find(iter -> first) -> second;
				#pragma omp parallel for shared(temp, iter)
				for(int k = 0; k < WEEK_ARRAY_DIM; k++) {
					temp -> weekAccidentsCounter[k] += iter -> second -> weekAccidentsCounter[k];
					temp -> weekLethal[k] += iter -> second -> weekLethal[k];
				}
				map.insert({iter -> first, temp});
			}
		}
	}
#endif
