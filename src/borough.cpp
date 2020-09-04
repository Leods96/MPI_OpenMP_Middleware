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

using namespace std;

/*
contains the information on the borough for each line
*/

borough_struct * getBorough(string token, unordered_map<string, borough_struct *> &borough_map) {
	if(borough_map.find(token) != borough_map.end())
		return borough_map.find(token) -> second;
	borough_struct * temp_borough = (borough_struct *) malloc (sizeof(borough_struct));
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
#ifdef _OPENMP
			int chunkSize = WEEK_ARRAY_DIM/omp_get_num_threads();
#endif
			#pragma omp parallel for schedule(dynamic, chunkSize) shared(temp_borough, b)
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
        int split_dimension = dim / omp_get_num_threads();
        int displacement = split_dimension * omp_get_thread_num();
        if(omp_get_thread_num() == omp_get_num_threads() - 1)
            split_dimension += dim - (split_dimension * omp_get_num_threads());
        borough_struct * starting_point = &f[displacement];   
        unordered_map<string, borough_struct *> borough_map_local;
        mergeBorough(starting_point, split_dimension, borough_map_local);
            //TODO remove for debug

        #pragma omp critical (debug)
        {
            cout << "Pre push_back - Sono " << omp_get_thread_num() << " map_dimension: " << borough_map_local.size() << endl;
		    set<string> bor;
		    for(auto i = borough_map_local.begin(); i != borough_map_local.end(); ++i) {
		    	bor.insert(i -> first);
		    }
		    for(auto s : bor)
		    	cout << s << endl;
		    cout << endl << endl;
        }  
        			//TILL here
        #pragma omp critical (push)
        {
            borough_map_array.push_back(borough_map_local);
        }
            //TODO remove for debug        
        #pragma omp critical (debug)
        {
            cout << "Post push_back - Sono " << omp_get_thread_num() << endl;

		    set<string> bor;
		    for(auto i : borough_map_array) {
		    	for(auto iter = i.begin(); iter != i.end(); ++iter) {
		    		bor.insert(iter -> first);
		    }
		    }
		    for(auto s : bor)
		    	cout << s << endl;
		    cout << endl << endl;
        }  
        			//TILL here
    }
    cout << "borough_map_array_size:  " << borough_map_array.size() << endl;    
    for(int distance = 1; distance < omp_get_num_threads(); distance *= 2) {
		#pragma omp parallel for 
		for(int i = 0; i < omp_get_num_threads() - distance; i += 2 * distance) {
			mergeBorough(borough_map_array[i+distance], borough_map_array[i]);
		}
	}
	            //TODO remove for debug        
            cout << "End function " << endl;

		    set<string> bor;
		    	for(auto iter = borough_map_array[0].begin(); iter != borough_map_array[0].end(); ++iter) {
		    		bor.insert(iter -> first);
		    }
		    for(auto s : bor)
		    	cout << s << endl;
		    cout << endl << endl;  
        			//TILL here
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
				int chunkSize = WEEK_ARRAY_DIM/omp_get_num_threads();
				#pragma omp parallel for schedule(dynamic, chunkSize) shared(temp, iter)
				for(int k = 0; k < WEEK_ARRAY_DIM; k++) {
					temp -> weekAccidentsCounter[k] += iter -> second -> weekAccidentsCounter[k];
					temp -> weekLethal[k] += iter -> second -> weekLethal[k];
				}
				map.insert({iter -> first, temp});
			}
		}
	}
#endif