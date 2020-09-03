#include <string>
#include <unordered_map>
#include <string.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif


#include "factor.h"

using namespace std;

factor_struct * getFactor(string token, unordered_map<string, factor_struct *> &factor_map) {
	if(factor_map.find(token) != factor_map.end())
		return factor_map.find(token) -> second;
	factor_struct * temp_factor = (factor_struct *) malloc (sizeof(factor_struct));
	factor_map.insert({token, temp_factor});
	strcpy(temp_factor -> name, token.c_str());
	temp_factor -> accidentsNumber = 0;
	temp_factor -> lethalAccidentsNumber = 0;
	temp_factor -> deathsNumber = 0;
	return temp_factor;
}

void mergeFactor(factor_struct * f, int dim, unordered_map<string, factor_struct *> &factor_map) {
    //TODO skippare la porzione della root e diminuire le iterazioni
	//string token = "";
	factor_map.clear();
    for(int distance = 1; distance < dim; distance *= 2)
        for(int i = 0; i < dim - distance ; i += distance * 2) {
            string token(f[i].name);
            //for (int j = 0; j < NAME_DIM && f[i].name[j] != '\0'; j++) {
            //	token = token + f[i].name[j];
            //}
            if(factor_map.find(token) == factor_map.end())
                factor_map.insert({token, &f[i]});
            else {
                factor_struct * temp_factor = factor_map.find(token) -> second;
                temp_factor -> accidentsNumber += f[i].accidentsNumber;
                temp_factor -> lethalAccidentsNumber += f[i].lethalAccidentsNumber;
                temp_factor -> deathsNumber += f[i].deathsNumber;
            }
            //token = "";
        }
}

#ifdef _OPENMP
void mergeFactorRecursive(factor_struct * f, int dim, unordered_map<string, factor_struct *> &factor_map) {
    vector<unordered_map<string, factor_struct *>> factor_map_array;
    #pragma omp parallel
    {
        int split_dimension = dim / omp_get_num_threads();
        if(omp_get_thread_num() == omp_get_num_threads() - 1)
            split_dimension += dim - (split_dimension * omp_get_num_threads());
        factor_struct * starting_point = &(f[(split_dimension) * omp_get_thread_num()]);
        unordered_map<string, factor_struct *> factor_map_local;
        mergeFactor(starting_point, split_dimension, factor_map_local);
        #pragma omp critical
        {
            factor_map_array.push_back(factor_map_local);
        }
    }
    for(int distance = 1; distance < omp_get_num_threads(); distance *= 2) {
		#pragma omp parallel for 
		for(int i = 0; i < omp_get_num_threads() - distance; i += 2 * distance) {
			mergeFactor(factor_map_array[i+distance], factor_map_array[i]);
		}
	}
	factor_map = factor_map_array[0];
}
#endif

#ifdef _OPENMP
	void mergeFactor(unordered_map<string, factor_struct *> &map_to_be_merged, unordered_map<string, factor_struct *> &map) {
		for (auto iter = map_to_be_merged.begin(); iter != map_to_be_merged.end(); ++iter) {
			if(map.find(iter -> first) == map.end())
				map.insert({iter -> first, iter -> second});
			else {
				factor_struct * temp = map.find(iter -> first) -> second;
				temp -> accidentsNumber += iter -> second -> accidentsNumber;
				temp -> lethalAccidentsNumber += iter -> second -> lethalAccidentsNumber;
				temp -> deathsNumber += temp -> deathsNumber;
			}
		}
	}
#endif