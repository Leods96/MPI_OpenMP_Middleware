#ifndef FACTOR_H
#define FACTOR_H

#include "constant.h"
#include <string>
#include <unordered_map>

typedef struct {
	char name[NAME_DIM];
	int accidentsNumber;
	int lethalAccidentsNumber;
	int deathsNumber;
} factor_struct;

factor_struct * getFactor(std::string token, std::unordered_map<std::string, factor_struct *> &factor_map);

void mergeFactor(factor_struct * f, int dim, std::unordered_map<std::string, factor_struct *> &factor_map);

#ifdef _OPENMP
void mergeFactor(std::unordered_map<std::string, factor_struct *> &map_to_be_merged, std::unordered_map<std::string, factor_struct *> &map);

void mergeFactorRecursive(factor_struct * f, int dim, std::unordered_map<std::string, factor_struct *> &factor_map);

#endif



#endif