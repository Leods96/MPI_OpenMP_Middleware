#ifndef BOROUGH_H
#define BOROUGH_H

#include "constant.h"
#include <string>
#include <unordered_map>

typedef struct {
	char name[NAME_DIM];
	int weekAccidentsCounter[WEEK_ARRAY_DIM];
	int weekLethal[WEEK_ARRAY_DIM];
} borough_struct;

borough_struct * getBorough(std::string token, std::unordered_map<std::string, borough_struct *> &borough_map);

void mergeBorough(borough_struct * b, int dim, std::unordered_map<std::string, borough_struct *> &borough_map);

#ifdef _OPENMP
void mergeBorough(std::unordered_map<std::string, borough_struct *> &map_to_be_merged, std::unordered_map<std::string, borough_struct *> &map);
#endif

#endif