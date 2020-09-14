#ifndef FILE_MANAGEMENT_H
#define FILE_MANAGEMENT_H

#include <string>
#include <unordered_map>
#include "factor.h"
#include "borough.h"

int parseLine(std::string line, std::unordered_map<std::string, borough_struct *> &borough_map, std::unordered_map<std::string, factor_struct *> &factor_map, int weekLethalCounter[]);

#endif