#ifndef DATE_H
#define DATE_H

#define YEAR_DIM 54
const int daysInYear[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};

/*
contains the information on the date for each line
*/
typedef struct {
	int month;
	int day;
	int year;
	int index;
} date_struct;

date_struct * computeDate(std::string token);

#endif