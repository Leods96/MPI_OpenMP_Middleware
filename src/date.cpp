#include <string>
#include <stdlib.h>
#include <sstream>
#include <stdio.h>
#include "date.h"
#include "constant.h"

date_struct * computeDate(std::string token) {
	date_struct * date = new date_struct();
    std::stringstream ss(token);
    getline(ss, token, '/');
	date -> month = stoi(token);
    getline(ss, token, '/');
	date -> day = stoi(token);
    getline(ss, token, '/');
	date -> year = stoi(token);
	int numberOfDays = daysInYear[date -> month - 1] + date -> day;	//week computed starting from 1 gen 2012, 1 week = from sunday to saturday
	if (date -> month > 1 && (date -> year % 400 == 0 || (date -> year % 4 == 0 && date -> year % 100 != 0))) { //Bisestile
		numberOfDays++;
	}
	int yeardiff = date -> year - BASE_YEAR; //difference between our year and base
	int bisestileyears = 0;
	int shift = 0;
	if(yeardiff != 0) {
		bisestileyears = (yeardiff - 1) / 4 + 1; //+ 1 because the 2012(the base) is bisestile each four years there is the bisetile and the year after there is a shift of 2 day that why diff - 1
		shift = ((yeardiff - bisestileyears) + (bisestileyears * 2)) % 7;
	}
	if(numberOfDays < 7 - shift) //if we are in the first week
		date -> index = yeardiff * YEAR_DIM;
	else { //otherwise the computation is done over the rest of the days without the first week
		numberOfDays -= 7 - shift;
		date -> index = yeardiff * YEAR_DIM + (int)(numberOfDays / 7) + 1;
	}
	if(yeardiff < 0 || date->index >= WEEK_ARRAY_DIM)
		date->index = -1;
	return date;
}