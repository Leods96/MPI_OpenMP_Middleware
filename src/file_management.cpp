
#include <stdlib.h>
#include <sstream>
#include <string>
#include <unordered_map>
#include "date.h"
#include <set>
#include "borough.h"
#include "factor.h"

using namespace std;

void parseLine(string line, unordered_map<string, borough_struct *> &borough_map, unordered_map<string, factor_struct *> &factor_map, int weekLethalCounter[]) {
	stringstream ss(line);
	string token;
	string borough = "";
	int i = 0, deaths = 0;
	set<string> factors;
	date_struct * date;

    while(getline(ss, token, ',') && i < 23) {
    	if(token[0] == '"') {
    		while(token[token.length()-1] != '"')
    			getline(ss, token, ',');
    	} else if(token.compare("") != 0)
	    	switch(i) {
	    		case 0: //date
	    			date = computeDate(token);
	    			break;
				case 2: //borough
					borough = token;
	    			break;
				case 11:	//death
				case 13:
				case 15:
					deaths += stoi(token);
					break;
				case 17: //death + borough update
					deaths += stoi(token);
					if(deaths)							
						weekLethalCounter[date -> index] ++; 
					break;
				case 18: //factor1
				case 19: //factor2
				case 20: //factor3
				case 21: //factor4
				case 22: //factor5
					factors.insert(token);
					break;
				default:
				//do nothing
					break;
	    	}
        i++;
    }

    if(borough.compare("") != 0) {
   		borough_struct * temp_borough = getBorough(borough, borough_map);
   		temp_borough -> weekAccidentsCounter[date -> index] ++;
   		if(deaths)
   			temp_borough -> weekLethal[date -> index] ++;
    }

    factor_struct * temp_factor;
    for(auto f : factors) {
    	temp_factor = getFactor(f, factor_map);
		temp_factor -> accidentsNumber++;
		if(deaths) {
			temp_factor -> lethalAccidentsNumber++;
			temp_factor -> deathsNumber += deaths;
		}
    }
}