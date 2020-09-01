#include <cstring>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <sstream>
#include <unordered_map>
#include <chrono>
#include <mpi.h>

#define WEEK_ARRAY_DIM 324 
#define YEAR_DIM 54
#define NAME_DIM 50

using namespace std;

const int daysInYear[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
/*	
each year is composed at most by 54 week, each week start from sun to sat, 
if one year start from saturday the first week is composed by only one day and the first week is burned 
*/
int weekLethalCounter[WEEK_ARRAY_DIM]; 

/*
contains the information on the date for each line
*/
typedef struct {
	int month;
	int day;
	int year;
	int index;
} date_struct;
date_struct* date = (date_struct*) malloc (sizeof(date_struct));

/*
contains the information on the contributing factor for each line
*/
typedef struct {
	char name[NAME_DIM];
	int accidentsNumber;
	int lethalAccidentsNumber;
	int deathsNumber;
} factor_struct;
unordered_map<string, factor_struct *> factor_map; //Second query structure

/*
contains the information on the borough for each line
*/
typedef struct {
	char name[NAME_DIM];
	int weekAccidentsCounter[WEEK_ARRAY_DIM];
	int weekLethal[WEEK_ARRAY_DIM];
} borough_struct;
unordered_map<string, borough_struct *> borough_map; //Third query structure

void computeDate(string token) {
    stringstream ss(token);
    getline(ss, token, '/');
	date -> month = stoi(token);
    getline(ss, token, '/');
	date -> day = stoi(token);
    getline(ss, token, '/');
	date -> year = stoi(token);
	//week computed starting from 1 gen 2012, 1 week = from sunday to saturday
	int numberOfDays = daysInYear[date -> month - 1] + date -> day;
	if (date -> month > 1 && (date -> year % 400 == 0 || (date -> year % 4 == 0 && date -> year % 100 != 0))) { //Bisestile
		numberOfDays++;
	}
	int yeardiff = date -> year - 2012; //difference between our year and 2012 that is the base
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
}

borough_struct * getBorough(string token) {
	if(borough_map.find(token) != borough_map.end())
		return borough_map.find(token) -> second;
	borough_struct * temp_borough = (borough_struct *) malloc (sizeof(borough_struct));
	borough_map.insert({token, temp_borough});
	strcpy(temp_borough -> name, token.c_str());
	for(int i = 0; i < WEEK_ARRAY_DIM; i++) {
		temp_borough -> weekAccidentsCounter[i] = 0;
		temp_borough -> weekLethal[i] = 0;
	}
	return temp_borough;
}

factor_struct * getFactor(string token) {
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

void mergeFactor(factor_struct * f, int dim) {
	cout << "start merge" << endl;
	//TODO skippare la porzione della root e diminuire le iterazioni
	string token = "";
	factor_map.clear();
		cout << "clean" << endl;

	for(int i = 0; i < dim ; i++) {
		//token = f[i].name;
		for (int j = 0; j < NAME_DIM && f[i].name[j] != '\0'; j++) {
			token = token + f[i].name[j];
		}
		cout << "Name = " << token << endl;

		if(factor_map.find(token) == factor_map.end())
			factor_map.insert({token, &f[i]});
		else {
			factor_struct * temp_factor = factor_map.find(token) -> second;
			temp_factor -> accidentsNumber += f[i].accidentsNumber;
			temp_factor -> lethalAccidentsNumber += f[i].lethalAccidentsNumber;
			temp_factor -> deathsNumber += f[i].deathsNumber;
			factor_map.insert({token, temp_factor});
		}
		token = "";
	}
	cout << "End merge" << endl;

}

void mergeBorough(borough_struct * b, int dim) {
	cout << "Start merge" << endl;
	string token = "";
	borough_map.clear();
	for(int i = 0; i < dim; i++) {
		token = b[i].name;
		if(borough_map.find(token) == borough_map.end())
			borough_map.insert({token, &b[i]});
		else {
			borough_struct * temp_borough = borough_map.find(token) -> second;
			for(int k = 0; k < WEEK_ARRAY_DIM; k++) {
				temp_borough -> weekAccidentsCounter[k] += b[i].weekAccidentsCounter[k];
				temp_borough -> weekLethal[k] += b[i].weekLethal[k];
			}
		}	
		token = "";	
	}
}

template <typename T>
void getArrayFromMap(unordered_map<string, T *> map, T array[]) {
	int i = 0;
	for(auto iter = map.begin(); iter != map.end(); ++iter) {
		/*strcpy(array[i].name, iter -> second -> name);
		array[i].accidentsNumber = iter -> second -> accidentsNumber;
		array[i].lethalAccidentsNumber = iter -> second -> lethalAccidentsNumber;
		array[i].deathsNumber = iter -> second -> deathsNumber;*/
		memcpy(&array[i], iter -> second, sizeof(T));
		i++;
	}
} 

void parseLine(string line) {
	stringstream ss(line);
	string token;
	int i = 0;
	factor_struct * temp_factor[5];
	borough_struct * temp_borough = NULL;
	int deaths = 0;

	//TODO if the same accident have more than 1 contrib fact the data are multiplied -> check if it is ok
	//it is possible to put a bool into the struct and update it only of it is false...at the and set to false
    while(getline(ss, token, ',') && i < 23) {
    	if(token[0] == '"') {
    		while(token[token.length()-1] != '"')
    			getline(ss, token, ',');
    	} else if(token.compare("") != 0)
	    	switch(i) {
	    		case 0: //date
	    			computeDate(token);
	    			//cout << "date: " << token << " - day: " << date -> day << " month: " << date -> month << " year: " << date -> year << " index: " << date -> index << endl;
	    			break;
				case 2: //borough
					temp_borough = getBorough(token);
	    			break;
				case 11:
				case 13:
				case 15:	//death
					deaths += stoi(token);
					break;
				case 17: //death + borough update
					deaths += stoi(token);
					if (temp_borough != NULL) 
						temp_borough -> weekAccidentsCounter[date -> index] ++;
					if(deaths) {
						if(temp_borough != NULL)
							temp_borough -> weekLethal[date -> index] ++;
						weekLethalCounter[date -> index] ++; 
					}
					break;
				case 18: //factor1
					temp_factor[0] = getFactor(token);
					temp_factor[0] -> accidentsNumber ++;
					if(deaths) {
						temp_factor[0] -> lethalAccidentsNumber++;
						temp_factor[0] -> deathsNumber += deaths;
					}
					break;
				case 19: //factor2
					temp_factor[1] = getFactor(token);
					temp_factor[1] -> accidentsNumber ++;
					if(deaths) {
						temp_factor[1] -> lethalAccidentsNumber++;
						temp_factor[1] -> deathsNumber += deaths;
					}
					break;
				case 20: //factor3
					temp_factor[2] = getFactor(token);
					temp_factor[2] -> accidentsNumber ++;
					if(deaths) {
						temp_factor[2] -> lethalAccidentsNumber++;
						temp_factor[2] -> deathsNumber += deaths;
					}
					break;
				case 21: //factor4
					temp_factor[3] = getFactor(token);
					temp_factor[3] -> accidentsNumber ++;
					if(deaths) {
						temp_factor[3] -> lethalAccidentsNumber++;
						temp_factor[3] -> deathsNumber += deaths;
					}
					break;
				case 22: //factor5
					temp_factor[4] = getFactor(token);
					temp_factor[4] -> accidentsNumber ++;
					if(deaths) {
						temp_factor[4] -> lethalAccidentsNumber++;
						temp_factor[4] -> deathsNumber += deaths;
					}
					break;
				default:
				//do nothing
					break;
	    	}
	    	i++;
    }
}

void printFirstQuery(int result[WEEK_ARRAY_DIM]) {
	cout << "----------------------------First query----------------------------" << endl << endl;
	for(int i = 0; i < WEEK_ARRAY_DIM; i++) {
		if(i % YEAR_DIM == 0)
			cout << 2012 + (i / YEAR_DIM) << " :" << endl;
		if(result[i] != 0)
			cout << "\t - week " << (i % YEAR_DIM) + 1 << " : " << result[i] << "death" << endl;
	}
}

void printSecondQuery() {
	cout << "----------------------------Second query----------------------------" << endl << endl;
	for(auto iter = factor_map.begin(); iter != factor_map.end(); ++iter) {
		auto p = iter -> second;
		cout << endl << iter -> first << " - #Accidents = " << p -> accidentsNumber << " - Average = " 
		<< (double)p -> lethalAccidentsNumber / (double)p -> accidentsNumber * 100 << endl;
	}
}

void printThirdQuery() {
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
	int rank, size, i = 0, filesize, starting_offset, limit;
    streampos begin, end;
	string line;

	MPI_Datatype boroughtype, factortype;
	MPI_Datatype typef[4] = {MPI_CHAR, MPI_INT, MPI_INT, MPI_INT}, typeb[3] = {MPI_CHAR, MPI_INT, MPI_INT};
	int blocklenf[4] = {NAME_DIM, 1, 1, 1}, blocklenb[3] = {NAME_DIM, WEEK_ARRAY_DIM, WEEK_ARRAY_DIM};
	MPI_Aint dispf[4], dispb[3];
	MPI_Status stat;

    factor_struct * f = (factor_struct *) malloc (sizeof(factor_struct));
    borough_struct * b = (borough_struct *) malloc (sizeof(borough_struct));

    dispf[0] = reinterpret_cast<std::uintptr_t>(f -> name) - reinterpret_cast<std::uintptr_t>(f);
	dispf[1] = reinterpret_cast<std::uintptr_t>(&f -> accidentsNumber) - reinterpret_cast<std::uintptr_t>(f);
    dispf[2] = reinterpret_cast<std::uintptr_t>(&f -> lethalAccidentsNumber) - reinterpret_cast<std::uintptr_t>(f);
    dispf[3] = reinterpret_cast<std::uintptr_t>(&f -> deathsNumber) - reinterpret_cast<std::uintptr_t>(f);

    dispb[0] = reinterpret_cast<std::uintptr_t>(b -> name) - reinterpret_cast<std::uintptr_t>(b);
    dispb[1] = reinterpret_cast<std::uintptr_t>(&b -> weekAccidentsCounter) - reinterpret_cast<std::uintptr_t>(b);
    dispb[2] = reinterpret_cast<std::uintptr_t>(&b -> weekLethal) - reinterpret_cast<std::uintptr_t>(b);

    MPI_Init(&argc, &argv);

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

    begin = file.tellg();
    file.seekg(0, ios::end);
    end = file.tellg();
    filesize = end-begin;

    starting_offset = filesize/size * rank;
    limit = starting_offset + filesize/size;

    file.seekg(starting_offset);

	getline(file,line); //positioning on the beginnin of the next line, 0 skip the header

    while(getline(file,line)) { //parse the file since it is finished 
    	parseLine(line);
    	if(file.tellg() > limit) //if we reach the limit -> exit
    		break;
    }


	file.close();
//Gather for first query
    int first_query_res_buf[WEEK_ARRAY_DIM];
    MPI_Reduce(weekLethalCounter, first_query_res_buf, WEEK_ARRAY_DIM, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
//Gatherv for second query
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
    	mergeFactor(factor_recv_buf, recv_buf_dim);
//Gatherv for third query
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
    	mergeBorough(borough_recv_buf, recv_buf_dim);
//Print of the results
	MPI_Finalize();
	if(rank == 0) {
	    auto stop = chrono::high_resolution_clock::now();

	    printFirstQuery(first_query_res_buf);
	    printSecondQuery();
	    printThirdQuery();

	    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	    cout << "Execution time = " << duration.count() << " milliseconds" << endl;
	}
	return 0;
}
