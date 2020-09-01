#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <sstream>

using namespace std;

int main() {
	ifstream file; 
	file.open("../../../NYPD_Motor_Vehicle_Collisions.csv");
    if (!file.is_open()) {
        cout << "Can not open file" << endl;
        return -1;
    }

    string line;
    getline(file,line);
    stringstream ss(line);
   	string token;

    for(int i = 0; getline(ss, token, ','); i++) {
    	cout << i << " - " << token << endl;
    }

    int comma = 0;
    int comma_30 = 0;
    int comma_29 = 0;
    int comma_28 = 0;
    for(int i = 1; getline(file,line); i++) {
		for (int j = 0; (j = line.find(',', j)) != string::npos; j++) {
			comma++;
		}
	    if(comma == 29)
	    	comma_29++; 
	    else if(comma == 28)
	    	comma_28++;
	    else {
	    	comma_30++;
	    	cout << "30 commas at line " << i << endl;
	    }
	    	
	    comma = 0;
    }

    cout << "Lines with 30 commas = " << comma_30 << endl;
    cout << "Lines with 29 commas = " << comma_29 << endl;
    cout << "Lines with 28 commas = " << comma_28 << endl;

    return 0;
}


