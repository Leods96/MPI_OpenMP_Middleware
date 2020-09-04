#include <cstring>
#include <iostream>
#include <fstream>

using namespace std;



int main() {
	string inputfilename = "../../../NYPD_Motor_Vehicle_Collisions.csv"; 
	string outputfilename = "cloned_file.csv";
	string line;
	ofstream outputFile;
	ofstream fs;
	ifstream inputfile; 
	int multiplier = 10;

	outputFile.open(outputfilename);

   	for (int i = 0; i < multiplier; i++) {
   		inputfile.open(inputfilename);
   		getline(inputfile, line);
   		if(i == 0)
	   		outputFile << line << endl;  	
	   	while(getline(inputfile, line)) 
	   		outputFile << line << endl;  
	   	inputfile.close();
    }
   	
   	outputFile.close();
   	inputfile.close();   	

}