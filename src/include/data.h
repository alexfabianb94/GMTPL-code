//: data.h

#ifndef DATA_H
#define DATA_H

#include <ilcplex/ilocplex.h>
#include <set>
#include <vector>
#include <ctime>

class Data 
{
private:
    void definirDimParam();
    bool verifyString(std::string data_string, std::string read_string);
    IloEnv ambiente;

public:
    IloInt n, deposito, numero_clusters;
    IloNum factor;
	IloNumArray xcoord;
	IloNumArray ycoord;
    std::string instance_name;
	std::vector <std::string> color;
	std::vector <std::set<IloInt>> clusters;
	IloArray <IloNumArray> c;
	IloArray <IloNumArray> t;
    Data(std::string ruta_archivo, IloEnv &env);
    void showData();
};

#endif
