//: solver_multi.h

#ifndef SOLVER_MULTI_M_H
#define SOLVER_MULTI_M_H

#include <ilcplex/ilocplex.h>
#include "data.h"

class SolverMulti {
private:

    IloEnv ambiente;
    IloCplex cplex;
	IloModel modelo;
    IloArray <IloBoolVarArray> x;
    IloArray <IloBoolVarArray> y;
    IloArray <IloArray <IloIntVarArray>> f;
    bool solver_success = false;
    Data datos;
    int TimeLimit = 5000;
    double start_time_det = 0;
    double start_time_cplex = 0;
    std::string model_name = "Multi-commodity";

public:
        void solve_model(std::string file_name, IloInt multiplicador, IloInt Max_Z1, IloInt Max_Z2);
        SolverMulti (IloEnv & env, Data & data);
        void showSolution();
        void graphSolution(std::string file_name);
        void saveResult(std::string file_name, IloInt multiplicador);
};

#endif

