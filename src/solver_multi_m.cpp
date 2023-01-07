//: solver.cpp

#include <ilcplex/ilocplex.h>
#include <set>
#include <vector>
#include <ctime>
#include "include/data.h"
#include "include/solver_multi_m.h"


SolverMMulti::SolverMMulti(IloEnv & env, Data & data): 
		cplex(env), modelo(env), x(env), y (env), f(env), datos(data), ambiente(env){
}

void SolverMMulti::solve_model(std::string file_name, IloInt multiplicador, IloInt Max_Z1, IloInt Max_Z2)
{
	//Inicialización del las clases del software

  
    x = IloArray <IloBoolVarArray> (ambiente, datos.n);
    y = IloArray <IloBoolVarArray> (ambiente, datos.n);
	f = IloArray <IloArray <IloIntVarArray>> (ambiente, datos.n);
	for (int i = 0; i < datos.n; i++)
	{
		x[i] = IloBoolVarArray(ambiente, datos.n);
		y[i] = IloBoolVarArray(ambiente, datos.n);
		f[i] = IloArray <IloIntVarArray>(ambiente, datos.n);
		for (int j = 0; j < datos.n; j++)
		{
			f[i][j] = IloIntVarArray(ambiente, datos.n, 0, datos.n);
		}
	}

	

	IloExpr Z1(ambiente);
	for (int i = 0; i < datos.n; i++)
	{
		for (int j = 0; j < datos.n; j++)
		{
			if (i != j){
				Z1 += datos.c[i][j] * x[i][j];
			}
		}
	}
	

    IloExpr Z2(ambiente);
	for (int i = 0; i < datos.n; i++)
	{
		for (int j = 0; j < datos.n; j++)
		{
			
			if (i != j){
				Z2 += 3*datos.t[i][j] * y[i][j];
			}
		}
	}
	for (int i = 0; i < datos.n; i++)
	{
		for (int j = 0; j < datos.n; j++)
		{
			for (int k = 0; k < datos.n; k++)
			{
				if (i != j && i != k){
					Z2 += datos.t[i][j] * f[i][j][k];
				}
			}
		}
	}

	for (int p = 0; p < datos.numero_clusters; p++)
	{
		for (std::set<IloInt>::iterator k = datos.clusters[p].begin(); k != datos.clusters[p].end(); k++){
			for (int j = 0; j < datos.n; j++)
			{
				if (datos.clusters[p].find(j) == datos.clusters[p].end() && *k != j){
					Z2 += datos.t[*k][j] * f[*k][j][*k];
				}
			}
		}
	}

	for (int p = 0; p < datos.numero_clusters; p++)
	{
		for (std::set<IloInt>::iterator i = datos.clusters[p].begin(); i != datos.clusters[p].end(); i++)
		{
			for (std::set<IloInt>::iterator j = datos.clusters[p].begin(); j != datos.clusters[p].end(); j++)
			{

				if (*i != *j){
					Z2 += datos.t[*i][*j] * x[*i][*j];}
			}
		}
	}





	//modelo.add(IloMinimize(ambiente, mult Z1 + Z2));
	modelo.add(IloMinimize(ambiente, multiplicador*Z1 + (1-multiplicador)*Z2));
    
    if (Max_Z1 != 10000000000)
		modelo.add(Z1 <= Max_Z1);
	if (Max_Z2 != 10000000000)
    	modelo.add(Z2 <= Max_Z2);
    
	Z1.end();
    Z2.end();
	//ambiente.out()  << endl << "Funcion objetivo ingresada correctamente" << endl << endl;

	for (int p = 0; p < datos.numero_clusters; p++)
	{
		for (std::set<IloInt>::iterator i = datos.clusters[p].begin(); i != datos.clusters[p].end(); i++){
			IloExpr suma(ambiente);
			for (std::set<IloInt>::iterator j = datos.clusters[p].begin(); j != datos.clusters[p].end(); j++)
			{
				suma += y[*i][*j];
			}
			modelo.add(suma == 1);
			suma.end();


		}
	}

	//ambiente.out()  << "Restriccion 1 ingresada correctamente" << endl;

	for (int j = 0; j < datos.n; j++)
	{
		IloExpr suma(ambiente);
		for (int i = 0; i < datos.n; i++)
		{
			if (i != j){
				suma += x[i][j];
			}
		}
		modelo.add(y[j][j] == suma);
		suma.end();
	}

	//ambiente.out()  << "Restriccion 2 ingresada correctamente" << endl;




	for (int j = 0; j < datos.n; j++)
	{
		IloExpr suma_izq(ambiente);
		IloExpr suma_der(ambiente);

		for (int i = 0; i < datos.n; i++)
		{
			if (i != j){
				suma_izq += x[i][j];
			}
		}

		for (int k = 0; k < datos.n; k++)
		{
			if (k != j){
				suma_der += x[j][k];
			}
		}


		modelo.add(suma_izq == suma_der);
		suma_izq.end();
		suma_der.end();
	}

	//ambiente.out()  << "Restriccion 3 ingresada correctamente" << endl;

	for (int p = 0; p < datos.numero_clusters; p++)
	{
		IloExpr suma(ambiente);
		for (std::set<IloInt>::iterator j = datos.clusters[p].begin(); j != datos.clusters[p].end(); j++)
		{
			for (int i = 0; i < datos.n; i++)
			{
				if (datos.clusters[p].find(i) == datos.clusters[p].end() && i != *j){
					suma += x[i][*j];
				}
			}
			
		}
		modelo.add(suma == 1);
		suma.end();
	}

	//ambiente.out()  << "Restriccion 4 ingresada correctamente" << endl;

	for (int p = 0; p < datos.numero_clusters; p++)
	{
		for (std::set<IloInt>::iterator i = datos.clusters[p].begin(); i != datos.clusters[p].end(); i++){
			for (std::set<IloInt>::iterator j = datos.clusters[p].begin(); j != datos.clusters[p].end(); j++)
			{
				if (*i != *j){
					modelo.add(y[*i][*j] <= y[*j][*j]);
				}
			}
		}
	}

	//ambiente.out()  << "Restriccion 8 ingresada correctamente" << endl;



	for (int j = 0; j < datos.n; j++){
		for (int k = 0; k < datos.n; k++){

			if (j != datos.deposito && k != datos.deposito && j != k){
				IloExpr suma_izq(ambiente);
				IloExpr suma_der(ambiente);

				for (int i = 0; i < datos.n; i++)
				{
					if (i != j && i != datos.deposito){
						suma_izq += f[i][j][k];
					}
				}

				for (int h = 0; h < datos.n; h++)
				{
					if (h != j){
						suma_der += f[j][h][k];
					}
				}


				modelo.add(suma_izq == suma_der);
				suma_izq.end();
				suma_der.end();

			}
		}


	}


	//ambiente.out()  << "Restriccion 5 ingresada correctamente" << endl;



	for (int p = 0; p < datos.numero_clusters; p++)
	{
		for (std::set<IloInt>::iterator j = datos.clusters[p].begin(); j != datos.clusters[p].end(); j++)
		{
			for (std::set<IloInt>::iterator k = datos.clusters[p].begin(); k != datos.clusters[p].end(); k++)
			{
				if (*k != *j){
					modelo.add(f[*k][*j][*k] <= y[*k][*j] + x[*k][*j]);}
			}
		}
	}

	for (int i = 0; i < datos.n; i ++)
	{
		for (int j = 0; j < datos.n; j++)
		{
			for (int k = 0; k < datos.n; k ++)
			{
				if (i != k){
					modelo.add(f[i][j][k] <= x[i][j]);}
			}
		}
	}



	//ambiente.out()  << "Restriccion 6 ingresada correctamente" << endl;


	for (int p = 0; p < datos.numero_clusters; p++)
	{
		for (std::set<IloInt>::iterator k = datos.clusters[p].begin(); k != datos.clusters[p].end(); k++){
			for (int j = 0; j < datos.n; j++)
			{
				if (datos.clusters[p].find(j) == datos.clusters[p].end() && *k != j){
					modelo.add(f[*k][j][*k] <= x[*k][j]);
				}
			}
		}
	}

	//ambiente.out()  << "Restriccion 7 ingresada correctamente" << endl;

	
	for (int k = 0; k < datos.n; k++)
	{
		if (k != datos.deposito)
		{
			IloExpr suma(ambiente);
			for (int j = 0; j < datos.n; j++)
			{
				if (k != j){
					suma += f[k][j][k];
				}
			}
			modelo.add(suma == 1);
			suma.end();
		}
	}

	for (int k = 0; k < datos.n; k++)
	{
		if (k != datos.deposito)
		{
			IloExpr suma(ambiente);
			for (int i = 0; i < datos.n; i++)
			{
				if (i != datos.deposito){
					suma += f[i][datos.deposito][k];
				}
			}
			modelo.add(suma == 1);
			suma.end();
		}
	}

	//ambiente.out()  << "Restriccion 8 ingresada correctamente" << endl;
	
	

	clock_t tiempo_inicio = clock();
	cplex.extract(modelo);
	cplex.setParam(IloCplex::TiLim, TimeLimit);
	cplex.setParam(IloCplex::WorkMem, 24631);
	start_time_det = cplex.getDetTime();
	start_time_cplex = cplex.getCplexTime();
	if(cplex.solve())
		solver_success = true;
	clock_t tiempo = clock() - tiempo_inicio;

	//cplex.exportModel("modelo_2.lp");
	//cplex.writeSolution("solucion.sol");

	saveResult(file_name, multiplicador);
}


void SolverMMulti::showSolution()
{
	if(solver_success)
	{
		IloNum Valor_Z1 = 0;
		IloNum Valor_Z2 = 0;
		
		ambiente.out() << "\nValores de X:\n";
		for (int i = 0; i < datos.n; i++)
		{
			for (int j = 0; j < datos.n; j++)
			{
				if (i != j){
					if (cplex.getValue(x[i][j]) >= 0.98)
					{

						ambiente.out() << "\nX[" << i + 1 << "][" << j + 1 << "] = " << cplex.getValue(x[i][j]);
						Valor_Z1 += cplex.getValue(x[i][j]) * datos.c[i][j];
					}
				}

			}
		}

		ambiente.out() << "\n\nValores de Y:\n";

		for (int i = 0; i < datos.n; i++)
		{
			for (int j = 0; j < datos.n; j++)
			{
					try
					{
						if (cplex.getValue(y[i][j]) >= 0.98)
						{
						
							ambiente.out() << "\nY[" << i + 1 << "][" << j + 1 << "] = " << cplex.getValue(y[i][j]);
							
							
						}
					}
					catch (...){}

			}
		}

		ambiente.out() << "\n\nValores de F:\n";
		for (int i = 0; i < datos.n; i++)
		{
			for (int j = 0; j < datos.n; j++)
			{
				for (int k = 0; k < datos.n; k++)
				{
					if (i != j){
						try
						{
							if (cplex.getValue(f[i][j][k]) >= 0.98)
							{
								ambiente.out() << "\nf[" << i + 1 << "][" << j + 1 << "][" << k + 1 << "] = " << cplex.getValue(f[i][j][k]);
							}
						}
						catch(...){}
					}
				}
			}
		}

		for (int i = 0; i < datos.n; i++)
		{
			for (int j = 0; j < datos.n; j++)
			{
				
				if (i != j){
					try
					{
						Valor_Z2 += 3*datos.t[i][j] * cplex.getValue(y[i][j]);
					}
					catch(...){}
				}
			}
		}
		for (int i = 0; i < datos.n; i++)
		{
			for (int j = 0; j < datos.n; j++)
			{
				for (int k = 0; k < datos.n; k++)
				{
					if (i != j && i != k){
						try
						{
							Valor_Z2 += datos.t[i][j] * cplex.getValue(f[i][j][k]);
						}
						catch(...){}
					}
				}
			}
		}

		for (int p = 0; p < datos.numero_clusters; p++)
		{
			for (std::set<IloInt>::iterator k = datos.clusters[p].begin(); k != datos.clusters[p].end(); k++){
				for (int j = 0; j < datos.n; j++)
				{
					if (datos.clusters[p].find(j) == datos.clusters[p].end() && *k != j){
						try
						{
							Valor_Z2 += datos.t[*k][j] * cplex.getValue(f[*k][j][*k]);
						}
						catch(...){}
					}
				}
			}
		}

		for (int p = 0; p < datos.numero_clusters; p++)
		{
			for (std::set<IloInt>::iterator i = datos.clusters[p].begin(); i != datos.clusters[p].end(); i++)
			{
				for (std::set<IloInt>::iterator j = datos.clusters[p].begin(); j != datos.clusters[p].end(); j++)
				{

					if (*i != *j){
						try
						{
							Valor_Z2 += datos.t[*i][*j] * cplex.getValue(x[*i][*j]);
						}
						catch(...){}
							
					}
				}
			}
		}

		ambiente.out() << "\nZ: " << cplex.getBestObjValue() <<"\n";
		ambiente.out() << "Valor Z1: " << Valor_Z1 << "\n";
		ambiente.out() << "Valor Z2: " << Valor_Z2 << "\n";

	}
	



}

void SolverMMulti::graphSolution(std::string file_name)
{
	if(solver_success)
	{
		std::ofstream grafo(file_name, std::ios::out);
		if (!grafo.is_open())
			throw std::runtime_error("Archivo " + file_name + " no abierto");

		grafo << "graph\n[\n\thierarchic\t1\n\tdirected\t1\n";
		for (int i = 0; i < datos.n; i++)
		{
			grafo << "\tnode\n\t[\n\t\tid\t"<< i + 1 
				<< "\n\t\tgraphics\n\t\t[\n\t\t\tx\t" << datos.xcoord[i]
				<< "\n\t\t\ty\t" << datos.ycoord[i]
				<< "\n\t\t\tw\t" << 50
				<< "\n\t\t\th\t" << 50
				<< "\n\t\t\ttype\t\x22" << "ellipse"
				<< "\x22\n\t\t\tfill\t\x22" << datos.color[i]
				<< "\x22\n\t\t]\n\t\tLabelGraphics\n\t\t[\n\t\t\ttext\t\x22" << i + 1
				<< "\x22\n\t\t\tfontSize\t" << 20
				<< "\n\t\t]\n\t]\n";
		}


		for (int i = 0; i < datos.n; i++)
		{
			for (int j = 0; j < datos.n; j++)
			{
				if (i != j){
					if (cplex.getValue(x[i][j]) >= 0.98)
					{
						grafo << "\tedge\n\t[\n\t\tsource\t" << i + 1
							<< "\n\t\ttarget\t" << j + 1
							<< "\n\t\tgraphics\n\t\t[\n\t\t\tfill\t\x22" << "#000000"
							<< "\x22\n\t\t\ttargetArrow\t\x22" << "standard"
							<< "\x22\n\t\t]\n\t]\n";

					}
				}

			}
		}

		for (int i = 0; i < datos.n; i++)
		{
			for (int j = 0; j < datos.n; j++)
			{
				if (i != j){
					try
					{
						if (cplex.getValue(y[i][j]) >= 0.98)
						{
							grafo << "\tedge\n\t[\n\t\tsource\t" << i + 1
								<< "\n\t\ttarget\t" << j + 1
								<< "\n\t\tgraphics\n\t\t[\n\t\t\tstyle\t\x22" << "dashed"
								<< "\x22\n\t\t\tfill\t\x22" << "#000000"
								<< "\x22\n\t\t\ttargetArrow\t\x22" << "standard"
								<< "\x22\n\t\t]\n\t]\n";

						}
					}
					catch (...){}
				}
			}
		}
	}
}

void SolverMMulti::saveResult(std::string file_name, IloInt multiplicador)
{
	std::ofstream salidas(file_name, std::ios::out | std::ios::app);
	if (!salidas.is_open())
		throw std::runtime_error("Archivo " + file_name + " no abierto");
	
	IloNum Valor_Z1 = 0;
	IloNum Valor_Z2 = 0;
	
	for (int i = 0; i < datos.n; i++)
	{
		for (int j = 0; j < datos.n; j++)
		{
			if (i != j){
				if (cplex.getValue(x[i][j]) >= 0.98)
				{
					Valor_Z1 += cplex.getValue(x[i][j]) * datos.c[i][j];
				}
			}
		}
	}


	for (int i = 0; i < datos.n; i++)
		{
			for (int j = 0; j < datos.n; j++)
			{
				
				if (i != j){
					try
					{
						Valor_Z2 += 3*datos.t[i][j] * cplex.getValue(y[i][j]);
					}
					catch(...){}
				}
			}
		}
		for (int i = 0; i < datos.n; i++)
		{
			for (int j = 0; j < datos.n; j++)
			{
				for (int k = 0; k < datos.n; k++)
				{
					if (i != j && i != k){
						try
						{
							Valor_Z2 += datos.t[i][j] * cplex.getValue(f[i][j][k]);
						}
						catch(...){}
					}
				}
			}
		}

		for (int p = 0; p < datos.numero_clusters; p++)
		{
			for (std::set<IloInt>::iterator k = datos.clusters[p].begin(); k != datos.clusters[p].end(); k++){
				for (int j = 0; j < datos.n; j++)
				{
					if (datos.clusters[p].find(j) == datos.clusters[p].end() && *k != j){
						try
						{
							Valor_Z2 += datos.t[*k][j] * cplex.getValue(f[*k][j][*k]);
						}
						catch(...){}
					}
				}
			}
		}

		for (int p = 0; p < datos.numero_clusters; p++)
		{
			for (std::set<IloInt>::iterator i = datos.clusters[p].begin(); i != datos.clusters[p].end(); i++)
			{
				for (std::set<IloInt>::iterator j = datos.clusters[p].begin(); j != datos.clusters[p].end(); j++)
				{

					if (*i != *j){
						try
						{
							Valor_Z2 += datos.t[*i][*j] * cplex.getValue(x[*i][*j]);
						}
						catch(...){}
							
					}
				}
			}
		}

	
	std::string obj = "";
	if (multiplicador == 1)
		obj = "Z1";
	if (multiplicador == 0)
		obj = "Z2";


	salidas << model_name << "," 
			<< datos.instance_name << "," 
			<< datos.n << "," 
			<< datos.numero_clusters << "," 
			<< datos.factor << "," 
			<< cplex.getStatus() << "," 
			<< obj << ","
			<< Valor_Z1 << "," 
			<< Valor_Z2 << "," 
			<< cplex.getObjValue() << "," 
			<< cplex.getBestObjValue() << "," 
			<< 100 * (cplex.getObjValue() - cplex.getBestObjValue()) 
						/ cplex.getObjValue() << "," 
			<< cplex.getNbinVars() << ","
        	<< cplex.getNintVars() << ","
			<< cplex.getNrows() << ","
        	<< cplex.getNcols() << ","
			<< cplex.getNnodes() << ","
			<< cplex.getNnodesLeft() << ","
			<< cplex.getTime() << ","
			<< cplex.getCplexTime() << "," 
			<< cplex.getCplexTime() - start_time_cplex << ","
			<< cplex.getDetTime() << ","
			<< cplex.getDetTime() - start_time_det << "\n";
			//<< "," << tiempo 
	salidas.close();

}

