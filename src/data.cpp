//: data.cpp

#include <ilcplex/ilocplex.h>
#include <set>
#include <vector>
#include "include/data.h"


void Data::showData()
{
	
	ambiente.out() << "n: " << n;
	ambiente.out() << "\ndeposito: " << deposito;
	ambiente.out() << "\nclusters: " << numero_clusters;
	ambiente.out() << "\nfactor: " << factor;

	ambiente.out() << "\n\n";

	for (int nodo = 0; nodo < n; nodo++)
	{
		ambiente.out() << "\nxcoord[" << nodo << "]= " << xcoord[nodo];;
	}

	ambiente.out() << "\n\n";

	for (int nodo = 0; nodo < n; nodo++)
	{
		ambiente.out() << "\nycoord[" << nodo << "]= " << ycoord[nodo];
	}

	ambiente.out() << "\n\n";

	for (int nodo = 0; nodo < n; nodo++)
	{
		ambiente.out() << "\ncolor[" << nodo << "]= " << color[nodo];
	}

	ambiente.out() << "\n\n";

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			ambiente.out() << "\nc[" << i << ", " << j << "]= " << c[i][j];
		}
	}

	ambiente.out() << "\n\n";

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			ambiente.out() << "\nt[" << i << ", " << j << "]= " << t[i][j];
		}
	}

	ambiente.out() << "\n\n";

	for (int i = 0; i < clusters.size(); i++)
	{
		ambiente.out() << "\nNp[" << i << "]: ";
		for (std::set<IloInt>::iterator iterador = clusters[i].begin(); iterador != clusters[i].end(); iterador++)
		{
			ambiente.out() << *iterador << ", ";
		}
	}

	ambiente.out() << "\n\n";

}

void Data::definirDimParam()
{
	xcoord = IloNumArray(ambiente, n);
	ycoord = IloNumArray(ambiente, n);
	c = IloArray <IloNumArray>(ambiente, n);
	t = IloArray <IloNumArray>(ambiente, n);

	for (int i = 0; i < n; i++){
		c[i] = IloNumArray(ambiente, n);
		t[i] = IloNumArray(ambiente, n);
	}
}

bool Data::verifyString(std::string data_string, std::string read_string)
{
	if (data_string != read_string)
		throw std::runtime_error("Error al leer parametro " + data_string);
}

Data::Data (std::string ruta_archivo, IloEnv & env)
{
	ambiente = env;
	instance_name = ruta_archivo;
	std::ifstream archivo(instance_name, std::ios::in);
	if (!archivo.is_open())
		throw std::runtime_error("Archivo " + instance_name + " no encontrado");
	
	std::string verificador;
	IloInt M = 100000;

	archivo >> verificador;
	verifyString("n:", verificador);
	archivo >> n;

	archivo >> verificador;
	verifyString("deposito:", verificador);
	archivo >> deposito;
	deposito = deposito - 1;

	archivo >> verificador;
	verifyString("clusters:", verificador);
	archivo >> numero_clusters;

	archivo >> verificador;
	verifyString("factor:", verificador);
	archivo >> factor;
	definirDimParam();

	archivo >> verificador;
	verifyString("xcoord[n]:", verificador);
	for (int i = 0; i < n; i++){
		archivo >> xcoord[i];
	}

	archivo >> verificador;
	verifyString("ycoord[n]:", verificador);
	for (int i = 0; i < n; i++){
		archivo >> ycoord[i];
	}

	archivo >> verificador;
	verifyString("color[n]:", verificador);
	std::string color_aux;
	for (int i = 0; i < n; i++)
	{
		archivo >> color_aux;
		color.push_back(color_aux);
	}

	for (int i = 0; i < numero_clusters; i++)
	{
		std::set<IloInt> aux;
		archivo >> verificador;
		verifyString("Np[" + std::to_string(i + 1) + "]:", verificador);
		std::string valor;
		for (int j = 0; j < n; j++)
		{
			archivo >> valor;
			if (valor != ";"){
				aux.insert(stoi(valor) - 1);
			}
			else {
				clusters.push_back(aux);
				aux.clear();
				break;
			}
		}
	}

	archivo >> verificador;
	verifyString("c[n,n]:", verificador);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j)
			{
				c[i][j] = M;
			}
			else
			{
				archivo >> c[i][j];
			}
		}
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			t[i][j] = c[i][j] * factor;
		}
	}

	archivo >> verificador;
	if (!archivo.eof())
			std::cerr << "Cuidado, revisar final de " << ruta_archivo << "\n";

	archivo.close();
}
