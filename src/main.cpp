//: Cplex.cpp

#include <ilcplex/ilocplex.h>
#include "include/data.h"
#include "include/solver_flow.h"
#include "include/solver_multi.h"
#include "include/solver_mtz.h"
#include "include/solver_mtz_m.h"
#include "include/solver_multi_m.h"

#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int argc, char **argv)
{

	IloEnv ambiente;
	std::string instance = "parametros.dat";
	std::string save_file = "results.csv";
	std::string obj = "Z1";
	std::string graph = "grafo.gml";
	bool show_parameter = false;
	bool show_solution = true;
	bool graph_bool = false;
	bool solve_problem = true;
	IloInt mult = 1;
	IloNum Max_Z1 = 10000000000;
	IloNum Max_Z2 = 10000000000;
		try
		{
			po::options_description desc{"Opciones"};
			desc.add_options()
			("help,h", "Pantalla de ayuda")
				("instance,i", po::value<std::string>(), "Instancia del modelo")
				("save-file,s", po::value<std::string>(), "Archivo csv de guardado")
				("obj,o", po::value<std::string>(), "Funcion objetivo a optimizar")
				("graph,g", po::value<std::string>(), "Realiza el grafo de la solucion en archivo indicado")
				("Max-Z1", po::value<IloNum>(), "Z1 Máximo")
				("Max-Z2", po::value<IloNum>(), "Z2 Máximo")
				("show-parameters", "Mostrar los parametros del modelo")
				("hide-solution", "Oculta la solucion")
				("no-solve", "Desactiva resolucion del problema");

			po::variables_map vm;
			po::store(
				po::parse_command_line(argc, argv, desc), vm);
			po::notify(vm);

			if (vm.count("help")){
				std::cout << desc << '\n';
				return EXIT_SUCCESS;
			}
			else {
				if (vm.count("instance")){
					instance = vm["instance"].as<std::string>();
					std::cout << "instance: " << instance << '\n';
				}
				if (vm.count("save-file")){
					save_file = vm["save-file"].as<std::string>().c_str();
					std::cout << "save-file: " << save_file << '\n';
				}
				if (vm.count("Max-Z1")){
					Max_Z1 = vm["Max-Z1"].as<IloNum>();
					std::cout << "Max_Z1: " << Max_Z1 << '\n';
				}
				if (vm.count("Max-Z2")){
					Max_Z2 = vm["Max-Z2"].as<IloNum>();
					std::cout << "Max_Z2: " << Max_Z2 << '\n';
				}
				if (vm.count("obj")){
					obj = vm["obj"].as<std::string>().c_str();
					if (obj == "Z2")
						mult = 0;
					if (obj != "Z1" && obj != "Z2")
						std::cerr << "Error en lectura de Fo. Se cargará Z1 \n";
					std::cout << "obj: " << obj << '\n';
				}
				if (vm.count("graph")){
					graph = vm["graph"].as<std::string>().c_str();
					graph_bool = true;
					std::cout << "graph: " << graph << '\n';
				}
				if (vm.count("show-parameters")){
					show_parameter = true;
				}
				if (vm.count("hide-solution")){
					show_solution = false;
				}
				if (vm.count("no-solve")){
					solve_problem = false;
				}
			}
			
		}
		catch (const po::error &ex)
		{
			std::cerr << ex.what() << '\n';
			return EXIT_FAILURE;
		}

	try
	{
		Data datos(instance, ambiente);
		//SolverFlow solver(ambiente, datos);
		//SolverMulti solver(ambiente, datos);
		//SolverMTZ solver(ambiente, datos);
		SolverMMTZ solver(ambiente, datos);
		//SolverMMulti solver(ambiente, datos);
		if(show_parameter)
			datos.showData();
		if(solve_problem)
			solver.solve_model(save_file, mult, Max_Z1, Max_Z2);
		if (show_solution)
			solver.showSolution();
		if (graph_bool)
			solver.graphSolution(graph);
	}
	catch(const std::exception& e)
	{
		std::cerr << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
		
	ambiente.end();
	return EXIT_SUCCESS;
}


