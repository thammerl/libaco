#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <tclap/CmdLine.h>
#include "decomp.h"

static std::string filepath;
static unsigned int ants = 10;
static unsigned int iterations = 5;
static double alpha = 1.0;
static double beta = 1.0;
static double rho = 0.1;
static unsigned int heuristic = 0;
static unsigned int graph_type = 0;
static bool print_tour_flag = false;
static double time_limit = 300;

static void parse_options(int argc, char *argv[]) {
  TCLAP::CmdLine cmd("Ant Colony Optimization for Tree Decomposition", ' ', "0.1");
  TCLAP::ValueArg<unsigned int> ants_arg ("m", "ants", "number of ants", false, ants, "integer");
  TCLAP::ValueArg<unsigned int> iterations_arg ("i", "iterations", "number of iterations", false, iterations, "positive integer");
  TCLAP::ValueArg<double> alpha_arg ("a", "alpha", "alpha (influence of pheromone trails)", false, alpha, "double");
  TCLAP::ValueArg<double> beta_arg("b", "beta", "beta (influence of heuristic information)", false, beta, "double");
  TCLAP::ValueArg<double> rho_arg("e", "rho", "pheromone trail evaporation rate", false, rho, "double");
  std::vector<unsigned int> allowed;
  allowed.push_back(0);
  allowed.push_back(1);
  TCLAP::ValuesConstraint<unsigned int> allowed_values( allowed );
  TCLAP::ValueArg<unsigned int> graph_type_arg("g", "graph", "0: AdjacencyMatrix 1: AdjacencyList", false, graph_type, &allowed_values);
  TCLAP::ValueArg<unsigned int> heuristic_arg("j", "heuristic", "0: min_degree 1: min_fill", false, heuristic, &allowed_values);
  TCLAP::ValueArg<std::string> filepath_arg("f", "file", "path to the graph file", true, "", "filepath");
  TCLAP::SwitchArg print_tour_arg("o", "printord", "print best elimination ordering in iteration");
  TCLAP::ValueArg<double> time_limit_arg("t", "time", "terminate after n seconds (after last iteration is finished)", false, time_limit, "double");
  cmd.add(ants_arg);
  cmd.add(iterations_arg);
  cmd.add(alpha_arg);
  cmd.add(beta_arg);
  cmd.add(rho_arg);
  cmd.add(graph_type_arg);
  cmd.add(heuristic_arg);
  cmd.add(filepath_arg);
  cmd.add(print_tour_arg);
  cmd.add(time_limit_arg);
  cmd.parse(argc, argv);
  ants = ants_arg.getValue();
  iterations = iterations_arg.getValue();
  alpha = alpha_arg.getValue();
  beta = beta_arg.getValue();
  rho = rho_arg.getValue();
  heuristic = heuristic_arg.getValue();
  graph_type = graph_type_arg.getValue();
  filepath = filepath_arg.getValue();
  print_tour_flag = print_tour_arg.getValue();
  time_limit = time_limit_arg.getValue();
}

heuristicf get_heuristic_function() {
  switch(heuristic) {
    case 0:
      return Heuristic::min_degree;
      break;
    case 1:
      return Heuristic::min_fill;
      break;
  }
}

template <class T> OptimizationProblem *get_optimization_problem(heuristicf heuristic_function) {
  OptimizationProblem *op;
  T &graph1 = (T &) Parser::parse_dimacs<T>(filepath.c_str());
  op = new DecompProblem<T>(&graph1, heuristic_function);
  return op;
}

void print_tour(std::vector<unsigned int> tour) {
  for(unsigned int i=0;i<tour.size();i++) {
    std::cout << tour[i] << ((i == (tour.size()-1)) ? "" : ",");
  }
}

double timer() {
  static bool initialized_time = false;
  static clock_t time;
  if(!initialized_time) {
    time = clock();
    initialized_time = true;
    return 0.0;
  } else {
    clock_t time_diff = clock() - time;
    double elapsed_time = time_diff * 1.0 / CLOCKS_PER_SEC;
    return elapsed_time;
  }
}

int main(int argc, char *argv[]) {
  try {
    parse_options(argc, argv);
  } catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    exit(0);
  }

  heuristicf heuristic_function = get_heuristic_function();

  OptimizationProblem *op;
  try {
    switch(graph_type) {
      case 0:
        op = get_optimization_problem<AdjacencyMatrixGraph>(heuristic_function);
        break;
      case 1:
        op = get_optimization_problem<AdjacencyListGraph>(heuristic_function);
        break;
    }
  }
  catch(FileNotFoundException e) {
    std::cerr << "error: could not open " << e.what()  << std::endl;
    exit(0);
  }

  SimpleAntColony *colony = new SimpleAntColony(op, ants, alpha, beta, rho);
  std::cout << "iter\ttime\tbest\tbest_it";
  std::cout << (print_tour_flag ? "\tordering" : "");
  std::cout << std::endl;
  timer();
  for(int i=0;i<iterations && timer() < time_limit;i++) {
    colony->run();
    std::cout << (i+1) << "\t" << timer() << "\t" << colony->get_best_tour_length() << "\t";
    std::cout << colony->get_best_tour_length_in_iteration();
    if(print_tour_flag) {
      std::cout << "\t";
      print_tour(colony->get_best_tour_in_iteration());
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << "best\tordering" << std::endl;
  std::cout << colony->get_best_tour_length() << "\t";
  print_tour(colony->get_best_tour());
  std::cout << std::endl;
  delete colony;
}
