#include <iostream>
#include <fstream>
#include <acotreewidth/decomp.h>

DecompLocalSearch::DecompLocalSearch(std::vector<unsigned int> initial_solution, EvaluationFunction &eval_func, Neighbourhood &neighbourhood) : LocalSearch(initial_solution, eval_func, neighbourhood) {
}

void DecompLocalSearch::search_neighbourhood() {
  const std::vector<unsigned int> &solution = neighbourhood_->next_neighbour_solution();
  double quality = eval_func_->eval_solution(solution);
  if(quality > best_so_far_quality_) {
    best_so_far_solution_ = solution;
    best_so_far_quality_ = quality;
  }
  neighbourhood_->set_solution(best_so_far_solution_);
}

double Heuristic::min_degree(const Graph &graph, unsigned int vertex) {
  return 1.0 / (graph.get_degree(vertex) + 1);
}

double Heuristic::min_fill(const Graph &graph, unsigned int vertex) {
  std::vector<unsigned int> neighbours = graph.get_neighbours(vertex);
  unsigned int edges_to_fill = 0;
  for(unsigned int i=0;i<neighbours.size();i++) {
    for(unsigned int j=i;j<neighbours.size();j++) {
      if(i!=j && !graph.is_edge(i, j)) {
        edges_to_fill += 1;
      }
    }
  }
  return 1.0 / (edges_to_fill + 1);
}

FileNotFoundException::FileNotFoundException(const char *filepath) {
  filepath_ = filepath;
}

const char *FileNotFoundException::what() const throw() {
  return filepath_;
}

HyperGraph &Parser::parse_hypertreelib(const char *filepath) throw(FileNotFoundException) {
  HyperGraph *graph = new HyperGraph();
  std::ifstream file(filepath);
  char buf[1024*8];

  if(!file) {
    throw FileNotFoundException(filepath);
  }

  std::map<std::string, unsigned int> vertex_ids;
  int edge_nr = 0;
  int vertex_nr = 0;
  while(file.good()) {
    std::string line;
    std::string atom;
    std::size_t found;
    file.getline(buf, 1024*8);
    line = std::string(buf);
    //strip whitespaces
    found = line.find(" ");
    while(found != std::string::npos) {
      line = line.replace(int(found), 1, "");
      found = line.find(" ");
    }

    if (line.length() > 0 && line[0] != '%' && line[0] != '\r' && line[0] != '<') {
      //strip ')' character
      found = line.find(")");
      line.replace(int(found), 1, "");
      //strip '.' character
      found = line.find(".");
      if(found != std::string::npos) {
        line.replace(int(found), 1, "");
      }

      found = line.find("(");
      atom = line.substr(0, int(found)); 
      line.replace(0, int(found)+1, "");
      found = line.find(",");
      std::vector<unsigned int> vars;
      while(found != std::string::npos) {
        std::string var = line.substr(0, int(found));
        if(vertex_ids.find(var) == vertex_ids.end()) {
          vars.push_back(vertex_nr);
          graph->set_vertex_label(vertex_nr, var);
          vertex_ids[var] = vertex_nr;
          vertex_nr++;
        } else {
          vars.push_back(vertex_ids[var]);
        }
        line.replace(0, int(found) + 1, "");
        found = line.find(",");
      }
      graph->add_hyperedge(edge_nr, vars);
      graph->set_edge_label(edge_nr, atom);
      edge_nr++;
    }
  }

  file.close();
  return *graph;
}
