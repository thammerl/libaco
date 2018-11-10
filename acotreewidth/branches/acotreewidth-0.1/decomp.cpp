#include "decomp.h"
#include <iostream>

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

