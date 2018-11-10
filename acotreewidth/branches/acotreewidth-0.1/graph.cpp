#include "graph.h"

AdjacencyMatrixGraph::AdjacencyMatrixGraph(unsigned int vertices) : Matrix<unsigned short int>(vertices, 0) {
}

AdjacencyMatrixGraph::AdjacencyMatrixGraph(const AdjacencyMatrixGraph &graph) : Matrix<unsigned short int>(graph) {
}

AdjacencyMatrixGraph &AdjacencyMatrixGraph::operator=(const AdjacencyMatrixGraph &graph) {
  Matrix<unsigned short int>::operator=(graph);
  return *this;
}

unsigned int AdjacencyMatrixGraph::number_of_vertices() const {
  return matrix_->size();
}

std::vector<unsigned int> AdjacencyMatrixGraph::get_neighbours(unsigned int vertex) const {
  std::vector<unsigned int> neighbours;
  for(int i=0;i<number_of_vertices();i++) {
    if((*matrix_)[vertex][i]) {
      neighbours.push_back(i);
    }
  }
  return neighbours;
}

void AdjacencyMatrixGraph::add_edge(unsigned int v, unsigned int w) {
  (*matrix_)[v][w] = 1;
  (*matrix_)[w][v] = 1;
}

bool AdjacencyMatrixGraph::is_edge(unsigned int v, unsigned int w) const {
  return (*matrix_)[v][w];
}

void AdjacencyMatrixGraph::remove_edge(unsigned int v, unsigned int w) {
  (*matrix_)[v][w] = 0;
  (*matrix_)[w][v] = 0;
}

unsigned int AdjacencyMatrixGraph::get_degree(unsigned int vertex) const {
  unsigned int degree=0;
  for(unsigned int i=0;i<(*matrix_).size();i++) {
    if((*matrix_)[vertex][i]) {
      degree+=1;
    }
  }
  return degree;
}

AdjacencyListGraph::AdjacencyListGraph(unsigned int vertices, unsigned short int default_value) {
  vertices_ = new std::vector< AdjacencyMap<unsigned short int> >(vertices, AdjacencyMap<unsigned short int>(default_value));
}

AdjacencyListGraph::~AdjacencyListGraph() {
  delete vertices_;
}

AdjacencyListGraph::AdjacencyListGraph(const AdjacencyListGraph &graph) {
  vertices_ = new std::vector< AdjacencyMap<unsigned short int> >(graph.number_of_vertices(), AdjacencyMap<unsigned short int>(0));
  *vertices_ = *graph.vertices_;
}

AdjacencyListGraph &AdjacencyListGraph::operator=(const AdjacencyListGraph &graph) {
  *vertices_ = *graph.vertices_;
  return *this;
}

void AdjacencyListGraph::add_edge(unsigned int v, unsigned int w) {
  (*vertices_)[v].add_neighbour(w, 1);
  (*vertices_)[w].add_neighbour(v, 1);
}

bool AdjacencyListGraph::is_edge(unsigned int v, unsigned int w) const {
  return (*vertices_)[v][w];
}

void AdjacencyListGraph::remove_edge(unsigned int v, unsigned int w) {
  (*vertices_)[v].remove_neighbour(w);
  (*vertices_)[w].remove_neighbour(v);
}

unsigned int AdjacencyListGraph::number_of_vertices() const {
  return vertices_->size();
}

std::vector<unsigned int> AdjacencyListGraph::get_neighbours(unsigned int vertex) const {
  return (*vertices_)[vertex].get_neighbours();
}

unsigned int AdjacencyListGraph::get_degree(unsigned int vertex) const {
  return (*vertices_)[vertex].size();
}
