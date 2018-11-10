#include <map>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <exception>
#include "ants.h"
#include "util.h"

typedef double (*heuristicf)(const Graph &graph, unsigned int vertex);

namespace Heuristic {
  double min_degree(const Graph &graph, unsigned int vertex);
  double min_fill(const Graph &graph, unsigned int vertex);
}

template <class T> class DecompProblem : virtual public OptimizationProblem {
  private:
    T *graph_;
    T *elim_graph_;
    std::map<unsigned int,bool> visited_vertices_;
    heuristicf heuristic_; 
  public:
    DecompProblem(T *graph, heuristicf heuristic=Heuristic::min_degree) {
      graph_ = graph;
      elim_graph_ = new T(*graph);
      heuristic_ = heuristic;
    }

    ~DecompProblem() {
      delete graph_;
      delete elim_graph_;
    }

    unsigned int get_max_tour_size() {
      return graph_->number_of_vertices();
    }

    unsigned int number_of_vertices() {
      return graph_->number_of_vertices();
    }

    std::map<unsigned int,double> get_feasible_start_vertices() {
      std::map<unsigned int,double> vertices;
      for(unsigned int i=0;i<graph_->number_of_vertices();i++) {
        vertices[i] = heuristic_(*graph_, i);
      }
      return vertices;
    }

    std::map<unsigned int,double> get_feasible_neighbours(unsigned int vertex) {
      std::map<unsigned int,double> vertices;
      for(unsigned int i=0;i<graph_->number_of_vertices();i++) {
        if (!visited_vertices_[i]) {
          vertices[i] = heuristic_(*elim_graph_, i);
        }
      }
      return vertices;
    }

    double eval_tour(Tour &tour) {
      return compute_width(tour);
    }

    double pheromone_update(Tour &tour) {
      return 1.0 / tour.get_length();
    }

    void added_vertex_to_tour(unsigned int vertex) {
      eliminate(*elim_graph_, vertex);
      visited_vertices_[vertex] = true;
    }

    bool is_tour_complete(Tour &tour) {
      return tour.size() == graph_->number_of_vertices();
    }

    void cleanup() {
      *elim_graph_ = *graph_;
      visited_vertices_.clear();
    }

    void eliminate(Graph &graph, unsigned int vertex) {
      std::vector<unsigned int> neighbours = graph.get_neighbours(vertex);
      for(int i=0;i<neighbours.size();i++) {
        for(int j=i;j<neighbours.size();j++) {
          if(i!=j) {
            graph.add_edge(neighbours[i],neighbours[j]);
          }
        }
        graph.remove_edge(vertex,neighbours[i]);
      }
    }

    unsigned int compute_width(Tour &tour) {
      unsigned int width = 0;
      T graph(*graph_);
      for(unsigned int i=0;i<tour.size();i++) {
        unsigned int w = graph.get_neighbours(tour[i]).size();
        eliminate(graph, tour[i]);
        if(w > width) {
          width = w;
        }
      }
      return width;
    }

};

namespace Parser {
  template <class T> Graph &parse_dimacs(const char *filepath) throw(FileNotFoundException) {
    Graph *graph;
    int number_of_vertices;
    int vertex_a, vertex_b;
    char problem[5];
    char buf[1024];
    char flag;
    std::ifstream file(filepath);

    if(!file) {
      throw FileNotFoundException(filepath);
    }

    while(file.good()) {
      file >> flag;
      switch(flag) {
        case 'p':
          file >> problem;
          if (!strcmp(problem, "edge")) {
            file >> number_of_vertices;
            graph = new T(number_of_vertices);
          }
          break;
        case 'e':
          file >> vertex_a;
          file >> vertex_b;
          graph->add_edge(vertex_a-1, vertex_b-1);
          break;
        default:
          file.getline(buf, 1024);
          break;
      }
    }
    file.close();
    return *graph;
  }

  template <class T> Graph &parse_hypertreelib(const char *filepath) {
  }
}
