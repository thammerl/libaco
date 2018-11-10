#include <map>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <exception>
#include <map>
#include <stack>
#include <list>
#include <libaco/ants.h>
#include <libaco/util.h>
#include <liblocalsearch/localsearch.h>

template <class T> class EliminationGraph : public T {
  private:
    std::vector<bool> *eliminated_vertices_;
    std::vector< std::vector<unsigned int> > *neighbour_lists_;
    std::vector< std::vector<unsigned int> > *neighbour_lists_length_;
    std::stack<unsigned int> eliminations_;
  public:
    EliminationGraph(const T &graph) : T(graph) {
      eliminated_vertices_ = new std::vector<bool>(graph.number_of_vertices(), false);
      neighbour_lists_ = new std::vector< std::vector<unsigned int> >(graph.number_of_vertices(), std::vector<unsigned int>(graph.number_of_vertices(), 0));
      neighbour_lists_length_ = new std::vector< std::vector<unsigned int> >(graph.number_of_vertices()+1, std::vector<unsigned int>(graph.number_of_vertices(), 0));

      for(unsigned int i=0;i<graph.number_of_vertices();i++) {
        std::vector<unsigned int> neighbours = graph.get_neighbours(i);
        (*neighbour_lists_length_)[i][0] = neighbours.size();
        for(unsigned j=0;j<neighbours.size();j++) {
          (*neighbour_lists_)[i][j] = neighbours[j];
        }
      }
    }

    ~EliminationGraph() {
      delete eliminated_vertices_;
      delete neighbour_lists_;
      delete neighbour_lists_length_;
    }

    unsigned int number_of_eliminations() {
      return eliminations_.size();
    }

    void eliminate(unsigned int vertex) {
      std::vector<unsigned int> neighbours = T::get_neighbours(vertex);
      
      for(int k=0;k<T::number_of_vertices();k++) {
        (*neighbour_lists_length_)[k][eliminations_.size()+1] = (*neighbour_lists_length_)[k][eliminations_.size()];
      }

      for(int i=0;i<neighbours.size();i++) {
        for(int j=i+1;j<neighbours.size();j++) {
          if(!T::is_edge(neighbours[i], neighbours[j])) {
            T::add_edge(neighbours[i], neighbours[j]);
            (*neighbour_lists_)[neighbours[i]][(*neighbour_lists_length_)[neighbours[i]][eliminations_.size()+1]] = neighbours[j];
            (*neighbour_lists_length_)[neighbours[i]][eliminations_.size()+1]++;
            (*neighbour_lists_)[neighbours[j]][(*neighbour_lists_length_)[neighbours[j]][eliminations_.size()+1]] = neighbours[i];
            (*neighbour_lists_length_)[neighbours[j]][eliminations_.size()+1]++;
          }
        }
        T::remove_edge(vertex,neighbours[i]);
      }

      (*eliminated_vertices_)[vertex] = true;
      eliminations_.push(vertex);
    }

    void rollback(unsigned int eliminations) {
      unsigned int eliminations_old = eliminations_.size();
      while(eliminations_.size() > eliminations) {
        unsigned int vertex = eliminations_.top();
        eliminations_.pop();
        (*eliminated_vertices_)[vertex] = false;
      }

      for(unsigned int v=0;v<T::number_of_vertices();v++) {
        unsigned int neighbours_before = (*neighbour_lists_length_)[v][eliminations];
        unsigned int neighbours_after = (*neighbour_lists_length_)[v][eliminations_old];
        for(unsigned int i=0;i<neighbours_after;i++) {
          if(i < neighbours_before) {
            if(!(*eliminated_vertices_)[(*neighbour_lists_)[v][i]]) {
              T::add_edge(v, (*neighbour_lists_)[v][i]);
            }
          } else {
            T::remove_edge(v, (*neighbour_lists_)[v][i]);
          }
        }
      }
    }
};

typedef double (*heuristicf)(const Graph &graph, unsigned int vertex);

namespace Heuristic {
  double min_degree(const Graph &graph, unsigned int vertex);
  double min_fill(const Graph &graph, unsigned int vertex);
}

template <class T> class DecompProblem : public OptimizationProblem, public EvaluationFunction {
  protected:
    T *graph_;
    EliminationGraph<T> *elim_graph_;
    EliminationGraph<T> *eval_elim_graph_;
    std::map<unsigned int,bool> visited_vertices_;
    std::vector<unsigned int> vertex_width_;
    heuristicf heuristic_; 
  public:
    DecompProblem(T *graph, heuristicf heuristic=Heuristic::min_degree) {
      graph_ = graph;
      elim_graph_ = new EliminationGraph<T>(*graph);
      eval_elim_graph_ = new EliminationGraph<T>(*graph);
      vertex_width_.reserve(graph->number_of_vertices());
      heuristic_ = heuristic;
    }

    ~DecompProblem() {
      delete graph_;
      delete elim_graph_;
      delete eval_elim_graph_;
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

    double eval_tour(const std::vector<unsigned int> &tour) {
      double width = compute_width(tour);
      return width;
    }

    double pheromone_update(unsigned int v, double tour_length) {
      return (vertex_width_[v] / tour_length) * (1.0 / tour_length);
    }

    void added_vertex_to_tour(unsigned int vertex) {
      vertex_width_[vertex] = elim_graph_->get_neighbours(vertex).size();
      elim_graph_->eliminate(vertex);
      visited_vertices_[vertex] = true;
    }

    bool is_tour_complete(const std::vector<unsigned int> &tour) {
      return tour.size() == graph_->number_of_vertices();
    }

    double eval_solution(const std::vector<unsigned int> &solution) {
      double quality = 1.0 / eval_tour(solution);
      return quality;
    }

    std::vector<unsigned int> apply_local_search(const std::vector<unsigned int> &tour) {
      TwoOptNeighbourhood neighbourhood;
      HillClimbing climbing(tour, *this, neighbourhood);
      climbing.search_neighbourhood();
      return climbing.get_best_so_far_solution();
      //return tour;
    }

    void cleanup() {
      elim_graph_->rollback(0);
      visited_vertices_.clear();
    }

    virtual unsigned int compute_width(const std::vector<unsigned int> &tour) = 0;
};

template <class T> class TreeDecompProblem : public DecompProblem<T> {
  public:
    TreeDecompProblem(T *graph, heuristicf heuristic=Heuristic::min_degree) : DecompProblem<T>(graph, heuristic) {
    }

    unsigned int compute_width(const std::vector<unsigned int> &tour) {
      static std::vector<unsigned int> prev_tour;
      unsigned int width = 0;

      unsigned int rollback = 0;
      if(prev_tour.size() == tour.size()) {
        while(rollback < tour.size() && rollback < DecompProblem<T>::eval_elim_graph_->number_of_eliminations()) {
          if(tour[rollback] == prev_tour[rollback]) {
            rollback++;
          } else {
            break;
          }
        }
      }

      DecompProblem<T>::eval_elim_graph_->rollback(rollback);
      for(unsigned int i=0;i<tour.size();i++) {
        unsigned int w = DecompProblem<T>::eval_elim_graph_->get_neighbours(tour[i]).size();
        DecompProblem<T>::eval_elim_graph_->eliminate(tour[i]);
        if(w > width) {
          width = w;
        }

        if(tour.size()-i <= width) {
          break;
        }
      }

      prev_tour = tour;
      return width;
    }
};

template <class T> class HyperTreeDecompProblem : public DecompProblem<T> {
  private:
    HyperGraph *hypergraph_;
  public:
    HyperTreeDecompProblem(HyperGraph *hypergraph, heuristicf heuristic=Heuristic::min_degree) : DecompProblem<T>(&hypergraph->get_primal_graph<T>(), heuristic) {
      hypergraph_ = hypergraph;
    }

    ~HyperTreeDecompProblem() {
      delete hypergraph_;
    }

    unsigned int compute_width(const std::vector<unsigned int> &tour) {
      unsigned int width = 0;
      EliminationGraph<T> graph(*DecompProblem<T>::graph_);
      for(unsigned int i=0;i<tour.size();i++) {
        std::vector<unsigned int> vertices = graph.get_neighbours(tour[i]);
        vertices.push_back(tour[i]);
        // set covering
        unsigned int w = compute_greedy_hyperedge_covering(vertices).size();
        graph.eliminate(tour[i]);
        if(w > width) {
          width = w;
        }
      }
      return width;
    }

    std::vector<unsigned int> compute_greedy_hyperedge_covering(std::vector<unsigned int> vertices) {
      std::vector<unsigned int> edges;
      std::map< unsigned int, std::vector<unsigned int> > edges_vertices;
      while(!vertices.empty()) {
        for(unsigned int i=0;i<vertices.size();i++) {
          unsigned int vertex = vertices[i];
          std::vector<unsigned int> vertex_edges = hypergraph_->get_edges_for_vertex(vertex);
          for(unsigned int j=0;j<vertex_edges.size();j++) {
            edges_vertices[vertex_edges[j]].push_back(vertex);
          }
        }

        //get edge that covers a maximum of vertices
        unsigned int maximum_edge;
        std::vector<unsigned int> vertices_covered;
        for(std::map< unsigned int, std::vector<unsigned int> >::iterator it=edges_vertices.begin();it!=edges_vertices.end();it++) {
          if((*it).second.size() > vertices_covered.size()) {
            maximum_edge = (*it).first;
            vertices_covered = (*it).second;
          }
        }

        //remove vertices covered by edge from vertices
        for(unsigned int k=0;k<vertices_covered.size();k++) {
          for(std::vector<unsigned int>::iterator it2=vertices.begin();it2!=vertices.end();it2++) {
            if(*it2 == vertices_covered[k]) {
              vertices.erase(it2);
              break;
            }
          }
        }

        edges.push_back(maximum_edge);
      }
      return edges;
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

  HyperGraph &parse_hypertreelib(const char *filepath) throw(FileNotFoundException);
}
