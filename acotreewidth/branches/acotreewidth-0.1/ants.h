#include <string>
#include <vector>
#include <list>
#include <map>
#include "graph.h"

class PheromoneMatrix : public Matrix<double> {
  private:
    double evaporation_rate_;
  public:
    PheromoneMatrix(int vertices, double evaporation_rate);
    void evaporate(unsigned int v, unsigned int w);
    void evaporate_all();
    unsigned int size();
};

class Tour {
  private:
    std::vector<unsigned int> *vertices_;
    unsigned int capacity_;
    double length_;
  public:
    Tour(unsigned int vertices);
    ~Tour();
    double get_length();
    void set_length(double length);
    std::vector<unsigned int> get_vertices();
    Tour &add_vertex(unsigned int vertex);
    unsigned int &operator[](const unsigned int vertex) const;
    unsigned int size() const;
    unsigned int capacity();
    void clear();
    Tour &operator=(const Tour &t);
    bool operator<(const Tour &t);
};

class OptimizationProblem {
  public:
    virtual unsigned int get_max_tour_size() = 0;
    virtual unsigned int number_of_vertices() = 0;
    virtual std::map<unsigned int,double> get_feasible_start_vertices() = 0;
    virtual std::map<unsigned int,double> get_feasible_neighbours(unsigned int vertex) = 0;
    virtual double eval_tour(Tour &tour) = 0;
    virtual double pheromone_update(Tour &tour) = 0;
    virtual void added_vertex_to_tour(unsigned int vertex) = 0;
    virtual bool is_tour_complete(Tour &tour) = 0;
    virtual void cleanup() = 0;
};

class Ant {
  protected:
    struct MultiMapComp {
      bool operator()(const double &p1, const double &p2) const;
    };
    Tour *tour;
    virtual void local_pheromone_update() = 0;
    void update_tour_length(OptimizationProblem &op);
    void add_vertex_to_tour(OptimizationProblem &op, unsigned int vertex);
    std::multimap<double,unsigned int,MultiMapComp> get_feasible_vertices(OptimizationProblem &op, PheromoneMatrix &pheromones, double alpha, double beta);
    unsigned int current_vertex();
    unsigned int choose_next_vertex_with_likelihood(std::multimap<double,unsigned int,MultiMapComp> probabilities);
  public:
    Ant(unsigned int vertices);
    Ant(const Ant &ant);
    Ant &operator=(const Ant &ant);
    bool operator<(const Ant &ant);
    ~Ant();
    double get_tour_length();
    std::vector<unsigned int> get_vertices();
    void reset();
    void apply_local_search(OptimizationProblem &op);
    virtual void construct_solution(OptimizationProblem &op, PheromoneMatrix &pheromones, double alpha, double beta) = 0;
    virtual void offline_pheromone_update(OptimizationProblem &op, PheromoneMatrix &pheromones) = 0;
};

class SimpleAnt : virtual public Ant {
  protected:
    void local_pheromone_update();
  public:
    SimpleAnt(unsigned int vertices);
    SimpleAnt(const SimpleAnt &ant);
    SimpleAnt &operator=(const SimpleAnt &ant);
    void construct_solution(OptimizationProblem &op, PheromoneMatrix &pheromones, double alpha, double beta);
    void offline_pheromone_update(OptimizationProblem &op, PheromoneMatrix &pheromones);
};

template<class T> class AntColony {
  private:
    T *best_so_far_;
    T *best_iteration_;
    double alpha_;
    double beta_;

    void construct_ants_solutions() {
      for(typename std::list<T>::iterator it=ants_->begin();it!=ants_->end();it++) {
        T &ant = (*it);
        ant.construct_solution(*problem_, *pheromones_, alpha_, beta_);
        problem_->cleanup();
      }
    }

    void apply_local_search() {
      for(typename std::list<T>::iterator it=ants_->begin();it!=ants_->end();it++) {
        T &ant = (*it);
        ant.apply_local_search(*problem_);
      }
    }

    void reset_ants() {
      for(typename std::list<T>::iterator it=ants_->begin();it!=ants_->end();it++) {
        it->reset();
      }
    }

    void update_best_tours() {
      ants_->sort();
      (*best_iteration_) = ants_->front();
      if(best_so_far_->get_tour_length() > best_iteration_->get_tour_length()) {
        (*best_so_far_) = (*best_iteration_);
      }
    }

    virtual void update_pheromones() = 0;
  protected:
    PheromoneMatrix *pheromones_;
    std::list<T> *ants_;
    OptimizationProblem *problem_;
  public:
    AntColony(OptimizationProblem *op, int number_of_ants=25, double alpha=1, double beta=1, double evaporation_rate=0.1) {
      problem_ = op;
      ants_ = new std::list<T>(number_of_ants, T(problem_->get_max_tour_size()));
      pheromones_ = new PheromoneMatrix(problem_->get_max_tour_size()+1, evaporation_rate);
      alpha_ = alpha;
      beta_ = beta;
      best_so_far_ = new T(problem_->get_max_tour_size());
      best_iteration_ = new T(problem_->get_max_tour_size());
    }

    ~AntColony() {
      delete problem_;
      delete ants_;
      delete pheromones_;
      delete best_so_far_;
      delete best_iteration_;
    };

    void run() {
      construct_ants_solutions();
      apply_local_search();
      update_pheromones();
      update_best_tours();
      reset_ants();
    }

    std::vector<unsigned int> get_best_tour() {
      return best_so_far_->get_vertices();
    }

    std::vector<unsigned int> get_best_tour_in_iteration() {
      return best_iteration_->get_vertices();
    }

    double get_best_tour_length() {
      return best_so_far_->get_tour_length();
    }

    double get_best_tour_length_in_iteration() {
      return best_iteration_->get_tour_length();
    }
};

class SimpleAntColony : public AntColony<SimpleAnt> {
  public:
    SimpleAntColony(OptimizationProblem *op, int number_of_ants=10, double alpha=1, double beta=100, double evaporation_rate=0.1);
  protected:
    void update_pheromones();
};
