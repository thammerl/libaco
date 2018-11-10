#include <limits>
#include <iostream>
#include <list>
#include <cmath>
#include "ants.h"
#include "util.h"

PheromoneMatrix::PheromoneMatrix(int vertices, double evaporation_rate) : Matrix<double>(vertices, 1.0) {
  evaporation_rate_ = evaporation_rate;
}

void PheromoneMatrix::evaporate(unsigned int v, unsigned int w) {
  (*matrix_)[v][w] *= 1 - evaporation_rate_;
}

void PheromoneMatrix::evaporate_all() {
  for(int i=0;i<matrix_->size();i++) {
    for(int j=0;j<matrix_->size();j++) {
      evaporate(i,j);
    }
  }
}

unsigned int PheromoneMatrix::size() {
  return matrix_->size();
}

Tour::Tour(unsigned int vertices) {
  vertices_ = new std::vector<unsigned int>();
  capacity_ = vertices;
  length_ = UINT_MAX;
}

Tour::~Tour() {
  delete vertices_;
}

double Tour::get_length() {
  return length_;
}

std::vector<unsigned int> Tour::get_vertices() {
  return *vertices_;
}

void Tour::set_length(double length) {
  length_ = length;
}

Tour &Tour::add_vertex(unsigned int vertex) {
  vertices_->push_back(vertex);
  return *this;
}

unsigned int Tour::size() const {
  return vertices_->size();
}

unsigned int Tour::capacity() {
  return capacity_;
}

unsigned int &Tour::operator[](const unsigned int vertex) const {
  return (*vertices_)[vertex];
}

void Tour::clear() {
  vertices_->clear();
  length_ = std::numeric_limits<double>::max();
}

Tour &Tour::operator=(const Tour &t) {
  this->vertices_->resize(t.size());
  for(unsigned int i=0;i<this->size();i++) {
    (*this)[i] = t[i];
  }
  this->length_ = t.length_;
  this->capacity_ = t.capacity_;
}

bool Tour::operator<(const Tour &t) {
  return this->length_ < t.length_;
}

Ant::Ant(unsigned int vertices) {
  tour = new Tour(vertices);
}

Ant::Ant(const Ant &ant) {
  tour = new Tour(ant.tour->capacity());
  for(unsigned int i=0;i<ant.tour->size();i++) {
    (*tour)[i] = (*ant.tour)[i];
  }
}

Ant &Ant::operator=(const Ant &ant) {
  (*tour) = (*ant.tour);
}

bool Ant::operator<(const Ant &ant) {
  return (*this->tour) < (*ant.tour);
}

Ant::~Ant() {
  delete tour;
}

double Ant::get_tour_length() {
  return tour->get_length();
}

std::vector<unsigned int> Ant::get_vertices() {
  return tour->get_vertices();
}

void Ant::add_vertex_to_tour(OptimizationProblem &op, unsigned int vertex) {
  tour->add_vertex(vertex);
  op.added_vertex_to_tour(vertex);
}

unsigned int Ant::current_vertex() {
  if(tour->size() > 0) {
    return (*tour)[tour->size()-1];
  } else {
    return -1;
  }
}

std::multimap<double,unsigned int,Ant::MultiMapComp> Ant::get_feasible_vertices(OptimizationProblem &op, PheromoneMatrix &pheromones, double alpha, double beta) {
  std::map<unsigned int,double> vertices;
  if(current_vertex() == -1) {
    vertices = op.get_feasible_start_vertices();
  } else {
    vertices = op.get_feasible_neighbours(current_vertex());
  }
  double base = 0.0;
  std::map<unsigned int,double>::iterator it;
  for(it=vertices.begin();it!=vertices.end();it++) {
    unsigned int vertex = (*it).first;
    double heuristic_value = (*it).second;
    if(current_vertex() == -1) {
      base += pow(pheromones[pheromones.size()-1][vertex], alpha) * pow(heuristic_value, beta);
    } else {
      base += pow(pheromones[current_vertex()][vertex], alpha) * pow(heuristic_value, beta);
    }
  }

  std::multimap<double,unsigned int,MultiMapComp> probabilities;
  for(it=vertices.begin();it!=vertices.end();it++) {
    //std::cout << pheromones[pheromones.size()-1][(*it).first] << std::endl;
    unsigned int vertex = (*it).first;
    double heuristic_value = (*it).second;
    if(current_vertex() == -1) {
      (*it).second = pow(pheromones[pheromones.size()-1][vertex], alpha) * pow(heuristic_value, beta) / base;
    }
    else {
      (*it).second = pow(pheromones[current_vertex()][vertex], alpha) * pow(heuristic_value, beta) / base;
    }
    //std::cout << "vertex: " << (*it).first << " heuristic: " << heuristic_value << " probability: " << (*it).second << std::endl;
    probabilities.insert(std::pair<double,unsigned int>((*it).second, (*it).first));
  }
  //std::cout << "***************************************" << std::endl;
  return probabilities;
}

void Ant::update_tour_length(OptimizationProblem &op) {
  tour->set_length(op.eval_tour(*tour));
}


bool Ant::MultiMapComp::operator()(const double &p1,const double &p2) const {
  return p1 > p2;
}

unsigned int Ant::choose_next_vertex_with_likelihood(std::multimap<double,unsigned int,MultiMapComp> probabilities) {
  unsigned int vertex = -1;
  unsigned int random_value = Util::random_number(RAND_MAX);
  double tmp = 0;
  for(std::multimap<double,unsigned int,MultiMapComp>::iterator it=probabilities.begin();it!=probabilities.end();) {
    double probability = (*it).first;
    std::vector<unsigned int> vertices;
    while(it!=probabilities.end() && (*it).first == probability) {
      vertices.push_back((*it).second);
      it++;
    }
    if(tmp+probability*vertices.size()*RAND_MAX > random_value) {
      unsigned int v = Util::random_number(vertices.size());
      vertex = vertices[v];
      break;
    }
    tmp += probability*vertices.size()*RAND_MAX;
  }
  return vertex;
}

void Ant::apply_local_search(OptimizationProblem &op) {
  //TODO: local search
}

void Ant::reset() {
  tour->clear();
}

SimpleAnt::SimpleAnt(unsigned int vertices) : Ant(vertices) {
}

SimpleAnt::SimpleAnt(const SimpleAnt &ant) : Ant(ant) {
}

SimpleAnt &SimpleAnt::operator=(const SimpleAnt &ant) {
  Ant::operator=(ant);
  return (*this);
}

void SimpleAnt::construct_solution(OptimizationProblem &op, PheromoneMatrix &pheromones, double alpha, double beta) {
  while(!op.is_tour_complete(*tour)) {
    std::multimap<double,unsigned int,MultiMapComp> vertices = get_feasible_vertices(op, pheromones, alpha, beta);
    unsigned int vertex = choose_next_vertex_with_likelihood(vertices);
    add_vertex_to_tour(op, vertex);
  }
  update_tour_length(op);
  /*for (int i=0;i<tour->size();i++) {
    std::cout << (*tour)[i] << ",";
  }
  std::cout << "width: " << tour->get_length() << std::endl;*/
  local_pheromone_update();
}

void SimpleAnt::offline_pheromone_update(OptimizationProblem &op, PheromoneMatrix &pheromones) {
  for(unsigned int i=0;i<tour->size();i++) {
    if(i==0) {
      pheromones[tour->size()][(*tour)[i]]+=op.pheromone_update(*tour);
    } else {
      pheromones[(*tour)[i-1]][(*tour)[i]]+=op.pheromone_update(*tour);
    }
  }
}

void SimpleAnt::local_pheromone_update() {
  // no local update by SimpleAnt
}

SimpleAntColony::SimpleAntColony(OptimizationProblem *op, int number_of_ants, double alpha, double beta, double evaporation_rate) : AntColony<SimpleAnt>(op, number_of_ants, alpha, beta, evaporation_rate) {
}

void SimpleAntColony::update_pheromones() {
  pheromones_->evaporate_all();
  for(std::list<SimpleAnt>::iterator it=ants_->begin();it!=ants_->end();it++) {
    SimpleAnt &ant = (*it);
    ant.offline_pheromone_update(*problem_, *pheromones_);
  }
}
