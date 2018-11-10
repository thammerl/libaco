#include <vector>
#include <map>

template <class T> class Row {
  private:
    std::vector<T> *values_;
  public:
    Row(const Row &nh) {
      values_ = new std::vector<T>(*(nh.values_));
    }

    Row(unsigned int vertices, const T &initialValue) {
      values_ = new std::vector<T>(vertices, initialValue);
    }

    ~Row() {
      delete values_;
    }

    Row &operator=(const Row &nh) {
      *this->values_ = *nh.values_;
      return *this;
    }

    T &operator[](unsigned int i) {
      return (*values_)[i];
    }
};

template <class T> class Matrix {
  protected:
    std::vector< Row<T> > *matrix_;
  public:
    Matrix(unsigned int vertices, const T &initialValue) {
      matrix_ = new std::vector< Row<T> >(vertices, Row<T> (vertices, initialValue));
    };

    ~Matrix() {
      delete matrix_;
    }

    Matrix(const Matrix &matrix) {
      matrix_ = new std::vector < Row<T> >(matrix.matrix_->size(), Row<T>(matrix.matrix_->size(), 0));
      *matrix_ = *(matrix.matrix_);
    }

    Row<T> &operator[](unsigned int i) {
      return (*matrix_)[i];
    }

    Matrix<T> &operator=(const Matrix<T> &matrix) {
      *this->matrix_ = *matrix.matrix_;
      return *this;
    }
};

class Graph {
  public:
    virtual void add_edge(unsigned int v, unsigned int w) = 0;
    virtual bool is_edge(unsigned int v, unsigned int w) const = 0;
    virtual void remove_edge(unsigned int v, unsigned int w) = 0;
    virtual unsigned int number_of_vertices() const = 0;
    virtual std::vector<unsigned int> get_neighbours(unsigned int vertex) const = 0;
    virtual unsigned int get_degree(unsigned int vertex) const = 0;
};

class AdjacencyMatrixGraph : public Graph, public Matrix<unsigned short int> {
  public:
    AdjacencyMatrixGraph(unsigned int vertices);
    AdjacencyMatrixGraph(const AdjacencyMatrixGraph &graph);
    AdjacencyMatrixGraph &operator=(const AdjacencyMatrixGraph &graph);
    void add_edge(unsigned int v, unsigned int w);
    bool is_edge(unsigned int v, unsigned int w) const;
    void remove_edge(unsigned int v, unsigned int w);
    unsigned int number_of_vertices() const;
    std::vector<unsigned int> get_neighbours(unsigned int vertex) const;
    unsigned int get_degree(unsigned int vertex) const;
};

template <class T> class AdjacencyMap {
  private:
    std::map<unsigned int, T> values_;
    T default_value_;
  public:
    AdjacencyMap(T default_value) {
      default_value_ = default_value;
    }

    void add_neighbour(unsigned int vertex, T value) {
      values_[vertex] = value;
    }

    void remove_neighbour(unsigned int vertex) {
      values_.erase(vertex);
    }

    AdjacencyMap &operator=(const AdjacencyMap &map) {
      values_ = map.values_;
    }

    T &operator[](unsigned int vertex) {
      typename std::map<unsigned int, T>::iterator it = values_.find(vertex);
      if (it != values_.end()) {
        return (*it).second;
      } else {
        return default_value_;
      }
    }

    unsigned int size() {
      return values_.size();
    }

    std::vector<unsigned int> get_neighbours() {
      std::vector<unsigned int> neighbours;
      typename std::map<unsigned int, T>::iterator it;
      for(it=values_.begin();it!=values_.end();it++) {
        neighbours.push_back((*it).first);
      }
      return neighbours;
    }
};

class AdjacencyListGraph : public Graph {
  private:
    std::vector< AdjacencyMap<unsigned short int> > *vertices_;
  public:
    AdjacencyListGraph(unsigned int vertices, unsigned short int default_value=0);
    ~AdjacencyListGraph();
    AdjacencyListGraph(const AdjacencyListGraph &graph);
    AdjacencyListGraph &operator=(const AdjacencyListGraph &graph);
    void add_edge(unsigned int v, unsigned int w);
    bool is_edge(unsigned int v, unsigned int w) const;
    void remove_edge(unsigned int v, unsigned int w);
    unsigned int number_of_vertices() const;
    std::vector<unsigned int> get_neighbours(unsigned int vertex) const;
    unsigned int get_degree(unsigned int vertex) const;
};
