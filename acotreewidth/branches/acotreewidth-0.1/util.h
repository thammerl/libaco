#include <cstdlib>
#include <ctime>
#include <exception>

class FileNotFoundException : public std::exception {
  private:
    const char *filepath_;
  public: 
    FileNotFoundException(const char *filepath);
    const char *what() const throw();
};

namespace Util {
  unsigned int random_number(unsigned int range=RAND_MAX);
}
