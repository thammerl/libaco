#include "util.h"

FileNotFoundException::FileNotFoundException(const char *filepath) {
  filepath_ = filepath;
}

const char *FileNotFoundException::what() const throw() {
  return filepath_;
}

static bool seeded = false;
unsigned int Util::random_number(unsigned int range) {
  if(!seeded) {
    srand(time(0));
    seeded = true;
  }
  return (rand() % range);
}
