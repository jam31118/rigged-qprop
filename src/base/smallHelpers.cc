#include "smallHelpers.hh"

FILE* fopen_with_check(std::string fname, const char* attr) {
  FILE* file_ptr = fopen(fname.c_str(), attr);
  if (file_ptr==NULL) {
    cerr << fname << " can not be opened!" << endl;
    exit(-1);
  };
  return file_ptr;
};
