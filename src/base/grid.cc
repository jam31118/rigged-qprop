#include <grid.h>

long num_of_basis(long qprop_dim, long ell_grid_size) {
  switch (qprop_dim) {
    case 34: return ell_grid_size;
    case 44: return ell_grid_size * ell_grid_size;
    default: return -1;
  }
}
