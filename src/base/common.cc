#include "common.hh"

double get_grid_size(parameterListe para_ini, parameterListe para_prop, parameterListe para_tsurff) {

  double beyond_R_distance_temp;
  try { beyond_R_distance_temp = para_tsurff.getDouble("beyond-R"); }
  catch (std::exception&) { beyond_R_distance_temp = 0.0; }
  const double beyond_R_distance = beyond_R_distance_temp;
  double grid_size = para_tsurff.getDouble("R-tsurff") + beyond_R_distance + para_prop.getDouble("imag-width");

  return grid_size;
}
