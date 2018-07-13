#include "common.hh"

double get_grid_size(parameterListe para_ini, parameterListe para_prop, parameterListe para_tsurff) {

  double beyond_R_distance_temp;
  try { beyond_R_distance_temp = para_tsurff.getDouble("beyond-R"); }
  catch (std::exception&) { beyond_R_distance_temp = 0.0; }
  const double beyond_R_distance = beyond_R_distance_temp;
  double grid_size = para_tsurff.getDouble("R-tsurff") + beyond_R_distance + para_prop.getDouble("imag-width");

  return grid_size;
}


double get_effpot_alpha(parameterListe para_ini) {
  
  double alpha;
  try { alpha = para_ini.getDouble("effpot-alpha"); }
  catch (std::exception&) { alpha = 0.0; } // default value - if alpha = 0.0 means no effective potential but just become coulumb potential.

  return alpha;

}
