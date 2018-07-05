#include "print-imagpot.hh"

int print_imagpot(int argc, char **argv) {
  
  parameterListe para_ini("initial.param");
  parameterListe para_prop("propagate.param");
  parameterListe para_tsurff("tsurff.param");
  
  // Add a parameter which determines radial distance between R-tsurff and imag-pot
  double beyond_R_distance_temp;
  try { beyond_R_distance_temp = para_tsurff.getDouble("beyond-R"); }
  catch (std::exception&) { beyond_R_distance_temp = 0.0; }
  const double beyond_R_distance = beyond_R_distance_temp;
  double grid_size = para_tsurff.getDouble("R-tsurff") + beyond_R_distance + para_prop.getDouble("imag-width");

  const double delta_r = para_ini.getDouble("delta-r");

  // Construct a grid object
  grid g_prop;
  g_prop.set_dim(para_ini.getLong("qprop-dim"));
  g_prop.set_ngps(long(grid_size/delta_r), para_prop.getLong("ell-grid-size"), 1); 
  g_prop.set_delt(delta_r);
  g_prop.set_offs(0, 0, 0);


  const long imag_potential_width=long(para_prop.getDouble("imag-width")/delta_r);
  imagpot imaginarypot(imag_potential_width);

  string imagpot_data_file_name("imagpot.dat");
  FILE *imagpot_data_file = fopen_with_check(imagpot_data_file_name, "w");

  for (long rho_index = 0; rho_index < g_prop.ngps_x(); rho_index++) {
    fprintf(imagpot_data_file, "%18ld %15.10e\n", rho_index, imaginarypot(rho_index, 0, 0, 0, g_prop));
  }

  fclose(imagpot_data_file);

  return 0;
}
