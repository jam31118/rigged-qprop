#include "printer.hh"

int print_imagpot(int argc, char **argv) {
  
  parameterListe para_ini("initial.param");
  parameterListe para_prop("propagate.param");
  parameterListe para_tsurff("tsurff.param");
  
  // Add a parameter which determines radial distance between R-tsurff and imag-pot
  double grid_size = get_grid_size(para_ini, para_prop, para_tsurff);
  const double delta_r = para_ini.getDouble("delta-r");

  // Construct a grid object
  grid g_prop;
  g_prop.set_dim(para_ini.getLong("qprop-dim"));
  g_prop.set_ngps(long(grid_size/delta_r), para_prop.getLong("ell-grid-size"), 1); 
  g_prop.set_delt(delta_r);
  g_prop.set_offs(0, 0, 0);


  const long imag_potential_width=long(para_prop.getDouble("imag-width")/delta_r);
  imagpot imaginarypot(imag_potential_width);

  string imagpot_data_file_name("imagpot.bin");
  FILE *imagpot_data_file = fopen_with_check(imagpot_data_file_name, "wb");

  double imagpot_value;
  for (long rho_index = 0; rho_index < g_prop.ngps_x(); rho_index++) {
    imagpot_value = imaginarypot(rho_index, 0, 0, 0, g_prop);
    fwrite(&imagpot_value, sizeof(double), 1, imagpot_data_file);
  }

  fclose(imagpot_data_file);

  return 0;
}


int print_scalarpot(int argc, char **argv) {

  parameterListe para_ini("initial.param");
  parameterListe para_prop("propagate.param");
  parameterListe para_tsurff("tsurff.param");
  
  // Add a parameter which determines radial distance between R-tsurff and imag-pot
  double grid_size = get_grid_size(para_ini, para_prop, para_tsurff);
  const double delta_r = para_ini.getDouble("delta-r");

  // Construct a grid object
  grid g_prop;
  g_prop.set_dim(para_ini.getLong("qprop-dim"));
  g_prop.set_ngps(long(grid_size/delta_r), para_prop.getLong("ell-grid-size"), 1); 
  g_prop.set_delt(delta_r);
  g_prop.set_offs(0, 0, 0);
  
  scalarpot atomic_potential(para_ini.getDouble("nuclear-charge"), para_ini.getDouble("pot-cutoff"), get_effpot_alpha(para_ini));

  string scalarpot_data_file_name("scalarpot.bin");
  FILE *scalarpot_data_file = fopen_with_check(scalarpot_data_file_name, "wb");

  long rho_index;
  double rho_value;
  double scalarpot_value;
  for (rho_index = 0; rho_index < g_prop.ngps_x(); rho_index++) {
    rho_value = (rho_index + 1) * delta_r;
    scalarpot_value = atomic_potential(rho_value, 0, 0, 0, 0);
    fwrite(&scalarpot_value, sizeof(double), 1, scalarpot_data_file);
  }

  fclose(scalarpot_data_file);

  return 0;
  
}


