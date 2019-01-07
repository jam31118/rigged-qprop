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


std::string get_default_parameter_file_path(void) { 
  char *qprop_home_path_c_str = getenv(qprop_home_env_var_name);
  if ( qprop_home_path_c_str == NULL ) { 
    std::cerr << "[ERROR] " << qprop_home_env_var_name << " does not exist\n"; 
    exit(-1);
  }
  std::string qprop_home_path(qprop_home_path_c_str);
  std::string path_to_default_parameter_file = qprop_home_path + "/" + std::string(path_from_qprop_home_to_default_param_file);
//  std::cout << "[ LOG ] default_param_file: " << path_to_default_parameter_file << std::endl;
  if ( access( path_to_default_parameter_file.c_str(), F_OK ) != -1 ) {
    return path_to_default_parameter_file;
  } else {
    std::cerr << "[ERROR] The given file doesn't exist: " << path_to_default_parameter_file << std::endl;
    exit(-1);
  }
}


parameterListe get_default_parameter_list_object(void) {
//  std::string path_to_default_parameter_file = get_default_parameter_file_path();
  parameterListe default_parameter_list(get_default_parameter_file_path());
  return default_parameter_list;
}


