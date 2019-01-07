#ifndef _COMMON_HH_
#define _COMMON_HH_

#include <iostream>
#include <string>
#include <cstdlib>
#include <unistd.h>
#include "parameter.hh"

//#define CONFIG_CHAR_MAX_NUM 500

//
//const struct config {
//  char initial_wf_file_name[CONFIG_CHAR_MAX_NUM] = "";
//}

const char * const qprop_home_env_var_name = "QPROP_HOME";
const char * const path_from_qprop_home_to_default_param_file = "./src/base/default.param";

double get_grid_size(parameterListe para_ini, parameterListe para_prop, parameterListe para_tsurff);

double get_effpot_alpha(parameterListe para_ini);

std::string get_default_parameter_file_path(void);

parameterListe get_default_parameter_list_object(void);



#endif // _COMMON_HH_
