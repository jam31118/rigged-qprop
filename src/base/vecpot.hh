#include <iostream>

#include <parameter.hh>

#include "potentials.hh"

using std::endl;
using std::cout;

struct vecpot_param {
  double omega, E0, num_cycles, phase_pi;
};




void print_vecpot_param(struct vecpot_param vparam) {
    cout << "omega: " << vparam.omega << endl;
    cout << "E0: " << vparam.E0 << endl;
    cout << "num-cycles: " << vparam.num_cycles << endl;
    cout << "phase / pi " << vparam.phase_pi << endl;
}

bool if_exist_get_vecpot_param ( parameterListe& para_prop, char direction, 
    int vecpot_index, struct vecpot_param *vp) {


  // construct suffix for each parameter name in a form: '-i-j'
  string suffix = string("-") + to_string(direction) 
    + string("-") + to_string(vecpot_index);

  cout << "[ LOG ] Processing vecpot param with suffix: " << suffix << endl;

  bool vecpot_set_exist = true;

  try {
    vp->omega = para_prop.getDouble(string("omega")+suffix);
    vp->E0 = para_prop.getDouble(string("max-electric-field")+suffix);
    vp->num_cycles = para_prop.getDouble(string("num-cycles")+suffix);
    vp->phase_pi = para_prop.getDouble(string("phase-pi")+suffix);
  } catch (std::exception&) {
    vecpot_set_exist = false;
  }
  
  if (vecpot_set_exist) { print_vecpot_param(*vp); }

  return vecpot_set_exist;
}

int get_number_of_vecpot_param_set(parameterListe& para_prop, char direction) {
  cout << "[ LOG ] Getting number of vecpot param set\n";
  int number_of_vecpot_param_set = 0;
  struct vecpot_param buf;
  int vecpot_param_set_index = 1;
  while (if_exist_get_vecpot_param(para_prop, direction, vecpot_param_set_index, &buf)) {
    number_of_vecpot_param_set++;
    vecpot_param_set_index++;
  }
  return number_of_vecpot_param_set;
}

bool vecpot_param_name_is_in_original_form(parameterListe& para) {
  bool is_in_original_form = true;
  try { para.getDouble("omega"); }
  catch (std::exception&) { is_in_original_form = false; }
  return is_in_original_form;
}

int vecpot_param_with_original_param_34 ( parameterListe& para, struct vecpot_param *vp ) {
  
  try {
    vp->omega = para.getDouble("omega");
    vp->E0 = para.getDouble("max-electric-field");
    vp->num_cycles = para.getDouble("num-cycles");
  } catch (std::exception&) {
    std::cerr << "[ERROR] Some vecpot potential\
      parameters in original form is missing." << endl;
    std::cerr << "[ERROR] Check all 'omega', 'max-electric-field', 'num-cycles' "
      << "are in corresponding parameter file." << endl;
    return 1;
  }

  try { vp->phase_pi = para.getDouble("phase") / M_PI; }
  catch (std::exception&) { vp->phase_pi = 0.0; }  // default value
  
  return 0;
}

int construct_vecpot_with_original_param_34 ( parameterListe& para,
    vecpot& vecpot_x, vecpot& vecpot_y, vecpot& vecpot_z ) {

  struct vecpot_param vp;
  if (vecpot_param_with_original_param_34(para, &vp)) {
    std::cerr << "[ERROR] Failed to construct vecpot with param set in original form.\n";
    return 1; }

  vecpot_x = vecpot(vp.omega, 1.0, 0.0, 0.0);
  vecpot_y = vecpot(vp.omega, 1.0, 0.0, 0.0);
  vecpot_z = vecpot(vp.omega, vp.num_cycles, vp.E0, vp.phase_pi * M_PI);

  return 0;
}


int if_exist_get_superposed_vecpot(parameterListe& para, char direction, 
    superposed_vecpot& svp, long num_of_vecpot) {

  vecpot *vp_arr = new vecpot[num_of_vecpot];
  bool vparam_exist = false;
  struct vecpot_param vparam;
  long vecpot_index;
  long vecpot_index_in_param_name;
  for (vecpot_index=0; vecpot_index < num_of_vecpot; vecpot_index++) {
    vecpot_index_in_param_name = vecpot_index + 1;
    vparam_exist = if_exist_get_vecpot_param(para, direction, 
        vecpot_index_in_param_name, &vparam);
    if ( ! vparam_exist ) {
      std::cerr << "[ERROR] vecpot param set with direction " << direction 
        << " and index " << vecpot_index_in_param_name << endl;
      return 1;
    }
    vp_arr[vecpot_index] = vecpot(vparam.omega, vparam.num_cycles, vparam.E0, vparam.phase_pi * M_PI);
  }
  svp = superposed_vecpot(vp_arr, num_of_vecpot);
  return 0;
}

void print_vecpot(vecpot& vp, const char *tag) {
  cout << "[ LOG ] ---------------------------------" << endl; 
  cout << "[ LOG ] vecpot with tag: " << tag << endl; 
  cout << "[ LOG ] - omega = " << vp.get_omega() << endl;
  cout << "[ LOG ] - E0 = " << vp.get_E0() << endl;
  cout << "[ LOG ] - num_cycles " << vp.get_num_cycles() << endl;
  cout << "[ LOG ] - phase " << vp.get_phase() << endl;
  cout << "[ LOG ] ---------------------------------" << endl; 
}

void print_superposed_vecpot(superposed_vecpot& svp, const char *tag) {
  long vp_index;
  vecpot *vp_arr = svp.get_vecpot_arr();
  cout << "[ LOG ] **********************************" << endl;
  for (vp_index = 0; vp_index < svp.get_num_of_vecpot(); vp_index++) {
    print_vecpot(vp_arr[vp_index], tag);
  }
  cout << "[ LOG ] **********************************" << endl;
}



int construct_superposed_vecpot_at_direction(parameterListe& para,
    char direction, superposed_vecpot& svp) {
  int number_of_vecpot_param_set = get_number_of_vecpot_param_set(para, direction);
  int return_code = if_exist_get_superposed_vecpot(para, 
      direction, svp, number_of_vecpot_param_set);
  if ( return_code != 0 ) { 
    std::cerr << "[ERROR] Failed to get vecpot for direction "
        << direction << " and number of vecpot " 
        << number_of_vecpot_param_set << endl;
    return 1;
  }
  return 0;
}


//int construct_vecpot ( long dim, parameterListe& para_prop, 
//    vecpot& vecpot_x, vecpot& vecpot_y, vecpot& vecpot_z ) {
int construct_vecpot ( long dim, parameterListe& para_prop, 
    superposed_vecpot& vecpot_x, superposed_vecpot& vecpot_y, superposed_vecpot& vecpot_z ) {

  if (dim==34) {
    //char direction = 'z';
    cout << "[ LOG ] in constructing clause\n";
    
    if ( vecpot_param_name_is_in_original_form(para_prop) ) {
      cout << "[ LOG ] vecpot param names are in original form\n";
      vecpot vpx, vpy, vpz;
      construct_vecpot_with_original_param_34(para_prop, vpx, vpy, vpz);
      print_vecpot(vpz, "vpz");
      vecpot *vpx_arr = new vecpot[1]; vpx_arr[0] = vpx;
      vecpot *vpy_arr = new vecpot[1]; vpy_arr[0] = vpy;
      vecpot *vpz_arr = new vecpot[1]; vpz_arr[0] = vpz;
      vecpot_x = superposed_vecpot(vpx_arr, 1);
      vecpot_y = superposed_vecpot(vpy_arr, 1);
      vecpot_z = superposed_vecpot(vpz_arr, 1);
    }

    else {
      cout << "[ LOG ] vecpot param names are in '-i-j' form\n";

      struct vecpot_param vecpot_param_single;

      vecpot *vpx = new vecpot(0.1, 1.0, 0.0, 0.0);
      vecpot_x = superposed_vecpot(vpx, 1);
      vecpot *vpy = new vecpot(0.1, 1.0, 0.0, 0.0);
      vecpot_y = superposed_vecpot(vpy, 1);
      
      // get vecpot_z
      if ( construct_superposed_vecpot_at_direction(
            para_prop, 'z', vecpot_z) ) { return 1; }

//      int number_of_vecpot_param_set_z = get_number_of_vecpot_param_set(para_prop, 'z');
//      cout << "[ LOG ] number of vecpot_param set in z : "
//        << number_of_vecpot_param_set_z << endl;
//      int return_code = 1;
//      return_code = if_exist_get_superposed_vecpot(para_prop, direction, 
//          vecpot_z,number_of_vecpot_param_set_z);
//      if ( return_code != 0 ) { std::cerr << "[ERROR] Failed to get vecpot for direction "
//        << direction << " and number of vecpot " << number_of_vecpot_param_set_z << endl;
//      }
      
    }

  }

  else if (dim == 44) {
    // get vecpot_x
    if ( construct_superposed_vecpot_at_direction(
          para_prop, 'x', vecpot_x) ) { return 1; }
    // get vecpot_y
    if ( construct_superposed_vecpot_at_direction(
          para_prop, 'y', vecpot_y) ) { return 1; }
    // get vecpot_z
    vecpot *vpz = new vecpot(0.1, 1.0, 0.0, 0.0);
    vecpot_z = superposed_vecpot(vpz, 1);
    
  } else {
    std::cerr << "[ERROR] Unexpected propagation mode (dimension): "
      << dim << endl;
  }
}

