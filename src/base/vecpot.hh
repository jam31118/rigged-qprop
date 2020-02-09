#ifndef _VECPOT_HH_
#define _VECPOT_HH_


#include <iostream>
#include <cstring>
#include <string>

#include <parameter.hh>

#include "potentials.hh"


struct vecpot_param {
  double omega, E0, num_cycles, phase_pi, start_time;

  vecpot_param(): start_time(0.0) {};

  vecpot_param(
      double omega_in, double E0_in, double num_cycles_in, double phase_pi_in, 
      double start_time_in=0.0):
    omega(omega_in), num_cycles(num_cycles_in), phase_pi(phase_pi_in),
    start_time(start_time_in) {};

};

const char valid_directions[] = "xyz";


void print_vecpot_param(struct vecpot_param vparam);

bool if_exist_get_vecpot_param ( parameterListe& para_prop, char direction, int vecpot_index, struct vecpot_param *vp);

int get_number_of_vecpot_param_set(parameterListe& para_prop, char direction);

bool vecpot_param_name_is_in_original_form(parameterListe& para);

int vecpot_param_with_original_param_34 ( parameterListe& para, struct vecpot_param *vp );

int construct_vecpot_with_original_param_34 ( parameterListe& para, vecpot& vecpot_x, vecpot& vecpot_y, vecpot& vecpot_z );

int if_exist_get_superposed_vecpot(parameterListe& para, char direction, superposed_vecpot& svp, long num_of_vecpot);

void print_vecpot(vecpot& vp, const char *tag);

void print_superposed_vecpot(superposed_vecpot& svp, const char *tag);

int construct_superposed_vecpot_at_direction(parameterListe& para, char direction, superposed_vecpot& svp);

int construct_vecpot ( long dim, parameterListe& para_prop, superposed_vecpot& vecpot_x, superposed_vecpot& vecpot_y, superposed_vecpot& vecpot_z );




#endif // _VECPOT_HH_
