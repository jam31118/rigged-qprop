#include "eval-tsurff.h"

using std::cout;
using std::string;

int eval_tsurff(int argc, char **argv) {

  parameterListe para_ini("initial.param");
  parameterListe para_prop("propagate.param");
  parameterListe para_tsurff("tsurff.param");

  long dimens = para_ini.getLong("qprop-dim");

  superposed_vecpot vecpot_x, vecpot_y, vecpot_z;
  construct_vecpot(dimens, para_prop, vecpot_x, vecpot_y, vecpot_z);

  if ( dimens == 34 ) {
    print_superposed_vecpot(vecpot_z, "superposed_vecpot_z");
    cout << "vecpot_z duration: " << vecpot_z.get_duration() << endl;
  } else if ( dimens == 44 ) {
    print_superposed_vecpot(vecpot_x, "superposed_vecpot_x");
    cout << "vecpot_x duration: " << vecpot_x.get_duration() << endl;
    print_superposed_vecpot(vecpot_y, "superposed_vecpot_y");
    cout << "vecpot_y duration: " << vecpot_y.get_duration() << endl;
  }

  tsurffSpectrum<superposed_vecpot, superposed_vecpot, superposed_vecpot>
    tsurff(para_ini, para_prop, para_tsurff, vecpot_x, vecpot_y, vecpot_z);

  tsurff.time_integration();
  tsurff.polar_spectrum();
  tsurff.print_partial_amplitudes();

  return 0;
};

