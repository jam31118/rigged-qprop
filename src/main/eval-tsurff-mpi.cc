#include <memory>
#include <complex>
#include <string>
#include <iostream>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

typedef std::complex<double> cplxd;
typedef std::unique_ptr<cplxd[]> cplxd_ptr;

#include <tsurffSpectrum.hh>
#include <potentials.hh>

#include <vecpot.hh>

using std::cout;
using std::string;

int main(int argc, char **argv) {

#ifdef HAVE_MPI
  int ierr = MPI_Init(&argc, &argv);
  cout << "mpi init returns " << ierr << endl;
#endif

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

#ifdef HAVE_MPI
  int ierrfin = MPI_Finalize();
  cout << "mpi finalize returns " << ierrfin << endl;
#endif

  return 0;
};
