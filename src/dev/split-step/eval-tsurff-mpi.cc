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
  const double n_cyc=para_prop.getDouble("num-cycles");
  vecpot vecpot_x(para_prop.getDouble("omega"), n_cyc, 0.0, 0.0);
  vecpot vecpot_y(para_prop.getDouble("omega"), n_cyc, 0.0, 0.0);
  vecpot vecpot_z(para_prop.getDouble("omega"), n_cyc, para_prop.getDouble("max-electric-field"), 0.0);
  tsurffSpectrum<vecpot, vecpot, vecpot> tsurff(para_ini, para_prop, para_tsurff, vecpot_x, vecpot_y, vecpot_z);
  tsurff.time_integration();
  // tsurff.print_int_dt_psi();
  // tsurff.print_wigner(0);
  // tsurff.print_Plms();
  // tsurff.print_bessel();
  tsurff.polar_spectrum();
  tsurff.print_partial_amplitudes();
#ifdef HAVE_MPI
  int ierrfin = MPI_Finalize();
  cout << "mpi finalize returns " << ierrfin << endl;
#endif
  return 0;
};
