#ifndef wavefunction_h
#define wavefunction_h wavefunction_h
#include<iostream>
#include<assert.h>
#include<complex>
#include<fstream>
#include<cmath>
#include<cstdlib>
#include<string>

using namespace std;

void f(double, double[], double[]);

typedef std::complex<double> cplxd;

class fluid;
class hamop;
class grid;

class wavefunction
{
  public:
    wavefunction(long x=0) : wf_dim(x), start(new cplxd[x]) { }

    wavefunction(const wavefunction& v) 
      {
	wf_dim  = v.wf_dim;
	start = new cplxd[wf_dim];
	for (long i = 0; i < wf_dim; i++)
	  start[i] = v.start[i];
      }

    virtual ~wavefunction() { delete [] start;}

    //
    // Functions
    //
    void apply(const wavefunction &a, const wavefunction &b, 
	       const wavefunction &c, const wavefunction &psi);
    cplxd* begin() { return start;}
    const cplxd* begin() const { return start;}
    
    void calculate_staticpot(grid g, hamop hamil);

    wavefunction calculate_Theta(grid g, const fluid &degeneracies, const fluid &ms);
    fluid calculate_radial_density(grid g, const fluid &degeneracies);

    wavefunction conjugate();

    void dump_to_file(grid g, FILE* os, int stepwidth, double fact);
    void dump_xvector_to_file(long ngpsx,  FILE* os, int stepwidth);
    void dump_to_file_sh(grid g, FILE* os, int stepwidth);
    int  dump_to_file_sh(grid g, FILE* os, int stepwidth, int iv);
    void dump_expect_z(grid g, FILE* os, fluid &degeneracies, const fluid &ms);



    void embed_as_x(grid g, long yindex, long zindex, const wavefunction &v);
    void embed_as_x(grid g, long l, long m, long i, const wavefunction &v);
    void embed_as_xy(grid g, long zindex, const wavefunction &v);
    cplxd energy(double time, grid g, hamop hamil, int me,
      const wavefunction &staticpot, double charge);

    wavefunction extract_x(grid g, long other_one, long other_two);
    wavefunction extract_y(grid g, long other_one, long other_two);
    wavefunction extract_z(grid g, long other_one, long other_two);
    wavefunction extract_xy(grid g, long zindex);
    wavefunction expwf(cplxd a);
    wavefunction expwf(double a);
    cplxd expect_z(grid g, fluid &degeneracies, const fluid &ms);
    cplxd expect_z(grid g);
    cplxd expect_cycl_pol_plus(grid g);
    cplxd expect_cycl_pol_minus(grid g);
    cplxd accel_z(grid g, int m, double charge);

    cplxd expect_rr(grid g, fluid &degeneracies);

    int init(long isize);
    void init(grid g, int inittype, double width, fluid &ells);
    void init_rlm(grid g, int inittype, double w, fluid &ells, fluid &ms);
    void init(grid g, int inittype, FILE* file, int ooi);
    int init(grid g, FILE* file, int ooi, int iv);



    void regrid(grid g, grid g_small, const wavefunction &v);
    void regrid_and_rebin(grid g, grid g_small, const wavefunction &v);

    void select_single_orbital(grid g, grid g_small, int orb_of_interest,
      double ell, const wavefunction &v);


    double totalenergy_single_part(grid g, const wavefunction &orb_energs, const fluid &degeneracies);
   
    void nullify();
    double norm(grid g);
    void normalize(grid g, const fluid &ms);    
    void normalize(grid g);    
    wavefunction project(grid g, grid gorig, wavefunction &orig);
    cplxd project(grid g, grid gorig, wavefunction &orig, long zindex);
    
    int load(FILE*, int);

    cplxd* end()   { return start + wf_dim;}
    const cplxd* end()   const { return start + wf_dim;}

    void propagate(cplxd timestep, 
		   double time, 
		   grid g, 
		   hamop hamil, 
		   int me, 
		   const wavefunction &staticpot, 
		   const fluid &wf_one, 
		   const wavefunction &wf_two, 
		   const wavefunction &wf_three, 
		   const fluid &ms, 
		   double charge,
		   int propornot[]);
    
    void propagate(cplxd timestep, 
		   double time, 
		   grid g, 
		   hamop hamil, 
		   int me, 
		   const wavefunction &staticpot, 
		   int m, 
		   double charge);
      
    wavefunction orb_fieldenergies(double  time, grid g, hamop hamil, int me, const fluid &ms);

    wavefunction orbital_energies(double time, grid g, hamop hamil, int me,  const wavefunction &staticpot, double charge);
    cplxd orbital_energy(double time, grid g, hamop hamil, int me,  const wavefunction &staticpot, double charge, long orb_no);


    fluid orbital_norms(grid g);

    wavefunction orbital_hartrees(double time, grid g, int me, 
				  const fluid &wf_one);
    cplxd orbital_hartree(double time, grid g, int me, 
				  const fluid &wf_one, long orb_no);

    void realific();

    void solve(const wavefunction &a, const wavefunction &b, 
	       const wavefunction &c, const wavefunction &psi, long dimens);
    void solve_du(const wavefunction &a, const wavefunction &b, 
               const wavefunction &c, const wavefunction &psi, long dimens);
    void solve_toep(double a, double b, double b_upperleft, double b_lowerright,
               double c, const wavefunction &psi, long dimens);
    wavefunction sqrtwf();
    wavefunction sqrtrealwf();

    wavefunction mult_diag_with_diag(const wavefunction &a);

    wavefunction invert();
        long  wf_size() const {return wf_dim;}
    //
    // Operators
    //

    wavefunction& operator *= (double z); 
    wavefunction& operator *= (cplxd z); 
    cplxd&  operator[](long index)
    {
      //      assert(index >= 0  &&  index < wf_dim);
      return start[index];
    }

    const cplxd&  operator[](long index) const   
    {   
      //      assert(index >= 0  &&  index < wf_dim);
      return start[index];
    }

    wavefunction& operator=(const wavefunction&);

  private:
    long  wf_dim;
    cplxd   *start; 

    void do_muller_ell(cplxd timestep, 
		       double time, 
		       grid g, 
		       hamop hamil, 
		       const wavefunction &staticpot, 
		       int me, 
		       double charge,
		       int m);
      
    void do_muller_ellm(cplxd timestep, 
				  double time, 
				  grid g, 
				  hamop hamil, 
				  const wavefunction &staticpot, 
				  int me, 
				  double charge);

};

ostream& operator<<(ostream& os, const wavefunction& v);
//istream& operator>>(istream& is, wavefunction& v);
cplxd operator * (const wavefunction &v, const wavefunction &w );
wavefunction operator * (double z, const wavefunction &v);
wavefunction operator * (cplxd z, const wavefunction &v);
wavefunction operator * (const wavefunction &v, double z);
wavefunction operator / (const wavefunction &v, double z);
wavefunction operator * (const wavefunction &v, cplxd z);
wavefunction operator + (const wavefunction &v, const wavefunction &vv);
wavefunction operator - (const wavefunction &v, const wavefunction &vv);
wavefunction operator + (const wavefunction &v, const fluid &vv);
wavefunction operator + (const fluid &v, const wavefunction &vv);


#endif // wavefunction_h
