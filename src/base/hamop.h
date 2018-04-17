#ifndef hamop_h
#define hamop_h hamop_h
#include<functional>
#include<complex>
#include<iostream>
#include<cmath>
#include<cstdlib>
#include<grid.h>

// using namespace std;

class wavefunction;
class fluid;

class hamop
{
 public:
  hamop(grid g, 
	std::function<double(double, int)> fpx,
	std::function<double(double, int)> fpy,
	std::function<double(double, int)> fpz,
	std::function<double(double, double, double, double, int)> fpsx,
	std::function<double(double, double, double, double, int)> fpsy,
	std::function<double(double, double, double, double, int)> fpsz,
	std::function<double(long, long, long, double, grid)> fpi,
	std::function<double(double, int)> fpf
       );

  hamop(void);

 int init(grid g,
          std::function<double(double, int)> fpx,
	  std::function<double(double, int)> fpy,
	  std::function<double(double, int)> fpz,
	  std::function<double(double, double, double, double, int)> fpsx,
	  std::function<double(double, double, double, double, int)> fpsy,
	  std::function<double(double, double, double, double, int)> fpsz,
	  std::function<double(long, long, long, double, grid)> fpi,
	  std::function<double(double, int)> fpf
	  );

 int init(std::function<double(double, int)> fpx,
	  std::function<double(double, int)> fpy,
	  std::function<double(double, int)> fpz,
	  std::function<double(double, double, double, double, int)> fpsx,
	  std::function<double(double, double, double, double, int)> fpsy,
	  std::function<double(double, double, double, double, int)> fpsz,
	  std::function<double(long, long, long, double, grid)> fpi,
	  std::function<double(double, int)> fpf
	  );

  double vecpot_x(double time, int me);  
  double vecpot_y(double time, int me);  
  double vecpot_z(double time, int me);
  double scalarpot(double x, double y, double z, double time, int me);
  double scalarpotx(double x, double y, double z, double time, int me);
  double scalarpoty(double x, double y, double z, double time, int me);
  double scalarpotz(double x, double y, double z, double time, int me);
  double imagpot(long xindex, long yindex, long zindex, double time, grid g);
  double field(double time, int me);

 private:

  double delta_x, delta_y, delta_z;

  std::function<double(double, int)> hamopvecpotx;
  std::function<double(double, int)> hamopvecpoty;
  std::function<double(double, int)> hamopvecpotz;
  std::function<double(double, double, double, double, int)> hamopscalarpotx;
  std::function<double(double, double, double, double, int)> hamopscalarpoty;
  std::function<double(double, double, double, double, int)> hamopscalarpotz;
  std::function<double(long, long, long, double, grid)> hamopimagpot;
  std::function<double(double, int)> hamopfield;
};



#endif // hamop_h
