#include<hamop.h>
#include<functional>
#include<cmath>
#include<wavefunction.h>
#include<fluid.h>

// using std::isnan;

hamop::hamop(grid g, 
	     std::function<double(double, int)> fpx,
	     std::function<double(double, int)> fpy,
	     std::function<double(double, int)> fpz,
	     std::function<double(double, double, double, double, int)> fpsx,
	     std::function<double(double, double, double, double, int)> fpsy,
	     std::function<double(double, double, double, double, int)> fpsz,
	     std::function<double(long, long, long, double, grid)> fpi,
	     std::function<double(double, int)> fpf
	    )
{
  delta_x=g.delt_x();
  delta_y=g.delt_y();
  delta_z=g.delt_z();

  hamopvecpotx          = fpx;
  hamopvecpoty          = fpy;
  hamopvecpotz          = fpz;
  hamopscalarpotx       = fpsx;
  hamopscalarpoty       = fpsy;
  hamopscalarpotz       = fpsz;
  hamopimagpot          = fpi;
  hamopfield            = fpf;  
}


hamop::hamop(void)
{
  delta_x=0.0;
  delta_y=0.0;
  delta_z=0.0;

  // hamopvecpotx          = NULL;
  // hamopvecpoty          = NULL;
  // hamopvecpotz          = NULL;
  // hamopscalarpotx       = NULL;
  // hamopscalarpoty       = NULL;
  // hamopscalarpotz       = NULL;
  hamopimagpot          = NULL;
  hamopfield            = NULL;
}

//
//
//
int hamop::init(
		std::function<double(double, int)> fpx,
		std::function<double(double, int)> fpy,
		std::function<double(double, int)> fpz,
		std::function<double(double, double, double, double, int)> fpsx,
		std::function<double(double, double, double, double, int)> fpsy,
		std::function<double(double, double, double, double, int)> fpsz,
		std::function<double(long, long, long, double, grid)> fpi, 
		std::function<double(double, int)> fpf
		)
{
  delta_x=0.0;
  delta_y=0.0;
  delta_z=0.0;

  hamopvecpotx          = fpx;
  hamopvecpoty          = fpy;
  hamopvecpotz          = fpz;
  hamopscalarpotx       = fpsx;
  hamopscalarpoty       = fpsy;
  hamopscalarpotz       = fpsz;
  hamopimagpot          = fpi;
  hamopfield            = fpf;

}

int hamop::init(grid g,
		std::function<double(double, int)> fpx,
		std::function<double(double, int)> fpy,
		std::function<double(double, int)> fpz,
		std::function<double(double, double, double, double, int)> fpsx,
		std::function<double(double, double, double, double, int)> fpsy,
		std::function<double(double, double, double, double, int)> fpsz,
		std::function<double(long, long, long, double, grid)> fpi, 
		std::function<double(double, int)> fpf
		)
{
  delta_x=g.delt_x();
  delta_y=g.delt_y();
  delta_z=g.delt_z();

  hamopvecpotx          = fpx;
  hamopvecpoty          = fpy;
  hamopvecpotz          = fpz;
  hamopscalarpotx       = fpsx;
  hamopscalarpoty       = fpsy;
  hamopscalarpotz       = fpsz;
  hamopimagpot          = fpi;
  hamopfield            = fpf;

}


/*! \fn double hamop::vecpot_x(double time, int me)
    returns the vector potential \f$ A_x(t) \f$ .
    
    \param time  time ()
    \param me    a parameter to choose the intensity. 
    
    \remarks Parameter \a me is possibly obsolete. 

*/

double hamop::vecpot_x(double time, int me)
{
   double result;
   result = hamopvecpotx(time, me);
   if (std::isnan(result)) fprintf(stderr, "err: hamop::vecpot_x() is nan!");
   return result;
};


//
//
//
/*! \fn double hamop::vecpot_y(double time, int me)
    returns the vector potential \f$ A_y(t) \f$ .
    
    \param time  time ()
    \param me    a parameter to choose the intensity. 
    
    \remarks Parameter \a me is possibly obsolete. 

*/
double hamop::vecpot_y(double time, int me)
{
   double result;
   result = hamopvecpoty(time, me);
   if (std::isnan(result)) fprintf(stderr, "err: hamop::vecpot_y() is nan!");
   return result;
};

//
//
//
/*! \fn double hamop::vecpot_z(double time, int me)
    returns the vector potential \f$ A_z(t) \f$ .
    
    \param time  time ()
    \param me    a parameter to choose the intensity. 
    
    \remarks Parameter \a me is possibly obsolete. 

*/
double hamop::vecpot_z(double time, int me)
{
   double result;
   result = hamopvecpotz(time, me);
   if (std::isnan(result)) fprintf(stderr, "err: hamop::vecpot_z() is nan!");
   return result;
};

/*! \fn double hamop::scalarpot(double x, double y, double z,
                                double time, int me)
    returns the scalar potential \f$ U(r,t) \f$.
    				
    \param x     radial coordinate in the current implementation
    \param y     dummy
    \param z     dummy
    \param time  time 
    \param me    a parameter to choose the potential.
    
    \remarks Parameter \a me is possibly obsolete.
*/

double hamop::scalarpot(double x, double y, double z, double time, int me)
{
   double result;
   result = hamopscalarpotx(x, y, z, time, me)
              + hamopscalarpoty(x, y, z, time, me)
	      + hamopscalarpotz(x, y, z, time, me);
   if (std::isnan(result)) fprintf(stderr, "err: hamop::scalarpot() is nan!");
   return result;
}  

/*! \fn double hamop::scalarpotx(double x, double y, double z, 
               double time, int me)
	       
    returns the scalar potential \f$ U_x(r,t) \f$. The 
    function is possibly obsolete.
    				
    \param x     radial coordinate in the current implementation
    \param y     dummy
    \param z     dummy
    \param time  time 
    \param me    a parameter to choose the potential.
    
    \remarks Parameter \a me is possibly obsolete.
*/
double hamop::scalarpotx(double x, double y, double z, double time, int me)
{
  double result;
  result = hamopscalarpotx(x, y, z, time, me);
  if (std::isnan(result)) fprintf(stderr, "err: hamop::scalarpotx() is nan!");
  return result;
};

/*! \fn double hamop::scalarpoty(double x, double y, double z, 
               double time, int me)
	       
    returns the scalar potential \f$ U_y(r,t) \f$. The 
    function is possibly obsolete.
    				
    \param x     radial coordinate in the current implementation
    \param y     dummy
    \param z     dummy
    \param time  time 
    \param me    a parameter to choose the potential.
    
    \remarks Parameter \a me is possibly obsolete.
*/
double hamop::scalarpoty(double x, double y, double z, double time, int me)
{
  double result;
  result = hamopscalarpoty(x, y, z, time, me);
  if (std::isnan(result)) fprintf(stderr, "err: hamop::scalarpoty() is nan!");
  return result;

};

/*! \fn double hamop::scalarpotz(double x, double y, double z, 
               double time, int me)
	       
    returns the scalar potential \f$ U_z(r,t) \f$. The 
    function is possibly obsolete.
    				
    \param x     radial coordinate in the current implementation
    \param y     dummy
    \param z     dummy
    \param time  time 
    \param me    a parameter to choose the potential.
    
    \remarks Parameter \a me is possibly obsolete.
*/
double hamop::scalarpotz(double x, double y, double z, double time, int me)
{
  double result;
  result = hamopscalarpotz(x, y, z, time, me);
  if (std::isnan(result)) fprintf(stderr, "err: hamop::scalarpotz() is nan!");
  return result;

};


/*! \fn double hamop::imagpot(long xindex,long yindex,long zindex,
                              double time, grid g)
    returns the imaginary part of scalar potential \f$U(r)\f$.
    This potential have to effectively absorb the electron flux at the border.
    				
    \param xindex     index of radial distance in grid \a g.
    \param yindex     dummy
    \param zindex     dummy
    \param time       time 
    \param me         a parameter to choose the potential.
    
*/

double hamop::imagpot(long xindex, long yindex, long zindex, double time, grid g)
{
  double result;
  result = hamopimagpot(xindex,yindex,zindex,time,g);
  if (std::isnan(result)) fprintf(stderr, "err: hamop::imagpot() is nan!");
  return result;
}  

/*! \fn double hamop::field(double time, int me)
    returns the time dependent scalar potential \f$ U(t) \f$.
    
    \param time  time 
    \param me    a parameter to choose the intensity. 
    
    \remarks Parameter \a me is possibly obsolete. 
*/
double hamop::field(double time, int me)
{
  double result;
  result = hamopfield(time, me);
  if (std::isnan(result)) fprintf(stderr, "err: hamop::field() is nan!");
  return result;
}
