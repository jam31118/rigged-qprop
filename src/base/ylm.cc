#include<ylm.h>

/*! \fn complex Ylm(long l, long m, double theta, double phi)
    evaluation of Ylm from Racah
    Expansion yields polynomial in sin(theta)
    using 5.2.(17) from Varshalovich et al. (1988).
*/

cplxd ylm2(long l, long m, double theta, double phi)
{
  cplxd ylm;
  long L, s;
  double fact;

  if (l<0) {fprintf(stderr, "err: Ylm2: l<0\n"); exit(-69);}
  if (l<fabs(m)) {fprintf(stderr, "err: Ylm2: l<|m|\n"); exit(-70);}
      
  if  (!((l-m)&1)) {L=l; fact=1.0;} else {L=l-1; fact=cos(theta);};

  if (sin(theta) != 0.0)
  {
    for(s=(long)(fabs(m)); s<=L; s=s+2)
    {
       ylm = ylm + pow(-1.0,(s+m)/2.0) * factorial(l+s) / 
	    (factorial((s+m)/2)*factorial((s-m)/2)*pow(2.0,double(s))) *
            (factorial((L+m)/2)*factorial((L-m)/2)) / 
	    (factorial((L+s)/2)*factorial((L-s)/2)) *
             pow(sin(theta),double(s-fabs(m)));
    }

    ylm = ylm * sqrt((l+l+1)/(4*M_PI*factorial(l+m)*factorial(l-m))) *  
    pow(sin(theta),double(fabs(m))) * fact * exp(cplxd(0.0, m*phi));
  }
  else if (m == 0) 
  {
    ylm = sqrt((l+l+1)/M_PI) * fact / 2;
  } 
  else 
  {
    ylm = 0.0; 
  }

  return(ylm);
}




cplxd ylm(long l, long m, double theta, double phi)
{
  double ylm;  
  double ylm0,ylm1;
  long lindex;

  if (m != 0) {fprintf(stderr, "err: Ylm: not implemented yet for m <> 0 \n"); exit(-71);}
  if (l<0) {fprintf(stderr, "err: Ylm: l<0\n"); exit(-69);}
  if (l<fabs(m)) {fprintf(stderr, "err: Ylm: l<|m|\n"); exit(-70);}
  

  if (l==0) 
    {
      ylm=0.5/sqrt(M_PI);
    }
  else
    {
      if (l==1)
	{
	  ylm=0.5*sqrt(3.0/M_PI)*cos(theta);
	}
      else
	{
	  ylm0=0.5/sqrt(M_PI);
	  ylm1=0.5*sqrt(3.0/M_PI)*cos(theta);
	  lindex=2;
	  do 
	    {
	      ylm=sqrt((2.0*lindex-1.0)*(2.0*lindex+1.0))/(1.0*lindex)*cos(theta)*ylm1 - (lindex-1.0)/(1.0*lindex)*sqrt((2.0*lindex+1.0)/(2.0*lindex-3.0))*ylm0;
	      ylm0=ylm1;
	      ylm1=ylm;
	      lindex++;
	    }
	  while (lindex<=l);
	};
    };



  return cplxd(ylm,0.0);
}
