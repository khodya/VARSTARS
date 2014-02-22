// This program is modeling a light curve for a HD52571 system
//
// Model:
// 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_NUM_CASES 10

double pi = 3.1415926; 
double eps = 1e-7;

/* model variables:
   r    radius of the orbit
   rad1 visible radius of a bigger star
   rad2 visible radius of a smaller star
   b1   decrease intensity factor for a bigger star
   b2   decrease intensity factor for a smaller star 
   T    orbital period */
   
double r = 1;
double rad1 = 1;
double rad2 = 0.9;
double b1 = 1;
double b2 = 0.9;
double T = 1.6;
double m1 = 1.1;
double m2 = 0.9;

/* Statistics of stars configuration cases 
   cases: 
   0 - far away
   1 - contiguous
   2 - overlap, acute angle
   3 - overlap, right angle
   4 - overlap, obtuse angle
*/
int cases_count[MAX_NUM_CASES];

/* This function will return a visible distance between stars
   consider an inclination angle.

   Input parameters:
   double dist -- real distance between stars (in orbit plane)
   double incl -- an inclination angle in radians
   double pa   -- positional angle in radians  */
double vdist(double incl, double pa)
{
  // only for (0, 0) coord center
  return  r * sqrt( cos(pa) * cos(pa) \
		    + cos(incl) * cos(incl) * sin(pa) * sin(pa) );
}

/* This function will return a visible square projection. */
/* Input parameters:
   t time
   a visible distance between stars 
   incl inclination angle
*/
double v_circle_square(double t, double a, double incl)
{
  // current positional angle of a secondary star 
  double phi = 0;
  phi = 2 * pi / T * t;

  // visible distance
  double vd = vdist(incl, phi);

  // upper edge
  double upper_edge = rad1 + rad2;

  // find min of radii
  double min = rad1;
  double max = rad2;
  int is_rad1_max = 0;

  if ( rad2 < rad1 )
    {
      min = rad2;
      max = rad1;
      is_rad1_max = 1;
    }

  // lower_edge
  double lower_edge = max - min;

  // stars are contiguous to each other
  if ( abs( upper_edge - vd ) <= eps )
    {
      return  pi * ( rad1 * rad1 + rad2 * rad2 );
    }

  // stars are far away from each other
  if ( abs( upper_edge - vd ) > eps )
    {
      return  pi * ( rad1 * rad1 + rad2 * rad2 );
    }

  // stars are overlapping ( includes subcases )


  return 0;
}

int main( int arc, const char **argv )
{ 
  printf("vdist tests:\n");
  printf("%f\n", vdist(0, pi / 2));
  printf("%f\n", vdist(pi / 3, pi / 2));
  printf("%f\n", vdist(pi / 4, pi / 2));

 
}
