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
   
double r = 5;
double rad1 = 2;
double rad2 = 1.75;
double b1 = 1;
double b2 = 0.9;
double T = 1.6;
double m1 = 1.1;
double m2 = 0.9;

/* Statistics of stars configuration cases:
   0 - far away ( a > r1 + r2 )
   1 - contiguous ( a == r1 + r2 )
   2 - overlap, acute arc ( a > max radius && a < r1 + r2)
   3 - overlap, right arc ( a < max radius )
   4 - overlap, obtuse arc ( a < max radius )
   5 - full eclipse
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
double v_circle_square(double t, double incl)
{
  // square terms
  double s = 0;   // sum of circles
  s = pi * ( pow(rad1, 2) + pow(rad2, 2) );

  // current positional angle of a secondary star 
  double phi = 0;
  phi = 2 * pi / T * t;

  // visible distance
  double vd = vdist(incl, phi);
  printf("vd : %f\n", vd);
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
  if ( fabs( upper_edge - vd ) <= eps )
    {
      printf("1\n");
      return  s;
    }

  // stars are far away from each other
  if ( vd > upper_edge )
    {
      printf("2\n");
      return  s;
    }

  // full eclipse (case 5)
  if ( is_rad1_max &&  vd + rad2 - rad1 <= eps )
    {
      return pi * pow(rad1, 2);
    }
       

  // stars are overlapping ( includes subcases )
  double psi1 = 0;
  double psi2 = 0;
  double arg1 = 0;
  double arg2 = 0;

  printf("vd: %f\n", vd);
  arg1 = (pow(vd, 2) + pow(rad1, 2) - pow(rad2, 2)) / 2.0 / vd / rad1;
  arg2 = (pow(vd, 2) + pow(rad2, 2) - pow(rad1, 2)) / 2.0 / vd / rad2;

  psi1 = acos(arg1);
  psi2 = acos(arg2);

  printf("psi1: %f\n", psi1 * 180 / pi);
  printf("psi2: %f\n", psi2 * 180 / pi);

  // if rad1 is max
  if ( is_rad1_max )
    {
      // square terms


      
      double segm1 = pow(rad1, 2) * ( 2.0l * psi1 - sin(2 * psi1) ) / 2.0;

      double segm2 = pow(rad2, 2) * ( 2.0l * psi2 - sin(2 * psi2) ) / 2.0;

      double semisquare = pi * pow(rad2 ,2) / 2.0;

      double obtuse_arc = pow(rad2, 2) * ( 2.0l * psi2 ) / 2.0;

      double triangle  = pow(rad2, 2) * sin( -2.0l * psi2) / 2.0;

      printf("s: %f\n", s);
      printf("segm1: %f\n", segm1);
      printf("segm2: %f\n", segm2);
      printf("semisquare: %f\n", semisquare);
      printf("obtuse arc: %f\n", obtuse_arc);
      printf("triangle: %f\n", triangle);
      
	  
      // case 2 (acute arc)
      if ( psi2 - pi / 2.0l <= eps )
	return s - segm1 - segm2;


      // case 3 (right arc)
      if ( fabs(1.0l * psi2 - pi / 2.0l ) <= eps )
	{
	  
	  return s - semisquare - segm1;    
	}
      // case 4 (obtuse arc)
      return s - segm1 - obtuse_arc - triangle;	
    }

  return 0;
}

void tests()
{
  printf("tests:\n");

  printf("eps: %f\n", eps);

  printf("vdist:\n");

  printf("%f\n", vdist(0, pi / 2));
  printf("%f\n", vdist(pi / 3, pi / 2));
  printf("%f\n", vdist(pi / 4, pi / 2));
  printf("%f\n", vdist(1.332797, pi / 2));

  printf("v_circle_square:\n");

  printf("%f\n", v_circle_square(0, pi / 2) );
}

void run()
{
  int N = 100;
  double a = 0;
  double b = T;
  double h = (b - a) / (N - 1);
  double t = a;
  
  FILE *f;
  f = fopen("output.txt", "w");

  for (int i = 0; i < N; i++ )
    {
      fprintf(f, "%f %f\n", t, v_circle_square(t, pi / 2));
      t = t + h;
    }

	fclose(f);
}

int main( int arc, const char **argv )
{ 
  tests();
  run();
}
