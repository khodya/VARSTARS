#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// ellipse parameters
double a = 1.0;
double b = 0.5;
double c = 0.5;


/* function: get_square
   input parameters:
   phi - rotation angle

*/
double get_square(double phi)
{
  // visible semi-axises of ellipse
  double va;
  double vb;

  // find major semi-axis
  // psi angle - radius-vector angle of maximum projection
  double psi = 0;
  if ( fabs( phi - M_PI / 2) > 1e-3  )
    {
      psi  = atan( 1.0L * b / a * tan( phi ));
      va = sqrt( pow( a * cos( psi ), 2) +
		 pow( b * sin( psi ), 2));
    }
  else
    {
     
      va = b;
    }

  // find minor semi-axis
  vb = b;

  //debug:
  printf("%lf %lf %lf %lf %lf\n",
	 phi * 180 / M_PI,
	 psi * 180 / M_PI,
	 b / a * tan( phi ),
	 va,
	 M_PI * va * vb);

  return M_PI * va * vb;
}



void run()
{
  int N = 90;
  double left = 0;
  double right = M_PI / 2;
  double step = (right - left) / N;
  double alpha;

  FILE *f;
  f = fopen("ellipse.txt", "w");
  for (int i=0; i < N + 1; i++)
    {
      fprintf(f, "%lf %lf\n",
	      alpha,
	      get_square(alpha));
      alpha += step;
    }
  
  fclose(f);
}

 void tests()
 {
   printf("circle square: %lf\n", M_PI * b * b);
 }

int main( int arc, const char **argv )
{ 
  //  run();
  tests();
  return 0;
}
