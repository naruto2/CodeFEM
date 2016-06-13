#include <iostream>


main() {
  double x = 0.0;
  double y = 0.0;
  double h = 0.5;

  for ( x = 0.0; x <= 1.0; x += h )
    for ( y = 0.0; y <= 1.0; y += h ) {
      if ( x == 0.0 && y == 0.0 ) { printf("%f %f v0\n", x, y ); continue; }
      if ( y == 0.0 && x == 1.0 ) { printf("%f %f v1\n", x, y ); continue; }
      if ( y == 1.0 && x == 1.0 ) { printf("%f %f v2\n", x, y ); continue; }
      if ( y == 1.0 && x == 0.0 ) { printf("%f %f v3\n", x, y ); continue; }

      if ( y == 0.0 ) { printf("%f %f e0\n", x, y ); continue; }
      if ( x == 1.0 ) { printf("%f %f e1\n", x, y ); continue; }
      if ( y == 1.0 ) { printf("%f %f e2\n", x, y ); continue; }
      if ( x == 0.0 ) { printf("%f %f e3\n", x, y ); continue; }

      printf("%f %f\n", x, y ); 
   }

  printf("%f %f\n", h/2.0, h/2.0 );
  printf("%f %f\n", 1.0-h/2.0,  1.0-h/2.0 );
}
