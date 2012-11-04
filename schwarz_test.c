#include <stdio.h>
#include <stdlib.h>
#include "schwarz.h"

int main() {

  const double Tsurf = 255;
  const double T[] = {100
                ,100
                ,100
                ,100
                ,100
                ,100
                ,100
                ,100
                ,100
                ,100
  };

  const int nlyr = 10;
  const int nlev = nlyr+1;
 
  const double dmu = 0.001;

  double *eup=calloc(nlev, sizeof(double));
  double *edn=calloc(nlev, sizeof(double));
  double *enet=calloc(nlyr, sizeof(double));

  int error = -1000;
  int ilev, ilyr;

  //schwarzschild(const double tau, const double *T, const int nlev, const double Ts, double *edn, double *eup);
  error = schwarzschild(1.0, T, nlev, Tsurf, edn, eup);

  /* E-netto */
  for(ilyr=0; ilyr<nlyr; ilyr++) {
    enet[ilyr] = eup[ilyr+1] + edn[ilyr] - eup[ilyr] - edn[ilyr+1];
  }

  if(error == EXIT_SUCCESS) {
    // output of eup, L, edn and enet
    for (ilev=0; ilev < nlev; ilev++) {
      printf (  "--%d-- ilev, eup = %f,edn = %f \t--%d--\n" ,ilev,  eup[ilev], edn[ilev], ilev);
      if(ilev < nlyr)
        printf ("  %d   ilyr, enetto = %f             \t  %d\n", ilev, enet[ilev], ilev);
    }
  } else {
    printf(" HEEEELP, ERROR %d\n", error);
  }

  return EXIT_SUCCESS;
}
