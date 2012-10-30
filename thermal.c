#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double planck (double T) {
  double B=(5.67E-8)*pow(T,4)/M_PI;
  return B;
}

double Lup (double tau, const double Lbelow, const double Tlyr) {
  double Lup = Lbelow * exp(-tau) + planck(Tlyr)*(1.0-exp(-tau));
  //printf ("Lup(): L = %f, tau=%f, Lbelow=%f, Tlyr=%f\n", tau, Lup, Lbelow, Tlyr);
  return Lup;
}

int main() {
  const double Tsurf = 255;
  const double T[] = {108.801737
		,148.944353
		,172.361770
		,189.763422
		,203.897599
		,215.935755
		,226.497577
		,235.954833
		,244.549807
		,252.450005
  };

  const int nlyr = 10;
  const int nlev = nlyr+1;
  const double dmu = 0.01;
  const double dtau =1.000/nlyr;

  double *L=calloc(nlev, sizeof(double));
  double *Eup=calloc(nlyr, sizeof(double));

  int ilyr;
  double mu;
  
  for(ilyr=0; ilyr < nlyr; ilyr++) {
    Eup[ilyr] = 0;
  }

  //surface L
  L[nlev-1] = planck(Tsurf);

  /*
   * UP
   * ilev = 1, 2, ..., nlev-1 (=nlyr)
   * ilyr = 1, 2, ..., nlyr-1
   *
   * Wegen up: ilev -1 = ilyr
   *
   * L[ilev] = L[ilev+1] exp(-dtau)  +  B(T[ilyr+1]) ( 1 - exp(-dtau) )
   * quasi  L[ilyr] = L[ilyr+1] exp(-dtau)  +  B(T[ilyr]) (1 - exp(-dtau))
   *
   */
  for (ilyr=nlyr-1;ilyr>0; ilyr--) {
    // integration ueber raumwinkel -> 2*PI und Integration ueber mu
    for (mu=dmu; mu <= 1; mu += dmu) {
      L[ilyr]=Lup(dtau/mu, L[ilyr+1], T[ilyr]);
      Eup[ilyr] += L[ilyr]*mu*dmu;
    }
    Eup[ilyr] *= 2*M_PI;
    printf ("ilyr = %d, Eup = %f, Lbelow = %f\n" ,ilyr,  Eup[ilyr], L[ilyr]);
  }

  //printf("Die, die, die! %d %f\n", nlev+nlyr,L[nlev+nlyr+500]); 
  return 0;
}
