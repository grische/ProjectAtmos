#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <plplot/plplot.h>

double planck (double T) {
  double B=(5.67E-8)*pow(T,4)/M_PI;
  return B;
}

/*
void Lup (double *L, int nlyr, int ilyr, double dtau, const double *T, double Tsurf) {
  L[ilyr]=planck(Tsurf)*exp(-dtau*(nlyr-ilyr))+
    planck(T[ilyr])*(1.0-exp(-dtau*(nlyr-ilyr)));
  
  printf ("L[%d] = %f , T[%d] = %f, plank(T[%d]) = %f\n", ilyr, L[ilyr], ilyr, T[ilyr], ilyr, planck(T[ilyr]));  
}
*/

double Lup (double tau, const double T, double Tsurf) {
  double Lup = planck(Tsurf)*exp(-tau) + planck(T)*(1.0-exp(-tau));
  printf ("Lup(): tau=%f, L = %f\n", tau, Lup);
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
  const double dtau =0.00000001/nlyr;

  double *L=calloc(nlyr, sizeof(double));
  double *Eup=calloc(nlyr, sizeof(double));
  int ilyr;
  //  int jlyr;
  double mu;
  const double dmu = 0.01;
  Eup [0]=0 ;

  //  for (ilyr=0;ilyr<nlyr; ilyr++) {
    for (mu=dmu; mu <= 1; mu += dmu) {
      L[0]=Lup(dtau*(nlyr-0)/mu, T[0], Tsurf);
      Eup[0] += L[0]*mu*dmu;
    }
    Eup[0] *= 2*M_PI;
    //  }
    printf ("Eup = %f\n" , Eup[0]);
  


  return 0;
}
