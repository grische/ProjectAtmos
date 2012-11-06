#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SIGMA 5.67E-8   /* W / (m^2 K^4) */
#define c1 3.74E-16 /* W / m^2 */
#define c2 1.44E-2 /* m K */


double boltzmann (const double T) {
  return SIGMA * pow(T, 4);
}

double boltzmann_overpi (const double T) {
  return boltzmann(T)/M_PI;
}

double planck (const double T, const double lamda) {
  return c1/(pow(lamda,5)*(exp(c2/(lamda*T))-1.0));
}

double TEmission (double tau, const double Lbelow, const double Tlyr) {
    double TEmission = Lbelow * exp(-tau) + boltzmann_overpi(Tlyr)*(1.0-exp(-tau));
    //printf ("lup(): L = %f, tau=%f, Lbelow=%f, Tlyr=%f\n", tau, lup, Lbelow, Tlyr);
    return TEmission;
}

int schwarzschild(const double* deltatau, const double *T, const int nlev, const double Ts, double *edn, double *eup) {
  
  const double dmu = 0.01;
  //  const double dtau = tau/(nlev-1);

  double *lup=calloc(nlev, sizeof(double));
  double *ldn=calloc(nlev, sizeof(double));

  int ilev;
  double mu;
  
  /* Remove next for-loop for optimization, it is not needed */
  for(ilev=0; ilev < nlev; ilev++) {
    eup[ilev] = 0;
    edn[ilev] = 0;
  }

  /*
   * UP
   * ilev = 1, 2, ..., nlev-1 (=nlyr)
   * ilyr = 1, 2, ..., nlyr-1
   *
   * Wegen up: ilev -1 = ilyr
   *
   * lup[ilev] = lup[ilev+1] exp(-dtau)  +  B(T[ilyr+1]) ( 1 - exp(-dtau) )
   * quasi  lup[ilyr] = lup[ilyr+1] exp(-dtau)  +  B(T[ilyr]) (1 - exp(-dtau))
   *
   */
  //Calculation of eup from surface
  lup[nlev-1] = boltzmann_overpi(Ts);
  eup[nlev-1] = lup[nlev-1]*M_PI;

  for (mu=dmu/2; mu <= 1; mu += dmu) {
    for (ilev=nlev-2;ilev>=0; ilev--) {
    // integration ueber raumwinkel -> 2*PI und Integration ueber mu
      lup[ilev]=TEmission(deltatau[ilev]/mu, lup[ilev+1], T[ilev]);
      eup[ilev] += lup[ilev]*mu*dmu*2.0*M_PI;
    }
  }
  
  /*
   * DOWN
   * 
   * ldn[ilev] = ldn[ilev-1] exp(-dtau) + B(T[ilyr-1]) ( 1-exp(-dtau))
   * quasi ldn[ilyr] = ldn[ilyr-1] exp(-dtau) + B(T[ilyr]) (1-exp(-dtau))
   *
   */
  ldn[0] = 0;
  edn[0] = 0;
  
  for (mu=dmu/2; mu <= 1; mu += dmu) {
    for (ilev=1; ilev<nlev; ilev++) {
      
      ldn[ilev]= TEmission(deltatau[ilev-1]/mu, ldn[ilev-1], T[ilev-1]); 
      edn[ilev] += ldn[ilev]*mu*dmu*2.0*M_PI;
    }
  }

  // output of eup, L, edn
  // for (ilev=0; ilev < nlev; ilev++) {
  //  printf (  "--%d-- ilev, eup = %f, Lbelow = %f,edn = %f \t--%d--\n" ,ilev,  Eup[ilev], lup[ilev], edn[ilev], ilev);
  // }

  return EXIT_SUCCESS;
}
