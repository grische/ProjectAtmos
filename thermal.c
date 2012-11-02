#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double planck (double T) {
  double B=(5.67E-8)*pow(T,4)/M_PI;
  return B;
}

double FLup (double tau, const double Lbelow, const double Tlyr) {
  double FLup = Lbelow * exp(-tau) + planck(Tlyr)*(1.0-exp(-tau));
  //printf ("Lup(): L = %f, tau=%f, Lbelow=%f, Tlyr=%f\n", tau, Lup, Lbelow, Tlyr);
  return FLup;
}

double FLdow (double tau, const double Labove, const double Tlyr) {
  double FLdow = Labove * exp(-tau) + planck(Tlyr)*(1-exp(-tau));
  return FLdow;


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
  const double dmu = 0.001;
  const double dtau =1.000/nlyr;

  double *Lup=calloc(nlev, sizeof(double));
  double *Ldow=calloc(nlev, sizeof(double));
  double *Eup=calloc(nlev, sizeof(double));
  double *Edow=calloc(nlev, sizeof(double));
  double *Enet=calloc(nlyr, sizeof(double));



  int ilyr;
  int ilev;
  double mu;
  
  /* Remove next for-loop for optimization, it is not needed */
  for(ilyr=0; ilyr < nlyr; ilyr++) {
    Eup[ilyr] = 0;
    Edow[ilyr] = 0;
  }

  /*
   * UP
   * ilev = 1, 2, ..., nlev-1 (=nlyr)
   * ilyr = 1, 2, ..., nlyr-1
   *
   * Wegen up: ilev -1 = ilyr
   *
   * Lup[ilev] = Lup[ilev+1] exp(-dtau)  +  B(T[ilyr+1]) ( 1 - exp(-dtau) )
   * quasi  Lup[ilyr] = Lup[ilyr+1] exp(-dtau)  +  B(T[ilyr]) (1 - exp(-dtau))
   *
   */
  //Calculation of Eup from surface
  Lup[nlev-1] = planck(Tsurf);
  Eup[nlev-1] = Lup[nlev-1]*M_PI;

  for (mu=dmu; mu <= 1; mu += dmu) {
    for (ilev=nlev-2;ilev>=0; ilev--) {
    // integration ueber raumwinkel -> 2*PI und Integration ueber mu
      Lup[ilev]=FLup(dtau/mu, Lup[ilev+1], T[ilev]);
      Eup[ilev] += Lup[ilev]*mu*dmu*2.0*M_PI;
    }
  }
  
  /*
   * DOWN
   * 
   * Ldow[ilev] = Ldow[ilev-1] exp(-dtau) + B(T[ilyr-1]) ( 1-exp(-dtau))
   * quasi Ldow[ilyr] = Ldow[ilyr-1] exp(-dtau) + B(T[ilyr]) (1-exp(-dtau))
   *
   */
  Ldow[0] = 0;
  Edow[0] = 0;
  
  for (mu=dmu; mu <= 1; mu += dmu) {
    for (ilev=1; ilev<nlev; ilev++) {
      
      Ldow[ilev]= FLdow(dtau/mu, Ldow[ilev-1], T[ilev-1]); 
      Edow[ilev] += Ldow[ilev]*mu*dmu*2.0*M_PI;
    }
  }

  // Net Energy per layer: Enet

  for (ilyr=0; ilyr<nlyr; ilyr++) {
    Enet[ilyr]=Eup[ilyr+1]+Edow[ilyr]-Eup[ilyr]-Edow[ilyr+1];
  }

  // output of Eup, L, Edow and Enet
  for (ilev=0; ilev < nlev; ilev++) {
    printf (  "--%d-- ilev, Eup = %f, Lbelow = %f,Edow = %f \t--%d--\n" ,ilev,  Eup[ilev], Lup[ilev], Edow[ilev], ilev);
    if(ilev < nlyr)
      printf ("  %d   ilyr, Enetto = %f                      \t\t\t  %d\n", ilev, Enet[ilev], ilev);
  }

  //printf("Die, die, die! %d %f\n", nlev+nlyr,L[nlev+nlyr+500]); 
  return 0;
}
