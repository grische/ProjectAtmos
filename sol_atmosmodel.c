#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <plplot/plplot.h>
#include "schwarz.h"
#include "mstrnx/ascii.h"
#include "mstrnx/rodents_solar.h"


//--------------------- Define constants ----------------------//

#define FALSE 0
#define TRUE 1

#define R 8.31447     /* J/(mol K) */
#define M 0.0289644   /* kg / mol */
#define g 9.8065      /* m/s^2 */

#define PSURF 1000     /* hPa */

#define TIME_MAX 3600 * 24 * 360 * 10 /* seconds */

//------------------- Implement sortfunc ----------------------//

static int sortfunc (const void *a, const void *b) {             /* Defining sortfunc */
  
  if(*((double*)a)<*((double*)b)) {
    return 1;
  }
  if(*((double*)a)>*((double*)b)) {
    return -1;
  }  
  return 0;
  
}


//----------------------Plot routine -------------------------//


void plotall(int nlyr, double* T, double* plyr, double* z, double* deltaTday) {
  /* Plot T against p */

  pladv(1);     /* select subpage 1  */
  plvsta();     /* standard viewport */
  plclear();    /* clear subpage     */
  plcol0 (15);  /* color black       */

  plwind( 0, 400, PSURF, 0 );  /* xmin, xmax, ymin, ymax */
  plbox( "bcnst", 100, 0, "bcnst", 150.0, 0 );
  pllab ("temperature [K]", "p [hPa]", "");  /* axis labels     */

  plcol0 (14);                         /* color blue  */
  plline (nlyr, T, plyr);  /* plot temperature profile  */

  plcol0 (15);                        /* color black */

  /* Plot T against z */

  pladv(3);     /* select subpage 1  */
  plvsta();     /* standard viewport */
  plclear();    /* clear subpage     */
  plcol0 (15);  /* color black       */

  plwind( 0, 400, 0, 40000 );  /* xmin, xmax, ymin, ymax */
  plbox( "bcnst", 100, 0, "bcnst", 5000.0, 0 );
  pllab ("temperature [K]", "z [m]", "");  /* axis labels     */

  plcol0 (14);                         /* color blue  */
  plline (nlyr, T, z);  /* plot temperature profile  */

  plcol0 (15);                        /* color black */

  /* Plot Heating rate against p */

  pladv(2);     /* select subpage 1  */
  plvsta();     /* standard viewport */
  plclear();    /* clear subpage     */
  plcol0 (15);  /* color black       */

  plwind( -20, 20, PSURF, 0 );  /* xmin, xmax, ymin, ymax */
  plbox( "bcnst", 2, 0, "bcnst", 150.0, 0 );
  pllab ("Heating Rate [T/day]", "p [hPa]", "");  /* axis labels     */

  plcol0 (12);                         /* color blue  */
  plline (nlyr, deltaTday, plyr);  /* plot temperature profile  */

  plcol0 (15);                        /* color black */


  /* Plot Heating rate against z */

  pladv(4);     /* select subpage 1  */
  plvsta();     /* standard viewport */
  plclear();    /* clear subpage     */
  plcol0 (15);  /* color black       */


  plwind( -20, 20, 0, 40000 );  /* xmin, xmax, ymin, ymax */
  plbox( "bcnst", 2, 0, "bcnst", 5000.0, 0 );
  pllab ("Heating Rate [T/day]", "z [m]", "");  /* axis labels     */

  plcol0 (12);                         /* color blue  */
  plline (nlyr, deltaTday, z);  /* plot temperature profile  */

  plcol0 (15);                        /* color black */
  
}



//------------------------- Defining variables and integers ----------------//

int main() {                                                    
  
  // http://cdiac.ornl.gov/pns/current_ghg.html
  const double c_co2 = 390;
  const double c_n2o = 0.32;
  const double c_co = 0.1;
  const double c_ch4 = 1.775;
  const double c_o2 = 20.9e+4;
  
  
  const double cp=1004;  /* J/kg K */
  const double deltat=3600.0*4.0;  /* s */
  const double Ra=287; /* J/kg K */

  int timesteps = 0;
  int ilev=0;
  int ilyr=0;
  const int nlyr=20;  /* Number of Layers */
  const int nlev=nlyr+1;
  int instabil=FALSE;

  int status;
  
  
  

  const double deltap=PSURF/nlyr;
  const double kappa=Ra/cp;
  
  double *p=calloc(nlev,sizeof(double));
  double *T=calloc(nlyr,sizeof(double));
  double *theta=calloc(nlev,sizeof(double));
  double *plyr=calloc(nlyr,sizeof(double));
  double *z=calloc(nlev,sizeof(double));
  double *zlyr=calloc(nlyr,sizeof(double));
  double *deltazlyr=calloc(nlyr,sizeof(double));
  double *deltazlev=calloc(nlev,sizeof(double));
  double *h2o=calloc(nlyr,sizeof(double));
  double *o3=calloc(nlyr,sizeof(double));
  
  double ***dtaumol;
  double **wgt;
  double *wavelength;
  int nbnd=0;
  int *nch;
  int iv;
  int iq;
  
  double *tmplev=calloc(nlev, sizeof(double)); /* temporary vector with length lev */
  double *tmplyr=calloc(nlyr, sizeof(double)); /* temporary vector with length lyr */
  
  double *eup=calloc(nlev, sizeof(double));
  double *edn=calloc(nlev, sizeof(double));
  double *euptmp=calloc(nlev, sizeof(double));
  double *edntmp=calloc(nlev, sizeof(double));
  double *edirtmp=calloc(nlev, sizeof(double));
  double *enet=calloc(nlyr, sizeof(double));

  double *deltaTday=calloc(nlyr, sizeof(double));
  double *deltaT=calloc(nlyr, sizeof(double));
  
  //------------- Define variables for Rodents---------------//
  
  const double gr_albedo=0.3; /* 1 */
  
  const int mu_counterlimit = (int)(24.0 * 3600 / deltat);
  double* const mu_0 = calloc(mu_counterlimit, sizeof(double));
  int mu_counter;
  for(mu_counter=0; mu_counter < mu_counterlimit; mu_counter++) {
    mu_0[mu_counter] = cos( 2.0*M_PI / mu_counterlimit * (mu_counter + mu_counter+1) / 2.0 );
    printf("mu_0[%d] = %e\n", mu_counter, mu_0[mu_counter]);
  }
  mu_counter = 0;
  
  double* S_0;
  double* omega_0=calloc(nlyr,sizeof(double));
  double* gassy=calloc(nlyr,sizeof(double));
  double* f=calloc(nlyr,sizeof(double));
  
  
  for(ilyr=0; ilyr<nlyr; ilyr++) {
   omega_0[ilyr] = 0.f;
   gassy[ilyr] = 0; //0.85f;
   f[ilyr] = 0; //0.8f;
  }
  
  {
    double *temp1, *temp2, *temp3;
    int temp4;
      status=read_4c_file("mstrnx.data/solar.dat", &temp1, &temp2, &temp3, &S_0, &temp4);
      if (status !=0) {
      printf("Error reading Solar.dat\n");
      return EXIT_FAILURE;
      }
      
      free(temp1);
      free(temp2);
      free(temp3);
  } 
  
  //--------------------- Plot color -------------------------//

  plscolbg (255, 255, 255);   /* background color white */
  plscol0  (15, 0, 0, 0);     /* set color 15 to black  */
  plsdev ("xwin"); /* if not called, the user is asked! */
  plinit ();
  
  plssub(2,2);

  
  //--------------------- Calculate pressure for layers and levels -------------//

  for(ilev=0; ilev<nlev; ilev++) {                              /* Calculation of the Pressure at the Levels p[ilev] */
    p[ilev]=PSURF*ilev/(nlev-1);
  }
    
  for(ilyr=0; ilyr<nlyr; ilyr++) {                            /* Calculation of the Pressure in the Layers*/
    plyr[ilyr]= 0.5*(p[ilyr]+p[ilyr+1]);
  }
   
  
  //----------------------- Reading afglus ------------------------//
  
  {
    
    double *ztemp;
    double *ptemp;
    double *Ttemp;
    double *h2otemp;
    double *o3temp;

    int status;
    int ntemp;
    int itemp;
    unsigned int mintemp = 0;
    unsigned int mintemp2 = 0;
    
    status = read_5c_file("mstrnx.data/afglus.atm", &ztemp, &ptemp, &Ttemp, &h2otemp, &o3temp, &ntemp);
    if (status !=0) {
      printf("Error reading Temperature profile\n");
      return EXIT_FAILURE;
    }

    /*
      for (itemp=0; itemp<ntemp; itemp++) {
      printf("itemp = %d, ztemp = %f, ptemp = %f, Ttemp = %f, h2otemp = %f, o3temp = %f\n", itemp, ztemp[itemp], ptemp[itemp], Ttemp[itemp], h2otemp[itemp], o3temp[itemp]);
      }
    */  
    
    //-------------------- Interpolate T, h2o and o3 to layers --------------//
    
    for (ilyr=0; ilyr<nlyr; ilyr++) {
      // printf("ilev = %d, searching for %f\n", ilyr, plyr[ilyr]);
      for (itemp=0; itemp<ntemp; itemp++) {
	if ( abs(ptemp[mintemp]- plyr[ilyr])>abs(ptemp[itemp]-plyr[ilyr])) {
	  mintemp=itemp;
	}
      }
      if (ptemp[mintemp]<plyr[ilyr]) {
	mintemp2 = mintemp+1;
      }
      else {
	mintemp2 = mintemp-1;
      }
      // printf("\tmintemp  = %d,\tp = %f\n", mintemp, ptemp[mintemp]);
      // printf("\tmintemp2 = %d,\tp = %f\n", mintemp2, ptemp[mintemp2]);

      const double weight2 = abs(ptemp[mintemp]-plyr[ilyr]);
      const double weight1 = abs(ptemp[mintemp2]-plyr[ilyr]);
      const double norm = weight1 + weight2;
      
      /* T in levels?!?! */
      T[ilyr] = Ttemp[mintemp]*weight1 + Ttemp[mintemp2]*weight2;
      T[ilyr] /= norm;
      
      h2o[ilyr] = h2otemp[mintemp]*weight1 + h2otemp[mintemp2]*weight2;
      h2o[ilyr] /= norm;

      o3[ilyr] = o3temp[mintemp]*weight1 + o3temp[mintemp2]*weight2;
      o3[ilyr] /= norm;
    }
  }
  
  
  for (ilyr=0; ilyr<nlyr; ilyr++)  {
    printf("ilyr = %d, plyr = %f,\tT = %f,\th2o = %f,\to3 = %f\n", ilyr, plyr[ilyr], T[ilyr], h2o[ilyr], o3[ilyr]);
  }
  
  /* 
     for (ilyr=0; ilyr<nlyr; ilyr++)  {
     printf("T[%d] = %f\n", ilyr, T[ilyr]);
     }
     for (ilyr=0; ilyr<nlyr; ilyr++)  {
     printf("o3[%d] = %f\n", ilyr, o3[ilyr]);
     }
     for (ilyr=0; ilyr<nlyr; ilyr++)  {
     printf("h2o[%d] = %f\n", ilyr, h2o[ilyr]);
     }
  */
  
   
  
  
  //--------------------------------- Starting while-loop ----------------------------------//
  //----------------------------------------------------------------------------------------//
  
  while (timesteps*deltat<TIME_MAX) {
    // printf("\nNew time %d: T = %f\n", (int)(timesteps*deltat), T[nlyr-1]);
    timesteps++;
  
    for(ilev=0; ilev<nlev; ilev++) {                    /* Reseting edn, eup and enet */
      edn[ilev] = 0;
      eup[ilev] = 0;
      if(ilev < nlyr) //enet is defined for ilyr only
        enet[ilev]=0;
    }
    
    
    //--------------------- Calculate z ----------------------//
    
    for (ilyr=0;ilyr<nlyr; ilyr++) {                         /* z for layers */
      deltazlyr[ilyr]=(Ra*T[ilyr]*deltap)/(plyr[ilyr]*g);  
    }
      
    zlyr[nlyr-1]=0; //check this one time
    for (ilyr=nlyr-2; ilyr >= 0; ilyr--) {
      zlyr[ilyr]=zlyr[ilyr+1]+deltazlyr[ilyr];
    }
      
     
    for (ilev=0;ilev<nlev-1; ilev++) {                      /* z for levels */
      deltazlev[ilev]=(Ra*T[ilev]*deltap)/(plyr[ilev]*g);  
    }
      
    z[nlev-1]=0;
    for (ilev=nlev-2; ilev >= 0; ilev--) {
      z[ilev]=z[ilev+1]+deltazlev[ilev];
     
    }
      
      
      
    //--------------------- Use k-distribution ---------------------//
      
    /* mu_0 is defined for a 8-part cycle of the earth, i.e. a timestep of 3 hours! */
    if(++mu_counter == mu_counterlimit)
      mu_counter = 0;

    //calculate it every two days!
    if((timesteps-1) % 10 == 0)
      //warning, wavelength is in nm!
      status = ck_mstrnx (z, plyr, T, h2o, o3, nlev, /*lyrflag*/ 1, c_co2, c_n2o, c_co, c_ch4, c_o2, &dtaumol, &wgt, &wavelength, &nbnd, &nch);
    
    for(iv=0; iv<nbnd; iv++) {
      //printf("iv = %d, wavelength = %e, S0 = %e\n", iv, wavelength[iv], S_0[iv]);
	  
      
      //------------- Use schwarzschild and sum edn eup --------//
	
      
      for(iq=0; iq<nch[iv]; iq++) {
	schwarzschild2(dtaumol[iv][iq], T, nlev, T[nlyr-1], edntmp, euptmp, wavelength[iv]*1e-6,wavelength[iv+1]*1e-6, tmplev, tmplyr);	

	for(ilev=0; ilev<nlev; ilev++) {
	  edn[ilev] += edntmp[ilev]*wgt[iv][iq];
	  eup[ilev] += euptmp[ilev]*wgt[iv][iq];
	  //printf("ilev = %d, iv = %d of %d, iq = %d of %d, wgt = %e, edntmp = %e, euptmp = %e\n", ilev, iv, nbnd, iq, nch[iv], wgt[iv][iq], edntmp[ilev], euptmp[ilev]);
	}
	
	//------------------- Include Rodents---------------------//  
	if(mu_0[mu_counter] > 0) {
	  rodents_solar(nlyr, dtaumol[iv][iq], omega_0, gassy, f, S_0[iv], mu_0[mu_counter], gr_albedo, edntmp, euptmp, edirtmp);
	  
	  for(ilev=0; ilev<nlev; ilev++) {
	    edn[ilev] += (edntmp[ilev]+edirtmp[ilev])*wgt[iv][iq];
	    eup[ilev] += euptmp[ilev]*wgt[iv][iq];
	    //printf("ilev = %d, iv = %d of %d, iq = %d of %d, wgt = %e, edntmp = %e, euptmp = %e, edirtmp = %e\n", ilev, iv, nbnd, iq, nch[iv], wgt[iv][iq], edntmp[ilev], euptmp[ilev], edirtmp[ilev]);
	  }
	}
      }
    }
     
    //---------------- Calculate enet -----------------------// 
    
    for(ilyr=0; ilyr<nlyr-1; ilyr++) {
      enet[ilyr] = eup[ilyr+1] + edn[ilyr] - eup[ilyr] - edn[ilyr+1];
     
    }

    for (ilev=0; ilev<nlev; ilev++) {
      //printf("ilev = %d, eup = %f, edn = %f\n", ilev, eup[ilev], edn[ilev]);
    }

    //---------------- Calculate temperature gain per layer -----------//
    
    for (ilyr=0; ilyr<nlyr-1; ilyr++)  {
      deltaT[ilyr]=(enet[ilyr]*g*deltat)/((deltap*100.0)*cp);
      T[ilyr] +=deltaT[ilyr];
      //printf("dT[ilyr] = %f\n", ilyr, deltaT[ilyr]);
    }
  
    //----------------Calculate Temperature gain for ground layer ---------//
    
    deltaT[nlyr-1]=((edn[nlyr-1]-eup[nlyr-1])*g*deltat)/((deltap*100.0)*cp);
    T[nlyr-1] += deltaT[nlyr-1];


    //---------- Convert T to theta ------------//
    
    for (ilyr=0; ilyr<nlyr; ilyr++) {                           
      theta[ilyr]=pow(PSURF/plyr[ilyr], kappa)*T[ilyr];
      //printf("%d %f\n", ilyr, theta[ilyr]);
    }
    
    
    //------------ Convection -------------------//
    
    instabil = FALSE;
    
    for (ilyr=0; ilyr<nlyr;ilyr++) {                            /* Testing for Instability */
      if (theta[ilyr+1]>theta[ilyr]) {
	instabil = TRUE;
	break;
      }
    }

    if (instabil) {                                             /* Convection - Sorting of Layers according to theta */
      qsort (theta, nlyr, sizeof(double), sortfunc);
      //  printf ("Konvektion :-)\n");
    }
    
    for (ilyr=0; ilyr<nlyr; ilyr++) {                           /* Conversion from theta to T */
      T[ilyr] = theta[ilyr] /  (pow(PSURF/plyr[ilyr], kappa));

      //printf("%d %f\n", ilyr, theta[ilyr]);
    }

    
    //-------------------- Print every x timestep ------------//
    
    if(timesteps % (int)(10) == 0) {
      printf ("timestep %d\n", timesteps);

      /* Printing time in readable format
       * In order to remove days, hours or minutes
       * just add // in front of the appropriate line
       */
      int time = timesteps*deltat;
      printf("time ");
      printf("%dd = ", (int)(time / 3600 / 24));
      printf("%dh = ", (int)(time / 3600 ));
      printf("%dmin = ", (int)(time / 60 ));
      printf("%ds\n", time);

      printf(" Tsurf=%f\n", T[nlyr-1]);

      for (ilyr=0; ilyr<nlyr; ilyr++) {
	deltaTday[ilyr]=deltaT[ilyr]*86.4/2;
      }

      plotall(nlyr, T, plyr, zlyr, deltaTday);

      /* for (ilyr=0; ilyr<nlyr; ilyr++) { */
      /*   printf("ilyr %d, z=%f,  plyr=%f,theta=%f, T=%f\n", ilyr, z[ilyr], plyr[ilyr],theta[ilyr], T[ilyr]); */
      /* } */

      for (ilyr=0; ilyr<nlyr; ilyr++) {
      	printf("p%f, edn%f, eup%f, enet%f\n", p[ilyr], edn[ilyr], eup[ilyr], enet[ilyr]);
      }
      printf("p%f, edn%f, eup%f\n", p[ilyr], edn[ilyr], eup[ilyr]);
    }
   
  }                                     

  free(tmplyr);
  free(tmplev);
  //----------------- End of while-loop---------------------//
  
  printf("\nTime %d: T = %f\n", (int)(timesteps*deltat), T[nlyr-1]);
  for (ilyr=0; ilyr<nlyr; ilyr++) {
    printf("%d %f\n", ilyr, T[ilyr]);
  }

  printf("\nTsurf=%f\n", T[nlyr-1]);

  return 0;
}
