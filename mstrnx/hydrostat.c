#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hydrostat.h"

#define EPSILON 1.0E-6

/****************************************************************/
/* provide a hydrostatic profile with a given lapse rate Gamma, */
/* surface temperature T0, surface pressure p0; the z profile   */
/* is equidistant between 0 and zmax [m] with nlev levels;      */
/* denscol is the density integrated over each layer [cm**-2]   */
/****************************************************************/

void hydrostatic_profile (int nlev, double zmax, double p0, double T0, double Gamma, int verbose,
			  double **z, double **p, double **T, double **dens, double **denscol)
{ 
  int ilyr=0, ilev=0;
  double dz=0;


  *z       = calloc (nlev-1, sizeof(double));
  *p       = calloc (nlev-1, sizeof(double));
  *T       = calloc (nlev-1, sizeof(double));
  *dens    = calloc (nlev-1, sizeof(double));
  *denscol = calloc (nlev-1, sizeof(double));

  dz = zmax / (double) (nlev-1);

  /* pressure, density */
  if (verbose) {
    fprintf (stderr, "--------------------- hydrostatic profile --------------------\n");
    fprintf (stderr, "     z[km]  p[hPa]     T[K]  dens [1/cm**3]  denscol [1/cm**2]\n");
  }

  for (ilev=0;ilev<nlev;ilev++) {
    (*z)[ilev] = zmax - dz * (double) ilev;
    (*T)[ilev] = T0 - Gamma*(*z)[ilev];
    (*p)[ilev] = p0*pow(1.0-Gamma*(*z)[ilev]/T0,G_EARTH/R_AIR/Gamma);  // hydrostatic pressure profile
    
    (*dens)[ilev] = (*p)[ilev]*100.0/K_BOLTZMANN/(*T)[ilev]/1.0E6;
  }
	    
  /* integrate density over layer */
  for (ilyr=0;ilyr<nlev-1;ilyr++)
    calc_denscol_hydrostatic_layer ((*p)[ilyr+1], (*p)[ilyr], &((*denscol)[ilyr]));
  

  if (verbose) {
    fprintf (stderr, "%10.2f %7.2f %8.3f    %9.6e\n", (*z)[0]/1000.0, (*p)[0], (*T)[0], (*dens)[0]);
    for (ilyr=0;ilyr<nlev-1;ilyr++)
      fprintf (stderr, "%10.2f %7.2f %8.3f    %9.6e       %9.6e\n", (*z)[ilyr+1]/1000.0, (*p)[ilyr+1], (*T)[ilyr+1], (*dens)[ilyr+1], (*denscol)[ilyr]);
    
    fprintf (stderr, "--------------------------------------------------------------\n\n");
  }    
}




/****************************************************************/
/* calculate denscol (if c==NULL) or cdenscol for a profile     */
/* with given pressure, assuming that the concentration         */
/* c [mol/mol] varies linearely with p within each layer        */       
/****************************************************************/

void calc_denscol_hydrostatic_profile (int nlev, double *p, double *T, double *c, int verbose,
				       double **cdenscol)
{ 
  int ilyr=0;
  double p1=0, p2=0, c1=0, c2=0;

  double *denscol  = calloc (nlev-1, sizeof(double));

  if (verbose)
    if (c==NULL)
      fprintf (stderr, "c==NULL - assuming unity concentration\n");

  *cdenscol = calloc (nlev-1, sizeof(double));

  /* pressure, density */
  if (verbose) {
    if (c!=NULL) {
      fprintf (stderr, "--------------------------------------------- hydrostatic profile ---------------------------------------------\n");
      fprintf (stderr, " p[hPa]     T[K]  c[mol/mol]  dens [1/cm**3]  cdens [1/cm**3]  denscol [1/cm**2]  cdenscol [1/cm**2]\n");
    }
    else {
      fprintf (stderr, "--------------------------------------------- hydrostatic profile ---------------------------------------------\n");
      fprintf (stderr, " p[hPa]     T[K]  dens [1/cm**3]  denscol [1/cm**2]\n");
    }
  }

  /* integrate density over layer */
  for (ilyr=0;ilyr<nlev-1;ilyr++) {
    p1=p[ilyr+1];   /* lower level pressure [hPa] */
    p2=p[ilyr];     /* upper level pressure [hPa] */

    if (c!=NULL) {
      c1=c[ilyr+1];   /* lower level concentration [mol/mol] */
      c2=c[ilyr];     /* upper level concentration [mol/mol] */
    }

    calc_denscol_hydrostatic_layer (p1, p2, &(denscol[ilyr]));
    if (c!=NULL) 
      calc_denscol_hydrostatic_layer_linear_conc (p1, p2, c1, c2, &((*cdenscol)[ilyr]));
    else 
      (*cdenscol)[ilyr]=denscol[ilyr];
  }

  if (verbose) {
    if (c!=NULL) {
      fprintf (stderr, "%7.2f %8.3f  %6.4e    %9.6e     %9.6e\n", p[0], T[0], c[0], p[0]*100.0/K_BOLTZMANN/T[0]/1e6, p[0]*100.0/K_BOLTZMANN/T[0]/1e6 * c[0]);
      for (ilyr=0;ilyr<nlev-1;ilyr++)
	fprintf (stderr, "%7.2f %8.3f  %6.4e    %9.6e     %9.6e       %9.6e        %9.6e\n", 
		 p[ilyr+1], T[ilyr+1], c[ilyr+1], 
		 p[ilyr+1]*100.0/K_BOLTZMANN/T[ilyr+1]/1e6, p[ilyr+1]*100.0/K_BOLTZMANN/T[ilyr+1]/1e6*c[ilyr+1], 
		 denscol[ilyr], (*cdenscol)[ilyr]);
    }
    else {
      fprintf (stderr,   "%7.2f %8.3f    %9.6e\n", p[0], T[0], p[0]*100.0/K_BOLTZMANN/T[0]/1e6);
      for (ilyr=0;ilyr<nlev-1;ilyr++)
	fprintf (stderr, "%7.2f %8.3f    %9.6e       %9.6e\n", 
		 p[ilyr+1], T[ilyr+1],
		 p[ilyr+1]*100.0/K_BOLTZMANN/T[ilyr+1]/1e6, 
		 denscol[ilyr]);
    }
  }

  if (verbose)
    fprintf (stderr, "---------------------------------------------------------------------------------------------------------------\n");

  free (denscol);
}	    



/* Integrate air density over an atmospheric layer, assuming that the pressure */
/* is hydrostatic and the temperature varies linearely over the layer;         */
/* the solution is exact for this case; specify p[hPa] and T[K] to obtain the  */
/* integrated column in cm^-2                                                  */

void calc_denscol_hydrostatic_layer (double p1, double p2, double *denscol)
{
  /* factor 100.0 for hPa -> Pa and 1/1.0E4 for 1/m2 -> 1/cm2 */
  *denscol = R_AIR/K_BOLTZMANN/G_EARTH*100.0*(p1-p2)/1.0E4;  
}


/* Integrate trace gas density over a layer between z1 and z2,                     */
/* assuming that the pressure is hydrostatic and the temperature and concentration */
/* vary linearely over the layer; the solution is exact for this case; specify     */
/* z[m], p[hPa], T[K] and c[mol/mol] to obtain the integrated columns in cm^-2     */

void calc_denscol_hydrostatic_layer_linear_conc (double p1, double p2, double c1, double c2, 
						 double *cdenscol)
{
  /* factor 100.0 for hPa -> Pa and 1/1E4 for 1/m2 -> 1/cm2 */
  *cdenscol = R_AIR/K_BOLTZMANN/G_EARTH*100.0*(p1-p2)*(c1+c2)/2.0/1E4;
}
