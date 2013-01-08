#ifndef __hydrostat_h
#define __hydrostat_h 1

#if defined (__cplusplus)
extern "C" {
#endif

#define K_BOLTZMANN 1.3806503E-23
#define R_AIR 287.0
#define G_EARTH 9.81


void hydrostatic_profile (int nlev, double zmax, double p0, double T0, double Gamma, int verbose,
			  double **z, double **p, double **T, double **dens, double **denscol);

void calc_denscol_hydrostatic_profile (int nlev, double *p, double *T, double *c, int verbose,
				       double **cdenscol);

void calc_denscol_hydrostatic_layer (double p1, double p2, double *denscol);
void calc_denscol_hydrostatic_layer_linear_conc (double p1, double p2, double c1, double c2, 
						 double *cdenscol);


#if defined (__cplusplus)
}
#endif

#endif /* __hydrostat_h */

