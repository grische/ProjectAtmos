#ifndef SCHWARZ_H_INCLUDED
#define SCHWARZ_H_INCLUDED

double boltzmann(const double T);
double boltzmann_overpi(const double T);
double planck(const double T, const double lamda);
int schwarzschild(const double *deltatau, const double *T, const int nlev, const double Ts, double *edn, double *eup, const double lambda);

#endif
