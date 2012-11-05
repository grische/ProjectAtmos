#ifndef SCHWARZ_H_INCLUDED
#define SCHWARZ_H_INCLUDED

double boltzmann(const double T);
double boltzmann_overpi(const double T);
int schwarzschild(const double tau, const double *T, const int nlev, const double Ts, double *edn, double *eup);

#endif
