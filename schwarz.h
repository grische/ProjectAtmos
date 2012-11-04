#ifndef SCHWARZ_H_INCLUDED
#define SCHWARZ_H_INCLUDED

double boltzman(const double T);
double boltzman_overpi(const double T);
int schwarzschild(const double tau, const double *T, const int nlev, const double Ts, double *edn, double *eup);

#endif
