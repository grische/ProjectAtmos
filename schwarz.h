#ifndef SCHWARZ_H_INCLUDED
#define SCHWARZ_H_INCLUDED

double boltzmann(const double T);
double boltzmann_overpi(const double T);
double planck(const double T, const double lambda, const double lambda2);
float plkint_(const double *wvllo, const double *wvlhi, const double *T); //comes from plkint.f
int schwarzschild(const double *deltatau, const double *T, const int nlev, const double Ts, double *edn, double *eup, const double lambda);
int schwarzschild2(const double *deltatau, const double *T, const int nlev, const double Ts, double *edn, double *eup, const double lambdalow, const double lambdahigh);

#endif
