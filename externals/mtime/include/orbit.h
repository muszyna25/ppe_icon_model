#ifndef ORBIT_H
#define ORBIT_H

void orbit_vsop87(double jde, double *ra, double *dec, double *dis, double *gha);

void orbit_kepler(double pvetim, double *pdisse, double *pdec, double *pra);

void earth_position(double t, double *l, double *b, double *r);

#endif
