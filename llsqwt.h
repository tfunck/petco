/*******************************************************************************

  llsqwt.h

  Version:
  2002-09-27 VO
  2002-10-28 VO

*******************************************************************************/
#ifndef _LLSQWT_H
#define _LLSQWT_H
/******************************************************************************/
int LLSQWT_TEST;
/******************************************************************************/
extern int llsqwt(
  double *x, double *y, int n, double *wx, double *wy, double tol, double *w,
  double *ic, double *slope, double *nwss, double *sic, double *sslope,
  double *cx, double *cy
);
extern int best_llsqwt(
  double *x, double *y, double *wx, double *wy, int nr, int min_nr, int mode,
  double *slope, double *ic, double *nwss, double *sslope, double *sic,
  double *cx, double *cy, int *bnr
);
/******************************************************************************/
extern int llsqperp(double *x, double *y, int nr,
  double *slope, double *ic, double *ssd
);
extern int llsqperp3(double *x, double *y, int nr,
  double *slope, double *ic, double *ssd
);
/******************************************************************************/
extern int quadratic(double a, double b, double c, double *m1, double *m2);
/******************************************************************************/
extern int medianline(double *x, double *y, int nr, double *slope, double *ic);
/******************************************************************************/
#endif

