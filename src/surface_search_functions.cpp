/*IMa2p 2009-2015 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, and Arun Sethuraman */

/* this file is numerical recipes stuff with some small modifications */
#undef GLOBVARS
#include "imamp.hpp"              // may not really be needed in this particular file

/*********** LOCAL STUFF **********/

#define NRANSI
#define NR_END 1
#define FREE_ARG char*

/* local function prototypes */

static double sign (double a, double b);

/* local functions */


static double
sign (double a, double b)
{
  if (b > 0.0)
    return fabs (a);

  else
    return (-fabs (a));
}


/************ GLOBAL FUNCTIONS ******************/

/* NR stuff for finding minimum using golden section */

#define gold            1.618034
#define glimit          100.0
#define tiny            1.0e-20
#define r               0.61803399

/* modified mnbrak() and golden() from NR code.  used only for the one dimensional optimization */
void
mnbrakmod (int ndim, int firsttree, int lasttree, double *ax, double *bx,
           double *cx, double *fa, double *fb, double *fc,
           double (*func) (int, int, int, double, int), int ifneeded)
{

  /* Programs using routine MNBRAK must supply an EXTERNAL
     FUNCTION func(x:REAL):REAL FOR which a minimum is TO be found */
  /* I converted this to a function and included code to check for
     a failure to find a braket */
  double ulim, u, rr, q, fu, dum;
  *fa = func (ndim, firsttree, lasttree, *ax, ifneeded);
  *fb = func (ndim, firsttree, lasttree, *bx, ifneeded);
  if (*fb > *fa)
  {
    dum = *ax;
    *ax = *bx;
    *bx = dum;
    dum = *fb;
    *fb = *fa;
    *fa = dum;
  }
  *cx = fabs (*bx + gold * (*bx - *ax));
  *fc = func (ndim, firsttree, lasttree, *cx, ifneeded);
_L1:
  /* need to add some excapes for floating point problems */
  if (*fb >= *fc && *fb > -DBL_MAX && !(*fb == 0 && *fc == 0))
  {
    rr = (*bx - *ax) * (*fb - *fc);
    q = (*bx - *cx) * (*fb - *fa);
    u = (double) *bx - ((*bx - *cx) * q - (*bx - *ax) * rr) /
      (2.0 * sign (FMAX (fabs (q - rr), (float) tiny), q - rr));
    ulim = *bx + glimit * (*cx - *bx);
    if ((*bx - u) * (u - *cx) > 0.0)
    {
      fu = func (ndim, firsttree, lasttree, u, ifneeded);
      if (fu < *fc)
      {
        *ax = *bx;
        *fa = *fb;
        *bx = u;
        *fb = fu;
        goto _L1;
      }
      if (fu > *fb)
      {
        *cx = u;
        *fc = fu;
        goto _L1;
      }
      u = *cx + gold * (*cx - *bx);
      fu = func (ndim, firsttree, lasttree, u, ifneeded);
    }
    else if ((*cx - u) * (u - ulim) > 0.0)
    {
      fu = func (ndim, firsttree, lasttree, u, ifneeded);
      if (fu < *fc)
      {
        *bx = *cx;
        *cx = u;
        u = *cx + gold * (*cx - *bx);
        *fb = *fc;
        *fc = fu;
        fu = func (ndim, firsttree, lasttree, u, ifneeded);
      }
    }
    else if ((u - ulim) * (ulim - *cx) >= 0.0)
    {
      u = ulim;
      fu = func (ndim, firsttree, lasttree, u, ifneeded);
    }
    else
    {
      u = *cx + gold * (*cx - *bx);
      fu = func (ndim, firsttree, lasttree, u, ifneeded);
    }
    *ax = *bx;
    *bx = *cx;
    *cx = u;
    *fa = *fb;
    *fb = *fc;
    *fc = fu;
    goto _L1;
  }
  if (*fb < -DBL_MAX)
    printf (" minus infinity in mnbrakmod \n");
}


#undef gold
#undef glimit
#undef tiny
double
goldenmod (int ndim, int firsttree, int lasttree, double ax, double bx,
           double cx, double tol_, double *xmin, double (*func) (int, int,
                                                                 int, double, int), int ifneeded)
{

  /* Programs using routine GOLDEN must supply an EXTERNAL
     FUNCTION func(x:REAL):REAL whose minimum is TO be found. */
  double Result, f1, f2, c, x0, x1, x2, x3;
  c = 1.0 - r;
  x0 = ax;
  x3 = cx;
  if (fabs (cx - bx) > fabs (bx - ax))
  {
    x1 = bx;
    x2 = bx + c * (cx - bx);
  }
  else
  {
    x2 = bx;
    x1 = bx - c * (bx - ax);
  }
  f1 = func (ndim, firsttree, lasttree, x1, ifneeded);
  f2 = func (ndim, firsttree, lasttree, x2, ifneeded);
  while (fabs (x3 - x0) > tol_ * (fabs (x1) + fabs (x2)) && (x3 > -DBL_MAX))
  {
    if (f2 < f1)
    {
      x0 = x1;
      x1 = x2;
      x2 = r * x1 + c * x3;
      f1 = f2;
      f2 = func (ndim, firsttree, lasttree, x2, ifneeded);
      continue;
    }
    x3 = x2;
    x2 = x1;
    x1 = r * x2 + c * x0;
    f2 = f1;
    f1 = func (ndim, firsttree, lasttree, x1, ifneeded);
  }
  if (f1 < f2)
  {
    Result = f1;
    *xmin = x1;
  }
  else
  {
    Result = f2;
    *xmin = x2;
  }
  if (x3 < -DBL_MAX)
    printf (" minus infinity in goldenmod \n");
  return Result;
}
