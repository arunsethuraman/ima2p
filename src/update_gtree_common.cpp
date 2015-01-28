/*IMa2p 2009-2015 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, and Arun Sethuraman */

#undef GLOBVARS
#include "imamp.hpp"
#include "update_gtree_common.hpp"
#include "xtrapbits.hpp"

/* This file:
  Includes local functions with three locals for declaration:
  declared locally - see static declarations below
  declared in update_gree_common.h
	used by  this file,  update_gtree.c and update_t_NW.c
  declared globally in imamp.h */

#define M_LNPI 1.14472988584940017414342735135
#ifndef M_LN2
#define M_LN2 0.69314718055994530941723212146
#endif /* M_LN2 */
//extern double logfact[100 * ABSMIGMAX + 1];  JH  3_1_2010 was duplicated
extern double loghalffact[50 * ABSMIGMAX];
extern int BITNUMBERTRUE[256];

#define  MIGCLOSEFRAC  0.9  // when both edges updating, chance of zero or 1 migration event in last migration period

#define gsl_sf_lnfact(n) logfact[(n)]
#define IMA_math_lnfact_half(v,n,quotient,remainder) \
  (quotient) = (n) / 2; \
  (remainder) = (n) % 2; \
  if ((remainder) == 0) \
    { \
      (v) = logfact[(quotient)]; \
    } \
  else \
    { \
      (v) = M_LNPI / 2 \
          + logfact[2 * (quotient) + 2] \
          - logfact[(quotient) + 1] \
          - (2 * (quotient) + 2) * M_LN2; \
    }

// variables declared in update_gtree
extern struct edgemiginfo oldedgemig;
extern struct edgemiginfo oldsismig;
extern struct edgemiginfo newedgemig;
extern struct edgemiginfo newsismig;


/* declarded in update_gtree.h */
/* extern int mrootdrop; not used here */
/* extern double lmrootdrop; not used here */
/* extern struct genealogy_weights holdgweight_updategenealogy; not used here */
/* extern struct genealogy_weights holdallgweight_updategenealogy; not used here */
/* extern struct probcalc holdallpcalc_updategenealogy; not used here */
extern int rootmove;            /* declared in update_gtree.c */
static int holddownA[MAXLINKED];
static int medgedrop;           /* NO USAGE: SHOULD BE DELETED */
static double lmedgedrop;       /* NO USAGE: SHOULD BE DELETED */
static double holdsisdlikeA[MAXLINKED];
struct genealogy holdgtree;
static struct migstruct mlist[ABSMIGMAX];     // used in getmprob() mlist very wasteful of space - could use dynamic memory and checkmig() 

//variables local to this file 
#define ESTARTLENGTH  20        // fairly arbitrary starting size of sgEvent and sgEvent_indx
#define EADDLENGTH  10          // fairly aribtrary amount to add when needed

static int elength;
static struct gtreeevent *sgEvent; 
static unsigned long *sgEvent_indx;
static struct edge *copyedge;
static double *nnminus1;

/* prototype of local functions */
static void IMA_initmemory_edgemiginfo (struct edgemiginfo *em);
static double integrate_coalescent_term (int cc, double fc, double hcc,
                                         double max, double min);
static double integrate_migration_term (int cm, double fm, double max,
                                        double min);
static double integrate_migration_term_expo_prior (int cm, double fm,
                                                   double exmean);
static double integrate_splitrate_term (double f, double max, double min);
static int simmpath (int ci, struct edgemiginfo *edgem, int period, int numm,
                     int lastm, double timein, double upt, int pop,
                     int constrainpop);

/********* LOCAL FUNCTIONS **************/

void
IMA_initmemory_edgemiginfo (struct edgemiginfo *em)
{
  if (modeloptions[NOMIGRATION] == 0)
  {
    em->mtimeavail = static_cast<double *> (malloc (npops * sizeof (double)));
    em->mp = static_cast<int *> (malloc (npops * sizeof (int)));
    em->sisid = -1;
    em->li = -1;
  }
  else
  {
    em->mtimeavail = static_cast<double *> (malloc (npops * sizeof (double)));
    em->mp = static_cast<int *> (malloc (npops * sizeof (int)));

    em->sisid = -1;
    em->li = -1;
  }
  return;
}

double
integrate_coalescent_term (int cc, double fc, double hcc, double max,
                           double min)
{                               /* note fc includes inheritance scalar, i.e. it is fc/h  - see treeweight  
                                   the hcc term is cc*log(h) */
  double p;
  double a, b, c, d;
  double ug,lg,ugalt, fullg;
  if (cc > 0)
  {
    assert (fc > 0);
    if (min == 0)
    {
      ug = uppergamma ((int) cc - 1, 2 * fc / max);
      if (cc > 1)
      {
        fullg = logfact[cc-2];
/* if uppergamma is too small or too close to the full integral  (fullg) then use lowergamma */
        if (fullg - ug < 1e-15 || fullg - ug  > 7.0978271289338397e+02)
        {
          lg  = lowergamma ((int) cc - 1, 2 * fc / max);
          if (fullg > lg)
          {
            LogDiff(ugalt,fullg,lg);
            if (fabs(ugalt-ug) > 1e-10)
              ug=ugalt;
          }
        }
      }
      p = ug + LOG2 - hcc + (1 - cc) * log (fc);
    }
    else
    {
#ifdef MORESTABLE
      a = uppergamma ((int) cc - 1, 2 * fc / max);
      b = uppergamma ((int) cc - 1, 2 * fc / min);
      if (a <= b)
        IM_err (IMERR_LOGDIFF, " pos1 in integrate_coalescent_term() cc %d, fc %lf, hcc %lf, max %lf, min %lf, a %lf  b %lf",cc,fc,hcc,max,min,a,b); 
      LogDiff (p, a, b);
      p += (LOG2 - hcc + (1 - cc) * log (fc));
#else
      p = log (exp (uppergamma ((int) cc - 1, 2 * fc / max)) 
               - exp (uppergamma ((int) cc - 1, 2 * fc / min))) 
          + LOG2 - hcc + (1 - cc) * log (fc);
#endif /* MORESTABLE */
    }
  }
  else
  {
    if (2 * fc / max > 0)       /* cc == 0 */
    {
      if (min == 0)
      {
//  fc values coming out very large.  this deals with it,  but why so large ????
#ifdef MORESTABLE
        a = log (max) - 2.0 * fc / max;
        b = LOG2 + log (fc) + uppergamma (0, 2.0 * fc / max);
        if (a <= b)
          IM_err (IMERR_LOGDIFF, " pos2 in integrate_coalescent_term() cc %d, fc %lf, hcc %lf, max %lf, min %lf, a %lf  b %lf",cc,fc,hcc,max,min,a,b); 
        LogDiff (p, a, b);
#else
        p = log (max * exp (-2 * fc / max) - 2 * fc * exp (uppergamma (0, 2 * fc / max)));
#endif /* MORESTABLE */
      }
      else
      {
#ifdef MORESTABLE
        a = uppergamma (0, 2 * fc / max);
        b = uppergamma (0, 2 * fc / min);
        if (a <= b)
          IM_err (IMERR_LOGDIFF, " pos3 in integrate_coalescent_term() cc %d, fc %lf, hcc %lf, max %lf, min %lf, a %lf  b %lf",cc,fc,hcc,max,min,a,b); 
        LogDiff (c, a, b);
        c += LOG2 + log (fc);
        a = log (max) - 2.0 * fc / max;
        b = log (min) - 2.0 * fc / min;
        if (a <= b)
          IM_err (IMERR_LOGDIFF, " pos4 in integrate_coalescent_term() cc %d, fc %lf, hcc %lf, max %lf, min %lf, a %lf  b %lf",cc,fc,hcc,max,min,a,b); 
        LogDiff (d, a, b);
        if (d <= c)
          IM_err (IMERR_LOGDIFF, " pos5 in integrate_coalescent_term() cc %d, fc %lf, hcc %lf, max %lf, min %lf, d %lf  c %lf",cc,fc,hcc,max,min,d,c); 
        LogDiff (p, d, c);
#else
        p = log (max * exp (-2 * fc / max) - min * exp (-2 * fc / min)
                 - 2 * fc * (exp (uppergamma (0, 2 * fc / max)) -
                             exp (uppergamma (0, 2 * fc / min))));
#endif /* MORESTABLE */
      }
    }
    else
    {
      p = log (max - min);
    }
  }

  return p;
}                               /* integrate_coalescent_term */

double
integrate_migration_term (int cm, double fm, double max, double min)
/* checked pretty extensively against numerical integration  9/12/2010 */
{
  double p;
  double a, b, c;

  double ug,lg, lgalt, fullg;
  if (cm > 0)
  {
    assert (fm > 0);
    if (min == 0)
    {
      lg =lowergamma ((int) cm + 1, fm * max);
      fullg = logfact[cm];
      /* if lowergamma is too small or too close to the full integral  (fullg) then use uppergamms */
      if (fullg - lg < 1e-15 || fullg - lg  > 7.0978271289338397e+02)
      {
        ug = uppergamma ((int) cm + 1, fm * max);
        if (fullg > ug)
        {
          LogDiff(lgalt,fullg,ug);
          if (fabs(lgalt-lg) > 1e-12)
            lg=lgalt;
        }
      }
      p = (-1 - cm) * log (fm) + lg;
    }
    else
    {

#ifdef MORESTABLE
      a = uppergamma (cm + 1, fm * min);
      b = uppergamma (cm + 1, fm * max);
      if (a<=b)
          IM_err (IMERR_LOGDIFF, " pos 1 in integrate_migration_term() cm %d, fm %lf, max %lf, min %lf, a %lf  b %lf",cm,fm,max,min,a,b); 
      LogDiff (c, a, b);
      p = (-1 - cm) * log (fm) + c;
#else
      p = (-1 - cm) * log (fm) 
          + log (exp (uppergamma (cm + 1, fm * min)) 
                 - exp (uppergamma (cm + 1, fm * max)));
#endif /* MORESTABLE */
    }
  }
  else
  {
    if (fm > MPRIORMIN)                 /* cm == 0   use a cutoff, because sometimes the value that comes in is very low, when it should be zero */
    {
      if (min == 0)
      {
        if (max == MPRIORMIN)
        {
          p = 0;
        }
        else
        {
#ifdef MORESTABLE
          a = 0.0;
          b = -fm * max;
          if (a<=b)
            IM_err (IMERR_LOGDIFF, " pos 2 in integrate_migration_term() cm %d, fm %lf, max %lf, min %lf, a %lf  b %lf",cm,fm,max,min,a,b); 
          LogDiff (c, a, b);
          p = c - log (fm);
#else
          p = log ((1 - exp (-fm * max)) / fm);
#endif /* MORESTABLE */
        }
      }
      else
      {
#ifdef MORESTABLE
        a = -fm * min;
        b = -fm * max;
        if (a<=b)
         IM_err (IMERR_LOGDIFF, " pos 3 in integrate_migration_term() cm %d, fm %lf, max %lf, min %lf, a %lf  b %lf",cm,fm,max,min,a,b); 
        LogDiff (c, a, b);
        p = c - log (fm);
#else
        p = log ((exp (-fm * min) - exp (-fm * max)) / fm);
#endif /* MORESTABLE */
      }
    }
    else
    {
      p = log (max - min);
    }
  }

//      assert(p > -1e200 && p < 1e200);
  return p;
}                               /* integrate_migration_term */

double
integrate_migration_term_expo_prior (int cm, double fm, double exmean)
{
  double p;
    p = -log (exmean) +( -(cm + 1) * log (fm + 1.0/exmean)) + logfact[cm];
  return p;
}                               /* integrate_migration_term_expo_prior */

double
integrate_splitrate_term (double f, double max, double min)
{
  if (min > 0)
  {
    return (2 - npops) * log (f) + logfact[npops] +
      log (exp (uppergamma (npops - 2, f / max)) -
           exp (uppergamma (npops - 2, f / min)));
  }
  else
  {
    return (2 - npops) * log (f) + logfact[npops] + uppergamma (npops - 2,
                                                                f / max);
  }
}                               //integrate_splitrate_term




/* simulate the migration path along the moved edge
this gets called starting from the top, period by period  
'period' is the current period,  numm is the number of migration events so far on this edge 
if it must end up in a particular population,  then constrainpop is that population */
int
simmpath (int ci, struct edgemiginfo *edgem, int period, int numm, int lastm,
          double timein, double upt, int pop, int constrainpop)
{
  int i, lastpop, startm;
  struct edge *gtree;
                    /* CR 110715.1 */
  int dupCheck;     /* flag turns off dup migration time check of mig events */ 
  int migIndex;     /* index used to look for duplicate migration times */  

  assert (numm > 0);
  startm = lastm + 1;
  lastm = lastm + numm;
  gtree = C[ci]->G[edgem->li].gtree;
  checkmig (lastm + 1, &gtree[edgem->edgeid].mig, &gtree[edgem->edgeid].cmm);

  do
  { 
    for (i = startm; i <= lastm; i++)
        edgem->mig[i].mt = upt + uniform () * timein;
    edgem->mig[i].mt = -1;
    dupCheck=0;     
    if (numm > 1)
    {
      hpsortmig (&edgem->mig[startm] - 1, numm);
        
      /* CR 110715.1
       * look for duplicate migration times in sorted event list.
       * This solves a charateristic of the Mersennes Twister 
       * random number generator in which identical random numbers may be
       * returned from the random number sequence in a very small number
       * of calls.  With some seeds it was noted as small as within 4 calls.
       */
      for (migIndex = startm; migIndex < lastm; ++migIndex)
      {
        if  (edgem->mig[migIndex].mt != edgem->mig[migIndex + 1].mt)
        {   
           continue;
        }
        else
        {  /* if duplicate time found, a new migration path must be simulated */
           dupCheck=1;
           break;
        }
      }
    }
    else
    {   
      dupCheck=0;
    }
  } while (dupCheck == 1); /* when no dup times found, exit loop   */

  lastpop = pop;
  for (i = startm; i <= lastm; i++)
  {
    if (constrainpop >= 0 && i >= lastm-1) 
    {
      if (i==lastm-1)
        edgem->mig[i].mp =  picktopop2 (lastpop, C[ci]->plist[period], npops - period,
                  constrainpop);
      else //i==lastm
        edgem->mig[i].mp = constrainpop;
    }
    else
      edgem->mig[i].mp = picktopop (lastpop, C[ci]->plist[period], npops - period);

    lastpop = edgem->mig[i].mp;
  }
  return lastm;
} // simmpath

/*
  lastpop = pop;
  if (constrainpop < 0)
  {
    for (i = startm; i <= lastm; i++)
    {
      edgem->mig[i].mp =
        picktopop (lastpop, C[ci]->plist[period], npops - period);
      lastpop = edgem->mig[i].mp;
    }
  }
  else 
  {
    if (numm >= 2)
    {
      i = startm;
      while (i < lastm - 1)
      {
        edgem->mig[i].mp =
          picktopop (lastpop, C[ci]->plist[period], npops - period);
        lastpop = edgem->mig[i].mp;
        i++;
      }
      edgem->mig[lastm - 1].mp =
        picktopop2 (lastpop, C[ci]->plist[period], npops - period,
                    constrainpop);
      edgem->mig[lastm].mp = constrainpop;
    }
    else
    {
      if (numm == 1)
        edgem->mig[lastm].mp = constrainpop;
    }
  }
  return lastm;  */
//}                               /* simmpath */



/********** FUNCTIONS DECLARED IN UPDATE_GTREE_COMMON.c***********/


/* calcmrate()
the count of migration events divided by the length of the relevant portion of the genealogy yields the 
migration rate on the genealogy, per unit time 

the migration rate prior is in units of migration events per mutation event 
to convert this to a relevant prior for the locus,  need to get the mutation rate scalar for this locus

8/08  removed use of mutation rate scalar in this
had thought it made sense for a given locus,  but now think not - essentially a bug
*/

/* 9/25/08  updated this
revised genealogy updating so that the current migration rate is based on the current number of migration events and the 
current length of the branch that is being updated 
reasoned that this might work better than using the rate that occurs for the entitre genealogy - e.g. help to avoid promoting
correlations and improve mixing  */

/* 7/21/10  after a great deal of messing with this,  settled on this function
also, calls to this function should be based on single edges,  not estiamtes of migration over the entire tree */ 

double
calcmrate (int mc, double mt)
{
  assert (mc >= 0);
  if (mt <= 0.0)
    return 1.0;
  if (mc == 0)
  {
    if (mt < 1)
      return 0.1; 
    else
      return 0.1 / mt;
  }
  else
  {
    if (mt < 1)
      return (double) mc;
    else
      return ((double) mc) / mt;
  }
}                               //calcmrate


/* SANGCHUL: Wed Dec 10 10:42:24 EST 2008
 * We could have used a simple memset function if members of structure
 * edgemiginfo are not dynamically allocated. We have to replace memset
 * function with function IMA_reset_edgemiginfo in order to reset structure
 * edgemiginfo.
 * member li must be initialized once. It should be not be changed in function
 * IMA_reset_edgemiginfo. Global edgemiginfo should have -1 value of li, which
 * is used in function checkmig. */
void
IMA_reset_edgemiginfo (struct edgemiginfo *em)
{
  em->edgeid = -1;
  em->sisid = -1;
  em->b = -1;
  em->e = -1;
  em->pop = -1;
  em->temppop = -1;
  em->fpop = -1;
  em->upt = -1.0;
  em->dnt = -1.0;
  assert (em->mtimeavail != NULL);
  memset (em->mtimeavail, 0, npops * sizeof (double));
  em->mtall = 0.0;              /* why this was -1.0 */
  assert (em->mp != NULL);
  memset (em->mp, 0, npops * sizeof (int));
  em->mpall = 0;                /* why this was -1 */
  assert (em->mig != NULL);
  em->mig[0].mt = -1.0;
  return;
}


void
init_gtreecommon (void)         // initialize copy edge
{
  int i, li, largestsamp;
  copyedge = static_cast<struct edge *> (malloc (3 * (sizeof (struct edge))));
  for (i = 0; i < 3; i++)
  {         // very wasteful of space - should use dynamic memory and checkmig()
    copyedge[i].mig = static_cast<struct migstruct *> 
                    (malloc (ABSMIGMAX * sizeof (struct migstruct)));   
  }
  if (somestepwise)
    for (i = 0; i < 3; i++)
    {
      copyedge[i].A = static_cast<int *> (calloc (MAXLINKED, sizeof (int)));
      copyedge[i].dlikeA = static_cast<double *> 
                            (calloc (MAXLINKED, sizeof (double)));
    }
  for (i = 0, li = 0; li < nloci; li++)
    if (i < L[li].numgenes)
      i = L[li].numgenes;
  largestsamp = i;
  nnminus1 = static_cast<double *> 
             (malloc ((largestsamp + 1) * sizeof (double)));
  nnminus1[0] = 0;
  for (i = 1; i <= largestsamp; i++)
  {
    nnminus1[i] = (double) (i) * ((double) i - 1);
  }
  IMA_initmemory_edgemiginfo (&newedgemig);
  IMA_initmemory_edgemiginfo (&newsismig);
  IMA_initmemory_edgemiginfo (&oldedgemig);
  IMA_initmemory_edgemiginfo (&oldsismig);
}                               /* init_gtreecommon */


void
free_gtreecommon (void)
{
  int i;
  for (i = 0; i < 3; i++)
  {
    XFREE (copyedge[i].mig);
    if (somestepwise)
    {
      XFREE (copyedge[i].A);
      XFREE (copyedge[i].dlikeA);
    }
  }
  XFREE (copyedge);
  XFREE (nnminus1);
  return;
}                               //free_gtreecommon


void
init_holdgtree (struct genealogy *g, int numgenes)
{
  int i;
  int numlines = 2 * numgenes - 1;
  g->gtree = static_cast<struct edge *> (calloc ((size_t) numlines, (sizeof (struct edge))));
  for (i = 0; i < numlines; i++)
  {
    g->gtree[i].mig = static_cast<struct migstruct *> (malloc (MIGINC * sizeof (struct migstruct)));
    g->gtree[i].mig[0].mt = -1;
    g->gtree[i].cmm = MIGINC;
    g->gtree[i].up[0] = -1;
    g->gtree[i].up[1] = -1;
    g->gtree[i].down = -1;
    g->gtree[i].time = 0;
    g->gtree[i].mut = -1;
    g->gtree[i].pop = -1;
    if (somestepwise)
    {
      g->gtree[i].A = static_cast<int *> (calloc (MAXLINKED, sizeof (int)));
      g->gtree[i].dlikeA = static_cast<double *> 
                            (calloc (MAXLINKED, sizeof (double)));
    }
  }
  init_genealogy_weights (&(g->gweight));
}                               //init_holdgtree

void
free_holdgtree (struct genealogy *g, int numgenes)
{
  int i, numlines = 2 * numgenes - 1;
  for (i = 0; i < numlines; i++)
  {
    XFREE (g->gtree[i].mig);
    if (somestepwise)
    {
      XFREE (g->gtree[i].A);
      XFREE (g->gtree[i].dlikeA);
    }
  }
  free_genealogy_weights (&g->gweight);
  XFREE (g->gtree);
}                               //free_holdgtree

/* FIXME: DELETE VARIABLES: lmedgedrop, medgedrop!? */
void
storeoldedges (int ci, int li, int edge, int sisedge, int downedge)
{
  int i;
  double uptime;
  struct edge *gtree = C[ci]->G[li].gtree;
  copyedge[0].down = gtree[edge].down;
  i = -1;

  do
  {
    i++;
    if (!(ABSMIGMAX > i))
    {
      IM_err (IMERR_TOOMANYMIG, "step %d: locus [%d] edge [%d] mig %d > %d",
              step, li, edge, ABSMIGMAX, i);
    }
    copyedge[0].mig[i] = gtree[edge].mig[i];
  } while (copyedge[0].mig[i].mt > -0.5);
  copyedge[0].cmm = gtree[edge].cmm;
  medgedrop = i;
  if (edge < L[li].numgenes)
  {
    uptime = 0;
  }
  else
  {
    uptime = gtree[gtree[edge].up[0]].time;
  }
  if (uptime < C[ci]->tvals[lastperiodnumber - 1])
  {
    if (gtree[edge].time < C[ci]->tvals[lastperiodnumber - 1])
    {
      lmedgedrop = gtree[edge].time - uptime;
    }
    else
    {
      lmedgedrop = C[ci]->tvals[lastperiodnumber - 1] - uptime;
    }
  }
  else
  {
    lmedgedrop = 0;
    assert (medgedrop == 0);
  }
  copyedge[0].time = gtree[edge].time;
  copyedge[0].pop = gtree[edge].pop;
  copyedge[1].down = gtree[sisedge].down;
  i = -1;
  do
  {
    i++;
    if (!(ABSMIGMAX > i))
    {
      IM_err (IMERR_TOOMANYMIG, "step %d: locus [%d] edge [%d] mig %d > %d",
              step, li, sisedge, ABSMIGMAX, i);
    }
    copyedge[1].mig[i] = gtree[sisedge].mig[i];
  } while (copyedge[1].mig[i].mt > -0.5);
  copyedge[1].cmm = gtree[sisedge].cmm;
  copyedge[1].time = gtree[sisedge].time;
  copyedge[1].pop = gtree[sisedge].pop;
  copyedge[2].down = gtree[downedge].down;
  i = -1;
  do
  {
    i++;
    if (!(ABSMIGMAX > i))
    {
      IM_err (IMERR_TOOMANYMIG, "step %d: locus [%d] edge [%d] mig %d > %d",
              step, li, downedge, ABSMIGMAX, i);
    }
    copyedge[2].mig[i] = gtree[downedge].mig[i];
  } while (copyedge[2].mig[i].mt > -0.5);
  copyedge[2].cmm = gtree[downedge].cmm;
  copyedge[2].time = gtree[downedge].time;
  copyedge[2].pop = gtree[downedge].pop;
  copyedge[0].up[0] = gtree[edge].up[0];
  copyedge[0].up[1] = gtree[edge].up[1];
  copyedge[1].up[0] = gtree[sisedge].up[0];
  copyedge[1].up[1] = gtree[sisedge].up[1];
  copyedge[2].up[0] = gtree[downedge].up[0];
  copyedge[2].up[1] = gtree[downedge].up[1];
  if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
    storeAinfo (li, gtree, edge, sisedge, downedge);
}                               /* storoldeges */


/* set the gtree back to the way it was */
void
restoreedges (int ci, int li, int edge, int sisedge, int downedge,
              int newsisedge)
/*all this can be optimized some*/
{
  int i, j, ai, down;
  struct edge *gtree = C[ci]->G[li].gtree;
  if (newsisedge != sisedge)
  {
    down = gtree[downedge].down;
    if (down != -1)
    {
      if (gtree[down].up[0] == downedge)
        gtree[down].up[0] = newsisedge;
      else
        gtree[down].up[1] = newsisedge;
    }
    else
    {
      C[ci]->G[li].root = newsisedge;
      C[ci]->G[li].roottime = gtree[gtree[newsisedge].up[0]].time;
      assert (C[ci]->G[li].roottime <= TIMEMAX);
    }
    gtree[newsisedge].down = down;
    if (down != -1)
    {
      i = 0;
      while (gtree[newsisedge].mig[i].mt > -0.5)
        i++;
      j = -1;

      do
      {
        j++;
        checkmig (i + j, &(gtree[newsisedge].mig), &(gtree[newsisedge].cmm));
        gtree[newsisedge].mig[i + j] = gtree[downedge].mig[j];
      } while (gtree[downedge].mig[j].mt > -0.5);
    }
    else
    {
      gtree[newsisedge].mig[0].mt = -1;
      if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
        for (ai = (L[li].model == JOINT_IS_SW); ai < L[li].nlinked; ai++)
          gtree[newsisedge].dlikeA[ai] = 0;
    }
    gtree[newsisedge].time = gtree[downedge].time;
  }
  gtree[edge].down = copyedge[0].down;
  i = -1;

  do
  {
    i++;
    checkmig (i, &(gtree[edge].mig), &(gtree[edge].cmm));
    gtree[edge].mig[i] = copyedge[0].mig[i];
  } while (gtree[edge].mig[i].mt > -0.5);
  gtree[edge].time = copyedge[0].time;
  gtree[edge].pop = copyedge[0].pop;
  down = gtree[sisedge].down;
  gtree[sisedge].down = copyedge[1].down;
  if (down != -1)
  {
    if (gtree[down].up[0] == sisedge)
      gtree[down].up[0] = downedge;

    else
      gtree[down].up[1] = downedge;
  }
  i = -1;

  do
  {
    i++;
    checkmig (i, &(gtree[sisedge].mig), &(gtree[sisedge].cmm));
    gtree[sisedge].mig[i] = copyedge[1].mig[i];
  } while (gtree[sisedge].mig[i].mt > -0.5);
  gtree[sisedge].time = copyedge[1].time;
  gtree[sisedge].pop = copyedge[1].pop;
  gtree[downedge].down = copyedge[2].down;
  i = -1;

  do
  {
    i++;
    checkmig (i, &(gtree[downedge].mig), &(gtree[downedge].cmm));
    gtree[downedge].mig[i] = copyedge[2].mig[i];
  } while (gtree[downedge].mig[i].mt > -0.5);
  gtree[downedge].time = copyedge[2].time;
  gtree[downedge].pop = copyedge[2].pop;
  gtree[downedge].up[0] = copyedge[2].up[0];
  gtree[downedge].up[1] = copyedge[2].up[1];
  if (gtree[downedge].down == -1)
  {
    C[ci]->G[li].roottime = gtree[gtree[downedge].up[0]].time;
    assert (C[ci]->G[li].roottime <= TIMEMAX);
    C[ci]->G[li].root = downedge;
  }
  if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
    for (ai = (L[li].model == JOINT_IS_SW); ai < L[li].nlinked; ai++)
    {
      gtree[edge].A[ai] = copyedge[0].A[ai];
      gtree[sisedge].A[ai] = copyedge[1].A[ai];
      gtree[downedge].A[ai] = copyedge[2].A[ai];
      gtree[edge].dlikeA[ai] = copyedge[0].dlikeA[ai];
      gtree[sisedge].dlikeA[ai] = copyedge[1].dlikeA[ai];
      gtree[downedge].dlikeA[ai] = copyedge[2].dlikeA[ai];
      if (holdsisdlikeA[ai] != 0)
        gtree[newsisedge].dlikeA[ai] = holdsisdlikeA[ai];
    }
}                               /* restoreedges  */

// changes made 5/27/08  so that this could be called directly from update_genealogy_covar() 
// 12/17/09 JH  is this now not necessary now that update_genealogy_covar() is gone ?? 
/* info for edge goes in copyedge[0],
   info for sisedge goes in copyedge[1]
   info for downedge goes in copyedge[2] */
void
storeAinfo (int li, struct edge *gtree, int edge, int sisedge,
            int downedge)
{
  int ai;

  for (ai = (L[li].model == JOINT_IS_SW); ai < L[li].nlinked; ai++)
  {
    copyedge[0].A[ai] = gtree[edge].A[ai];
    copyedge[0].dlikeA[ai] = gtree[edge].dlikeA[ai];
    copyedge[1].A[ai] = gtree[sisedge].A[ai];
    copyedge[1].dlikeA[ai] = gtree[sisedge].dlikeA[ai];
    copyedge[2].A[ai] = gtree[downedge].A[ai];
    if (gtree[downedge].down != -1)   // added this 5/27/08 
    {
      copyedge[2].dlikeA[ai] = gtree[downedge].dlikeA[ai];
      holddownA[ai] = gtree[gtree[downedge].down].A[ai];      // changed this 5/27/08
    }
    else
    {
      copyedge[2].down = -1;  // inserted this 5/27/08
      holddownA[ai] = -1;
      copyedge[2].dlikeA[ai] = 0;
    }
  }
  return;
}


double
getmprob(int ci, struct edgemiginfo *edgem,
          struct edgemiginfo *sisem,struct edgemiginfo *oldedgem,
          struct edgemiginfo *oldsisem)
{
  double tempp, n, d, t;
  double r;
  int periodi, cm[2], pop[2], topop, popc;
  int lastm_2_pop;
  double pathc;
  int ii, lastmigrationperiod;
  struct edgemiginfo *mm, *oldmm;

  tempp = 0;
  lastmigrationperiod = IMIN(edgem->e, lastperiodnumber-1);
  if (sisem->mtall <= 0)       // only deal with edgem 
  {
    cm[0] = 0;
    //for (periodi = edgem->b; periodi <= edgem->e - 1; periodi++)
    for (periodi = edgem->b; periodi <= edgem->e ; periodi++)
      if ((periodi < lastmigrationperiod) || (periodi == lastmigrationperiod && edgem->e == lastperiodnumber))
      {
        assert (edgem->mtimeavail[periodi] > 0.0);
        r = calcmrate(oldedgem->mp[periodi],oldedgem->mtimeavail[periodi])* edgem->mtimeavail[periodi];//9/2/10
        tempp += edgem->mp[periodi] * log(r/(edgem->mtimeavail[periodi] * (npops - (periodi + 1)))) - r;
        assert (tempp != 0);
        cm[0] += edgem->mp[periodi];

      }
    if (edgem->e < lastperiodnumber)
    {

      /* if assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1  do not enter here */
      r = calcmrate(oldedgem->mp[edgem->e],oldedgem->mtimeavail[edgem->e])* edgem->mtimeavail[edgem->e];//9/2/10
      if (edgem->e == lastperiodnumber - 1 
          && (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0 || npops == 2))     // only 2 pops in last period
      {
        if (ODD (edgem->mp[edgem->e]))
        {
          tempp += edgem->mp[edgem->e] * log(r/edgem->mtimeavail[edgem->e]) - mylogsinh (r);
        }
        else
        {
          tempp += edgem->mp[edgem->e] * log(r/edgem->mtimeavail[edgem->e]) - mylogcosh(r);
        }
      }
      else                      // 3 or more pops in the last period 
      {
        if (cm[0] == 0)
          pop[0] = edgem->pop;
        else
          pop[0] = edgem->mig[cm[0] - 1].mp;  // the last migration before this last period 
        while (C[ci]->poptree[pop[0]].e <= edgem->e)
          pop[0] = C[ci]->poptree[pop[0]].down;
        topop = edgem->fpop;
        popc = (npops - edgem->e - 1);
        if (pop[0] == topop)
        {
            d = log (1 - r * exp (-r));
        }
        else           
        {
          d = log ((1 -  exp (-r)));
        }
        switch (edgem->mp[edgem->e])
        {
        case 0: assert(pop[0]==topop);
                n = -r; 
                break;
        case 1: assert(pop[0] != topop);
                n = log(r/edgem->mtimeavail[edgem->e]) - r;
                break;
        default:
                {
                  assert(edgem->mp[edgem->e] >=2);
                  if (edgem->mp[edgem->e]==2)
                    lastm_2_pop = pop[0];  
                  else
                    lastm_2_pop = edgem->mig[edgem->mpall-3].mp;  
/*assert(checkpop(ci,edgem->e, lastm_2_pop)); */
                  if (lastm_2_pop == topop)
                    pathc = -log((double) popc); 
                  else
                    pathc = -log((double) popc-1); 
                  n = edgem->mp[edgem->e] * log(static_cast<double>(r)/edgem->mtimeavail[edgem->e]) - r + (2- edgem->mp[edgem->e]) * log(static_cast<double>(popc)) + pathc;
                }
        }
        tempp += n - d;
      }
    }
  }
  else  // both edges 
  {
    if (edgem->mtall > 0 && sisem->mtall > 0 && edgem->e < lastperiodnumber)
    {
      for (ii=0;ii<2;ii++)  // probability of simulating fpop
      {
        mm = (ii==0)? edgem : sisem;
        oldmm = (ii==0)? oldedgem : oldsisem;
        pop[ii] = mm->pop;
        if (mm->mpall > 0 && mm->e > 0)
        {
          t = C[ci]->tvals[mm->e-1];
          periodi= -1;
          while (mm->mig[periodi+1].mt >= 0 && mm->mig[periodi+1].mt < t)
          {
            periodi++;
          }
          if (periodi >= 0)
            pop[ii] = mm->mig[periodi].mp;
        }
        if (mm->e > 0)
          while (C[ci]->poptree[pop[ii]].e <= mm->e)
            pop[ii] = C[ci]->poptree[pop[ii]].down;
      }
      if (pop[0] == pop[1])
      {
        if (pop[0] == edgem->fpop)
          tempp = log(MIGCLOSEFRAC);
        else
          tempp = log((1.0-MIGCLOSEFRAC) / (double) (npops-edgem->e-1));
      }
      else
      {
        if (npops - edgem->e == 2)
          tempp = log(0.5);
        else
        {
          if (edgem->fpop == pop[0] || edgem->fpop == pop[1])
            tempp = log(0.5* MIGCLOSEFRAC);
          else
            tempp = log((1.0 - MIGCLOSEFRAC)/(double) (npops-edgem->e-2));
        }
      }
    } 
    else
      tempp = 0.0;
    lastmigrationperiod = IMIN(edgem->e, lastperiodnumber-1);
    for (ii = 0;ii<2;ii++)
    {
      mm = (ii==0)? edgem : sisem;
      oldmm = (ii==0)? oldedgem : oldsisem;
      if (mm->mtall > 0)
      {
        cm[ii] = 0;
        for (periodi = mm->b; periodi <= mm->e ; periodi++)
          if ((periodi < lastmigrationperiod) ||
              (periodi == lastmigrationperiod && mm->e == lastperiodnumber))
        {
          assert (mm->mtimeavail[periodi] > 0.0);
          r = calcmrate(oldmm->mp[periodi],oldmm->mtimeavail[periodi])* mm->mtimeavail[periodi];//9/2/10
          tempp += mm->mp[periodi] * log(r/(mm->mtimeavail[periodi] * (npops - (periodi + 1)))) - r;
          assert (tempp != 0);
          cm[ii] += mm->mp[periodi];

        }

        if (mm->e < lastperiodnumber)
        {
          r = calcmrate(oldmm->mp[mm->e],oldmm->mtimeavail[mm->e])* mm->mtimeavail[mm->e];//9/2/10
          /* if assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1  do not enter here */
          if (mm->e == lastperiodnumber - 1 
              && (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0 || npops == 2))     // only 2 pops in last period
          {
            assert (mm->mtimeavail[mm->e] > 0);
            if (ODD (mm->mp[mm->e]))
            {
              tempp += mm->mp[mm->e] * log(r/mm->mtimeavail[mm->e]) - mylogsinh (r);
            }
            else
            {
              tempp += mm->mp[mm->e] * log(r/mm->mtimeavail[mm->e]) - mylogcosh(r);
            }
          }
          else                      // 3 or more pops in the last period 
          {
            if (cm[ii] == 0)
              pop[ii] = mm->pop;
            else
              pop[ii] = mm->mig[cm[ii] - 1].mp;  // the last migration before this last period 
            while (C[ci]->poptree[pop[ii]].e <= mm->e)
              pop[ii] = C[ci]->poptree[pop[ii]].down;
            topop = mm->fpop;
            popc = (npops - mm->e - 1);
            if (pop[ii] == topop)
            {
              d = log (1 - r * exp (-r));
            }
            else           
            {
              d = log ((1 -  exp (-r)));
            }
            switch (mm->mp[mm->e])
            {
            case 0: assert(pop[ii]==topop);
                    n = -r;
                    break;
            case 1: assert(pop[ii] != topop);
                    n = log(r/mm->mtimeavail[mm->e]) - r; 
                    break;
            default:
                    assert(mm->mp[mm->e] >=2);
                    if (mm->mp[mm->e]==2)
                      lastm_2_pop = pop[ii];  
                    else
                      lastm_2_pop = mm->mig[mm->mpall-3].mp;  
                    if (lastm_2_pop == topop)
                      pathc = -log((double) popc); 
                    else
                      pathc = -log((double) popc-1); 
                    n = mm->mp[mm->e] * log(static_cast<double>(r)/mm->mtimeavail[mm->e]) - r + (2- mm->mp[mm->e]) * log(static_cast<double>(popc)) + pathc;
                    break;
              }
            tempp += n - d;
          }
        }
      }
    }
  }
  assert (tempp > -1e200 && tempp < 1e200);
  return tempp;
}                               /* getmprob*/



/*void summiginfo(struct edgemiginfo *e, struct edgemiginfo *s)
{
  struct migstruct mig[ABSMIGMAX];
  int i,j;
  s->upt += e->upt;
  s->dnt += e->dnt;
  s->fpop += e->fpop;
  s->b += e->b;
  s->e += e->e;
  s->mpall += e->mpall;
  s->mtall += e->mtall;
  for (i = e->b; i <= e->e; i++)
  {
    s->mtimeavail[i] +=    e->mtimeavail[i];
    for (j=0;j<10;j++)
    {
      s->mig[j].mt += e->mig[j].mt;
      s->mig[j].mp += e->mig[j].mp;
    }
  }
}*/

/* fillmiginfoperiods() */
/* called by fillmiginfo  and by addmigration.  Fills in info on time periods spanned by an edge*/
/* b and e are the time periods during with the beginning and ends of the edge occur, respectively */
void
fillmiginfoperiods (int ci, struct edgemiginfo *em)
{
  int i;
  //double tcheck = 0; debug
  em->b = 0;
  while (em->upt > C[ci]->tvals[em->b])
    em->b++;
  em->e = em->b;
  while (em->dnt > C[ci]->tvals[em->e])
    em->e++;
  if (em->e == em->b)
  {
    if (em->b == lastperiodnumber)
      em->mtimeavail[em->b] = 0;
    else
    {
      em->mtimeavail[em->b] = em->dnt - em->upt;
      assert(em->mtimeavail[em->b] > 0);
    }
  }
  else
  {
    em->mtimeavail[em->b] = C[ci]->tvals[em->b] - em->upt;
    assert(em->mtimeavail[em->b] > 0);
    if (em->e == lastperiodnumber)
      em->mtimeavail[em->e] = 0;
    else
    {
      em->mtimeavail[em->e] = em->dnt - C[ci]->tvals[em->e - 1];
      assert(em->mtimeavail[em->e] > 0);
    }
    for (i = em->b + 1; i < em->e; i++)
    {
      em->mtimeavail[i] = C[ci]->tvals[i] - C[ci]->tvals[i - 1];
      assert(em->mtimeavail[i] > 0);
    }
  }
  if (em->b < lastperiodnumber)
  {
    em->mtall = DMIN (C[ci]->tvals[lastperiodnumber - 1], em->dnt) - em->upt;
    /*for (i = em->b; i <= em->e; i++)  //debugging
    {
      if (i < lastperiodnumber)
      {
        assert(em->mtimeavail[i] > 0);
    //    tcheck += em->mtimeavail[i]; debug
      }
    } */
    //assert(fabs(tcheck - em->mtall) < 1e-14); debug
  }
  else
  {
    em->mtall = 0;
  }


  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
  {
    assert (em->b == 0 && em->e == 0);
  }
  return;
}                               /*fillmiginfoperiods */

/* fillmiginfo() 
for updating edges and migration it is necessary to gather the info about the original edge 
and the original sister edge into one place to make sure it is easily available 
initialize oldedgemig and oldsismig*/

void
fillmiginfo (int ci, int li, struct edge *gtree, int edge, int sisedge)
{
  int i, j;
  double uptime;

/* REMOVED!
  memset (&oldedgemig, 0, sizeof (struct edgemiginfo));
  memset (&oldsismig, 0, sizeof (struct edgemiginfo));
*/
  IMA_reset_edgemiginfo (&oldedgemig);
  IMA_reset_edgemiginfo (&oldsismig);

  oldedgemig.edgeid = edge;
  oldsismig.edgeid = sisedge;
  oldedgemig.li = li;
  if (edge < L[li].numgenes)
    uptime = 0;
  else
    uptime = gtree[gtree[edge].up[0]].time;
  oldedgemig.pop = gtree[edge].pop;
  oldedgemig.dnt = gtree[edge].time;
  oldedgemig.upt = uptime;
  oldedgemig.fpop = gtree[gtree[edge].down].pop;
  fillmiginfoperiods (ci, &oldedgemig);
  oldedgemig.mig[0].mt = -1;
  for (j = 0; j <= oldedgemig.b; j++)
    oldedgemig.mp[j] = 0;
  oldedgemig.mpall = 0;
  i = 0;
  j = oldedgemig.b;
  while (gtree[edge].mig[i].mt > -0.5)
  {
    while (gtree[edge].mig[i].mt > C[ci]->tvals[j])
    {
      j++;
      oldedgemig.mp[j] = 0;
    }
    oldedgemig.mp[j]++;
    oldedgemig.mpall++;
    oldedgemig.mig[i] = gtree[edge].mig[i];
    i++;
  }
  oldedgemig.mig[i].mt = -1;
  if (sisedge >= 0)
  {
    oldsismig.li= li;
    // assert (gtree[sisedge].down == C[ci]->G[li].root); not true if this is called from t updating 
    if (sisedge < L[li].numgenes)
      uptime = 0;
    else
      uptime = gtree[gtree[sisedge].up[0]].time;
    oldsismig.pop = gtree[sisedge].pop;
    oldsismig.dnt = gtree[sisedge].time;
    oldsismig.upt = uptime;
    oldsismig.fpop = gtree[gtree[sisedge].down].pop;
    fillmiginfoperiods (ci, &oldsismig);
    oldsismig.mig[0].mt = -1;
    for (j = 0; j <= oldsismig.b; j++)
      oldsismig.mp[j] = 0;
    oldsismig.mpall = 0;
    i = 0;
    j = oldsismig.b;
    while (gtree[sisedge].mig[i].mt > -0.5)
    {
      while (gtree[sisedge].mig[i].mt > C[ci]->tvals[j])
      {
        j++;
        oldsismig.mp[j] = 0;
      }
      oldsismig.mp[j]++;
      oldsismig.mpall++;
      oldsismig.mig[i] = gtree[sisedge].mig[i];
      i++;
    }
    oldsismig.mig[i].mt = -1;
  }
}                               /* fillmiginfo */


// copy the migration info in newedgemig and newsismig  to the genealogy
void
copynewmig_to_gtree (int ci, int li)
{
  int i, pop, downperiod;
  struct edge *gtree = C[ci]->G[li].gtree;
  i = 0;
  while (newedgemig.mig[i].mt > 0)
  {
    checkmig (i, 
              &(gtree[newedgemig.edgeid].mig), 
              &(gtree[newedgemig.edgeid].cmm));
    gtree[newedgemig.edgeid].mig[i] = newedgemig.mig[i];
    i++;
  }
  gtree[newedgemig.edgeid].mig[i].mt = -1;
  if (newsismig.edgeid >= 0)
  {
    i = 0;

    //pop = newsismig.pop;
    while (newsismig.mig[i].mt > 0)
    {
      checkmig (i, &(gtree[newsismig.edgeid].mig),
                &(gtree[newsismig.edgeid].cmm));
      gtree[newsismig.edgeid].mig[i] = newsismig.mig[i];
      i++;
    }
    gtree[newsismig.edgeid].mig[i].mt = -1;

    // set top of down edge - seems to be necessary when root moves 
    if (i > 0)
      pop = gtree[newsismig.edgeid].mig[i - 1].mp;
    else
      pop = gtree[newsismig.edgeid].pop;
    downperiod = findperiod (ci, gtree[newsismig.edgeid].time);
    while (C[ci]->poptree[pop].e <= downperiod && C[ci]->poptree[pop].e != -1)
      pop = C[ci]->poptree[pop].down;
    gtree[gtree[newsismig.edgeid].down].pop = pop;
  }
}                               /* copynewmig_to_gtree */


/* store a few basic tree statistics */
void
storegenealogystats (int ci, int li, int mode)
{
  static double holdlength, holdtlength;
  static double holdroottime;
  static int holdroot;
  static int holdmig;
  if (mode == 0)
  {
    holdlength = C[ci]->G[li].length;
    holdtlength = C[ci]->G[li].tlength;
    holdroottime = C[ci]->G[li].roottime;
    holdroot = C[ci]->G[li].root;
    holdmig = C[ci]->G[li].mignum;
  }
  else
  {
    C[ci]->G[li].length = holdlength;
    C[ci]->G[li].tlength = holdtlength;
    C[ci]->G[li].mignum = holdmig;
    C[ci]->G[li].roottime = holdroottime;
    C[ci]->G[li].root = holdroot;
  }
}                               /* storegenealogystats */


/* return the population to which migration happens, chosen at random from those in plist*/
int
picktopop (int nowpop, int plist[], int numpops)
{
  int topop;

  do
  {
    topop = randposint (numpops);
  } while (plist[topop] == nowpop);

  assert (topop < numpops);
  return plist[topop];
}                               // picktopop


/* return the population to which migration happens, chosen at random from those in plist*/
/* this is same as picktopop() except that this will specifically avoid population notother */
int
picktopop2 (int nowpop, int plist[], int numpops, int notother)
{
  int topop;
  do
  {
    topop = randposint (numpops);
  } while (plist[topop] == nowpop || plist[topop] == notother);

  assert (topop < numpops);
  return plist[topop];
}                               //picktopop2



/* adds migration for single edges,  updates temppop (population id of the edge at the time under consideration) as needed */
int
mwork_single_edge (int ci, struct edgemiginfo *edgem,struct edgemiginfo *oldedgem, int lastmigperiod)  
{
  int periodi, lastm, mpall;
  double r;
  double timestartperiod;


  assert (lastmigperiod < lastperiodnumber);
  if (lastmigperiod < edgem->b)
    return 0;
  lastm = -1;
  mpall = 0;
  timestartperiod = edgem->upt;
  for (periodi = edgem->b;periodi <= lastmigperiod; periodi++)
  {
    assert(periodi < lastperiodnumber);
    while (C[ci]->poptree[edgem->temppop].e <= periodi && C[ci]->poptree[edgem->temppop].e != -1)
      edgem->temppop = C[ci]->poptree[edgem->temppop].down;
    r = calcmrate(oldedgem->mp[periodi],oldedgem->mtimeavail[periodi])*edgem->mtimeavail[periodi]; //9/2/10
    if (periodi < edgem->e)           // can do migration to any available pop
    {
      if ((edgem->mp[periodi] = poisson (r, -1)) > 0)
        lastm =
          simmpath (ci, edgem, periodi, edgem->mp[periodi], lastm,
                    edgem->mtimeavail[periodi], timestartperiod, edgem->temppop,-1);
    }
    else
    {
      if (npops - periodi == 2)
      {
        if (edgem->temppop == edgem->fpop)      //even
        {
          edgem->mp[periodi] = poisson (r, 0);
        }
        else                    // odd
        {
          edgem->mp[periodi] = poisson (r, 1);
        }
        if (edgem->mp[periodi] > 0)
          lastm =
            simmpath (ci, edgem, periodi, edgem->mp[periodi], lastm,
                      edgem->mtimeavail[periodi], timestartperiod, edgem->temppop, -1);
      }
      else
      {
        if (edgem->temppop == edgem->fpop)      //cannot be just one migration event
        {
          edgem->mp[periodi] = poisson (r, 3);
        }
        else                    // cannot be zero migration events 
        {
          edgem->mp[periodi] = poisson (r, 2);
        }
        if (edgem->mp[periodi] > 0)
          lastm =
            simmpath (ci, edgem, periodi, edgem->mp[periodi], lastm,
                      edgem->mtimeavail[periodi], timestartperiod,
                      edgem->temppop, edgem->fpop);
      }
    }
    if (edgem->mp[periodi] > 0)
    {
      mpall += edgem->mp[periodi];
      edgem->temppop = edgem->mig[lastm].mp;
    }
    timestartperiod = C[ci]->tvals[periodi];
  }
  if (edgem->mtimeavail[periodi] > 0)  
  {
    if (C[ci]->poptree[edgem->temppop].e == periodi)
      edgem->temppop = C[ci]->poptree[edgem->temppop].down;
  }
  return mpall;
}                               /* mwork_single_edge */


void
mwork_two_edges(int ci, struct edgemiginfo *edgem, struct edgemiginfo *sisem, struct edgemiginfo *oldedgem, 
    struct edgemiginfo *oldsisem,int lastmigperiod, int* empall, int* smpall) 
{
  int periodi,lastperiodi, lastm[2], mpall[2];
  int b, ii;
  double r;
  double timestartperiod[2];
  struct edgemiginfo *mm, *oldmm;

  assert (lastmigperiod < lastperiodnumber);
  assert(lastmigperiod >= edgem->b);
  lastm[0] = lastm[1] = -1;
  mpall[0] = mpall[1] = 0;
  timestartperiod[0] = edgem->upt;
  timestartperiod[1] = sisem->upt;
  b = IMIN(edgem->b,sisem->b);
  /*  two portions of both edges that come before lastmigperiod */
  if (edgem->e == lastperiodnumber)
    lastperiodi = lastmigperiod;
  else
    lastperiodi = lastmigperiod - 1;
  for (periodi = b; periodi <= lastperiodi; periodi++)
  {
    assert(periodi < lastperiodnumber);
    for (ii=0;ii<2;ii++)
    {
      mm = (ii==0) ? edgem : sisem;
      oldmm = (ii==0) ? oldedgem : oldsisem;
      if (mm->b <= periodi)
      {
        while (C[ci]->poptree[mm->temppop].e <= periodi && C[ci]->poptree[mm->temppop].e != -1)
          mm->temppop = C[ci]->poptree[mm->temppop].down;
        r = calcmrate(oldmm->mp[periodi],oldmm->mtimeavail[periodi])*mm->mtimeavail[periodi];//9/2/10
        if ((mm->mp[periodi] = poisson (r, -1)) > 0)
          lastm[ii] = simmpath (ci, mm, periodi, mm->mp[periodi], lastm[ii],
                     mm->mtimeavail[periodi], timestartperiod[ii], mm->temppop,-1);
        if (mm->mp[periodi] > 0)
        {
          mpall[ii] += mm->mp[periodi];
          mm->temppop = mm->mig[lastm[ii]].mp;
        }
        timestartperiod[ii] = C[ci]->tvals[periodi];
      }
    }
    
  }
  if (periodi == lastperiodnumber)
  {
    edgem->fpop = sisem->fpop = C[ci]->rootpop;
  }
  else
  {
    assert(edgem->e < lastperiodnumber);
  /*  two portions of both edges that are in lastmigperiod */
  /* first determine the population state of the edges where they join */
    if (C[ci]->poptree[edgem->temppop].e == periodi)
      edgem->temppop = C[ci]->poptree[edgem->temppop].down;
    if (C[ci]->poptree[sisem->temppop].e == periodi)
      sisem->temppop = C[ci]->poptree[sisem->temppop].down;
    if (edgem->temppop == sisem->temppop)
    {
      if (uniform() < MIGCLOSEFRAC) 
      {
        edgem->fpop = sisem->fpop = edgem->temppop;
      }
      else
      {
        edgem->fpop = sisem->fpop = picktopop(edgem->temppop,C[ci]->plist[periodi],npops - periodi);
      }
    }
    else
    {
      if (npops - periodi == 2)
      {
        if (uniform() < 0.5) 
          edgem->fpop = sisem->fpop = edgem->temppop;
        else
          edgem->fpop = sisem->fpop = sisem->temppop;
      }
      else
      {
        if (uniform() < MIGCLOSEFRAC) 
        {
          if (uniform() < 0.5) 
            edgem->fpop = sisem->fpop = edgem->temppop;
          else
            edgem->fpop = sisem->fpop = sisem->temppop;
        }
        else // pick something besides edgem->temppop and sisem->temppop
        {
          edgem->fpop = sisem->fpop = picktopop2(edgem->temppop,C[ci]->plist[periodi],npops - periodi,sisem->temppop);
        }
      }
    }
    for (ii=0;ii<2;ii++)
    {
      mm = (ii==0) ? edgem : sisem;
      oldmm = (ii==0) ? oldedgem : oldsisem;
      r = calcmrate(oldmm->mp[periodi],oldmm->mtimeavail[periodi])*mm->mtimeavail[periodi];//9/2/10
      if (npops - periodi == 2)
      {
        if (mm->temppop == mm->fpop)      //even
        {
          mm->mp[periodi] = poisson (r, 0);
        }
        else                    // odd
        {
          mm->mp[periodi] = poisson (r, 1);
        }
        if (mm->mp[periodi] > 0)
          lastm[ii] =
            simmpath (ci, mm, periodi, mm->mp[periodi], lastm[ii],
                      mm->mtimeavail[periodi], timestartperiod[ii], mm->temppop, -1);
      }
      else
      {
        if (mm->temppop == mm->fpop)      //cannot be just one migration event
        {
          mm->mp[periodi] = poisson (r, 3);
        }
        else                    // cannot be zero migration events 
        {
          mm->mp[periodi] = poisson (r, 2);
        }
        if (mm->mp[periodi] > 0)
          lastm[ii] =
            simmpath (ci, mm, periodi, mm->mp[periodi], lastm[ii],
                      mm->mtimeavail[periodi], timestartperiod[ii],
                      mm->temppop, mm->fpop);
      }
      if (mm->mp[periodi] > 0)
      {
        mpall[ii] += mm->mp[periodi];
        mm->temppop = mm->mig[lastm[ii]].mp;
      }
      if (mm->mtimeavail[periodi] > 0)  
      {
        if (C[ci]->poptree[mm->temppop].e == periodi)
          mm->temppop = C[ci]->poptree[mm->temppop].down;
      }
    }
  }
  *empall = mpall[0];
  *smpall = mpall[1];
}                               /* mwork_two_edges*/



/************ GLOBAL FUNCTIONS **************/

// JH 4/1/2011  added these functions that Vitor wrote to trap errors in realloc()
// VS - designed new reallocSgv function to call when sgEvent are reallocated
// this could help to understand if the program crash because of these errors
static void *reallocSgv(struct gtreeevent *sge, int n) {

	// create a new pointer that will save the realloc of the old one
	struct gtreeevent *newSge = static_cast <gtreeevent *> (realloc(sge, n*sizeof(struct gtreeevent)));

	// check that nmig is different from zero
	if(n<1) {
		printf("ERROR with realloc: n is smaller than 1");
	}

	// If the realloc returns null need to free the mig pointer
	if((newSge==NULL) && (sge != NULL)) {
		free(sge); // if the realloc fails - this may happen when nmig is zero!
	}

	return newSge;
}

// VS - designed new reallocSgv function to call when sgEvent_indx are reallocated
// or other unsigned long pointers
static void *reallocSgI(unsigned long *sgeI, int n) {
	
	// create a new pointer that will save the realloc of the old one
	unsigned long *newSgeI = static_cast <unsigned long *>(realloc(sgeI, n*sizeof(unsigned long)));

	// check that nmig is different from zero
	if(n<1) {
		printf("ERROR with reallocSgI: n is smaller than 1");
	}

	// If the realloc returns null need to free the mig pointer
	if((newSgeI==NULL) && (sgeI != NULL)) {
		free(sgeI); // if the realloc fails - this may happen when nmig is zero!
	}

	return newSgeI;
}



/* return the period in the population tree of a particular time point */
int
findperiod (int ci, double t)
{
  int k = 0;
  assert (t >= 0);
  while (k < lastperiodnumber && C[ci]->tvals[k] <= t)
    k++;
  return k;
}

/* return the population that the edge is in at time.  If time is a population splitting time, it returns
the population state in the period above that split*/
int
nowedgepop (int ci, struct edge *gtree, double ptime)
{
  int pop, j;

  pop = gtree->pop;
  j = 0;
  while (gtree->mig[j].mt > -1 && gtree->mig[j].mt < ptime)
  {
    pop = gtree->mig[j].mp;
    j++;
  }
  while (ptime > C[ci]->poptree[pop].time && pop != -1)
  {
    pop = C[ci]->poptree[pop].down;
  }
  assert (pop >= 0 && pop < numtreepops);
  return pop;
}

void
init_treeweight (void)
{
  elength = ESTARTLENGTH;
  sgEvent = static_cast<struct gtreeevent *> 
            (malloc (elength * sizeof (struct gtreeevent)));
  sgEvent_indx = static_cast<unsigned long *> 
                 (malloc (elength * sizeof (unsigned long)));
}                               // init_treeweight

void
free_treeweight (void)
{
  XFREE (sgEvent);
  XFREE (sgEvent_indx);
}                               // free_treeweight


/* updated this on 11/3/08  to avoid most of the mallocing that was being done
created static dynamic arrays  sgEvent and sgEvent_indx  and just realloc them if they are not big enough
Also dropped the use of memset - just not needed */

void
treeweight (int ci, int li)
/* events are migrations, m = 1, or coalescent, c = 0 or population splitting, t = -1
sgEvent contains information for events. event i is contained in e+i and includes:
  the time till the next event
  the event that happens at that time cmt (0 or 1 or -1)
  the population in which the event happens 
if that event is migration, the pop is the one from which the migrant is leaving (going back in time)
the number of individuals in each population, just before the event happens.

accumulates info in G->gweight  
information is accumulated in each period k,  even if a population spans multiple periods. 
The info on coalescent events will have to be summed over periods for those populations that 
span multiple periods.  will have to be summed for population s

G->gweight must have been previously initialized to zero by a call to setzero_genealogy_weights()
*/
{
  int i, ii, j, jj, k, n[2 * MAXPOPS - 1], nsum, ncount;
  int ip, jp;
  int nowpop;
  double t, timeinterval, fmtemp, lasttime;
  int ec;
  double h2term, hlog;
  struct genealogy *G = &(C[ci]->G[li]);
  struct edge *gtree = G->gtree;
  struct genealogy_weights *gweight = &(G->gweight);
  int ng = L[li].numgenes;
  double lastsplitt, timeadd;
  static struct gtreeevent *sgp;
  

  ncount = ng - 1;              // there must be ng-1 coalescent events
  G->mignum = 0;
  for (i = 0; i < 2 * ng - 1; i++)
  {
    j = 0;
    while (gtree[i].mig[j].mt > -0.5)
      j++;
    G->mignum += j;             
  } // count migration events  = must this be done here?  seems necessary but maybe only the first time

  ncount += G->mignum;

  /* add an 'event' for every population tree split time that is younger than the root of the genealogy */
  ncount += findperiod (ci, G->roottime);

  if (ncount + 2 >= elength) // allocate more memory
  {
    sgEvent = static_cast<struct gtreeevent *> 
              (realloc (sgEvent, (ncount + EADDLENGTH) * 
                        sizeof (struct gtreeevent)));
    assert (sgEvent !=  NULL);
    sgEvent_indx = static_cast<unsigned long *> 
                    (realloc (sgEvent_indx, (ncount + EADDLENGTH) * 
                              sizeof (unsigned long)));
    assert (sgEvent_indx !=  NULL);
    elength = ncount + EADDLENGTH;
  }

  ec = 0;
  sgp = sgEvent;
  for (i = 0; i < 2 * ng - 1; i++)
  {
    nowpop = gtree[i].pop;
    assert (nowpop < numtreepops);
    if (i < ng)
    {
      t = 0;
    }
    else
    {
      // every internal branch begins in a population by a coalescent event, we know the time and the pop
      t = gtree[gtree[i].up[0]].time;
      assert (i >= ng);
      sgp->time = t;
      sgp->cmt = 0;
      sgp->periodi = findperiod (ci, t);
      sgp->pop = nowpop;
      ec += 1;
      sgp++;
    }
    j = 0;
    while (gtree[i].mig[j].mt > -0.5)
    {
      //assert(gtree[i].mig[j].mt >= t); 
      t = gtree[i].mig[j].mt;
      sgp->time = t;
      sgp->periodi = findperiod (ci, t);

      /* could change population by passing into the next period */
      while (C[ci]->poptree[nowpop].e <= sgp->periodi
             && C[ci]->poptree[nowpop].e != -1)
      {
        nowpop = C[ci]->poptree[nowpop].down;
      }
      assert (ISELEMENT (nowpop, C[0]->periodset[sgp->periodi]));
      sgp->pop = nowpop;
      sgp->topop = gtree[i].mig[j].mp;
      assert (nowpop != sgp->topop);
      nowpop = sgp->topop;
      assert (sgp->topop >= 0 && sgp->topop < numtreepops);
      sgp->cmt = 1;        // a migration event
      ec += 1;
      sgp++;
      j++;
    }
  }

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0)
  {
    k = findperiod (ci, G->roottime);
    for (i = 0; i < k; i++)
    {
      sgp->time = C[ci]->tvals[i];
      sgp->periodi = i;
      sgp->cmt = -1;
      ec += 1;
      sgp++;
    }
  }
  assert (ncount == ec); /* FIXME: island model is checked here? */
  indexx ((unsigned long) ec, sgEvent - 1, sgEvent_indx - 1);  // sorts an index of locations in ec by time;-1 because NR routines use arrays that begin at 1

  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    for (i = 0; i < npops; i++)
      n[i] = 0;
    for (i = 0; i < ng; i++)
    {
      n[gtree[i].pop]++;
    }
  }
  else
  {
    for (i = 0; i < npops; i++)
      n[i] = L[li].samppop[i];  // use regular num indexing  (not sets)
  }

  for (i = npops; i < 2 * npops - 1; i++)
    n[i] = 0;
  nsum = L[li].numgenes;
  lasttime = 0;
  G->length = G->tlength = 0;   // does not really belong here, but this is an easy place to measure tlength, the total length of genealogy
  h2term = 1 / (2 * L[li].hval);
  lastsplitt = C[ci]->tvals[lastperiodnumber - 1];
  k = 0;                        //period
  for (j = 0; j < ec; j++)
  {
    sgp = sgEvent + (sgEvent_indx[j] - 1); // -1 because NR sort routines use arrays that begin at 1
    timeinterval = sgp->time - lasttime;
    assert (timeinterval >= 0);
    timeadd = nsum * timeinterval;
    G->length += timeadd;
    if (sgp->time < lastsplitt)
    {
      G->tlength += timeadd;
    }
    else
    {
      if (lasttime < lastsplitt)
        G->tlength += nsum * (lastsplitt - lasttime);
    }
    lasttime = sgp->time;
    for (ii = 0; ii < npops - k; ii++)
    {
      if (npops > 1)
      {
        ip = C[ci]->plist[k][ii];
      }
      else
      {
        ip = 0;
      }
      gweight->fc[k][ii] += nnminus1[n[ip]] * timeinterval * h2term;
      assert(gweight->fc[k][ii] >= 0);
      assert (gweight->fc[k][ii] < DBL_MAX);
      if (!modeloptions[NOMIGRATION] && k < lastperiodnumber)
      {
        fmtemp = n[ip] * timeinterval;
        for (jj = 0; jj < npops - k; jj++)
          if (jj != ii)
          {
            gweight->fm[k][ii][jj] += fmtemp;
          }
      }
    }
    switch (sgp->cmt)
    {
    case 0:
      if (npops > 1)
      {
        assert (ISELEMENT (sgp->pop, C[ci]->periodset[findperiod (ci, lasttime)]));
        ip = sgp->pop;
        assert (n[ip] > 1);
        ii = 0;
        // this while() may be slow, could have some kind of lookup 
        while (C[ci]->plist[k][ii] != ip)
        {
          ii++;
        }
      }
      else
      {
        ip = 0;
        ii = 0;
      }
      gweight->cc[k][ii]++;
      n[ip]--;
      nsum--;
      break;
    case 1:
      assert (npops > 1);
      ip = sgp->pop;
      jp = sgp->topop;
      ii = jj = 0;
      // this while() may be slow, could have some kind of lookup 
      while (C[ci]->plist[k][ii] != ip)
        ii++;
      assert (ii < npops - k);
      // this while() may be slow, could have some kind of lookup 
      while (C[ci]->plist[k][jj] != jp)
        jj++;
      assert (jj < npops - k);
      gweight->mc[k][ii][jj]++;
      assert (n[ip] > 0);
      n[ip]--;
      n[jp]++;
      break;
    case -1:
      assert (npops > 1);
      assert (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0);
      assert (sgp->time >= C[ci]->tvals[k]);
      assert (nsum > 1);
      assert (k < numsplittimes);
      k++;
      assert(k<=npops);
      n[C[ci]->addpop[k]] =
        n[C[ci]->droppops[k][0]] + n[C[ci]->droppops[k][1]];
      n[C[ci]->droppops[k][0]] = n[C[ci]->droppops[k][1]] = 0;
      break;
    }
  }
  hlog = log (L[li].hval);
  if (hlog != 0.0)  //added 9/13/2010
  {
    for (k = 0; k < numsplittimes + 1; k++)
      for (i = 0; i < npops - k; i++)
        gweight->hcc[k][i] += hlog * gweight->cc[k][i];
  }
  
  assert (nsum == 1);
}                               /*treeweight */

#undef ESTARTLENGTH
#undef EADDLENGTH



/* integrate_tree_prob - determine integrated probability of the genealogy  */
/* this function does not know about loci,  all the information needed is contained in 
  *gweight,  regardless of whether it is for one genealogy (called from updategenealogy) or all genealogies (called
  from changet ) 
  regardless of whether it is all genealogies or a single genealogy the caculations results are put into C[ci]->allpcalc*/

void
integrate_tree_prob (int ci,
                     struct genealogy_weights *gweight,
                     struct genealogy_weights *holdgweight,
                     struct probcalc *pcalc,
                     struct probcalc *holdpcalc, double *holdt)
{
  double psum;
  int i, j;
  int c, holdc;
  double hc, f, holdf;
  int checkm;
  int tempc = 0, holdtempc = 0;
  /* check to see if there are any migration events where there should not be under the model */
  if (modeloptions[NOMIGRATION] == 0 && nomigrationchecklist.n > 0)
  {
    i = 0;
    checkm = 1;
    while (checkm && i < nomigrationchecklist.n)
    {
      checkm = gweight->mc[nomigrationchecklist.p[i]]
                          [nomigrationchecklist.r[i]]
                          [nomigrationchecklist.c[i]] == 0;
      i++;
    }
  }
  else
  {
    checkm = 1;
  }

  if (!checkm)
  {
    psum = -MYDBL_MAX;          // very large neg value to ensure rejection of update if migration events occur where they should not under the model
  // not entirely sure if this is the best way to reject updates 
  }
  else
  {
    psum = 0;
    for (i = 0; i < numpopsizeparams; i++)
    {
      c = holdc = 0;
      f = holdf = hc = 0.0;
      for (j = 0; j < itheta[i].wp.n; j++)
      {
        c += gweight->cc[itheta[i].wp.p[j]][itheta[i].wp.r[j]];
        holdc += holdgweight->cc[itheta[i].wp.p[j]][itheta[i].wp.r[j]];
        f += gweight->fc[itheta[i].wp.p[j]][itheta[i].wp.r[j]];
        holdf += holdgweight->fc[itheta[i].wp.p[j]][itheta[i].wp.r[j]];
        hc += gweight->hcc[itheta[i].wp.p[j]][itheta[i].wp.r[j]];
      }
      tempc += c;
      holdtempc += c;
      if (c == holdc && f == holdf)
      {
        pcalc->qintegrate[i] = holdpcalc->qintegrate[i];
      }
      else
      {
        pcalc->qintegrate[i] =
          integrate_coalescent_term (c, f, hc, itheta[i].pr.max,
                                     itheta[i].pr.min);
      }

      psum += pcalc->qintegrate[i];
      //     assert(psum > -1e200 && psum < 1e200);
    }
assert(tempc == total_numgenes-nloci);
assert(holdtempc==tempc);

    if (!modeloptions[NOMIGRATION])
    {
      for (i = 0; i < nummigrateparams; i++)
      {
        c = holdc = 0;
        f = holdf = 0.0;
        for (j = 0; j < imig[i].wp.n; j++)
        {
          c += gweight->mc[imig[i].wp.p[j]][imig[i].wp.r[j]][imig[i].wp.c[j]];
          holdc +=
            holdgweight->mc[imig[i].wp.p[j]][imig[i].wp.r[j]][imig[i].
                                                              wp.c[j]];
          f += gweight->fm[imig[i].wp.p[j]][imig[i].wp.r[j]][imig[i].wp.c[j]];
          holdf +=
            holdgweight->fm[imig[i].wp.p[j]][imig[i].wp.r[j]][imig[i].
                                                              wp.c[j]];
        }
        if (c == holdc && f == holdf)
        {
          pcalc->mintegrate[i] = holdpcalc->mintegrate[i];
        }
        else
        {

          if (modeloptions[EXPOMIGRATIONPRIOR])
            pcalc->mintegrate[i] =
              integrate_migration_term_expo_prior (c, f, imig[i].pr.mean);
          else
            pcalc->mintegrate[i] =
              integrate_migration_term (c, f, imig[i].pr.max, imig[i].pr.min);
        }
        psum += pcalc->mintegrate[i];
        //     assert(psum > -1e200 && psum < 1e200);
      }
    }
    //assert(psum > -1e200 && psum < 1e200);
  }
  //assert(psum > -1e200 && psum < 1e200);
  pcalc->probg = psum;
}                               /* integrate_tree_prob */

// same as integrate_tree_prob but does not have checks to see if values match current ones */
void
initialize_integrate_tree_prob (int ci,
                                struct genealogy_weights *gweight,
                                struct probcalc *pcalc)
{
  double psum;
  int i, j;
  int c;
  double hc, f;
  int checkm;
  /* check to see if there are any migration events where there should not be under the model */
  if (modeloptions[NOMIGRATION] == 0 && nomigrationchecklist.n > 0)
  {
    i = 0;
    checkm = 1;
    while (checkm && i < nomigrationchecklist.n)
    {
      checkm = gweight->mc[nomigrationchecklist.p[i]]
                          [nomigrationchecklist.r[i]]
                          [nomigrationchecklist.c[i]] == 0;
      i++;
    }
  }
  else
  {
    checkm = 1;
  }

  if (!checkm)
  {
    psum = -MYDBL_MAX;          // very large neg value to ensure rejection of update if migration events occur where they should not under the model
  }
  else
  {
    psum = 0;
    for (i = 0; i < numpopsizeparams; i++)
    {
      c = 0;
      f = hc = 0.0;
      for (j = 0; j < itheta[i].wp.n; j++)
      {
        c += gweight->cc[itheta[i].wp.p[j]][itheta[i].wp.r[j]];
        f += gweight->fc[itheta[i].wp.p[j]][itheta[i].wp.r[j]];
        hc += gweight->hcc[itheta[i].wp.p[j]][itheta[i].wp.r[j]];
      }
      pcalc->qintegrate[i] =
        integrate_coalescent_term (c, f, hc, itheta[i].pr.max,
                                   itheta[i].pr.min);
      psum += pcalc->qintegrate[i];
      assert (psum > -1e200 && psum < 1e200);
    }
    if (!modeloptions[NOMIGRATION])
    {
      for (i = 0; i < nummigrateparams; i++)
      {
        c = 0;
        f = 0.0;
        for (j = 0; j < imig[i].wp.n; j++)
        {
          c += gweight->mc[imig[i].wp.p[j]][imig[i].wp.r[j]][imig[i].wp.c[j]];
          f += gweight->fm[imig[i].wp.p[j]][imig[i].wp.r[j]][imig[i].wp.c[j]];
        }
        if (modeloptions[EXPOMIGRATIONPRIOR])
        {
          pcalc->mintegrate[i] =
            integrate_migration_term_expo_prior (c, f, imig[i].pr.mean);
        }
        else
        {
          pcalc->mintegrate[i] =
            integrate_migration_term (c, f, imig[i].pr.max, imig[i].pr.min);
        }
        psum += pcalc->mintegrate[i];
      }
    }
    assert (psum > -1e200 && psum < 1e200);
  }
  pcalc->probg = psum;
}                               /* initialize_integrate_tree_prob */



/* for HKY model */
void
copyfraclike (int ci, int li)
{
  int i, j, k;
  struct edge *gtree = C[ci]->G[li].gtree;
  int ng = L[li].numgenes;
  for (i = ng; i < 2 * ng - 1; i++)
  {
    if (gtree[i].hkyi.newfrac[0][0] != -1)
    {
      for (j = 0; j < L[li].numsites; j++)
        for (k = 0; k < 4; k++)
          gtree[i].hkyi.frac[j][k] = gtree[i].hkyi.newfrac[j][k];
    }
  }
}

/* for HKY model */
void
storescalefactors (int ci, int li)
{
  int i, k;
  struct edge *gtree = C[ci]->G[li].gtree;
  for (i = L[li].numgenes; i < 2 * L[li].numgenes - 1; i++)
    for (k = 0; k < L[li].numsites; k++)
      gtree[i].hkyi.oldscalefactor[k] = gtree[i].hkyi.scalefactor[k];
}


/* for HKY model */
void
restorescalefactors (int ci, int li)
{
  int i, k;
  struct edge *gtree = C[ci]->G[li].gtree;
  for (i = L[li].numgenes; i < 2 * L[li].numgenes - 1; i++)
  {
    for (k = 0; k < L[li].numsites; k++)
    {
      gtree[i].hkyi.scalefactor[k] = gtree[i].hkyi.oldscalefactor[k];
    }
  }
}


/*  finishSWupdateA() returns the likelihood associated with the allele states on the parts of the gtree that have changed
   it returns the difference in likelihoods associated with these parts.  thus the overall likelihood of data can just be
   updated by the value returned from this function. 

	when  finishSWupdate() is entered the gtree is not complete, as the allele state of the new downedge is not known
    must pick an allele state for the new downedge:
    the new downedge has the same number as the old downedge, so it is just called downedge hereafter

    at the start the current allele value for downedge is oldA - this however comes from a different part 
       of the gtree from the old downedge that was erased
   
    consider - the point at the base of edge is connected to 3 other nodes (2 up and 1 down) unless
    edge is the root in which case there are only 2 ups
    
    start w/ a old value of A at the top of downedge = oldA

    determine the mean (weighted by the length of the branch to the node) of the difference between the A values
    of the connecting nodes and oldA 
    Treat this mean value as the parameter of a geometric distribution.
    pick an rv that is in the necessary direction (+ or -) and add this to oldA  this becomes newA. 
    
    the allele value for downedge becomes newA. the new genealogy is now complete.
    
    now consider the reverse update from new to old A

    again, calculate the weighted mean difference (for the old genealogy) between adjacent allele states and newA. 
    pick a geometric rv 

    Consider the update as this - When the gtree jumped from one genealogy to another, this caused oldA to be set to newA

        
    Aterm = log(Pr(oldstate in node | T* -> T)/Pr(new state in node | T -> T*)) is the ratio of two geometric probabilities. 
	Aterm is part of the Hastings term, the proposal ratio.
    This is passed to the calling function along with the new likelihood. */
double
finishSWupdateA (int ci, int li, int ai, int edge,
                 int downedge, int sisedge, int newsisedge,
                 double u, double *Aterm)
{
  int i, j, newA, oldA;
  int d, dA, e[3];
  double geotermnew, geotermold;
  double upt, t[3];
  int wsumdiff;
  double oldlikeadj = 0, likeadj = 0;
  struct edge *gtree = C[ci]->G[li].gtree;
  oldA = gtree[downedge].A[ai];
  if (newsisedge != sisedge)
    holdsisdlikeA[ai] = gtree[newsisedge].dlikeA[ai];
  else
    holdsisdlikeA[ai] = 0;


  /* NOTE: Function storeoldedge must be called to set copyedge. Sang Chul has 
   * used finishSWupdateA without calling function storeoldedge. */
  oldlikeadj =
    copyedge[0].dlikeA[ai] + copyedge[1].dlikeA[ai] +
    copyedge[2].dlikeA[ai] + holdsisdlikeA[ai];
  e[0] = edge;
  e[1] = newsisedge;
  e[2] = gtree[e[0]].down;
  for (i = 0; i < 3; i++)
  {
    if (e[i] >= L[li].numgenes)
    {
      upt = gtree[gtree[e[i]].up[0]].time;
      assert (upt == gtree[gtree[e[i]].up[1]].time);
    }
    else
    {
      upt = 0;
    }

    if (gtree[e[i]].down != -1)
    {
      t[i] = gtree[e[i]].time - upt;
      assert (t[i] >= 0.0);
    }
  }
  assert (gtree[e[0]].time == gtree[e[1]].time);
  wsumdiff = 0;
  for (i = 0, j = 0; i < 3; i++)
  {
    if (gtree[e[i]].down != -1)
    {
      if (i < 2)
        d = gtree[e[i]].A[ai] - oldA;
      else
        d = gtree[gtree[e[i]].down].A[ai] - oldA;
      wsumdiff += abs (d);
      j++;
    }
  }
  geotermnew = j / ((double) (wsumdiff + j));
  assert (geotermnew > 0);
  if (geotermnew > 0.95)
    geotermnew = 0.95;
  dA = geometric (geotermnew) - 1;
  if (bitran () /*uniform() < 0.5 */ )
    dA = -dA;
  if (dA >= 0)
    newA = IMIN (L[li].maxA[ai], oldA + dA);
  else
    newA = IMAX (L[li].minA[ai], oldA + dA);
  gtree[downedge].A[ai] = newA;
  assert (newA >= L[li].minA[ai]);
  dA = newA - oldA;             /* difference between new value, and what it would be based simply on weighted mean */
  if (gtree[sisedge].down != -1 && sisedge != newsisedge)
  {
    if (sisedge >= L[li].numgenes)
    {
      upt = gtree[gtree[sisedge].up[0]].time;
    }
    else
    {
      upt = 0;
    }

    /* set dlikeA for the old sisedge (which is now continuous with the old downedge */
    d = gtree[sisedge].A[ai] - gtree[gtree[sisedge].down].A[ai];
    assert (d <= L[li].maxA[ai]);
    gtree[sisedge].dlikeA[ai] = -(gtree[sisedge].time - upt) * u + log (bessi (d, (gtree[sisedge].time - upt) * u));
  }
  else
  {
    gtree[sisedge].dlikeA[ai] = 0;
  }

  likeadj = gtree[sisedge].dlikeA[ai];
  for (i = 0, j = 0; i < 3; i++)
  {
    if (gtree[e[i]].down != -1)
    {
      if (i < 2)
        d = gtree[e[i]].A[ai] - newA;
      else
        d = gtree[gtree[e[i]].down].A[ai] - newA;
      gtree[e[i]].dlikeA[ai] = -(t[i] * u) + log (bessi (d, t[i] * u));
      likeadj += gtree[e[i]].dlikeA[ai];
    }
  }
  wsumdiff = 0;
  j = 0;
  for (i = 0; i < 2; i++)
  {
    d = copyedge[i].A[ai] - newA;
    wsumdiff += abs (d);
    j++;
  }
  if (copyedge[2].down != -1)
  {
    d = holddownA[ai] - newA;
    wsumdiff += abs (d);
    j++;
  }
  geotermold = j / ((double) (wsumdiff + j));
  if (geotermold > 0.95)
    geotermold = 0.95;

  /* FIXME: Tue Jan 27 17:35:50 EST 2009
   * SANGCHUL thought that we could return 0 without computing anything. Now
   * that I think that we may have to compute parts of the above to propose a
   * new state of alleles at internal nodes. We may place this before computing
   * likelihood because we do not need compute likelihood when we run without
   * data. In the same reason, we may want to rearrange the code of function
   * likelihoodSW. */
  if (calcoptions[DONTCALCLIKELIHOODMUTATION])
  {
    *Aterm = 0;
    return 0;
  }
  else
  {
    *Aterm =
      (abs (dA) * log (1 - geotermold) + log (geotermold)) -
      (abs (dA) * log (1 - geotermnew) + log (geotermnew));
    return (likeadj - oldlikeadj);
  }
}                               /* finishSWupdateA */


#if 0
/*  CR 110825.1
 *  updateA() not used and no longer compiled into code
 */

double
updateA (int ci, int li, int ai, double u, int *count)
/* implements a geometric distribution based on the weighted mean difference in allele sizes surrounding the value in the 
node to be updated. branch lengths are used for the weighting */
/* would it be good to randomize the order in which nodes are updated ?  */
{
  int i, upA[2], downA, newA, oldA;
  int d;
  int up[2], upup[2], down;
  double pdg, like, dlikeup[2], dlikedown;
  double tup[2], tdown, weightsum, metropolishastingsterm;
  double wsumdiff, geotermnew, geotermold, Aterm;
  struct edge *gtree = C[ci]->G[li].gtree;
  int ng = L[li].numgenes;
  *count = 0;
  if (calcoptions[DONTCALCLIKELIHOODMUTATION])
  {
    *count = (int) (updateAfrac * (ng - 1));    // don't do anything as it can have no effect 
    return 0;
  }
  for (i = ng; i < 2 * ng - 1; i++)
  {
    /* don't bother with all nodes  just do uprop of them */
    if (uniform () < updateAfrac)
    {
      L[li].A_rec[ai].upinf->tries++;
      up[0] = gtree[i].up[0];
      up[1] = gtree[i].up[1];
      upup[0] = gtree[up[0]].up[0];
      upup[1] = gtree[up[1]].up[0];
      down = gtree[i].down;
      if (upup[0] == -1)
        tup[0] = gtree[up[0]].time;
      else
        tup[0] = gtree[up[0]].time - gtree[upup[0]].time;
      if (upup[1] == -1)
        tup[1] = gtree[up[1]].time;
      else
        tup[1] = gtree[up[1]].time - gtree[upup[1]].time;
      oldA = gtree[i].A[ai];
      upA[0] = gtree[gtree[i].up[0]].A[ai];
      upA[1] = gtree[gtree[i].up[1]].A[ai];
      pdg = gtree[gtree[i].up[0]].dlikeA[ai] + gtree[gtree[i].up[1]].dlikeA[ai];

      wsumdiff = abs (oldA - upA[0]) / tup[0] + abs (oldA - upA[1]) / tup[1];
      weightsum = 1 / tup[0] + 1 / tup[1];
      if (down != -1)
      {
        tdown = gtree[i].time - gtree[up[0]].time;
        downA = gtree[gtree[i].down].A[ai];
        wsumdiff += abs (oldA - downA) / tdown;
        weightsum += 1 / tdown;
        pdg += gtree[i].dlikeA[ai];
      }
      else
      {
        downA = -1;
      }
      geotermnew = weightsum / (wsumdiff + weightsum);
      assert (geotermnew > 0);
      if (geotermnew > 0.95)
        geotermnew = 0.95;
      d = geometric (geotermnew) - 1;
      assert (d >= 0);
      if (bitran () /*uniform() < 0.5 */ )
        newA = IMAX (L[li].minA[ai], oldA - d);
      else
        newA = IMIN (L[li].maxA[ai], oldA + d);

      if (newA != oldA)
      {
        wsumdiff = abs (newA - upA[0]) / tup[0] + abs (newA - upA[1]) / tup[1];
        if (down != -1)
        {
          wsumdiff += abs (newA - downA) / tdown;
        }
        geotermold = weightsum / (wsumdiff + weightsum);
        assert (geotermold > 0);
        if (geotermold > 0.95)
          geotermold = 0.95;
        d = newA - upA[0];
        like = dlikeup[0] = -(tup[0] * u) + log (bessi (d, tup[0] * u));
        d = newA - upA[1];
        like += dlikeup[1] = -(tup[1] * u) + log (bessi (d, tup[1] * u));
        if (down != -1)
        {
          d = newA - downA;
          like += dlikedown = -(tdown * u) + log (bessi (d, tdown * u));
        }
        Aterm = abs (newA - oldA) * log ((1 - geotermold) / (1 - geotermnew)) + log (geotermold / geotermnew);

        /* 5/19/2011 JH adding thermodynamic integration */
        /* this use of beta is not affected by whether or not the probability of the genealogy is included, 
          since it is not present in this MH term */
        metropolishastingsterm = exp (beta[ci] * (like - pdg) + Aterm);
        if (metropolishastingsterm >= 1.0
            || metropolishastingsterm > uniform ())
        {
          gtree[i].A[ai] = newA;
          if (down != -1)
            gtree[i].dlikeA[ai] = dlikedown;
          gtree[up[0]].dlikeA[ai] = dlikeup[0];
          gtree[up[1]].dlikeA[ai] = dlikeup[1];
          //*count = *count + 1;
          *count = *count + 1;
          L[li].A_rec[ai].upinf->accp++;
        }
      }
    }
  }
  like = 0;
  for (i = 0; i < 2 * ng - 1; i++)
  {
    if (gtree[i].down != -1)
    {
      like += gtree[i].dlikeA[ai];
    }
  }
  return like;
}                               /*updateA */

#endif  // #if 0

/* getnewt() these window widths just seem to work,  no method to them */
/*whichupdate:
0 : NW
1:  RY1
*/
#define MAXWIN 10
// 9/10/10  trimmed down getnewt()
double
getnewt (int timeperiod, double t_u_prior, double t_d_prior, double oldt, int whichupdate)
{
  double twin, newt;
  double U;


  twin =  (t_d_prior - t_u_prior) / (log (static_cast<double>(nloci) + 1) * (npops - timeperiod));
  U = uniform ();
  newt = (oldt - twin / 2) + U * twin;
  if (newt >= t_d_prior)
    newt = 2.0 * t_d_prior - newt;
  else
  {
    if (newt <= t_u_prior)
      newt = 2.0 * t_u_prior - newt;
  }
  assert (newt < t_d_prior && newt > t_u_prior);
  return newt;
}                               //getnewt()


