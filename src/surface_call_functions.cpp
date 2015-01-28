/*IMa2p 2009-2015 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, and Arun Sethuraman */

/* funccall.c   functions associated w/ the surface, and that make calls optimization routines */
#undef GLOBVARS
#include "imamp.hpp"

#define LOG_10_2  0.30102999566398119521
#define OCUTOFF  10
#define LOWVAL 1e-200
#define LOWLOG  -1e200

/*********** LOCAL STUFF **********/


/* function prototypes */

static double marginp (int param, int firsttree, int lasttree, double x, int dummy);       //calculate marginal probability
static double marginbis (double (*func) (double, double, int, int), double x1, double x2, double yadjust, int pi);       // find a value associated w/ a certain marginal probability


/****** LOCAL FUNCTIONS *********/


#define OFFSCALEVAL 1
double
marginp (int param, int firsttree, int lasttree, double x, int dummy)
{
  int ei, p;
  double hval, prob, temp, max, min, sumtemp = 0;
  double meani;

  dummy = 0; /* a dummy value is set to remove compiler complaints. */
  if (param < numpopsizeparams)
  {
    max = itheta[param].pr.max;
    min = itheta[param].pr.min;
  }
  else//param < numpopsizeparams + nummigrateparams
  {
    max = imig[param - numpopsizeparams].pr.max;
    min = imig[param - numpopsizeparams].pr.min;
    if (modeloptions[EXPOMIGRATIONPRIOR])
      meani = 1.0 / imig[param - numpopsizeparams].pr.mean;
  }
  if (x < min || x > max)
    return OFFSCALEVAL;
    //return 0;

  for (ei = firsttree; ei < lasttree; ei++)
  {
    if (param < numpopsizeparams)
    {
      p = param;
      hval = gsampinf[ei][gsamp_hccp + p];
      temp =
        -gsampinf[ei][gsamp_qip + p] + gsampinf[ei][gsamp_ccp +
                                                    p] * (LOG2 -
                                                          log (x)) -
        hval - 2 * gsampinf[ei][gsamp_fcp + p] / x;
      sumtemp += exp (temp);
    }
    else //param < numpopsizeparams + nummigrateparams
    {
      assert (param < numpopsizeparams + nummigrateparams);
      p = param - numpopsizeparams;
      if (modeloptions[EXPOMIGRATIONPRIOR])
        temp = -gsampinf[ei][gsamp_mip + p] + log (meani) - x * meani +
          INTEGERROUND (gsampinf[ei][gsamp_mcp + p]) *
          log (x) - gsampinf[ei][gsamp_fmp + p] * x;
      else
        temp = -gsampinf[ei][gsamp_mip + p] +
          INTEGERROUND (gsampinf[ei][gsamp_mcp + p]) *
          log (x) - gsampinf[ei][gsamp_fmp + p] * x;
      sumtemp += exp (temp);
    }
  }
  sumtemp /= (lasttree - firsttree + (firsttree == 0)); 
  prob = sumtemp;
  return -prob;                 /* negative because a minimization routine is used */
}                               /* marginp */
#undef OFFSCALEVAL

#define JMAX 40
#define BISTOL  1e-4
double
marginbis (double (*func) (double, double, int, int), double x1, double x2,
           double yadjust, int pi)
{
  int j;
  double dx, f, fmid, xmid, rtb;
  f = (*func) (x1, yadjust, pi, 1);
  fmid = (*func) (x2, yadjust, pi, 1);
  if (f * fmid >= 0.0)
    return DBL_MIN;             /* not found does not appear to be a point corresponding to 95% limit */
  rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
  for (j = 1; j <= JMAX; j++)
  {
    fmid = (*func) (xmid = rtb + (dx *= 0.5), yadjust, pi, 1);
    if (fmid <= 0.0)
      rtb = xmid;
    if (fabs (dx) < BISTOL || fmid == 0.0)
      return rtb;
  }
  return DBL_MAX;               /* Too many bisections in marginbis - cannot find root */
}


/********** GLOBAL FUNCTIONS ***********/

/* margincalc does the same thing as marginp, but gets called by different functions
unlike marginp, margincalc automatically uses all trees 
can be used for direct calculation
also can be used for root finding. 
also by using a nonzero value of jadust can be used to find the value of x that is
associated with a particular value of y (i.e. the likelihood) 
if logi==1  return the logarithm 
 */

double
margincalc (double x, double yadjust, int pi, int logi)
{
  int ei, p;
  double hval, sum, temp,meani;
  sum = 0;

  /* calculate meani for migration params only */
  if (modeloptions[EXPOMIGRATIONPRIOR] && 
        (pi >= numpopsizeparams) && /* cr 110907.1 pi must be at least equal */ 
                                    /* to numpopsizeparams */
        (pi < numpopsizeparams + nummigrateparams)) 
  {
    meani = 1.0 / imig[pi - numpopsizeparams].pr.mean;
  }


  for (ei = 0; ei < genealogiessaved; ei++)
  {
    if (pi < numpopsizeparams)
    { /*  for popsize params only  */
      p = pi;
      hval = gsampinf[ei][gsamp_hccp + p];
      temp =
        -gsampinf[ei][gsamp_qip + p] +
        INTEGERROUND (gsampinf[ei][gsamp_ccp + p]) * (LOG2 - log (x)) -
        hval - 2 * gsampinf[ei][gsamp_fcp + p] / x;
      sum += exp (temp);
    }
    else //pi < numpopsizeparams + nummigrateparams
    {  /*  for migration params only  */
      p = pi - numpopsizeparams;

      if (modeloptions[EXPOMIGRATIONPRIOR])
        sum += exp (-gsampinf[ei][gsamp_mip + p] + log (meani) - x * meani +
                    INTEGERROUND (gsampinf[ei][gsamp_mcp + p]) *
                    log (x) - gsampinf[ei][gsamp_fmp + p] * x);
      else
        sum += exp (-gsampinf[ei][gsamp_mip + p] +
                    INTEGERROUND (gsampinf[ei][gsamp_mcp + p]) *
                    log (x) - gsampinf[ei][gsamp_fmp + p] * x);
    }
  }
  sum /= genealogiessaved;

  if (logi)
  {
    if (sum <= 0)
      sum = LOWLOG;
    else
      sum = log (sum);
  }
  sum -= yadjust;
  return sum;
}                               /* marginalcalc */


#define SEARCHSTARTFRAC  4      // fraction of position in parameter range to start at
void
marginalopt (int firsttree, int lasttree, double *mlval, double *peakloc)
{
  int i;
  double ftol, ax, bx, cx, fa, fb, fc, xmax, ml;
  double axt, bxt, cxt, prior;
  double max0, min0, max1, min1;
  double (*func) (int, int, int, double, int);
  ftol = 1e-7;
  func = marginp;
  for (i = 0;  i < numpopsizeparams + nummigrateparams; i++)
  {
    if (i < numpopsizeparams)
    {
      prior = itheta[i].pr.max;
    }
    else
    {
      prior = imig[i - numpopsizeparams].pr.max;
    }

    ax = prior;
    bx = prior / 2;
    mnbrakmod (i, firsttree, lasttree, &ax, &bx, &cx, &fa, &fb, &fc, func,0); // why is this passing i??  shouldn't it be the number of parameters?
    axt = ax;
    bxt = bx;
    cxt = cx;
    bx = prior / 2;
    ax = MINPARAMVAL;
    mnbrakmod (i, firsttree, lasttree, &ax, &bx, &cx, &fa, &fb, &fc, func,0); // why is this passing i??  shouldn't it be the number of parameters?
    if (axt < bxt && axt < cxt)
      min0 = axt;
    if (bxt < axt && bxt < cxt)
      min0 = bxt;
    if (cxt < axt && cxt < bxt)
      min0 = cxt;
    if (axt > bxt && axt > cxt)
    {
      if (axt > prior)
        axt = prior;
      max0 = axt;
    }
    if (bxt > axt && bxt > cxt)
    {
      if (bxt > prior)
        bxt = prior;
      max0 = bxt;
    }
    if (cxt > axt && cxt > bxt)
    {
      if (cxt > prior)
        cxt = prior;
      max0 = cxt;
    }
    if (ax < bx && ax < cx)
      min1 = ax;
    if (bx < ax && bx < cx)
      min1 = bx;
    if (cx < ax && cx < bx)
      min1 = cx;
    if (ax > bx && ax > cx)
      max1 = ax;
    {
      if (ax > prior)
        ax = prior;
      max1 = ax;
    }
    if (bx > ax && bx > cx)
    {
      if (bx > prior)
        bx = prior;
      max1 = bx;
    }
    if (cx > ax && cx > bx)
    {
      if (cx > prior)
        cx = prior;
      max1 = cx;
    }
    if (max0 <= min1 || max1 <= min0)
    {
      peakloc[i] = -1;
    }
    else
    {
      ml = -goldenmod (i, firsttree, lasttree, ax, bx, cx, ftol, &xmax, func,0); // why is this passing i??  shouldn't it be the number of parameters?
      mlval[i] = ml;
      peakloc[i] = xmax;
    }
  }
}                               /* marginalopt */


#define  DOWN95  1.92
double
margin95 (double mlval[], double peakloc[], int pi, int UL)
{
  double x1, x2, x, yadjust;
  if (UL == 0)                  // lower
  {
    x1 = MINPARAMVAL;
    x2 = peakloc[pi];
  }
  else                          // upper
  {
    x1 = peakloc[pi];
    if (pi < numpopsizeparams)
    {
      x2 = itheta[pi].pr.max;
    }
    else //pi < numpopsizeparams + nummigrateparams
    {
      x2 = imig[pi - numpopsizeparams].pr.max;
    }
  }
  yadjust = log (mlval[pi]) - DOWN95;
  x = marginbis (margincalc, x1, x2, yadjust, pi);
  return x;
}                               /* margin95 */

#define NUMTREEINT 2            // consider NUMTREEINT batches of trees for finding peaks,  one way to check for convergence

/* findmarginpeaks() 
	1) find marginal peaks for the main model  - save points in peakloc
	2) find 95% confidence limits on marginal peak locations 
		these calls ultimately go to the function marginp() which determines the marginal function value
    3)  does an LLR test on migration parameters using as a test distribution a distributino that is 50% 0 
    and 50% x^2_1df  This distribution has values
    2.70554  at p=0.05   The ratio of probabilities (as opposed to twice the log ratio) is 3.86813
    5.41189	  at p = 0.01  the ratio of prbabilities is 14.9685
    9.54954	 at p = 0.001  the ratio of probabilities is 118.483

*/

/* findmarginpeaks is a complex function that builds a table of parameter peak locations
estimated from the marginal posterior density
also does LLR tests */ 

void
findmarginpeaks (FILE * outfile, float *holdpeakloc)
{
  int i, j, k, ii, ilo, ihi, iihi, iilo, p;
  double temp, prior;
  double **mlval, **peakloc;
  double **popmigmlval, **popmigpeakloc;
  double *migtest, *popmigtest, *temptest;
  int *mpop, *mterm;
  char **popmigstr;
  int nmi, tempmpop, thetai, found, mi;
  int firsttree, lasttree;
  int printerrorfootnote = 0;
  int printsigfootnote = 0;
  int nummigprint;
                            // c++ will append null terminator during init
  char sig[4][4] = { "ns", "*", "**", "***" };
  char llrstring[20];
  double maxp, max0p;

  p = numpopsizeparams + nummigrateparams;
  mlval = orig2d_alloc2Ddouble (NUMTREEINT + 1, p);
  peakloc = orig2d_alloc2Ddouble (NUMTREEINT + 1, p);
  if (outputoptions[POPMIGPARAMHIST])
  {
    popmigmlval = orig2d_alloc2Ddouble (NUMTREEINT + 1, nummigrateparams);
    popmigpeakloc = orig2d_alloc2Ddouble (NUMTREEINT + 1, nummigrateparams);
    popmigtest = static_cast<double *> 
                (malloc (nummigrateparams * sizeof (double)));
    mpop = static_cast<int *> (malloc (nummigrateparams * sizeof (int)));
    mterm = static_cast<int *> (malloc (nummigrateparams * sizeof (int)));
    popmigstr = static_cast<char **> 
                (malloc (nummigrateparams * sizeof (char *)));
    for (i=0;i<nummigrateparams;i++)
      popmigstr[i] = static_cast<char *> (malloc(PARAMSTRLEN *sizeof(char)));
    nmi = 0;
    if (modeloptions[PARAMETERSBYPERIOD])
    {
      for (k = 0; k < lastperiodnumber; k++)
      {
        for (i = 0; i < npops - k; i++)
        {
          tempmpop = C[0]->plist[k][i];
          thetai = 0;
          found = 0;
          while (!found && thetai < numpopsizeparams)
          {
            found = (k == atoi (&itheta[thetai].str[1]) && tempmpop == atoi (&itheta[thetai].str[3]));
            if (!found)
              thetai++;
          }
          assert (thetai < numpopsizeparams);
          for (mi = 0; mi < nummigrateparams; mi++)
          {
            found = 0;
            if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])
            {
              found = (k == atoi (&imig[mi].str[1])
                       && (tempmpop == atoi (&imig[mi].str[3])
                           || tempmpop == atoi (&imig[mi].str[6])));
            }
            else
            {
              found = (k == atoi (&imig[mi].str[1])
                       && tempmpop == atoi (&imig[mi].str[3]));
            }
            if (found)
            {
              mpop[nmi] = tempmpop;
              mterm[nmi] = mi;
              sprintf (popmigstr[nmi], "%d,2N%d", k, tempmpop);
              strcat (popmigstr[nmi], "M");
              if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])
              {
                strcat (popmigstr[nmi], &imig[mi].str[3]);
              }
              else
              {
                strcat (popmigstr[nmi], &imig[mi].str[1]);
              }
              nmi++;
            }
          }
        }
      }
    }
    else
    {
      for (i = 0; i < numtreepops - 1; i++)
      {
        thetai = i;
        for (mi = 0; mi < nummigrateparams; mi++)
        {
          if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])
          {
            found = ((thetai == atoi (&imig[mi].str[1])) || (thetai == atoi (&imig[mi].str[4])));
          }
          else
          {
            found = thetai == atoi (&imig[mi].str[1]);
          }
          if (found)
          {
            mpop[nmi] = thetai;
            mterm[nmi] = mi;
            sprintf (popmigstr[nmi], "2N%d", thetai);
            strcat (popmigstr[nmi], "M");
            strcat (popmigstr[nmi], &imig[mi].str[1]);
            nmi++;
          }
        }
      }
    }
  }
  migtest = static_cast<double *> (malloc (nummigrateparams * sizeof (double)));
  FP "\nMarginal Peak Locations and Probabilities\n");
  FP "=========================================\n");
  FP "  peak locations are estimated using a peak finding algorithm, which may\n");
  FP "  fail if the curve has multiple peaks. All peaks can also be found, and\n"); 
  FP "  related analyses conducted, by plotting the histograms\n\n");
  printf ("Finding marginal peaks and probabilities \n");
  genealogiessaved = IMIN (MAXGENEALOGIESTOSAVE - 1, genealogiessaved);   // trap cases when too saving too many trees is attempted
  if (genealogiessaved <= 10)
  {
    FP " TOO FEW TREES SAVED - MARGINAL VALUES NOT FOUND \n");
    for (i = 0; i < p; i++)
      holdpeakloc[i] = -1;
  }
  else
  {
    for (firsttree = 0, lasttree = (int) genealogiessaved / NUMTREEINT, j = 0;
         j < NUMTREEINT; j++)
    {
      marginalopt (firsttree, lasttree, mlval[j], peakloc[j]);
      if (outputoptions[POPMIGPARAMHIST])
        marginalopt_popmig (firsttree, lasttree, popmigmlval[j], popmigpeakloc[j],mpop, mterm);
      firsttree = lasttree + 1;
      lasttree += (int) genealogiessaved / NUMTREEINT;
      if (lasttree > genealogiessaved)
        lasttree = genealogiessaved;
    }
    marginalopt (0, genealogiessaved, mlval[NUMTREEINT], peakloc[NUMTREEINT]);
    for (i = 0; i < nummigrateparams; i++)
    {
      maxp = -marginp (i + numpopsizeparams, 0, genealogiessaved, peakloc[NUMTREEINT][i + numpopsizeparams], 0);
      max0p = -marginp (i + numpopsizeparams, 0, genealogiessaved, MINPARAMVAL, 0);
      migtest[i] = 2 * log(maxp/max0p);
    }

    if (outputoptions[POPMIGPARAMHIST])
    {
      marginalopt_popmig (0, genealogiessaved, popmigmlval[NUMTREEINT], popmigpeakloc[NUMTREEINT], mpop, mterm);
      for (i = 0; i < nummigrateparams; i++)
      {
        maxp = -marginpopmig (mterm[i], 0, genealogiessaved, popmigpeakloc[NUMTREEINT][i],mpop[i]);
        max0p = -marginpopmig (mterm[i], 0, genealogiessaved, MINPARAMVAL,mpop[i]);
        popmigtest[i] = 2 * log(maxp/max0p);
      }
    }

    k = 0;
    for (; k <= numsplittimes; k++)
    {
      FP "\n");
      FP "Period %d\n", k);
      FP "--------\n");
      if (k < lastperiodnumber)       // set the upper bound on the ii loop,  no migration parameters if k==lastperiodnumber
      {
        iilo = 1;
        if (outputoptions[POPMIGPARAMHIST])
          iihi = 3;
        else
          iihi = 2;
      }
      else
      {
        iilo = 1;
        iihi = 1;
      }
      for (ii = iilo; ii <= iihi; ii++)        // small loop ii==0 for population size parameters, ii==1 for migration parameters
      {
        switch (ii)
        {
        case 1:                //population size
          FP "Population Size Parameters\n Param:");
          for (i = 0; i < numpopsizeparams; i++)
            if (itheta[i].b == k)
              FP "\t %s\tP", itheta[i].str);
          FP "\n");
          ilo = 0;
          ihi = numpopsizeparams;
          break;
        case 2:                // migration 
          nummigprint = 0;
          for (i = 0; i < nummigrateparams; i++)
            nummigprint += imig[i].b == k;
          if (nummigprint)
          {
            FP "Migration Rate Parameters\n Param:");
            for (i = 0; i < nummigrateparams; i++)
              if (imig[i].pr.max > MPRIORMIN)
                if (imig[i].b == k)
                  FP "\t %s\tP", imig[i].str);
            FP "\n");
          }
          ilo = numpopsizeparams;
          ihi = numpopsizeparams + nummigrateparams;
          temptest = migtest;
          break;
        case 3:                // 2NM 
          nummigprint = 0;
          for (i = 0; i < nummigrateparams; i++)
            nummigprint += imig[i].b == k;
          if (nummigprint)
          {
            FP "Population Migration (2NM) Terms\n Term:");
            for (i = 0; i < nummigrateparams; i++)
              if (imig[mterm[i]].pr.max > MPRIORMIN)
                if (imig[mterm[i]].b == k)
                  FP "\t %s\tP", popmigstr[i]);
            FP "\n");
          }
          ilo = numpopsizeparams;
          ihi = numpopsizeparams + nummigrateparams;
          temptest = popmigtest;
          break;
        }
        for (j = 0; j <= NUMTREEINT; j++)
          if (ii == 0 || (ii == 2 && nummigprint > 0) || (ii == 3 && nummigprint > 0) || ii == 1)
          {
            if (j < NUMTREEINT)
              FP " Set%d", j);
            else
              FP " All");
            for (i = ilo; i < ihi; i++)
              if ((ii == 0)
                  || (ii == 1 && itheta[i - ilo].b == k)
                  || (ii == 2 && imig[i - ilo].b == k)
                  || (ii==3 && imig[mterm[i - ilo]].b == k)
                  )

              {
                if (ii < 3)
                {
                  if (peakloc[j][i] >= 0)
                  {
                    FP "\t%7.3lf\t%7.3lf", peakloc[j][i], mlval[j][i]);
                  }
                  else
                  {
                    FP "\terror*\t");
                    printerrorfootnote = 1;
                  }
                }
                else
                {
                  if (popmigpeakloc[j][i-ilo] >= 0)
                  {
                    FP "\t%7.3lf\t%7.3lf", popmigpeakloc[j][i-ilo], popmigmlval[j][i-ilo]);
                  }
                  else
                  {
                    FP "\terror*\t");
                    printerrorfootnote = 1;
                  }
                }
              }
            FP "\n");
          }
        if (ii == 0 || ii == 1 || ((ii == 2|| ii==3)&& nummigprint > 0))
        {
          if (ii < 3)
          {
            FP " LR95%%Lo");
            for (i = ilo; i < ihi; i++)
            {
              if ((ii == 0)
                  || (ii == 1 && itheta[i - ilo].b == k)
                  || (ii == 2 && imig[i - ilo].b == k))

              {
                if (peakloc[NUMTREEINT][i] >= 0)
                {
                  temp =
                    margin95 (mlval[NUMTREEINT], peakloc[NUMTREEINT], i, 0);
                  if (temp <= 0)
                  {
                    FP "\t<min\t");
                  }
                  else
                  {
                    if (temp >= DBL_MAX || temp <= DBL_MIN)
                      FP "\tna\t");
                    else
                      FP "\t%7.3lf\t", temp);
                  }
                }
                else
                {
                  FP "\t\t");
                }
              }
            }
            FP "\n");
            FP " LR95%%Hi");
            for (i = ilo; i < ihi; i++)
            {
              if ((ii == 0)
                  || (ii == 1 && itheta[i - ilo].b == k)
                  || (ii == 2 && imig[i - ilo].b == k))

              {
                if (ii == 1)
                  prior = itheta[i - ilo].pr.max;
                else
                  prior = imig[i - ilo].pr.max;
                if (peakloc[NUMTREEINT][i] >= 0)
                {
                  temp = margin95 (mlval[NUMTREEINT], peakloc[NUMTREEINT], i, 1);
                  if (temp >= prior || temp <= DBL_MIN)
                  {
                    FP "\t>max\t");
                  }
                  else
                  {
                    if (temp <= DBL_MIN)
                      FP "\tna\t");
                    else
                      FP "\t%7.3lf\t", temp);
                  }
                }
                else
                {
                  FP "\t\t");
                }
              }
            }
            FP "\n");
          }
          if (ii == 2 || ii == 3)
          {
            FP " LLRtest ");
            for (i = 0; i < nummigrateparams; i++)
              if ((ii==2 && imig[i].b == k && imig[i].pr.max > MPRIORMIN ) ||
                  (ii==3 && imig[mterm[i]].b == k && imig[mterm[i]].pr.max > MPRIORMIN ))
              {
                if (fabs (temptest[i]) > 1e6)
                  sprintf (llrstring, "bad value\t");
                else
                {
                  if (temptest[i] > 9.54954)
                  {
                    printsigfootnote = 1;
                    sprintf (llrstring, "%7.3lf%s\t", temptest[i], sig[3]);
                  }
                  else if (temptest[i] > 5.41189)
                    sprintf (llrstring, "%7.3lf%s\t", temptest[i], sig[2]);
                  else if (temptest[i] > 2.70554)
                    sprintf (llrstring, "%7.3lf%s\t", temptest[i], sig[1]);
                  else
                    sprintf (&llrstring[0], "%7.3lf%s\t", temptest[i], sig[0]);
                }
                FP "%s", llrstring);
              }
            FP "\n");
          }
          if (ii == 1)
          {
            FP " LastPeriod");
            for (i = 0; i < numpopsizeparams; i++)
              if (itheta[i].b == k)
                FP "\t%d\t", itheta[i].e);
            FP "\n");
          }
          if (ii == 2 || ii==3)
          {
            FP " LastPeriod");
            for (i = 0; i < nummigrateparams; i++)
            {
              if (ii==2 && imig[i].b == k && imig[i].pr.max > MPRIORMIN )
                FP "\t%d\t", imig[i].e);
              if (ii==3 && imig[mterm[i]].b == k && imig[mterm[i]].pr.max > MPRIORMIN )
                FP "\t%d\t", imig[mterm[i]].e);
            }
            FP "\n");
          }
        }
      }
    }
    for (i = 0; i < p; i++)
      holdpeakloc[i] = (float)
        (peakloc[NUMTREEINT][i] < 0) ? (double) MINPARAMVAL : peakloc[NUMTREEINT][i];
  }
  if (printerrorfootnote==1)
    FP "*  peak not found possibly due to multiple peaks (check plot of marginal density) \n");
  if (printsigfootnote == 1)
  {
     FP " migration rate likelihood ratio test - see Nielsen and Wakeley (2001)\n migration significance levels :  * p < 0.05;   **  p < 0.01,   *** p < 0.001\n");
  }
  
  FP "\n");


  orig2d_free2D ((void **) mlval, NUMTREEINT + 1);
  orig2d_free2D ((void **) peakloc, NUMTREEINT + 1);
  XFREE (migtest);
  if (outputoptions[POPMIGPARAMHIST])
  {
    orig2d_free2D ((void **) popmigmlval, NUMTREEINT + 1);
    orig2d_free2D ((void **) popmigpeakloc, NUMTREEINT + 1);    
    XFREE (popmigtest);
    XFREE(mpop);
    XFREE(mterm);
    for (i=0;i<nummigrateparams;i++)
      XFREE(popmigstr[i]);
    XFREE(popmigstr);
  }
}                               /* findmarginpeaks */


