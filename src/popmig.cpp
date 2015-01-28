/*IMa2p 2009-2015 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, and Arun Sethuraman */

#undef GLOBVARS
#include "imamp.hpp"


/* calculate the probability density of the product 2NM  */
double
calc_popmig (int thetai, int mi, double x, int prob_or_like)
{
  int ei, cc, mc;
  int ccp, fcp, hcp, mcp, fmp, qip, mip;
  double sum, temp, temp1, temp2, fc, fm, hc, qintg, mintg, mmax, qmax;
  double a,b;
  ccp = gsamp_ccp + thetai;
  fcp = gsamp_fcp + thetai;
  hcp = gsamp_hccp + thetai;
  mcp = gsamp_mcp + mi;
  fmp = gsamp_fmp + mi;
  qip = gsamp_qip + thetai;
  mip = gsamp_mip + mi;
  mmax = imig[mi].pr.max;
  qmax = itheta[thetai].pr.max;

  for (sum = 0, ei = 0; ei < genealogiessaved; ei++)

  {
    cc = (int) gsampinf[ei][ccp];
    fc = (double) gsampinf[ei][fcp];
    hc = (double) gsampinf[ei][hcp];
    mc = (int) gsampinf[ei][mcp];
    fm = (double) gsampinf[ei][fmp];
    qintg = (double) gsampinf[ei][qip];
    mintg = (double) gsampinf[ei][mip];
    if (fc == 0 && cc == 0 && fm > 0)
    {
      temp1 = LOG2 - (mc * log (fm)) - hc - qintg - mintg;
      temp2 = log (exp (uppergamma (mc, 2 * fm * x / qmax)) - exp (uppergamma (mc, mmax * fm)));
    }
    else
    {
      if (fm == 0 && mc == 0 && fc > 0)
      {
        temp1 = LOG2 - (cc * log (fc)) - hc - qintg - mintg;
        temp2 = log (exp (uppergamma (cc, 2 * fc / qmax)) - exp (uppergamma (cc, fc * mmax / x)));
      }
      else
      {
        if (fc == 0 && cc == 0 && mc == 0 && fm == 0)
        {
          temp1 = log (2 * log (mmax * qmax / (2 * x))) - hc - qintg - mintg;;
          temp2 = 0;
        }
        else// fc>0, cc>0, mc> 0, fm> 0
        {
          temp1 = LOG2 + (mc * log (x)) - ((cc + mc) * log (fc + fm * x)) - hc - qintg - mintg;
          //temp2 = log ( exp(uppergamma (cc + mc, 2 * (fc + fm * x) / qmax)) - exp(uppergamma (cc + mc, mmax * (fm + fc / x))) );
          a = uppergamma (cc + mc, 2 * (fc + fm * x) / qmax);
          b = uppergamma (cc + mc, mmax * (fm + fc / x));
          /* the difference between two upppergamma()s is the negative of */
          /* the difference between two lowergamma()s.  However numbers   */
          /* for a or b often come out to be very large and equal under   */
          /* under uppergamma  or underlowergramma, in which case it is   */
          /* necessary to trap this and use the other function            */
          if (a==b)
          {
            b = lowergamma (cc + mc, 2 * (fc + fm * x) / qmax);
            a = lowergamma (cc + mc, mmax * (fm + fc / x));
          }
		  if (a>b)
          {
			LogDiff(temp2,a,b);
          }
		  else
          {
            temp1 = temp2 = 0.0;
          }
        }
      }
    }
    if ((temp1 + temp2 < 700 ) && (temp1 + temp2 > -700 )) // skip things that cannot be exped
      sum += exp (temp1 + temp2);
  }
  sum /= genealogiessaved;
  if (prob_or_like)
  {
    temp = 2 * (log (qmax) + log (mmax) - log (2 * x)) / (qmax * mmax);
    sum /= temp;
  }
  return sum;
}                               //calc_popmig

/* calculate the probability density of the product 2NM when migration has an exponential prior */
#define OCUTOFF  10
double
calc_pop_expomig (int thetai, int mi, double x, int prob_or_like)
{
  int ei, cc, mc;
  int ccp, fcp, hcp, mcp, fmp, qip, mip;
  double sum, temp, temp1, temp2, temp3,fc, fm, hc, qintg, mintg, mmean, qmax;
  int zadj, maxz = -1000000000;
  double acumm = 0;
  ccp = gsamp_ccp + thetai;
  fcp = gsamp_fcp + thetai;
  hcp = gsamp_hccp + thetai;
  mcp = gsamp_mcp + mi;
  fmp = gsamp_fmp + mi;
  qip = gsamp_qip + thetai;
  mip = gsamp_mip + mi;
  mmean = imig[mi].pr.mean;
  qmax = itheta[thetai].pr.max;

  for (sum = 0, ei = 0; ei < genealogiessaved; ei++)

  {
    cc = (int) gsampinf[ei][ccp];
    fc = (double) gsampinf[ei][fcp];
    hc = (double) gsampinf[ei][hcp];
    mc = (int) gsampinf[ei][mcp];
    fm = (double) gsampinf[ei][fmp];
    qintg = (double) gsampinf[ei][qip];
    mintg = (double) gsampinf[ei][mip];
    temp1 = x + fc *mmean + fm * mmean * x;
    temp2 = LOG2 - hc - log(mmean) - qintg - mintg;
    temp3 = 2 * temp1/(mmean * qmax);

    if (mc == 0 && cc == 0 )
    {
      temp3 = uppergamma(0,temp3);
      temp2 += temp3;
    }
    else
    {
      if (cc == 0) // mc > 0 
      {
        temp3 = uppergamma(mc,temp3);
        temp2 += temp3 - mc* log(fm + 1/mmean);
      }
      else
      {
        if (mc == 0) //cc > 0
        {
          temp3 = uppergamma(cc,temp3);
          temp2 += temp3 + -cc *log(fc + x*(fm + 1/mmean));; 
        }
        else
        {
          temp3 = uppergamma(mc+cc,temp3);
          temp2 += temp3 - cc*log(x) + (cc+mc)*log(mmean*x/temp1); 
        }
      }
    }
    eexp (temp2, &eexpsum[ei].m, &eexpsum[ei].z);
    if (eexpsum[ei].z > maxz)
      maxz = eexpsum[ei].z;
  }
  for (ei = 0; ei < genealogiessaved; ei++)
  {
    zadj = eexpsum[ei].z - (maxz - OCUTOFF);
    eexpsum[ei].m *= pow (10.0, (double) zadj);
    acumm += eexpsum[ei].m;
  }
  sum = log (acumm) + (maxz - OCUTOFF) * LOG10;
  sum -= log((double) genealogiessaved);
  sum = exp(sum);
  if (prob_or_like)
  {
    temp = 2 * exp(uppergamma(0,2*x/(mmean*qmax)))/ (qmax * mmean);
    sum /= temp;
  }
  return sum;
}                               //calc_pop_expomig
#undef OCUTOFF  

#define OFFSCALEVAL 1.0

/*  do calculations for finding the peak of the marginal density for a 2NM term
    this is similar to calc_popmig().  It is called by marginalopt_popmig()*/

double
marginpopmig (int mi, int firsttree, int lasttree, double x, int thetai)
{
  int ei, cc, mc;
  int ccp, fcp, hcp, mcp, fmp, qip, mip;
  double sum, temp1, temp2, fc, fm, hc, qintg, mintg, mmax, qmax;
  double a,b;
  double max, min;

  max = itheta[thetai].pr.max*imig[mi].pr.max/2.0;
  min = 0;
  if (x < min || x > max)
    return OFFSCALEVAL;
  
  ccp = gsamp_ccp + thetai;
  fcp = gsamp_fcp + thetai;
  hcp = gsamp_hccp + thetai;
  mcp = gsamp_mcp + mi;
  fmp = gsamp_fmp + mi;
  qip = gsamp_qip + thetai;
  mip = gsamp_mip + mi;
  mmax = imig[mi].pr.max;
  qmax = itheta[thetai].pr.max;

  for (sum = 0,ei = firsttree; ei < lasttree; ei++)

  {
    cc = (int) gsampinf[ei][ccp];
    fc = (double) gsampinf[ei][fcp];
    hc = (double) gsampinf[ei][hcp];
    mc = (int) gsampinf[ei][mcp];
    fm = (double) gsampinf[ei][fmp];
    qintg = (double) gsampinf[ei][qip];
    mintg = (double) gsampinf[ei][mip];
    if (fc == 0 && cc == 0 && fm > 0)
    {
      temp1 = LOG2 - (mc * log (fm)) - hc - qintg - mintg;
      temp2 = log (exp (uppergamma (mc, 2 * fm * x / qmax)) - exp (uppergamma (mc, mmax * fm)));
    }
    else
    {
      if (fm == 0 && mc == 0 && fc > 0)
      {
        temp1 = LOG2 - (cc * log (fc)) - hc - qintg - mintg;
        temp2 = log (exp (uppergamma (cc, 2 * fc / qmax)) - exp (uppergamma (cc, fc * mmax / x)));
      }
      else
      {
        if (fc == 0 && cc == 0 && mc == 0 && fm == 0)
        {
          temp1 = log (2 * log (mmax * qmax / (2 * x))) - hc - qintg - mintg;;
          temp2 = 0;
        }
        else
        {
          temp1 = LOG2 + (mc * log (x)) - ((cc + mc) * log (fc + fm * x)) - hc - qintg - mintg;
          a = uppergamma (cc + mc, 2 * (fc + fm * x) / qmax);
          b = uppergamma (cc + mc, mmax * (fm + fc / x));
          /* the difference between two upppergamma()s is the negative of */
          /* the difference between two lowergamma()s.  However numbers   */
          /* for a or b often come out to be very large and equal under   */
          /* under uppergamma  or underlowergramma, in which case it is   */
          /* necessary to trap this and use the other function            */
          if (a==b)
          {
            b = lowergamma (cc + mc, 2 * (fc + fm * x) / qmax);
            a = lowergamma (cc + mc, mmax * (fm + fc / x));
          }
		  if (a>b)
          {
			LogDiff(temp2,a,b);
          }
		  else
          {
            temp1 = temp2 = 0.0;
          }

          //temp2 = log (exp (uppergamma (cc + mc, 2 * (fc + fm * x) / qmax)) - exp (uppergamma (cc + mc, mmax * (fm + fc / x))));
        }
      }
    }
    if ((temp1 + temp2 < 700 ) && (temp1 + temp2 > -700 )) // skip things that cannot be exped
      sum += exp (temp1 + temp2);
  }
  sum /= (lasttree - firsttree + (firsttree == 0)); 
  return -sum; 
}                               /* marginpopmig */


/* calculate the probability density of the product 2NM when migration has an exponential prior */
#define OCUTOFF  10
double
marginpop_expomig (int mi, int firsttree, int lasttree, double x, int thetai)
{
  int ei, cc, mc;
  int ccp, fcp, hcp, mcp, fmp, qip, mip;
  double sum, temp1, temp2, temp3,fc, fm, hc, qintg, mintg, mmean, qmax;
  int zadj, maxz = -1000000000;
  double acumm = 0;
  double max, min;


  min = 0;
  max = EXPOMIGPLOTSCALE * imig[mi].pr.mean;
  if (x < min || x > max)
    return OFFSCALEVAL;

  ccp = gsamp_ccp + thetai;
  fcp = gsamp_fcp + thetai;
  hcp = gsamp_hccp + thetai;
  mcp = gsamp_mcp + mi;
  fmp = gsamp_fmp + mi;
  qip = gsamp_qip + thetai;
  mip = gsamp_mip + mi;
  mmean = imig[mi].pr.mean;
  qmax = itheta[thetai].pr.max;

  for (sum = 0, ei = firsttree; ei < lasttree; ei++)

  {
    cc = (int) gsampinf[ei][ccp];
    fc = (double) gsampinf[ei][fcp];
    hc = (double) gsampinf[ei][hcp];
    mc = (int) gsampinf[ei][mcp];
    fm = (double) gsampinf[ei][fmp];
    qintg = (double) gsampinf[ei][qip];
    mintg = (double) gsampinf[ei][mip];
    temp1 = x + fc *mmean + fm * mmean * x;
    temp2 = LOG2 - hc - log(mmean) - qintg - mintg;
    temp3 = 2 * temp1/(mmean * qmax);

    if (mc == 0 && cc == 0 )
    {
      temp3 = uppergamma(0,temp3);
      temp2 += temp3;
    }
    else
    {
      if (cc == 0) // mc > 0 
      {
        temp3 = uppergamma(mc,temp3);
        temp2 += temp3 - mc* log(fm + 1/mmean);
      }
      else
      {
        if (mc == 0) //cc > 0
        {
          temp3 = uppergamma(cc,temp3);
          temp2 += temp3 + -cc *log(fc + x*(fm + 1/mmean));; 
        }
        else
        {
          temp3 = uppergamma(mc+cc,temp3);
          temp2 += temp3 - cc*log(x) + (cc+mc)*log(mmean*x/temp1); 
        }
      }
    }
    eexp (temp2, &eexpsum[ei].m, &eexpsum[ei].z);
    if (eexpsum[ei].z > maxz)
      maxz = eexpsum[ei].z;
  }
  for (ei = firsttree; ei < lasttree; ei++)
  {
    zadj = eexpsum[ei].z - (maxz - OCUTOFF);
    eexpsum[ei].m *= pow (10.0, (double) zadj);
    acumm += eexpsum[ei].m;
  }
  sum = log (acumm) + (maxz - OCUTOFF) * LOG10;
  sum -= log(static_cast<double>(lasttree) - firsttree + (firsttree == 0));
  sum = exp(sum);
  return -sum;
}                               //marginpop_expomig
#undef OCUTOFF  
#undef OFFSCALEVAL

#define SEARCHSTARTFRAC  4      // fraction of position in parameter range to start at
/* find the peak of the marginal density for a 2NM term.  This is similar to the marginalopt() function
  this calls marginpopmig() */ 

#define POPMIG_MARGINBINS 100
void
marginalopt_popmig (int firsttree, int lasttree, double *mlval, double *peakloc, int *mpop, int *mterm)
{
  int i,j;
  double ftol, ax, bx, cx, fa, xmax, ml;
  double upperbound;
  double maxf;
  int maxj;
  double (*func) (int, int, int, double, int);
  ftol = 1e-7;

  if (modeloptions[EXPOMIGRATIONPRIOR])
    func = marginpop_expomig;
  else
    func = marginpopmig;
  for (i = 0; i < nummigrateparams; i++)
  {
    if (modeloptions[EXPOMIGRATIONPRIOR])
      upperbound = EXPOMIGPLOTSCALE * imig[mterm[i]].pr.mean;
    else
      upperbound = itheta[mpop[i]].pr.max * imig[mterm[i]].pr.max/2.0;
/* 8/31/09  JH
  it turns out these curves often have multiple peaks, and mnbrak just was not working then.
  replace use of mnbrak with function calls over POPMIG_MARGINBINS to find an interval with 
  a peak.  This is quite crude */
    ax = MINPARAMVAL;
    maxf = 1e100;
    maxj = -1;
    for (j = 0;j< POPMIG_MARGINBINS;j++)
    {
      ax = MINPARAMVAL + j*(upperbound /(double) POPMIG_MARGINBINS);
      fa = func (mterm[i], firsttree, lasttree, ax, mpop[i]);
//printf("j  %d  ax  %lf  fa  %lf \n",j, ax, fa);
      if (fa < maxf)
      {
        maxf = fa;
        maxj = j;
      }
    }
//printf("i %d  maxf  %lf   maxj  %d\n",i,maxf, maxj);
    if (maxj == 0)
    {
      ax  =  MINPARAMVAL;
      cx = MINPARAMVAL + upperbound /(double) POPMIG_MARGINBINS;
      bx = (ax + cx)/2.0;
    }
    else
    {
      if (maxj == POPMIG_MARGINBINS -1)
      {
        ax = MINPARAMVAL + ((double) POPMIG_MARGINBINS - 1)*(upperbound /(double) POPMIG_MARGINBINS);
        bx = MINPARAMVAL + ((double) POPMIG_MARGINBINS - 2)*(upperbound /(double) POPMIG_MARGINBINS);
        cx = upperbound;
      }
      else
      {
        bx = MINPARAMVAL + ((double) maxj - 1)*(upperbound /(double) POPMIG_MARGINBINS);
        cx = MINPARAMVAL + ((double) maxj + 1)*(upperbound /(double) POPMIG_MARGINBINS);
        ax = MINPARAMVAL + ((double) maxj)*(upperbound /(double) POPMIG_MARGINBINS);
      }
    }
//printf("ax %lf  bx %lf  cx %lf\n",ax,bx,cx);
    //{
      ml = -goldenmod (mterm[i], firsttree, lasttree, ax, bx, cx, ftol, &xmax, func,mpop[i]);
      mlval[i] = ml;
      peakloc[i] = xmax;
    //}
//printf("ml  %lf  xmax %lf\n",ml, xmax);
  }
}                               /* marginalopt_popmig */
