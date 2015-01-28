/*IMa2p 2009-2015 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, and Arun Sethuraman */
#undef GLOBVARS
#include "imamp.hpp"
/* calculate marginal likelihood */


/*
implement thermodynamic integration over a large number of intervals

we want the marginal likelihood under the model

p(D) 

for p_Bi(D|G)

where Bi is a heating value i,  for a total of  j values   0<i<j-1

Let L_i  be the mean of p_Bi(D|G) sampled over the course of the run

L_i = Sum[p_Bi(D|G)]/k    for k samples

Then 

p(d) = Sum[ (Bi-B(i-1)) (L_i + L_(i-1))/2, {i, 1, j-1}] 

this is trapezoidal rule 

Also record the harmonic mean  - have to use eexp()

p(d) = 1/ [ SUM[ 1/(p(D|G)]/k ]

*/
#define LOG_10_2  0.30102999566398119521
#define OCUTOFF  10

double thermosum[MAXCHAINS]; 



void initmarginlikecalc()
{
  int i;
  //harmonicsum = 0.0;
 // harmonicsum_eexp = 0.0;
  for (i=0;i<numchains;i++)
    thermosum[i]  = 0.0;
}

/* harmonic mean calculations,  has to deal with using eexp() */

void summarginlikecalc(void)
{
  int i;

  /*  10/6/2011  commented out some old code for calculation harmonic mean   
    this was moved to harmonicmarginlikecalc() 
  static int ei = 0;
  int tempz;
  int i, zadj;
  double tempm;

  //harmonicsum += 1.0/exp(C[0]->allpcalc.pdg); simple form has floating point problems
  eexp(C[0]->allpcalc.pdg,&tempm,&tempz);
  if (ei < HARMONICMEANCHECK)
  {
    harmonicsump[ei].m = 1.0/tempm;
    harmonicsump[ei].z = -tempz;
    if (harmonicsump[ei].z > maxz)
      maxz = harmonicsump[ei].z;
  }
  else
  {
    if (ei == HARMONICMEANCHECK)
    {
      for (i = 0; i < ei; i++)
      {
        zadj = harmonicsump[i].z - (maxz - OCUTOFF);
        harmonicsum_eexp += harmonicsump[i].m * pow (10.0, (double) zadj);
      }
    }
    zadj = -tempz - (maxz - OCUTOFF);
    harmonicsum_eexp += (1.0/tempm) * pow (10.0, (double) zadj);
  }
  ei++; */
  for (i=0;i<numchains;i++)
    thermosum[i] += C[i]->allpcalc.pdg; 
}

double harmonicmarginlikecalc(int k)
{
  double hmlog;
  int gi, zadj;
  int pdgp;
  int tempz;
  double tempm;
  struct extendnum *harmonicsump;
  int maxz = 0;
  double  harmonicsum_eexp = 0.0;

  harmonicsump = (struct extendnum *) malloc ((size_t) ((genealogiessaved + 1) * sizeof (struct extendnum)));

  pdgp = 4*numpopsizeparams + 3* nummigrateparams;  // position in gsampinf[gi] that holds pdg
  for (gi = 0; gi < genealogiessaved; gi++)
  {
    eexp(gsampinf[gi][pdgp],&tempm,&tempz);
    harmonicsump[gi].m = 1.0/tempm;
    harmonicsump[gi].z = -tempz;
    if (harmonicsump[gi].z > maxz)
      maxz = harmonicsump[gi].z;
  }
  for (gi = 0; gi < genealogiessaved; gi++)
  {
    zadj = harmonicsump[gi].z - (maxz - OCUTOFF);
    harmonicsum_eexp += harmonicsump[gi].m * pow (10.0, (double) zadj);
  }
  hmlog = -(log (harmonicsum_eexp) - log( (double) genealogiessaved) + (maxz - OCUTOFF) * LOG10);
  XFREE(harmonicsump);
  return hmlog;
}
/* made this change based on audes recommendation 9/16/2011 */ 
double thermomarginlikecalc(int k)
{
  int i;
  double width, sum; 
  double mcalc = 0.0;
  /* trapezoid rule 
  for (i= (numchains - 2); i>= 0;i--)
  {
    if (beta[i+1]==0.0)
      mcalc += (beta[i+1]- beta[i]) * (thermosum[i]/k);
    mcalc += (beta[i+1]- beta[i]) * ((thermosum[i]/k + thermosum[i+1]/k)/2);
  } */
  
  /* Simpson's rule  - must ensure previously that the number of intervals is even*/
  /* so the number of chains must  be odd.
    but since we are counting from zero,  the value of numchains must  even */
  width = 1.0/ (float) (numchains -1);
  sum = 0.0;
  for (i = 0;i<=numchains-1;i+=2)
  {
    if (i != (numchains - 1))  // for chain[numchains-1] use a value of 0.0 because under thermodynamic integration the mean likelihood at Beta=0 is 0.0 
      sum += 4.0 * (thermosum[i]/k);
  }
  for (i = 1;i<=numchains-2;i+=2)
  {
    sum += 2.0 * (thermosum[i]/k);
  }
  mcalc = width * sum/3.0; 
  return mcalc;
}


