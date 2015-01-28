/*IMa2p 2009-2015 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, and Arun Sethuraman */
#undef GLOBVARS
#include "imamp.hpp"
#include "updateassignment.hpp"


/* stuff for writing and reading the state of the mcmc to a file 
to write the state of the mcmc state space to a file
call writemcmc(int ci)
to load it from a file 
call readmcmc(int ci)
after calling readme  should probably call treeweight() and integrate_tree_prob () to do the main likelihood and probability math
e.g. for infinite sites data from chain 0
the program IMa() has similar functions
can write an mcmc statespace in IMa() to a file and then read it in IMamp() to compare programs (for 2 population models)

codes for atype 
0  for integer
1  for long
2  for float
3  for double
4  for character
5  for unsigned long
*/


static void awrite (FILE * mcffile, const char *name, int atype, int iu, void *a);
static void aread (FILE * mcffile, const char *name, int atype, int iucheck,
                   void *a);
static void init_p (void);
/**** LOCAL FUNCTIONS *****/
void
awrite (FILE * mcffile, const char *name, int atype, int iu, void *a)
{
  int *ip;
  long *lip;
  float *fp;
  double *dp;
  char *cp;
  unsigned long *ulp;
  int i;

  fprintf (mcffile, "%s %d %d ", name, atype, iu);
  switch (atype)
  {
  case 0:
    for (i = 0, ip = static_cast<int *> (a); i < iu; i++)
      fprintf (mcffile, "%d ", *(ip + i));
    break;
  case 1:
    for (i = 0, lip = static_cast<long *> (a); i < iu; i++)
      fprintf (mcffile, "%ld ", *(lip + i));
    break;
  case 2:
    for (i = 0, fp = static_cast<float *> (a); i < iu; i++)
      fprintf (mcffile, "%g ", *(fp + i));
    break;
  case 3:
    for (i = 0, dp = static_cast<double *> (a); i < iu; i++)
      fprintf (mcffile, "%.10lg ", *(dp + i));
    break;
  case 4:
    for (i = 0, cp = static_cast<char *> (a); i < iu; i++)
      fprintf (mcffile, "%c", *(cp + i));
    if (*(cp + i) != '\0')
      fprintf (mcffile, "%c", '\0');
    break;
  case 5:
    for (i = 0, ulp = static_cast<unsigned long *> (a); i < iu; i++)
      fprintf (mcffile, "%lu ", *(ulp + i));
    break;
  }
  fprintf (mcffile, "\n");
}                               /* arraywrite */


void
aread (FILE * mcffile, const char *name, int atype, int iucheck, void *a)
{
  int *ip;
  long *lip;
  float *fp;
  double *dp;
  char *cp;
  unsigned long *ulp;
  int typecheck, i, iu;
  char namecheck[100];

  iucheck = 0; /* To remove warning: unused parameter. */
  fscanf (mcffile, "%99s %d %d ", namecheck, &typecheck, &iu);
  if (strcmp (name, namecheck) != 0)
  {
    IM_err(IMERR_MCFREADFAIL,"variable names do not match: %s  <> %s ", name, namecheck);
  }
  if (typecheck != atype)
  {
    IM_err(IMERR_MCFREADFAIL,"variable types do not match: %d  <> %d ", atype, typecheck);
  }
  switch (atype)
  {
  case 0:
    for (i = 0, ip = static_cast<int *> (a); i < iu; i++, ip++)
      fscanf (mcffile, "%d ", ip);
    break;
  case 1:
    for (i = 0, lip = static_cast<long *> (a); i < iu; i++, lip++)
      fscanf (mcffile, "%ld ", lip);
    break;
  case 2:
    for (i = 0, fp = static_cast<float *> (a); i < iu; i++, fp++)
      fscanf (mcffile, "%g ", fp);
    break;
  case 3:
    for (i = 0, dp = static_cast<double *> (a); i < iu; i++, dp++)
      fscanf (mcffile, "%lg ", dp);
    break;
  case 4:
    for (i = 0, cp = static_cast<char *> (a); i < iu; i++, cp++)
      fscanf (mcffile, "%c", cp);
    *cp = '\0';
    break;
  case 5:
    for (i = 0, ulp = static_cast<unsigned long *> (a); i < iu; i++, ulp++)
      fscanf (mcffile, "%lu ", ulp);
    break;
  }
  return;
}                               /* arrayread */

void
init_p (void)
{
  int ci, li, ai;
  for (ci = 0; ci < numchains; ci++)
  {
    setzero_genealogy_weights (&C[ci]->allgweight);
    C[ci]->allpcalc.pdg = 0;
    for (li = 0; li < nloci; li++)
    {
      setzero_genealogy_weights (&C[ci]->G[li].gweight);
      treeweight (ci, li);
      switch (L[li].model)
      {
      case INFINITESITES:
        C[ci]->G[li].pdg = C[ci]->G[li].pdg_a[0] =
          likelihoodIS (ci, li, C[ci]->G[li].uvals[0]);
        break;
      case HKY:
        if (assignmentoptions[JCMODEL] == 1)
        {
          C[ci]->G[li].pdg = C[ci]->G[li].pdg_a[0] =
            likelihoodJC (ci, li, C[ci]->G[li].uvals[0]);
        }
        else
        {
          C[ci]->G[li].pdg = C[ci]->G[li].pdg_a[0] =
            likelihoodHKY (ci, li, C[ci]->G[li].uvals[0],
                           C[ci]->G[li].kappaval, -1, -1, -1, -1);
          copyfraclike (ci, li);
          storescalefactors (ci, li);
        }
        break;
      case STEPWISE:
        somestepwise = 1;
        C[ci]->G[li].pdg = 0;
        for (ai = 0; ai < L[li].nlinked; ai++)
        {
          C[ci]->G[li].pdg_a[ai] = likelihoodSW (ci, li, ai, C[ci]->G[li].uvals[ai], 1.0);
          C[ci]->G[li].pdg += C[ci]->G[li].pdg_a[ai];
        }
        break;
      case JOINT_IS_SW:
        somestepwise = 1;
        // JH  fixed bug here on 6/1/09  these two lines should not have been included here 
        //makeJOINT_IS_SW (ci, li);
        //treeweight (ci, li);
        C[ci]->G[li].pdg = C[ci]->G[li].pdg_a[0] =
          likelihoodIS (ci, li, C[ci]->G[li].uvals[0]);
        for (ai = 1; ai < L[li].nlinked; ai++)
        {
          C[ci]->G[li].pdg_a[ai] = likelihoodSW (ci, li, ai, C[ci]->G[li].uvals[ai], 1.0);
          C[ci]->G[li].pdg += C[ci]->G[li].pdg_a[ai];
        }
        break;
      }
      sum_treeinfo (&(C[ci]->allgweight), &(C[ci]->G[li].gweight));
      C[ci]->allpcalc.pdg += C[ci]->G[li].pdg;
    }
    initialize_integrate_tree_prob (ci, &(C[ci]->allgweight),
                                    &C[ci]->allpcalc);
  }

}                               // init_p

/***********GLOBAL FUNCTIONS **********/
#ifdef aa
#undef aa
#endif /* aa */
#define aa   awrite (mcffile,
/* the order of things written to the dck file is nearly arbitrary.  file can be large */


void
writemcf (char mcffilename[])
{
  int i, j, li, ci;
  FILE *mcffile;

  if ((mcffile = fopen (mcffilename, "w")) == NULL)
  {
    IM_err(IMERR_CREATEFILEFAIL,"Error creating mcffile: %s", mcffilename);
  }
  for (ci = 0; ci < numchains; ci++)
  {
    // tvalues 
    for (i = 0; i < numsplittimes; i++)
      aa "tvalue", 3, 1, &(C[ci]->tvals[i]));

    //mutation scalars 
    for (li = 0; li < nloci; li++)
    {
      if (nloci > 1 || L[li].nlinked > 0)
      {
        for (i = 0; i < L[li].nlinked; i++)
          aa "uvalue", 3, 1, &(C[ci]->G[li].uvals[i]));
      }
      if (L[li].model == HKY)
        aa "kappavalue", 3, 1, &(C[ci]->G[li].kappaval));
      aa "pi[4]", 3, 4, &(C[ci]->G[li].pi[0]));
      for (i = 0; i < L[li].numlines; i++)
      {
        aa "up[2]", 0, 2, &(C[ci]->G[li].gtree[i].up[0]));
        aa "down", 0, 1, &(C[ci]->G[li].gtree[i].down));
        aa "mut", 0, 1, &(C[ci]->G[li].gtree[i].mut));
        aa "pop", 0, 1, &(C[ci]->G[li].gtree[i].pop));
        if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
        {
          aa "A[]", 0, L[li].nlinked, &(C[ci]->G[li].gtree[i].A[0]));
          aa "dlikeA[]", 3, L[li].nlinked, &(C[ci]->G[li].gtree[i].dlikeA[0]));
        }
        aa "time", 3, 1, &(C[ci]->G[li].gtree[i].time));
        //aa "cmm",0,1,&(C[ci]->G[li].gtree[i].cmm));
        j = 0;
        do
        {
          aa "mig[].mt", 3, 1, &(C[ci]->G[li].gtree[i].mig[j].mt));
          if (C[ci]->G[li].gtree[i].mig[j].mt > 0)
            aa "mig[].mp", 0, 1, &(C[ci]->G[li].gtree[i].mig[j].mp));
          j++;
        }
        while (C[ci]->G[li].gtree[i].mig[j - 1].mt > 0);
        if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
        {
          aa "asn", 0, 1, &(C[ci]->G[li].gtree[i].pop)); 
        }

        if (L[li].model == HKY && i >= L[li].numgenes)
        {
          /* don't think need to save scalefactors, as they get recalculated by makefrac 
             aa "C[ci]->G[li].gtree[i].scalefactor[0]",3,L[li].numsites,&(C[ci]->G[li].gtree[i].scalefactor[0]));
             aa "C[ci]->G[li].gtree[i].oldscalefactor[0]",3,L[li].numsites,&(C[ci]->G[li].gtree[i].oldscalefactor[0]));  */
          for (j = 0; j < L[li].numsites; j++)
          {
            aa "C[ci]->G[li].gtree[i].hkyi.frac[j]", 3, 4, &(C[ci]->G[li].gtree[i].hkyi.frac[j][0]));
            aa "C[ci]->G[li].gtree[i].hkyi.newfrac[j]", 3, 4, &(C[ci]->G[li].gtree[i].hkyi.newfrac[j][0]));
          }
        }
      }
      /* should not be necessary as these get initialized anyway
         if (L[li].model == HKY || L[li].model == INFINITESITES
         || L[li].model == JOINT_IS_SW)
         {
         for (i = 0; i < L[li].numgenes; i++)
         {
         aa "L[li].seq[i]", 0, L[li].numsites,
         &(L[li].seq[i][0]));
         }
         }
         if (L[li].model == HKY)
         {
         aa "L[li].mult", 0, L[li].numsites,
         &(L[li].mult[0]));
         }
         if (L[li].model == INFINITESITES)
         {
         aa "L[li].badsite", 0, L[li].numbases,
         &(L[li].badsite[0]));
         } */
    }
  }
  fclose (mcffile);
}                               /* writemcf */


#ifdef aa
#undef aa
#endif /* aa */
#define aa  aread(mcffile,
/*readmcf() is very similar to writemcf(), except it includes some mallocs() and a small number of initializations at the end*/

/* this reads the mcffile
if the number of chains in the current run (call this cc for now) is different than the number of chains in the mcffile (call this cf) then:
- if cc <= cf,  then the first cc chains in the mcffile are loaded
- if cc > cf,  then the first cf chains in the mcffile are loaded into the first cf positions in C[]
	then the mcffile is closed and reopened and the chains are reloaded in to the next available positions
	this keeps getting done until all the cc positions in C[] have been loaded  */
/* need to revise this so that if cc > cf,  then all chains above cf are copied from the chain cf
  this way we are less likely to have a chain that is much better than its beta value, 
  at least compared to what happens when we start copying low number chains into highly heated positions */ 
void readmcf (char mcffilename[])
{
  int i, j, li, ci, lastci = -1;
  double uptime;
  FILE *mcffile;
  char checkeofc;
  int largetimeflag;

  if ((mcffile = fopen (mcffilename, "r")) == NULL)
  {
    IM_err(IMERR_READFILEOPENFAIL,"Error opening mcffile: %s", mcffilename);
  }


  for (ci = 0; ci < numchains; ci++)
  {
// tvalues 
    largetimeflag = 0;
    for (i = 0; i < numsplittimes; i++)
    {
      aa "tvalue", 3, 1, &(C[ci]->tvals[i]));
      largetimeflag = largetimeflag || (C[ci]->tvals[i] > T[i].pr.max);
      //assert(C[ci]->tvals[i] < T[i].pr.max);
      C[ci]->poptree[C[ci]->droppops[i + 1][0]].time =
        C[ci]->poptree[C[ci]->droppops[i + 1][1]].time = C[ci]->tvals[i];
    }

    //mutation scalars 
    for (li = 0; li < nloci; li++)
    {
      if (nloci > 1 || L[li].nlinked > 0)
      {
        for (i = 0; i < L[li].nlinked; i++)
          aa "uvalue", 3, 1, &(C[ci]->G[li].uvals[i]));
      }
      if (L[li].model == HKY)
        aa "kappavalue", 3, 1, &(C[ci]->G[li].kappaval));
      aa "pi[4]", 3, 4, &(C[ci]->G[li].pi[0]));
      C[ci]->G[li].mignum = 0;
      C[ci]->G[li].tlength = 0;
      C[ci]->G[li].length = 0;
      for (i = 0; i < L[li].numlines; i++)
      {
        aa "up[2]", 0, 2, &(C[ci]->G[li].gtree[i].up[0]));
        aa "down", 0, 1, &(C[ci]->G[li].gtree[i].down));
        if (C[ci]->G[li].gtree[i].down == -1)
        {
          C[ci]->G[li].root = i;
        }

        aa "mut", 0, 1, &(C[ci]->G[li].gtree[i].mut));
        aa "pop", 0, 1, &(C[ci]->G[li].gtree[i].pop));
        if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
        {
          aa "A[]", 0, L[li].nlinked, &(C[ci]->G[li].gtree[i].A[0]));
          aa "dlikeA[]", 3, L[li].nlinked,
            &(C[ci]->G[li].gtree[i].dlikeA[0]));
        }
        aa "time", 3, 1, &(C[ci]->G[li].gtree[i].time));
        if (i != C[ci]->G[li].root)
        {
          if (i < L[li].numgenes)
            uptime = 0;
          else
            uptime = C[ci]->G[li].gtree[C[ci]->G[li].gtree[i].up[0]].time;

          C[ci]->G[li].tlength += C[ci]->G[li].gtree[i].time - uptime;
          if (uptime < C[ci]->tvals[npops - 1])
            C[ci]->G[li].tlength +=
              DMIN (C[ci]->G[li].gtree[i].time - uptime,
                    C[ci]->tvals[npops - 1]) - uptime;
        }

        //aa "cmm",0,1,&(C[ci]->G[li].gtree[i].cmm));  don't read this, but reset current value if neeed
        j = 0;
        do
        {
          aa "mig[].mt", 3, 1, &(C[ci]->G[li].gtree[i].mig[j].mt));
          if (C[ci]->G[li].gtree[i].mig[j].mt > 0)
          {
            aa "mig[].mp", 0, 1, &(C[ci]->G[li].gtree[i].mig[j].mp));
            C[ci]->G[li].mignum++;
          }
          j++;
          checkmig (j, &(C[ci]->G[li].gtree[i].mig),
                    &(C[ci]->G[li].gtree[i].cmm));
        } while (C[ci]->G[li].gtree[i].mig[j - 1].mt > 0);
        if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
        {
          aa "asn", 0, 1, &(C[ci]->G[li].gtree[i].pop)); 
        }

        if (L[li].model == HKY && i >= L[li].numgenes)
        {
          /* don't think need to save scalefactors, as they get recalculated by makefrac 
             aa "C[ci]->G[li].gtree[i].scalefactor[0]",3,L[li].numsites,&(C[ci]->G[li].gtree[i].scalefactor[0]));
             aa "C[ci]->G[li].gtree[i].oldscalefactor[0]",3,L[li].numsites,&(C[ci]->G[li].gtree[i].oldscalefactor[0]));  */
          for (j = 0; j < L[li].numsites; j++)
          {
            aa "C[ci]->G[li].gtree[i].hkyi.frac[j]", 3, 4,
              &(C[ci]->G[li].gtree[i].hkyi.frac[j][0]));
            aa "C[ci]->G[li].gtree[i].hkyi.newfrac[j]", 3, 4,
              &(C[ci]->G[li].gtree[i].hkyi.newfrac[j][0]));
          }
        }
      }
      C[ci]->G[li].roottime = C[ci]->G[li].gtree[C[ci]->G[li].gtree[C[ci]->G[li].root].up[0]].time;
    }
    if ((checkeofc = (char) fgetc (mcffile)) == EOF)
    {
      if (ci == lastci)  
        IM_err (IMERR_MCFSPLITTIMEPROB,
                " can't load trees, probably because of multiple instances of splittime time conflict with t prior");
      if (ci < numchains - 1)
      {
        f_close (mcffile);
        if ((mcffile = fopen (mcffilename, "r")) == NULL)
        {
           IM_err(IMERR_READFILEOPENFAIL,"Error reopening mcffile: %s", mcffilename);
        }
        lastci = ci;
      }
    }
    else
    {
      ungetc (checkeofc, mcffile);
    }
    if (largetimeflag)  // skip that locus
      ci--;
  }
  fclose (mcffile);
  init_p ();
}                               // readmcf
