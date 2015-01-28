/*IMa2p 2009-2015 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, and Arun Sethuraman */
#undef GLOBVARS
#include "imamp.hpp"
#include "updateassignment.hpp"
#include "update_gtree_common.hpp"

/* 1/8/09 rewrote this using new data structures   */
/* Most initializations  */

/*********** LOCAL STUFF **********/

extern double pi[MAXLOCI][4];
extern int urri[2 * MAXLOCI][2 * MAXLOCI];      // used mostly in update_mc_params.c
extern double urrlow[2 * MAXLOCI][2 * MAXLOCI], urrhi[2 * MAXLOCI][2 * MAXLOCI];        // used mostly in update_mc_params.c

static double **uvals;
static char startpoptreestring[POPTREESTRINGLENGTHMAX];
static double geomeanvar;
static int **numsitesIS;        // maxpops * maxloci   temporarily holds the number of polymorphic sites in loci with infinite sites model
static double uterm;
static double *popsizeprior_fromfile,**mprior_fromfile;
static int maxpossiblemigrateparams;
/* used if modeloptions[ANCESTRALPOPSIZESSETTOLARGESTSAMPLED]==1 
  varrank[] holds the sum of 4Nu estimates across loci for the sampled populations - fill in set_uvals()
  varranktree[i]:
    if i< npops,  then it is just the sampled population number and its corresponding value of varrank[i]
    else  it holds the sampled population number for which the population size parameter for this ancestral population will be the same as

*/
static double varrank[MAXPOPS] = {0.0};
static struct {
  int samppopnum;
  double var; 
}  varranktree[2*MAXPOPS - 1];


// 1/12/09  why is this stuff on struct edgemiginfo here  ??  what does it do??

/* We initialize these four in function [[IMA_initmemory_aftersetpopntree]]. */
/* global variables of update_gtree */
extern struct edgemiginfo oldedgemig;
extern struct edgemiginfo oldsismig;
extern struct edgemiginfo newedgemig;
extern struct edgemiginfo newsismig;

/* prototypes */
static void addpopsizeidtext(int *fpstri, char fpstr[]);
static void fillvarranktree(int nowpop);
static int numvarHKY (int li, int b, int e);
static int checkaresis (int period, int i, int j);
static void setimigarray (int ii);
static void setup_iparams ();
static void set_nomigrationchecklist ();
static void setuinfo (double summut);
static double set_uvals (void);
static void set_x (struct value_record *v, int isint);
static void init_value_record (struct value_record *v, int isint);
static void init_g_rec (int li);
static void init_a_rec (int li);
static void init_lpgpd_v (void);
/* 8/26/2011 */
void fillm_times_rec(int j, int i, int pi, int pj);
static void init_migration_counts (void);
//static void fillmrec(int j, int i, int pi, int pj);
//static void init_migration_counts_times (void);

static void init_mutation_scalar_rec (int li);
static void start_setup_L (char infilename[], int *fpstri, char fpstr[], int currentid);
static void add_priorinfo_to_output(char priorfilename[],int *fpstri, char fpstr[]);
static void start_setup_C (void);
static void finish_setup_C ();
static void set_tvalues (int ci);
static void setup_T ();
static void finish_setup_L (void);
static void reportparamcounts(int *fpstri, char fpstr[]);
/******** LOCAL FUNCTIONS ************/

/* used if modeloptions[ANCESTRALPOPSIZESSETTOLARGESTSAMPLED]==1  */
void addpopsizeidtext(int *fpstri, char fpstr[])
{
  int i;
  SP"\nAncestral Population Size Parameter Simplification\n");
  SP"--------------------------------------------------\n");
  SP"\nEach ancestral population size set to be the same as the size its largest descendant sampled population\n");
  SP"Sampled Populations and Polymorphism Summaries\n"); 
  SP"\tPop#\tpolysum\n");
  for (i=0;i<npops;i++)
    SP"\t%d\t%.3lf\n",i, varrank[i]);
  SP"\nAncestral Populations Connected to which Sampled Populations\n");
  SP"\tAncPop#\tSampPop#\n");
  for (i=npops;i<numtreepops;i++)
    SP"\t%d\t%d\n",i,varranktree[i].samppopnum);
  SP"\n");
} //addpopsizeidtext(int *fpstri, char fpstr[])

void fillvarranktree(int nowpop)
{
  if (nowpop < npops)
  {
    varranktree[nowpop].samppopnum = nowpop;
    varranktree[nowpop].var = varrank[nowpop];
  }
  else
  {
    fillvarranktree(C[0]->poptree[nowpop].up[0]);
    fillvarranktree(C[0]->poptree[nowpop].up[1]);
    if (varranktree[C[0]->poptree[nowpop].up[0]].var < varranktree[C[0]->poptree[nowpop].up[1]].var)
    {
      varranktree[nowpop] = varranktree[C[0]->poptree[nowpop].up[1]];
    }
    else
    {
      varranktree[nowpop] = varranktree[C[0]->poptree[nowpop].up[0]];
    }
  }

}// fillvarranktree,  recursive 

int
numvarHKY (int li, int b, int e)
{
  int i, j, tot = 0;
  for (i = 0; i < L[li].numsites; i++)

  {
    j = b + 1;
    while (j < e && L[li].seq[b][i] == L[li].seq[j][i])
      j++;
    if (j < L[li].numgenes)
      tot++;
  }
  return tot;
}

int
checkaresis (int period, int i, int j)  // returns 1 if populations i and j are both in period, and they are sister populations
{
  int aresis, k;
  aresis = 0;
  if (ISELEMENT (C[0]->plist[period][i], C[0]->periodset[period])
      && ISELEMENT (C[0]->plist[period][j], C[0]->periodset[period]))
  {
    for (k = period + 1; k < npops; k++)
    {
      if ((C[0]->droppops[k][0] == C[0]->plist[period][i]
           && C[0]->droppops[k][1] == C[0]->plist[period][j])
          || (C[0]->droppops[k][0] == C[0]->plist[period][j]
              && C[0]->droppops[k][1] == C[0]->plist[period][i]))
        aresis = 1;
    }
  }
  return aresis;
}                               //checkaresis 

void
setimigarray (int ii)
{
  int jj;

  if (modeloptions[EXPOMIGRATIONPRIOR])
  {
    if (!calcoptions[USEPRIORFILE])
      imig[ii].pr.mean = mprior;
    else
    {
      if (modeloptions[ONEMIGRATIONPARAMETER])
        imig[ii].pr.mean =  max_m_from_priorfile;
      else
        imig[ii].pr.mean = imig[ii].pr.max; //this was previously set in setup_iparams()
    }
    imig[ii].pr.max = EXPOMIGPLOTSCALE * imig[ii].pr.mean;
    imig[ii].pr.min = 0;
  }
  else
  {
    if (!calcoptions[USEPRIORFILE])
      imig[ii].pr.max = DMAX (mprior, MPRIORMIN); 
    else
    {
      if (modeloptions[ONEMIGRATIONPARAMETER])
        imig[ii].pr.max =  max_m_from_priorfile;
      //else this was previously set in setup_iparams()
    }
    imig[ii].pr.min = 0;
    imig[ii].pr.mean = (imig[ii].pr.max - imig[ii].pr.min) / 2.0;
  }
  imig[ii].xy = static_cast<plotpoint *> 
                (calloc (GRIDSIZE, sizeof (struct plotpoint)));
  for (jj = 0; jj < GRIDSIZE; jj++)
  {
    imig[ii].xy[jj].x =
      imig[ii].pr.min +
      ((jj + 0.5) * (imig[ii].pr.max - imig[ii].pr.min)) / GRIDSIZE;
  }
  imig[ii].wp.n = 0;

}                               // setimigarray

/* initialize the splitting rate, population size and migration instances of struct i_param. 
   This function is complicated with lots of sections because of many user options */
void
setup_iparams (void)
{
  int i, j, k, ci = 0, ii, jj, kk, pj, ni, mi, mcheck;
  char temps[PARAMSTRLEN];

  /* set up the population size parameters */
  if (modeloptions[PARAMETERSBYPERIOD])
  {
    numpopsizeparams = (npops * (npops + 1)) / 2;       // every population in every period gets a parameter
  }
  else if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
  {
    numpopsizeparams = npops;
  }
  else if (modeloptions[ANCESTRALPOPSIZESSETTOLARGESTSAMPLED] == 1)
  {
    numpopsizeparams = npops;
    fillvarranktree(numtreepops-1);
  }
  else
  {
    numpopsizeparams = numtreepops;     // every distinct population gets a parameter
  }

  itheta = static_cast<i_param *> 
                (malloc (numpopsizeparams * sizeof (struct i_param)));
  for (i = 0; i < numpopsizeparams; i++)
  {

    if (calcoptions[USEPRIORFILE])
      itheta[i].pr.max = popsizeprior_fromfile[i];
    else
      itheta[i].pr.max = thetaprior;
    itheta[i].pr.min = 0;
    itheta[i].xy = static_cast<plotpoint *> (calloc (GRIDSIZE, sizeof (struct plotpoint)));
    for (j = 0; j < GRIDSIZE; j++)
    {
      itheta[i].xy[j].x =
        itheta[i].pr.min +
        ((j + 0.5) * (itheta[i].pr.max - itheta[i].pr.min)) / GRIDSIZE;
    }
  }
  if (modeloptions[PARAMETERSBYPERIOD])
  {
    for (i = 0, ii = 0, k = npops; i <= lastperiodnumber; i++, k--)
    {
      for (j = 0; j < k; j++)
      {
        pj = C[0]->plist[i][j];
        sprintf (itheta[ii].str, "q%d,%d", i, pj);
        itheta[ii].b = C[0]->poptree[pj].b;
        itheta[ii].e = C[0]->poptree[pj].e;
        itheta[ii].wp.n = 1;
        itheta[ii].wp.p = static_cast<int *> 
                    (malloc (itheta[ii].wp.n * sizeof (int)));
        itheta[ii].wp.r = static_cast<int *> 
                    (malloc (itheta[ii].wp.n * sizeof (int)));
        *itheta[ii].wp.p = i;
        *itheta[ii].wp.r = j;
        ii++;
      }
    }
  }
  else 
  {
    for (ii = 0; ii < numpopsizeparams; ii++)
    {
      sprintf (itheta[ii].str, "q%d", ii);
      itheta[ii].b = C[0]->poptree[ii].b;
      itheta[ii].e = C[0]->poptree[ii].e;
      itheta[ii].wp.n = 0;
      if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
      {
        assert (modeloptions[ANCESTRALPOPSIZESSETTOLARGESTSAMPLED] == 0);
        itheta[ii].wp.n = 1;
        if (itheta[ii].wp.n > 0)
        {
          itheta[ii].wp.p = static_cast<int *> 
                    (malloc (itheta[ii].wp.n * sizeof (int)));
          itheta[ii].wp.r = static_cast<int *> 
                    (malloc (itheta[ii].wp.n * sizeof (int)));
        }
        ni = 0;
        for (i = 0; i < lastperiodnumber; i++)
          for (j = 0; j < npops; j++)
          {
            if (ii == C[0]->plist[i][j])
            {
              itheta[ii].wp.p[ni] = i;
              itheta[ii].wp.r[ni] = j;
              ni++;
            }
          }
        assert (ni == itheta[ii].wp.n);

      }
      else 
      {
        if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0
            && modeloptions[ANCESTRALPOPSIZESSETTOLARGESTSAMPLED] == 0)
        {
          for (i = 0, k = npops; i <= lastperiodnumber; i++, k--)
            for (j = 0; j < k; j++)
              itheta[ii].wp.n += (ii == C[0]->plist[i][j]);
          if (itheta[ii].wp.n > 0)
          {
            itheta[ii].wp.p = static_cast<int *> 
                        (malloc (itheta[ii].wp.n * sizeof (int)));
            itheta[ii].wp.r = static_cast<int *> 
                        (malloc (itheta[ii].wp.n * sizeof (int)));
          }
          ni = 0;
          for (i = 0, k = npops; i <= lastperiodnumber; i++, k--)
            for (j = 0; j < k; j++)
            {
              if (ii == C[0]->plist[i][j])
              {
                itheta[ii].wp.p[ni] = i;
                itheta[ii].wp.r[ni] = j;
                ni++;
              }
            }
          assert (ni == itheta[ii].wp.n);
        }
        else if (modeloptions[ANCESTRALPOPSIZESSETTOLARGESTSAMPLED] == 1)
        {
          /* need to know what sampled population is  */
          for (i = 0, k = npops; i <= lastperiodnumber; i++, k--)
            for (j = 0; j < k; j++)
              itheta[ii].wp.n += (ii == C[0]->plist[i][j] || ii==varranktree[C[0]->plist[i][j]].samppopnum);
          if (itheta[ii].wp.n > 0)
          {
            itheta[ii].wp.p = static_cast<int *> 
                            (malloc (itheta[ii].wp.n * sizeof (int)));
            itheta[ii].wp.r = static_cast<int *> 
                            (malloc (itheta[ii].wp.n * sizeof (int)));
          }
          ni = 0;
          for (i = 0, k = npops; i <= lastperiodnumber; i++, k--)
            for (j = 0; j < k; j++)
            {
              if (ii == C[0]->plist[i][j] || ii==varranktree[C[0]->plist[i][j]].samppopnum)
              {
                itheta[ii].wp.p[ni] = i;
                itheta[ii].wp.r[ni] = j;
                ni++;
              }
            }
          assert (ni == itheta[ii].wp.n);
        }
      }
    }
  }
  /* set up the migration parameters 
   Three, somewhat repetitive loops here:
    1. count how many migration parameters are needed. 
    2. build  most of the imig[] structures, including str, b, e; start to build wp part of imig[]
    3.  finish building wp parts of imig[]

     there is a complex series of checks to ensure that the intended model is being followed
     loop through periods
     check order:
     MIGRATIONBETWEENSAMPLED  - migration only between sampled populations
     PARAMETERSBYPERIOD  - every population size and migration parameter applies for only 1 period
     NOMIGBETWEENNONSISTERS  - set migration to zero between non-sister populations
     SINGLEMIGRATIONBOTHDIRECTIONS
     if ONEMIGRATIONPARAMETER applies  then there is only one migration rate,  whenever there is migration in the model
   */
   
  if (!modeloptions[PARAMETERSBYPERIOD])
    maxpossiblemigrateparams = 2*(npops-1)*(npops-1);
  /*count how many migration parameters are needed, use mi 
    set up a standard sequence of checks that get repeated when building imig[] */
  mi = 0;
  if (!modeloptions[NOMIGRATION] && npops > 1)
  {
    if (modeloptions[ONEMIGRATIONPARAMETER])
    {
      mi = 1;
      goto outsidefirstloop; // get out of this big condition section
    }
    for (k = 0; k < lastperiodnumber; k++)
      for (i = 0; i < npops - k - 1; i++)
        for (j = i + 1; j < npops - k; j++) if ( k == 0 ||
              ( modeloptions[MIGRATIONBETWEENSAMPLED]==0 &&
               (modeloptions[PARAMETERSBYPERIOD]==1 ||
                ( !ISELEMENT (C[ci]->plist[k][i], C[ci]->periodset[k - 1]) ||
                  !ISELEMENT (C[ci]->plist[k][j], C[ci]->periodset[k - 1])  
                  ) ) ) )
          /* tricky condition.  proceed to set mcheck if:
            where are considering only sampled populations (k==0)  or
            (not considering just sampled populations  AND
            EITHER we are considering each parameter to apply only to one period
                  OR one of the two populations being considered was not in the previous (k-1) period (meaning it is a new ancestral pop) */
          {
          /* if mcheck ends up as 1, we have one more parameter to add
             if it is also the case that migration is not the same in both directions,  then we have two parameters to add
             
             mcheck can get set to 0 if we require two populations to be sisters, and they are not
             or if the prior given for that pair of populations is zero 
             */
            if (modeloptions[NOMIGBETWEENNONSISTERS])
              mcheck = checkaresis (k, i, j);
            else
            {
              if (calcoptions[USEPRIORFILE]==0 || (calcoptions[USEPRIORFILE] && mprior_fromfile[C[ci]->plist[k][i]][C[ci]->plist[k][j]] > MINPARAMVAL))
                mcheck = 1;
              else
                mcheck = 0;
            }
            if (mcheck)
            {
              mi++;
              mi += modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS]==0;
            }
          }
outsidefirstloop: 
    assert (!assignmentoptions[POPULATIONASSIGNMENTINFINITE] || mi == npops * (npops - 1));

    /* now do another loop, much like the previous one, 
      but now start to build the imig[] structures, including str, b, e
      start to build the wp part of imig[] and determine wp.n*/
    nummigrateparams = mi;
    imig = static_cast<i_param *> 
                (malloc (nummigrateparams * sizeof (struct i_param)));
    if (modeloptions[ONEMIGRATIONPARAMETER])
    {

      mi = 0;
      setimigarray (mi);
      imig[mi].wp.n = 0;
      imig[mi].b = 0;
      imig[mi].e = lastperiodnumber - 1;
      sprintf (imig[mi].str, "m_all");
    }
    else
      mi = -1;
    for (k = 0; k < lastperiodnumber; k++)
      for (i = 0; i < npops - k - 1; i++)
        for (j = i + 1; j < npops - k; j++) if ((k == 0) ||
              (!modeloptions[MIGRATIONBETWEENSAMPLED] && (modeloptions[PARAMETERSBYPERIOD] ||
                ((!ISELEMENT(C[ci]->plist[k][i], C[ci]->periodset[k - 1]))
                 ||  (!ISELEMENT(C[ci]->plist[k][j], C[ci]->periodset[k - 1]))))))
          {
            if (modeloptions[NOMIGBETWEENNONSISTERS])
              mcheck = checkaresis (k, i, j);
            else
            {
               if (!calcoptions[USEPRIORFILE] || (calcoptions[USEPRIORFILE] && mprior_fromfile[C[ci]->plist[k][i]][C[ci]->plist[k][j]] > MINPARAMVAL))
                mcheck = 1;
               else
                 mcheck = 0;
            }
            if (mcheck)
            {
              if (!modeloptions[ONEMIGRATIONPARAMETER])
              {
                mi++;
                if (calcoptions[USEPRIORFILE])
                {
                    imig[mi].pr.max = mprior_fromfile[C[ci]->plist[k][i]][C[ci]->plist[k][j]];
                    if (imig[mi].pr.max <= MINPARAMVAL)
                      IM_err(IMERR_PRIORFILEVALS,"migration rate set too low from %d to %d",i,j);

                }
                setimigarray (mi);
                imig[mi].b = k;
                imig[mi].e = k;
                imig[mi].wp.n = 1;
              }
              else
                imig[mi].wp.n++;
              imig[mi].wp.n += modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS]; // can't also have ONEMIGRATIONPARAMETER
              if (modeloptions[PARAMETERSBYPERIOD] && !modeloptions[ONEMIGRATIONPARAMETER])
              {
                sprintf (imig[mi].str, "m%d,", k);
                if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])
                  sprintf (temps, "%d<>%d", C[ci]->plist[k][i], C[ci]->plist[k][j]);
                else
                  sprintf (temps, "%d>%d", C[ci]->plist[k][i], C[ci]->plist[k][j]);
                strcat (imig[mi].str, temps);
              }
              else
              {
                kk = k + 1;
                while (kk < lastperiodnumber
                       && ISELEMENT (C[ci]->plist[k][i],
                                     C[ci]->periodset[kk])
                       && ISELEMENT (C[ci]->plist[k][j],
                                     C[ci]->periodset[kk]))
                {
                  if (!modeloptions[ONEMIGRATIONPARAMETER])
                    imig[mi].e = kk;
                  kk++;
                  imig[mi].wp.n++;
                  imig[mi].wp.n += modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS]; // can't also have ONEMIGRATIONPARAMETER
                }
                if (!modeloptions[ONEMIGRATIONPARAMETER])
                {
                  if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])  // can't also have ONEMIGRATIONPARAMETER
                    sprintf (imig[mi].str, "m%d<>%d", C[ci]->plist[k][i], C[ci]->plist[k][j]);
                  else
                    sprintf (imig[mi].str, "m%d>%d", C[ci]->plist[k][i], C[ci]->plist[k][j]);
                }
              }
              if (!modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])
              {
                if (!modeloptions[ONEMIGRATIONPARAMETER])
                {
                  mi++;
                  if (calcoptions[USEPRIORFILE])
                  {
                    imig[mi].pr.max = mprior_fromfile[C[ci]->plist[k][j]][C[ci]->plist[k][i]];
                    if (imig[mi].pr.max <= MINPARAMVAL)
                      IM_err(IMERR_PRIORFILEVALS,"migration rate set too low from %d to %d",i,j);
                  }
                  setimigarray (mi);
                  imig[mi].b = k;
                  imig[mi].e = k;
                  imig[mi].wp.n = 1;
                }
                else
                  imig[mi].wp.n++;
                if (modeloptions[PARAMETERSBYPERIOD] && !modeloptions[ONEMIGRATIONPARAMETER])
                {
                  sprintf (imig[mi].str, "m%d,%d>%d", k, C[ci]->plist[k][j], C[ci]->plist[k][i]);
                }
                else
                {
                  kk = k + 1;
                  while (kk < lastperiodnumber &&
                         ISELEMENT (C[ci]->plist[k][i],
                                    C[ci]->periodset[kk])
                         && ISELEMENT (C[ci]->plist[k][j],
                                       C[ci]->periodset[kk]))
                  {
                    if (!modeloptions[ONEMIGRATIONPARAMETER])
                      imig[mi].e = kk;
                    kk++;
                    imig[mi].wp.n++;
                  }
                  if (!modeloptions[ONEMIGRATIONPARAMETER])
                    sprintf (imig[mi].str, "m%d>%d", C[ci]->plist[k][j], C[ci]->plist[k][i]);
                }
              }
            }
          }
    
    /* initialize the parts of wp */
    for (mi = 0; mi < nummigrateparams; mi++)
    {
      if (imig[mi].wp.n > 0)
      {
        imig[mi].wp.p = static_cast<int *> 
                    (malloc (imig[mi].wp.n * sizeof (int)));
        imig[mi].wp.r = static_cast<int *> 
                    (malloc (imig[mi].wp.n * sizeof (int)));
        imig[mi].wp.c = static_cast<int *> 
                    (malloc (imig[mi].wp.n * sizeof (int)));
      }
    }
    /* now do another loop, much like the previous two. this time finish building the wp parts of imig[]
      using wp.n that was determined in the previous loop */
    if (modeloptions[ONEMIGRATIONPARAMETER])
    {
      mi = 0;
      ni=-1;
    }
    else
      mi = -1;
    for (k = 0; k < lastperiodnumber; k++)
      for (i = 0; i < npops - k - 1; i++)
        for (j = i + 1; j < npops - k; j++) if ((k == 0) || 
              (!modeloptions[MIGRATIONBETWEENSAMPLED] && (modeloptions[PARAMETERSBYPERIOD] ||
                ((!ISELEMENT (C[ci]->plist[k][i], C[ci]->periodset[k - 1]))
                 || (!ISELEMENT(C[ci]->plist[k][j], C[ci]->periodset[k - 1]))))))
          {
            if (modeloptions[NOMIGBETWEENNONSISTERS])
              mcheck = checkaresis (k, i, j);
            else
            {
               if (!calcoptions[USEPRIORFILE] || (calcoptions[USEPRIORFILE] && mprior_fromfile[C[ci]->plist[k][i]][C[ci]->plist[k][j]] > MINPARAMVAL))
                mcheck = 1;
               else
                 mcheck = 0;
            }
            if (mcheck)
            {
               if (!modeloptions[ONEMIGRATIONPARAMETER])
               {
                mi++;
                ni = 0;
               }
               else
                 ni++;
              imig[mi].wp.p[ni] = k;
              imig[mi].wp.r[ni] = i;
              imig[mi].wp.c[ni] = j;
              assert (ni < imig[mi].wp.n);
              if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])   // can't also have ONEMIGRATIONPARAMETER
              {
                ni++;
                imig[mi].wp.p[ni] = k;
                imig[mi].wp.r[ni] = j;
                imig[mi].wp.c[ni] = i;
                assert (ni < imig[mi].wp.n);
              }
              if (!modeloptions[PARAMETERSBYPERIOD])
              {
                kk = k + 1;
                while (kk < lastperiodnumber &&
                       ISELEMENT (C[ci]->plist[k][i],
                                  C[ci]->periodset[kk])
                       && ISELEMENT (C[ci]->plist[k][j],
                                     C[ci]->periodset[kk]))
                {
                  ii = 0;
                  while (C[0]->plist[kk][ii] != C[0]->plist[k][i])
                    ii++;
                  assert (ii < npops - kk);
                  jj = 0;
                  while (C[0]->plist[kk][jj] != C[0]->plist[k][j])
                    jj++;
                  assert (jj < npops - kk);
                  ni++;
                  imig[mi].wp.p[ni] = kk;
                  imig[mi].wp.r[ni] = ii;
                  imig[mi].wp.c[ni] = jj;
                  assert (ni < imig[mi].wp.n);
                  if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])   // can't also have ONEMIGRATIONPARAMETER
                  {
                    ni++;
                    imig[mi].wp.p[ni] = kk;
                    imig[mi].wp.r[ni] = jj;
                    imig[mi].wp.c[ni] = ii;
                    assert (ni < imig[mi].wp.n);
                  }
                  kk++;
                }
              }
              if (!modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS]) 
              {
                if (!modeloptions[ONEMIGRATIONPARAMETER])
                {
                  mi++;
                  ni = 0;
                }
                else
                  ni++;
                imig[mi].wp.p[ni] = k;
                imig[mi].wp.r[ni] = j;
                imig[mi].wp.c[ni] = i;
                assert (ni < imig[mi].wp.n);
                if (!modeloptions[PARAMETERSBYPERIOD])
                {
                  kk = k + 1;
                  while (kk < lastperiodnumber &&
                         ISELEMENT (C[ci]->plist[k][i],
                                    C[ci]->periodset[kk])
                         && ISELEMENT (C[ci]->plist[k][j],
                                       C[ci]->periodset[kk]))
                  {
                    ii = 0;
                    while (C[0]->plist[kk][ii] != C[0]->plist[k][i])
                      ii++;
                    assert (ii < npops - kk);
                    jj = 0;
                    while (C[0]->plist[kk][jj] != C[0]->plist[k][j])
                      jj++;
                    assert (jj < npops - kk);
                    ni++;
                    imig[mi].wp.p[ni] = kk;
                    imig[mi].wp.r[ni] = jj;
                    imig[mi].wp.c[ni] = ii;
                    kk++;
                    assert (ni < imig[mi].wp.n);
                  }
                }
              }
            }
          }                           // initialize wp
  }

/*
	type     |  # values  |  cumulative total at end
    cc	       numpopsizeparams   numpopsizeparams
	fc	       numpopsizeparams   2*numpopsizeparams
	hcc	       numpopsizeparams   3*numpopsizeparams
	mc         nummigrateparams   3*numpopsizeparams + nummigrateparams
	fm         nummigrateparams   3*numpopsizeparams + 2*nummigrateparams
	qintegrate numpopsizeparams   4*numpopsizeparams + 2*nummigrateparams
	mintegrate nummigrateparams   4*numpopsizeparams + 3*nummigrateparams
	pdg             1        4*numpopsizeparams + 3*nummigrateparams +  1 
	probg           1        4*numpopsizeparams + 3*nummigrateparams +  2
	t          numsplittimes 4*numpopsizeparams + 3*nummigrateparams +  numsplittimes + 2
*/

  /* set values of position markers for gsampinf */
  
  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
  {
    assert (numpopsizeparams == npops);
    if (modeloptions[NOMIGRATION] == 0)
    {
      assert (nummigrateparams == npops * (npops - 1));
    }
  }

  gsamp_ccp = 0;
  gsamp_fcp = gsamp_ccp + numpopsizeparams;
  gsamp_hccp = gsamp_fcp + numpopsizeparams;
  gsamp_mcp = gsamp_hccp + numpopsizeparams;
  gsamp_fmp = gsamp_mcp + nummigrateparams;
  gsamp_qip = gsamp_fmp + nummigrateparams;
  gsamp_mip = gsamp_qip + numpopsizeparams;
  gsamp_pdgp = gsamp_mip + nummigrateparams;
  gsamp_probgp =  gsamp_pdgp + 1;
  gsamp_tp = gsamp_probgp + 1;

  /*// to redo the 48 locus IMa runs of 2006
itheta[0].pr.max = 3.82;
itheta[1].pr.max = 0.87;
itheta[2].pr.max = 1.53;
imig[0].pr.max = 6.0;
imig[1].pr.max = 3.0; */
}                               /* setup_iparams */


/* nomigrationchecklist is three arrays (p for period, r for row, and c for column), each of a set length (nomigrationchecklist.n )

nomigrationchecklist is used to check gweight->mc[][][]
gweight->mc[k][i][j] is the number of migration events in the genealogies in period k 
from population C[ci]->plist[k][i] to population C[ci]->plist[k][j]

So the ith element of the arrays in nomigrationchecklist is used to check the value of 
gweight->mc[nomigrationchecklist.p[i]][nomigrationchecklist.r[i]][nomigrationchecklist.c[i]]
(this is done in update_gree_common.c).  If that value if not zero the update gets rejected. 

building nomigrationchecklist requires identifying which pairs of populations, 
C[ci]->plist[k][i] and C[ci]->plist[k][j], in each period k,  should not be
exchanging migrants.

These pairs are those that meet any of the following:
1. modeloptions[MIGRATIONBETWEENSAMPLED]==1 and the populations are not both sampled populations
2. modeloptions[NOMIGBETWEENNONSISTERS]==1  and the populations are not sisters
3. they have had their priors set to zero in a priorfile

These correspond to the following:
1. modeloptions[MIGRATIONBETWEENSAMPLED]==1 && C[ci]->plist[k][i] >= npops && C[ci]->plist[k][j] >= npops
2. modeloptions[NOMIGBETWEENNONSISTERS]==1 && checkaresis (k, i, j)==1
3. calcoptions[USEPRIORFILE]==1 && mprior_fromfile[C[ci]->plist[k][i]][C[ci]->plist[k][j]] <= MINPARAMVAL


To set this up we must go through two loops,  first to cound how many (i.e. set nomigrationchecklist.n)
Then after initializing the arrays in nomigrationchecklist, we must do another similar loop to fill them up. 
*/

/* this should not be needed if modeloptions[NOMIGRATION]==1 */
void
set_nomigrationchecklist ()
{
  int n, i, j, k, ci = 0;

  nomigrationchecklist.n = 0;
  nomigrationchecklist.p = NULL;
  nomigrationchecklist.r = NULL;
  nomigrationchecklist.c = NULL;
  for (k = 0; k < lastperiodnumber; k++)
    for (i = 0; i < npops - k - 1; i++)
      for (j = i + 1; j < npops - k; j++)
      {
        if ((modeloptions[MIGRATIONBETWEENSAMPLED]==1 && (C[ci]->plist[k][i] >= npops || C[ci]->plist[k][j] >= npops)) ||
            (modeloptions[NOMIGBETWEENNONSISTERS]==1 && checkaresis (k, i, j)==0) ||
            (calcoptions[USEPRIORFILE]==1 && mprior_fromfile[C[ci]->plist[k][i]][C[ci]->plist[k][j]] <= MINPARAMVAL))
        {
          nomigrationchecklist.n+=2;
        }
      }
    if (nomigrationchecklist.n > 0)
    {
      nomigrationchecklist.p = static_cast<int *> 
                        (malloc (nomigrationchecklist.n * sizeof (int)));
      nomigrationchecklist.r = static_cast<int *> 
                        (malloc (nomigrationchecklist.n * sizeof (int)));
      nomigrationchecklist.c = static_cast<int *> 
                        (malloc (nomigrationchecklist.n * sizeof (int)));
    }
    for (n = -1, k = 0; k < lastperiodnumber; k++)
      for (i = 0; i < npops - k - 1; i++)
        for (j = i + 1; j < npops - k; j++)
        {
          if ((modeloptions[MIGRATIONBETWEENSAMPLED]==1 && (C[ci]->plist[k][i] >= npops || C[ci]->plist[k][j] >= npops)) ||
            (modeloptions[NOMIGBETWEENNONSISTERS]==1 && checkaresis (k, i, j)==0) ||
            (calcoptions[USEPRIORFILE]==1 && mprior_fromfile[C[ci]->plist[k][i]][C[ci]->plist[k][j]] <= MINPARAMVAL))
          {
              n++;
              nomigrationchecklist.p[n] = k;
              nomigrationchecklist.r[n] = i;
              nomigrationchecklist.c[n] = j;
              n++;
              nomigrationchecklist.p[n] = k;
              nomigrationchecklist.r[n] = j;
              nomigrationchecklist.c[n] = i;
          }
        }
  assert(nomigrationchecklist.n == n+1);
  return;
}                               // set_nomigrationchecklist

#define MAXPRIORSCALETRY 10000000       // max number of times to try getting starting mutation rates compatible with priors
void
setuinfo (double summut)
{
  int i, li, ui, uj, k;
  int numuprior = 0;
  double priorscaletry = 0;
  double U, r;
  int doneu, rcheck, numupair, *upriorpairlist1, *upriorpairlist2;
  double prodcheck, maxr, newr, d;

/*  ul contains a list of locations by locus and within locus of all mutation rate scalars
    uii is the position in that list
	in set_mcparam_values() the priors on mutation rate scalars are set to standard pos (max) and neg (min) values (they are on a log scale)
	   and are stored in the mcinf.pr   (e.g. C[0]->G[0].u[0].mcinf.pr)
	if the input file has prior ranges on mutation yets,  these are read into uperyear.pr  (e.g. C[ci]->G[li].u[ui].uperyear.pr)
	The ratios of these uperyear values among loci must be taken to reset the values of mcinf.pr   
*/
  for (i = 0, li = 0; li < nloci; li++)
    for (ui = 0; ui < L[li].nlinked; ui++)
    {
      ul[i].l = li;
      ul[i].u = ui;
      L[li].uii[ui] = i;
      i++;
    }

  /* reset the mutation rate values to have geometric mean of 1 */
  for (prodcheck = 1, ui = 0; ui < nurates; ui++)
  {
    C[0]->G[ul[ui].l].uvals[ul[ui].u] =
      exp (uvals[ul[ui].l][ul[ui].u] - summut / nurates);
    prodcheck *= C[0]->G[ul[ui].l].uvals[ul[ui].u];
  }
  if (fabs (log (prodcheck)) > 1e-5)
    IM_err(IMERR_MUTSCALEPRODUCTFAIL,"product of mutation scalars not close to or equal to 1: %lf",prodcheck); 
  for (ui = 0; ui < nurates; ui++)
  {
    if (L[ul[ui].l].uperyear_prior[ul[ui].u].min != 0)
      numuprior++;
  }

  if (calcoptions[MUTATIONPRIORRANGE])
  {
    assert (numuprior > 1);
    numupair = numuprior * (numuprior - 1) / 2;
    upriorpairlist1 = static_cast<int *> 
                    (calloc ((size_t) numupair, sizeof (int)));
    upriorpairlist2 = static_cast<int *> 
                    (calloc ((size_t) numupair, sizeof (int)));
    k = 0;
    for (ui = 0; ui < nurates - 1; ui++)
    {
      for (uj = ui + 1; uj < nurates; uj++)
      {
        if (L[ul[ui].l].uperyear_prior[ul[ui].u].min != 0
            && L[ul[uj].l].uperyear_prior[ul[uj].u].min != 0)

        {
          upriorpairlist1[k] = ui;
          upriorpairlist2[k] = uj;
          k++;
        }
      }
    }

    /* urri[i][j] has a 0 if neither i nor j has a prior. 1 if i has a prior and j does not,  -1 if i does not have a pior and j does  */
    for (ui = 0; ui < nurates; ui++)
    {
      for (uj = 0; uj < nurates; uj++)
      {
        urri[ui][uj] = 0;
        if (ui != uj
            && L[ul[ui].l].uperyear_prior[ul[ui].u].min != 0
            && L[ul[uj].l].uperyear_prior[ul[uj].u].min != 0)
        {
          urrlow[ui][uj] =
            log (L[ul[ui].l].uperyear_prior[ul[ui].u].min /
                 L[ul[uj].l].uperyear_prior[ul[uj].u].max);
          urrhi[ui][uj] =
            log (L[ul[ui].l].uperyear_prior[ul[ui].u].max /
                 L[ul[uj].l].uperyear_prior[ul[uj].u].min);
          urri[ui][uj] = urri[uj][ui] = 2;
        }
        else
        {
          if (ui != uj
              && L[ul[ui].l].uperyear_prior[ul[ui].u].min != 0
              && L[ul[uj].l].uperyear_prior[ul[uj].u].min == 0)
            urri[ui][uj] = 1;
          if (ui != uj
              && L[ul[uj].l].uperyear_prior[ul[uj].u].min != 0
              && L[ul[ui].l].uperyear_prior[ul[ui].u].min == 0)
            urri[ui][uj] = -1;
        }
      }
    }
    /* need to set all of the uscalers so that their ratios are in the ranges defined by the priors */
    /* it is possible that suiteable sets of scalars will not be able to be found */
    maxr = 3 * L[0].u_rec[0].pr.max;

    do
    {
      doneu = 1;
      for (i = 0; i < numupair; i++)

      {
        ui = upriorpairlist1[i];
        uj = upriorpairlist2[i];
        r =
          log (C[0]->G[ul[ui].l].uvals[ul[ui].u] /
               C[0]->G[ul[uj].l].uvals[ul[uj].u]);
        rcheck = (r >= urrlow[ui][uj] && r <= urrhi[ui][uj]);
        doneu = doneu && rcheck;
        if (!rcheck)
        {
          do
          {
            U = uniform ();
            if (U > 0.5)
              newr = r + maxr * (2.0 * U - 1.0);

            else
              newr = r - maxr * U * 2.0;
            if (newr > maxr)
              newr = 2.0 * maxr - newr;

            else if (newr < -maxr)
              newr = 2.0 * (-maxr) - newr;
            d = exp ((newr - r) / 2);
          } while ((newr <= urrlow[ui][uj] || newr >= urrhi[ui][uj]));
          C[0]->G[ul[ui].l].uvals[ul[ui].u] *= d;
          C[0]->G[ul[uj].l].uvals[ul[uj].u] /= d;
        }
      }
      priorscaletry++;
    } while (!doneu && priorscaletry < MAXPRIORSCALETRY);

    if (priorscaletry >= MAXPRIORSCALETRY)
      IM_err(IMERR_MUTSCALARPRIORRANGEFAIL,"More than %d failed attempts at finding a valid set of mutation scalars, prior ranges probably too restrictive",MAXPRIORSCALETRY);
    for (prodcheck = 1, ui = 0; ui < nurates; ui++)
      prodcheck *= C[0]->G[ul[ui].l].uvals[ul[ui].u];
      //prodcheck *= C[0]->G[ul[ui].l].u[ul[ui].u].uperyear.val;

    if (fabs (log (prodcheck)) > 1e-5)
    {
      IM_err(IMERR_MUTSCALEPRODUCTFAIL,"product of mutation scalars not close to or equal to 1: %lf",prodcheck); 
    }
    XFREE (upriorpairlist1);
    XFREE (upriorpairlist2);
  }
}                               /* setuinfo */

double
set_uvals (void)
  /* set mutation rate scalars */
  /* get a relative mutation rate for a locus by summing up estimates of 4Nu for all the populations
     the geometric mean of these relative values are then used to set the starting mutation rates */
{
  int j, k, li, ui;
  int b, e;
  double w, w2;
  double temp, setutemp, sumq;
  int npopstemp;
  

  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    npopstemp = npops + 1;
  }
  else
  {
    npopstemp = npops;
  }


  if (assignmentoptions[POPULATIONASSIGNMENT] == 0)
  {
    /*************************************************************************/
  setutemp = 0;
  for (li = 0; li < nloci; li++)
  {
    if (L[li].model == INFINITESITES
        || L[li].model == JOINT_IS_SW 
        || L[li].model == HKY)
    {
      sumq = 0;
      b = 0;
      e = L[li].samppop[0] - 1;
      for (k = 0; k < npopstemp; k++)
      {
        if (e >= b)
        {
          for (j = 1, w = 0.0; j <= e; j++)
            w += 1 / (double) j;
          if (w == 0)
            w = 1;
          if (L[li].model == HKY)
            temp = numvarHKY (li, b, e) / w;
          else
            temp = numsitesIS[li][k] / w;
        }
        else
        {
          temp = 0;
        }
        if (modeloptions[ANCESTRALPOPSIZESSETTOLARGESTSAMPLED] && k<npops)
        {
          varrank[k] += temp;
        }
        b = e + 1;
        if (k + 1 == npopstemp)
        {
          e = e + L[li].numgenesunknown;
        }
        else
        {
          e = e + L[li].samppop[k + 1];
        }
        sumq += temp;

      }
      sumq = DMAX (0.1, sumq);  // nominal low values in case of zero variation 
      uvals[li][0] = log (sumq);
      setutemp += uvals[li][0];
    }

    if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
    {
      somestepwise = 1;
      if (L[li].model == JOINT_IS_SW)
        ui = 1;
      else
        ui = 0;

      for (; ui < L[li].nlinked; ui++)
      {
        b = 0;
        e = L[li].samppop[0] - 1;
        sumq = 0;
        for (k = 0; k < npopstemp; k++)
        {
          if (e >= b)
          {
            for (j = b, w = 0, w2 = 0; j <= e; j++)
            {
              w += L[li].A[ui][j];
              w2 += SQR ((double) L[li].A[ui][j]);
            }
            if (k == npops)
            {
              temp =
                (double) 2 *(w2 -
                             SQR (w) / L[li].numgenesunknown) /
                (L[li].numgenesunknown);
            }
            else
            {
              temp =
                (double) 2 *(w2 -
                             SQR (w) / L[li].samppop[k]) / (L[li].samppop[k]);
            }
            sumq += temp;
          }
          else
          {
            temp = 0;
          }
          if (modeloptions[ANCESTRALPOPSIZESSETTOLARGESTSAMPLED] && k<npops)
          {
            varrank[k] += temp;
          }
          b = e + 1;
          if (k + 1 == npopstemp)
          {
            e = e + L[li].numgenesunknown;
          }
          else
          {
            e = e + L[li].samppop[k + 1];
          }
        }
        sumq = DMAX (0.1, sumq);        // nominal low values in case of zero variation 
        uvals[li][ui] = log (sumq);
        setutemp += uvals[li][ui];
      }
    }
  }

  // geomeanvar is used to set the prior on thetas, in setup_iparams(), depending on command line options 
  /* sum is the sum of the logs, across loci, of the maximal amount of variation found among the sampled populations */
  /* the scalar of thetaprior is multiplied times the geometric mean of the maximal estimates of variation across loci */
  geomeanvar = exp (setutemp / nurates);

    /*************************************************************************/
  }
  /* SANGCHUL: Sat May 16 16:48:10 EDT 2009
   * For assignment option, starting u values are all 1.
   */
  else if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    for (li = 0; li < nloci; li++)
    {
      if (L[li].model == INFINITESITES
          || L[li].model == JOINT_IS_SW 
          || L[li].model == HKY)
      {
        uvals[li][0] = 0.0;
      }

      if (L[li].model == STEPWISE 
          || L[li].model == JOINT_IS_SW)
      {
        somestepwise = 1;
        if (L[li].model == JOINT_IS_SW)
        {
          ui = 1;
        }
        else
        {
          ui = 0;
        }
        for (; ui < L[li].nlinked; ui++)
        {
          uvals[li][ui] = 0.0;
        }
      }
    }
    geomeanvar = 1.0;
    setutemp = 0.0;
  }
  return setutemp;
}                               /*set_uvals */

void
set_x (struct value_record *v, int isint)
{
  int i;
  for (i = 0; i < GRIDSIZE; i++)
  {
    if (v->do_logplot)
    {
      v->xy[i].x = exp (v->plotrange.min + (i + 0.5) * v->plotrescale); // the same scale is used for all mutation rate scalars 
    }
    else
    {
      if (isint)
        v->xy[i].x = (double) i;
      else
        v->xy[i].x = v->plotrange.min + ((i + 0.5) * (v->plotrange.max * v->plotrescale - v->plotrange.min)) / GRIDSIZE;
    }
  }
}                               // set_x 

void
init_value_record (struct value_record *v, int isint)
{
  if (v->do_xyplot)
  {
    v->xy = static_cast<plotpoint *> 
                (calloc (GRIDSIZE, sizeof (struct plotpoint)));
    set_x (v, isint);
  }

  if (v->do_trend)
    v->trend = static_cast<double *> (calloc (TRENDDIM, sizeof (double)));

  v->beforemin = v->aftermax = 0;
}                               // init_value_Record

void
init_g_rec (int li)
{
  int num_g_update_types = IM_UPDATE_GENEALOGY_NUMBER;
  L[li].g_rec = static_cast<chainstate_record_updates_and_values *> 
                (malloc (sizeof (struct chainstate_record_updates_and_values)));
  sprintf (L[li].g_rec->str, "gtree_%d", li);
  L[li].g_rec->num_uptypes = num_g_update_types;
  L[li].g_rec->upnames = static_cast<strnl *> 
                (malloc (L[li].g_rec->num_uptypes * sizeof (strnl)));
  sprintf (L[li].g_rec->upnames[IM_UPDATE_GENEALOGY_ANY],      "branch     ");
  sprintf (L[li].g_rec->upnames[IM_UPDATE_GENEALOGY_TOPOLOGY], "topology   ");
  sprintf (L[li].g_rec->upnames[IM_UPDATE_GENEALOGY_TMRCA],    "tmrca      ");
  L[li].g_rec->upinf = static_cast<update_rate_calc *> 
        (calloc ((size_t) L[li].g_rec->num_uptypes, 
        sizeof (struct update_rate_calc)));
  L[li].g_rec->num_vals = 1;
  L[li].g_rec->v = static_cast<value_record *> 
            (malloc (L[li].g_rec->num_vals * sizeof (struct value_record)));
  sprintf (L[li].g_rec->v->str, "g%d_tmrca", li);
  sprintf (L[li].g_rec->v->strshort, "tmrca%d", li);
  L[li].g_rec->v->plotrange.min = 0;
  L[li].g_rec->v->do_autoc = 1;
  L[li].g_rec->v->do_xyplot = outputoptions[PRINTTMRCA];
  L[li].g_rec->v->do_trend = 1;
  L[li].g_rec->v->plotrescale = 1.0;
  L[li].g_rec->v->do_logplot = 0;
  

}                               // init_g_rec

void
init_a_rec (int li)
{
  int num_a_update_types = IM_UPDATE_ASSIGNMENT_NUMBER;

  L[li].a_rec = static_cast<chainstate_record_updates_and_values *> 
            (malloc (sizeof (struct chainstate_record_updates_and_values)));
  sprintf (L[li].a_rec->str, "asn_%d", li);
  L[li].a_rec->num_uptypes = num_a_update_types;
  L[li].a_rec->upnames = static_cast<strnl *> 
                    (malloc (num_a_update_types * sizeof (strnl)));

  sprintf (L[li].a_rec->upnames[IM_UPDATE_ASSIGNMENT_RELABEL], "relabel");
  sprintf (L[li].a_rec->upnames[IM_UPDATE_ASSIGNMENT_BF], "bf");

  L[li].a_rec->upinf = static_cast<update_rate_calc *> 
                (calloc ((size_t) L[li].a_rec->num_uptypes, 
                          sizeof (struct update_rate_calc)));
  L[li].a_rec->num_vals = 1;
  L[li].a_rec->v = static_cast<value_record *> 
            (malloc (L[li].a_rec->num_vals * sizeof (struct value_record)));
  sprintf (L[li].a_rec->v->str, "a%2d_asn", li);
  strncpy (L[li].a_rec->v->strshort, L[li].a_rec->v->str, PARAMSTRLENSHORT - 1);
  L[li].a_rec->v->plotrange.min = 0;
  L[li].a_rec->v->plotrange.max = 10;
  L[li].a_rec->v->do_autoc = 1;
  L[li].a_rec->v->do_xyplot = 0;
  L[li].a_rec->v->do_trend = 1;
  L[li].a_rec->v->plotrescale = 1.0;
  L[li].a_rec->v->do_logplot = 0;
  init_value_record (L[li].a_rec->v, 0);
  return;
}                               /* init_a_rec */

void
init_a_rec_multichain ()
{
  int ci;
  Cupinf = static_cast<im_chainstate_updateonly *> 
            (malloc (numchains * sizeof (im_chainstate_updateonly)));

  for (ci = 0; ci < numchains; ci++)
  {
    sprintf (Cupinf[ci].str, "ac_%d", ci);
    Cupinf[ci].num_uptypes = IM_UPDATE_ASSIGNMENT_NUMBER;
    Cupinf[ci].upnames = static_cast<strnl *> 
                (malloc (IM_UPDATE_ASSIGNMENT_NUMBER * sizeof (strnl)));
    Cupinf[ci].upinf = static_cast<update_rate_calc *> 
                (calloc (IM_UPDATE_ASSIGNMENT_NUMBER, 
                         sizeof (struct update_rate_calc)));
    sprintf (Cupinf[ci].upnames[IM_UPDATE_ASSIGNMENT_RELABEL], "relabel");
    sprintf (Cupinf[ci].upnames[IM_UPDATE_ASSIGNMENT_BF], "bf");
  }
  return;
}


void
init_lpgpd_v (void)
{
  lpgpd_v = static_cast<value_record *> (malloc (sizeof (struct value_record)));
  sprintf (lpgpd_v->str, "Log[P(G)+P(D|G)]");
  sprintf (lpgpd_v->strshort, "Log[P]");
  lpgpd_v->plotrange.min = lpgpd_v->plotrange.max = 0;
  lpgpd_v->do_xyplot = 0;
  lpgpd_v->do_logplot = 0;
  lpgpd_v->do_trend = 1;
  lpgpd_v->do_autoc = 1;
  lpgpd_v->plotrescale = 1.0;
  init_value_record (lpgpd_v, 0);
  return;
}                               // init_lpgpd_rec
/* 8/26/2011 */
void fillm_times_rec(int j, int i, int pi, int pj)
{
  char stemp[6];
  sprintf(stemp,"m%d>%d",pi,pj);
  strcpy (migration_counts[j][i].str,stemp);
  migration_counts[j][i].plotrange.max = (double) GRIDSIZE - 1;       // used for counts so the grid position is the count #
  migration_counts[j][i].plotrange.min = 0;
}

/*void fillmrec(int j, int i, int pi, int pj)
{
  char stemp[6];
  sprintf(stemp,"%d>%d",pi,pj);
  strcpy (migration_counts_times[j][i].str, stemp);
  strcat (migration_counts_times[j][i].str, "_t");
  migration_counts_times[j][i].plotrange = T[numsplittimes - 1].pr; // same range as for splitting times
  strcpy (migration_counts_times[j][i+1].str,stemp);
  strcat (migration_counts_times[j][i+1].str, "_#");
  migration_counts_times[j][i+1].plotrange.max = (double) GRIDSIZE - 1;       // used for counts so the grid position is the count #
  migration_counts_times[j][i+1].plotrange.min = 0;
}*/

void
init_migration_counts (void)
{
  int i, j, k, pi, pj, ji, jj, mrows;//, nummigdirs;
  mrows = nloci + (nloci > 1); 
  nummigdirs = 2*(npops-1)*(npops-1);
  /*if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])
    nummigdirs = 2*nummigrateparams;  
  else
    nummigdirs = nummigrateparams;   */ 
  migration_counts = static_cast<value_record **> (malloc (mrows * sizeof (struct value_record *)));
  for (j = 0; j < mrows; j++)
    migration_counts[j] = static_cast<value_record *>
      (malloc (nummigdirs * sizeof (struct value_record)));
  for (j = 0; j < mrows; j++)
  {
    i=0;
    for (k= 0; k < lastperiodnumber; k++)
    {
      for (ji=0; ji< npops-k;ji++) 
      {
          pi = C[0]->plist[k][ji];
          if (k==0)
          {
            for (jj=ji+1; jj< npops-k;jj++) 
            {
              pj = C[0]->plist[k][jj];
              if (pi != pj)
              {
                fillm_times_rec(j,i,pi,pj);
                fillm_times_rec(j,i+1,pj,pi);
                i+=2;
}
            }
          }
          else
          {
            pj = C[0]->addpop[k];
            if (pi != pj)
            {
              fillm_times_rec(j,i,pi,pj);
              fillm_times_rec(j,i+1,pj,pi);
              i+=2;
            }
          }
      }
    }
    assert(i==nummigdirs);
  }
  for (j = 0; j < mrows; j++)
    for (i = 0; i < nummigdirs; i++)
    {
      migration_counts[j][i].do_xyplot = 1;
      migration_counts[j][i].do_logplot = 0;
      migration_counts[j][i].do_trend = 0;
      migration_counts[j][i].do_autoc = 0;
      migration_counts[j][i].plotrescale = 1.0;
      init_value_record (&migration_counts[j][i], 1);
    }
}                               //init_migration_counts

/*
void
init_migration_counts_times (void)
{
  int i, j, k, pi, pj, ji, jj, numhists, mrows, nummigdirs;
  mrows = nloci + (nloci > 1); 
  nummigdirs = 2*(npops-1)*(npops-1);
  numhists = 2 * nummigdirs;
  migration_counts_times = static_cast<value_record **> 
                        (malloc (mrows * sizeof (struct value_record *)));
  for (j = 0; j < mrows; j++)
    migration_counts_times[j] = static_cast<value_record *> 
                        (malloc (numhists * sizeof (struct value_record)));
  for (j = 0; j < mrows; j++)
  {
    i=0;
    for (k= 0; k < lastperiodnumber; k++)
    {
      for (ji=0; ji< npops-k;ji++) 
      {
          pi = C[0]->plist[k][ji];
          if (k==0)
          {
            for (jj=ji+1; jj< npops-k;jj++) 
            {
              pj = C[0]->plist[k][jj];
              if (pi != pj)
              {
                fillmrec(j,i,pi,pj);
                fillmrec(j,i+2,pj,pi);
                i+=4;
              }
            }
          }
          else
          {
            pj = C[0]->addpop[k];
            if (pi != pj)
            {
              fillmrec(j,i,pi,pj);
              fillmrec(j,i+2,pj,pi);
              i+=4;
            }
          }
      }
    }
    assert(i==numhists);
  }
  for (j = 0; j < mrows; j++)
    for (i = 0; i < numhists; i++)
    {
      migration_counts_times[j][i].do_xyplot = 1;
      migration_counts_times[j][i].do_logplot = 0;
      migration_counts_times[j][i].do_trend = 0;
      migration_counts_times[j][i].do_autoc = 0;
      migration_counts_times[j][i].plotrescale = 1.0;
      if (ODD (i))
        init_value_record (&migration_counts_times[j][i], 1);
      else
        init_value_record (&migration_counts_times[j][i], 0);
    }
}                               //init_migration_counts_times
*/


void
init_mutation_scalar_rec (int li)
{
  int i, ui, ai;
  int num_u_update_types = 1;   // number of different types of updates for mutation rate scalars and kappa values
  int num_A_update_types = 1;   // number of different types of updates for A states 
  double uscale;

  // initialize u_rec
  L[li].u_rec = static_cast<chainstate_record_updates_and_values *> 
            (malloc (L[li].nlinked *
                     sizeof (struct chainstate_record_updates_and_values)));
  // set priors and windows
  for (ui = 0; ui < L[li].nlinked; ui++)
  {
    L[li].u_rec[ui].pr.max = log (UMAX);        // 10000
    L[li].u_rec[ui].pr.min = -log (UMAX);
    L[li].u_rec[ui].win = L[li].u_rec[ui].pr.max / nloci;
  }
  uscale = (2 * log (UMAX)) / GRIDSIZE; // 10000
  // name the mutation rate scalars     
  if (L[li].model == STEPWISE)
  {
    for (ai = 0; ai < L[li].nlinked; ai++)
      sprintf (L[li].u_rec[ai].str, "%dSW%d", li, ai);
  }
  else
  {
    sprintf (L[li].u_rec[0].str, "%du ", li);
    if (L[li].model == JOINT_IS_SW)
      for (ai = 1; ai < L[li].nlinked; ai++)
        sprintf (L[li].u_rec[ai].str, "%dSW%d", li, ai);
  }
  for (ui = 0; ui < L[li].nlinked; ui++)
  {
    L[li].u_rec[ui].num_uptypes = num_u_update_types;
    L[li].u_rec[ui].upnames = static_cast<strnl *> 
                    (malloc (L[li].u_rec[ui].num_uptypes * sizeof (strnl)));
    for (i = 0; i < num_u_update_types; i++)
    {
      sprintf (L[li].u_rec[ui].upnames[i], "scalar update");
    }
    L[li].u_rec[ui].upinf = static_cast<update_rate_calc *> 
            (calloc ((size_t) L[li].u_rec[ui].num_uptypes, 
                     sizeof (struct update_rate_calc)));
    L[li].u_rec[ui].num_vals = 1;
    L[li].u_rec[ui].v = static_cast<value_record *> 
            (malloc (L[li].u_rec[ui].num_vals * sizeof (struct value_record)));
    for (i = 0; i < L[li].u_rec[ui].num_vals; i++)
    {
      strcpy (L[li].u_rec[ui].v[i].str, L[li].u_rec[ui].str);   // will this ever get used ? 
      L[li].u_rec[ui].v[i].do_xyplot = 1;
      L[li].u_rec[ui].v[i].plotrescale = uscale;
      L[li].u_rec[ui].v[i].do_logplot = 1;
      L[li].u_rec[ui].v[i].do_autoc = 0;
      L[li].u_rec[ui].v[i].do_trend = 1;
      L[li].u_rec[ui].v[i].plotrange = L[li].u_rec[ui].pr;
      init_value_record (&L[li].u_rec[ui].v[i], 0);
    }
  }

  // do kappa_rec 
  if (L[li].model == HKY)
  {
    L[li].kappa_rec = static_cast<chainstate_record_updates_and_values *> 
                (malloc (sizeof (struct chainstate_record_updates_and_values)));
    L[li].kappa_rec->pr.max = KAPPAMAX;
    L[li].kappa_rec->pr.min = 0.0;
    L[li].kappa_rec->win = 2.0;
    sprintf (L[li].kappa_rec->str, "%d_Ka", li);
    L[li].kappa_rec->num_uptypes = num_u_update_types;
    L[li].kappa_rec->upnames = static_cast<strnl *> 
                            (malloc (num_u_update_types * sizeof (strnl)));
    for (i = 0; i < num_u_update_types; i++)
    {
      sprintf (L[li].kappa_rec->upnames[i], "kappa update");
    }
    L[li].kappa_rec->upinf = static_cast<update_rate_calc *> 
                    (calloc ((size_t) num_u_update_types, 
                             sizeof (struct update_rate_calc)));
    L[li].kappa_rec->num_vals = 1;
    L[li].kappa_rec->v = static_cast<value_record *> 
            (malloc (L[li].kappa_rec->num_vals * sizeof (struct value_record)));
    strcpy (L[li].kappa_rec->v->str, L[li].kappa_rec->str);
    strncpy (L[li].kappa_rec->v->strshort, L[li].kappa_rec->str,
             PARAMSTRLENSHORT - 1);
    L[li].kappa_rec->v->do_xyplot = 1;
    L[li].kappa_rec->v->do_logplot = 0;
    L[li].kappa_rec->v->do_trend = 1;
    L[li].kappa_rec->v->plotrange = L[li].kappa_rec->pr;
    L[li].kappa_rec->v->plotrescale = 1.0;
    L[li].kappa_rec->v->do_autoc = 0;

    init_value_record (L[li].kappa_rec->v, 0);
  }

  // do A_rec 
  // A_rec not used as of sometime in 2010, A gets enough updates when updating genealogy
  if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
  {
    L[li].A_rec = static_cast<chainstate_record_updates_and_values *> 
                (malloc (L[li].nlinked * 
                         sizeof (struct chainstate_record_updates_and_values)));
    for (ai = 0; ai < L[li].nlinked; ai++)
    {
      if (L[li].umodel[ai] == STEPWISE)
        sprintf (L[li].A_rec[ai].str, "%dSW%d", li, ai);
      else
        sprintf (L[li].A_rec[ai].str, "NOTE STEPWISE");
      L[li].A_rec[ai].num_uptypes = num_A_update_types;
      L[li].A_rec[ai].upnames = static_cast<strnl *> 
                            (malloc (num_A_update_types * sizeof (strnl)));
      for (i = 0; i < num_A_update_types; i++)
      {
        sprintf (L[li].A_rec[ai].upnames[i], "STR update");
      }
      L[li].A_rec[ai].upinf = static_cast<update_rate_calc *> 
                        (calloc ((size_t) num_A_update_types, 
                         sizeof (struct update_rate_calc)));
      L[li].A_rec[ai].num_vals = 0;
    }
  }
}                               // init_mutation_scalar_rec

void add_priorinfo_to_output(char priorfilename[],int *fpstri, char fpstr[])
{
  int i;
  SP "\nParameter Priors\n");
  SP "-----------------\n");
  if (!calcoptions[USEPRIORFILE])
  {
    SP "  Population size parameters maximum value : %.4lf \n", thetaprior);
    if (modeloptions[EXPOMIGRATIONPRIOR])
      SP "  Migration rate parameters exponential distribution mean : %.4lf \n",
        mprior);
    else
      SP "  Migration rate parameters maximum value: %.4lf \n", mprior);
    if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0)
    {
      SP "  Splitting time : %.4lf\n", tprior);
    }
    SP "\n");
  }
  else
  {
    SP"  Prior distribution terms given in file: %s\n\n",priorfilename);
    SP"  Splitting time maximum values \n");
    SP"\tPeriod\tPrior maximum value\n");
    for (i=0;i<numsplittimes;i++)
      SP"\t%s\t%.3lf\n",T[i].str,T[i].pr.max);
    SP"  Population size parameters maximum values \n");
    SP"\tPopulation\tPrior maximum value\n");
    for (i=0;i<numtreepops;i++)
      SP"\t%s\t%.3lf\n",itheta[i].str,itheta[i].pr.max);
    if (modeloptions[EXPOMIGRATIONPRIOR]==0)
    {
      SP"  Migration parameters maximum values \n");
      SP"\tMigration rate\tMaximum value\n");
    }
    else
    {
      SP"  Migration parameter prior means (exponential priors)\n");
      SP"\tMigration rate\tPrior mean value\n");
    }
    for (i=0;i<nummigrateparams;i++)
      SP"\t%s\t%.3lf\n",imig[i].str,imig[i].pr.max);
  }
}
void
start_setup_L (char infilename[], int *fpstri, char fpstr[], int currentid)
{
  int li, i;

  /* get the number of loci and the number of populations from the top of the datafile */
  read_datafile_top_lines (infilename, fpstri, fpstr, startpoptreestring);
  /* setup a temporary struture to record how much variation there is,  used for picking starting values of mutation scalars */
  numsitesIS = static_cast<int **> (malloc (nloci * sizeof (int *)));
                                // rows in this matrix are malloced in setup_L
  uvals = static_cast<double **> (malloc (nloci * sizeof (double *)));   
  for (i = 0; i < nloci; i++)
  {
    numsitesIS[i] = static_cast<int *> 
                        (calloc ((size_t) npops + 1, sizeof (int)));
  }
  readdata (infilename, fpstri, fpstr, numsitesIS, currentid);
  for (li = 0; li < nloci; li++)
  {
    init_mutation_scalar_rec (li);
    init_g_rec (li);
    uvals[li] = static_cast<double *> 
                    (malloc (L[li].nlinked * sizeof (double)));
    if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
    {
      init_a_rec (li);
    }
  }
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    init_a_rec_multichain ();
  }
  uterm = set_uvals ();
  if (modeloptions[ADDGHOSTPOP])
  {
    assert (assignmentoptions[POPULATIONASSIGNMENT] == 0);                          
    if (npops + 1 > MAXPOPS)
    {
      IM_err (IMERR_INPUTFILEINVALID, 
              "ghost population makes number of population (%d) greater than MAXPOPS [%d]", 
              npops, MAXPOPS);
    }
     for (li = 0; li < nloci; li++)
     {
      L[li].samppop[npops] = 0;
     }
    npops++;
    numtreepops += 2;
    lastperiodnumber++;
    numsplittimes++;
   }
  gi_largestngenes = 0;
  gi_largestnumsites = 0;
  for (li = 0; li < nloci; li++)
  {
    if (gi_largestngenes < L[li].numgenes)
    {
      gi_largestngenes = L[li].numgenes;
    }
    if (gi_largestnumsites < L[li].numsites)
    {
      gi_largestnumsites = L[li].numsites;
    }
  }
 return;
}                               //start_setup_L

void
start_setup_C (void)
{
  int ci, li;

  C = static_cast<chain **>             //points to an array of chains
            (malloc (numchains * sizeof (struct chain *)));     
  for (ci = 0; ci < numchains; ci++)
  {
    C[ci] = static_cast<chain *> (malloc (sizeof (struct chain)));
    sprintf (C[ci]->name, "%d", ci);
    if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
    {
      C[ci]->nasn = static_cast<int *> (malloc (npops * sizeof (int)));
    }
  }
  for (ci = 0; ci < numchains; ci++)
  {
    init_genealogy_weights (&C[ci]->allgweight);

    /* CR: 110516.2
     * change malloc to calloc so that memory would be initialized to
     * a known state before use.
     */
    C[ci]->G = static_cast<genealogy *> 
                (calloc ((size_t) nloci, sizeof (struct genealogy)));
    for (li = 0; li < nloci; li++)
    {
      init_genealogy_weights (&(C[ci]->G[li].gweight));
    }
    setup_poptree (ci, startpoptreestring);
    set_tvalues (ci);
  }
}                               // start_setup_C


void
finish_setup_C ()
{
  int ci, li, i, ai;
  int nosimmigration;
  init_treeweight ();
  init_gtreecommon ();
  if (modeloptions[NOMIGRATION]==0)
  {
    set_nomigrationchecklist ();
  }
  for (ci = 0; ci < numchains; ci++)
  {
    for (li = 0; li < nloci; li++)
    {
      C[ci]->G[li].gtree = static_cast<edge *> 
                        (malloc (L[li].numlines * sizeof (struct edge)));
      C[ci]->G[li].uvals = static_cast<double *> 
                        (malloc (L[li].nlinked * sizeof (double)));
      C[ci]->G[li].pdg_a = static_cast<double *> 
                        (malloc (L[li].nlinked * sizeof (double)));
    }
  }
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    /* Individual names are set here. */
    imaAsnSet ();
  }

  for (ci = 0; ci < numchains; ci++)
  {
    init_probcalc (&(C[ci]->allpcalc));
    for (li = 0; li < nloci; li++)
    {
      for (i = 0; i < L[li].numlines; i++)
      {
        C[ci]->G[li].gtree[i].mig = static_cast<migstruct *> 
                            (malloc (MIGINC * sizeof (struct migstruct)));
        C[ci]->G[li].gtree[i].cmm = MIGINC;
        C[ci]->G[li].gtree[i].mig[0].mt = -1;
        C[ci]->G[li].gtree[i].up[0] = -1;
        C[ci]->G[li].gtree[i].up[1] = -1;
        C[ci]->G[li].gtree[i].down = -1;
        C[ci]->G[li].gtree[i].time = 0;
        C[ci]->G[li].gtree[i].mut = -1;
        C[ci]->G[li].gtree[i].pop = -1;
        C[ci]->G[li].gtree[i].ei = i;
        C[ci]->G[li].gtree[i].exist = 'T';

        if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
        {
          /* CR: 110516.1
           * change malloc to calloc so that memory would be initialized to
           * a known state before use.
           */
          C[ci]->G[li].gtree[i].A = static_cast<int *> 
                            (calloc ((size_t) L[li].nlinked, sizeof (int)));
          C[ci]->G[li].gtree[i].dlikeA = static_cast<double *> 
                            (calloc ((size_t) L[li].nlinked, sizeof (double)));
        }

        if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
        {
          if (L[li].model == INFINITESITES)
          {
            C[ci]->G[li].gtree[i].seq = static_cast<int *> 
                                (malloc (L[li].numsites * sizeof (int)));
            if (i < L[li].numgenes)
              {
                memcpy (C[ci]->G[li].gtree[i].seq, L[li].seq[i], L[li].numsites * sizeof (int));
              }
          }
        }
      }
      if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
      {
        C[ci]->G[li].mut = static_cast<int *> 
                        (malloc (L[li].numsites * (sizeof (int))));
      }
    }
    /* we need member ei to be ready for initial assignment */
    imaAsnInitAssign (ci);

    if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
    {
      for (i = 0; i < npops; i++)
      {
        assert (C[ci]->poptree[i].b == 0);
        assert (C[ci]->poptree[i].e == 1);
        C[ci]->poptree[i].time = TIMEMAX;
      }
    }
    else
    {
      for (i = 0; i < 2 * npops - 1; i++)
      {
        if (C[ci]->poptree[i].e == -1)
          C[ci]->poptree[i].time = TIMEMAX;
        else
          C[ci]->poptree[i].time = C[ci]->tvals[C[ci]->poptree[i].e - 1];
      }
    }

    if (ci == 0)                //&& nurates > 1)
    {
      setuinfo (uterm);
    }

    for (li = 0; li < nloci; li++)
    {
      if (ci > 0 && nurates > 1)
      {
        for (ai = 0; ai < L[li].nlinked; ai++)
          C[ci]->G[li].uvals[ai] = C[0]->G[li].uvals[ai];
      }
      if (nurates == 1)
      {
        C[ci]->G[li].uvals[0] = 1.0;
      }

      nosimmigration = !(modeloptions[PARAMETERSBYPERIOD] || nummigrateparams == maxpossiblemigrateparams);
      
      switch (L[li].model)
      {
      case INFINITESITES:
        makeIS (ci, li,nosimmigration);
        treeweight (ci, li);
        C[ci]->G[li].pdg = C[ci]->G[li].pdg_a[0] = likelihoodIS (ci, li, C[ci]->G[li].uvals[0]);
        break;
      case HKY:
        if (assignmentoptions[JCMODEL] == 1)
        {
          makeHKY (ci, li,nosimmigration);
          treeweight (ci, li);
          C[ci]->G[li].pdg_a[0] = likelihoodJC (ci, li, C[ci]->G[li].uvals[0]);
          C[ci]->G[li].pdg = C[ci]->G[li].pdg_a[0];
        }
        else
        {
          for (i = 0; i < 4; i++)
          {
            C[ci]->G[li].pi[i] = pi[li][i];
          }
          C[ci]->G[li].kappaval = 2.0;    // starting kappa value
          makeHKY (ci, li,nosimmigration);
          treeweight (ci, li);
          C[ci]->G[li].pdg_a[0] = likelihoodHKY (ci, li, C[ci]->G[li].uvals[0], 
                                                 C[ci]->G[li].kappaval, 
                                                 -1, -1, -1, -1);
          C[ci]->G[li].pdg = C[ci]->G[li].pdg_a[0];
          copyfraclike (ci, li);
          storescalefactors (ci, li);
        }
        break;
      case STEPWISE:
        somestepwise = 1;

        for (ai = 0; ai < L[li].nlinked; ai++)
          for (i = 0; i < L[li].numgenes; i++)
            C[ci]->G[li].gtree[i].A[ai] = L[li].A[ai][i];
        makeSW (ci, li,nosimmigration);
        treeweight (ci, li);
        C[ci]->G[li].pdg = 0;
        for (ai = 0; ai < L[li].nlinked; ai++)
        {
          C[ci]->G[li].pdg_a[ai] = likelihoodSW (ci, li, ai, C[ci]->G[li].uvals[ai], 1.0);
          C[ci]->G[li].pdg += C[ci]->G[li].pdg_a[ai];
        }

        break;
      case JOINT_IS_SW:
        for (ai = 1; ai < L[li].nlinked; ai++)
          for (i = 0; i < L[li].numgenes; i++)
            C[ci]->G[li].gtree[i].A[ai] = L[li].A[ai][i];
        somestepwise = 1;
        makeJOINT_IS_SW (ci, li,nosimmigration);
        treeweight (ci, li);
        C[ci]->G[li].pdg = C[ci]->G[li].pdg_a[0] =
          likelihoodIS (ci, li, C[ci]->G[li].uvals[0]);
        for (ai = 1; ai < L[li].nlinked; ai++)
        {
          C[ci]->G[li].pdg_a[ai] = likelihoodSW (ci, li, ai, C[ci]->G[li].uvals[ai], 1.0);
          C[ci]->G[li].pdg += C[ci]->G[li].pdg_a[ai];
        }
        break;
      }
     
      C[ci]->allpcalc.pdg += C[ci]->G[li].pdg;
      sum_treeinfo (&(C[ci]->allgweight), &(C[ci]->G[li].gweight));
    }
    initialize_integrate_tree_prob (ci, &C[ci]->allgweight, &C[ci]->allpcalc);
  }
  return;
}                               // finish_setup_C


/* set random times for poptree */
  /* first pick trandom times over interval from 0 to tprior 
     use broken stick model - simple dirichlet distribution
   */
/* to fix values:
in imamp.h  turn on these undefs
#undef DO_RY1UPDATE
#undef DO_RY1UPDATE
and at the end of set_tvalues()
use code like 
 C[ci]->tvals[0]=0.1;C[ci]->tvals[1]=0.2;C[ci]->tvals[2]=0.4;C[ci]->tvals[3]=0.8; */

void
set_tvalues (int ci)
{
  int i;
  double sum;
  double times[MAXPOPS];

  if (npops == 1)
  {
    C[ci]->tvals = static_cast<double *> (malloc (sizeof (double)));
    C[ci]->tvals[0] = TIMEMAX;
  }
  else if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
  {
    C[ci]->tvals = static_cast<double *> (malloc (sizeof (double)));
    C[ci]->tvals[0] = TIMEMAX;
  }
  else
  {
    C[ci]->tvals = static_cast<double *> 
                    (malloc ((lastperiodnumber + 1) * sizeof (double)));
    for (sum = 0, i = 0; i < npops; i++)
    {
      times[i] = expo (1.0);
      sum += times[i];
    }

    for (i = 0; i < lastperiodnumber; i++)
    {
      times[i] *= (T[i].pr.max - T[i].pr.min) / sum;
      times[i] += T[i].pr.min;
      if (i > 0)
        times[i] += times[i - 1];
      C[ci]->tvals[i] = times[i];
      assert (C[ci]->tvals[i] < T[i].pr.max && C[ci]->tvals[i] > T[i].pr.min);
    }
    C[ci]->tvals[i] = TIMEMAX;
  }
#ifndef DO_NWUPDATE
#ifndef DO_RY1UPDATE
 if (lastperiodnumber ==2  & T[0].pr.max == 5.0)
 {
   C[ci]->tvals[0]=1.0;C[ci]->tvals[1]=4.0;
 }
 else
 {
  for (i=0;i<lastperiodnumber;i++)
    C[ci]->tvals[i] = (i+1.0)/(lastperiodnumber+1.0) * T[i].pr.max;
 }
 //C[ci]->tvals[0] = 0.01;
 //C[ci]->tvals[1] = 0.02;
#endif
#endif  
  return;
}                               /* set_tvalues */

/* JH 2/18/10 inserted this */
#define MAXTRANGE_IF_ISLANDMODEL  100.0  // have to pick something 
void
setup_T ()
{
  int i, li;
  /* JH 2/18/10  moved this here */ 
  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1  || npops == 1)
  {
    T = NULL;
    tprior = 0.0;
    for (li=0;li<nloci;li++)
    {
      L[li].g_rec->v->plotrange.max = MAXTRANGE_IF_ISLANDMODEL;
      init_value_record (L[li].g_rec->v, 0);
    }
  }
  else
  {
    T = static_cast<chainstate_record_updates_and_values *> 
            (calloc ((size_t) (npops - 1), 
             sizeof (struct chainstate_record_updates_and_values)));
    for (i = 0; i < lastperiodnumber; i++)
    {
      if (!calcoptions[USEPRIORFILE])
        tperiodpriors[i] = tprior;
      T[i].pr.max = tperiodpriors[i];
      T[i].pr.min = 0;
    }
    for (i = 0; i < lastperiodnumber; i++)
    {
      sprintf (T[i].str, "t%d", i);
      T[i].num_uptypes = IM_UPDATE_TIME_NUMBER;
      T[i].upnames = static_cast<strnl *> 
                    (malloc (T[i].num_uptypes * sizeof (strnl)));
#ifdef DO_NWUPDATE
      if (assignmentoptions[POPULATIONASSIGNMENTBF] == 0)
      {
        sprintf (T[i].upnames[IM_UPDATE_TIME_NW], "NielsenWakeley");
      }
      else
      {
        sprintf (T[i].upnames[IM_UPDATE_TIME_NW], "RannalaYang");
      }
#endif
#ifdef DO_RY1UPDATE
      sprintf (T[i].upnames[IM_UPDATE_TIME_RY1], "RannalaYang");
#endif
      T[i].upinf = static_cast<update_rate_calc *> 
                (calloc ((size_t) T[i].num_uptypes, 
                         sizeof (struct update_rate_calc)));
      T[i].num_vals = 1;
      T[i].v = static_cast<value_record *> 
                (malloc (T[i].num_vals * sizeof (struct value_record)));
      strcpy (T[i].v->str, T[i].str);
      strncpy (T[i].v->strshort, T[i].v->str, PARAMSTRLENSHORT - 1);
      T[i].v->do_xyplot = 1;
      T[i].v->do_trend = 1;
      T[i].v->do_logplot = 0;
      T[i].v->do_autoc = 1;
      T[i].v->plotrange.max = tperiodpriors[i];
      T[i].v->plotrange.min = 0;
      T[i].v->plotrescale = 1.0;
      init_value_record (T[i].v, 0);
      /* 2/18/10 JH  inserted this,  needed for plotting TMRCAs */
      if (i==0)
        for (li=0;li<nloci;li++)
        {
          L[li].g_rec->v->plotrange.max = tperiodpriors[lastperiodnumber-1] * TIMEPRIORMULTIPLIER;
          init_value_record (L[li].g_rec->v, 0);
        }
    }
  }

}                               //setup_T



void
finish_setup_L (void)
{
  int li;
  for (li = 0; li < nloci; li++)
  {
    XFREE (numsitesIS[li]);
    XFREE (uvals[li]);
  }
  XFREE (numsitesIS);
  XFREE (uvals);
}                               // finish_setup_L

void reportparamcounts(int *fpstri, char fpstr[])
{
  SP"\nParameter Counts\n");
  SP"----------------\n");
  SP"   Population sizes : %d\n",numpopsizeparams);
  SP"   Migration rates  : %d\n",nummigrateparams);
  SP"   Parameters in the MCMC simulation\n");
  SP"      Splitting times : %d\n",numsplittimes);
  SP"      Mutation scalars: %d\n",nurates); 
  SP"      HKY Kappa (ti/tv) ratios: %d\n",nkappas); 
} // reportparamcounts(void)

/**********  GLOBAL FUNCTIONS  *******/

void
setup (char infilename[], int *fpstri, char fpstr[], char priorfilename[], int currentid)
{
  start_setup_L (infilename, fpstri, fpstr, currentid);

  if (calcoptions[USEPRIORFILE])
  {
    popsizeprior_fromfile = static_cast<double *> 
                        (malloc (numtreepops * sizeof(double)));
    mprior_fromfile = alt2d_alloc2Ddouble(numtreepops,numtreepops);
    readpriorfile(priorfilename,popsizeprior_fromfile,mprior_fromfile);
  }
  
  setup_T ();
  start_setup_C ();
  setup_iparams ();

  if (modeloptions[ANCESTRALPOPSIZESSETTOLARGESTSAMPLED])
  {
    addpopsizeidtext(fpstri,fpstr);
  }
  finish_setup_C ();
  reportparamcounts(fpstri, fpstr);
  finish_setup_L ();            // somethings in L need info from T , free up numsitesIS and uvals

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0
      && npops > 1)
  {
    init_t_NW ();
    init_t_RY ();
  }
  init_updategenealogy ();
  init_lpgpd_v ();
  if (outputoptions[MIGRATEHIST])
    /* 8/26/2011 */
    init_migration_counts ();
    //init_migration_counts_times ();

  add_priorinfo_to_output(priorfilename,fpstri, fpstr);
  init_autoc_pointers ();
  if (calcoptions[USEPRIORFILE])
  {
    XFREE(popsizeprior_fromfile);
    alt2d_free2D(mprior_fromfile);
  }
  /* 5/19/2011 JH adding thermodynamic integration */
  if (calcoptions[CALCMARGINALLIKELIHOOD])
    initmarginlikecalc();
}                               // setup() 
