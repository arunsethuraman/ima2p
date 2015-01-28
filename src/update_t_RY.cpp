/*IMa2p 2009-2015 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, and Arun Sethuraman */
#undef GLOBVARS

#include "imamp.hpp"
#include "update_gtree_common.hpp"
#include "updateassignment.hpp"


/*********** LOCAL STUFF **********/

static struct genealogy_weights holdallgweight_t_RY;
static struct genealogy_weights holdgweight_t_RY[MAXLOCI];
static struct probcalc holdallpcalc_t_RY;
static int largestsamp;
static int **skipflag;

static double beforesplit (int tnode, double oldt, double newt, double tau_u, double ptime);
static double aftersplit (int tnode, double oldt, double newt, double tau_d, double ptime);
static void storegenealogystats_all_loci (int ci, int mode);

static double forwardRY3 (double ptime, /* double oldt, */double r, double t_u_prior);
static double backwardRY3 (double ptime, /* double newt, */double r,double t_u_prior);
/********* LOCAL FUNCTIONS **************/

void
storegenealogystats_all_loci (int ci, int mode)
{
  static double holdlength[MAXLOCI], holdtlength[MAXLOCI];
  static double holdroottime[MAXLOCI];
  static int holdroot[MAXLOCI];
  static int holdmig[MAXLOCI];
  int li;
  if (mode == 0)
  {
    for (li = 0; li < nloci; li++)
    {
      holdlength[li] = C[ci]->G[li].length;
      holdtlength[li] = C[ci]->G[li].tlength;
      holdroottime[li] = C[ci]->G[li].roottime;
      holdroot[li] = C[ci]->G[li].root;
      holdmig[li] = C[ci]->G[li].mignum;
    }
  }
  else
  {
    for (li = 0; li < nloci; li++)
    {
      C[ci]->G[li].length = holdlength[li];
      C[ci]->G[li].tlength = holdtlength[li];
      C[ci]->G[li].mignum = holdmig[li];
      C[ci]->G[li].roottime = holdroottime[li];
      C[ci]->G[li].root = holdroot[li];
    }
  }
  return;
}                               // storegenealogystats  

double
aftersplit (int tnode, double oldt, double newt, double tau_d, double ptime)
{
  if (tnode == lastperiodnumber - 1)
  {
     return ptime + newt - oldt;
  }
  else
  {
    return tau_d - (tau_d - newt) * (tau_d - ptime) / (tau_d - oldt);
  }
}

double
beforesplit (int tnode, double oldt, double newt, double tau_u, double ptime)
{
  if (tnode == 0)
  {
    return ptime * newt / oldt;
  }
  else
  {
    return tau_u + (ptime - tau_u) * (newt - tau_u) / (oldt - tau_u);
  }
}


/*************GLOBAL FUNCTIONS ******************/


void
init_t_RY (void)
{
  int li, j;
  init_genealogy_weights (&holdallgweight_t_RY);
  for (li = 0; li < nloci; li++)
    init_genealogy_weights (&holdgweight_t_RY[li]);
  init_probcalc (&holdallpcalc_t_RY);
  for (largestsamp = 0, j = 0; j < nloci; j++)
    if (largestsamp < L[j].numlines)
      largestsamp = L[j].numlines;
  skipflag = alloc2Dint (nloci, 2 * largestsamp - 1);
}                               // init_changet_RY


void
free_t_RY (void)
{
  int li;
  free_genealogy_weights (&holdallgweight_t_RY);
  for (li = 0; li < nloci; li++)
  {
    free_genealogy_weights (&holdgweight_t_RY[li]);
  }
  free_probcalc (&holdallpcalc_t_RY);
  orig2d_free2D ((void **) skipflag, nloci);

}                               // free_changet_RY


/* 
Notes on changet_RY()  implements updating of Rannala and Yang (2003)   

This application is pretty much the same as theirs - changing times on a species tree that contains a gene tree. 
The big difference is that IM includes migration.   This means that we have to count migration events and change 
migration times in the same was as we change coalescent times and count how many get changed. 

R&Y also only change coalescent times in populations affected by a changed t.  But because we have migration
there is more entanglement between populations.  It seems best to change all times within an interval that is 
affected by a changing t. 


in R&Y usage
u (upper) means older, deeper in the genealogy
l (lower)  means younger, more recent in the genealogy

In R&Y 
Tu next older split time
Tl  next more recent split time
T - current split time
T*  - proposed time 
t - time of some event 
t* - time of that event after update

If  Tu> t > T
t* = Tu - (Tu-t) (Tu-t*)/(Tu-T)

If Tl< t<= T
t* = Tl + (t-tl) (T*-Tl)/(T-Tl)

nl = number of nodes with Tl< t < T
ml = number of migration events with Tl< t < T

nu = number of nodes with Tu> t > T
mu = number of migration events with Tu> t > T

MH criteria  
p(X|G*,t*)p(G*,t*) (Tu-T*)^(nu+mu) (T*-Tl)^(nl+ml)
---------------------------------------------------
p(X|G,t)p(G,t)      (Tu-T)^(nu+mu) (T-Tl)^(nl+ml)


but this causes confusion with jhey usage in which 
u means upper - more recent. 

so here we use  u  for upper (meaning more recent)
use d for down  (older)

tau current split time
tau*  new value 
tau_d - older node time (i.e. time of next oldest node - deeper in the genealogy)
tau_u - more recent node time (i.e. time of next youngest node - more recent in time)
tau_d  > tau > tau_u

if tau is the first node,  then tau_u = 0. 
if tau is the last node, then tau_d = infinity

for an event at time t where tau_u < t < tau_d

if t > tau  see aftersplit()  t is older than the current split time 
t* = tau_d - (tau_d - tau*)(tau_d - t)/(tau_d-tau)  (A7)

if t <= tau  see beforesplit()
t* = tau_u + (tau* - tau_u)(t - tau_u)/(tau - tau_u) (A8)

m is the number of events moved using A7,  n is the number moved using A8

then Hastings term for the update is:
 tau_d - tau*      tau* - tau_u 
(------------)^m  (------------)^n
 tau-u - tau        tau - tau_u

 For IM,  we use the same except m and n include both includes migation and coalescent events

For edges where downtime < tau_d || uptime > tau_d  are not involved in the update
For any edge that spends zero time in either splitpop  or the ancestral pop, during the tau_u/tau_d interval
it is  possible to not update the coalescent time or migration times of 

The difficulty is that it is possible that an uptime for one edge gets moved because the times on one of its daughter edges got moved. 
This means that for checking about skipping an edge, because it is not involved in any
population associated with the splittin time update
we need to check entire groups of branches that descend from 
an edge that crosses the tau_u boundary. 

use a scoring system for edges  
-1 to ignore because out of bounds above
-2 to ignore because out of bounds below
0 to  deal with, edge is relevant
1  to ignore because not in splitpops or ancestral pop

set up a recursive function  
for an edge that crosses the tau_d line,  check to see if that edge
and all descendent edges that to not cross tau_l  are in the 
populations affected by the splitting time update
If all of those edges are not relevant then 
they all get a skipflag value of 1
If any one of them does spend any time in any of the 
populations involved int the population split
then all of the edges have their skipflag value 
set to 0  

*/

/* let u refer to the more recent time  and d to the older time  */
int
changet_RY1 (int ci, int timeperiod)    // after Rannala and Yang (2003)  - rubberband method
{
  double metropolishastingsterm, newt, oldt;
  double pdgnew[MAXLOCI + MAXLINKED], pdgnewsum, pdgoldsum, probgnewsum,
    temppdg;
  double t_u_hterm, t_d_hterm, tpw;
  int li, i, j, ecd, ecu, emd, emu, ai, ui;
  double U;
  struct genealogy *G;
  struct edge *gtree;
  double t_d, t_u, t_u_prior, t_d_prior;
  double holdt[MAXPERIODS];


  if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
  {
    assertgenealogy (ci);
  }

  t_u = (timeperiod == 0) ? 0 : C[ci]->tvals[timeperiod - 1];
  t_d =
    (timeperiod ==
     (lastperiodnumber - 1)) ? TIMEMAX : C[ci]->tvals[timeperiod + 1];
  t_d_prior = DMIN (T[timeperiod].pr.max, t_d);
  t_u_prior = DMAX (T[timeperiod].pr.min, t_u);
  oldt = C[ci]->tvals[timeperiod];
  newt = getnewt (timeperiod, t_u_prior, t_d_prior, oldt, 1);
  
  t_u_hterm = (newt - t_u) / (oldt - t_u);
  if (timeperiod == lastperiodnumber - 1)
  {
    t_d_hterm = 1;
  }
  else
  {
    t_d_hterm = (t_d - newt) / (t_d - oldt);
  }

  copy_treeinfo (&holdallgweight_t_RY, &C[ci]->allgweight);  // try turning this off and forcing all recalculations
  copy_probcalc (&holdallpcalc_t_RY, &C[ci]->allpcalc);
  for (i = 0; i < lastperiodnumber; i++)
    holdt[i] = C[ci]->tvals[i];


  pdgoldsum = C[ci]->allpcalc.pdg;
  setzero_genealogy_weights (&C[ci]->allgweight);
  ecd = ecu = emd = emu = 0;
  pdgnewsum = 0;
  probgnewsum = 0;
  storegenealogystats_all_loci (ci, 0);
  C[ci]->tvals[timeperiod] = newt;
  for (i = 0; i < nurates; i++)
    pdgnew[i] = 0;
  for (li = 0; li < nloci; li++)
  {
    G = &(C[ci]->G[li]);
    gtree = G->gtree;
    copy_treeinfo (&holdgweight_t_RY[li], &G->gweight);
    for (i = 0; i < L[li].numlines; i++)
    {
      if (gtree[i].down != -1)
      {
        if (gtree[i].time <= oldt && gtree[i].time > t_u)

        {
          //assert (skipflag[li][i] == 0);turn off 9/19/10
          gtree[i].time =
            beforesplit (timeperiod, oldt, newt, /* t_d, */ t_u, gtree[i].time);
          assert (gtree[i].time != newt);
          ecu++;
        }
        else
        {
          if (gtree[i].time > oldt && gtree[i].time < t_d)
          {
           // assert (skipflag[li][i] == 0); turn off 9/19/10
            gtree[i].time =
              aftersplit (timeperiod, oldt, newt, t_d, /* t_u, */ gtree[i].time);
            assert (gtree[i].time != newt);
            ecd++;
          }
          //else  do not change the time
        }
        j = 0;
        while (gtree[i].mig[j].mt > -0.5)
        {
          assert (gtree[i].mig[j].mt < C[0]->tvals[lastperiodnumber]);
          if (gtree[i].mig[j].mt <= oldt && gtree[i].mig[j].mt > t_u)
          {
            gtree[i].mig[j].mt =
              beforesplit (timeperiod, oldt, newt, /* t_d, */ t_u,
                           gtree[i].mig[j].mt);
            emu++;
          }
          else
          {
            assert (oldt < C[0]->tvals[lastperiodnumber]);
            if (gtree[i].mig[j].mt > oldt && gtree[i].mig[j].mt < t_d)
            {
              gtree[i].mig[j].mt =
                aftersplit (timeperiod, oldt, newt, t_d, /* t_u, */
                            gtree[i].mig[j].mt);
              emd++;
            }
            // else no need to change the time
          }
          j++;
        }
      }
    }
    if (G->roottime <= oldt && G->roottime > t_u
        /* && skipflag[li][G->root] == 0 turn off 9/19/10*/)
      G->roottime =
        beforesplit (timeperiod, oldt, newt, /* t_d, */ t_u, G->roottime);
    else if (G->roottime > oldt && G->roottime < t_d
            /* && skipflag[li][G->root] == 0 turn off 9/19/10*/)
      G->roottime =
        aftersplit (timeperiod, oldt, newt, t_d, /* t_u, */ G->roottime);
    setzero_genealogy_weights (&G->gweight);
        
    treeweight (ci, li);

    sum_treeinfo (&C[ci]->allgweight, &G->gweight);
    ai = 0;
    ui = L[li].uii[ai];

    switch (L[li].model)
    {
      assert (pdgnew[ui] == 0);
    case HKY:
      if (assignmentoptions[JCMODEL] == 1)
      {
        temppdg = pdgnew[ui] =
          likelihoodJC (ci, li, G->uvals[0]);
      }
      else
      {
        temppdg = pdgnew[ui] =
          likelihoodHKY (ci, li, G->uvals[0], G->kappaval, -1, -1, -1, -1);
      }
      break;
    case INFINITESITES:
      temppdg = pdgnew[ui] = likelihoodIS (ci, li, G->uvals[0]);
      break;
    case STEPWISE:
      temppdg = 0;
      for (; ai < L[li].nlinked; ai++)
      {
        ui = L[li].uii[ai];
        assert (pdgnew[ui] == 0);
        pdgnew[ui] = likelihoodSW (ci, li, ai, G->uvals[ai], 1.0);
        temppdg += pdgnew[ui];
      }
      break;
    case JOINT_IS_SW:
      temppdg = pdgnew[ui] = likelihoodIS (ci, li, G->uvals[0]);
      for (ai = 1; ai < L[li].nlinked; ai++)
      {
        ui = L[li].uii[ai];
        assert (pdgnew[ui] == 0);
        pdgnew[ui] = likelihoodSW (ci, li, ai, G->uvals[ai], 1.0);
        temppdg += pdgnew[ui];
      }
      break;
    }
    pdgnewsum += temppdg;
  }

  assert (!ODD (ecd));
  assert (!ODD (ecu));
  ecd /= 2;
  ecu /= 2;
  integrate_tree_prob (ci, &C[ci]->allgweight, &holdallgweight_t_RY,
                       &C[ci]->allpcalc, &holdallpcalc_t_RY, &holdt[0]);   // try enforcing full cacullation
  tpw = gbeta * (pdgnewsum - pdgoldsum);
/* 5/19/2011 JH adding thermodynamic integration */
  if (calcoptions[CALCMARGINALLIKELIHOOD])
  {
    metropolishastingsterm = beta[ci] * tpw + (C[ci]->allpcalc.probg - holdallpcalc_t_RY.probg) + (ecd + emd) * log (t_d_hterm) + (ecu + emu) * log (t_u_hterm);
  }
  else
  {
  tpw += C[ci]->allpcalc.probg - holdallpcalc_t_RY.probg;
    metropolishastingsterm = beta[ci] * tpw + (ecd + emd) * log (t_d_hterm) +   (ecu + emu) * log (t_u_hterm);
  }
  //assert(metropolishastingsterm >= -1e200 && metropolishastingsterm < 1e200);
  U = log (uniform ());
  if (U < DMIN(1.0, metropolishastingsterm))  //9/13/2010 
  //if (metropolishastingsterm >= 0.0 || metropolishastingsterm > U)
  {
    for (li = 0; li < nloci; li++)
    {
      C[ci]->G[li].pdg = 0;
      for (ai = 0; ai < L[li].nlinked; ai++)
      {
        C[ci]->G[li].pdg_a[ai] = pdgnew[L[li].uii[ai]];
        C[ci]->G[li].pdg += C[ci]->G[li].pdg_a[ai];
      }
      if (L[li].model == HKY)
      {
        storescalefactors (ci, li);
        copyfraclike (ci, li);
      }
    }
    C[ci]->allpcalc.pdg = pdgnewsum;
    C[ci]->poptree[C[ci]->droppops[timeperiod + 1][0]].time =
      C[ci]->poptree[C[ci]->droppops[timeperiod + 1][1]].time = newt;

    if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
    {
      assertgenealogy (ci);
    }
    return 1;
  }
  else
  {
    copy_treeinfo (&C[ci]->allgweight, &holdallgweight_t_RY);
    copy_probcalc (&C[ci]->allpcalc, &holdallpcalc_t_RY);
    assert (pdgoldsum == C[ci]->allpcalc.pdg);
    C[ci]->tvals[timeperiod] = oldt;
    for (li = 0; li < nloci; li++)
    {
      G = &(C[ci]->G[li]);
      gtree = G->gtree;
      storegenealogystats_all_loci (ci, 1);
      copy_treeinfo (&G->gweight, &holdgweight_t_RY[li]);
      for (i = 0; i < L[li].numlines; i++)
      {
        if (gtree[i].down != -1)
        {
          if (gtree[i].time <= newt && gtree[i].time > t_u)
          {
           // assert (skipflag[li][i] == 0); turned off 9/19/10
            gtree[i].time =
              beforesplit (timeperiod, newt, oldt, /* t_d, */ t_u, gtree[i].time);
            //cecu++;
          }

          else
          {
            if (gtree[i].time > newt && gtree[i].time < t_d)
            {
             //assert (skipflag[li][i] == 0); turned off 9/19/10
              gtree[i].time =
                aftersplit (timeperiod, newt, oldt, t_d, /* t_u, */ gtree[i].time);
              //cecl++;
            }
          }
          j = 0;
          while (gtree[i].mig[j].mt > -0.5)
          {
            if (gtree[i].mig[j].mt <= newt && gtree[i].mig[j].mt > t_u)
            {
              gtree[i].mig[j].mt =
                beforesplit (timeperiod, newt, oldt, /* t_d, */
                             t_u, gtree[i].mig[j].mt);
              //cemu++;
            }
            else if (gtree[i].mig[j].mt > newt && gtree[i].mig[j].mt < t_d)
            {
              gtree[i].mig[j].mt =
                aftersplit (timeperiod, newt, oldt, t_d, /* t_u, */
                            gtree[i].mig[j].mt);
              //ceml++;
            }
            j++;
          }
        }
      }
//        assert(fabs(C[ci]->G[li].gtree[  C[ci]->G[li].gtree[C[ci]->G[li].root].up[0]].time - C[ci]->G[li].roottime) < 1e-8);    
    }
    /*    assert(ecu==cecu/2);
       assert(ecd==cecl/2);
       assert(emu==cemu);
       assert(emd==ceml); */
    for (li = 0; li < nloci; li++)
    {
      if (L[li].model == HKY)
        restorescalefactors (ci, li);
      /* have to reset the dlikeA values in the genealogies for stepwise model */
      if (L[li].model == STEPWISE)
        for (ai = 0; ai < L[li].nlinked; ai++)
          likelihoodSW (ci, li, ai, C[ci]->G[li].uvals[ai], 1.0);
      if (L[li].model == JOINT_IS_SW)
        for (ai = 1; ai < L[li].nlinked; ai++)
          likelihoodSW (ci, li, ai, C[ci]->G[li].uvals[ai], 1.0);
      // assert(fabs(C[ci]->G[li].gtree[  C[ci]->G[li].gtree[C[ci]->G[li].root].up[0]].time - C[ci]->G[li].roottime) < 1e-8);    
    }
    if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
    {
      assertgenealogy (ci);
    }
    return 0;
  }
}                               /* changet_RY1 */

