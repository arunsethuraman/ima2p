/*IMa2p 2009-2015 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, and Arun Sethuraman */
#undef GLOBVARS
#include "imamp.hpp"
#include "update_gtree_common.hpp"
#include "update_gtree.hpp"
#include "updateassignment.hpp"

/*********** LOCAL STUFF **********/

struct edgemiginfo oldedgemig;
struct edgemiginfo oldsismig;
struct edgemiginfo newedgemig;
struct edgemiginfo newsismig;

/* declarded in update_gtree_common.h  not used here 
struct edge *copyedge; 
extern int holddownA[MAXLINKED];
extern int medgedrop;
extern double lmedgedrop;
extern double holdsisdlikeA[MAXLINKED];
extern struct genealogy holdgtree;
*/

static int mrootdrop;
static double lmrootdrop;
static struct genealogy_weights holdgweight_updategenealogy;
static struct genealogy_weights holdallgweight_updategenealogy;
static struct probcalc holdallpcalc_updategenealogy;
int rootmove;                   /* used in update_gtree.c and update_gtree_covar.c */ // 12/16/09 JH note that update_gree_covar.c no longer used
static double mlist[2 * ABSMIGMAX];     // used in mwork_single_edge() mlist very wasteful of space - could use dynamic memory and checkmig() 
static int mplist[2 * ABSMIGMAX];

/* find the time when two populations join */
double
findjointime (int ci, int slidepop, int sispop, double edgeuptime,
              double sisuptime)
{
  int edgeperiod, sisperiod;
  double jointime;
  edgeperiod = findperiod (ci, edgeuptime);
  sisperiod = findperiod (ci, sisuptime);
  while (edgeperiod < sisperiod)

  {
    edgeperiod++;
    if (slidepop == C[ci]->droppops[edgeperiod][0]
        || slidepop == C[ci]->droppops[edgeperiod][1])
      slidepop = C[ci]->poptree[slidepop].down;
  }

  while (sisperiod < edgeperiod)
  {
    sisperiod++;
    if (sispop == C[ci]->droppops[sisperiod][0]
        || sispop == C[ci]->droppops[sisperiod][1])
      sispop = C[ci]->poptree[sispop].down;
  }

  // at this point edgeperiod == sisperiod 
  while (slidepop != sispop)
  {
    edgeperiod++;
    if (slidepop == C[ci]->droppops[edgeperiod][0]
        || slidepop == C[ci]->droppops[edgeperiod][1])
      slidepop = C[ci]->poptree[slidepop].down;
    if (sispop == C[ci]->droppops[edgeperiod][0]
        || sispop == C[ci]->droppops[edgeperiod][1])
      sispop = C[ci]->poptree[sispop].down;
  }

  if (edgeperiod == 0)
    jointime = 0;
  else
    jointime = C[ci]->tvals[edgeperiod - 1];
  return jointime;
}                               //findjointime

void
slider_nomigration (int ci, int li, int slidingedge, int *sis,
                    double *timepoint, double *slidedist)
/* this is very similar to slider, but with a an extra if/else on the slides up for the case of zero migration */
/* timepoint points at C[ci]->G[li].gtree[*slidingedge].time and is the current position of the sliding point, slidedist is the distance it must move 
with multiple populations,  and no migration, the sliding edge can only be in the population it started in, or an ancestral population
so the same goes for the point on the edge on which the slide is currently at. 
At the beginning of a slide the point is necessarily valid  - All down slides are valid for any distance
the upper limit at any point is the 
maximum of the top of the sliding edge and the time at which the sliding edge and the sister edge are in different populations  
if distance is negative  move up
	determine the upper limit ( MAX(top of edge, split time of sis and edge,  time of node of sis) )
	if upper limit is at a node,  pick left or right  at random
		call slider and continue up
	else reflect  (switch sign on remaining distance)
		call slider and and move down
else  move down
	if reach a node,  pick down or up at random
		if down,  call slider and continue down
		if up,  switch sign on remaining distance
			call slider and move up
kinds of upper limits
 top of sliding edge
 top of sister edge (a node) 
 beginning of period when the population of the edge and the sister branch come together into the same population 
*/
{
  double edgeuptime, sisuptime, popjointime;
  struct edge *gtree = C[ci]->G[li].gtree;
  int slidepop, sispop;
  if (*slidedist < 0)
  {
    /* go up */
    *slidedist = -*slidedist;
    if (slidingedge < L[li].numgenes)
      edgeuptime = 0;
    else
      edgeuptime = gtree[gtree[slidingedge].up[0]].time;

    if (*sis < L[li].numgenes)
      sisuptime = 0;
    else
      sisuptime = gtree[gtree[*sis].up[0]].time;

    assert (*timepoint >= edgeuptime);
    slidepop = gtree[slidingedge].pop;
    sispop = gtree[*sis].pop;
    if (slidepop != sispop)
      popjointime = findjointime (ci, slidepop, sispop, edgeuptime, sisuptime);
    else
      popjointime = 0;

    if (popjointime > edgeuptime && popjointime > sisuptime)
    {
      if (*slidedist < *timepoint - popjointime)
        /* slide up and stop,  sis remains the same */
      {
        *timepoint -= *slidedist;
        *slidedist = 0;
        assert (*timepoint > popjointime);
        return;
      }
      else
      {

        /* slide up and reflect, sis remains the same, leave slidedist positive so slidingedge goes down with next call to slider */
        *slidedist -= *timepoint - popjointime;
        *timepoint = popjointime;
        assert (*slidedist > 0);
        slider_nomigration (ci, li, slidingedge, sis, timepoint, slidedist);
        return;
      }
    }
    else
    {
      if (sisuptime == 0 || edgeuptime >= sisuptime)
      {
        if (*slidedist < *timepoint - edgeuptime)

          /* slide up and stop,  sis remains the same */
        {
          *timepoint -= *slidedist;
          *slidedist = 0;
          assert (*timepoint > edgeuptime);
          return;
        }
        else
        {
          /* slide up and reflect, sis remains the same, leave slidedist positive so slidingedge goes down with next call to slider */
          *slidedist -= *timepoint - edgeuptime;
          *timepoint = edgeuptime;
          assert (*slidedist > 0);
          slider_nomigration (ci, li, slidingedge, sis, timepoint, slidedist);
          return;
        }
      }
      else
      {
        /* edgeuptime is less than sis up time, and thus slidingedge can reach a node */
        if (*slidedist < *timepoint - sisuptime)

          /* slide up and stop,  sis remains the same */
        {
          *timepoint -= *slidedist;
          *slidedist = 0;
          assert (*timepoint > sisuptime);
          return;
        }
        else
        {

          /* slide up and reach a node, pick one side at random and recurse */
          *slidedist -= *timepoint - sisuptime;
          *timepoint = sisuptime;
          if (bitran () /*uniform() < 0.5 */ )
          {
            *sis = gtree[*sis].up[0];
          }
          else
          {
            *sis = gtree[*sis].up[1];
          }

          /* reset slidedist to negative, so slidingedge continues up the gtree in next call to slider */
          *slidedist = -*slidedist;
          assert (*slidedist < 0);
          slider_nomigration (ci, li, slidingedge, sis, timepoint, slidedist);
          return;
        }
      }
    }
  }
  else
  {
    /* go down */
    if (gtree[*sis].down == -1 || *timepoint + *slidedist < gtree[*sis].time)
    {

      /* if sis is the root, or distance is less than to next down node, just slide down */
      *timepoint += *slidedist;
      if (*timepoint >= TIMEMAX)
        *timepoint = TIMEMAX;
      *slidedist = 0;
      return;
    }
    else
    {

      /* a down node is reached */
      *slidedist -= gtree[*sis].time - *timepoint;
      *timepoint = gtree[*sis].time;
      if (bitran ())
      {
        /* begin to slide down the down node */
        *sis = gtree[*sis].down;
        slider_nomigration (ci, li, slidingedge, sis, timepoint, slidedist);
        return;
      }
      else
      {
        /* begin to slide up the sis  */
        if (gtree[gtree[*sis].down].up[0] == *sis)
        {
          *sis = gtree[gtree[*sis].down].up[1];
        }
        else
        {
          *sis = gtree[gtree[*sis].down].up[0];
        }
        *slidedist = -*slidedist;
        slider_nomigration (ci, li, slidingedge, sis, timepoint, slidedist);
        return;
      }
    }
  }
}                               /* slider_nomigration */

void
slider (int ci, int li, int slidingedge, int *sis, double *timepoint,
        double *slidedist)
/* this is not ready for case of multiple populations and zero migration */
/* timepoint points at C[ci]->G[li].gtree[*slidingedge].time and is the current position of the sliding point, 
slidedist is the distance it must move 
do not restructure the gtree. just figure out when and on which branch timepoint ends up on, 
this will be the new sisterbranch
use recursion */
{
  double uplimit;
  struct edge *gtree = C[ci]->G[li].gtree;
  if (*slidedist < 0)
  {
    /* go up */
    *slidedist = -*slidedist;
    if (slidingedge < L[li].numgenes)
      uplimit = 0;
    else
      uplimit = gtree[gtree[slidingedge].up[0]].time;
    assert (*timepoint >= uplimit);

    /* if uplimit >= sis up time  - slidingedge cannot reach a node */
    if (gtree[*sis].up[0] == -1 || uplimit >= gtree[gtree[*sis].up[0]].time)
    {
      if (*slidedist < *timepoint - uplimit)
        /* slide up and stop,  sis remains the same */
      {
        *timepoint -= *slidedist;
        *slidedist = 0;
        assert (*timepoint > uplimit);
        return;
      }
      else
      {
        /* slide up and reflect, sis remains the same, leave slidedist positive so slidingedge goes down with next call to slider */
        *slidedist -= *timepoint - uplimit;
        *timepoint = uplimit;
        assert (*slidedist > 0);
        slider (ci, li, slidingedge, sis, timepoint, slidedist);
        return;
      }
    }
    else
    {
      /* uplimit is less than sis up time, and thus slidingedge can reach a node */
      if (*slidedist < *timepoint - gtree[gtree[*sis].up[0]].time)
        /* slide up and stop,  sis remains the same */
      {
        *timepoint -= *slidedist;
        *slidedist = 0;
        assert (*timepoint > gtree[gtree[*sis].up[0]].time);
        return;
      }
      else
      {
        /* slide up and reach a node, pick one side at random and recurse */
        *slidedist -= *timepoint - gtree[gtree[*sis].up[0]].time;
        *timepoint = gtree[gtree[*sis].up[0]].time;
        if (bitran () /*uniform() < 0.5 */ )
        {
          *sis = gtree[*sis].up[0];
        }
        else
        {
          *sis = gtree[*sis].up[1];
        }

        /* reset slidedist to negative, so slidingedge continues up the gtree in next call to slider */
        *slidedist = -*slidedist;
        assert (*slidedist < 0);
        slider (ci, li, slidingedge, sis, timepoint, slidedist);
        return;
      }
    }
  }
  else
  {
    /* go down */
    if (gtree[*sis].down == -1 || *timepoint + *slidedist < gtree[*sis].time)
    {
      /* if sis is the root, or distance is less than to next down node, just slide down */
      *timepoint += *slidedist;
      if (*timepoint >= TIMEMAX)
        *timepoint = TIMEMAX;
      *slidedist = 0;
      return;
    }
    else
    {
      /* a down node is reached */
      *slidedist -= gtree[*sis].time - *timepoint;
      *timepoint = gtree[*sis].time;
      if (bitran ())
      {
        /* begin to slide down the down node */
        *sis = gtree[*sis].down;
        slider (ci, li, slidingedge, sis, timepoint, slidedist);
        return;
      }
      else
      {
        /* begin to slide up the sis  */
        if (gtree[gtree[*sis].down].up[0] == *sis)
        {
          *sis = gtree[gtree[*sis].down].up[1];
        }
        else
        {
          *sis = gtree[gtree[*sis].down].up[0];
        }
        *slidedist = -*slidedist;
        slider (ci, li, slidingedge, sis, timepoint, slidedist);
        return;
      }
    }
  }
}                               /* slider */

void
joinsisdown (int ci, int li, int sis, int *tmrcachange)
{

  /* extend sis, and XFREE up the down edge */
  int i, j, ai, downdown, down;
  double uptime;
  struct edge *gtree = C[ci]->G[li].gtree;
  down = gtree[sis].down;
  i = 0;
  while (gtree[sis].mig[i].mt > -0.5)
    i++;
  j = -1;

  do
  {
    j++;
    checkmig (i + 1, &(gtree[sis].mig), &(gtree[sis].cmm));
    gtree[sis].mig[i] = gtree[down].mig[j];
    i++;
  } while (gtree[down].mig[j].mt > -0.5);
  gtree[sis].time = gtree[down].time;

  /* set the up to which sis now connects */
  gtree[sis].down = gtree[down].down;
  downdown = gtree[sis].down;
  if (downdown != -1)
  {
    rootmove = 0;
    if (gtree[downdown].up[0] == down)
      gtree[downdown].up[0] = sis;
    else
      gtree[downdown].up[1] = sis;
    mrootdrop = 0;
    lmrootdrop = 0;
  }
  else
  {
    rootmove = 1;
    *tmrcachange += 1;
/* figure out total time and number of migrants being dropped */
    i = -1;
    do
    {
      i++;
    } while (gtree[sis].mig[i].mt > -0.5);

    /* mrootdrop and lmrootdrop do not seem to do anything. We may want to
     * delete two variables? */
    mrootdrop = i; 
    if (sis < L[li].numgenes)
    {
      uptime = 0;
    }
    else
    {
      uptime = gtree[gtree[sis].up[0]].time;
    }

    if (uptime < C[ci]->tvals[lastperiodnumber - 1])
    {
      if (C[ci]->G[li].roottime < C[ci]->tvals[lastperiodnumber - 1])
        lmrootdrop = C[ci]->G[li].roottime - uptime;
      else
        lmrootdrop = C[ci]->tvals[lastperiodnumber - 1] - uptime;
    }
    else
    {
      lmrootdrop = 0;
    }

    C[ci]->G[li].root = sis;
    gtree[sis].down = -1;
    gtree[sis].time = TIMEMAX;
    if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
      for (ai = (L[li].model == JOINT_IS_SW); ai < L[li].nlinked; ai++)
        gtree[sis].dlikeA[ai] = 0;
    gtree[sis].mig[0].mt = -1;
    C[ci]->G[li].roottime = uptime;
  }
}                               /* joinsisdown */

void
splitsisdown (int ci, int li, int slidingedge, int down, int newsis)
{

  /* split newsis into two parts, and make a new down edge out of the lower part */
  int i, j, downdown, nowpop;
  double curt;
  struct edge *gtree = C[ci]->G[li].gtree;
  curt = gtree[slidingedge].time;
  gtree[down].time = gtree[newsis].time;
  gtree[newsis].time = curt;

  /* set the up  of the edge to which down now connects, depends on whether newsis is the root */
  downdown = gtree[newsis].down;
  if (downdown != -1)
  {
    if (gtree[downdown].up[0] == newsis)
      gtree[downdown].up[0] = down;
    else
      gtree[downdown].up[1] = down;
  }
  else
  {
    /* newsis is the current root so the root must move down */
    C[ci]->G[li].root = down;
    C[ci]->G[li].roottime = curt;
    rootmove = 1;
    if (C[ci]->G[li].roottime > TIMEMAX)
      IM_err(IMERR_ROOTTIMEMAXFAIL, "roottime greater than TIMEMAX, chain: %d,locus: %d, roottime %lf, TIMEMAX %lf",ci,li,C[ci]->G[li].roottime,TIMEMAX);
    gtree[down].mig[0].mt = -1;
    if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
      for (i = (L[li].model == JOINT_IS_SW); i < L[li].nlinked; i++)
        gtree[down].dlikeA[i] = 0;
  }
  gtree[down].down = downdown;

  /* divide the migration events along newsis into upper part for newsis and lower part for down */
  /* this might have bugs setting the population of gtree[down] */
  i = 0;
  while (gtree[newsis].mig[i].mt > -0.5 && gtree[newsis].mig[i].mt < curt)
    i++;
  if (i > 0)
    nowpop = gtree[newsis].mig[i - 1].mp;
  else
    nowpop = gtree[newsis].pop;
  j = findperiod (ci, curt);
  /* SANGCHUL: Thu Oct  8 14:14:31 EDT 2009
   * We could use
   * nowpop = saC.popndown[nowpop][j];
   * instead of using the followng while statement. 
   * */
  while (C[ci]->poptree[nowpop].e <= j && C[ci]->poptree[nowpop].e != -1)
    nowpop = C[ci]->poptree[nowpop].down;
  gtree[down].pop = nowpop;
  j = 0;
  if (downdown != -1)
  {
    while (gtree[newsis].mig[j + i].mt > -0.5)
    {
      checkmig (j, &(gtree[down].mig), &(gtree[down].cmm));
      gtree[down].mig[j] = gtree[newsis].mig[j + i];
      assert (nowpop != gtree[down].mig[j].mp);
      nowpop = gtree[down].mig[j].mp;
      j++;
    }
  }
  gtree[newsis].mig[i].mt = -1;
  gtree[down].mig[j].mt = -1;
  gtree[newsis].down = gtree[slidingedge].down = down;
  gtree[down].up[0] = newsis;
  gtree[down].up[1] = slidingedge;
  return;
}                               /* splitsisdown */


/* called by addmigration(),  does most of the migration either by adding events or by calling mwork_single_edge() 
    mwork_single_edge is called for the simpler cases, 
    mwork_two_edges is called when two edges need updating
     */

void  getm (int ci,struct edgemiginfo *edgem,struct edgemiginfo *sisem, struct edgemiginfo *oldedgem,struct edgemiginfo *oldsisem)
{
  int lastmigperiod;
//  int i, ii;

/*for (ii = 0;ii<2;ii++)
for (i=0;i<5;i++)
{
  checkptype[ii][i] = -100;
  checkrr[ii][i] = 0;
}  //8_30_10 */

  lastmigperiod = IMIN(edgem->e,lastperiodnumber-1);
  if ( sisem->edgeid == -1  /* sisem->mtall <= 0*/)  // no sister edge, or sister edge not in a period where migration can occur
  {
    assert(sisem->mtimeavail[0] == 0);
    edgem->mpall = mwork_single_edge (ci, edgem,oldedgem, lastmigperiod); 
  }
  else
  {
    assert (edgem->e == sisem->e);
    if (edgem->mtall <= 0)   // edge has no length in a period with migration,  so just do sister edge
    {
      assert(edgem->mtimeavail[0] == 0);
      sisem->mpall = mwork_single_edge (ci, sisem,oldsisem, lastmigperiod);
    }
    else
    {  // both edge and sis have length in periods with migration
      mwork_two_edges(ci, edgem, sisem,oldedgem, oldsisem, lastmigperiod, &edgem->mpall, &sisem->mpall);

    }
  }
}  //getm

/* add migration to edge, and to its sister if edge connects to the root, 
 * return the log of the hastings ratio of update probabilities 
 *
 * add migration events to edge that just slid.  Also if it slid down the 
 * root node, and moved the root, then migration events may need to be 
 * added to the sister branch as well 
 * 
 * CR 111006.1 10/6/2011 JH   removed oldmigcount, oltlength, newmigcount 
 * and newtlength these had once been used for calculating the migration 
 * weight but this is now done using the  edgemiginfo structures passed 
 * to and from getmprob()
 *
 * Additional note: local vars mtime and mcount also removed.
 */
/* 10/6/2011 JH   removed oldmigcount, oltlength, newmigcount and newtlength
  these had once been used for calculating the migration weight
  but this is now done using the  edgemiginfo structures passed to and from getmprob()*/

double
addmigration (int ci, int li)
{
  int newsis, edge;
  double weight;
  double mproposenum, mproposedenom, temp;
  struct edge *gtree = C[ci]->G[li].gtree;

  assert (C[ci]->G[li].mignum >= 0 && C[ci]->G[li].tlength > 0);
  IMA_reset_edgemiginfo (&newedgemig);
  IMA_reset_edgemiginfo (&newsismig);
  newedgemig.edgeid = edge = oldedgemig.edgeid;
  newedgemig.li = li;
  if (edge < L[li].numgenes)
    newedgemig.upt = 0;
  else
    newedgemig.upt = gtree[gtree[edge].up[0]].time;
  newedgemig.fpop = gtree[gtree[edge].down].pop;
  newedgemig.pop = newedgemig.temppop = gtree[edge].pop;
  newedgemig.dnt = gtree[edge].time;
  newedgemig.mig[0].mt = -1;
  fillmiginfoperiods (ci, &newedgemig);
  if (gtree[edge].down == C[ci]->G[li].root)    /* simulate migration on the sister branch as well */
  {
    //IMA_reset_edgemiginfo (&newsismig);
    newedgemig.fpop = -1;       //pop unknown, as edge must be determined by migration 
    if (gtree[gtree[edge].down].up[0] == edge)
      newsis = gtree[gtree[edge].down].up[1];
    else
      newsis = gtree[gtree[edge].down].up[0];
    if (newsis < L[li].numgenes)
      newsismig.upt = 0;
    else
      newsismig.upt = gtree[gtree[newsis].up[0]].time;
    newsismig.edgeid = newsis;
    newsismig.li = li;
    newsismig.fpop = -1;        //pop unknown, as edge must be determined by migration 
    newsismig.pop = newsismig.temppop = gtree[newsis].pop;
    newsismig.dnt = gtree[newsis].time;
    newsismig.mig[0].mt = -1;
    fillmiginfoperiods (ci, &newsismig);
  }
  else
  {
    newedgemig.fpop = gtree[gtree[edge].down].pop;
    newsismig.edgeid = -1;
    newsismig.mtall = 0; 
    newsismig.b = newsismig.e = -1;// no second edge to deal with
  }

  assert((newsismig.mtall > 0 && newsismig.edgeid >= 0) || newsismig.mtall == 0);


  getm (ci, &newedgemig, &newsismig, &oldedgemig, &oldsismig);

  if (newedgemig.fpop == -1)
  {
    assert(newsismig.edgeid != -1);
    assert(newedgemig.e == lastperiodnumber);
    newedgemig.fpop  = newsismig.fpop  = C[ci]->G[li].root;

  }
  if (gtree[newedgemig.edgeid].down == C[ci]->G[li].root)
    gtree[C[ci]->G[li].root].pop = newedgemig.fpop;


    temp = getmprob (ci, &newedgemig, &newsismig, &oldedgemig, &oldsismig);

  mproposedenom = temp;
  assert (temp > -1e200 && temp < 1e200);

  /* calculate probability of reverse update    */
  temp = getmprob (ci, &oldedgemig, &oldsismig, &newedgemig, &newsismig);
  mproposenum = temp;
  assert (temp > -1e200 && temp < 1e200);
  weight = mproposenum - mproposedenom;
  return weight;
}                               /* addmigration */


/********GLOBAL FUNCTIONS *******/

void
init_updategenealogy (void)
{
  init_genealogy_weights (&holdallgweight_updategenealogy);
  init_genealogy_weights (&holdgweight_updategenealogy);
  init_probcalc (&holdallpcalc_updategenealogy);
}                               // init_updategenealogy


void
free_updategenealogy (void)
{
  free_genealogy_weights (&holdallgweight_updategenealogy);
  free_genealogy_weights (&holdgweight_updategenealogy);
  free_probcalc (&holdallpcalc_updategenealogy);
}                               //free_updategenealogy

#define SLIDESTDVMAX 20   // why this ??

/* steps in picking a new genealogy 
- pick an edge, the bottom of which will slide
- save all info for that edge, sis and the down edge that will be freed up
- join the sis and down edges, freeing up an edge
	if down was the root, then sis becomes the root
-set aside the number for the down edge - this is the freed up edge and will get used again later
-do the sliding and pick a new sis and a location - but do not change the gtree. 
-split the newsis edge into sis and down edges, and divide the migration events accordingly
-connect the original edge to the new spot
-add migration events to this edge
- calculate the probability of the update, in forwrad and reverse directions
- calculate the total Metropolis Hastings terms for genealogy update,  accept or reject 
*/
/* for gtreeprint calls,  use callingsource = 0 */


/* 9/25/08  updated this
revised genealogy updating so that the current migration rate is based on the current number of migration events and the 
current length of the branch that is being updated 
reasoned that this might work better than using the rate that occurs for the entitre genealogy - e.g. help to avoid promoting
correlations and improve mixing  */

/* 7/15/2010  JH revised updategenealogy() and addmigration() and descendant functions to fix a problem associated with the 
Hastings term calculation of the migration path */ 


//prune this   9/21/2010

/* CR 111006.1 10/6/2011 JH   removed local variables mpart and tlengthpart
 * these had once been used for calculating the migration weight
 * but this is now done using the  edgemiginfo structures passed to
 * and from getmprob() in addmigration()
 */ 
int
updategenealogy (int ci, int li, int *topolchange, int *tmrcachange)
{
  int ai, i;
  int edge, oldsis, newsis, freededge, accp;
  double newpdg, newpdg_a[MAXLINKED];
  double migweight, metropolishastingsterm, U;
  double tpw;
  double Aterm[MAXLINKED], Atermsum;
  double slidedist;
  double slideweight, holdslidedist, slidestdv;
  struct genealogy *G = &(C[ci]->G[li]);
  struct edge *gtree = G->gtree;
  int rejectIS;
  double like;
  double holdt[MAXPERIODS];

  if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
  {
    assertgenealogyloc (ci, li);
  }
// initialize and make copies structures that hold quantities for calculating prob of genealogy
  copy_treeinfo (&holdgweight_updategenealogy, &G->gweight);
  copy_treeinfo (&holdallgweight_updategenealogy, &C[ci]->allgweight);
  // store summary stats of the genealogy
  storegenealogystats (ci, li, 0);
  for (i = 0; i < lastperiodnumber; i++)
    holdt[i] = C[ci]->tvals[i];  // JH is this actually necessary ? 
  *tmrcachange = 0;
  *topolchange = 0;
  // Atermsum only used for Stepwise mutation model
  Atermsum = 0;

/* pick an edge, identify freedup edge (the down edge) and the sister edge */
  do
  {
    edge = randposint (L[li].numlines);
  } while (gtree[edge].down == -1);
  freededge = gtree[edge].down;
  if ((oldsis = gtree[freededge].up[0]) == edge)
    oldsis = gtree[freededge].up[1];

  /* copy information on the edge,  and if it connects to the root, then the sister edge as well */
  if (gtree[edge].down == G->root)
  {
    fillmiginfo (ci, li, gtree, edge, oldsis);
  }
  else
  {
    fillmiginfo (ci, li, gtree, edge, -1);
  }

  /* store information on the genealogy before changing it */
  storeoldedges (ci, li, edge, oldsis, freededge);
// remove any migrations  from the slidingedge 
  gtree[edge].mig[0].mt = -1;
/* slide edge, pick a distance and slide it  */
  
/* PROBLEM it is possible for this to generate too big a slide distance if the sample size is large and there are a lot of very short edges 
  this causes the recursion in slider to crash*/ 
  slidestdv = DMIN (SLIDESTDVMAX, G->roottime/3 );
  holdslidedist = slidedist = normdev (0.0, slidestdv); 

// join the sister and the down branches at the point where edge used to connect, this frees up the down branch 
  joinsisdown (ci, li, oldsis, tmrcachange);

                            // use when debugging with gtreeprint()
  //gtree[edge].down = -1;  // not necessary but makes gtreeprint() output easier to read for intermediate stages
  //gtree[freededge].down = -1; // not necessary but makes gtreeprint() output easier to read for intermediate stages

// do the slide and identify the new sister branch and where new connection point for the edge is 
  newsis = oldsis;
  if (modeloptions[NOMIGRATION])
    slider_nomigration (ci, li, edge, &newsis, &(gtree[edge].time),&slidedist);
  else
    slider (ci, li, edge, &newsis, &(gtree[edge].time), &slidedist);
  *topolchange += (oldsis != newsis);

// now separate the new sister branch into a shorter sis branch and a down branch 
  splitsisdown (ci, li, edge, freededge, newsis);

  if (rootmove)
  {
    slideweight = -log (normprob (0.0, slidestdv, holdslidedist));
    slidestdv = DMIN (SLIDESTDVMAX, G->roottime / 3); 
    slideweight += log (normprob (0.0, slidestdv, holdslidedist)); 
  }
  else
  {
    slideweight = 0;
  }

// add migration events 
  if (modeloptions[NOMIGRATION])
  {
    migweight = 0;
  }
  else
  {
    migweight = addmigration (ci, li);
  }

  // copy the migration info in newedgemig and newsismig  to the genealogy
  copynewmig_to_gtree (ci, li);

// determine all the weights needed for calculating the probability of the genealogy
  setzero_genealogy_weights (&G->gweight);

  treeweight (ci, li);
  sum_subtract_treeinfo (&C[ci]->allgweight, &G->gweight,
                         &holdgweight_updategenealogy);

/* calculate P(D|G)  for new genealogy */
  rejectIS = 0;                 /* use this to catech when P(D|G) for IS model is zero */
  newpdg = 0;

  switch (L[li].model)
  {
  case HKY:
    if (assignmentoptions[JCMODEL] == 1)
    {
      newpdg_a[0] = likelihoodJC (ci, li, G->uvals[0]);
      newpdg = newpdg_a[0];
    }
    else
    {
      newpdg = newpdg_a[0] =
        likelihoodHKY (ci, li, G->uvals[0], G->kappaval, edge,
                       freededge, oldsis, newsis);
    }
    break;
  case INFINITESITES:
    newpdg = newpdg_a[0] = like = likelihoodIS (ci, li, G->uvals[0]);
    rejectIS = (like == REJECTINFINITESITESCONSTANT);
    break;
  case STEPWISE:
    {
      for (ai = 0, newpdg = 0; ai < L[li].nlinked; ai++)
      {
        newpdg_a[ai] =
          G->pdg_a[ai] + finishSWupdateA (ci, li, ai, edge, freededge,
                                          oldsis, newsis,
                                          G->uvals[ai], &Aterm[ai]);

        newpdg += newpdg_a[ai];
        Atermsum += Aterm[ai];
      }

//            checklikelihoodSW(ci, li,G->u[ai].mcinf.val);  
      break;
    }
  case JOINT_IS_SW:
    newpdg = newpdg_a[0] = likelihoodIS (ci, li, G->uvals[0]);
    rejectIS = (newpdg == REJECTINFINITESITESCONSTANT);
    for (ai = 1; ai < L[li].nlinked; ai++)
    {
      newpdg_a[ai] =
        G->pdg_a[ai] + finishSWupdateA (ci, li, ai, edge, freededge,
                                        oldsis, newsis,
                                        G->uvals[ai], &Aterm[ai]);
      newpdg += newpdg_a[ai];
      Atermsum += Aterm[ai];
    }

    //checklikelihoodSW(ci, li,Q[ci]->us[li]);  
    break;
  }
  accp = 0;

/* final weight calculation */
/* tpw is the ratio of new and old prior probability of the genealogies.  It is actually the ratio of the total across all loci,  but
since only genealogy li is being changed at the present time,  the ratio works out to just be the ratio for genealogy li */
  copy_probcalc (&holdallpcalc_updategenealogy, &C[ci]->allpcalc);
  /* Find all internal node sequences and mutations of a full genealogy. */
  if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
  {
    if (L[li].model == INFINITESITES)
    {
      accp = IMA_genealogy_findIntSeq (ci, li);
      if (accp == 0)
      {
        rejectIS = 1;
      }
      accp = 0;
    }
  }
  if (rejectIS == 0)
  {
    // the metropolis term includes p(D|G) and p(G),  
    tpw = -C[ci]->allpcalc.probg;
    integrate_tree_prob (ci, &C[ci]->allgweight,
                         &holdallgweight_updategenealogy, &C[ci]->allpcalc,
                         &holdallpcalc_updategenealogy, &holdt[0]);
    tpw += C[ci]->allpcalc.probg;
    /* 5/19/2011 JH adding thermodynamic integration */
    if (calcoptions[CALCMARGINALLIKELIHOOD])
    {
      metropolishastingsterm = exp (beta[ci] * gbeta*(newpdg - G->pdg) + tpw + migweight + slideweight + Atermsum);
    }
    else
    {
    metropolishastingsterm = tpw + gbeta*(newpdg - G->pdg);
      metropolishastingsterm = exp (beta[ci] * metropolishastingsterm + migweight + slideweight + Atermsum);
    }
    U = uniform ();
    if (U < DMIN(1.0, metropolishastingsterm))  //9/13/2010 
    {
      /* accept the update */
      C[ci]->allpcalc.pdg -= G->pdg;
      C[ci]->allpcalc.pdg += newpdg;
      G->pdg = newpdg;
      for (ai = 0; ai < L[li].nlinked; ai++)
        G->pdg_a[ai] = newpdg_a[ai];
      if (L[li].model == HKY)
      {
        copyfraclike (ci, li);
        storescalefactors (ci, li);
      }
      accp = 1;
    }
  }
  /* reject the update */
  if (accp == 0)
  {
    // put the edges back 
    restoreedges (ci, li, edge, oldsis, freededge, newsis);

    // copy summary stats back
    storegenealogystats (ci, li, 1);

    // reset HKY terms
    if (L[li].model == HKY)
      restorescalefactors (ci, li);
    // copy back all the weights and results associated with calculating the probability of the genealogy 
    copy_probcalc (&C[ci]->allpcalc, &holdallpcalc_updategenealogy);
    copy_treeinfo (&C[ci]->allgweight, &holdallgweight_updategenealogy);
    copy_treeinfo (&G->gweight, &holdgweight_updategenealogy);
    *topolchange = 0;
    *tmrcachange = 0;
  }
  
/* do updates at nodes for stepwise loci, regardless of whether slide update was accepted.  This could go somewhere else  */
  
  return accp;
}   /*updategenealogy */



#undef SLIDESTDVMAX
