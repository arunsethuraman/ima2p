/*IMa2p 2009-2015 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, and Arun Sethuraman */

#ifdef _DEBUG

/* used only when debugging  and when the gtreeprint(ci, li) is called 
generates a file called gtreeprint.out - can be very large if gtreeprint() is in a loop - careful */

#undef GLOBVARS
#include "imamp.hpp"

/*********** LOCAL STUFF **********/

/* treevals only used in gtreeprint() and poptreeprint() for debugging */
struct treevals
{
  int nodenum;
  int c_or_m;
  int up1;
  int up2;
  int down;
  double utime;
  double dtime;
  double timei;
  int migcount;
  int startpop;
  int finpop;
  int mut;
  int A[MAXLINKED];
  double dlikeA[MAXLINKED];
  int flag;
} *savetree[MAXLOCI];

/* local function prototypes */
static void shellgtreevals (struct treevals *lptr, int length);

/* local functions */

void
shellgtreevals (struct treevals *lptr, int length)
{
  double aln = 1.442695022, tiny = 1.0e-5;
  struct treevals t;
  static int nn, m, lognb2, i, j, k, l;
  lognb2 = (int) floor (log (static_cast<double>(length)) * aln + tiny);
  m = length;
  for (nn = 1; nn <= lognb2; nn++)
  {
    m = m / 2;
    k = length - m;
    for (j = 0; j <= k - 1; j++)
    {
      i = j;
    reloop:l = i + m;
      if ((lptr + l)->dtime < (lptr + i)->dtime)
      {
        t = *(lptr + i);
        *(lptr + i) = *(lptr + l);
        *(lptr + l) = t;
        i = i - m;
        if (i >= 0)
          goto reloop;
      }
    }
  }
}                               /* shellfithet */


/********** GLOBAL FUNCTIONS ***********/

/* callsource not that useful,  commented out on 4/1/2008 */
void
gtreeprint (int ci, int li/*, int step , int callsource */ )
/* use this in debugging mode to see what a genealogy looks like for a particular parameter set */
{
  int p, i, j, k, site, node, up1, up2, nm, totalm;
  int startpop, periodi;
  double upt;
  struct treevals *tvcptr;
  struct treevals *tvmptr;
  char sc[MAXPERIODS][5];
  int checkpt[MAXPERIODS];
  FILE *treeprintfile;
  char treefilename[] = "gtreeprint.out";
  char sourcestring[100];
  struct genealogy *G = &(C[ci]->G[li]);
  struct edge *gtree = G->gtree;
  int fromi, toi;

/*	switch (callsource)
		{
		case 8 : strcpy(sourcestring,"make_JOINT_IS_SW()"); break;
		case 7 : strcpy(sourcestring,"makeSW()"); break;
		case 6 : strcpy(sourcestring,"setup_chains()"); break;
		case 5 : strcpy(sourcestring,"addmigration()"); break;
		case 3 : strcpy(sourcestring,"makeIS()"); break;
		case 4 : strcpy(sourcestring,"main()"); break;
		case 0 : strcpy(sourcestring,"updategenealogy()"); break;
		case 1 : strcpy(sourcestring,"treeweight()"); break;
		case 2 : strcpy(sourcestring,"integrate_tree_prob()"); break;
		default : strcpy(sourcestring,"UNKNOWN");
		} */
  for (i = 0; i < lastperiodnumber; i++)
  {
    sprintf (sc[i], "t%d>", i);
    checkpt[i] = 1;
  }
  totalm = 0;
  for (i = 0; i < 2 * L[li].numgenes - 1; i++)
  {
    j = 0;
    while (gtree[i].mig[j].mt > -0.5)
      j++;
    totalm += j;
  }
  tvcptr = static_cast<treevals *>
	  (malloc ((2 * L[li].numgenes - 1) * (sizeof (struct treevals))));
  tvmptr = static_cast<treevals *>
	  (malloc (totalm * sizeof (struct treevals)));
  for (i = 0; i < 2 * L[li].numgenes - 1; i++)
    tvcptr[i].mut = 0;
  for (site = 0; site < L[li].numsites; site++)
  {
    for (j = 0; j < L[li].numgenes; j++)
      gtree[j].mut = L[li].seq[j][site];
    for (j = L[li].numgenes; j < L[li].numlines; j++)
      gtree[j].mut = -1;
    for (j = 0; j < L[li].numgenes; j++)
      labelgtree (ci, li, j);
    for (node = L[li].numgenes; node < L[li].numlines; node++)
    {
      i = node - L[li].numgenes;
      up1 = gtree[node].up[0];
      up2 = gtree[node].up[1];
      if (gtree[node].down != -1)
        if (i != gtree[gtree[node].down].up[0]
            && i != gtree[gtree[node].down].up[1])
          tvcptr[node].flag = 2;
      if (gtree[node].up[0] != -1)
        if (i != gtree[gtree[node].up[0]].down)
          tvcptr[node].flag = 3;
      if (gtree[node].up[1] != -1)
        if (i != gtree[gtree[node].up[1]].down)
          tvcptr[node].flag = 4;
      if ((gtree[up1].mut == 0
           && gtree[up2].mut == 1)
          || (gtree[up1].mut == 1 && gtree[up2].mut == 0))
      {
        if (gtree[node].down == -1)     /* root - not clear where mutation is , put it on left */
        {
          tvcptr[up1].mut++;
        }
        else
        {
          if (gtree[gtree[node].down].mut == gtree[up1].mut)
          {
            tvcptr[up2].mut++;
          }
          else
          {
            tvcptr[up1].mut++;
          }
        }
      }
    }
  }

  for (i = 0, k = 0; i < 2 * L[li].numgenes - 1; i++)
  {
    tvcptr[i].nodenum = i;
    tvcptr[i].c_or_m = 0;
    tvcptr[i].flag = 0;
    tvcptr[i].up1 = gtree[i].up[0];
    tvcptr[i].up2 = gtree[i].up[1];
    tvcptr[i].down = gtree[i].down;
    tvcptr[i].dtime = gtree[i].time;
    if (i >= L[li].numgenes)
      upt = gtree[gtree[i].up[0]].time;
    else
      upt = 0;
    tvcptr[i].utime = upt;
    tvcptr[i].timei = gtree[i].time - upt;
    tvcptr[i].startpop = gtree[i].pop;
    j = 0;
    while (gtree[i].mig[j].mt > -1)
      j++;
    tvcptr[i].migcount = j;
    if (j > 0)
      tvcptr[i].finpop = gtree[i].mig[j - 1].mp;
    else
      tvcptr[i].finpop = tvcptr[i].startpop;
    while (tvcptr[i].dtime > C[ci]->poptree[tvcptr[i].finpop].time)
    {
      tvcptr[i].finpop = C[ci]->poptree[tvcptr[i].finpop].down;
    }
    if (tvcptr[i].finpop >= numtreepops)
      tvcptr[i].flag = 5;
    if (L[li].model == STEPWISE)
    {
      tvcptr[i].A[0] = gtree[i].A[0];
      tvcptr[i].dlikeA[0] = gtree[i].dlikeA[0];
    }
    if (L[li].model == JOINT_IS_SW)
    {
      tvcptr[i].A[0] = gtree[i].A[1];
      tvcptr[i].dlikeA[0] = gtree[i].dlikeA[1];
    }
    nm = j;
    if (nm > 0)
      for (startpop = gtree[i].pop, j = 0; j < nm; j++)
      {
        tvmptr[k].nodenum = i;
        tvmptr[k].c_or_m = 1;
        tvmptr[k].dtime = gtree[i].mig[j].mt;
        periodi = findperiod (ci, tvmptr[k].dtime);

        /* could change population by passing into the next period */
        while (C[ci]->poptree[startpop].e <= periodi
               && C[ci]->poptree[startpop].e != -1)
          startpop = C[ci]->poptree[startpop].down;
        tvmptr[k].startpop = startpop;
        tvmptr[k].finpop = gtree[i].mig[j].mp;
        startpop = tvmptr[k].finpop;
        k++;
      }
  }
  shellgtreevals (tvcptr, L[li].numgenes);
  shellgtreevals (tvcptr + L[li].numgenes, L[li].numgenes - 1);
  if (!(treeprintfile = fopen (treefilename, "a")))
  {
    IM_err(IMERR_APPENDFILEFAIL,"Error opening tree file for appending");
  }
  fprintf (treeprintfile,
           "chain: %d Locus: %d  Step: %d  Calling function: %s\n", ci, li,
           step, sourcestring);
  fprintf (treeprintfile, "   split times:");
  for (i = 0; i < lastperiodnumber; i++)
    fprintf (treeprintfile, "  %8.4f", C[ci]->tvals[i]);
  fprintf (treeprintfile, "\n");
  fprintf (treeprintfile, "  current p(D|G): %.5f current p(G) %.5f\n", C[ci]->G[li].pdg , C[ci]->allpcalc.pdg);
  fprintf (treeprintfile,
           "\tNode#\tup1\tup2\tdown\tdtime\tutime\ttimei\tpop\tfinpop\t#mig\tflag");
  if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
    fprintf (treeprintfile, "\tA\tld\n");
  else
    fprintf (treeprintfile, "\n");
  p = 0;
  for (i = 0; i < 2 * L[li].numgenes - 1; i++)
  {
    while (tvcptr[i].dtime > C[ci]->tvals[p] && checkpt[p])

    {
      fprintf (treeprintfile,
               "%s========================================================================================\n",
               sc[p]);
      checkpt[p] = 0;
      p++;
    }
    //  assert(tvcptr[i].down == -1 || gtree[tvcptr[i].down].pop == tvcptr[i].finpop);
    fprintf (treeprintfile,
             "\t%3d\t%3d\t%3d\t%3d\t%7.4f\t%7.4f\t%7.4f\t%3d\t%3d\t%d\t%d",
             tvcptr[i].nodenum, tvcptr[i].up1, tvcptr[i].up2,
             tvcptr[i].down, tvcptr[i].dtime, tvcptr[i].utime,
             tvcptr[i].timei, tvcptr[i].startpop, tvcptr[i].finpop,
             tvcptr[i].migcount, tvcptr[i].flag);
    if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
      fprintf (treeprintfile, "\t%3d\t%7.5f\n", tvcptr[i].A[0],
               tvcptr[i].dlikeA[0]);
    else
      fprintf (treeprintfile, "\n");
    for (j = 0; j < totalm; j++)
    {
      if (tvmptr[j].nodenum == tvcptr[i].nodenum)
        fprintf (treeprintfile, "\t%3d\t\t\t\t%7.4f\t\t\t%3d\t%3d\n",
                 tvmptr[j].nodenum, tvmptr[j].dtime,
                 tvmptr[j].startpop, tvmptr[j].finpop);
    }
  }
  fprintf (treeprintfile, "\n");
  fprintf (treeprintfile, "Migration counts :\n");
  for (i = 0; i < nummigrateparams; i++)
  {
    fprintf (treeprintfile, "  %s\t", imig[i].str);
    fromi = atoi (&imig[i].str[1]);
    toi = atoi (&imig[i].str[3]);
    for (j = 0, k = 0; j < totalm; j++)
    {
      k += (fromi == tvmptr[j].startpop && toi == tvmptr[j].finpop);
    }
    fprintf (treeprintfile, "%d\n", k);
  }
  fprintf (treeprintfile, "\n");
  f_close (treeprintfile);
  treeprintfile = NULL;
  XFREE (tvcptr);
  XFREE (tvmptr);
}                               /* gtreeprint */

void
poptreeprint (int ci/*, int step*/)
/* use this in debugging mode to see what a genealogy looks like for a particular parameter set */
{
  int i;
  double upt;
  struct treevals *tvptr;
  FILE *treeprintfile;
  char treefilename[] = "poptreeprint.out";
  tvptr = static_cast<treevals *>
	  (malloc ((2 * npops - 1) * (sizeof (struct treevals))));
  for (i = 0; i < 2 * npops - 1; i++)
  {
    tvptr[i].nodenum = i;
    tvptr[i].up1 = C[ci]->poptree[i].up[0];
    tvptr[i].up2 = C[ci]->poptree[i].up[1];
    tvptr[i].down = C[ci]->poptree[i].down;
    tvptr[i].dtime = C[ci]->poptree[i].time;
    if (i >= npops)
      upt = C[ci]->poptree[C[ci]->poptree[i].up[0]].time;
    else
      upt = 0;
    tvptr[i].utime = upt;
    tvptr[i].timei = C[ci]->poptree[i].time - upt;
  }
  shellgtreevals (tvptr, 2 * npops - 1);
  if (!(treeprintfile = fopen (treefilename, "a")))
  {
    IM_err(IMERR_APPENDFILEFAIL,"Error opening tree file for appending");
  }
  fprintf (treeprintfile, "Chain: %d  Step: %d\n", ci, step);
  fprintf (treeprintfile, "   split times:");
  for (i = 0; i < lastperiodnumber; i++)
    fprintf (treeprintfile, "  %8.4f", C[ci]->tvals[i]);
  fprintf (treeprintfile, "\n");
  fprintf (treeprintfile, "Node#\tup1\tup2\tdown\tdtime\tutime\ttimei\n");
  for (i = 0; i < 2 * npops - 1; i++)
  {
    fprintf (treeprintfile, "\t%3d\t%3d\t%3d\t%3d\t%7.4f\t%7.4f\t%7.4f\n",
             tvptr[i].nodenum, tvptr[i].up1, tvptr[i].up2,
             tvptr[i].down, tvptr[i].dtime, tvptr[i].utime, tvptr[i].timei);
  }
  fprintf (treeprintfile, "\n");
  f_close (treeprintfile);
  treeprintfile = NULL;
  XFREE (tvptr);
}                               /* poptreeprint */

#endif
