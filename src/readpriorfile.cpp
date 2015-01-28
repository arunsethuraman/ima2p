/*IMa2p 2009-2015 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, and Arun Sethuraman */
#undef GLOBVARS
#include "imamp.hpp"

/* 10/3/09  work in progress
  adding code for reading a file that contains priors of parameters */ 


/* format of the prior file 

a title line 
then one or more lines beginning with # 

then popstring
then popstring with popsize priors
then popstring with t priors
then migration prior matrix 


*/

static int ancestralpopnums[2*MAXPOPS - 1]; 
static double *pp;
static char *treestringspot;
static int numpopnodes;
static struct popedge *temppoptree;
static void readprior_parenth0 (void);
static void readprior_parenth (int mode, int tempcurrent, int startparenth);
static void readprior_poptreeread (int mode, char *poptreestring);

#define MAXPRIORTEXTLINE 500

/* readprior_parenth0()
come in on first opening parenthesis
the ith opening parenthesis is associated with a number after its corresponding closing parenthesis

readprior_parenth0() will go through parentheses from left to right,
when it comes to a close parenthesis it records the ancestral node number for that pairing

ancestral node numbers range from npops to 2*npops-2 and proceed from lowest
to highest in order of what time they occured

the ancestral node number is recorded in the array ancestralpopnums[]

the position in the array is the count of which parenthesis pair has just closed, plus npops

in other words if it is the 0'th parenthesis pair (i.e. the first one that opened, meaning it is the
outermost pair),  then the ancestral node number is recorded in ancestralpopnums[npops]

If it is the ith pair that has closed, it is recorded in ancestralpopnums[npops + i]

This seemed to be the only way to get the correct labeling of internal nodes. 

Then when the function readprior_parenth() is called,  the correct times and populatino sizes can be associated 
with these ancestral populations. 
*/

void
readprior_parenth0 (void)
{
  int itemp;
  char *ne;
  int psetlist[MAXPOPS], nextlistspot, popennum;
  nextlistspot = 0;
  popennum = 0;
  ne = treestringspot;
  while (*ne != '\0')
  {
    if (*ne == '(')
    {
      psetlist[nextlistspot] = popennum;
      nextlistspot++;
      popennum++;
      ne++;
    }
    else
    {
      if (*ne == ')')
      {
        ne += 2;
        itemp = strtol (ne, &ne, 10);
        ancestralpopnums[npops + psetlist[nextlistspot - 1]] = itemp;
        nextlistspot--;
      }
      else
      {
        ne++;
      }
    }
  }
}                               /* parenth0 */


void
readprior_parenth (int mode, int tempcurrent, int startparenth)
/* recursive - reenters for every open parentheses */
{
  int itemp, current;
  static int nextnode;
  static int periodi;
  double val;
  char *ne;

  if (startparenth == 1)
  {
    nextnode = -1;
    periodi = 0;
  }
  current = ancestralpopnums[tempcurrent];

  treestringspot++;
  while (isspace (*(treestringspot + 1))) // why  + 1 ? 
    treestringspot++;

  /* next could be:
     - a number  (a simple upnode) read it and get ':' and float number after it
     - an open parenthesis (a complex upnode - call parenth)
     - a comma  skip it
     - a close parenthesis - close the node 
   */
  do
  {
    if (isdigit (*treestringspot))
    {
      /*itemp = atoi (treestringspot);
      treestringspot++; */
      itemp = strtol(treestringspot,&ne,10);							// read the id of population
      treestringspot=ne;
      treestringspot++; /* skip colon */
      val = strtod(treestringspot,&ne);		 // read time or size of population
      treestringspot=ne;
      if (mode ==1)
      {
        //assert(current == itemp);
        temppoptree[itemp].time = val;
        temppoptree[current].up[temppoptree[current].numup] = itemp;
        temppoptree[current].numup++;
        temppoptree[itemp].down = current;
      }
      if (mode==2) 
      {
        pp[itemp] = val;
      }
    }
    if (*treestringspot == ',')
      treestringspot++;
    if (*treestringspot == '(')
    {
      if (nextnode == -1)
        nextnode = npops + 1;
      else
        nextnode++;
      if (mode==1)
      {
        temppoptree[ancestralpopnums[nextnode]].down = current;
        temppoptree[current].up[temppoptree[current].numup] =  ancestralpopnums[nextnode];
        temppoptree[current].numup++;
      }
      readprior_parenth (mode, nextnode, 0);
    }
  } while (*treestringspot != ')');
  treestringspot++;             /* skip parentheses and colon*/
  if (*treestringspot == ':')
  {
    treestringspot++;
    itemp = strtol(treestringspot,&ne,10); // read the id of population
    treestringspot=ne;
    treestringspot++; /* skip colon */
    val = strtod(treestringspot,&ne);		 // read time or size of population
    treestringspot=ne;
    if (itemp < npops)
      IM_err (IMERR_POPTREESTRINGFAIL,
              " wrong number of ancestral populations indicated. string %s, step %d",
              treestringspot, step);
    periodi = itemp - npops;
    if (mode==1)
    {
      assert(current == itemp);
      temppoptree[current].time = val;
      temppoptree[current].b = periodi + 1;
      temppoptree[temppoptree[current].up[0]].e =
      temppoptree[temppoptree[current].up[1]].e = periodi + 1;
    }
    if (mode  == 2)
    {
      pp[itemp] = val;
    }

  }
  else
  { // is it possible to get here ? 
    temppoptree[current].b = periodi + 1;
    temppoptree[temppoptree[current].up[0]].e =
      temppoptree[temppoptree[current].up[1]].e = periodi + 1;
    periodi++;
  }
  if (temppoptree[current].down != -1)
  {
    numpopnodes++;
    current = temppoptree[current].down;
  }
  else
  {
    periodi++;
    temppoptree[current].e = -1;
  }
}                               /* readprior_parenth */

void
readprior_poptreeread (int mode, char *poptreestring)
{
  int i, j;

  /* read in the tree string until enough parentheses are found */
  /* pcount counts parentheses '(' is +1 ')' is -1 repeat until 0 */
  treestringspot = poptreestring;
  if (mode==0)
  {
    numpopnodes = 0;
    for (i = 0; i < npops; i++)
    {
      temppoptree[i].b = 0;
      temppoptree[i].numup = 0;
      temppoptree[i].up = static_cast<int *> (malloc (2 * sizeof (int)));
      for (j = 0; j < 2; j++)
        temppoptree[i].up[j] = -1;
      temppoptree[i].down = -1;
    }
    for (; i < numtreepops; i++)
    {
      temppoptree[i].numup = 0;
      temppoptree[i].up = static_cast<int *> (malloc (2 * sizeof (int)));
      for (j = 0; j < 2; j++)
        temppoptree[i].up[j] = -1;
      temppoptree[i].down = -1;
    }
    temppoptree[npops].down = -1;
    readprior_parenth0 ();
  }
  else
    readprior_parenth (mode, npops, 1);
  return;
}  /* end readprior_poptreeread */


void readpriorfile(char priorfilename[],double *popsizepriorvals, double **mpriorvals)
{
  FILE *priorfile; 
  char  *chpt;
  int i,j;
  char *priortextline;
  char *treetext;
  if ((priorfile = fopen (priorfilename, "r")) == NULL)
  {
    IM_err(IMERR_READFILEOPENFAIL,"Error opening file with prior values: %s", priorfilename);
  }
  pp = popsizepriorvals;
  priortextline = static_cast<char *> (malloc(MAXPRIORTEXTLINE*sizeof(char)));
  temppoptree = static_cast<popedge *> (malloc (numtreepops * sizeof (struct popedge)));
  for (i=0;i<3;i++)
  {
    while (fgets(priortextline,MAXPRIORTEXTLINE,priorfile)!= NULL && priortextline[0]=='#') {};
    treetext = &(priortextline[0]);
    readprior_poptreeread(i,treetext);
  }
  for (i=npops;i<numtreepops;i++)
  {
    tperiodpriors[i-npops]=temppoptree[temppoptree[i].up[0]].time;
    if (temppoptree[temppoptree[i].up[0]].time != temppoptree[temppoptree[i].up[1]].time)
      IM_err(IMERR_PRIORFILEVALS,"two split times specified in tree string in prior file are not equal: %lf, %lf",temppoptree[temppoptree[i].up[0]].time,temppoptree[temppoptree[i].up[1]].time);
  }
  for (i=0;i<numsplittimes-1;i++)
    if (tperiodpriors[i] > tperiodpriors[i+1])
      IM_err(IMERR_PRIORFILEVALS,"earlier max split time greater than max for older split:%lf, %lf", tperiodpriors[i],tperiodpriors[i+1]);
  // read in migration rates
  max_m_from_priorfile = -1.0;
  for (i=0;i<numtreepops;i++)									
	  {
      while (fgets(priortextline,MAXPRIORTEXTLINE,priorfile)!= NULL && priortextline[0]=='#') {};
	  chpt = &priortextline[0];
	  for (j=0;j<numtreepops;j++)
      {
	    mpriorvals[i][j] = strtod(chpt,&chpt);
        if (mpriorvals[i][j] > max_m_from_priorfile )
          max_m_from_priorfile  = mpriorvals[i][j];
        if (i>j && ( (mpriorvals[i][j]==0.0 && mpriorvals[j][i] > 0.0) ||(mpriorvals[i][j] > 0.0 && mpriorvals[j][i] == 0.0)))
          IM_err(IMERR_PRIORFILEVALS,"reciprocal migration rates both not zero or both not greater than zero: %d->%d %lf; %d->%d, %lf", i,j,mpriorvals[i][j],j,i,mpriorvals[j][i]);
       }
	  }

  for (i=0; i < numtreepops; i++)
    XFREE(temppoptree[i].up);
  XFREE(temppoptree);
  XFREE(priortextline);
  fclose(priorfile);
} // readpriorfile
   

