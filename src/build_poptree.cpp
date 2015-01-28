/*IMa2p 2009-2015 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, and Arun Sethuraman */
#undef GLOBVARS
#include "imamp.hpp"

/*********** LOCAL STUFF **********/
static int numpopnodes;
static char *treestringspot;
static int pos;
static int ancestralpopnums[2 * MAXPOPS - 1];
static void parenth0 (int current);
static void parenth (int ci, int tempcurrent, int startparenth);
static void fillplist (int ci);
static void simpoptree (int ci);
static void poptreeread (int ci, char *poptreestring);
static void poptreewrite (int ci, char *buildstr);

/* popstring primary format:
for npops,  the populations are numbered 0 thru npops-1
format includes branching pattern and node order in time
every node is flanked by parentheses, and every pair of items within a node is separated by a comma
every closed parentheses is followed by a colon and a node sequence number
node sequence numbers are ordered from most recent (npops) oldest,  2*(npops - 1)

In other words, the external nodes (sampled populations) are numbered
from 0 to npops - 1,  and the internal nodes are numbered from npops to 2*(npops - 1)
*/

/* parenth0() identifies the order of nodes as they will be reached by parenth()
then it associates these with the correct period and ancestral population size numbers
and puts these into ancestralpopnums[] to be used by parenth()  


How it works:
come in on first opening parenthesis
the ith opening parenthesis is associated with a number after its corresponding closing parenthesis

parenth0() will go through parentheses from left to right,
when it comes to a close parenthesis it records the ancestral node number for that pairing

ancestral node numbers range from numpops to 2*numpops-2 and proceed from lowest
to highest in order of what time they occured

the ancestral node number is recorded in the array ancestralpopnums[]

the position in the array is the count of which parenthesis pair has just closed, plus numpops

in other words if it is the 0'th parenthesis pair (i.e. the first one that opened, meaning it is the
outermost pair),  then the ancestral node number is recorded in ancestralpopnums[numpops]

If it is the ith pair that has closed, it is recorded in ancestralpopnums[numpops + i]

This seemed to be the only way to get the correct labeling of internal nodes. 

Then when the function parenth() is called,  the correct times and populatino sizes can be associated 
with these ancestral populations. 
*/

/**** LOCAL FUNCTIONS *****/

void
parenth0 (int current)
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

/* parenth()  reads the tree topology and puts it into poptree structure */

void
parenth (int ci, int tempcurrent, int startparenth)
/* recursive - reenters for every open parentheses */
{
  int itemp, current, i;
  static int nextnode;
  static int periodi;

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
      itemp = atoi (treestringspot);
      treestringspot++;
      C[ci]->poptree[current].up[C[ci]->poptree[current].numup] = itemp;
      C[ci]->poptree[current].numup++;
      C[ci]->poptree[current].time = 0;
      C[ci]->poptree[itemp].down = current;
    }
    if (*treestringspot == ',')
      treestringspot++;
    if (*treestringspot == '(')
    {
      if (nextnode == -1)
        nextnode = npops + 1;
      else
        nextnode++;
      C[ci]->poptree[ancestralpopnums[nextnode]].down = current;
      C[ci]->poptree[current].up[C[ci]->poptree[current].numup] =
        ancestralpopnums[nextnode];
      C[ci]->poptree[current].numup++;
      C[ci]->poptree[current].time = 0;
      parenth (ci, nextnode, 0);
    }
  } while (*treestringspot != ')');
  treestringspot++;             /* skip parentheses */
  if (*treestringspot == ':')
  {
    treestringspot++;
    i = atoi (treestringspot);
    if (i < npops)
      IM_err (IMERR_POPTREESTRINGFAIL,
              " wrong number of ancestral populations indicated. string %s, step %d",
              treestringspot, step);
    assert (i >= npops);
    periodi = i - npops;
    C[ci]->poptree[current].b = periodi + 1;
    C[ci]->poptree[C[ci]->poptree[current].up[0]].e =
      C[ci]->poptree[C[ci]->poptree[current].up[1]].e = periodi + 1;
    if (i >= 10)
      treestringspot += 2;
    else
      treestringspot++;
  }
  else
  { // is it possible to get here? 
    C[ci]->poptree[current].b = periodi + 1;
    C[ci]->poptree[C[ci]->poptree[current].up[0]].e =
      C[ci]->poptree[C[ci]->poptree[current].up[1]].e = periodi + 1;
    periodi++;
  }
  if (C[ci]->poptree[current].down != -1)
  {
    numpopnodes++;
    current = C[ci]->poptree[current].down;
  }
  else
  {
    periodi++;
    C[ci]->poptree[current].e = -1;
    C[ci]->rootpop = current;
  }
}                               /* parenth */


/* find root,  move up, recurse, build the string, makes strings 
 * that contain node order information 
 *
 * This function is used for debugging only.  CR 110825.1  
 */
void
makepoptreestring (int ci, int curpop, char *buildstr)
{
  int i;
  char ss[5];
  if (curpop == -1)
  {
    i = 0;
    while (C[ci]->poptree[i].down != -1)
      i++;
    curpop = i;
    pos = 0;
  }
  ss[0] = '\0';
  sprintf (ss, "(\0");
  strinsert (buildstr, ss, pos);
  pos++;
  for (i = 0; i < C[ci]->poptree[curpop].numup; i++)
  {
    if (C[ci]->poptree[curpop].up[i] < npops)
    {
      ss[0] = '\0';
      sprintf (ss, "%d\0", C[ci]->poptree[curpop].up[i]);
      strinsert (buildstr, ss, pos);
      pos++;
    }
    else
    {
      makepoptreestring (ci, C[ci]->poptree[curpop].up[i],
                         buildstr /*, pos */ );
    }
    if (i < C[ci]->poptree[curpop].numup - 1)
    {
      ss[0] = '\0';
      sprintf (ss, ",\0");
      strinsert (buildstr, ss, pos);
      pos++;
    }
    else
    {
      ss[0] = '\0';
      sprintf (ss, "):%d\0", C[ci]->poptree[C[ci]->poptree[curpop].up[0]].b);
      strinsert (buildstr, ss, pos);
      pos += 3;
    }
  }
}                               /* makepoptreestring */


/* rewrite() rewrites the treestring in a standard order
    swivels nodes,  if both have node sequence values, the one with the lower node sequence value (periodi[]) goes on the left
    if only one has a node sequence value,  it goes on the right
	when neither has a node sequence value, the one with the lowest node number go on the left */
/* this is simply sorting for a pair.  To handle multifurcations, must put in proper sorting */
/* it actually works,  checked on simulated random trees on 4/17/07 */
/*rewrite if for trees that have node sequence values, it is based on an older version rewrite for trees without node sequence values */
/* works recursively */

void rewritecheckchar(char c)
{
  if ((isdigit(c) || c=='(' || c==',' || c==')' || c==':') == 0)
    IM_err (IMERR_POPTREESTRINGFAIL," something wrong in formatting of population string in input file");
}

void rewrite (char *substr)
{
  int slengths[MAXPOPS];
  int pcount, subpos, subcount;
  char holdsubs[MAXPOPS][POPTREESTRINGLENGTHMAX];
  int firstint[MAXPOPS];
  int i, j, k;
  int periodi[MAXPOPS];
  pos = 1;
  subpos = pos;
  subcount = 0;
  pcount = 0;
  slengths[subcount] = 0;

  do
  {
    if (substr[pos] == '(')
      pcount++;
    if (substr[pos] == ')')
      pcount--;
    pos++;
    slengths[subcount]++;
    if (pcount == 0)
    {
      if (slengths[subcount] > 1)
      {
        pos++;
        i = atoi (&substr[pos]);
        periodi[subcount] = i;
        if (i >= 10)
        {
          pos += 2;
          slengths[subcount] += 3;
        }
        else
        {
          pos++;
          slengths[subcount] += 2;
        }
      }
      else
      {
        periodi[subcount] = -1;
      }
      assert (!(pos < subpos));
      strncpy (holdsubs[subcount], &substr[subpos], (size_t) (pos - subpos));
      holdsubs[subcount][slengths[subcount]] = '\0';
      i = 0;
      while (!isdigit (holdsubs[subcount][i]))
        i++;
      firstint[subcount] = atoi (&holdsubs[subcount][i]);
      subcount++;
      slengths[subcount] = 0;
      if (substr[pos] == ',')
      {
        pos++;
      }
      subpos = pos;
    }
  } while (pos < (int) strlen (substr));
  if ((periodi[0] > periodi[1] && periodi[0] >= 0 && periodi[1] >= 0)
      || (periodi[0] >= 0 && periodi[1] < 0))
  {
    substr[0] = '(';
    j = slengths[1];
    for (i = 1, k = 0; i <= j; i++, k++)
    { 
      rewritecheckchar(holdsubs[1][k]);
      substr[i] = holdsubs[1][k];
    }
    subpos = 1;
    substr[i] = '\0';
    if (slengths[1] > 2)
      rewrite (&substr[subpos]);
    substr[i] = ',';
    i++;
    subpos = i;
    j += 1 + slengths[0];
    for (k = 0; i <= j; i++, k++)
    {
      rewritecheckchar(holdsubs[0][k]);
      substr[i] = holdsubs[0][k];
    }
    substr[i] = '\0';
    if (slengths[0] > 2)
      rewrite (&substr[subpos]);
    substr[i] = ')';
  }
  else
  {
    if (firstint[0] > firstint[1] && periodi[0] < 0 && periodi[1] < 0)
    {
      substr[0] = '(';
      j = slengths[1];
      for (i = 1, k = 0; i <= j; i++, k++)
      {
        rewritecheckchar(holdsubs[1][k]);
        substr[i] = holdsubs[1][k];
      }
      subpos = 1;
      if (slengths[1] > 2)
        rewrite (&substr[subpos]);
      substr[i] = ',';
      i++;
      subpos = i;
      j += 1 + slengths[0];
      for (k = 0; i <= j; i++, k++)
      {
        rewritecheckchar(holdsubs[0][k]);
        substr[i] = holdsubs[0][k];
      }
      if (slengths[0] > 2)
        rewrite (&substr[subpos]);
      substr[i] = ')';
    }
    else
    {
      substr[0] = '(';
      subpos = 1;
      substr[slengths[0] + 1] = '\0';
      if (slengths[0] > 2)
        rewrite (&substr[subpos]);
      substr[slengths[0] + 1] = ',';
      subpos = slengths[0] + 2;
      substr[slengths[0] + slengths[1] + 2] = '\0';
      if (slengths[1] > 2)
        rewrite (&substr[subpos]);
      substr[slengths[0] + slengths[1] + 2] = ')';
    }
  }
} /* rewrite */ ;

void
fillplist (int ci)
{
  int i, j, k;
  SET tempset;

  C[ci]->periodset[0] = EMPTYSET;
  for (i = 0; i < npops; i++)
    C[ci]->periodset[0] = UNION (C[ci]->periodset[0], SINGLESET (i));
  tempset = C[ci]->periodset[0];
  C[ci]->addpop[0] = C[ci]->droppops[0][0] = C[ci]->droppops[0][1] = -1;
  C[ci]->addpop[npops] = C[ci]->droppops[npops][0] = C[ci]->droppops[npops][1] = 0;
  for (i = 1; i < npops; i++)   // loop over periods
  {
    k = 0;
    for (j = 0; j < numtreepops; j++)
    {
      if (C[ci]->poptree[j].e == i)
      {
        if (k >= 2)
          IM_err (IMERR_POPTREESTRINGFAIL,
                  " wrong number ancestral populations indicated. step %d",
                  step);
        assert (ISELEMENT (j, tempset));
        tempset = SETREMOVE (tempset, j);       // remove j from the set
        C[ci]->droppops[i][k] = j;
        k++;
      }
      if (C[ci]->poptree[j].b == i)
      {
        tempset = SETADD (tempset, j);  // add j to the set
        C[ci]->addpop[i] = j;
      }
    }
    C[ci]->periodset[i] = tempset;
  }

  /* fill C[ci]->plist */
  if ((C[ci]->plist = static_cast<int **> (malloc (npops * sizeof (*C[ci]->plist)))) == NULL)
    IM_err (IMERR_MEM, "  plist malloc did not work.   npops %d, step %d",
            npops, step);
  for (i = 0; i < npops; i++)
  {
    if ((C[ci]->plist[i] =
        static_cast<int *> (malloc ((npops - i) * sizeof (*C[ci]->plist[i])))) == NULL)
      IM_err (IMERR_MEM,
              "  plist malloc did not work.   npops - i  %d, step %d",
              npops - i, step);
  }
  for (i = 0; i < npops; i++)
    C[ci]->plist[0][i] = i;
  for (i = 0; i < npops; i++)
  {
    j = 0;
    FORALL (k, C[ci]->periodset[i])
    {
      C[ci]->plist[i][j] = k;
      j++;
    }
  }
}                               /* fillplist */


/*
 * This function is used for debugging only.  CR 110825.1  
 */
void
simpoptree (int ci)
{
  int i, k1, k2, newpop, n, periodi;
  int list[2 * MAXPOPS - 1];
  for (i = 0; i < numtreepops; i++)
    C[ci]->poptree[i].up[0] = C[ci]->poptree[i].up[1] =
      C[ci]->poptree[i].down = -1;
  for (i = 0; i < npops; i++)
  {
    C[ci]->poptree[i].b = 0;
    list[i] = i;
  }
  for (periodi = 1, newpop = npops, n = npops; newpop < numtreepops;
       newpop++, n--, periodi++)
  {
    C[ci]->poptree[newpop].numup = 2;
    k1 = randposint (n);

    do
    {
      k2 = randposint (n);
    } while (k2 != k1);
    C[ci]->poptree[newpop].up[0] = list[k1];
    C[ci]->poptree[newpop].up[1] = list[k2];
    C[ci]->poptree[newpop].b = periodi;
    for (i = k1; i < n; i++)
      list[i] = list[i + 1];
    for (i = k2; i < n - 1; i++)
      list[i] = list[i + 1];
    list[n - 2] = newpop;
    C[ci]->poptree[C[ci]->poptree[newpop].up[0]].time =
      C[ci]->poptree[C[ci]->poptree[newpop].up[1]].time =
      C[ci]->tvals[periodi];
    C[ci]->poptree[C[ci]->poptree[newpop].up[0]].e =
      C[ci]->poptree[C[ci]->poptree[newpop].up[1]].e = periodi;
    C[ci]->poptree[C[ci]->poptree[newpop].up[0]].down =
      C[ci]->poptree[C[ci]->poptree[newpop].up[1]].down = newpop;
  }
  C[ci]->poptree[2 * npops - 2].down = -1;
  C[ci]->poptree[2 * npops - 2].time = TIMEMAX;
  fillplist (ci);
}                               /* simpoptree */

/* SANGCHUL: Wed Apr 15 16:46:44 EDT 2009
 * This is MacOSX specific error.
 * There is something that may be corrected for other debugger.
 * When Sang Chul tried to compile it with memory check option,
 * there is 0 memory access compilation error. This does not seem 
 * to happen when we just compile it without that option.
 */
void
poptreeread (int ci, char *poptreestring)
{
  int i, j;

  /* read in the tree string until enough parentheses are found */
  /* pcount counts parentheses '(' is +1 ')' is -1 repeat until 0 */
  numpopnodes = 0;
  for (i = 0; i < npops; i++)
  {
    C[ci]->poptree[i].b = 0;
    C[ci]->poptree[i].numup = 0;
    C[ci]->poptree[i].up = static_cast<int *> (malloc (2 * sizeof (int)));
    for (j = 0; j < 2; j++)
      C[ci]->poptree[i].up[j] = -1;
    C[ci]->poptree[i].down = -1;
  }
  for (; i < numtreepops; i++)
  {
    C[ci]->poptree[i].numup = 0;
    C[ci]->poptree[i].up = static_cast<int *> (malloc (2 * sizeof (int)));
    for (j = 0; j < 2; j++)
      C[ci]->poptree[i].up[j] = -1;
    C[ci]->poptree[i].down = -1;
  }
  rewrite (poptreestring);
  C[ci]->poptree[npops].down = -1;
  treestringspot = poptreestring;
  if (ci == 0)
    parenth0 (npops);
  parenth (ci, npops, 1);
  return;
}                               /* end treeread */


/*
 * This function is used for debugging only.  CR 110825.1  
 */
void
poptreewrite (int ci, char *buildstr)
{
  makepoptreestring (ci, -1, buildstr);

  //rewrite(buildstr);
  rewrite (buildstr);
}

/********** GLOBAL FUNCTIONS ***********/

/* for checking gtree simulation, writing and sorting 
struct p
		{
		char str[40];
		int  count;
		};
#define simtreenum 20000
void testtreewrite(char  startpoptreestring[])
	{
	FILE *outfile;
	char teststr[40];
	struct p popsimtest[simtreenum];
	int i, j, found, numsofar = 0;
outfile = fopen("testtree.out","w");
	for (i=0; i< simtreenum;i++)
		{
		teststr[0]=0;
		setup_poptree(0,startpoptreestring);
		poptreewrite(0,teststr );
		XFREE(C[0]->poptree);
		j = 0;
		found = 0;
		while (j < numsofar && found == 0)
			{
			if (strcmp(teststr,popsimtest[j].str)==0)
				{
				popsimtest[j].count++;
				found = 1;
				}
			j++;
			}
		if (found == 0)
			{
			strcpy(popsimtest[j].str, teststr);
			popsimtest[j].count = 1;
			numsofar++;
			}
		}
	fprintf(outfile,"numtrees: %d\n",numsofar);
	for (i=0; i< numsofar;i++)
		fprintf(outfile,"%d\t%s\n",popsimtest[i].count,popsimtest[i].str);
	fclose(outfile);
	}
*/
void
add_ghost_to_popstring (char poptreestring[])
{
  size_t i;
  int n;
  char stringBuf[POPTREESTRINGLENGTHMAX]; /* temp buffer to build newstring */
  char newstring[POPTREESTRINGLENGTHMAX];

  strcpy (newstring, "(");
  for (i=0;i<strlen(poptreestring);i++)
  {
    if (poptreestring[i] == ':')
    {
      n = atoi(&poptreestring[i+1]);
      sprintf(stringBuf,"%s%c%d",newstring,':',n+1);
      i++;
      if (n>=10)
        i++;
      strcpy(newstring, stringBuf);
    }
    else
    {
      strncat(newstring,&poptreestring[i],1);
    }
  }
  sprintf (poptreestring, "%s,%d):%d", newstring, npops, 2*npops);
}

void
setup_poptree (int ci, char startpoptreestring[])
{
  int i;
  if (npops == 1 )
  {
    C[ci]->poptree = static_cast<popedge *> (malloc (sizeof (struct popedge)));
    C[ci]->poptree[0].numup = 2;
    C[ci]->poptree[0].up = static_cast<int *> (malloc (2 * sizeof (int)));
    C[ci]->poptree[0].down = -1;
    C[ci]->poptree[0].time = TIMEMAX;
    C[ci]->poptree[0].b = 0;
    C[ci]->poptree[0].e = -1;
    C[ci]->poptree[0].up[0] = -1;
    C[ci]->poptree[0].up[1] = -1;
    //C[ci]->plist = NULL;
    fillplist (ci);
  }
  else if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
  {
    /* multiple branching population tree with a single split time */
    assert (numtreepops == npops + 1);
    C[ci]->poptree = static_cast<popedge *> (malloc ((npops + 1) * sizeof (struct popedge)));
    C[ci]->poptree[npops].numup = npops;
    C[ci]->poptree[npops].up = static_cast<int *> 
                (malloc (npops * sizeof (int)));
    C[ci]->poptree[npops].down = -1;
    C[ci]->poptree[npops].time = TIMEMAX;
    C[ci]->poptree[npops].b = 1;
    C[ci]->poptree[npops].e = -1;
    for (i = 0; i < npops; i++)
    {
      C[ci]->poptree[i].numup = 2;
      C[ci]->poptree[i].up = static_cast<int *> (malloc (2 * sizeof (int)));
      C[ci]->poptree[i].down = npops;
      C[ci]->poptree[i].time = TIMEMAX;
      C[ci]->poptree[i].b = 0;
      C[ci]->poptree[i].e = 1;
      C[ci]->poptree[i].up[0] = -1;
      C[ci]->poptree[i].up[1] = -1;
      C[ci]->poptree[npops].up[i] = i;
    }
    C[ci]->plist = static_cast<int **> (malloc (2 * sizeof (int *)));
    C[ci]->plist[0] = static_cast<int *> (malloc (npops * sizeof (int)));
    C[ci]->plist[1] = static_cast<int *> (malloc (sizeof (int)));
    for (i = 0; i < npops; i++)
      C[ci]->plist[0][i] = i;
    C[ci]->plist[1][0] = npops;

    C[ci]->periodset[0] = EMPTYSET;
    C[ci]->periodset[1] = EMPTYSET;
    for (i = 0; i < npops; i++)
    {
      C[ci]->periodset[0] = UNION (C[ci]->periodset[0], SINGLESET (i));
    }
    C[ci]->periodset[1] = UNION (C[ci]->periodset[1], SINGLESET (npops));
    C[ci]->rootpop = npops;
  }
  else
  {
    checktreestring (startpoptreestring);
    strcpy (C[ci]->poptreestring, startpoptreestring);

    //rewrite(startpoptreestring);  this will put the string in standard order,  not necessary now
    C[ci]->poptree = static_cast<popedge *> 
            (malloc (numtreepops * sizeof (struct popedge)));
    if (strlen (startpoptreestring) > 0)
    {
      poptreeread (ci, startpoptreestring);
      fillplist (ci);
    }
    else
    {
      // this was used for debugging poptreestring code 
      simpoptree (ci);
      poptreewrite (ci, C[ci]->poptreestring);
    }
  }


  return;
}                               /* setup_poptree */
