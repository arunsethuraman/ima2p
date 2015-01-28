/*IMa2p 2009-2015 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, and Arun Sethuraman */
#undef GLOBVARS
#include "imamp.hpp"
#include "updateassignment.hpp"

/* 

Calculating autocorrelations and ESS values:
---------------------------------------------
the general idea is to initialize things so that autocorrelations and ESS values can be calculated for a bunch 
of diverse things. 
At the start a set of three arrays of pointers are set up.  These point to the information needed, 
then they are passed to the function that does the calculation and the printing of values. 

void init_autoc_pointers(void):  initializes three arrays:
static struct autoc **autoc_pointer - an array of pointers to struct autoc which holds the information needed to do the calculations
static char **autoc_str_pointer - an array of pointers to strings that hold the names of things for which calculationsa 
are being done 
static double *autoc_vals - an array of numerical values that are being autocorrelated

See also 
static struct autoc **autoc_a_pointer = NULL;
static char **autoc_a_str_pointer = NULL;
static double *autoc_a_vals = NULL;
These are the same kinds of things, but for assignment related work


How to have something included in autocorrelation and ESS calculations and printouts:
-------------------------------------------------------------------------------------
The only things that should be needed, to add something new to the list of things to have 
autocorrelations and ESS values calculated is to have them included in the list of pointers in init_autoc_pointers() 
and to have the values to be included in the autocorrelation identifed in set_autoc_vals() 
Once things are set up,  the rest of the work should be done. 
So the only place changes should be needed is in init_autoc_pointers() and set_autoc_vals()


*/

static struct autoc **autoc_pointer = NULL;
static char **autoc_str_pointer = NULL;
static double *autoc_vals;
static int num_autoc = 0;

/* For ASSIGNMENT */
static struct autoc **autoc_a_pointer = NULL;
static char **autoc_a_str_pointer = NULL;
static double *autoc_a_vals = NULL;
static int num_autoc_a = 0;

/* these are the lag intervals overwhich autocorrelations are measured
 on 2/9/09  added AUTOCSTEPSCALAR as a multiplier of these to make 
 measurements over a much wider interval */

static int autoc_checkstep[AUTOCTERMS] =
  { 1, 10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000 };


/****** LOCAL FUNCTIONS ***********/
static void fillautoc (struct autoc *a, int n, double v);
static double calcautoc (struct autoc a);
static double printautocvalue (FILE * outto, struct autoc *ta);
static void integrate_autoc (FILE * outto, struct autoc ta[], double ac[]/*,int step */);
void printautoctable (FILE * outto, int numautoc, struct autoc **ac,
                      char **ac_str, const char *printacstring/*, int step*/);
void iautoc (struct autoc *a);  // set autoc values to zero
//AS: adding currentid to this 7/2/2014
void set_autoc_vals (double *pautoc_vals, int currentid);
void set_autoc_a_vals (double *pautoc_vals);
// struct autoc  used for recording values for measuring autocorrelations


/* fillauto() adds new values to the cov and var terms that are needed to calculate an autocorrelation
these are calculated using the value v, which is the current value of the variable being autocorrelated and a value
that was saved previously in position n of the the vals[] array */
void
fillautoc (struct autoc *a, int n, double v)
{
  a->cov.n++;
  a->cov.s += a->vals[n];
  a->cov.ss += v;
  a->cov.s2 += (a->vals[n] * v);
  a->var[0].s += a->vals[n];
  a->var[1].s += v;
  a->var[0].s2 += SQR (a->vals[n]);
  a->var[1].s2 += SQR (v);
}

/* calculated a correlation using the cov and var terms in a struct autoc */
double
calcautoc (struct autoc a)
{
  return (a.cov.s2 - (a.cov.s * a.cov.ss) / a.cov.n) / sqrt ((a.var[0].s2 - SQR (a.var[0].s) / a.cov.n) * (a.var[1].s2 - SQR (a.var[1].s) / a.cov.n));
}

double
printautocvalue (FILE * outto, struct autoc *ta)
{
  double ac;
  ac = calcautoc (*ta);
  if (ac > 1)
    ac = 1;
  fprintf (outto, "\t% .4f", ac);
  return ac;
}


/* integrate over autocorrelations values for values > AUTOCMIN */
#define AUTOCMIN  0.03
void
integrate_autoc (FILE * outto, struct autoc ta[], double ac[]/*, int step*/)
{
  double ess, temp1, temp2;
  long last_checkstep;
  int i, stopsum;
  ess = 0;
  stopsum = 0;
/* if AUTOCSTEPSCALAR is > 1 then we use the autocorrelation at step 0 and assume that it is equal to 1.
otherwise we use calculated autocorrelation measured for a lag of 1 */

  if (AUTOCSTEPSCALAR > 1)
    i = 0;
  else
    i = 1;
  for (; i < AUTOCTERMS; i++)
  {
    if (ta[i].cov.n > (AUTOCCUTOFF / AUTOCSTEPSCALAR))
    {
      if (ac[i] < AUTOCMIN)
        temp1 = 0;
      else
        temp1 = ac[i];
      if (AUTOCSTEPSCALAR > 1 && i == 0)
      {
        temp2 = 1.0;
        last_checkstep = 0;
      }
      else
      {
        if (ac[i - 1] < AUTOCMIN)
          temp2 = 0;
        else
          temp2 = ac[i - 1];
        last_checkstep = autoc_checkstep[i - 1];
      }
      if (stopsum == 0)
        ess += AUTOCSTEPSCALAR * (autoc_checkstep[i] - last_checkstep) * (temp1 + temp2) / 2;
      stopsum = (stopsum || temp1 == 0);
    }
  }
  if (stopsum)
    fprintf (outto, "\t%.0f", step / (1 + 2 * ess));
  else
    fprintf (outto, "\t< %.0f", step / (1 + 2 * ess));
}                               //integrate_autoc 



#define IFSTDOUTMAXSHOW 10
void
printautoctable (FILE * outto, int numautoc, struct autoc **ac, char **ac_str,
                 const char *printacstring/*, int step*/)
{
  int i, j;
  double **acvals;
  int numacprint;
  char numstr[20];

  /* 4/9/09  change output for stdout so it only does tables for L[P] and split times, not tmrcas 
   the total number of autoc tables for stdout will then be npops */
  //numacprint = (outto == stdout) ? IMIN (numautoc, IFSTDOUTMAXSHOW) : numautoc;
  if (strcmp (printacstring, "Population Assignment Autocorrelations and Effective Sample Size Estimates") == 0)
  {
    numacprint = (outto == stdout) ? IMIN (nloci, IFSTDOUTMAXSHOW) : numautoc;
  }
  else
  {
    numacprint = (outto == stdout) ? IMIN (numautoc, IFSTDOUTMAXSHOW) : numautoc;
  }
  acvals = orig2d_alloc2Ddouble (numacprint, AUTOCTERMS);
  fprintf (outto, "\n%s\n", printacstring);
  for (i = 0; i < (int) strlen (printacstring); i++)
    fprintf (outto, "-");
  fprintf (outto, "\n");
  fprintf (outto,
           "   # Steps Between Values and Autocorrelation Estimates \n");
  fprintf (outto, "\tSteps ");
  for (j = 0; j < numacprint; j++)
    fprintf (outto, "\t %s", ac_str[j]);
  fprintf (outto, "\n");
  for (i = 0; i < AUTOCTERMS; i++)
  {
    if (ac[0][i].cov.n > (AUTOCCUTOFF / AUTOCSTEPSCALAR))
    {
      if (((float) autoc_checkstep[i] * AUTOCSTEPSCALAR) > 1e4)
      {
        sprintf (&numstr[0], "%.1e",
                 (float) autoc_checkstep[i] * AUTOCSTEPSCALAR);
        fprintf (outto, "\t%s", shorten_e_num (&numstr[0]));
        //fprintf (outto, "\t%.1e", (float) autoc_checkstep[i]*AUTOCSTEPSCALAR);
      }
      else
      {
        fprintf (outto, "\t%d", autoc_checkstep[i] * AUTOCSTEPSCALAR);
      }
      for (j = 0; j < numacprint; j++)
      {
        acvals[j][i] = printautocvalue (outto, &ac[j][i]);
      }
      fprintf (outto, "\n");
    }
  }
  fprintf (outto, "\tESS");

/* integrate over autocorrelations values for values > 0.03 */

  for (j = 0; j < numacprint; j++)
    integrate_autoc (outto, ac[j], acvals[j]/*, step*/);
  fprintf (outto, "\n");
  orig2d_free2D ((void **) acvals, numacprint);
}                               //printautoctable

void
iautoc (struct autoc *a)        // initialize structure
{
  int i;
  for (i = 0; i < AUTOCTERMS; i++)
  {
    ieevent (&a[i].cov);
    ieevent (&a[i].var[0]);
    ieevent (&a[i].var[1]);
  }
}                               //iautoc 

/* this just puts the current values that are needed for the autorcorrelations all in one place */
void
set_autoc_vals (double *pautoc_vals, int currentid)
{
  int j, li, i = 0;
  int z = whichiscoldchain();
//  if ( z >= 0) {
  if (lpgpd_v->do_autoc == 1)
  {
    if (z >= 0 && currentid == 0) {	
    pautoc_vals[i] = C[z]->allpcalc.pdg  + C[z]->allpcalc.probg;
    i++;
    }
    #ifdef MPI_ENABLED
    if (z < 0 && currentid == 0) {
	double pautocrec = 0.0;
	MPI::COMM_WORLD.Recv(&pautocrec, 1, MPI::DOUBLE, MPI_ANY_SOURCE, 5656);
	pautoc_vals[i] = pautocrec;
	i++;
    }
    if (z >= 0 && currentid != 0) {
	double pautocrec = C[z]->allpcalc.pdg + C[z]->allpcalc.probg;
	MPI::COMM_WORLD.Send(&pautocrec, 1, MPI::DOUBLE, 0, 5656);
	i++;
    }
    #endif
	
  }
  for (j = 0; j < numsplittimes; j++)
  {
    if (T[j].v->do_autoc == 1)
    {
      if (z >= 0 && currentid == 0) {
	      pautoc_vals[i] = C[z]->tvals[j];
	      i++;
	}
	#ifdef MPI_ENABLED
	if (z < 0 && currentid == 0) {
		double pautocrec = 0.0;
		MPI::COMM_WORLD.Recv(&pautocrec, 1, MPI::DOUBLE, MPI_ANY_SOURCE, 6767);
		pautoc_vals[i] = pautocrec;
		i++;
	}
	if (z >= 0 && currentid != 0) {
		double pautocrec = C[z]->tvals[j];
		MPI::COMM_WORLD.Send(&pautocrec, 1, MPI::DOUBLE, 0, 6767);
		i++;
	}
	#endif
    }
  }
  for (li = 0; li < nloci; li++)
  {
    if (L[li].g_rec->v->do_autoc == 1)
    {
	if (z >= 0 && currentid == 0) {
	      pautoc_vals[i] = C[z]->G[li].roottime;
	      i++;
	}
	#ifdef MPI_ENABLED
	if (z < 0 && currentid == 0) {
		double pautocrec = 0.0;
		MPI::COMM_WORLD.Recv(&pautocrec, 1, MPI::DOUBLE, MPI_ANY_SOURCE, 7878);
		pautoc_vals[i] = pautocrec;
		i++;
	}
	if (z >= 0 && currentid != 0) {
		double pautocrec = C[z]->G[li].roottime;
		MPI::COMM_WORLD.Send(&pautocrec, 1, MPI::DOUBLE, 0, 7878);
		i++;
	}
	#endif
    }
  }
// }
}                               //set_autoc_vals

void
set_autoc_a_vals (double *pautoc_vals)
{
  int li;
  int i;

  i = 0;
  for (li = 0; li < nloci; li++)
  {
    if (L[li].a_rec->v->do_autoc == 1)
    {
      pautoc_vals[i] = imaAsnValue (0, li);
      i++;
    }
  }
}                               //set_autoc_vals

/*********** GLOBAL FUNCTIONS *************/


void
free_autoc_pointers (void)
{

  XFREE (autoc_pointer);
  XFREE (autoc_str_pointer);
  XFREE (autoc_vals);
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    XFREE (autoc_a_pointer);
    XFREE (autoc_a_str_pointer);
    XFREE (autoc_a_vals);
  }
}                               // free_autoc_pointers 

void
init_autoc_pointers (void)
{
  int i, j, li;
  if (autoc_pointer == NULL)
  {
    // count up the number of things to record ESS values 
    num_autoc += (lpgpd_v->do_autoc == 1);      // lpgpd
    if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0)
    {
      for (j = 0; j < numsplittimes; j++)
        num_autoc += (T[j].v->do_autoc == 1);
    }
    for (li = 0; li < nloci; li++)
      num_autoc += (L[li].g_rec->v->do_autoc == 1);     // tmrca values
    autoc_pointer = static_cast<autoc **> 
            (malloc (num_autoc * sizeof (struct autoc *)));
    autoc_str_pointer = 
            static_cast<char **> (malloc (num_autoc * sizeof (char *)));
    autoc_vals = static_cast<double *> (malloc (num_autoc * sizeof (double)));
    i = 0;
    if (lpgpd_v->do_autoc == 1)
    {
      autoc_pointer[i] = &(lpgpd_v->ac[0]);
      autoc_str_pointer[i] = lpgpd_v->strshort;
      i++;
    }
    if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0)
    {
      for (j = 0; j < numsplittimes; j++)
      {
        if (T[j].v->do_autoc == 1)
        {
          autoc_pointer[i] = &(T[j].v->ac[0]);
          autoc_str_pointer[i] = T[j].v->strshort;
          i++;
        }
      }
    }
    for (li = 0; li < nloci; li++)
      if (L[li].g_rec->v->do_autoc == 1)
      {
        autoc_pointer[i] = &(L[li].g_rec->v->ac[0]);
        autoc_str_pointer[i] = L[li].g_rec->v->strshort;
        i++;
      }
    assert (i == num_autoc);
    for (i = 0; i < num_autoc; i++)
      iautoc (autoc_pointer[i]);

    /* for ASSIGNMENT */
    if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
    {
      for (li = 0; li < nloci; li++)
        num_autoc_a += (L[li].a_rec->v->do_autoc == 1);
      autoc_a_pointer = static_cast<autoc **> 
            (malloc (num_autoc_a * sizeof (struct autoc *)));
      autoc_a_vals = static_cast<double *> 
            (malloc (num_autoc_a * sizeof (double)));
      autoc_a_str_pointer = static_cast<char **> 
            (malloc (num_autoc_a * sizeof (char *)));
      i = 0;
      for (li = 0; li < nloci; li++)
        if (L[li].a_rec->v->do_autoc == 1)
        {
          autoc_a_pointer[i] = L[li].a_rec->v->ac;
          autoc_a_str_pointer[i] = L[li].a_rec->v->strshort;
          i++;
        }
      for (i = 0; i < num_autoc_a; i++)
        iautoc (autoc_a_pointer[i]);
      assert (i == num_autoc_a);
    }

  }
  return;
}                               // init_autoc_pointers


/*called with start_autocorrelations == 1 at beginning and after burn 
AUTOCTERMS is a list of the different lag values in steps for which autocorrelations are recorded

nextstepcalc[AUTOCTERMS] for each lag value, the next step number at which another term is accumulated for the autocorrelation term for that lag
nextpossave[AUTOCTERMS] for each lag value, the position in the autoc_pointer[][].vals array that should get the next value to be recorded
nextposcalc[AUTOCTERMS] for each lag value, the position in the autoc_pointer[][].vals array that is involved in the next autocorrelation calculation
maxpos[AUTOCTERMS] for each lag value, the length of the autoc_pointer[][].vals array 
*/
void
checkautoc (int start_autocorrelations, int burndone, int burnsteps, int currentid)
{
  int i, j;
  int autoc_vals_recorded = 0;
  static int nextstepcalc[AUTOCTERMS],
    nextpossave[AUTOCTERMS], nextposcalc[AUTOCTERMS], maxpos[AUTOCTERMS];
  if (start_autocorrelations == 1)
  {
    for (i = 0; i < num_autoc; i++)
      iautoc (autoc_pointer[i]);

    /* For ASSIGNMENT */
    if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
    {
      for (i = 0; i < num_autoc_a; i++)
        iautoc (autoc_a_pointer[i]);
    }

    for (i = 0; i < AUTOCTERMS; i++)
    {
      nextstepcalc[i] = CHECKAUTOCWAIT + AUTOCINT * AUTOCSTEPSCALAR + autoc_checkstep[i] * AUTOCSTEPSCALAR + (burndone * burnsteps);
      nextpossave[i] = 0;
      nextposcalc[i] = 0;
      if (autoc_checkstep[i] <= AUTOCINT)
        maxpos[i] = 0;
      else
        maxpos[i] = (autoc_checkstep[i] / AUTOCINT) - 1;
      assert (maxpos[i] < AUTOCNEXTARRAYLENGTH);
    }

  }
  else
  {
    for (i = 0; i < AUTOCTERMS; i++)
    {
      if (step == nextstepcalc[i])
      {
       if (!autoc_vals_recorded)
       {
         set_autoc_vals (autoc_vals, currentid);
	///AS: Now if these values are on different processors, I need to broadcast them
	/*if (numprocesses > 1) {
		try {
			MPI::COMM_WORLD.Bcast(autoc_vals, num_autoc, MPI::DOUBLE, currentid);
		} catch (MPI::Exception e) {
			std::cout << "Bcast in autoc failed\n";
			MPI::COMM_WORLD.Abort(-1);
		}
	}*/
         if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
           set_autoc_a_vals (autoc_a_vals); 
         autoc_vals_recorded = 1;
       }
        for (j = 0; j < num_autoc; j++)
        {
          fillautoc (&autoc_pointer[j][i], nextposcalc[i], autoc_vals[j]);
        }

        /* For ASSIGNMENT */
        if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
        {
          set_autoc_a_vals (autoc_a_vals);      // must define this function or something like it
          for (j = 0; j < num_autoc_a; j++)
          {
            fillautoc (&autoc_a_pointer[j][i], nextposcalc[i], autoc_a_vals[j]);
            //fillautoc (&(L[li].a_rec->v->ac[i]), nextposcalc[i],IMA_assignment2value (0, li));
          }
        }

        nextposcalc[i]++;
        if (nextposcalc[i] > maxpos[i])
          nextposcalc[i] = 0;
        nextstepcalc[i] += AUTOCINT * AUTOCSTEPSCALAR;
      }
    }
    if ((long) AUTOCINT * AUTOCSTEPSCALAR * (long) (step / (AUTOCINT * AUTOCSTEPSCALAR)) == step)
    {
      // JH 9/30/09 bad bug here,  this next line should have been included but wasn't, also maybe two lines
      //  after that for assignment 
      if (!autoc_vals_recorded)
      {
        set_autoc_vals (autoc_vals, currentid);
        if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
          set_autoc_a_vals (autoc_a_vals);    
      }

      for (i = 0; i < AUTOCTERMS; i++)
      {
        for (j = 0; j < num_autoc; j++)
          autoc_pointer[j][i].vals[nextpossave[i]] = autoc_vals[j];
        /* For ASSIGNMENT */
        if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
        {
          for (j = 0; j < num_autoc_a; j++)
          {
            autoc_a_pointer[j][i].vals[nextpossave[i]] = autoc_a_vals[j];
            //L[li].a_rec->v->ac[i].vals[nextpossave[i]] = IMA_assignment2value (0, li);
          }
        }

        nextpossave[i]++;
        if (nextpossave[i] > maxpos[i])
          nextpossave[i] = 0;
      }
    }
  }
}                               /*checkautoc */


void
callprintautoctable (FILE * outto/*, int step */)
{

  if (autoc_pointer[0][0].cov.n > (AUTOCCUTOFF / AUTOCSTEPSCALAR))      // check one of the autocorrelation accumulators to see if enough counts have been made
  {
    printautoctable (outto, num_autoc, autoc_pointer, autoc_str_pointer,
                     "Autocorrelations and Effective Sample Size Estimates"/*, step*/);
    /* For ASSIGNMENT */
    if (assignmentoptions[POPULATIONASSIGNMENT] == 1
        && assignmentoptions[POPULATIONASSIGNMENTASSIGNED] == 0)
    {
      printautoctable (outto, num_autoc_a, autoc_a_pointer,
                       autoc_a_str_pointer,
                       "Population Assignment Autocorrelations and Effective Sample Size Estimates"/*,  step*/);
    }
  }
}                               // callprintautoctable 
