/*IMa2p 2009-2015 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, and Arun Sethuraman */

#define GLOBVARS
#include "imamp.hpp"
#include "updateassignment.hpp"
#ifndef WIN32
/* cr 110524.1  config.h part of Auto make process.  File will not be 
 * included if the NOCONFIG define is set in the makefile.  This change 
 * only affects build on platforms other than Windows.
 */
//#ifndef NOCONFIG
//#include <config.h>
//#endif /* NOCONFIG */
#endif

/* MPI definitions */
#ifdef MPI_ENABLED
#include <mpi.h>
#endif

#include <stdlib.h>

extern struct edgemiginfo oldedgemig;
extern struct edgemiginfo oldsismig;
extern struct edgemiginfo newedgemig;
extern struct edgemiginfo newsismig;

#define BURNTRENDSTARTDELAYDEFAULT  10000

/* Global variables => as of now, each process will store a copy of this */
/* AS: Should be more efficient to initialize just once on the head node, then broadcast this to all processes */
static FILE *outfile;
//static FILE *autocfile;
///AS: adding a dummy variable to carry the array pointed to by outfilename, so I can broadcast it
char *outfilename;
char outfile_name[FNSIZE];
char *trueassignment;
int gsampinflength; /* updateassignment.c needs to know this. */

/* This is the number of tries of updating assignment, which is equal to number
 * of individuals with unknown origin */
static int snupdatei;              

static char *loadfilebase;
static time_t starttime;
static time_t endtime;
static time_t chainstarttime;
static time_t timer;
static time_t lasttime;
static time_t remained_starttime;
static time_t remained_endtime;
static long burnduration, chainduration;
static int burndone;
static long int burnsteps;
static int burntrendstartdelay;
static int burndurationmode, cdurationmode;
static int genealogiestosave;
static int memforgenealogiessaved = 0;
static int savegenealogyint;
static int recordint;
static int printint;
static double generationtime;
static double scaleumeaninput = 0;
static int swaptries;
static int heatmode;
static char fpstr[50000];       // probably long enough
static int *fpstri;
static char oldoutfilename[FNSIZE];
static double hilocuslike[MAXLOCI];
static int numgenealogyfiles;
static FILE *genealogyinfosavefile;
static char genealogyinfosavefilename[FNSIZE];
static long int recordstep = 0;
static double hilike = -1e20, hiprob = -1e20;
static FILE *checkdonefile;
static double hval1, hval2;
static int trenddoublepoint;
static int trendspot = 0;
static int maxedoutgenealogysave = 0;
static char priorfilename[FNSIZE];
static char mcfwritefilename[FNSIZE];
static char mcfreadfilename[FNSIZE];
static char command_line[1001];
static char heatingterm_str[50], modeloptions_str[50], calcoptions_str[50], outputoptions_str[50];
static char *infilename;
///AS: Adding this so I can broadcast this
static char infile_name[FNSIZE];
static long seed_for_ran1;
static int **migcount, *migfrom, *migto;
static char migplotfilename[FNSIZE];
static char migrationnamefilename[FNSIZE];
static int migrationnamefrom,migrationnameto;
static FILE *migrationnamefile;
static FILE *migplotfile;
static char nestedmodelfilename[FNSIZE] = "\0";
static char defaultpriorfilename[14]= "imapriors.txt";
static int *swapper;
static int *swappee;
static int swapbetasonly = 1;

/*Local function prototypes  */
static void init_after_start_IMA ();
static void init_IMA ();
static void scan_commandline (int argc, char *argv[], int curr_id); ///AS: Adding curr_id
static void print_outputfile_info_string (void);
static void start (int argc, char *argv[], int currentid); ///AS: Adding currentid to all the files to be printed
static void qupdate (int curr_id, /*int *swapper, int *swappee, std::ofstream &f1, */int swapA, int swapB);
static void savegenealogyinfo (int current_id); ///AS: Adding current_id
static void reset_after_burn (int current_id); ///AS: Adding current_id
static void output_burntrendfile (int current_id); ///AS: Adding current_id/
static int run (int current_id); ///AS: Adding current_id
static void inctrend (int m, int t, struct value_record *v, double newval);
static void trend_reset (struct value_record *v, int nv);
static void trendrecord (int loadarrayj, int currentid);
static void recordval (struct value_record *v, double val);
static void record_migrations (int z, int currentid);
static void record (int currentid);
static void loadgenealogyvalues (void);
static void callasciicurves (void);
static void callasciitrend (FILE * outfile);
static void printoutput (int currentid);
static void intervaloutput (FILE * outto, int currentid);
static void free_ima_main_stuff ();
static void callprintacceptancerates (FILE * outto, int currentid);
static void printsteps (FILE * outto, double like);
static void check_to_record (int current_id); ///AS: Adding current_id
static void record_migration_names();
static void fillswaparrays (int *swapper, int *swappee);
//static int whichiscoldchain (void);

int main ( int argc, char *argv[]);


/* SANGCHUL: Mon Jan 12 20:32:41 EST 2009
 * All global variables can be initialized in function [[init_IMA]].
 * All variables must be finalized in function [[free_IMA]].
 * Sometimes we need to initialze some global variables after we call function
 * [[start]]. Function [[init_after_start_IMA]] does this job.
 * */
void
init_IMA ()
{
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    imaAsnInit ();
  }

  lastgenealogysaved = -1;
  return;
}

void
init_after_start_IMA ()
{
#ifdef DEBUG
  int ci;
#endif /* DEBUG */
  char Structurama[FNSIZE];

  if (assignmentoptions[PRINTSTRUCTURAMA] == 1)
    {
      strcpy (Structurama, outfilename);
      strcat (Structurama, ".in");
      IMA_convert_IM2Structurama (Structurama);
      IMA_output_structurama_bat (outfilename);
    }

  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    recordassignment_header (outfilename, 0, 0);
#ifdef DEBUG
    if (numchains > 1)
      {
        for (ci = 0; ci < numchains; ci++)
        {
          recordassignment_header (outfilename, ci, 1);
        }
      }
#endif /* DEBUG */
  }
  /* Function IMA_ninds returns the number of individuals with their label being
   * unknown. */
  /* snupdatei = IMA_nindsunknown (); */
  snupdatei = 1;
  return;
}


void
free_ima_main_stuff ()
{
  int i;
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    imaAsnFree ();
  }
  if (memforgenealogiessaved > 0)
  {
    for (i = 0; i < memforgenealogiessaved; i++)
    {
      XFREE (gsampinf[i]);
    }
    XFREE (gsampinf);
  }
  XFREE (fpstri);
  if (modeloptions[EXPOMIGRATIONPRIOR] || (runoptions[LOADRUN] && calcoptions[FINDJOINTPOSTERIOR]))
     XFREE(eexpsum);

  if (outputoptions[MIGRATEHIST])
  {
    for (i = 0; i < nloci + (nloci > 1); i++)
    {
      XFREE (migcount[i]);
    }
    XFREE (migcount);
    XFREE (migfrom);
    XFREE (migto);
  }

  if (trueassignment != NULL)
  {
    XFREE (trueassignment);
  }

  if (infilename != NULL)
  {
    XFREE (infilename);
  }

  if (outfilename != NULL)
  {
    XFREE (outfilename);
  }

  if (loadfilebase != NULL)
  {
    XFREE (loadfilebase);
  }
  XFREE (swapper);
  XFREE (swappee);

	
  freeanymemory ();

}                               //free_ima_main_stuff


#define GENERATIONTIMEDEFAULT   1
#define DEFAULTNUMGENEALOGIES 10000

void
scan_commandline (int argc, char *argv[], int curr_id)
{
  static int i, j, k;
  static char ch, ch1;
  static char pstr[256];
  /* option flags that, although values are assigned to the flags, the flag
   * is never used by the code:  Cp, Ep, Wp, Yp.  These variable flags
   * are being left in for completeness
   */
  int Ap, Bp, Cp, Dp, Ep, Fp, Gp, Hfp, Hnp, Hkp, Hap, Hbp, Ip, Jp, Lp, Mp, Op, Pp,
    Qp, Rp, Sp, Tp, Up, Vp, Wp, Yp, Zp;
  int Xp;
  double tempf;
  char *opt;
  const char *sep = ",:";
  int ioption;
  
  /* 5/19/2011 JH adding thermodynamic integration */
  int Kp = 0;

  time (&starttime);
  time (&remained_starttime);


  Ap = 0;                       /* Assignment options */
  Bp = 0;                       /* duration of burnin */
  Cp = 0;                       /* calculation options  - flag not used */
  Dp = 0;                       /* number of steps in between genealogy saves */
  Ep = 0;                       /* prior on splitting rate parameter -
                                        flag not used */
  Fp = 0;                       /* name of mcf file */
  Gp = 0;                       /* name of prior file */
  Hfp = 0;                      /* heating model */
  Hnp = 0;                      /* # of chains */
  Hkp = 0;                      /* # of swap attempts */
  Hap = 0;                      /* heat term1 */
  Hbp = 0;                      /* heat term2 */
  Ip = 0;                       /*input file */
  Jp = 0;                       /* used for programmer options */
  Lp = 0;                       /* duration of chain */
  Mp = 0;                       /* migration rate max */
  Op = 0;                       /* output file */
  Pp = 0;                       /* output options */
  Qp = 0;                       /* Theta max scalar */
  Rp = 0;                       /* run options */
  Sp = 0;                       /* random number seed */
  Tp = 0;                       /* Time maximum */
  Up = 0;                       /* generation time in years */
  Vp = 0;                       /* genealogy load file name base */
  Wp = 0;                       /* name of file with nested models -  
                                            flag not used */
  Yp = 0;                       /* mutation rate scalar for loci with mutation rates given in input file - for use with LOADRUN mode  - flag not used */
  Zp = 0;                       /* screen printout frequency */
  Xp = 0;                       /* True Assignment */
  trueassignment = NULL;
  if (curr_id == 0) {
	  printf ("executing program ...Scanning commandline happens only on the head node\n");
  }
  if (((argc == 2 && (char) toupper (argv[1][1]) == 'H') || argc == 1) && curr_id == 0)
  {
    printf ("IMa2p Program - copyright 2015 by Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, and Arun Sethuraman\n");
    printf ("Release date: %s\n",RELEASE_DATE);
    printf ("This program is run through a command line interface\n");
    printf ("To execute the program, type the program name followed by the necessary command line flags and options \n");
/*  assignment options not included in output for now 4/27/09  jhey 
    printf ("-a  Population assignment options: \n");
    printf ("    1  Invoke population assignment\n");
    printf ("    2  Invoke DNA barcoding\n");
    printf ("    3  Island model\n"); */
/*
    printf ("    0  Turn on check-points\n");
    printf ("    1  Invoke population assignment\n");
    printf ("    2  Relabel update\n");
    printf ("    3  Beerli and Felsenstein update\n");
    printf ("    4  Print info. for DNA Barcoding\n");
    printf ("    5  Island model\n");
    printf ("    6  Local assignment of genes\n");
    printf ("    7  Print input file for Structurama\n");
    printf ("    e.g, -a12 relabel update of assignment with population tree model\n");
    printf ("         -a13 Beerli and Felsenstein update of assignment with population tree model\n");
    printf ("         -a135 Beerli and Felsenstein update of assignment with island model\n");
    printf ("         -a012: turn on checking genealogy integrity additional to option 1,2\n");
    printf ("         -a124: print assignment proportion of a single unknown gene\n");
    printf ("         -a125: relabel update of assignment with island model\n");
    printf ("         -a127: print out the STRUCTURAMA input file additional to option 1,2\n");
    printf ("         -a126: local assignment of genes are allowed\n");
Island model must be in -j option.
1 Population Assignment, BF99 update, -a7 should be turned on
2 DNA Barcoding, NM05 update, -a4 should be turned on
-j6 is no migration. BF99 can be implemented for no migration case.
NM05 can be implemented for no migration case as well.
-a1 : Assignment, Print structurama output, BF99 or NM05 or both
-a2 : DNA Barcoding, Print info. for DNA Barcoding, BF99 or NM05 or both
-a1 -j6 : Assignment with no-migration
-a2 -j6 : DNA Barcoding with no-migration
-a13 : Assignment with island model
-a23 : DNA Barcoding with island model
*/
    printf ("-b  Duration of burn  (MCMC mode only)\n");
    printf ("    - If integer, the number of burnin steps \n");
    printf ("    - If floating point, the time in hours between writing of burntrend file\n");
    printf ("         run continues until file " "IMburn" " is no longer present\n");
    printf ("         in the directory, or if present, does not begin with 'y'\n");
    printf ("-c  Calculation options: \n");
    printf ("    0 Likelihood of data functions return a constant - posterior should equal prior \n");
    printf ("    1 Include ranges on mutation rates as priors on mutation rate scalars\n");
    printf ("    2 Joint posterior density calculations, for LLR tests of nested models use with -w (LOAD-GENEALOGY mode only)\n");
    printf ("    3 Get prior distribution terms from file (requires filename given with -g )\n");
    /* 5/19/2011 JH adding thermodynamic integration */
    /* as of 8/23/2011 this is still in progress,  do not show this help output in release versions 
    printf ("    4 Calculate the marginal likelihood, must specify -hn (odd number > 50 chains),  -hf, -ha, and -hb are ignored\n"); */
    printf ("-d  Number of steps between genealogy saving (MCMC mode only) (default 100)\n");
    printf ("-f  Name of file with saved Markov chain state generated in previous run - use with -r3\n");
    printf ("-g  Name of file with parameter priors  (requires -c3) default: '%s'\n",defaultpriorfilename);
    //printf("-jh  jh personal options \n");
    //printf("      0  alt data format for SW data -  one data line for each allele in each pop,  1st # is allele length, 2nd is # copies \n");
    //printf("      1  write migration
    //printf("      2  ghost island - used with a ghost population.  There is split time for the ancestral population, instead  the ghost and the ancestor of sampled populations are two islands with gene flow 
    printf ("-h  Heating terms (MCMC mode only): \n");
    printf ("  -hf Heating model: l linear (default); g geometric\n"); /* 8/23/2011  remove HTWOSTEP from use */
    printf ("  -hn Number of chains \n");
    printf ("  -hk Number of chain swap attempts per step (default = number of chains)\n");
    printf ("  -ha First heating parameter, effect depends on heating model \n");
    printf ("  -hb Second heating parameter, effect depends on heating model  \n");
    printf ("-i  Input file name (no spaces) \n");
    printf ("-j  Model options: \n");
    printf ("    1  Migration only between sister populations (no migration between non-sister populations)\n");
    printf ("    2  One migration parameter for each pair of populations (do not use with -p5)\n");
    printf ("    3  Migration only between sampled populations (ancestral populations have zero migration)\n");
    printf ("    4  Add a non-sampled ghost population to the model \n");
    printf ("    5  Separate population size and migration parameters in each period (lots of parameters) \n");
    printf ("    6  No migration in the model\n");
    printf ("    7  Migration prior follows exponential distribution with mean given by -m or in parameter prior file \n");
    printf ("    8  Each ancestral population size is assumed to be identical to that of their largest descendant population\n");
    printf ("    9  One single migration parameter for all pairs of populations (do not use with -p5)\n");
    printf ("-l  Run duration (default: %d genealogies sampled per locus):\n", DEFAULTNUMGENEALOGIES);
    printf ("     If in MCMC mode (i.e. not loading genealogies from a previous run) \n");
    printf ("       - If integer, the number of genealogies to save\n");
    printf ("         This value times -d value sets the # of steps in chain after burnin) \n");
    printf ("	    - If floating point, the time in hours between outputs. \n");
    printf ("         Run continues until file " "IMrun" " is no longer present\n");
    printf ("           in the directory, or if present, does not begin with 'y'\n");
    printf ("     If in load-genealogy mode (i.e. using -r0 to load genealogies from previous run)\n");
    printf ("       - Integer indicates number of genealogies to load from file(s) named with -r\n");
    printf ("-m  Migration prior value (maximum for uniform,  mean if exponential distribution is used \n");
    printf ("-o  Output file name (no spaces) default is 'outfile.txt' \n");
    printf ("-p  Output options: \n");
    printf ("    0 Turn off trend plots in outfile (default is to print trend plots)\n");
    printf ("    1 Turn off plots of marginal curves in outfile (default is to print marginal density plots)\n");
    printf ("    2 Print TMRCA histogram for each genealogy (MCMC mode only)\n");
    printf ("    3 Print histogram of parameters on demographic scales  (requires mutation rate(s) in data file)\n");
    printf ("    4 Print histogram of splitting times divided by prior (do not use with -j0 or when only 2 sampled populations\n");
    printf ("    5 Print estimates and histograms of population migration rate (2NM)\n");
    printf ("    6 Print pairwise probabilities that one parameter is greater than another \n");
    /* CR:110114.2  message text changed */
    printf ("    7 Print histograms of the number of migration events (MCMC mode only)\n");
    printf ("    8 Print joint estimate for splitting times (MCMC mode only, for models with 3, 4 or 5 populations)\n");
    printf ("-q  Maximum for population size parameters (4Nu) \n");
    printf ("-r  Run options \n");
    printf ("    0 LOAD-GENEALOGY Mode - load genealogies from previous run(s); also requires -v \n");
    printf ("    1 Do not save genealogies to a file (default saves sampled genealogies) \n");
    printf ("    2 Save the state of the Markov chain in a file - named with extension .mcf (MCMC mode only)\n");
    printf ("    3 Start run by loading a previously saved *.mcf file; requires -f (data and priors must be the same) \n");
    printf ("    4 Write all mutation related updates rates to stdout during the run (default is to suppress this)\n");
    printf ("    5 Print burntrend file at end of burnin period; use with -b followed by integer (MCMC mode only)\n");
    printf ("-s  Random number seed (default is taken from current time)\n");
    printf ("-t  Maximum time of population splitting\n");
    /*printf ("-t  Maximum time of population splitting ( do not use with -e) \n"); */
    printf ("-u  Generation time in years - for use with -p3 (default is %d) \n", GENERATIONTIMEDEFAULT);
    printf ("-v  Base name (no extension) of *.ti files with genealogy data  (requires use of -r0) \n");
    printf ("-w  Name of file with nested models to be tested (LOAD-GENEALOGY mode only), invokes -c2\n");
/*     printf ("-x  beta for raising the power to likelihood\n"); */
    printf ("-y  Mutation rate scalar for relevant loci - for use with -p3 \n");
    printf ("-z  Number of steps between screen output (default is %d) (MCMC mode only)\n",PRINTINTDEFAULT);
    

  #ifdef MPI_ENABLED
  MPI::Finalize();
  #endif
    exit (0);
  }
  else
  {
/*
command line circumstances:
all flags begin with '-'
-most flags are single letter flags
-some are double letter flags: h
-some flags are followed by a string or a character
others by a single number (int or float)
others by a string of integers

it is ok to have spaces between a flag and its values 

All flags are followed by at least something
no flag is followed by nothing 
*/
    strcpy (command_line, "");
    strcpy (heatingterm_str,"");
    strcpy (modeloptions_str,"");
    strcpy (calcoptions_str,"");
    strcpy (outputoptions_str,"");

    for (i = 1; i < argc; i++)
    {
      strcpy (pstr, argv[i]);
      strcat (command_line, " ");
      strcat (command_line, pstr);

      if (strlen (pstr) < 2)
        IM_err (IMERR_COMMANDLINEFORMAT, " one of the command line strings is too short: %s ",pstr);
      if (pstr[0] != '-')
        IM_err (IMERR_COMMANDLINEFORMAT, "command line flag not preceded by '-' : %s", pstr);
      ch = toupper (pstr[1]);
      ch1 = ' ';
      if (ch == 'H')
      {
        strcat (heatingterm_str, " ");
        strcat (heatingterm_str, pstr);
      }
      if (ch == 'J')
      {
        strcat (modeloptions_str, " ");
        strcat (modeloptions_str, pstr);
        ch1 = toupper (pstr[2]);
      }
      if (ch == 'C')
      {
        strcat (calcoptions_str, " ");
        strcat (calcoptions_str, pstr);
      }
      if (ch == 'P')
      {
        strcat (outputoptions_str, " ");
        strcat (outputoptions_str, pstr);
      }
      if (ch == 'H')
        ch1 = toupper (pstr[2]);
      if (strlen (argv[i]) == 2 || (i < argc - 1 && isdigit (argv[i + 1][0])))  // space separates flag from its number
      {
        i++;
        strcpy (pstr, argv[i]);
        strcat (command_line, " ");
        strcat (command_line, pstr);
      }
      else
      {
        if ((ch == 'H'))
          strdelete (pstr, 1, 3);
        else
          strdelete (pstr, 1, 2);
      }
      switch ((char) toupper (ch))
      {
      case 'A':
        j = (int) (strlen (pstr) - 1);
        if (comma_exists (pstr))
        {
          for (opt = strtok (pstr, sep); opt; opt = strtok (NULL, sep))
          {
            ioption = atoi (opt);
            if (ioption < 0 || ioption >= POPULATIONASSIGNMENT_NUMBER)
            {
              IM_err (IMERR_COMMANDLINEFORMAT, "option -a %s", pstr);
            }
            assignmentoptions[ioption] = 1;
          }
          Ap = 1;
        }
        else
        {
          while (j >= 0)
          {
            ioption = atoi (&pstr[j]);
            if (ioption < 0 || ioption >= POPULATIONASSIGNMENT_NUMBER)
            {
              IM_err (IMERR_COMMANDLINEFORMAT, "option -a %s", pstr);
            }
            assignmentoptions[ioption] = 1;
            pstr[j] = '\0';
            j--;
            Ap = 1;
          }
        }
        break;
      case 'B':
        tempf = atof (&pstr[0]);
        /* check to see if the value is floating point, in which case treat it as being in fractions of an hour  and convert to seconds */
        if (strchr (pstr, '.'))
        {
          burnduration = (int) (3600 * tempf);
          burndurationmode = TIMEINF;
	///AS: Should time be broadcast or should each processor have its own time?
          time (&lasttime);
          runoptions[PRINTBURNTREND] = 1;
        }
        else
        {
          burnduration = (int) tempf;
          burndurationmode = TIMESTEPS;
        }
	
	#ifdef MPI_ENABLED
	  try {
		MPI::COMM_WORLD.Bcast(&burnduration, 1, MPI::LONG, curr_id);
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	//std::cout << "Burn duration after broadcast is " << burnduration << "\n";
	try {
		MPI::COMM_WORLD.Bcast(&burndurationmode, 1, MPI::INT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
	//std::cout << "Burn duration mode after broacast is " << burndurationmode << "\n";
        Bp = 1;
        break;
      case 'C':
        j = (int) (strlen (pstr) - 1);
        while (j >= 0)
        {
          if (!isdigit (pstr[j]))
            IM_err (IMERR_COMMANDLINEFORMAT, "calculation option flag -c should be followed by a digit: %s",pstr);
          calcoptions[atoi (&pstr[j])] = 1;
          pstr[j] = '\0';
          j--;
        }
        break;
      case 'D':
        savegenealogyint = atoi (&pstr[0]);
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&savegenealogyint, 1, MPI::INT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
        Dp = 1;
        break;
      case 'F':
        strcpy (mcfreadfilename, pstr);
        put_spaces_in_filepaths(mcfreadfilename);
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&mcfreadfilename, 500*sizeof(char), MPI::CHAR, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
        Fp = 1;
        break;
      case 'G':
        strcpy (priorfilename, pstr);
        put_spaces_in_filepaths(priorfilename);
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&priorfilename, 500*sizeof(char), MPI::CHAR, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
        Gp = 1;
        break;
      case 'H':
        switch ((char) toupper (ch1))
        {
        case 'A':
          hval1 = atof (&pstr[0]);
          Hap = 1;
          break;
        case 'B':
          hval2 = atof (&pstr[0]);
          Hbp = 1;
          break;
        case 'N':
          numchains = atoi (&pstr[0]);
          Hnp = 1;
          break;
        case 'K':
          swaptries = atoi (&pstr[0]);
          Hkp = 1;
          break;
        case 'F':
          Hfp = 1;
          switch ((char) toupper (pstr[0]))
          {
         /* case 'T':
            heatmode = HTWOSTEP;
            break; */
          case 'G':
            heatmode = HGEOMETRIC;
            break;
          default:
            heatmode = HLINEAR;
            break;
          }
          break;
        default:
          IM_err (IMERR_COMMANDLINEHEATINGTERMS, "mistake in use of -h flag : %s", pstr);
        }
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&hval1, 1, MPI_FLOAT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	
	
	try {
		MPI::COMM_WORLD.Bcast(&hval2, 1, MPI_FLOAT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	
	try {
		MPI::COMM_WORLD.Bcast(&numchains, 1, MPI::INT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	
	try {
		MPI::COMM_WORLD.Bcast(&swaptries, 1, MPI::INT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	
	try {
		MPI::COMM_WORLD.Bcast(&heatmode, 1, MPI::INT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
        break;
      case 'I':
        infilename = strdup (pstr);
        put_spaces_in_filepaths(infilename);
	strncpy(infile_name, infilename, sizeof(infilename));	
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&infile_name, 500*sizeof(char), MPI::CHAR, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
        Ip = 1;
        break;
      case 'J':
        Jp = 0;
        if (!(toupper(ch1) == 'H'))
        {
          j = (int) (strlen (pstr) - 1);
          while (j >= 0)
          {
            if (!isdigit (pstr[j]))
              IM_err (IMERR_COMMANDLINEFORMAT, "model option flag -j should be followed by a digit: %s", pstr);
            modeloptions[atoi (&pstr[j])] = 1;
            pstr[j] = '\0';
            j--;
          }
        }
        else
        {
          j = (int) (strlen (pstr) - 1);
          if (isdigit(pstr[1])  && pstr[1] == '1')
          {
            migrationnamefrom = atoi (&pstr[2]);
            migrationnameto = atoi (&pstr[3]);
            strcpy(migrationnamefilename, &pstr[4]);
            put_spaces_in_filepaths(migrationnamefilename);
		///Bcast migrationnamefilename
            jheyoptions[WRITEMIGRATIONNAME] = 1;
          }
          else
          {
            while (j >= 0)
            {
              if (isdigit(pstr[j]))
              {
                jheyoptions[atoi (&pstr[j])] = 1;
              }
              pstr[j] = '\0';
              j--;
            }
          }
        };
        break;
      case 'L':
        tempf = atof (&pstr[0]);
        /* check to see if the value is floating point, in which case treat it as being in fractions of an hour  and convert to seconds */
        if (strchr (pstr, '.'))
        {
          chainduration = (int) (3600 * tempf);
          cdurationmode = TIMEINF;
          genealogiestosave = -1;
        }
        else
        {
          genealogiestosave = (int) tempf;
          cdurationmode = TIMESTEPS;
        }
        Lp = 1;
	#ifdef MPI_ENABLED
	if (numprocesses > 1) {
	try {
		//std::cout << "Broadcasting chainduration and it is " << chainduration << "\n";
		MPI::COMM_WORLD.Bcast(&chainduration, 1, MPI::LONG, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	
	try {
		//std::cout << "Broadcasting cdurationmode and it is " << cdurationmode << "\n";
		MPI::COMM_WORLD.Bcast(&cdurationmode, 1, MPI::INT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	
	try {
		MPI::COMM_WORLD.Bcast(&genealogiestosave, 1, MPI::INT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	}
	#endif

        break;
      case 'M':
        mprior = (double) atof (&pstr[0]);
        if (mprior == 0)
          modeloptions[NOMIGRATION] = 1;
        Mp = 1;
	#ifdef MPI_ENABLED	
	try {
		MPI::COMM_WORLD.Bcast(&mprior, 1, MPI_FLOAT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
        break;
      case 'O':
        outfilename = strdup (pstr);
        put_spaces_in_filepaths(outfilename);
	///AS: Has to be broadcast too?
        Op = 1;
        break;
      case 'P':
        Pp = 1;
        j = (int) (strlen (pstr) - 1);
        while (j >= 0)
        {
          if (!isdigit (pstr[j]))
          {
            IM_err (IMERR_COMMANDLINEFORMAT, "print option flag -p should be followed by a digit : %s",pstr);
          }
          k = atoi (&pstr[j]);
          outputoptions[k] = 1;
          pstr[j] = '\0';
          j--;
        }
        break;
      case 'Q':
        thetaprior = (double) atof (&pstr[0]);
        Qp = 1;
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&thetaprior, 1, MPI_FLOAT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
        break;
      case 'R':
        Rp = 1;
        j = (int) (strlen (pstr) - 1);
        while (j >= 0)
        {
          if (!isdigit (pstr[j]))
            IM_err (IMERR_COMMANDLINEFORMAT, "run option flag -r should be followed by a digit : %s ", pstr);
          k = atoi (&pstr[j]);
          runoptions[k] = 1;
          pstr[j] = '\0';
          j--;
        }
        break;
	///AS: This has to change - I need to set a different seed on each process
      case 'S':
        seed_for_ran1 = atoi (&pstr[0]);
	seed_for_ran1 = seed_for_ran1 * curr_id; ///Just to make sure that there's a different seed on each process
        if (!seed_for_ran1)
          seed_for_ran1 = curr_id;
	/*
	try {
		MPI::COMM_WORLD.Bcast(seed_for_ran1, 1, MPI::INT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}*/
        Sp = 1;
        break;
      case 'T':
        tprior = (double) atof (&pstr[0]);
        Tp = 1;
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&tprior, 1, MPI_FLOAT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
        break;
      case 'U':
        generationtime = atof (&pstr[0]);
        Up = 1;
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&generationtime, 1, MPI_FLOAT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
        break;
      case 'V':
        loadfilebase = strdup (pstr);
        put_spaces_in_filepaths(loadfilebase);
        Vp = 1;
        break;
      case 'W':
        strcpy (nestedmodelfilename, pstr);
        put_spaces_in_filepaths(nestedmodelfilename);
        calcoptions[FINDJOINTPOSTERIOR] = 1;
	///AS: nestedmodelfilename to be broadcast too?
        break;
      case 'Y':
        scaleumeaninput = atof (&pstr[0]);
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&scaleumeaninput, 1, MPI_FLOAT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
        break;
      case 'Z':
        printint = atoi (&pstr[0]);
        Zp = 1;
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&printint, 1, MPI::INT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
        break;
      case 'X':
        if (strchr (pstr, '.') == NULL)
        {
          trueassignment = strdup (pstr);
          gbeta = 1.0;
        }
        else
        {
          gbeta = atof (&pstr[0]);
          trueassignment = NULL;
        }
        Xp = 1;
	/*
	try {
		MPI::COMM_WORLD.Bcast(trueassignment, trueassignment.size(), MPI::CHAR, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	*/
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&gbeta, 1, MPI_FLOAT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
        break;
      default:
        IM_err (IMERR_COMMANDLINEFORMAT, &ch);
      }
    }
  }

//AS: Jared Knoblauch caught a bug here - if the numprocessors > 1, and heating parameters aren't specified
//the program hangs. Fixing this Thu Aug  4 13:40:59 EDT 2016
	#ifdef MPI_ENABLED
	if (numprocesses > 1) {
		if (Hnp == 0)
			IM_err (IMERR_COMMANDLINEHEATINGTERMS, "no -hn flag specified despite attempting to run on %d processors", numprocesses);
		if (Hnp == 1 && Hfp == 0)
			IM_err (IMERR_COMMANDLINEHEATINGTERMS, "no -heating mode specified");
	}
	#endif
  /* Check if command line options are compatible with single population. */
  if (infilename==0)
    IM_err (IMERR_READFILEOPENFAIL,  "pointer to input file not set,  check -i on command line");
  npops = imaInfileNpops (infilename);
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&npops, 1, MPI::INT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
  if (npops < 1 || npops > 10)
  {
    IM_err (IMERR_COMMANDLINEFORMAT, "Number of populations must be nonnegative and less than 10");
  }
  if (npops == 1)
  {
    if (modeloptions[NOMIGRATION] == 0)
      modeloptions[NOMIGRATION] = 1;
    if (Mp == 1)
    {
      IM_err (IMERR_COMMANDLINEFORMAT, "model option [single population] should not go with -m");
    }
    if (Tp == 1)
    {
      IM_err (IMERR_COMMANDLINEFORMAT, "model option [single population] should not go with -t");
    }
  }
  
///AS: Ignoring assignment options for now
  if (Ap == 1)
  {
    /* Allowed command line options for assignment work.
     * -a0, -a01, -a02, -a03, -a013, -023,
     * -a1, -a13
     * -a2, -a23
     * -a3,
     */
    if (assignmentoptions[POPULATIONASSIGNMENT] == 1
        || assignmentoptions[POPULATIONASSIGNMENTRELABEL] == 1
        || assignmentoptions[POPULATIONASSIGNMENTBF] == 1
        || assignmentoptions[POPULATIONASSIGNMENTASSIGNED] == 1
        || assignmentoptions[PRINTSTRUCTURAMA] == 1
        || assignmentoptions[JCMODEL] == 1
        || assignmentoptions[POPULATIONTREEWCH] == 1)
    {
      IM_errloc (AT, "option -a is being misused.");
    }

    if (assignmentoptions[POPULATIONASSIGNMENTASN] == 1
        && assignmentoptions[POPULATIONASSIGNMENTBAR] == 0)
    {
      assignmentoptions[POPULATIONASSIGNMENT] = 1;
      assignmentoptions[POPULATIONASSIGNMENTBF] = 1;
      assignmentoptions[PRINTSTRUCTURAMA] = 1;
    }
    else if (assignmentoptions[POPULATIONASSIGNMENTASN] == 0
             && assignmentoptions[POPULATIONASSIGNMENTBAR] == 1)
    {
      assignmentoptions[POPULATIONASSIGNMENT] = 1;
      assignmentoptions[POPULATIONASSIGNMENTRELABEL] = 1;
      assignmentoptions[POPULATIONASSIGNMENTASSIGNED] = 1;
    } 
    else if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
    {
      /* MCMC State check is allowed without assignment. */
    }
    else if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
    {
      /* Island model is allowed without assignment. */
    }
    else
    {
      IM_errloc (AT, "option -a  is being misused.");
    }
  }

  if (Xp == 0)
  {
    gbeta = 1.0;
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&gbeta, 1, MPI_FLOAT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
  }

  if (!Ip)
  {
    IM_err (IMERR_MISSINGCOMMANDINFO, " No data file given on command line");
  }
  if (!Lp && !runoptions[LOADRUN])
  {
    genealogiestosave = (int) DEFAULTNUMGENEALOGIES;
    cdurationmode = TIMESTEPS;
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&genealogiestosave, 1, MPI::INT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	
	try {
		MPI::COMM_WORLD.Bcast(&cdurationmode, 1, MPI::INT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
  }
  if (runoptions[LOADRUN] || modeloptions[NOMIGRATION])
    outputoptions[MIGRATEHIST] = 0;
  if (!Op)
    strcpy (outfilename, "outfile.txt");
#ifdef _DEBUG
  checkoutfileclosed (&outfile, outfilename);   // just make sure that outfilename does not name a file that is already opened
#endif
  if (runoptions[LOADRUN] && !Vp)
  {
    IM_err (IMERR_MISSINGCOMMANDINFO, " -r0 invoked without -v information, i.e. no base name for files containing genealogys was given on the command line");
  }
  if (!Hnp || runoptions[LOADRUN])
  {
    numchains = DEFCHAINS;      /* default value */
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&numchains, 1, MPI::INT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
  }
  else
  {
    /* 5/19/2011 JH adding thermodynamic integration */
    if (calcoptions[CALCMARGINALLIKELIHOOD])
    {
      if (runoptions[LOADRUN])
      {
        IM_err (IMERR_COMMANDLINEINCOMPAT, " Conflicting command line arguments, cannot estimate marginal likelihood in Load mode (-r0)");
      }
      heatmode = HEVEN;
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&heatmode, 1, MPI::INT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
    }
    else
    {
    if ((heatmode == HGEOMETRIC && numprocesses * numchains < 4)
          /*|| (heatmode == HTWOSTEP && numchains < 3) stopped using 8/23/2011 */ )
    {
      IM_err (IMERR_COMMANDLINEHEATINGTERMS, "too few chains specified in heating model");
    }
    if (!Hfp)
    {
      heatmode = HLINEAR;
      if (!Hap)
        hval1 = 0.05;           /* default value */
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&hval1, 1, MPI_FLOAT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
    }
    else
    {
      if (heatmode > HLINEAR)
      {
          /*if (heatmode == HTWOSTEP)  stopped using 8/23/2011
        {
          if (!Hap)
            hval1 = 0.05;
          if (!Hbp)
            hval2 = 2;
          }*/
        if (heatmode == HGEOMETRIC)
        {
          if (!Hap) {
            hval1 = 0.95;  // default value
	 #ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&hval1, 1, MPI_FLOAT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
	}
          if (!Hbp) {
            hval2 = 0.8;
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&hval2, 1, MPI_FLOAT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
	}
          // if (hval1 > 1.0)  //6/11/2010 JH  stopped this,  it turns out numbers slightly higher than 1 can be useful when the are 
            // a large number of chains
           // IM_err (IMERR_COMMANDLINEHEATINGTERMS, "ha commandline term is out of range, should be <= 1");
          if (hval1 > 1.1) //6/11/2010  JH  added this to avoid values much larger than 1
            IM_err (IMERR_COMMANDLINEHEATINGTERMS, "for geometric heating it is not useful to have the ha term be greater than 1.1");
          if (hval1 < 0.9)
            IM_err (IMERR_COMMANDLINEHEATINGTERMS, "for geometric heating it is not useful to have the ha term be less than 0.9");
          if (hval2 >= 1.0|| hval2 <= 0.0)
            IM_err (IMERR_COMMANDLINEHEATINGTERMS, "hb commandline term is out of range, should be < 1 and > 0)");
        }
      }
      else if (!Hap)
        hval1 = 0.05;           /* default value */
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&hval1, 1, MPI_FLOAT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
    }
  }
  }
  recordint = RECORDINTDEFAULT;
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&recordint, 1, MPI::INT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
  if (!Up) {
    generationtime = GENERATIONTIMEDEFAULT;
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&generationtime, 1, MPI::INT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
}
  if (!Dp) {
    savegenealogyint = SAVEGENEALOGYDEFAULT;
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&savegenealogyint, 1, MPI::INT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
	}
  if (!Zp) {
    printint = PRINTINTDEFAULT;
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&printint, 1, MPI::INT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
	}
  if (!Bp && !runoptions[LOADRUN])
  {
    IM_err (IMERR_MISSINGCOMMANDINFO,
            "No burn duration information given on command line (use -b)");
  }
  if (modeloptions[NOMIGRATION] && !Mp)
  {
    mprior = 0;
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&mprior, 1, MPI_FLOAT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
  }
  if (modeloptions[NOMIGRATION] && Mp && mprior > 0)
  {
    IM_err (IMERR_COMMANDLINEINCOMPAT, " Conflicting command line arguments, no migration set but migration prior > 0 : %lf",mprior);
  }
  if (Gp && !calcoptions[USEPRIORFILE])
  {
    calcoptions[USEPRIORFILE] = 1;
  }
  if (!Tp && !assignmentoptions[POPULATIONASSIGNMENTINFINITE]
      && npops > 1
      && !calcoptions[USEPRIORFILE])
  {
    IM_err (IMERR_MISSINGCOMMANDINFO,
            " No prior information provided for population splitting times (-t)");
  }

  if (!Mp && modeloptions[NOMIGRATION] != 1 && !calcoptions[USEPRIORFILE] && npops > 1)
  {
    IM_err (IMERR_MISSINGCOMMANDINFO,
            " No information provided for maximum value for migration parameter (-m)");
  }
  if (!Qp && !calcoptions[USEPRIORFILE])
  {
    IM_err (IMERR_MISSINGCOMMANDINFO,
            " No information provided for maximum value of 4Nu parameters");
  }

  if (runoptions[LOADMCSTATE] && !Fp)
  {
    IM_err (IMERR_MISSINGCOMMANDINFO,
            "  -r3 invoked without -f information, i.e. no filename given for markov chain state file,  for loading state of markov chain ");
  }
  if (modeloptions[PARAMETERSBYPERIOD] && modeloptions[ANCESTRALPOPSIZESSETTOLARGESTSAMPLED])
    IM_err (IMERR_COMMANDLINEINCOMPAT,  " model options in conflict: -%d with -%d",PARAMETERSBYPERIOD,ANCESTRALPOPSIZESSETTOLARGESTSAMPLED);
  if (!Hkp) // default is to swap a once step or once for every 10 chains each step whichver is greater
    //swaptries = numchains;
    swaptries = IMAX(1,numchains/10);
  else
  {
    if (numchains > 1 && swaptries > numchains * (numchains - 1) / 2)
      swaptries = numchains * (numchains - 1) / 2;
  }
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&swaptries, 1, MPI::INT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif

  if (runoptions[PRINTBURNTREND])
  {
    if (runoptions[LOADMCSTATE])
      burntrendstartdelay = 0;
    else
      burntrendstartdelay = BURNTRENDSTARTDELAYDEFAULT;
  }
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&burntrendstartdelay, 1, MPI::INT, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
  if (!Sp) {
    seed_for_ran1 = (long) time (NULL);
	
/*	try {
		MPI::COMM_WORLD.Bcast(seed_for_ran1, 1, MPI::LONG, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}*/
}
  if (strcmp (infilename, outfilename) == 0)
  {
    IM_err (IMERR_COMMANDLINEINCOMPAT, " Input and output file names are identical");
  }
	//std::cout << "Am I checking this at all??\n";
  if (cdurationmode == TIMESTEPS)
  {
    chainduration = genealogiestosave * savegenealogyint;
	//std::cout << "Broadcasting chain duration again here and it is " << chainduration << "\n";	
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&chainduration, 1, MPI::LONG, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
  }
  if (jheyoptions[WRITEMIGRATIONNAME])
  {
    migrationnamefile = fopen (migrationnamefilename, "w");
  }
  if (( modeloptions[NOMIGBETWEENNONSISTERS] ||
        modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS] ||
        modeloptions[MIGRATIONBETWEENSAMPLED] || 
        modeloptions[ADDGHOSTPOP] || 
        modeloptions[PARAMETERSBYPERIOD] ||
       /* modeloptions[EXPOMIGRATIONPRIOR] ||  */
        npops == 1 ||
        modeloptions[NOMIGRATION]) && calcoptions[USEPRIORFILE])
  {
    IM_err (IMERR_COMMANDLINEINCOMPAT,"incompatibility between a model option and the use of a file with parameter priors");
  }

  if (( modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS] ||
        npops == 1 ||
        modeloptions[NOMIGRATION]) &&  modeloptions[ONEMIGRATIONPARAMETER])
  {
    IM_err (IMERR_COMMANDLINEINCOMPAT,"incompatibility between a model option and the use of a single migration parameter");
  }
  if (calcoptions[FINDJOINTPOSTERIOR]==1 && runoptions[LOADRUN]==0)
  {
    IM_err (IMERR_COMMANDLINEINCOMPAT," finding joint posterior and tests of nested models requires L mode");
  }

  /* IMa stops its running if there is any conflicted options. */
  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
  {
    if (modeloptions[NOMIGRATION] == 1)
    {
      IM_err (IMERR_COMMANDLINEINCOMPAT, 
              "Island model must be allowed for migration events: do not use -j%d with -a%d", 
              NOMIGRATION, 
              POPULATIONASSIGNMENTINFINITE);
    }
    if ( modeloptions[NOMIGBETWEENNONSISTERS] == 1
        || modeloptions[PARAMETERSBYPERIOD] == 1)
    {
      IM_err (IMERR_COMMANDLINEINCOMPAT, 
              "Island model has only a single period: do not use -j%d, -j%d, or -j%d with -a%d", 
              NOMIGBETWEENNONSISTERS, PARAMETERSBYPERIOD, 
              POPULATIONASSIGNMENTINFINITE);
    }
    if (outputoptions[THISTDIVIDEBYPRIOR] == 1
        || outputoptions[PRINTJOINTTEST] == 1)
    {
      IM_err (IMERR_COMMANDLINEINCOMPAT, 
              "Island model has no split time: do not use -p%d or -p%d with -a%d", 
              THISTDIVIDEBYPRIOR, 
              PRINTJOINTTEST,
              POPULATIONASSIGNMENTINFINITE);
    }
    if (outputoptions[POPMIGPARAMHIST] && (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS] || modeloptions[ONEMIGRATIONPARAMETER]  ))
    {
       IM_err (IMERR_COMMANDLINEINCOMPAT, 
              "output option for 2NM, -p%d, cannot be used with model options -j%d or -j%d",POPMIGPARAMHIST,SINGLEMIGRATIONBOTHDIRECTIONS,ONEMIGRATIONPARAMETER);
    }
      
    if (Tp == 1)
    {
      IM_err (IMERR_COMMANDLINEINCOMPAT, 
              "Island model has no split time: do not use -t with -a%d", 
              POPULATIONASSIGNMENTINFINITE);
    }


  }
  if (calcoptions[USEPRIORFILE] && !Gp)
    strcpy(priorfilename,defaultpriorfilename);
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&priorfilename, 500, MPI::CHAR, curr_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	
 ///AS: Broadcast everything that was set in scan_commandline() function above
 try {
	MPI::COMM_WORLD.Bcast(outputoptions, 10, MPI::INT, curr_id);
	} catch (MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
 try {
	MPI::COMM_WORLD.Bcast(modeloptions, 10, MPI::INT, curr_id);
	} catch (MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
 try {
	MPI::COMM_WORLD.Bcast(calcoptions, 6, MPI::INT, curr_id);
	} catch (MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
	int swapbetasonly = 1;
	#ifdef MPI_ENABLED
	try {
		MPI::COMM_WORLD.Bcast(&swapbetasonly, 1, MPI::INT, curr_id);
	} catch (MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return;
	}
	#endif
///AS: Should I broadcast jheyoptions() too??
 
  return;
}                               // scan_commandline 


/*  this prints basic info to a string, fpstr, that later gets printed to the output file */
void
print_outputfile_info_string (void)
{
  fpstri = static_cast<int *> (malloc (sizeof (int)));
  *fpstri = 0;
  SP "IMa2p version 1.0 - Parallel Isolation with Migration Analysis  -  Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Arun Sethuraman 2015 \n");
  SP "Release date: %s\n\n",RELEASE_DATE);
  SP "\nINPUT AND STARTING INFORMATION \n");
  SP "================================\n");
  SP "\nCommand line string : %s \n", command_line);
  SP "  Input filename : %s \n", infilename);
  SP "  Output filename: %s \n", outfilename);
  SP "  Random number seed : %li \n", seed_for_ran1);
  SP "  Heating terms on command line : %s \n",heatingterm_str);
  SP "  Calculation options on command line : %s \n",calcoptions_str);
  SP "  Model options on command line : %s \n",modeloptions_str);
  SP "  Output options on command line : %s \n",outputoptions_str);
  if (calcoptions[DONTCALCLIKELIHOODMUTATION])
    SP "**NO DATA ** - Data likelihoods set to constant  posterior should equal prior \n");
  if (!runoptions[LOADRUN])
  {
    SP "- Run Duration - \n");
    switch (burndurationmode)
    {
    case TIMESTEPS:
      SP "     Burn period, # steps: %li \n", burnduration);
      break;
    case TIMEINF:
      SP "     Burn period, # seconds: %li (total burn duration depends on IMburn file status)\n", burnduration);
      break;
    };
    if (runoptions[PRINTBURNTREND])
    {
      SP "          -User option for printing trendline during, or at end of burnin period, invoked\n");
      SP "          -initial burn duration prior to beginning recording burntrend : %d steps\n", burntrendstartdelay);
    }

    switch (cdurationmode)
    {
    case TIMESTEPS:
      SP "     Record period, #saves: %d  #steps each: %li   total #steps: %li \n", genealogiestosave, chainduration / genealogiestosave, chainduration);
      break;
    case TIMEINF:
      SP "      Record period, # seconds per interval: %ld \n",
        chainduration);
      SP "      Run period indefinite in length - determined by user using 'IMrun' file\n");
      break;
    }
    
    

    SP "- Metropolis Coupling -\n");
    if (numprocesses > 1) {
	//SP "	Parallel Metropolis Couping implemented using %d processes \n", numprocesses);
	}
    if (numprocesses * numchains > 1)
    {
      SP "     Metropolis Coupling implemented using %d chains \n", numchains * numprocesses);
      switch (heatmode)
      {
      case HLINEAR:
        SP "     Linear Increment Model   term: %.3f\n", hval1);
        break;
     /* case HTWOSTEP:
        SP "     Twostep Increment Model   term1: %.3f  term2: %.3f\n",
          hval1, hval2);
        break; stopped using 8/23/2011 */
      case HGEOMETRIC:
        SP "     Geometric Increment Model   term1: %.3f  term2: %.3f\n",
          hval1, hval2);
        break;
      }
    }
    else
      SP "     None \n");
  }
  if (calcoptions[USEPRIORFILE])
    SP "Prior distribution terms loaded from file : %s \n",priorfilename);

  if (runoptions[LOADMCSTATE])
  {
    SP "Initial Markov chain state space loaded from file: %s\n",
      mcfreadfilename);
  }
  if (modeloptions[EXPOMIGRATIONPRIOR]==1)
  {
    SP"Exponential priors used for migration rate parameters\n");
  }
  if (runoptions[SAVEMCSTATEFILE])
  {
    SP "State of Markov chain saved to file : %s\n", mcfwritefilename);
  }
  if (!runoptions[DONTSAVEGENEALOGIES])
  {
    SP "All genealogy information used for surface estimation printed to file: %s\n", genealogyinfosavefilename);
  }
  if (outputoptions[DONTSAVEGENEALOGIES])
  {
    SP "All genealogy information used for surface estimation printed to file: %s\n", genealogyinfosavefilename);
  }
  if (outputoptions[MIGRATEHIST])
  {
    /* CR:110114.2  message text changed */
    SP "Distributions for migration event counts saved in file: %s%s\n",
        outfilename,".mpt");
  }

  if (outputoptions[PRINTTMRCA])
    SP "TMRCA  histograms printed \n");
}                               //  print_outputfile_info_string



/* AS - TODO: This function has to come under the MPI framework
 * only the head node should read the file and scan_commandline()
 * only head node needs to set up all the environment variables.
 * I am yet to figure out what this means for variables shared between processes as of Tue Mar 11 15:42:31 EDT 2014
 */
// sets up filenames and calls the main initialization functions 
void start (int argc, char *argv[], int currentid)
{
  int i;
  scan_commandline (argc, argv,currentid);
  /// AS - only head node performs scan_command_line - but this has been dealt with in scan_commandline
  ///Hereon, all functions/variables that are used have to be broadcast first
  //AS:TODO runoptions has to be broadcast
  if (runoptions[SAVEMCSTATEFILE])
  {
    //strcpy (mcfwritefilename, outfilename);
    //strcat (mcfwritefilename, ".mcf");
    sprintf(mcfwritefilename, "%s.mcf.%d", outfilename, currentid);  
//    strcat (mcfwritefilename, itoa(currentid, fn, 10));
  }
///AS:TODO runoptions has to be broadcast
  if (!runoptions[DONTSAVEGENEALOGIES])
  {
	if (currentid == 0) {
  //  strcpy (genealogyinfosavefilename, outfilename);
  //  strcat (genealogyinfosavefilename, ".ti");
  // strcat (genealogyinfosavefilename, currentid);
     sprintf(genealogyinfosavefilename, "%s.ti", outfilename);
	}
  }
///AS:TODO cdurationmode has to be broadcast
  if (cdurationmode == TIMEINF)
  {
  //  strcpy (oldoutfilename, outfilename);
  //  strcat (oldoutfilename, ".old");
  //  strcat (oldoutfilename, currentid);
      sprintf(oldoutfilename, "%s.old.%d", outfilename, currentid);
  }
  for (i = 0; i < MAXLOCI; i++)
    hilocuslike[i] = -1e20;
  print_outputfile_info_string ();
///AS:TODO this has to be taken out to mersenne twister
///AS: Also, I am setting seeds based on current process ID instead. Simpler!	
  setseeds (seed_for_ran1 * currentid);
  setlogfact ();
  // CR 120124.1 memory allocation in setup() must come after setheat()
  // AS: TODO: setheat function has to be modified such that it uses the process ID or currentid here
  setheat (hval1, hval2, heatmode, currentid);
  // AS: TODO: I think each processor should set up its own set of structures, to prevent sharing difficulties, pointer math
  // But this should also be set up such that all the processors can read the file at the same time, using MPI::File
  setup (infilename, fpstri, fpstr,priorfilename, currentid);
  // AS: TODO: Think this is fine
  checkautoc (1, 0, 0, currentid);
  // AS: TODO: Same here
  gsampinflength = calc_gsampinf_length ();
  // AS: TODO Looks like outputoptions has to be broadcast as well
  if (outputoptions[MIGRATEHIST])
  {
  //  strcpy (migplotfilename, outfilename);
  //  strcat (migplotfilename, ".mpt");
  //  strcat (migplotfilename, currentid);
      sprintf(migplotfilename, "%s.mpt.%d", outfilename, currentid);
  //AS: Why has the following code been commented out? If not for this, nummigdirs is never set!
  //AS: Correction - this is set as a global variable when calling setup in initialize file
  //AS: So wouldn't have to be broadcast unlike other parameters above
  /*  nummigdirs = 2*(npops-1)*(npops-1);
    if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])
      nummigdirs = 2*nummigrateparams;  
    else
      nummigdirs = nummigrateparams;  */
    migcount = static_cast<int **> (malloc ((nloci + (nloci > 1)) * sizeof (int *)));
  
    for (i = 0; i < nloci + (nloci > 1); i++)
    {
      migcount[i] = static_cast<int *> (malloc (nummigdirs * sizeof (int)));
    }
    migfrom = static_cast<int *> (malloc (nummigdirs * sizeof (int)));
    migto = static_cast<int *> (malloc (nummigdirs * sizeof (int)));
/* 8/26/2011 */
    for (i = 0; i < nummigdirs; i++)
    {
	//AS: TODO - migration_counts[][] has to be broadcast
      //migfrom[i]= atoi(migration_counts[0][i].str);
      migfrom[i]= atoi(strchr(migration_counts[0][i].str,'m')+1);
      migto[i] = atoi(strchr(migration_counts[0][i].str,'>')+1);
    }

/*    for (i = 0; i < nummigdirs; i++)
    {
      migfrom[i]= atoi(migration_counts_times[0][2*i].str);
      migto[i] = atoi(strchr(migration_counts_times[0][2*i].str,'>')+1);
    }  */

    } 
  // AS: TODO runoptions has to be broadcast - done
  if (runoptions[LOADMCSTATE])
  {
	//AS: previously, this was erroneously being written as outputfilename.mcf.currentid. This is wrong
	//Brought to my attention by Louis Plough on Mon Oct 19 12:18:15 EDT 2015
    sprintf(mcfreadfilename, "%s.mcf.%d", mcfreadfilename, currentid);  
    readmcf (mcfreadfilename);
  }
}                               /* start */



#define TUPDATEINC  0           //1            // do a single t parameter every TUPDATEINC steps
#define UUPDATEINC  4//9           // u parameters seem to mix well so do every 5 steps


void qupdate (int curr_id,/* int *swapper, int *swappee, std::ofstream &f1, */int swapA, int swapB)
{
  int i;
  int j, k, li, ci, ui;
  int changed;
  int qswapped = 0;
  int topolchange, tmrcachange;
  int periodpick, tupdatemethodpick;
  static int tui = 0;
  static int count_t_updatetypes = -1;
  static int NWi = -1, RY1i = -1;
  static int uui = 0;

  if (count_t_updatetypes == -1)
  {
    i = 0;
#ifdef DO_RY1UPDATE  
      i++;
      RY1i = i;
#endif
#ifdef DO_NWUPDATE
    if (modeloptions[NOMIGRATION] == 0)
    {
      i++;
      NWi = i;
    }
#endif
    count_t_updatetypes = i;
  }

  int z = whichiscoldchain();

  /* update genealogies */
  for (ci = 0; ci < numchains; ci++)
  {
    for (li = 0;li<nloci;li++)
    {
      if (ci == z)
      {
        L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_ANY].tries++;
        L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_TOPOLOGY].tries++;
        L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_TMRCA].tries++;
      }
        if (updategenealogy (ci, li, &topolchange, &tmrcachange))
        {
          if (ci == z)
          {
            L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_ANY].accp++;
            L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_TOPOLOGY].accp += (topolchange > 0);
            L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_TMRCA].accp += (tmrcachange > 0);
          }
        }
    }
  }

 /* update population assignment */
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    for (ci = 0; ci < numchains; ci++)
    {
      for (i = 0; i < snupdatei; i++)
      {
        if (assignmentoptions[POPULATIONASSIGNMENTRELABEL] == 1)
        {
          updateassignmentrelabel (ci);
        }
        if (assignmentoptions[POPULATIONASSIGNMENTBF] == 1)
        {
          updateassignmentbf (ci);
        }
      }
    }
    if (assignmentoptions[POPULATIONASSIGNMENTASSIGNED] == 1)
    {
      /* for DNA Barcoding */
      imaBarTick ();
    }
  } 

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0 && npops > 1  && tui == TUPDATEINC)
  {
    for (ci = 0; ci < numchains; ci++)
    {
      periodpick = randposint (numsplittimes);
      tupdatemethodpick = randposint (count_t_updatetypes)+1;    
      //if (modeloptions[NOMIGRATION])
//        tupdatemethodpick = 0;  // must do ry1 , NW does not work with no migration

      if (tupdatemethodpick == RY1i)
      {
#ifdef DO_RY1UPDATE
        changed = changet_RY1 (ci, periodpick);
        if (ci == z)
        {
          T[periodpick].upinf[IM_UPDATE_TIME_RY1].tries++;
          if (changed)
            T[periodpick].upinf[IM_UPDATE_TIME_RY1].accp++;
        }
#endif //DO_RY1UPDATE
      }
      if (tupdatemethodpick == NWi)
      {
#ifdef DO_NWUPDATE
        /* SANGCHUL: NW update does not see eye to eye with
         * assignment update at the moment. */
        if (assignmentoptions[POPULATIONASSIGNMENTBF] == 1)
        {
          changed = changet_RY1 (ci, periodpick);
        }
        else
        {
          changed = changet_NW (ci, periodpick);
        }
        if (ci == z)
        {
          T[periodpick].upinf[IM_UPDATE_TIME_NW].tries++;
          if (changed)
            T[periodpick].upinf[IM_UPDATE_TIME_NW].accp++;
        }
#endif /* DO_NWUPDATE */
      }
    }
    tui = 0;
  }
  else
  {
    tui++;
  }

  if (uui == UUPDATEINC)
  {
    if (nurates > 1)
      for (ci = 0; ci < numchains; ci++)
      {
        for (j = 0; j < (nurates - (nurates == 2)); j++)
        {
          ui = changeu (ci, j, &k);
          if (ci == z)

          {
            L[ul[j].l].u_rec[ul[j].u].upinf->tries++;
            L[ul[k].l].u_rec[ul[k].u].upinf->tries++;
            if (ui == 1)

            {
              L[ul[j].l].u_rec[ul[j].u].upinf->accp++;
              L[ul[k].l].u_rec[ul[k].u].upinf->accp++;
            }
          }
        }
      }
    else
    {
      if (nloci == 1 && L[0].model == HKY)      /* if there is just one HKY locus kappa needs updating on its own */
        for (ci = 0; ci < numchains; ci++)
          changekappa (ci);
    }
    uui = 0;
  }
  else
  {
    uui++;
  }
	/* If there is only 1 process */
  if (numchains > 1 && numprocesses == 1)
  {
    /* CR 110929.4 get rid of extraneous args in function call
     * although return val from swapchains() call is not used locally,it may
     * be /useful during debugging
     */
	//AS: This has to change - user can perhaps have the option of setting whether to swap pointers or heats?
	//AS: I'd assume both are fast anyway
	//AS: But in parallel, we can only swap heats, not pointers
	
  	swapbetasonly = 1;
    	qswapped = swapchains (swaptries, swapbetasonly, curr_id, heatmode);
	//f1 << "Step = " << step << "\n";
  }
	/* If we have more than 1 process, we are by default swapping only beta values
 * 	so no need to broadcast it*/
	//std::cout << "Rundone before this chain goes into swapping is " << rundone << "\n";
	#ifdef MPI_ENABLED
	if (numchains * numprocesses > 1)
	{
		//std::cout << "Attempting swap...\n";
		qswapped = swapchains_bwprocesses(/*swapper, swappee,*/curr_id, step, swaptries, swapbetasonly,
					chainduration, burnduration,/* f1, swapA, swapB*/ heatmode);
		if (swapflag == 0) {
			swapflag = 1;
		} else if (swapflag == 1) {
			swapflag = 0;
		}
	}
	#endif
	//AS: Combine all swap information before printing intervaloutput
	//
	//
 /*
	#ifdef MPI_ENABLED
	if (numprocesses > 1 && (step / ((int) printint * (int) printint)) == step && step > 0) {
		//std::cout << "Going into reduction in intervaloutput\n";
		for (int x = 0; x < numprocesses; x++) {
			for (int y = 0; y < numprocesses; y++) {
				try {
					MPI::COMM_WORLD.Reduce(&swaps_bwprocesses[x][y], &swaps_rec_bwprocesses[x][y], 1,
					MPI::INT, MPI::SUM, 0);
				} catch (MPI::Exception e) {
					std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
					MPI::COMM_WORLD.Abort(-1);
				}
			}
		}
		//AS: also have to send-receive all the swapcount variables into the swapcount_bwprocesses matrix
		if (curr_id == 0) {
			for (int x = 0; x < numchains; x++) {
				for (int y = 0; y < numchains; y++) {
					swapcount_bwprocesses[x][y] = swapcount[x][y];
				}
			}
		}	
		for (int x = 1; x < numprocesses; x++) {
			if (curr_id == x) {
				for (int y = 0; y < numchains; y++) {
					for (int z = 0; z < numchains; z++) {
						MPI::COMM_WORLD.Send(&swapcount[y][z], 1, MPI::INT, 0, 1234);
					}
				}
			}
			if (curr_id == 0) {
				for (int y = 0; y < numchains; y++) {
					for (int z = 0; z < numchains; z++) {
						MPI::COMM_WORLD.Recv(&swapcount_bwprocesses[x * numchains + y][x * numchains + z], 1, MPI::INT, x, 1234);
					}
				}
			}
		}
			
			

		}
		//std::cout << "Done with reduction in interval output\n";
	#endif
	*/
//	if (curr_id == 0) {		
//	std::cout << "Going into interval output on head node\n";
  intervaloutput (stdout, curr_id);
//	std::cout << "Done with interval output on head node\n";
//}
	#ifdef MPI_ENABLED
	if (numprocesses > 1) {
//	std::cout << "Going into barrier in intervaloutput, waiting for headnode to finish printing\n";
	MPI::COMM_WORLD.Barrier();
//	std::cout << "Done with barrier in intervaloutput!\n";
	}
	#endif
  if (step >= CHECKAUTOCWAIT)
    checkautoc (0, burndone, burnsteps, curr_id);

  return;
}                               /* qupdate */



void reset_after_burn (int current_id)
{
  int ci;
  int li, ui, i, j;
  time (&chainstarttime);
  burnsteps = step - 1;
  if (cdurationmode == TIMEINF)
  {
    time (&lasttime);
    time (&timer);
  }
  // reset acceptance rate accumulators
  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1
      || npops == 1)
  {
    /* no split times */
  }
  else
  {
    for (i = 0; i < lastperiodnumber; i++)
      for (j = 0; j < T[i].num_uptypes; j++)
        T[i].upinf[j].accp = T[i].upinf[j].tries = 0;
  }
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    for (ci = 0; ci < numchains; ci++)
    {
      for (j = 0; j < Cupinf[ci].num_uptypes; j++)
      {
        Cupinf[ci].upinf[j].accp = 0;
        Cupinf[ci].upinf[j].tries = 0;
      }
    }
  }

  for (li = 0; li < nloci; li++)
  {
    if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
    {
      for (j = 0; j < L[li].a_rec->num_uptypes; j++)
      {
        L[li].a_rec->upinf[j].accp = 0;
        L[li].a_rec->upinf[j].tries = 0;
      }
      if (assignmentoptions[POPULATIONASSIGNMENTASSIGNED] == 1)
      {
        imaBarReset ();
      }
    }
    for (ui = 0; ui < L[li].nlinked; ui++)
    {
      for (j = 0; j < L[li].u_rec[ui].num_uptypes; j++)
        L[li].u_rec[ui].upinf[j].accp = L[li].u_rec[ui].upinf[j].tries = 0;
      if (L[li].model == HKY)
        for (j = 0; j < L[li].kappa_rec->num_uptypes; j++)
          L[li].kappa_rec->upinf[j].accp = L[li].kappa_rec->upinf[j].tries =
            0;
      // A_rec not used as of sometime in 2010, A gets enough updates when updating genealogy
      if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
      {
        for (j = 0; j < L[li].A_rec[ui].num_uptypes; j++)
          L[li].A_rec[ui].upinf[j].accp = L[li].A_rec[ui].upinf[j].tries = 0;
      }
      for (j = 0; j < L[li].g_rec->num_uptypes; j++)
        L[li].g_rec->upinf[j].accp = L[li].g_rec->upinf[j].tries = 0;
    }
  }
/*updated the way checkautoc() initializes on 3/25/08 */
  checkautoc (1, burndone, burnsteps, current_id);
  if (!runoptions[DONTSAVEGENEALOGIES])
  {
	if (current_id == 0) {
    if ((genealogyinfosavefile = fopen (genealogyinfosavefilename, "w")) == NULL)
    {
      IM_err (IMERR_CREATEFILEFAIL, "Error creating file for holding genealogy information");
    }
    fprintf (genealogyinfosavefile,
             "-------------------------------------------\n\n");
    fprintf (genealogyinfosavefile, "Header for genealogy file:  %s\n\n",
             genealogyinfosavefilename);
    fprintf (genealogyinfosavefile,
             "-------------------------------------------\n\n");
    fprintf (genealogyinfosavefile, "%s\n", fpstr);
    fprintf (genealogyinfosavefile,
             "-------------------------------------------\n\n");
    fprintf (genealogyinfosavefile, "End of header for genealogy file:  %s\n\n",
             genealogyinfosavefilename);
    fprintf (genealogyinfosavefile,
             "-------------------------------------------\n\n");
    fprintf (genealogyinfosavefile, "VALUESSTART\n");
    f_close (genealogyinfosavefile);
	}
  }
  if (runoptions[SAVEMCSTATEFILE])
  {
    writemcf (mcfwritefilename);
  }
	///AS: This has to change as well - head node need not be the one that holds the cold chain!	
	///AS: so find cold chain and compute this and broadcast it
	int z = whichiscoldchain();
	  if ( z >= 0) {
		gloglikelihood = C[z]->allpcalc.pdg; 
	/*	try {
			MPI::COMM_WORLD.Bcast(&gloglikelihood, 1, MPI::DOUBLE, 0);
		} catch (MPI::Exception e) {
			std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
			MPI::COMM_WORLD.Abort(-1);
			return;
		}*/
			
	}
}                               /* reset_after_burn() */

void output_burntrendfile (int current_id)
{
  FILE *burntrendfile;
  char burntrendfilename[FNSIZE];
  assert (runoptions[PRINTBURNTREND]);
  //strcpy (burntrendfilename, outfilename);
  //strcat (burntrendfilename, ".burntrend.out");
	
  sprintf(burntrendfilename, "%s.burntrend.out.%d", outfilename, current_id);  
  if ((burntrendfile = fopen (burntrendfilename, "w")) == NULL)
  {
    IM_err (IMERR_CREATEFILEFAIL,
            "Error opening burntrend output file for writing");
  }
  printf ("\n\n========================\n");
  printf ("Printing Burn Trend File\n");
  printf ("========================\n\n");
  fprintf (burntrendfile,
           "Plots of Runtime Information and Parameter Trends during Burnin \n");
  fprintf (burntrendfile, "========================================\n");
  intervaloutput (burntrendfile, current_id);
  fprintf (burntrendfile, "========================================\n\n");
  fprintf (burntrendfile, "Current Step #: %d \n\n", step);
  if (trendspot > 1 && current_id == 0)
    callasciitrend (burntrendfile);
  else if (trendspot <= 1 && current_id == 0) {
    fprintf(burntrendfile, "burn period too short to plot trends,  trend recording begins at step %d \n",burntrendstartdelay); 
  }
  if (current_id > 0) {
	fprintf(burntrendfile, "Check burn trend file on head processor for trends!\n\n");
  }
  fclose (burntrendfile);
}                               // output_burntrendfile

#define CHECKINTERVALSTEPS 1000

/* run() determines the status of a run in terms of whether it is in the burnin phase or not
  and of whether the run should keep going. 
  
  if the burnin period has just ended,  some work is done

  run()  also  checks to see if it is time to print an output file */

int run (int current_id)
{
  static int checkinterval = 0;
  static int printburntrendstep = 0, burnrecordi = 1;
  int tempburndone;
  char ch;
//std::cout << "here burndone is " << burndone << "\n";
  if (burndone)
  {
    switch (cdurationmode)
    {
    case TIMESTEPS:{
//	std::cout << "Chain ID " << current_id << "Returning " << (genealogiessaved < genealogiestosave) << "\n";
    return (step < (chainduration + burnsteps)); ///AS: Is this correct? 
     /*if (rundone == 0) {
	return (1);
	} else {
	return (0);
	}*/
    // 	return (genealogiessaved < genealogiestosave);
	///AS: Ideally I want to stop running once all the necessary genealogies have been saved
      break;
	}
    case TIMEINF:
      if (checkinterval < CHECKINTERVALSTEPS)
      {
        checkinterval++;
        return (1);
      }
      else
      {
	//std::cout << "Node " << current_id << " has now entered this printing stuff\n";

        checkinterval = 0;
        if (maxedoutgenealogysave)
	{
		return(0);
	}
        time (&timer);
        if ((timer - lasttime) > chainduration)
        {
          if ((checkdonefile = fopen ("IMrun", "r")) == NULL)
          {
		std::cout << "IMrun is null so returning 0\n";
            return (0);
          }
          else
          {
            ch = (char) getc (checkdonefile);
            f_close (checkdonefile);
            if ((char) toupper (ch) != 'Y')
            {
		std::cout << "IMrun is present, but not Y inside so returning 0\n";
              return (0);
            }
            else
            {
		//if (current_id == 0) {
		//	std::cout << "Head node is here!!\n";
		//}
	//AS: combine all the swap variables before calling printoutput
	/*#ifdef MPI_ENABLED
	if (numprocesses > 1) {
		MPI::COMM_WORLD.Barrier();
	}
	//std::cout << "Going into reduction of swapcount variables on processor " << current_id << "\n";
	if (numprocesses > 1) {
		for (int x = 0; x < numprocesses; x++) {
			for (int y = 0; y < numprocesses; y++) {
				try {
					MPI::COMM_WORLD.Reduce(&swaps_bwprocesses[x][y], &swaps_rec_bwprocesses[x][y], 1,
					MPI::INT, MPI::SUM, 0);
				} catch (MPI::Exception e) {
					std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
					MPI::COMM_WORLD.Abort(-1);
				}
			}
		}
		
		for (int x = 0; x < numprocesses * numchains; x++) {
			for (int y = 0; y < numprocesses * numchains; y++) {
				try {
					MPI::COMM_WORLD.Reduce(&tempbasedswapcount[x][y], &tempbased_rec_swapcount[x][y], 1,
								MPI::INT, MPI::SUM, 0);
				} catch (MPI::Exception e) {
					std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
					MPI::COMM_WORLD.Abort(-1);
				}
			}
		}
		
		//AS: also have to send-receive all the swapcount variables into the swapcount_bwprocesses matrix
		if (current_id == 0) {
			for (int x = 0; x < numchains; x++) {
				for (int y = 0; y < numchains; y++) {
					swapcount_bwprocesses[x][y] = swapcount[x][y];
				}
			}
		}	
		for (int x = 1; x < numprocesses; x++) {
			if (current_id == x) {
				for (int y = 0; y < numchains; y++) {
					for (int z = 0; z < numchains; z++) {
						MPI::COMM_WORLD.Send(&swapcount[y][z], 1, MPI::INT, 0, 1234);
					}
				}
			}
			if (current_id == 0) {
				for (int y = 0; y < numchains; y++) {
					for (int z = 0; z < numchains; z++) {
						MPI::COMM_WORLD.Recv(&swapcount_bwprocesses[x * numchains + y][x * numchains + z], 1, MPI::INT, x, 1234);
					}
				}
			}
		}
			
			

		}
	//std::cout << "Done with reduction of swapcount variables\n";
	#endif*/

		if (current_id == 0) {
		//std::cout << "Going into save genealogies on head node!\n";
              if (!runoptions[DONTSAVEGENEALOGIES] && genealogiessaved > 0)
              {
                savegenealogyfile (genealogyinfosavefilename, genealogyinfosavefile, &lastgenealogysaved, gsampinflength);
              }
		}
              printoutput (current_id);
		//std::cout << "Done with printoutput on head node\n";
	
              //time (&lasttime); // start the clock again 
		#ifdef MPI_ENABLED
		if (numprocesses > 1) {
		MPI::COMM_WORLD.Barrier();
		}
		#endif
		time (&lasttime);
		//std::cout << "Out of barrier after printing stuff!\n";
            }
          }
        }
        return (1);
      }
      break;
    default:
	//std::cout << "Default returning 0\n";
      return (0);
      break;
	
    }
  }
  else
  {
    tempburndone = 0;
    switch (burndurationmode)
    {
    case TIMESTEPS:
      tempburndone = (step > burnduration);
      break;
    case TIMEINF:
      if (checkinterval < CHECKINTERVALSTEPS)
      {
        checkinterval++;
      }
      else
      {
        checkinterval = 0;
        time (&timer);
        tempburndone = (timer - lasttime) > burnduration;
      }
      break;
    default:
      return (0);
      break;
    }
    // plot burn trend if called for 
    if (runoptions[PRINTBURNTREND] && (tempburndone || step >= burntrendstartdelay))
    {
      if (burnrecordi == recordint)
      {
        trendrecord (-1, current_id);       // record values of parameters that are in mcmc
        burnrecordi = 1;
      }
      else
      {
        burnrecordi++;
      }
      if (tempburndone)
      {
	//AS: combine all swap variables before calling output_burntrendfile	
	/*#ifdef MPI_ENABLED
	if (numprocesses > 1) {
		for (int x = 0; x < numprocesses; x++) {
			for (int y = 0; y < numprocesses; y++) {
				try {
					MPI::COMM_WORLD.Reduce(&swaps_bwprocesses[x][y], &swaps_rec_bwprocesses[x][y], 1,
					MPI::INT, MPI::SUM, 0);
				} catch (MPI::Exception e) {
					std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
					MPI::COMM_WORLD.Abort(-1);
				}
			}
		}
		
		for (int x = 0; x < numprocesses * numchains; x++) {
			for (int y = 0; y < numprocesses * numchains; y++) {
				try {
					MPI::COMM_WORLD.Reduce(&tempbasedswapcount[x][y], &tempbased_rec_swapcount[x][y], 1,
								MPI::INT, MPI::SUM, 0);
				} catch (MPI::Exception e) {
					std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
					MPI::COMM_WORLD.Abort(-1);
				}
			}
		}
		
		//AS: also have to send-receive all the swapcount variables into the swapcount_bwprocesses matrix
		if (current_id == 0) {
			for (int x = 0; x < numchains; x++) {
				for (int y = 0; y < numchains; y++) {
					swapcount_bwprocesses[x][y] = swapcount[x][y];
				}
			}
		}	
		for (int x = 1; x < numprocesses; x++) {
			if (current_id == x) {
				for (int y = 0; y < numchains; y++) {
					for (int z = 0; z < numchains; z++) {
						MPI::COMM_WORLD.Send(&swapcount[y][z], 1, MPI::INT, 0, 1234);
					}
				}
			}
			if (current_id == 0) {
				for (int y = 0; y < numchains; y++) {
					for (int z = 0; z < numchains; z++) {
						MPI::COMM_WORLD.Recv(&swapcount_bwprocesses[x * numchains + y][x * numchains + z], 1, MPI::INT, x, 1234);
					}
				}
			}
		}
			
			

		}
	#endif*/
        output_burntrendfile (current_id);

        printburntrendstep = 1;
        /* now check to see if IMburn file is present */
        if (burndurationmode == TIMEINF)
        {
          if ((checkdonefile = fopen ("IMburn", "r")) != NULL)
          {
            ch = (char) getc (checkdonefile);
            f_close (checkdonefile);
            if ((char) toupper (ch) == 'Y')
            {
              tempburndone = 0;
              time (&lasttime); // start the clock again 
            }
          }
        }
      }
      else
      {
        printburntrendstep++;
      }
    }
    if (tempburndone)
    {
      burndone = 1;
      reset_after_burn (current_id);
    }
//std::cout << "Returning 1\n";
    return (1);
  }
}                               /* run */

#define  TRENDLASTPT (TRENDDIM - 1)
void inctrend (int m, int t, struct value_record *v, double newval)
{
  int j;
  for (j = m; j < TRENDLASTPT; j++)
  {
    /* assert (v->trend[j + 1] != 0); This asserts sometimes. */
    v->trend[j] = v->trend[j + 1];
  }
  v->trend[t] = newval;
  //assert(v != 0);
}

void trend_reset (struct value_record *v, int nv)
{
  int i, j;
  for (i = 0; i < nv; i++)
    for (j = 0; j < TRENDDIM; j++)
      v[i].trend[j] = 0;
}                               // trend_reset


/* Using trendrecord()
This function records trend lines
It works on instances of struct value_record 
the value_record must first be initialized (probably in initialize.c) 
Code for a particular value_record or array of value_records 
can be placed in trendrecord at two places:
  in the "if (burndone && reset == 0) " section
  and in the "if (recordinc == recordtrendinc)" section


explanation for how trendline data are accumulated:
 - movespot is the point at which values are deleted from the array, 
 - each new value is added to the end of the array (i.e. at trendspot)
 - all values from one position past movespot up to the end of the array are moved down one position
 - this erases the value at movespot and makes room for a new value at the end. 
 -each time the replacement point (movespot) reaches the end of the array the time 
	period between additions doubles 
values are added more slowly as the run proceeds.  
- this is because the time period doubles when movespot reaches the end 
- the values to the left of movespot have half the density (in time) of those to the right 
*/
//AS: adding currentid, so I can start recording values correctly on the head node
//2/7/2014

void trendrecord (int loadarrayj, int currentid)
{
  static int /*trendspot = 0,*/ recordtrendinc = 1, recordinc = 0, movespot =
    TRENDLASTPT;
  static int init = 1, reset = 0;
  int j, k, li, ui;

  if (burndone && reset == 0)   // reset all trend-related values after burnin
  {
    init = 1;
    reset = 1;
    if (lpgpd_v->do_trend)
    {
      trend_reset (lpgpd_v, 1);
    }
    for (li = 0; li < nloci; li++)
    {
      if (runoptions[LOADRUN] == 0)
      {
        if (L[li].g_rec->v->do_trend)
        {
          trend_reset (L[li].g_rec->v, 1);
        }
        for (ui = 0; ui < L[li].nlinked; ui++)
        {
          if (L[li].u_rec[ui].v->do_trend)
          {
            trend_reset (L[li].u_rec[ui].v, 1);
          }
        }
        if (L[li].model == HKY)
        {
          if (L[li].kappa_rec->v->do_trend)
          {
            trend_reset (L[li].kappa_rec->v, 1);
          }
        }
        if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
        {
          if (L[li].a_rec->v->do_trend == 1)
          {
            trend_reset (L[li].a_rec->v, 1);
          }
        }
      }
    }
    if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1
        || npops == 1)
    {
      /* no split time */
    }
    else
    {
      for (k = 0; k < lastperiodnumber; k++)
        if (T[k].v->do_trend)
          trend_reset (T[k].v, 1);
    }

/* ADD ADDITONAL trend_reset() calls here */
  }
  if (init == 1)
  {
    trendspot = 0;
    recordtrendinc = 1;
    recordinc = 0;
    movespot = TRENDLASTPT;
    init = 0;
  }
  recordinc++;
  if (recordinc == recordtrendinc)
  {
    if (runoptions[LOADRUN])
    {

      if (lpgpd_v->do_trend)
        inctrend (movespot, trendspot, lpgpd_v,
                  gsampinf[loadarrayj][gsamp_pdgp] +
                  gsampinf[loadarrayj][gsamp_probgp]);

      if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0
          && npops > 1)
      {
        for (j = 0; j < lastperiodnumber; j++)
          if (T[j].v->do_trend)
            inctrend (movespot, trendspot, T[j].v,
                      gsampinf[loadarrayj][gsamp_tp + j]);
      }
    }
    else
    {
      if (lpgpd_v->do_trend)
      {
	int z = whichiscoldchain();
	//Have to send/receive if the cold chain exists on some other processor
	//AS: 7/2/2014
	if (z >= 0 && currentid == 0) {
	        inctrend (movespot, trendspot, lpgpd_v,
                  C[z]->allpcalc.probg + C[z]->allpcalc.pdg);
	}
	#ifdef MPI_ENABLED
	if (z < 0 && currentid == 0) {
		double probrec = 0.0;
		MPI::COMM_WORLD.Recv(&probrec, 1, MPI::DOUBLE, MPI_ANY_SOURCE, 1313);
		inctrend (movespot, trendspot, lpgpd_v, probrec);
	}
	if (z >=0 && currentid != 0) {
		double probrec = C[z]->allpcalc.probg + C[z]->allpcalc.pdg;
		MPI::COMM_WORLD.Send(&probrec, 1, MPI::DOUBLE, 0, 1313);
	}
	#endif



	}

      if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1
          || npops == 1)
      {
        /* no split time */
      }
      else
      {
        for (j = 0; j < lastperiodnumber; j++) {
          if (T[j].v->do_trend) {
		int z = whichiscoldchain();
		if (z >= 0 && currentid == 0) {
	            inctrend (movespot, trendspot, T[j].v, C[z]->tvals[j]);
		}
		#ifdef MPI_ENABLED
		if (z < 0 && currentid == 0) {
			double probrec = 0.0;
			MPI::COMM_WORLD.Recv(&probrec, 1, MPI::DOUBLE, MPI_ANY_SOURCE, 1717);
			inctrend (movespot, trendspot, T[j].v, probrec);
		}
		if (z >= 0 && currentid != 0) {
			double probrec = C[z]->tvals[j];
			MPI::COMM_WORLD.Send(&probrec, 1, MPI::DOUBLE, 0, 1717);
		}
		#endif
		}
	}
      	}
    }
    for (li = 0; li < nloci; li++)
    {
      if (runoptions[LOADRUN] == 0)
      {
        for (ui = 0; ui < L[li].nlinked; ui++)
          if (L[li].u_rec[ui].v->do_trend) {
		int z = whichiscoldchain();
		if (z >= 0 && currentid == 0)  {
            inctrend (movespot, trendspot, L[li].u_rec[ui].v,
                      C[z]->G[li].uvals[ui]);
		}
		#ifdef MPI_ENABLED
		if (z < 0 && currentid == 0) {
			double probrec = 0.0;
			MPI::COMM_WORLD.Recv(&probrec, 1, MPI::DOUBLE, MPI_ANY_SOURCE, 1414);
			inctrend (movespot, trendspot, L[li].u_rec[ui].v, probrec);
		}
		if (z >= 0 && currentid != 0) {
			double probrec = C[z]->G[li].uvals[ui];
			MPI::COMM_WORLD.Send(&probrec, 1, MPI::DOUBLE, 0, 1414);
		}
		#endif
		}
	  }
        if (L[li].model == HKY)
          if (L[li].model == HKY) {
		int z = whichiscoldchain();
		if (z >= 0 && currentid == 0) {
	            inctrend (movespot, trendspot, L[li].kappa_rec->v,
                      C[z]->G[li].kappaval);
		}
		#ifdef MPI_ENABLED
		if (z < 0 && currentid == 0) {
			double probrec = 0.0;
			MPI::COMM_WORLD.Recv(&probrec, 1, MPI::DOUBLE, MPI_ANY_SOURCE, 1515);
			inctrend (movespot, trendspot, L[li].kappa_rec->v, probrec);
		}
		if (z >=0 && currentid != 0) {
			double probrec = C[z]->G[li].kappaval;
			MPI::COMM_WORLD.Send(&probrec, 1, MPI::DOUBLE, 0, 1515);
		}
		#endif

	  }
        if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
        {
          inctrend (movespot, trendspot, L[li].a_rec->v, 
                    imaAsnValue (0, li));
        }
      }
/* ADD ADDITONAL inctrend() calls here */
    if (movespot == TRENDLASTPT && trendspot == TRENDLASTPT)
    {
      movespot = 0;
      recordtrendinc *= 2;
    }
    else
    {
      movespot += (movespot < TRENDLASTPT);
    }
    trendspot += (trendspot < TRENDLASTPT);
    recordinc = 0;
  }
  trenddoublepoint = movespot;
}                               /* trendrecord */

/* calculates the bin number of an xy array of a value_record that a value falls in,  increments the count in that bin */
//AS: I need to change this for MPI version - every time I record something on the cold chain process,
//it also has to be updated on the head node, so printing can happen correctly

void recordval (struct value_record *v, double val)
{

  int k;
  double logval;
  if (v->do_logplot)
  {
    assert (!(val < 0.0));
    logval = log (val);
    k =
      (int) (GRIDSIZE * (logval + v->plotrange.max) /
             (2.0 * v->plotrange.max));
  }
  else
  {
    k = (int) (GRIDSIZE * val / v->plotrange.max);
  }
  
  if (k < 0)
  {
    v->beforemin++;
  }
  else if (k >= GRIDSIZE)
  {
    v->aftermax++;
  }
  else
  {
    assert (!(k < 0));          // FIXME: it's been crashing
    v->xy[k].y++;
  }
 
  return;
}

/* this is used to record the names of actual loci and gene copies that migrate, from to to
write these names to a file */
void record_migration_names(void)
 {
  int i, j, li;
  int from, to;
  struct edge *gtree;
  from = migrationnamefrom;
  to = migrationnameto;
  int z = whichiscoldchain();
  if (z >= 0) {
  for (li = 0; li < nloci; li++)
  {
    gtree = C[z]->G[li].gtree;
    for (i = 0; i < L[li].numlines; i++)
    {
      j = 0;
      while (gtree[i].mig[j].mt > 0)
      {
        if (from == nowedgepop (0, &gtree[i], gtree[i].mig[j].mt) && to == C[z]->G[li].gtree[i].mig[j].mp)
        {
          fprintf(migrationnamefile,"%s ",L[li].name);
          if (i< L[li].numgenes)
            fprintf(migrationnamefile,"%s ",L[li].gNames[i]);
          else
            fprintf(migrationnamefile,"internal ");
        }
        j++;
      }
    }
  }
  fprintf(migrationnamefile,"\n");
  } else {
	return;
 }
 } // record_migration_names



/* INSTRUCTIUONS to record a numerical value from the markov chain:
----------------------------------------------------
this works on instances of struct value_record

the value_record must be initiatlized (e.g. in initialize.c,  see e.g. init_g_rec)
this includes a call to init_value_record()

insert line(s) code into record()  below,  to make a call to recordval() 

*/

void record_migrations (int z, int currentid)
{
  int i, j, k, li, from, to, foundparam;
  struct edge *gtree;
	if (z >= 0) {

  for (j = 0; j < nloci + (nloci > 1); j++)
    for (i = 0; i < nummigdirs; i++)
    {
      migcount[j][i] = 0;
    }
  // 4/1/11  stopped recording distribution of migration times 
  for (li = 0; li < nloci; li++)
  {
    gtree = C[z]->G[li].gtree;
    for (i = 0; i < L[li].numlines; i++)
    {
      j = 0;
      while (gtree[i].mig[j].mt > 0)
      {
        from = nowedgepop (z, &gtree[i], gtree[i].mig[j].mt);
        to = C[z]->G[li].gtree[i].mig[j].mp;
        k = 0;
        foundparam = 0;
        while (k < nummigdirs && foundparam == 0)
        {
          if (from == migfrom[k] && to == migto[k])
            foundparam = 1;
          else
            k++;
        }
	//As: debug only
	//std::cout << "k is " << k << " and nummigdirs is " << nummigdirs << "\n";
        assert(k<nummigdirs);
        if (nloci == 1)
        {
          migcount[0][k]++;
          /* CR:110114.2  1 line removed */
        }
        else
        {
          migcount[li + 1][k]++;
          migcount[0][k]++;
          /* CR:110114.2  2 lines removed */
        }
        j++;
      }
    }
  }
  for (i = 0; i < nloci + (nloci > 1); i++)
    for (k = 0; k < nummigdirs; k++)
    {
      /* 8/26/2011 */
      recordval (&migration_counts [i][k], (double) migcount[i][k]);
      //recordval (&migration_counts_times[i][2 * k + 1], (double) migcount[i][k]);
    }
  }
}                               //record_migrations


void record (int currentid)
{
  int j, li, ui;
  struct genealogy *G;
  /* 5/19/2011 JH adding thermodynamic integration */
  int z = whichiscoldchain();

// if (z >= 0) {
  if (outputoptions[PRINTTMRCA])
    for (li = 0; li < nloci; li++)
    {
	//AS: 6/16/2014 - recording should happen on head node for me to print values correctly
	if (z >= 0 && currentid == 0) {
	      recordval (L[li].g_rec->v, C[z]->G[li].roottime);
	}
	#ifdef MPI_ENABLED
	if (z >=0 && currentid !=0) {
		double rtime = C[z]->G[li].roottime;
		MPI::COMM_WORLD.Send(&rtime, 1, MPI::DOUBLE, 0, 2323);
	}
	if (z < 0 && currentid == 0) {
		double rtime = 0.0;
		MPI::COMM_WORLD.Recv(&rtime, 1, MPI::DOUBLE, MPI_ANY_SOURCE, 2323);
		recordval (L[li].g_rec->v, rtime);
	}
	#endif
    }
  for (li = 0; li < nloci; li++)
  {
   if (z >= 0) {
    G = &(C[z]->G[li]);
   }


    for (ui = 0; ui < L[li].nlinked; ui++)
    {
	if (z >= 0 && currentid == 0) {
	      recordval (L[li].u_rec[ui].v, G->uvals[ui]);
	}
	#ifdef MPI_ENABLED
	if (z >=0 && currentid !=0) {
		double uvals = G->uvals[ui];
		MPI::COMM_WORLD.Send(&uvals, 1, MPI::DOUBLE, 0, 2424);
	}
	if (z < 0 && currentid == 0) {
		double uvals = 0.0;
		MPI::COMM_WORLD.Recv(&uvals, 1, MPI::DOUBLE, MPI_ANY_SOURCE, 2424);
		recordval (L[li].u_rec[ui].v, uvals);
	}
	#endif
    }

    if (L[li].model == HKY)
    {
	if (z >= 0 && currentid == 0) {
	      recordval (L[li].kappa_rec->v, G->kappaval);
	}
	#ifdef MPI_ENABLED
	if (numprocesses > 1) {
		if (z >= 0 && currentid != 0) {
			double kval = G->kappaval;
			MPI::COMM_WORLD.Send(&kval, 1, MPI::DOUBLE, 0, 2345);
		}
		if (z < 0 && currentid == 0) {
			double kval = 0.0;
			MPI::COMM_WORLD.Recv(&kval, 1, MPI::DOUBLE, MPI_ANY_SOURCE, 2345);
			recordval (L[li].kappa_rec->v, kval);
		}
	}
	#endif

    }
  }

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1
      || npops == 1)
  {
    /* we could have lastperiodnumber be 0 */
  }
  else
  {
    for (j = 0; j < lastperiodnumber; j++)
    {
	if ( z >= 0) {
	      assert (C[z]->tvals[j] > T[j].pr.min && C[z]->tvals[j] < T[j].pr.max);
	}
	if (z >= 0 && currentid == 0) {
	      recordval (T[j].v, C[z]->tvals[j]);
	}
	#ifdef MPI_ENABLED
	if (z >= 0 && currentid != 0) {
		double tvals = C[z]->tvals[j];
		MPI::COMM_WORLD.Send(&tvals, 1, MPI::DOUBLE, 0, 2525);
	}
	if (z < 0 && currentid == 0) {
		double tvals = 0.0;
		MPI::COMM_WORLD.Recv(&tvals, 1, MPI::DOUBLE, MPI_ANY_SOURCE, 2525);
		recordval (T[j].v, tvals);
	}
	#endif
    }
  }
  if (outputoptions[MIGRATEHIST])
  {
    record_migrations (z, currentid);
  }
 
  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
  {
    /* no split time */
  }
  else
  {
    if (npops >= 3 && npops <= 5 && outputoptions[PRINTJOINTTEST])
      setup_multi_t_arrays (z);
  }
//  }
  if (!outputoptions[DONTPRINTASCIITREND])
  {
    trendrecord (-1, currentid);
  }
   /* 5/19/2011 JH adding thermodynamic integration */
  if (calcoptions[CALCMARGINALLIKELIHOOD])
  {
    summarginlikecalc();
  } 
  
  return;
}                               /* record */


int
whichiscoldchain (void)
{
	int which = -1;
	for (int i = 0; i < numchains; i++) {
		if (beta[i] == 1.0) {
			which = i;
		}
	}
	//else return a -ve flag
	return which;
}

void savegenealogyinfo (int current_id)        // use floats to save space
{
  int i;
  #ifdef MPI_ENABLED
  MPI::Status status;
  MPI::Request request[1];
 #endif
  float *gsampinflocal;
//  request = new MPI::Request[1];

///Allocate memory for the full gsampinf only if I am on the head node
//else, I only need a local structure to be filled up, then I'll send that to the head node	
  if (genealogiessaved == 0)
  {
    if (cdurationmode == TIMESTEPS)
    {
	///AS: This is a problem - each process is going to now have its own set of genealogies
	//that are going to be saved, and the total number of genealogies should be = genealogiestosave
	//so we need to dynamically allocate this. This current method is wasteful.
	//Better way to do is to create a 2D vector gsampinf, keep push_back() into it dynamically
	//But keeping this for now as on Mon Mar 31 16:14:33 EDT 2014
      gsampinf = static_cast<float **> (malloc (genealogiestosave * sizeof (float *)));
      for (i = 0; i < genealogiestosave; i++)
        gsampinf[i] = static_cast<float *> (malloc (gsampinflength * sizeof (float)));
      memforgenealogiessaved = genealogiestosave;
    }
    else                        /* cdurationmode == TIMEINF */
    {
      gsampinf = static_cast<float **> (malloc (MAXGENEALOGIESTOSAVE * sizeof (float *)));
      for (i = 0; i < MAXGENEALOGIESTOSAVE; i++)
        gsampinf[i] = static_cast<float *> (malloc (gsampinflength * sizeof (float)));
      memforgenealogiessaved = MAXGENEALOGIESTOSAVE;
    }
  }
///AS: creating local gsampinf float array instead of length gsampinflength
//AS: ideally should do a check for cduration mode here as well and allocate gsampinflocal accordingly
	if (current_id != 0) {	
		gsampinflocal = static_cast<float *> (malloc (gsampinflength * sizeof (float)));
	}

  if (genealogiessaved >= (MAXGENEALOGIESTOSAVE - 1) && cdurationmode == TIMEINF)
  {
    printf (" maximum possible genealogies saved \n");
    maxedoutgenealogysave = 1;
	#ifdef MPI_ENABLED	
	try {
		MPI::COMM_WORLD.Bcast(&maxedoutgenealogysave, 1, MPI_INT, current_id);
	 
	} catch(MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
	}
	#endif

  }
  else
  {
	//save only cold chain info - this is going to be floating around, till finally
	//we can collate it into a single gsampinf for L mode
	//AS: But then again autocorrelations, etc are computed in runtime in M mode
	//so instead, I should broadcast and accumulate this on the head node.
	int z = whichiscoldchain();
	//AS: Debug only
	//std::cout << "Cold chain here is " << z << " on chain "<< current_id <<"\n";
	if (current_id == 0) {
		if ( z >= 0) {
		//	std::cout << "Saving genealogy here\n";
		/*	
		std::cout << "Here 1???";
		for (int q = 0; q < gsampinflength; q++) {
			std::cout << gsampinf[0][q] << "\t";
		}*/
		//std::cout << "\n";

		    savegsampinf (gsampinf[genealogiessaved], z);
			return;

		}
	} else if (current_id != 0) {
		//std::cout << "Here 2???";
		if ( z >= 0) {
			savegsampinf (gsampinflocal, z);
		//AS: Debug only
		//	std::cout << "Before sending gsampinflocal\n";
		//	for (int q = 0; q < gsampinflength; q++) {
		//		std::cout << gsampinflocal[q] << "\t";
		//	}
		//	std::cout << "\n";
		}
	}
	#ifdef MPI_ENABLED
	///if i am not on head node, send this info to head node
	if (current_id != 0 && numprocesses > 1 && z >= 0) {
		//std::cout << "Waiting to send...\n";
		//send my id first so head node knows where to receive from
		try {
			int tempcurrid = current_id;
			request[0] = MPI::COMM_WORLD.Isend(&tempcurrid, 1, MPI::INT, 0, 1212);
		} catch (MPI::Exception e) {
			std::cout << "Error in sending my id!\n";
			MPI::COMM_WORLD.Abort(-1);
		}
		try {
			MPI::Request::Waitall(1, &request[0]);
		} catch (MPI::Exception e) {
			std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
			std::cout << "Error sending my id!\n";
			MPI::COMM_WORLD.Abort(-1);
		}
		for (int v = 0; v < gsampinflength; v++) {
			try {
				request[0] = MPI::COMM_WORLD.Isend(&gsampinflocal[v], 1, MPI::FLOAT, 0, v*2);
			} catch (MPI::Exception e) {
				std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
				std::cout << "Error in sending genealogy\n";
				MPI::COMM_WORLD.Abort(-1);
			//genealogiessaved++;
			//	XFREE(gsampinflocal);
				return;
			}
			try {
				MPI::Request::Waitall(1, &request[0]);
			} catch (MPI::Exception e) {
				std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
				MPI::COMM_WORLD.Abort(-1);
				return;
			}

		}
		XFREE(gsampinflocal);
		
	} else if (current_id == 0 && numprocesses > 1 && z < 0) {
		//std::cout << "Waiting to receive...\n";
		//receive info from the other proces first about where to receive from
		int tempcurrid = 0;
		try {
			request[0] = MPI::COMM_WORLD.Irecv(&tempcurrid, 1, MPI::INT, MPI_ANY_SOURCE, 1212);
		} catch (MPI::Exception e) {
			std::cout << "Error in receiving the other person's id!\n";
			MPI::COMM_WORLD.Abort(-1);
		}	
		try {
			MPI::Request::Waitall(1, &request[0]);
		} catch (MPI::Exception e) {
			std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
			MPI::COMM_WORLD.Abort(-1);
			return;
		}
		for (int v = 0; v < gsampinflength; v++) {
			try {
				request[0] = MPI::COMM_WORLD.Irecv(&gsampinf[genealogiessaved][v], 1, MPI::FLOAT, tempcurrid, v*2);
			} catch (MPI::Exception e) {
				std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
				std::cout << "Error in receving genealogy\n";
				MPI::COMM_WORLD.Abort(-1);
				return;
			}
			try {
				MPI::Request::Waitall(1, &request[0]);
			} catch (MPI::Exception e) {
				std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
				return;
			}
		}
		//AS: Debug only
		//std::cout << "After receiving gsampinflocal\n";
		//for (int q = 0; q < gsampinflength; q++) {
		//	std::cout << gsampinf[genealogiessaved][q] << "\t";
		//}
		//std::cout << "\n";
	}
	#endif
	
	return;
  }
//	delete[] request;
}                               /* savegenealogyinfo */


/* SANGCHUL: This could be simplified considering its function.
 * */
void loadgenealogyvalues (void)
{
  char filenamewildcard[FNSIZE];
  char *ctp, *textline, *dataline, *c, tempc;
             /* CR 111019.1 increase from 12 to 20 handles larger data values */
  int charspervalue = 20; 
  FILE *sfile;
  int i, j, numgenealogies, totalnumgenealogies;
  /* int filefound, nofile; */
  char *defaultdir; 
  struct dirent *dir_entry;
  DIR *dp;
  int numfiles = 0;
  int numtoload, loadall, loaded, notloaded;
  size_t len_base;
  size_t l2;
  size_t len_defaultdir;
  char *basename;
  char *treefilename;

  float load_fraction;
  if (strlen (loadfilebase))
  {
    strcpy (filenamewildcard, loadfilebase);
  }
  else
  {
    strcpy (filenamewildcard, outfilename);
    strtrunc (filenamewildcard, (char) '.');
    strtrunc (filenamewildcard, (char) '-');
  }
//AS: adding a .0 to the end, since the head node is the one that collects all genealogies sampled
  strcat (filenamewildcard, "*.ti");
  SP "\nLOAD TREES (L) MODE INFORMATION\n");
  SP "============================================================================\n");
  SP "  Base filename for loading files with sampled genealogies: %s\n", filenamewildcard);
  SP "  Files loaded with sampled genealogies:\n");
  numgenealogyfiles = 0;
  textline = static_cast<char *> (malloc (300 * sizeof (char)));
  dataline = static_cast<char *> (malloc (gsampinflength * charspervalue * sizeof (char)));
  ctp = &textline[0];
  numgenealogies = totalnumgenealogies = 0;

  /* Find the default directory based on loadfilebase.
   * /this/directory/a.out -> defaultdir is /this/directory/
   * /a.out                -> /
   * a.out                 -> ./
   */
  imaDirBase (loadfilebase, &defaultdir); 
  len_defaultdir = strlen (defaultdir);
  
  if ((dp = opendir (defaultdir)) == NULL)
  {
    IM_err (IMERR_TIFILE, "cannot open directory: %s", defaultdir);
  }
  basename = strrchr (loadfilebase, '/');
  if (basename == NULL)
    {
      basename = strrchr (loadfilebase, '\\');
      if (basename == NULL)
        basename = loadfilebase;
      else
        basename++;
    }
  else
    {
      basename++;
    }
  len_base = strlen (basename);
  while ((dir_entry = readdir(dp)) != NULL)
  {
    if (!strncmp(dir_entry->d_name, basename, len_base)) 
    {
      l2 = strlen (dir_entry->d_name);
      if (!strcmp (&dir_entry->d_name[l2 - 3], ".ti"))
      {  
        numgenealogyfiles++; 
        /* We found one. */
        treefilename = static_cast<char *> (malloc ((len_defaultdir + l2 + 1) * sizeof (char)));
        sprintf (treefilename, "%s%s", defaultdir, dir_entry->d_name);

        /* Count the number of gene genealogies of the found file. */
        if ((sfile = fopen (treefilename, "r")) == NULL) 
        {
          IM_err (IMERR_TIFILE, " cannot open .ti file");
        }
        while (fgets (textline, 300, sfile)
               && strstr (textline, "VALUESSTART") == NULL && !feof (sfile));
        numgenealogies = 0;
        while ((tempc = fgetc (sfile)) != EOF)    // count lines
        {
          numgenealogies += (tempc == '\n');
        }
        if (numgenealogies < 1)
        {
          printf ("  *no genealogies loaded from file %s\n", dir_entry->d_name);
          SP "  *no genealogies loaded from file %s\n", dir_entry->d_name);
        }
        else
        {
          printf ("  loaded %d genealogies from genealogy file  %s\n", numgenealogies,
                  dir_entry->d_name);
          SP "  loaded %d genealogies from genealogy file  %s\n", numgenealogies,
            dir_entry->d_name);
        }
        fclose(sfile);
        totalnumgenealogies += numgenealogies;

        XFREE (treefilename);
      }
    }
  }
  closedir(dp);
  /* We have counted gene genealogies. */

  if (genealogiestosave > 0)
    numtoload = IMIN (totalnumgenealogies, genealogiestosave);
  else
    numtoload = totalnumgenealogies;
  numtoload = IMIN (numtoload, MAXGENEALOGIESTOSAVE);
  memforgenealogiessaved = numtoload;
  gsampinf = static_cast<float **> (malloc (numtoload * sizeof (float *)));
  loadall = numtoload >= totalnumgenealogies;
  load_fraction = (float) numtoload/ (float) totalnumgenealogies;
  loaded = 0;
  notloaded = 0;
  SP "\n  HEADER INFORMATION FROM FIRST GENEALOGY FILE\n");
  SP "  ---------------------------------------\n");
// now go through again and save the genealogies
  /* closedir (dp); */
  if ((dp = opendir (defaultdir)) == NULL)
  {
    IM_err (IMERR_TIFILE, "cannot open directory: %s", defaultdir);
  }
  numfiles = 0;

  while ((dir_entry = readdir(dp)) != NULL)
  {
    if (!strncmp(dir_entry->d_name, basename, len_base)) 
    {
      l2 = strlen (dir_entry->d_name);
      if (!strcmp (&dir_entry->d_name[l2 - 3], ".ti"))
      {  
        numfiles++;
        /* We found one. */
        treefilename = static_cast<char *> 
                (malloc ((len_defaultdir + l2 + 1) * sizeof (char)));
        sprintf (treefilename, "%s%s", defaultdir, dir_entry->d_name);

        /* Count the number of gene genealogies of the found file. */
        if ((sfile = fopen (treefilename, "r")) == NULL) 
        {
          IM_err (IMERR_TIFILE, " cannot open .ti file");
        }

        while (fgets (textline, 300, sfile) && strstr (textline, "VALUESSTART") == NULL && !feof (sfile))
          if (numfiles == 1)
          {
            SP "  ||%s", textline);
          }
        while (fgets (dataline, gsampinflength * charspervalue, sfile)
               && !feof (sfile))
        {
          if (loadall || (loaded == 0)
              ||  ((float) loaded/ (float) (loaded + notloaded)) <= load_fraction)
              //JH 7/24/09  change this to using a proportion((float) loaded / (float) notloaded) <= load_notload_ratio)
          {
            gsampinf[loaded] = static_cast<float *> 
                    (malloc (gsampinflength * sizeof (float)));
            c = dataline;
            for (i = 0; i < gsampinflength; i++)
            {
              sscanf (c, "%f", &gsampinf[loaded][i]);
              j = allwhitespace (c);
              if (j ==1 || j== -1)
                IM_err (IMERR_TIFILE, "Problem in .ti file %s,  too few values per genealogy, .ti file may have been generated with a different program",treefilename);
              c = nextwhite (c);
            }
            j = allwhitespace (c);
            if (j==0)
               IM_err (IMERR_TIFILE, "Problem in .ti file %s,  too many values per genealogy, .ti file may have been generated with a different program",treefilename);
            loaded++;
          }
          else
            notloaded++;
          /* gcounter++; */
        }
        fclose (sfile);

        XFREE (treefilename);
      }
    }
  }
  closedir(dp);

  SP "  END OF HEADER INFORMATION FROM FIRST GENEALOGY FILE\n");
  SP "  ----------------------------------------------\n\n");
  
  SP "  Number of files loaded : %d  total number of genealogies used: %d out of a total of: %d \n", numgenealogyfiles, loaded, totalnumgenealogies);
  printf("  Number of files loaded : %d  total number of genealogies used: %d out of a total of: %d \n", numgenealogyfiles, loaded, totalnumgenealogies);
  if (numgenealogies < 1)
    IM_err (IMERR_TIFILE, "  no genealogies loaded from .ti file(s)");
  /* closedir (dp); */
  genealogiessaved = loaded;
  //int z = whichiscoldchain();
//	if (z >= 0) {
  for (j = 0; j < genealogiessaved; j++)
  {
    // use full range of t, ignore t.pr.min > 0
    if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0
        && npops > 1)
    {
      for (i = 0; i < lastperiodnumber; i++)
      {
        recordval (T[i].v, gsampinf[j][gsamp_tp + i]);
      }
    }
    if (!outputoptions[DONTPRINTASCIITREND])
      trendrecord (j, 0);
  }
// }
  XFREE (textline);
  XFREE (dataline);
  XFREE (defaultdir); 
}                               /* sang chul's loadgenealogyvalues */

void printsteps (FILE * outto, double like)
{
  int ci;

  ci = 0;
  if (!burndone)
  {
    fprintf (outto, "=BURNIN-PERIOD===============================\n");
    fprintf (outto, "STEP # %d  p(D|G): %.3lf p(G): %.3lf\n", step,
             like, C[ci]->allpcalc.probg);
  }
  else
  {
    fprintf (outto, "=============================================\n");
    if (genealogiessaved > 0)
      fprintf (outto,
               "STEP # %ld # Genealogies Saved: %d p(D|G): %.1lf p(G): %.1f\n",
               step - burnsteps, genealogiessaved, like, C[ci]->allpcalc.probg);
    else
      fprintf (outto, "STEP # %ld  p(D|G): %.3lf p(G): %.3lf\n",
               step - burnsteps, like, C[ci]->allpcalc.probg);
  }
  return;
}

/* To print acceptance rates:
	reclist[] is an array of pointers to struct chainstate_record_updates_and_values
	set values of reclist[] to those structures for which you want to print acceptance rates
	call printacceptancerates ()
	*/

void callprintacceptancerates (FILE * outto, int currentid)
{
  int i, j, li;
  // length of this array must be fairly long, although it is technically possible to have MAXLOCI * MAXLINKED records,  but very unlikely
  struct chainstate_record_updates_and_values *reclist[MAXLOCI + MAXLINKED];
  struct chainstate_record_updates_and_values *reclist_saved[MAXLOCI + MAXLINKED];

// t values 
  for (i = 0; i < numsplittimes; i++)
    reclist_saved[i] = (T + i);
#ifdef DO_RY1UPDATE
#ifdef MPI_ENABLED
if (numprocesses > 1) {
int *y = static_cast<int *> (malloc (numsplittimes * sizeof (int)));
int *y_rec = static_cast<int *> (malloc (numsplittimes * sizeof (int)));
int *z = static_cast<int *> (malloc (numsplittimes * sizeof (int)));
int *z_rec = static_cast<int *> (malloc (numsplittimes * sizeof (int)));
 for (int x = 0; x < numsplittimes; x++) {
        y[x] = T[x].upinf[IM_UPDATE_TIME_RY1].tries;
        z[x] = T[x].upinf[IM_UPDATE_TIME_RY1].accp;
 }
        try {
                MPI::COMM_WORLD.Reduce(y, y_rec, numsplittimes, MPI::INT, MPI::SUM, 0);
        } catch (MPI::Exception e) {
                std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
                MPI::COMM_WORLD.Abort(-1);
        }
        try {
                MPI::COMM_WORLD.Reduce(z, z_rec, numsplittimes, MPI::INT, MPI::SUM, 0);
        } catch (MPI::Exception e) {
                std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
                MPI::COMM_WORLD.Abort(-1);
        }
        if (currentid == 0) {
        for (int x = 0; x < numsplittimes; x++) {
                T[x].upinf[IM_UPDATE_TIME_RY1].tries = y_rec[x];
                T[x].upinf[IM_UPDATE_TIME_RY1].accp = z_rec[x];
	}
	}
  XFREE(y);
  XFREE(y_rec);
  XFREE(z);
  XFREE(z_rec);
}
#endif
#endif

#ifdef DO_NWUPDATE
#ifdef MPI_ENABLED
if (numprocesses > 1) {
int *y = static_cast<int *> (malloc (numsplittimes * sizeof (int)));
int *y_rec = static_cast<int *> (malloc (numsplittimes * sizeof (int)));
int *z = static_cast<int *> (malloc (numsplittimes * sizeof (int)));
int *z_rec = static_cast<int *> (malloc (numsplittimes * sizeof (int)));
 for (int x = 0; x < numsplittimes; x++) {
        y[x] = T[x].upinf[IM_UPDATE_TIME_NW].tries;
        z[x] = T[x].upinf[IM_UPDATE_TIME_NW].accp;
 }
        try {
                MPI::COMM_WORLD.Reduce(y, y_rec, numsplittimes, MPI::INT, MPI::SUM, 0);
        } catch (MPI::Exception e) {
                std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
                MPI::COMM_WORLD.Abort(-1);
        }
        try {
                MPI::COMM_WORLD.Reduce(z, z_rec, numsplittimes, MPI::INT, MPI::SUM, 0);
        } catch (MPI::Exception e) {
                std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
                MPI::COMM_WORLD.Abort(-1);
        }
        if (currentid == 0) {
        for (int x = 0; x < numsplittimes; x++) {
                T[x].upinf[IM_UPDATE_TIME_NW].tries = y_rec[x];
                T[x].upinf[IM_UPDATE_TIME_NW].accp = z_rec[x];
        }
        }

 XFREE(y);
 XFREE(y_rec);
 XFREE(z);
 XFREE(z_rec);
}
#endif
#endif
  for (i = 0; i < numsplittimes; i++)
    reclist[i] = (T + i);

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1
      || npops == 1)
  {
    assert (numsplittimes == 0);
    /* no split time */
  }
  else
  {
    if (currentid == 0) 
	    printacceptancerates (outto, numsplittimes, reclist,
                          "Update Rates -- Population Splitting Times");
  }
#ifdef DO_RY1UPDATE
#ifdef MPI_ENABLED
        if (currentid == 0 && numprocesses > 1) {
                for (int j = 0; j < numsplittimes; j++) {
                        T[j].upinf[IM_UPDATE_TIME_RY1].tries = reclist_saved[j]->upinf[IM_UPDATE_TIME_RY1].tries;
                        T[j].upinf[IM_UPDATE_TIME_RY1].accp = reclist_saved[j]->upinf[IM_UPDATE_TIME_RY1].accp;
                }
        }
#endif
#endif


#ifdef DO_NWUPDATE
#ifdef MPI_ENABLED
        if (currentid == 0 && numprocesses > 1) {
                for (int j  = 0; j < numsplittimes; j++) {
                        T[j].upinf[IM_UPDATE_TIME_NW].tries = reclist_saved[j]->upinf[IM_UPDATE_TIME_NW].tries;
                        T[j].upinf[IM_UPDATE_TIME_NW].accp = reclist_saved[j]->upinf[IM_UPDATE_TIME_NW].accp;
                }
        }
#endif
#endif


// genealogy updates
  for (li = 0; li < nloci; li++)
    reclist_saved[li] = L[li].g_rec;

#ifdef MPI_ENABLED

if (numprocesses > 1) {

int *y = static_cast<int *> (malloc (nloci * sizeof (int)));
int *y_rec = static_cast<int *> (malloc (nloci * sizeof (int)));
int *z = static_cast<int *> (malloc (nloci * sizeof (int)));
int *z_rec = static_cast<int *> (malloc (nloci * sizeof (int)));
 for (int x = 0; x < nloci; x++) {
        y[x] = L[x].g_rec->upinf[IM_UPDATE_GENEALOGY_ANY].tries;
        z[x] = L[x].g_rec->upinf[IM_UPDATE_GENEALOGY_ANY].accp;
 }
        try {
                MPI::COMM_WORLD.Reduce(y, y_rec, nloci, MPI::INT, MPI::SUM, 0);
        } catch (MPI::Exception e) {
                std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
                MPI::COMM_WORLD.Abort(-1);
        }
        try {
                MPI::COMM_WORLD.Reduce(z, z_rec, nloci, MPI::INT, MPI::SUM, 0);
        } catch (MPI::Exception e) {
                std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
                MPI::COMM_WORLD.Abort(-1);
        }
        if (currentid == 0) {
        for (int x = 0; x < nloci; x++) {
                L[x].g_rec->upinf[IM_UPDATE_GENEALOGY_ANY].tries = y_rec[x];
                L[x].g_rec->upinf[IM_UPDATE_GENEALOGY_ANY].accp = z_rec[x];
        }
        }

for (int x = 0; x < nloci; x++) {
        y[x] = L[x].g_rec->upinf[IM_UPDATE_GENEALOGY_TOPOLOGY].tries;
        z[x] = L[x].g_rec->upinf[IM_UPDATE_GENEALOGY_TOPOLOGY].accp;
 }
        try {
                MPI::COMM_WORLD.Reduce(y, y_rec, nloci, MPI::INT, MPI::SUM, 0);
        } catch (MPI::Exception e) {
                std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
                MPI::COMM_WORLD.Abort(-1);
        }
        try {
                MPI::COMM_WORLD.Reduce(z, z_rec, nloci, MPI::INT, MPI::SUM, 0);
        } catch (MPI::Exception e) {
                std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
                MPI::COMM_WORLD.Abort(-1);
        }
        if (currentid == 0) {
        for (int x = 0; x < nloci; x++) {
                L[x].g_rec->upinf[IM_UPDATE_GENEALOGY_TOPOLOGY].tries = y_rec[x];
                L[x].g_rec->upinf[IM_UPDATE_GENEALOGY_TOPOLOGY].accp = z_rec[x];
        }
        }
for (int x = 0; x < nloci; x++) {
        y[x] = L[x].g_rec->upinf[IM_UPDATE_GENEALOGY_TMRCA].tries;
        z[x] = L[x].g_rec->upinf[IM_UPDATE_GENEALOGY_TMRCA].accp;
 }
        try {
                MPI::COMM_WORLD.Reduce(y, y_rec, nloci, MPI::INT, MPI::SUM, 0);
        } catch (MPI::Exception e) {
                std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
                MPI::COMM_WORLD.Abort(-1);
        }
        try {
                MPI::COMM_WORLD.Reduce(z, z_rec, nloci, MPI::INT, MPI::SUM, 0);
        } catch (MPI::Exception e) {
                std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
                MPI::COMM_WORLD.Abort(-1);
        }
        if (currentid == 0) {
        for (int x = 0; x < nloci; x++) {
                L[x].g_rec->upinf[IM_UPDATE_GENEALOGY_TMRCA].tries = y_rec[x];
                L[x].g_rec->upinf[IM_UPDATE_GENEALOGY_TMRCA].accp = z_rec[x];
        }
        }

  XFREE(y);
  XFREE(y_rec);
  XFREE(z);
  XFREE(z_rec);
}

#endif


for (li = 0; li < nloci; li++)
    reclist[li] = L[li].g_rec;

  if (currentid == 0) 
	  printacceptancerates (outto, nloci, reclist, "Update Rates -- Genealogies");

#ifdef MPI_ENABLED
if (numprocesses > 1 && currentid == 0) {
        for (li = 0; li < nloci; li++) {
                L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_ANY].tries = reclist_saved[li]->upinf[IM_UPDATE_GENEALOGY_ANY].tries;
                L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_TOPOLOGY].tries = reclist_saved[li]->upinf[IM_UPDATE_GENEALOGY_TOPOLOGY].tries;
                L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_TMRCA].tries = reclist_saved[li]->upinf[IM_UPDATE_GENEALOGY_TMRCA].tries;
                L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_ANY].accp = reclist_saved[li]->upinf[IM_UPDATE_GENEALOGY_ANY].accp;
                L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_TOPOLOGY].accp = reclist_saved[li]->upinf[IM_UPDATE_GENEALOGY_TOPOLOGY].accp;
                L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_TMRCA].accp = reclist_saved[li]->upinf[IM_UPDATE_GENEALOGY_TMRCA].accp;
        }
}
#endif

// mutation rate scalars
  if (nurates > 1
      && (runoptions[PRINTMUTATIONUPDATESTOSCREEN] || outto != stdout))
  {
#ifdef MPI_ENABLED
   if (numprocesses > 1 && currentid == 0) {
    for (i = 0, li = 0; li < nloci; li++)
      for (j = 0; j < L[li].nlinked; j++)
      {
        reclist_saved[i] = &L[li].u_rec[j];
        i++;
      }
   }
#endif
#ifdef MPI_ENABLED
if (numprocesses > 1) {
 for (li = 0;  li < nloci; li++) {
        int *y = static_cast<int *> (malloc (L[li].nlinked * sizeof (int)));
        int *y_rec = static_cast<int *> (malloc (L[li].nlinked * sizeof (int)));
        int *z = static_cast<int *> (malloc (L[li].nlinked * sizeof (int)));
        int *z_rec = static_cast<int *> (malloc (L[li].nlinked * sizeof (int)));

        for (j = 0; j < L[li].nlinked; j++) {
                y[j] = L[li].u_rec[j].upinf->tries;
                z[j] = L[li].u_rec[j].upinf->accp;
        }
        try {
                MPI::COMM_WORLD.Reduce(y, y_rec, L[li].nlinked, MPI::INT, MPI::SUM, 0);
        } catch (MPI::Exception e) {
                std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
                MPI::COMM_WORLD.Abort(-1);
        }
        try {
                MPI::COMM_WORLD.Reduce(z, z_rec, L[li].nlinked, MPI::INT, MPI::SUM, 0);
        } catch (MPI::Exception e) {
                std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
                MPI::COMM_WORLD.Abort(-1);
        }
        if (currentid == 0) {
                for (j = 0; j < L[li].nlinked; j++) {
                        L[li].u_rec[j].upinf->tries = y_rec[j];
                        L[li].u_rec[j].upinf->accp = z_rec[j];
                }
        }
        XFREE(y);
        XFREE(y_rec);
        XFREE(z);
        XFREE(z_rec);

        }
}
#endif
for (i = 0, li = 0; li < nloci; li++)
      for (j = 0; j < L[li].nlinked; j++)
      {
        reclist[i] = &L[li].u_rec[j];
        i++;
      }
// kappa values for HKY model
    if (currentid == 0)
	    printacceptancerates (outto, i, reclist,
                          "Update Rates -- Mutation Rate Scalars");
#ifdef MPI_ENABLED
if (numprocesses > 1 && currentid == 0) {
        for (li = 0; li < nloci; li++) {
                for (j = 0; j < L[li].nlinked; j++) {
                        L[li].u_rec[j].upinf->tries = reclist_saved[li]->upinf[j].tries;
                        L[li].u_rec[j].upinf->accp = reclist_saved[li]->upinf[j].accp;
                }
        }
}
#endif



    for (i = 0, li = 0; li < nloci; li++)
      for (j = 0; j < L[li].nlinked; j++)
      {
        if (L[li].umodel[j] == HKY)
        {
          reclist_saved[i] = L[li].kappa_rec;
          i++;
        }
      }

#ifdef MPI_ENABLED
if (numprocesses > 1) {
 for (li = 0;  li < nloci; li++) {
        int *y = static_cast<int *> (malloc (L[li].nlinked * sizeof (int)));
        int *y_rec = static_cast<int *> (malloc (L[li].nlinked * sizeof (int)));
        int *z = static_cast<int *> (malloc (L[li].nlinked * sizeof (int)));
        int *z_rec = static_cast<int *> (malloc (L[li].nlinked * sizeof (int)));
        for (j = 0; j < L[li].nlinked; j++) {
                if (L[li].umodel[j] == HKY) {
                        y[j] = L[li].kappa_rec[j].upinf->tries;
                        z[j] = L[li].kappa_rec[j].upinf->accp;
                }
        }
        try {
                MPI::COMM_WORLD.Reduce(y, y_rec, L[li].nlinked, MPI::INT, MPI::SUM, 0);
        } catch (MPI::Exception e) {
                std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
                MPI::COMM_WORLD.Abort(-1);
        }
        try {
                MPI::COMM_WORLD.Reduce(z, z_rec, L[li].nlinked, MPI::INT, MPI::SUM, 0);
        } catch (MPI::Exception e) {
                std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
                MPI::COMM_WORLD.Abort(-1);
        }
        if (currentid == 0) {
                for (j = 0; j < L[li].nlinked; j++) {
                        if (L[li].umodel[j] == HKY) {
                                L[li].kappa_rec[j].upinf->tries = y_rec[j];
                                L[li].kappa_rec[j].upinf->accp = z_rec[j];
                        }
                }
        }
        XFREE(y);
        XFREE(y_rec);
        XFREE(z);
        XFREE(z_rec);

        }
}
#endif

for (i = 0, li = 0; li < nloci; li++)
      for (j = 0; j < L[li].nlinked; j++)
      {
        if (L[li].umodel[j] == HKY)
        {
          reclist[i] = L[li].kappa_rec;
          i++;
        }
      }

    if (i > 0 && currentid == 0)
      printacceptancerates (outto, i, reclist,
                            "Update Rates -- HKY Model Kappa parameter");

#ifdef MPI_ENABLED
if (numprocesses > 1 && currentid == 0) {
        for (li = 0; li < nloci; li++) {
                for (j = 0; j < L[li].nlinked; j++) {
                        if (L[li].umodel[j] == HKY) {
                                L[li].kappa_rec[j].upinf->tries = reclist_saved[li]->upinf[j].tries;
                                L[li].kappa_rec[j].upinf->accp = reclist_saved[li]->upinf[j].accp;
                        }
                }
        }
}
#endif




// STR ancestral allele states 
// A_rec not used as of sometime in 2010, A gets enough updates when updating genealogy
//8/26/2011  turn this printing section off, as it only ever prints zeros when A updating is not used 
    /*
    for (i = 0, li = 0; li < nloci; li++)
      for (j = 0; j < L[li].nlinked; j++)
      {
        if (L[li].umodel[j] == STEPWISE)
        {
          reclist[i] = &L[li].A_rec[j];
          i++;
        }
      }
      
    if (i > 0)
    {
      printacceptancerates (outto, i, reclist,
                            "Update Rates -- STR Genealogy Allele States");
    } */
    }
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    for (li = 0; li < nloci; li++)
      reclist[li] = L[li].a_rec;

    if (assignmentoptions[POPULATIONASSIGNMENTLOCAL] == 1)
    {
      printacceptancerates (outto, nloci, reclist,
                            "Update Rates -- Assignment");
    }
    // assignment updating information 
    // set reclest
    // call printacceptancerates 
    //
    printacceptancerates_multichain (outto);
  }

  return;
}                               //callprintacceptancerates


/* set up arrays pointing to information to put in curve plots,  then call asciicurve
   some things that are plotted are based on struct value_record and others on struct i_param
   this is why we cannot simply call asciicurve() with a pointer to a single type of structure  */
void callasciicurves (void)
{
  struct plotpoint **curvexy;
  char **curvestr;
  int *curve_do_logplot;
  int numcurve = 0;
  int i, j, li, ui;
  int *nrecstep;
  //_CrtCheckMemory( );
// find out how many curves
  numcurve += numpopsizeparams;
  for (i = 0; i < nummigrateparams; i++)
    if (imig[i].pr.max > MPRIORMIN)
      numcurve++;
  numcurve += numsplittimes;
  if (runoptions[LOADRUN] == 0 && nurates > 1)
    numcurve += nurates;
  if (runoptions[LOADRUN] == 0)
    for (li = 0; li < nloci; li++)
      if (L[li].model == HKY)
        numcurve++;
// allocate
  curvexy = static_cast<plotpoint **> 
                (malloc (numcurve * sizeof (struct plotpoint *)));
  curvestr = static_cast<char **> (malloc (numcurve * sizeof (char *)));
  nrecstep = static_cast<int *> (malloc (numcurve * sizeof (int)));
  curve_do_logplot = static_cast<int *> (malloc (numcurve * sizeof (int)));
// assign
  j = 0;
  for (i = 0; i < numpopsizeparams; i++)
  {
    curvexy[j] = itheta[i].xy;
    curvestr[j] = &itheta[i].str[0];
    curve_do_logplot[j] = 0;
    nrecstep[j] = 1;
    j++;
  }
  for (i = 0; i < nummigrateparams; i++)
    if (imig[i].pr.max > MPRIORMIN)
    {
      curvexy[j] = imig[i].xy;
      curvestr[j] = &imig[i].str[0];
      curve_do_logplot[j] = 0;
      nrecstep[j] = 1;
      j++;
    }
  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1
      || npops == 1)
  {
    /* no split time */
  }
  else
  {
    for (i = 0; i < lastperiodnumber; i++)
    {
      curvexy[j] = T[i].v->xy;
      curvestr[j] = &T[i].v->str[0];
      curve_do_logplot[j] = T[i].v->do_logplot;
      nrecstep[j] = recordstep;
      j++;
    }
  }
  if (runoptions[LOADRUN] == 0 && nurates > 1)
  {
    for (li = 0; li < nloci; li++)
      for (ui = 0; ui < L[li].nlinked; ui++)
      {
        curvexy[j] = L[li].u_rec[ui].v->xy;
        curvestr[j] = &L[li].u_rec[ui].v->str[0];
        curve_do_logplot[j] = L[li].u_rec[ui].v->do_logplot;
        nrecstep[j] = recordstep;
        j++;
      }
  }
  if (runoptions[LOADRUN] == 0)
    for (li = 0; li < nloci; li++)
      if (L[li].model == HKY)
      {
        curvexy[j] = L[li].kappa_rec->v->xy;
        curvestr[j] = &L[li].kappa_rec->v->str[0];
        curve_do_logplot[j] = L[li].kappa_rec->v->do_logplot;
        nrecstep[j] = recordstep;
        j++;
      }
  assert (numcurve == j);
  for (j = 0; j < numcurve; j++)
    asciicurve (outfile, curvexy[j], curvestr[j], curve_do_logplot[j],
                nrecstep[j]);
//free
  XFREE (curvexy);
  XFREE (curvestr);
  XFREE (curve_do_logplot);
  XFREE (nrecstep);
}                               //callasciicurve 

// makes calls to asciitrend
void callasciitrend (FILE * outtofile)
{
  int i, li, ui;
  asciitrend (outtofile, lpgpd_v, trenddoublepoint, trendspot);

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1
      || npops == 1)
  {
    /* no split time */
  }
  else
  {
    for (i = 0; i < lastperiodnumber; i++)
      asciitrend (outtofile, T[i].v, trenddoublepoint, trendspot);
  }
  if (nurates > 1 && runoptions[LOADRUN] == 0)
  {
    for (li = 0; li < nloci; li++)
      for (ui = 0; ui < L[li].nlinked; ui++)
        asciitrend (outtofile, L[li].u_rec[ui].v, trenddoublepoint, trendspot);
  }
  for (li = 0; li < nloci; li++)
    if (L[li].model == HKY && runoptions[LOADRUN] == 0)
      asciitrend (outtofile, L[li].kappa_rec->v, trenddoublepoint, trendspot);

  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    for (li = 0; li < nloci; li++)
    {
      asciitrend (outtofile, L[li].a_rec->v, trenddoublepoint, trendspot);
    }
  }
  return;
}                               // callasciitrend 

void 
printoutput (int currentid)         // mostly calls functions in output.c
{
  int i;
  double seconds;
  int p;
  float *holdpeakloc;
  double multitpeak[MAXPOPS - 1];

if (currentid == 0) {
  if (runoptions[LOADRUN] == 0)
  {
    /* savegenealogyfile moved out of printoutput because it does not seem to belong to it. */
    if (cdurationmode == TIMEINF)
    {
      remove (oldoutfilename);
      rename (outfilename, oldoutfilename);
    }
  }

  if ((outfile = fopen (outfilename, "w")) == NULL)
  {
    IM_err (IMERR_CREATEFILEFAIL, "Error opening text file for writing");
  }
 
  printrunbasics (outfile, runoptions[LOADRUN], fpstr, burnsteps, recordint,
                  recordstep, savegenealogyint, endtime, starttime, hilike,
                  hiprob/*, step*/);
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    imaAsnPrintNumGenesPopn (outfile);
  }
 }
  // Bayes factor stuff of Sang Chul's fprintf (outfile, "Average loglikelihood: %lf\n", gloglikelihood);
  if (runoptions[LOADRUN] == 0)
  { //if (currentid == 0) {
    callprintacceptancerates (outfile, currentid);
    //}
    if (numprocesses * numchains > 1 && runoptions[LOADRUN] != 1) { 
	      printchaininfo (outfile, heatmode, hval1, hval2, currentid);
	}
	if (currentid == 0) {
   callprintautoctable (outfile/*, step*/);
    }
  }
if (currentid == 0) {
  if (outputoptions[PARAMGREATERTHAN])
  {
    print_greater_than_tests (outfile);
  }
  if (!modeloptions[EXPOMIGRATIONPRIOR] && calcoptions[DONTCALCLIKELIHOODMUTATION]==0)        //as of 11/19/09 have not yet done the math for case of migration with exponential prior
    print_means_variances_correlations (outfile);
  // init_surface_calc (); not used as of 8/24/09
/*  get marginal peaks */
  if (!calcoptions[DONTCALCLIKELIHOODMUTATION])
  {
    p = numpopsizeparams + nummigrateparams ;
    holdpeakloc = static_cast<float *> (malloc (p * sizeof (float)));
    printf ("surface calculations . . .\n");
    if (modeloptions[EXPOMIGRATIONPRIOR] || (runoptions[LOADRUN] && calcoptions[FINDJOINTPOSTERIOR]))
      eexpsum = (struct extendnum *) malloc ((size_t) ((genealogiessaved + 1) * sizeof (struct extendnum)));

    findmarginpeaks (outfile, holdpeakloc);

    if (runoptions[LOADRUN] && calcoptions[FINDJOINTPOSTERIOR])
    {
      closeopenout (&outfile, outfilename);
      /* CR 110921.1  Change outfile parameter type to (FILE **) */
      findjointpeaks(&outfile,outfilename,nestedmodelfilename,p);
    }

  }
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    imaAsnPrintNumGenesPopn (outfile);
  }
  else
  {
/* get joint splittime peak */
    if (npops >=3  && npops <= 5 && !runoptions[LOADRUN] && outputoptions[PRINTJOINTTEST])
    {
      return_joint_t (multitpeak);
      FP "\nEstimated joint splitting time from multi-dimensional histogram\n");
      FP "  number of bins per dimension %d\n", NUMTARRAYBINS);
      FP "  Posterior probability of estimated joint value of splitting time: %7.4lf\n", joint_t_prob (&multitpeak[0]));
      FP "---------------------------------------------------------------\n");
      for (i = 0; i < numsplittimes; i++)
        FP "   %s\t%.3lf\n", T[i].str, multitpeak[i]);
      FP "\n\n");
    }
  }
  printhistograms (outfile, recordstep, generationtime, scaleumeaninput);
  if (!outputoptions[DONTPRINTASCIICURVE])
  {
    FP "\n\nASCII Curves - Approximate Posterior Densities \n");
    FP "===================================================\n");
    callasciicurves ();
  }
  if (!outputoptions[DONTPRINTASCIITREND])
  {
    if (trendspot <= 1)
    {
      FP "burn period too short to plot trends,  trend recording begins at step %d \n",burntrendstartdelay); 
    }
    else
    {
      FP "\n\nASCII Plots of Parameter Trends \n");
      FP "===================================================\n");
      FP " - note points to the left of '!' on X axis have twice the density in time relative to points to the right\n\n");
      callasciitrend (outfile);
    }

  }
  time (&endtime);
  seconds = difftime (endtime, starttime);
  FP "Time Elapsed : %d hours, %d minutes, %d seconds \n\n",
    (int) seconds / (int) 3600,
    ((int) seconds / (int) 60) - ((int) 60 * ((int) seconds / (int) 3600)),
    (int) seconds - (int) 60 *((int) seconds / (int) 60));
  FP "\nEND OF OUTPUT\n");
  f_close (outfile);
  //free_surface_calc ();

  if (outputoptions[MIGRATEHIST] && nummigrateparams > 0)
  {

    if ((migplotfile = fopen (migplotfilename, "w")) == NULL)
    {
      IM_err (IMERR_CREATEFILEFAIL,
              "Error opening file for plotting migration amounts and times");
    }
    printmigrationhistograms (migplotfile, recordstep);
    f_close (migplotfile);
  }


  if (!calcoptions[DONTCALCLIKELIHOODMUTATION])
    XFREE (holdpeakloc);
  if (runoptions[SAVEMCSTATEFILE])
  {
    writemcf (mcfwritefilename);
  }
  if (jheyoptions[WRITEMIGRATIONNAME])
  {
    f_close(migrationnamefile);
    migrationnamefile = fopen (migrationnamefilename, "a");
  }
}
  return;
}                               /* printoutput */


void intervaloutput (FILE * outto, int currentid)
{
  int li;
  int j;
  double like;
  int ci;
  double multitpeak[MAXPOPS - 1];
  ci = 0;
  if (currentid == 0) {
  checkhighs (0, printint, &hilike, &hiprob, &like/*, step*/);
  }
  if (((step / (int) printint) * (int) printint == step && step > 0)
      || outto != stdout)
  {
   if (currentid == 0)
    printsteps (outto, like);

    callprintacceptancerates (outto, currentid);
   if (currentid == 0) {
    printcurrentvals (outto);
    callprintautoctable (outto /*, step*/);
    }
    if (numprocesses * numchains > 1)
      printchaininfo (outto, heatmode, hval1, hval2, currentid);
    if (currentid == 0) {
    if (genealogiessaved > 0)
    {
      time (&remained_endtime);
    }
    }
    /* For ASSIGNMENT */
    if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
    {
      printsteps (outto, like);
      imaAsnPrintNumGenesPopn (outto);
      for (ci = 0; ci < numchains; ci++)
      {
        for (j = 0; j < Cupinf[ci].num_uptypes; j++)
        {
          Cupinf[ci].upinf[j].accp = 0;
          Cupinf[ci].upinf[j].tries = 0;
        }
      }
      for (li = 0; li < nloci; li++)
      {
        for (j = 0; j < L[li].a_rec->num_uptypes; j++)
          L[li].a_rec->upinf[j].accp = L[li].a_rec->upinf[j].tries = 0;
/*
        if (assignmentoptions[POPULATIONASSIGNMENTASSIGNED] == 1)
        {
          IMA_assigned_reset ();
        }
*/
      }
    }

  }
if (currentid == 0) {
if (burndone)
  return_joint_t (multitpeak);
}
  return;
}                               /* intervaloutput */

// check if it is time to call record() and savegenealogyinfo(), and call if it is
void check_to_record (int current_id)
{
  static int i;
  static int j;
  static int init = 0;

  if (init == 0)
  {
    i = recordint;
    j = savegenealogyint;
    init = 1;
  }
  if (i == recordint)
  {
    record (current_id);                  // record values of parameters that are in mcmc
    recordstep++;
    i = 1;
  }
  else
  {
    i++;
  }

  if (j == savegenealogyint)
  {
    savegenealogyinfo (current_id);            // record values associated with genealogies
    if (jheyoptions[WRITEMIGRATIONNAME])
      record_migration_names();
    if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
    {
      recordassignmentloci (outfilename, 0);
    }
    genealogiessaved++;
 /*   if (genealogiessaved == genealogiestosave && current_id == 0) {
		std::cout << "genealogiessaved here in chain 0 is " << genealogiessaved << " and broadcasting now\n";
		try {
			MPI::COMM_WORLD.Bcast(&genealogiessaved, 1, MPI::INT, current_id);
		} catch (MPI::Exception e) {
			std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
			MPI::COMM_WORLD.Abort(-1);
			return;
		}
			
    }*/
    j = 1;
  }
  else
  {
    j++;
  }
}                               // check_to_record


void fillswaparrays(int *swapper, int *swappee)
{
	for (int i = 0; i < /*chainduration + burnduration +*/ 1; i++) {
		swapper[i] = rand() % numprocesses;
		swappee[i] = rand() % numprocesses;
		//MPI::COMM_WORLD.Bcast(&swapper[i], 1, MPI::INT, 0);
		//MPI::COMM_WORLD.Bcast(&swappee[i], 1, MPI::INT, 0);
	}
	return;
}

int main (int argc, char *argv[])
{
  #ifdef MPI_ENABLED
  MPI::Init(argc, argv);
  MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);
  MPI::Status status;
  #endif
  int currentid = 0;
  //AS: Debug only
  //std::ofstream f1;
  //char filename[15];
  /*extendor = 0;*/
  #ifdef MPI_ENABLED
  try {
	numprocesses = MPI::COMM_WORLD.Get_size();
	} catch (MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return 0;
	}
  try {
	currentid = MPI::COMM_WORLD.Get_rank();
	} catch (MPI::Exception e) {
		std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
		MPI::COMM_WORLD.Abort(-1);
		return 0;
	}
  #else
  numprocesses = 1;
  #endif
//AS: Debug only
//	if (currentid == 0) { 
//	sprintf(filename, "mctrace-%03d", currentid);
//	try {
//		f1.open(filename);
//	} catch (std::ofstream::failure e) {
//		std::cout << "Unable to open file " << filename << " for writing\n";
//		return (EXIT_FAILURE);
//	}
//	}
/*	sprintf(filename, "autoc-%03d", currentid);
	
  	if ((autocfile = fopen (filename, "w")) == NULL)
	  {
	    IM_err (IMERR_CREATEFILEFAIL, "Error opening autoc text file for writing");
	  }
*/

//if (numchains > 1 || currentid == 0) {
  init_IMA ();
  start (argc, argv, currentid);
	///AS: This has to change in case of chainduration and burnduration specified in minutes/hours
	//AS: Debug only
	//std::cout << "Chain duration " << chainduration << "\t Burn duration " << burnduration << "\n";
	//std::cout << "Genealogies to save " << genealogiestosave << "\n";	
	//AS: Old code - decided to pick swapper and swappee on the fly 7/2/2014
	
    	//swapper = static_cast<int *> (malloc ((/*chainduration + burnduration +*/ 1) * sizeof (int)));
	//swappee = static_cast<int *> (malloc ((/*chainduration + burnduration +*/ 1) * sizeof (int)));
	//fillswaparrays(swapper, swappee);

  if (runoptions[LOADRUN])
  {
	if (currentid == 0) {
	    loadgenealogyvalues ();
	    recordstep = genealogiessaved;
    if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
    {
      IMAgenealogyAlign (outfilename, genealogiessaved);
    }
	}
  printoutput (currentid);
	if (currentid == 0)
	  free_ima_main_stuff ();

  } else if (runoptions[LOADRUN] != 1){
  if (currentid == 0) {
  	printf ("Starting Markov chain.\n");
  }

  init_after_start_IMA ();
  step = 0;
  recordstep = 0;
//std::cout << "going into run...numchains is "<< numchains << "\n";  
//if (numchains > 1 || currentid == 0) {
while (run (currentid))
{
	//if (currentid == 0) {
	//AS: changing this as on 7/1/2014 to directly pick the chains that are swapping
	//AS: this is to avoid any bias in picking which chains are involved in swap
	//AS: we previously noted that there were more swaps within a process, than between
	//std::cout << "Here?\n";
	int swapA, swapB;
		if (numprocesses * numchains > 1) {
		swapA = rand() % (numprocesses * numchains);
		swapB = rand() % (numprocesses * numchains);
	
		while (swapA == swapB) {
			swapB = rand() % (numprocesses * numchains);
		}
		}
	//std::cout << "Step is: " << step << " and I am processor " << currentid << "\n";
      qupdate (currentid,/*swapper, swappee, f1,*/ swapA, swapB);
	//f1 << "Step is (qupdate done): " << step << " and I am Chain " << currentid << "\n";
    if (burndone)
    {
      check_to_record (currentid);
#ifdef COUNT_PRIOR_ASSIGNMENT
      if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
      {
        IMA_rgf_tickSaasn (0);
      }
#endif /* COUNT_PRIOR_ASSIGNMENT */
	int z = whichiscoldchain();
	if ( z >= 0 ) {
	      gloglikelihood += C[z]->allpcalc.pdg; 
	      gloglikelihood /= 2.0;
		/*if (numprocesses > 1) {
			std::cout << "Global log likelihood is " << gloglikelihood << " and broadcasting...\n";
		}*/
    	}
//AS: Debug only
//	f1 << "Chain ID is " << currentid << " Genealogies saved is" << genealogiessaved << "\n"; 

  }
	step++;
	//AS: Debug only
	//if (currentid == 0) {
	//	std::cout << "Current step is " << step << "\n";
	//}
}

	#ifdef MPI_ENABLED
	MPI::COMM_WORLD.Barrier();
	#endif
	rundone = 1;
//AS: Debug only
//	std::cout << "Chain ID is " << currentid << " Max genealogies to be saved is" << MAXGENEALOGIESTOSAVE << "\n";
//	std::cout << "rundone is " << rundone << "\n"; 
#ifdef COUNT_PRIOR_ASSIGNMENT
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    IMA_rgf_saveSaasn ();
  }
#endif /* COUNT_PRIOR_ASSIGNMENT */

  /* save genealogy info in *.ti file */
  if (!runoptions[DONTSAVEGENEALOGIES] && genealogiessaved > 0 && currentid == 0)
  {
//	std::cout << "Waiting here inside savegenealogyfile " << genealogiessaved << "\n"; 

    savegenealogyfile (genealogyinfosavefilename, genealogyinfosavefile, &lastgenealogysaved, gsampinflength);
  }
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1
      && assignmentoptions[POPULATIONASSIGNMENTASSIGNED] == 0)
  {
    IMAgenealogyAlign (outfilename, genealogiessaved);
  }

	//AS: combine all swap variables before calling printoutput on head node
	//AS: note that only the burntrend file printed on the head node, and the final output file (which is also printed on the head node)
	//will be able to print out the swap counts between successive chains.
	//AS: MPI check	
	/*#ifdef MPI_ENABLED
	if (numprocesses > 1) {
		for (int x = 0; x < numprocesses; x++) {
			for (int y = 0; y < numprocesses; y++) {
				try {
					MPI::COMM_WORLD.Reduce(&swaps_bwprocesses[x][y], &swaps_rec_bwprocesses[x][y], 1,
					MPI::INT, MPI::SUM, 0);
				} catch (MPI::Exception e) {
					std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
					MPI::COMM_WORLD.Abort(-1);
				}
			}
		}
		
		for (int x = 0; x < numprocesses * numchains; x++) {
			for (int y = 0; y < numprocesses * numchains; y++) {
				try {
					MPI::COMM_WORLD.Reduce(&tempbasedswapcount[x][y], &tempbased_rec_swapcount[x][y], 1,
								MPI::INT, MPI::SUM, 0);
				} catch (MPI::Exception e) {
					std::cout << e.Get_error_string() << e.Get_error_code() << "\n";
					MPI::COMM_WORLD.Abort(-1);
				}
			}
		}
		
		//AS: also have to send-receive all the swapcount variables into the swapcount_bwprocesses matrix
		if (currentid == 0) {
			for (int x = 0; x < numchains; x++) {
				for (int y = 0; y < numchains; y++) {
					swapcount_bwprocesses[x][y] = swapcount[x][y];
				}
			}
		}	
		for (int x = 1; x < numprocesses; x++) {
			if (currentid == x) {
				for (int y = 0; y < numchains; y++) {
					for (int z = 0; z < numchains; z++) {
						MPI::COMM_WORLD.Send(&swapcount[y][z], 1, MPI::INT, 0, 1234);
					}
				}
			}
			if (currentid == 0) {
				for (int y = 0; y < numchains; y++) {
					for (int z = 0; z < numchains; z++) {
						MPI::COMM_WORLD.Recv(&swapcount_bwprocesses[x * numchains + y][x * numchains + z], 1, MPI::INT, x, 1234);
					}
				}
			}
		}
			
			

		}
	#endif*/
		
		
		///AS: At the end of this, my new swapcount_bwprocesses matrix should have all counts added up
	
  //if (currentid == 0) {
	//AS: Debug only
	//std::cout << "Waiting for printing to finish up\n";
  printoutput (currentid);
  //}
// callprintautoctable(autocfile);

  if (jheyoptions[WRITEMIGRATIONNAME])
    f_close(migrationnamefile);
  free_ima_main_stuff ();
}
	//AS: Debug only
//   f_close(autocfile);
//	if (currentid == 0) {
//  f1.close();
//  }

  #ifdef MPI_ENABLED
  MPI::Finalize();
  #endif
  return 0;
}                               /* main */
