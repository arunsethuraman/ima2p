/*IMa2p 2009-2015 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, and Arun Sethuraman */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <time.h>
#include <ctype.h>
#include <assert.h>
#include <stdarg.h>
#include <iostream> 
#include <string> 

/* AS: uncomment this if compiling from source */
/* #ifndef MPI_ENABLED
 * #define MPI_ENABLED
 * #endif
 */


//AS: need config.h file for MPI definition
//AS: comment this out if compiling from source!!
#include "config.h"
#ifdef MPI_ENABLED
#include <mpi.h>
#endif
#include <fstream>
/* SANGCHUL: Wed Oct 21 15:08:02 EDT 2009
 * Use -DNDEBUG for a release version.
 * Use -DDEBUG for a debug version.
 * Note that _DEBUG macro is not standard. Visual C++ does define it when it
 * compiles source codes in debug mode. Compiling IMa2 with DEBUG macro is
 * preferred in linux and mac.
 */

#define RELEASE_DATE "January, 2015" 
//#define MPI_ENABLED

//#define SANITY_TEST   // ONLY used for running testbed 

#if defined(_DEBUG) && !defined(DEBUG)
#define DEBUG
#endif /* _DEBUG or DEBUG */

#ifdef  _MSC_VER

//#define NTDDI_VERSION  NTDDI_WIN7
#define NTDDI_VERSION 0x05000000   //NTDDI_WIN2K  
#define WINVER 0x0500   // _WIN32_WINNT_WIN2K
#define _WIN32_WINNT 0x0500  //_WIN32_WINNT_WIN2K

#ifndef WIN64 
// JH added 12/22/08 ,  this is a leak detection library,  only invoked in MVC++ debug mode and if vld is installed
//#include <vld.h> // a useful library for checking memory leaks under microsoft visual c++ 8.0 32bit, otherwise do not use this
#endif //WIN64
#endif /*  _MSC_VER */

#include "xtrapbits.hpp"
#include "utilities.hpp"

/* Microsoft compilers use _findfirst and _findnext to search directories for files.  
Most other compilers use a library called dirent.  It is not a formal standard
but is very widespread.  dirent does not have _findfirst or _findnext but instead 
uses functions called  opendir() closedir() readdir() and rewinddir().  
In order for this code to be portable,  I've included a file called sdirent that 
maps the usual dirent functions calls onto the microsoft functions.  If the compiler 
is not microsoft, then this file is not included, and the regular direct is used */

#ifdef _MSC_VER
#include <io.h>
//#include "msdirent.hpp"
#include <errno.h>              /* _findfirst and _findnext set errno iff they return -1 so these must be included */
//#elif  __GNUC__
#else /*  */
#include <dirent.h>
#endif /*  */

#ifndef _MSC_VER
#ifndef __forceinline
#define __forceinline __attribute__((__always_inline__)) inline
#endif /*  */
#endif /*  */

/* microsft visual studio compiler stuff */
/* disable deprecate warning  C4996 */
#ifdef _MSC_VER
#pragma warning( disable : 4996)
#endif /* _MSC_VER */

#define _CRT_SECURE_NO_WARNINGS

#ifdef _DEBUG
#ifdef _MSC_VER
#define MYDEBUG
# endif
#endif /*  */


/* By turning on conditional compiler definition MORESTABLE, 
 * we use LogDiff or LogSum2. Otherwise, we use the original form of equations.
 * Search codes for LogDiff and LogSum2 for example. */
#define MORESTABLE

/* splittime updateing  */
#define DO_NWUPDATE
//#undef DO_NWUPDATE    // use undef to turn this update off
#define DO_RY1UPDATE
//#undef DO_RY1UPDATE  // use undef to turn this update off

/***********************************************/
/******** SIMPLE DEFINITIONS MACROS ************/
/***********************************************/


/* CONSTANTS */
#define MAXLOCI  1000            //1000
#define MAXCHAINS 1000          //201
#define MAXGENES 1000           // maximum sample size for a locus  
#define MAXLENGENENAME 12       /* a gene name can be up to 10, plus 11-th can be @ for diploid */
#define LENGENENAME 10          /* a gene name can be up to 10 */
#define MAXPOPS  10             // MAXPOPS cannot exceed 10 because the treestring functions assume populations and nodes are represented by single integers
#define MAXPERIODS MAXPOPS+1
#define MAXTREEPOPS  (2*MAXPOPS - 1)
#define FNSIZE    500           // max file name length
#define POPTREESTRINGLENGTHMAX  100
#define NAMELENGTH 151          // max length of population names and locus names
#define MAXLINKED 15            // largest number of linked loci with the same genealogy - each neads its own mutation rate
#define TIMEMAX 1000000.0       // no branch can have a bottom time greater than this
#define MINPARAMVAL 0.0000001   //0.0001      // smallest parameter value for parameters that are in the MCMC
#define STARTMIGMAX 1000        /* max number per edge when trees are first built */
#define ADDMIGMAX 1000          /* max of number of migrations that can be added by poisson generator */
#define ABSMIGMAX  5000         /* absolute maximum # of migration events allowed  -  is 5000 too large? */
#define PARAMSTRLEN  12         /* length of string of parameter name */
#define PARAMSTRLENSHORT  7     /* length of string of parameter name */
#define UPDATELABELLEN 18       /* length of string for an update name */
#define MPRIORMIN 0.000001      /* small value for setting upper bound on migration to (near) zero */
#define BFMPRIORMIN 1.01        /* small value for setting upper bound on migration to (near) zero */
#define KAPPAMAX  100           /* maximum value of HKY parameter */
#define MIGINC 5                // 20 /* number of possbile migration events to add to branch migration array at a time *//* changed this from 20 to 5 on 4/3/08 */
#define DEFCHAINS  1            // default # chains
#define RECORDINTDEFAULT 10     // default # of steps between recording values of things - used to call record()
#define SAVEGENEALOGYDEFAULT 100  // default # of steps between recording information about the genealogies
#define MINSTRLENGTH  3         // minimum allowed number of STR repeats, so users don't use data that doesn't fit // JH changed to 3  11/29/2010
#define GRIDSIZE 1000           // # of bins in histogram
#define TRENDDIM 500            // number of points saved for the trendline plots
#define MARGIN2DGRIDSIZE 50     // # of bins on each axis in 2D histograms
#define PROFILEGRIDSIZE 50      // 100 //30  // # points along single dimension profile curve
#define MAXGENEALOGIESTOSAVE 300000 // the max to save in ram during a run - some bug,  setting this to 500,000 caseus memory management errors
#define PRINTINTDEFAULT 10000   // default # steps between writing to the screen
#define MAXLOADFILES 500        // max # of files with genealogy information that can be loaded
#define UMAX  10000.0           // highest value for u scalars  - can differ by UMAX^2 fold
#define HMAX  20                // highest value for h scalars, if HMAX=20 it gives a range of ratios for h scalars from 1/20  to 20  (i.e. they can differ by up to 400 fold)
#define DEFAULTBURNTRENDSTEP 10000      // default # steps between writing burn trend file,  if not given on command line
#define TIMEPRIORMULTIPLIER 10  // useful for plotting TMRCAs
#define EXPOMIGPLOTSCALE    20 //with exponential prior on migration, this number times the given prior mean value sets peak search interval and plot scale
#define REJECTINFINITESITESCONSTANT  -1000000000.0      // just some value not likelily to turn up in a calculation, indicates failure of IS model
#define NUMTARRAYBINS  100      // number of bins for multidimensional peak estimation of splitting times

/* autocorrelation estimation constants */
#define AUTOCTERMS 12           // the number of lag values for which autocorrelations are recorded
#define CHECKAUTOCWAIT 10000   /* 100000 */     // number of steps before recording autocorrelation values
#define AUTOCINT 1000           // interval between measurements for autocorrelations - cannot be changed easily
#define AUTOCNEXTARRAYLENGTH 1000       /*  = (largest value in autoc_checkstep[] divided by AUTOCINT) */
#define AUTOCCUTOFF 100  /*500 */       /* minimum number of measurements to have for autocorrelation before printing results */
#define AUTOCSTEPSCALAR 1//50//5       /* scale over which autocorrelation is measured, 1 means the scale is steps */



/* MACROS  */
#define  MYDBL_MAX DBL_MAX/1e10 // these avoid some over- under-flow issues, I think.  Not used much. May not be necessary
#define  MYDBL_MIN DBL_MIN * 1e10
#define IM_BESSI_MIN (-1e+100) /* very small value that would rejct any update with stepwise muation model */

#ifdef _MSC_VER
# define rint(x) floor((x) + 0.5)
#endif /*  _MSC_VER*/

#define POSROUND(a) (long) ((a)+0.5)    // simple rounding to positive integers
#define INTEGERROUND(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))        // simple rounding to integer
#define ODD(a) ( (((a) & 1) == 1) ? 1 : 0 )     /* for nonnegative integers , if odd returns 1, else 0 */
#define FP fprintf(outfile,     // handy way to avoid retyping the same thing again and again
#define SP *fpstri += sprintf(&fpstr[*fpstri],  // used to add text to the long string (called 'fpstr') that goes at the beginning of the output files
#define f_close(a)  fclose(a); (a) = NULL       //regular close() does not make the pointer null,  but it is useful to have this set to null
#define XFREE(p) do { free((p)); (p) = NULL; } while(0)

#define LOG10  2.3025850929940456840
#define LOG2 0.69314718055994530941723212146
#define LOG2HALF 0.34657359027997265470861606073   // half of LOG2, used in update_t_NW
#define LOG_DBL_MAX 7.0978271289338397e+02

/* these are supposedly bulletproof macros suggested by Melissa Hibnuz
 * but still, do not nest calls to these min and max macros */
static double sqrarg;

#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
static double dmaxarg1, dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))
static double dminarg1, dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ? (dminarg1) : (dminarg2))
static float maxarg1, maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
static float minarg1, minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ? (minarg1) : (minarg2))
static long lmaxarg1, lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ? (lmaxarg1) : (lmaxarg2))
static long lminarg1, lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ? (lminarg1) : (lminarg2))
static int imaxarg1, imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ? (imaxarg1) : (imaxarg2))
static int iminarg1, iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))

/**************/
/* SET MACROS */
/**************/
/* below are some routines for using sets ala pascal.
	I got this from "C reference manual" by Harbison p 175. useful for sets
	of small numbers from 0(inclusive) to the size of an unsigned long
	integer (32, exclusive). A number is in the set if the bit in the place
	position of that number (plus 1, because we permit the number 0 to be in
	the set) is a 1. If the bit is zero then it is not in the set.  Thus
	the numbers 0-31 can be in the set. */
/* any program using this must also declare typedef unsigned long set; */
/* changed to unsigned short to save memory - won't need more than 16 populations */
typedef unsigned int SET;

#define SET_BITS 32             // make sure this is correct for the compiler
#define EMPTYSET   ((SET) 0)    /* produces a SET of value 0 */
#define ISEMPTY(setj)  ((setj) == 0)    /*true if setj has no elements */
#define SINGLESET(i)  (((SET) 1) << (i))        /* makes a set of the integer i */
#define SETADD(setj,i)  ((setj) | SINGLESET (i))        /* setj with i added to it */
#define ISELEMENT(i,setj)   (SINGLESET((i))  & (setj) ) /*true if i in the setj */
#define INTERSECT(set1,set2)   ((set1) & (set2))
#define UNION(set1,set2)   ((set1) | (set2))
#define SETDIFF(set1,set2)   ((set1) ^ (set2))
#define SETREMOVE(setj,i)  ((setj) ^ SINGLESET(i))      /* setj with i removed from it , be sure that is is aleady in it */


#define FORALL(i,setj) \
  for ((i) = 0; (i) < SET_BITS; ++(i)) \
  if (ISELEMENT ((i), (setj)))

      /* forall : permits an action to be applied to all members of a set
         for example
         int k;
         forall(k,z) printf("%d ",k);   */

#define GETLOW(setj)  ( (setj) & -(setj) )
      /*returns the lowest member of the set, as a set 
         can be used in conjuntion with setdiff to reduce sets an item at a time */

/* Macro LogDiff approximately computes the value of v of the following
 * equation:
 * exp(v) = exp(a) - exp(b) 
 * where a > b. */
#define LogDiff(v,a,b) \
    if ((a) <= (b)) \
      IM_err (IMERR_LOGDIFF, " inside LogDiff() macro, a %lf  b %lf",(a),(b)); \
    if ((a) - (b) < 7.0978271289338397e+02) \
      (v) = (b) + log (exp ((a) - (b)) - 1.0); \
    else \
      (v) = (a)

/* Macro LogSum2 approximately computes the value of v of the following
 * equation:
 * exp(v) = exp(a) + exp(b) */
#define LogSum2(v,a,b) \
    if ((a) > (b)) \
      if ((a) - (b) < 7.0978271289338397e+02) \
        (v) = (b) + log (exp ((a) - (b)) + 1.0); \
      else \
        (v) = (a); \
    else \
      if ((b) - (a) < 7.0978271289338397e+02) \
        (v) = (a) + log (exp ((b) - (a)) + 1.0); \
      else \
        (v) = (b)


/******************************************/
/****** ENUMERATED THINGS ************/
/******************************************/
/*  jhey personal options using -jh flag
	SWINPUTOPTION   alternate way of inputing SW data when there is just one SW portion for a given genealogy 
    WRITEMIGRATIONNAME  -jh1filename  (i.e. the numberid symbol for WRITEMIGRATIONNAME is immediately followed by 
      two integers,  the from and to population numbers,  and then immediately followed by the filename.
      e.g. '-jh104migrate0to4.out')
    GHOSTISLAND  used with a ghost population.  There is split time for the ancestral population, instead
      the ghost and the ancestor of sampled populations are two islands with gene flow 
    */
enum
{ SWINPUTOPTION = 0,WRITEMIGRATIONNAME=1,GHOSTISLAND = 2,JHOPTIONNUMBER };

/* assignment related options 
 * Usage:
 *   -a1: assignment only
 *   -a2: DNA barcoding only
 *   -a1 -j6 : Assignment with no-migration
 *   -a2 -j6 : DNA Barcoding with no-migration
 *   -a13 : Assignment with island model
 *   -a23 : DNA Barcoding with island model
 *
 *   -a12: relabel update of assignment with population tree model
 *   -a012: turn on checking genealogy integrity additional to option 1,2
 *   -a124: print assignment proportion of a single unknown gene
 *   -a125: relabel update of assignment with island model
 *   -a127: print out the STRUCTURAMA input file additional to option 1,2
 * assignmentoptions[POPULATIONASSIGNMENT_NUMBER + 1];
 * POPULATIONASSIGNMENTASN
 * POPULATIONASSIGNMENTBAR
 *
 * POPULATIONASSIGNMENTCHECKPOINT: turns on genealogy integrity checks,
 * POPULATIONASSIGNMENT: turns on assignment,
 * POPULATIONASSIGNMENTRELABEL: relabel update
 * POPULATIONASSIGNMENTBF: Beerli and Felsenstein update
 * POPULATIONASSIGNMENTASSIGNED: for DNA Barcoding (only with single-locus data)
 * POPULATIONASSIGNMENTINFINITE: Island model
 * PRINTSTRUCTURAMA: generates an input file for Structurama */
enum
{
  POPULATIONASSIGNMENTCHECKPOINT     = 0,
  POPULATIONASSIGNMENTASN            = 1,
  POPULATIONASSIGNMENTBAR            = 2,
  POPULATIONASSIGNMENTINFINITE       = 3,
  POPULATIONASSIGNMENT               = 4,
  POPULATIONASSIGNMENTRELABEL        = 5,
  POPULATIONASSIGNMENTBF             = 6,
  POPULATIONASSIGNMENTASSIGNED       = 7,
  POPULATIONASSIGNMENTLOCAL          = 8,
  PRINTSTRUCTURAMA                   = 9,
  JCMODEL                            = 10,
  POPULATIONTREEWCH                  = 11,
  POPULATIONASSIGNMENT_NUMBER
};

/* calcoptions -c
    DONTCALCLIKELIHOODMUTATION - don't calculate p(Data|G)  if set to 1 then data likelihood functions return a constant
	MUTATIONPRIORRANGE - include pior ranges on mutation rate scalrs, as included in input file
    FINDJOINTPOSTERIOR - find joint posterior for full model and/or nested models loaded in nestedmodelfile
    USEPRIORFILE - read a file with prior's for each parameter 
*/
enum
{ 
  DONTCALCLIKELIHOODMUTATION = 0, 
  MUTATIONPRIORRANGE = 1,
  FINDJOINTPOSTERIOR = 2,
  USEPRIORFILE = 3,
  /* 5/19/2011 JH adding thermodynamic integration */
  CALCMARGINALLIKELIHOOD = 4,

  CALOPTIONSNUMBER = 5
  /*, FINDJOINTPEAK, PRINTPEAKFINDDETAILS */  
};

/*  modeloptions  using the -j flag
	NOMIGBETWEENNONSISTERS  - set migration to zero between non-sister populations
	SINGLEMIGRATIONBOTHDIRECTIONS
	MIGRATIONBETWEENSAMPLED  - migration only between sampled populations
	ADDGHOSTPOP - add a non-sampled ghost population to the model 
	PARAMETERSBYPERIOD  - every population size and migratino parameter applies for only 1 period
	NOMIGRATION  - no migration
	EXPOMIGRATIONPRIOR - exponential distribution on migration prior -m gives mean value 
    ANCESTRALPOPSIZESSETTOLARGESTSAMPLED - ancestral population size parameter is set to be the same as its largest sampled descendant population
	SINGLEPOPULATION - single population, no migration, no population tree, just a single theta      ---> replaced by npops == 1
    ONEMIGRATIONPARAMETER - just one migration parameter for the entire model  
	*/
enum
{ 
  NOMIGBETWEENNONSISTERS = 1, 
  SINGLEMIGRATIONBOTHDIRECTIONS,
  MIGRATIONBETWEENSAMPLED, 
  ADDGHOSTPOP, 
  PARAMETERSBYPERIOD,
  NOMIGRATION, 
  EXPOMIGRATIONPRIOR, 
  ANCESTRALPOPSIZESSETTOLARGESTSAMPLED,
  ONEMIGRATIONPARAMETER,
  MODELOPTIONSNUMBER
};


/* outputoptions -p
   DONTPRINTASCIITREND - don't print ascii trend lines
   DONTPRINTASCIICURVE - don't print ascii curves
   PRINTTMRCA - print a table of the TMRCA distribution 
   PRINTDEMOGHIST - print out distributions on demographic scales - requires mutation rates and generation times 
   THISTDIVIDEBYPRIOR - print histogram table of splittime values, divided by the prior - only used when splitting rate parameter is not used 
   POPMIGPARAMHIST - print the histograms and plots for the population migration rate parameters 
   PARAMGREATERTHAN - print probabilities of pairwise greater than probabilities 
   MIGRATEHIST - print out distributions of numbers and times of migration events
   PRINTJOINTTEST  - print the joint estimate for splitting times,  for more then 1 splitting time,  not in LOADRUN mode

 
*/
enum
{ 
  DONTPRINTASCIITREND, 
  DONTPRINTASCIICURVE, 
  PRINTTMRCA, 
  PRINTDEMOGHIST,
  THISTDIVIDEBYPRIOR, 
  POPMIGPARAMHIST, 
  PARAMGREATERTHAN,
  MIGRATEHIST,
  PRINTJOINTTEST,
  OUTPUTOPTIONSNUMBER
};


/* runoptions -r 
	LOADRUN  - load genealogy and t information saved in a previous run 
	DONTSAVEGENEALOGIES  - save gtree information for analytical integration  
 SAVEMCSTATEFILE - save the state of the markvon chain in a file
 LOADMCSTATE - load a previously saved markov chain state 
 PRINTMUTATIONUPDATESTOSCREEN - print all of the mutation related updates (scalars, kappas, A values) and values to stdout during the run 
    PRINTBURNTREND - print trend blot during burnin
*/
enum
{ 
  LOADRUN = 0,
  DONTSAVEGENEALOGIES,
  SAVEMCSTATEFILE,
  LOADMCSTATE,
  PRINTMUTATIONUPDATESTOSCREEN,
  PRINTBURNTREND
};

/* mutation models 
 * ---------------
 * INFINITESITES
 * HKY
 * STEPWISE    - one or more linked stepwise loci, no IS portion
 * JOINT_IS_SW - one IS portion and one or more STEPWISE portions
 * JC          - Jukes-Cantor model
 */
enum
{ 
  INFINITESITES, 
  HKY, 
  STEPWISE, 
  JOINT_IS_SW,
  JC,
  IM_MODEL_NUMBER
};

/*  different ways that the length of the run can be specified using the '-L' flag
	TIMESTEPS - time counted in # of steps
    TIMEINF  - run until first charcater in 'IMrun' is not 'y' 
	INDEFINITE - used for burnin */
enum
{ 
  TIMESTEPS, 
  TIMEINF, 
  INDEFINITE 
};

/* heating modes */
enum
{ 
  HLINEAR, 
  /*HTWOSTEP,  */
  HGEOMETRIC,
  /* 5/19/2011 JH adding thermodynamic integration */
  /* 8/23/2011  remove HTWOSTEP from use */
  HEVEN
};

/* kinds of acceptance of genealogy update */
enum
{
  IM_UPDATE_GENEALOGY_ANY,
  IM_UPDATE_GENEALOGY_TOPOLOGY,
  IM_UPDATE_GENEALOGY_TMRCA,
  /* IM_UPDATE_GENEALOGY_COVAR,  stopped using 12/15/09 JH */ 
  IM_UPDATE_GENEALOGY_NUMBER
};

/* update method of split time */
enum
{
#ifdef DO_NWUPDATE
  IM_UPDATE_TIME_NW,
#endif
#ifdef DO_RY1UPDATE
  IM_UPDATE_TIME_RY1,
#endif
  IM_UPDATE_TIME_NUMBER
};

/* update method of tree */
enum
{
  IM_UPDATE_TREE_WCH,
  IM_UPDATE_TREE_NUMBER
};

/* update method of assingment */
enum
{
  IM_UPDATE_ASSIGNMENT_RELABEL,
  IM_UPDATE_ASSIGNMENT_BF,
  IM_UPDATE_ASSIGNMENT_NUMBER
};

/* update method of assingment */
enum
{
  IM_ASSIGNRECORDINTERVAL = 10, /* interval to record assignment: -l# * 10 */
  IM_WINDOWSIZE_ASSIGNSPLITA = 20,      /* window size for changet_A: two periods / 20 */
  IM_WINDOWSIZE_ASSIGNSPLITAO = 10,     /* window size for changet_Ao: two periods / 10 */
  IM_NDELAYASSIGNMENTGENEALOGY = 200,   /* 200 for slide */
  IM_NDELAYASSIGNMENTSPLIT = 1,
  IM_SWD_BOTH = 0,
  IM_SWD_LEFT = 1,
  IM_SWD_RIGHT = 2,
  IM_SWD_NEITHER = 3,
  IM_COPYEDGE_EDGE = 0,
  IM_COPYEDGE_SISEDGE = 1,
  IM_COPYEDGE_DOWNEDGE = 2,
  IM_EDGE_LEFT = 0,
  IM_EDGE_RIGHT = 1,
  IM_EDGE_DOWN = 2,
  IM_MAXTRIES_RELABEL = 100
};

typedef char strn[PARAMSTRLEN];
typedef char strnl[UPDATELABELLEN];
/******************************************/
/***** STRUCTURES  ************/
/******************************************/

/* extendnum is used for calls to eexp() */
struct extendnum
{
  double m;
  int z;
}; 

/* eevent contains info needed to calculate the mean and variance, correlations and autocorrelations over course of run */
struct eevent
{
  double s;                     /* sum of times */
  double ss;                    /* in case need to sum a second variable */
  double s2;                    /* sum of squares of times */
  int n;                        /* number of events */
};

typedef struct eevent im_eevent;

/* for calculating autocorrelations */

struct autoc
{
  struct eevent cov;
  struct eevent var[2];
  double vals[AUTOCNEXTARRAYLENGTH];
};

typedef struct autoc im_autoc;

// used for calculating 90% HPD intervals 
struct hlists
{
  double v;
  double p;
};

typedef struct hlists im_hlists;

// terms needed for loci with HKY mutation model
struct hkyinfo
{
  double **frac;
  double **newfrac;
  double *scalefactor;
  double *oldscalefactor;
};

typedef struct hkyinfo im_hkyinfo;

/*Main data structures: edge, locus, parameter */
/* an edge is a branch in the genealogy
each edge gets a number.  the n sequences are numbers 0 thru n-1.  The remaining n-1 edges are numbered after that.
For edge i.
up[2] - contains the numbers of the edges to which edge i connects to  (i.e. the numbers of its descendants in the genealogy).
down - contains the number of the edge to which edge i connects down to  (i.e. the immediate ancestor in the genealogy
time - contains the time at the bottom of the edge, that is the time at the top of the ancestral edge (down). 
      the root edge has a negative value for time. 
*mig - pointer to arrya that contains the times of migration, and identity of pops migrated to, on edge i, these times are on the same scale as time, 
   and thus must be less than the population splitting time
cmm - current length of the memory available to hold migration times 
mut - used for labeling under INFINITESITES
A - an allele state for stepwise or other allelic model.  This is the state of the node at the top.
dlike - the likelihood of the distance between A, the allele state at the top of the node, and the allele state at the top
of the down edge 
pop - contains the population that the edge is in at its top, which may be different than 
    the one it is in at the bottom, depending on migration.
frac, newfrac, scalefactor, oldscalefactor - used for the HKY mutation model
*/
struct migstruct
{
  double mt; /* for the time                                                */
  int    mp; /* for the population a migration went to (in the coalescent); */
};

typedef struct migstruct im_migstruct;

/* struct edgmiginfo holds info about the genealogy edges being updated 
 * info is different than what is easily available in the genealogy itself
 * particularly in the arrays mtimeavail and mp which hold the times and number
 * of migration events in each period for the branch.
 * Comment on [[e]]: time period in which the edge ends  e >= b,  DO NOT CONFUSE
 * these values with poptree b and e  which refer to populations.
 *
 * When we move a branch, we could move either only the branch or the branch
 * with its sister. We may know whether a sister branch is involved in a branch
 * using member edgeid. If edgeid for a sister branch is negative, then we may
 * assume that the moving branch is not attached to a root.
 *
 * We may set the following members before using an edgemiginfo variable: li,
 * sisid, edgeid, upt, dnt, pop, temppop, fpop. By calling function 
 * [[fillmiginfoperiods]] we determine the values of b, e, mtimeavail, and 
 * mtall. We set mp, mpall, and mig by migration simulation. 
 * During migration simulations it is necessary to determine which population
 * the edge is currently in.  temppop is used for this */
struct edgemiginfo
{
  /* by manual */
  int li;
  int sisid;
  int edgeid;
  double upt;
  double dnt;
  int pop;
  int temppop;
  int fpop;
  /* by fillmiginfoperiods */
  int b;                        /* time period in which the edge begins */
  int e;
  double *mtimeavail;
  double mtall;
  /* by migration simulation or copying */
  int *mp;
  int mpall;
  struct migstruct mig[ABSMIGMAX];
};

typedef struct edgemiginfo im_edgemiginfo;

struct edge
{
  /* key edge members */
  int up[2];             /* daughter edge IDs: -1 for leaves       */
  int down;              /* down edge ID: -1 for root edge         */
  double time;           /* time at the top of edge                */
  struct migstruct *mig; /* migration events                       */
  int cmm;               /* current size of mig array              */
  /* supplementary members */
  int pop;               /* population the edge is in at its top   */
  int mut;               /* number of mutations on the edge        */
  int *A;
  double *dlikeA;
  struct hkyinfo hkyi;   /* only use if the mutation model for the */
                         /* locus is HKY                           */
  int nmig;              /* the number of migration events in mig  */
  int i;                 /* identifier of individual whose         */
  int ei;                /* index at gtree for saving purposes     */
                         /* changed from gi to ei  1/9/09 is this  */
                         /* used?                                  */
  char exist;            /* 'F' if edge is detached from a         */
                         /* genealogy                              */
  /* variables used for assignment */
  int *seq;              /* sequences at the top of edge           */
  /* Variable population tree due to Yong's */
  im_migstruct *fmig;    /* with fictitious migration events       */
  int fcmm;
  int fpop;
};

typedef struct edge im_edge;

/* struct gtreeevent  is used in treeweight() which sorts all of the events in a genealogy */
struct gtreeevent
{
  double time;
  int pop;     /* population the edge is in when the event happens        */
  int topop;   /* population the edge goes to if it is a migration event  */
  int cmt;     /* coalesce =0 or  migrate=1  or process reaches divT = -1 */
  int periodi; /* period the event is in                                  */
               /* noticed 2/16/09  could probably delete periodi as it    */
               /* seems not to get used                                   */
};

typedef struct gtreeevent im_gtreeevent;

struct upairlist
{
  int l;
  int u;
};

typedef struct upairlist im_upairlist;

struct priorvalues
{
  double min;
  double max;
  double mean;
};

typedef struct priorvalues im_priorrange;

struct plotpoint                /* (x axis are the possible values,  y axis is the counts) */
{
  double x;
  double y;
};

typedef struct plotpoint im_plotpoint;

struct update_rate_calc
{
  double accp;
  double tries;
};

typedef struct update_rate_calc im_update_rate_calc;

struct acceptancerate
{
  struct update_rate_calc *au_perturb;
};

typedef struct acceptancerate acceptancerate_t;

/* weightposition holds information on which values in the cc, fc, mc 
 * and fm arrays of genealogy_weights to sum for doing the integration */
struct weightposition
{
  int n;                        /* # of terms summed to calculate the weight                 */
  int *p;                       /* period - array of length n                                */
  int *r;                       /* row  - array of length n                                  */
  int *c;                       /* col  - array of length n (only use for migration weights) */
};

typedef struct weightposition im_weightposition;

struct i_param
{
  struct priorvalues pr;
  struct plotpoint *xy;         /* if used points to an array of gridize 
                                 * elements, each a plotpoint */
  strn str;
  int b;                        /* period when this parameter first appears */
  int e;                        /* period when this ends */
  struct weightposition wp;     /* info on where the weights are for 
                                 * calculating the integration */
};

typedef struct i_param im_i_param;

/* 1/8/09  new structure 
struct value_record  for recording stuff needed for posterior plots, trends and ESS
char  str[PARAMSTRLN];  name
struct plotpoint *xy;   recorded values that are can later be plotted,an array of gridize elements, each a plotpoint 
double *trend;    record of trend 
struct  autoc ac[AUTOCTERMS];  autocorrelation calculations
int beforemin;   number of values  below the prior
int aftermax;    number of values  above the prior
*/
struct value_record
{
  strnl str;
  strn strshort;
  int do_xyplot;
  int do_logplot;
  struct plotpoint *xy;
  struct priorvalues plotrange;
  double plotrescale;           // 1 unless needed for some reason
  int do_trend;
  double *trend;
  int do_autoc;
  struct autoc ac[AUTOCTERMS];
  int beforemin;
  int aftermax;
};

typedef struct value_record im_value_record;

/* 
1/9/09
 struct chainstate_record_updates_and_values
 replaces struct mc_param
 it is to be used in a way that allows mc_param  variables to be removed from the struct chain and from struct locus 
 it is intended to allow some simplification of code 

 this structure contains a bunch of stuff related to something that can be measured about the state of the Chain 0 
 examples include splitting times, mutation rate scalars,  assignment, measures of assignment,  features of genealogies etc etc 
 An instance of chainstate_record_updates_and_values can hold:
	info on priors and updating window width
	numbers of update types
	acceptance rates for each type
	pointers to another structure (value_record) that hold information on recorded values 

 components of chainstate_record_updates_and_values:
 char str[PARAMSTRLEN]  name of the thing being update 
 struct priorvalues pr;  
 double win    window width, may be used for updating
 int num_uptypes   number of different update types
 char **upnames   array of length num_uptypes of names of different update types
 struct update_rate_calc *upinf   - array of length num_uptypes of calculations of update rates 
 int num_vals    Number of different measurements made -  for some things this will be zero as often we just want to know about update rates, not values
 struct  value_record  *v    pointer to array of length num_vals 
*/
struct chainstate_record_updates_and_values
{
  strn str;
  struct priorvalues pr;
  double win;
  int num_uptypes;
  strnl *upnames;
  struct update_rate_calc *upinf;
  int num_vals;
  struct value_record *v;
};

typedef struct chainstate_record_updates_and_values
  im_chainstate_record_updates_and_values;

struct chainstate_updateonly
{
  strn str;
  int num_uptypes;
  strnl *upnames;
  struct update_rate_calc *upinf;
};

typedef struct chainstate_updateonly im_chainstate_updateonly;


struct popedge
{
  int numup;                    /* number of descendant nodes                 */
  int *up;                      /* the > 2 identities of the descendant nodes */
  int down;                     /* the down node                              */

  /* int ti; the time interval at the bottom of the branch    */
  double time;                  /* the time on the branch to the down node    */
  int b;                        /* time period in which the population starts */
  int e;                        /* time period that begins where the population 
                                 * ends, e > b    e == -1  if the population is 
                                 * the root */
};

typedef struct popedge im_popedge;

struct probcalc
{
  double *qintegrate;           /*  integration of theta terms                   */
  double *mintegrate;           /* integration of migration terms                */
  double pdg;                   /* probability of data given all genealogies     */
  double probg;                 /* prior prob of current genealogy  based on 
                                 * integration over parameter                    */
};

typedef struct probcalc im_probcalc;

/* genealogy_weights holds the information needed to do integrations 
 * of p(Genealogy) for the coalescent-related quantities (cc, hcc, fc),  
 * there is just an array with positions corresponding to population numbers 
 * for the migration-related quantitites there is a 3d array,  
 * with each layer being a period in the poptree, and within each period 
 * there is a 2D array of migration rates.   
 * The indexing of these does not follow population numbers,  
 * but rather follows the position of the population numbers that is 
 * given in plist */
struct genealogy_weights
{
  int **cc;                     /* coalescent counts for pops in period i            */
  double **hcc;                 /* inheritance weights for population in period i    */
  double **fc;                  /*coalescent weights for pops in period i            */
  int ***mc;                    /* migration counts between populations in period i, 
                                 * the positions in layer k of mc (i.e. mc[k] follow 
                                 * the same order and listing in plist[k]            */
  double ***fm;                 /* migration weights between populations in period i */
};

typedef struct genealogy_weights im_genealogy_weights;

/*  1/9/09  new struct locus 
*/

struct locus
{
  char name[NAMELENGTH];
  int pairs[MAXGENES];
  char gNames[MAXGENES][MAXLENGENENAME];        /* upto 100 gene names of at most 9-character string */
  int numgenesknown;
  int numgenesunknown;
  int ii[MAXGENES];  // what is this for ??  JH  7/9/09 - used for assignment, contains individual ids
  int numgenes;
  int samppop[MAXPOPS];
  int numlines;
  int model;                    /* the overall mutation model for the locus */
  double hval;
  int nlinked;                  /* # of linked portions  = 1 + # linked SW portions */
  int nAlinked;                 /* # of linked SW portions */
  int **A;                      /* points to microsat allele values, if nlinked > 1,  2D array */
  int maxA[MAXLINKED];          /* array of maximum allele size */
  int minA[MAXLINKED];          /* array of minimum allele size */
  int umodel[MAXLINKED];        /* array of locus specific model */
  int numbases;
  int numsites;
  int totsites;       /* added this 5/15/09,  fixed bad HKY bug, not sure why it was missing*/
  int **seq;
  int *mult;
  int *badsite;

  // mutation rate per year values and priors
  double uperyear_vals[MAXLINKED];
  struct priorvalues uperyear_prior[MAXLINKED];
  int uii[MAXLINKED];

  // records for mutation scalars 
  struct chainstate_record_updates_and_values *u_rec;
  struct chainstate_record_updates_and_values *kappa_rec;

  // records for microsat allele state updates 
  struct chainstate_record_updates_and_values *A_rec;

  // records for genealogy measurements 
  struct chainstate_record_updates_and_values *g_rec;

  struct chainstate_record_updates_and_values *a_rec;
};

typedef struct locus im_locus;

/* renamed form struct locus to struct genealogy  on 1/9/09
 moved a bunch of stuff out to the new struct locus

*gtree  pointer to the genealogy
root - edge number of root 
mignum - current # of migration events in gtree
roottime - time of root
length - length of gtree
tlength -  length of gtree for part with time < t
pdg - current total P(X|G)
pdg_a array of P(X|G) for each part of locus
*gtree  pointer to array of branch info
genealogy_weights - array of info that is saved when the genealogy is sampled
pi - used for HKY model
*/

struct genealogy
{
  /* genealogy stuff - different for each chain */
  struct edge *gtree;
  int root;             /* Should be updated in each genealogy update.                            */
  int mignum;           /* Is computed in treeweight.                                             */
  double roottime;      /* Should be updated in each genealogy update.                            */
  double length;        /* Is computed in treeweight.                                             */
  double tlength;       /* Is computed in treeweight.                                             */
  double pdg;           /* the total probability of the data given the genealogy (sum of pdg_a)   */
  double *pdg_a;        /* points to an array of of length nlinked  used for multiple linked loci */
  double *uvals;        /* points to an array of length nlinkded - the mutation scalar values     */
  double kappaval;      /* value of HKY parameter if needed */
  double pi[4];         /* for HKY model */
  double hilike;
  struct genealogy_weights gweight;
  double asn;
  int *mut;

  /* This is for computing Levine's likelihood. */
  /*UByteP *Pi;
  UByteP *N;
  UByteP *Nl;
  UByteP *Nr;
  UByteP **UV;
  UByteP **Ui;
  UByteP **Vi;
  UByteP **Ni;
  UByteP **Ri;
  UByteP **Zi;
  UByteP **URi; */
};

typedef struct genealogy im_genealogy;

struct chain
{
  char poptreestring[POPTREESTRINGLENGTHMAX];
  SET periodset[MAXPERIODS];
  int addpop[MAXPERIODS];
  int droppops[MAXPERIODS][2];
  int **plist;
  struct popedge *poptree;
  int rootpop;
  double *tvals;                // array of time values
  char name[10];                /* chain name or id */

  struct genealogy_weights allgweight;  // accumulated stuff used to calculate p(G) 
  struct probcalc allpcalc;     // prob(Data|G) and prob(G) stuff
  // removed  1/9/09
  //struct mc_param *t;           // array of split time info
  // added 1/9/09  just a simple array of time splittimes 
  struct genealogy *G;          // points to an array of genealogy 

  /* ASSIGNMENT: nasn - the number of individuals assigned to populations. */
  int *nasn;
};

typedef struct chain im_chain;

/******************************************/
/*******     GLOBAL VARIABLES     *********/
/******************************************/

/* SC: I would like to prefix global variable names with `g' later. */
 
/* if GLOBVARS is defined then gextern is ignored. This
causes the corresponding variable declarations to be definitions.
GLOBVARS is only defined in ima_main.c.  If GLOBVARS is not defined,
as is the case at the beginning of the other files, then gextern
gets replaced by extern and the variable declarations are simply
declarations */



#ifdef GLOBVARS
#define gextern
#else /*  */
#define gextern extern
#endif /*  */

gextern int lastgenealogysaved; 
gextern struct i_param *itheta;
gextern struct i_param *imig;
gextern double thetaprior, mprior, tprior;
gextern double tperiodpriors[MAXPOPS-1];
gextern struct chain **C;       //points to an array of pointers to chains 
//added 1/9/09 
gextern struct chainstate_record_updates_and_values *T;
//added 1/9/09 
gextern struct locus *L;
        // record for likelihood measurements
gextern struct value_record *lpgpd_v;
  // used if outputoptions[MIGRATEHIST]

/* 8/26/2011 */
gextern struct value_record **migration_counts;
//gextern struct value_record **migration_counts_times;

//added 3/17/09
gextern im_chainstate_updateonly *Cupinf;
gextern int step;
gextern int nurates;
gextern int nkappas; 
gextern int numchains;
//AS: adding numprocesses as a global
gextern int numprocesses;
//AS: adding swapflag
gextern int swapflag;
//AS: adding extendor
gextern int extendor;
gextern int nloci;
gextern int npops;
gextern int numtreepops;        /* number of distinct populations in the tree */
gextern int numpopsizeparams;   /* number of distinct population size parameters in the model */
gextern int nummigrateparams;   /* number of distinct migration rate parameters in the model */
gextern int nummigdirs;
gextern int jheyoptions[JHOPTIONNUMBER];
gextern int assignmentoptions[POPULATIONASSIGNMENT_NUMBER];
gextern int modeloptions[MODELOPTIONSNUMBER];
gextern int calcoptions[CALOPTIONSNUMBER];
gextern int outputoptions[OUTPUTOPTIONSNUMBER];
gextern int runoptions[PRINTBURNTREND + 1];
gextern double beta[MAXCHAINS];
//AS: adding an array to store other betas for printing from multiple processors
gextern double *allbetas;
//AS:: adding betaref1, betaref2
gextern double betaref1[MAXCHAINS];
gextern double betaref2[MAXCHAINS];
gextern double gbeta;
gextern double gloglikelihood;
gextern int genealogiessaved;
gextern int somestepwise;
gextern int countuprior;
gextern int counturateperyear;
gextern struct upairlist ul[2 * MAXLOCI];       /* listing of mutation rate scalars  not clear how big to make it,  because not clear how many portions each locus will have */
gextern int lastperiodnumber;   // just numsplittimes, but it comes up a lot 
gextern int numsplittimes;      // same as lastperiodnumber, but it comes up a lot
gextern float **gsampinf;
/* position markers for gsampinf */
gextern int gsamp_ccp;
gextern int gsamp_fcp;
gextern int gsamp_hccp;
gextern int gsamp_mcp;
gextern int gsamp_fmp;
gextern int gsamp_qip;
gextern int gsamp_mip;
gextern int gsamp_pdgp;
gextern int gsamp_probgp;
gextern int gsamp_tp;
//gextern struct update_rate_calc *NW_t_upinf;
gextern double logfact[100 * ABSMIGMAX + 1];
gextern struct weightposition nomigrationchecklist;
gextern struct extendnum *eexpsum; 
/* g for global and i for integer type */
gextern int gi_largestngenes; /* the largest number of genes among all loci */
gextern int gi_largestnumsites; /* the largest number of sites among all loci */
gextern int total_numgenes;  /* sum of sample sizes across all loci */ 
gextern double max_m_from_priorfile;

///AS: adding a flag variable to indicate if 
gextern int rundone;

///AS: adding a new swap counter for across processors. I keep adding things
//to this on each processor until the end of M mode, then do an all reduce to get the actual swap counts
//gextern int **swapcount_bwprocesses;

gextern int swapcount[MAXCHAINS][MAXCHAINS];
gextern int swapcount_bwprocesses[MAXCHAINS][MAXCHAINS];
gextern int swaps_bwprocesses[MAXCHAINS][MAXCHAINS];
gextern int swaps_rec_bwprocesses[MAXCHAINS][MAXCHAINS];
//AS: taking these out because the results are just too long for greater number of chains
gextern int tempbasedswapcount[MAXCHAINS][MAXCHAINS];
gextern int tempbased_rec_swapcount[MAXCHAINS][MAXCHAINS];

//gextern int numback, numforth, backcount, forthcount;


/*
#ifdef HPDBG
// define a file to hold heap info and dumps
gextern char  heapdebugfilename[50];
gextern FILE  *heapdebugfile;
gextern _CrtMemState memstate1, memstate2, memstate3;

#endif    */


/******************************************/
/***** GLOBAL FUNCTION PROTOTYPES *********/
/******************************************/

/**** GLOBAL FUNCTIONS IN autoc.c ****/
void init_autoc_pointers (void);
void free_autoc_pointers (void);
void checkautoc (int start_autocorrelations, int burndone, int burnsteps, int currentid);
void callprintautoctable (FILE * outto/*, int step*/);

/**** GLOBAL FUNCTIONS IN build_gtree.c ****/
void makeHKY (int ci, int li, int nosimmigration);
void makeIS (int ci, int li, int nosimmigration);
void makeJOINT_IS_SW (int ci, int li, int nosimmigration);
void makeSW (int ci, int li, int nosimmigration);

/**** GLOBAL FUNCTIONS IN build_poptree.c ****/
void testtreewrite (char startpoptreestring[]);
void add_ghost_to_popstring (char poptreestring[]);
void setup_poptree (int ci, char startpoptreestring[]);

/**** GLOBAL FUNCTIONS IN calc_prob_data.c ****/
void labelgtree (int ci, int li, int edge);
double likelihoodHKY (int ci, int li, double mutrate, double kappa, int e1,
                      int e2, int e3, int e4);

void calc_sumlogk (int ci, int li, double *sumlogk);
void free_sumlogk (void);
double likelihoodIS (int ci, int li, double mutrate);
double likelihoodSW (int ci, int li, int ai, double u, double tr);
void checklikelihoodSW (int ci, int li, int ai, double u);

/**** GLOBAL FUNCTIONS IN ginfo.c ****/
void init_genealogy_weights (struct genealogy_weights *gweight);
void setzero_genealogy_weights (struct genealogy_weights *gweight);
void free_genealogy_weights (struct genealogy_weights *gweight);
void init_probcalc (struct probcalc *pcalc);
void free_probcalc (struct probcalc *pcalc);
void copy_treeinfo (struct genealogy_weights *dest,
                    struct genealogy_weights *srce);
void copy_probcalc (struct probcalc *dest, struct probcalc *srce);
void sum_subtract_treeinfo (struct genealogy_weights *addup,
                            struct genealogy_weights *addtoplus,
                            struct genealogy_weights *addtominus);
int calc_gsampinf_length (void);
//AS: Added the index of cold chain to this as on Mon Mar 31 16:20:45 EDT 2014
void savegsampinf (float *g, int z);
void sum_treeinfo (struct genealogy_weights *addup,
                   struct genealogy_weights *addto);

/**** GLOBAL FUNCTIONS IN histograms.c ****/
void printhistograms (FILE * outfile, long int recordstep,
                      double generationtime, double scaleumeaninput);
void printmigrationhistograms (FILE * outfile, long int recordstep);

/**** GLOBAL FUNCTIONS IN initialize.c ****/
void setup (char infilename[], int *fpstri, char fpstr[], char priorfilename[], int currentid);

/**** GLOBAL FUNCTIONS IN output.c ****/
//void closeopenout (FILE ** outfile, char outfilename[]);
void closeopenout (FILE ** p_to_file, char fname[]);
void checkoutfileclosed (FILE ** outfile, char outfilename[]);

void printrunbasics (FILE * outfile, int loadrun, char fpstr[],
                     int burnsteps, int recordint, int recordstep,
                     int savegenealogyint, time_t endtime, time_t starttime,
                     double hilike, double hiprob/*, int step*/);
void checkhighs (int ci, int printint, double *hilike, double *hiprob,
                 double *like/*, int step*/);

//void asciitrend (FILE * outfile, struct value_record *v, int trenddoublepoint);
void asciitrend (FILE * outfile, struct value_record *v, int trenddoublepoint,int trendspot);
void asciicurve (FILE * outfile, struct plotpoint *a, char *qlabel,
                 int logscale, int recordstep);
void printacceptancerates (FILE * outto, int numrec,
                           struct chainstate_record_updates_and_values *rec[],
                           const char *printstring);
void printacceptancerates_multichain (FILE * outto);
void printcurrentvals (FILE * outto);

void savegenealogyfile (char genealogyinfosavefilename[], FILE * genealogyinfosavefile,
                   int *lastgenealogysavedvalue, int gsampinflength);
void preparehistogram (FILE * outfile, int mode, long int recordstep,
                       double scaleumeaninput, double generationtime);
void printmigratehist (FILE * outfile, int recordstep);
void printtmrcahist (FILE * outfile, int recordstep);
void print_means_variances_correlations (FILE * outfile);


/**** GLOBAL FUNCTIONS IN readata.c ****/
void read_datafile_top_lines (char infilename[], int *fpstri, char fpstr[],
                              char startpoptreestring[]);
void readdata (char infilename[], int *fpstri,
               char fpstr[], int **numsitesIS, int currentid);
int imaInfileNpops (const char *fn);

/**** GLOBAL FUNCTIONS IN surface_call_functions.c ****/
//void eexp (double x, double *m, int *z);
double margincalc (double x, double yadjust, int pi, int logi);
// double margincalceexp (double x, double yadjust, int pi, int logi); not used as of 8/25/09

void marginalopt (int firsttree, int lasttree, double *mlval,
                  double *peakloc);
double margin95 (double mlval[], double peakloc[], int pi, int UL);
void findmarginpeaks (FILE * outfile, float *holdpeakloc);


/**** GLOBAL FUNCTIONS IN surface_search_functions.c ****/
void mnbrakmod (int ndim, int firsttree, int lasttree, double *ax, double *bx,
                double *cx, double *fa, double *fb, double *fc,
                double (*func) (int, int, int, double, int), int ifneeded);
double goldenmod (int ndim, int firsttree, int lasttree, double ax,
                  double bx, double cx, double tol_, double *xmin,
                  double (*func) (int, int, int, double, int), int ifneeded);

/**** GLOBAL FUNCTIONS IN popmig.c ****/
double calc_popmig (int thetai, int mi, double x, int prob_or_like);
double calc_pop_expomig (int thetai, int mi, double x, int prob_or_like);
void marginalopt_popmig (int firsttree, int lasttree, double *mlval, double *peakloc, int *mpop, int *mterm);
double marginpopmig (int mi, int firsttree, int lasttree, double x, int thetai);

/**** GLOBAL FUNCTIONS IN gtint.c ****/
void print_greater_than_tests (FILE * outfile);

/**** GLOBAL FUNCTIONS IN swapchains.c ****/
void setheat (double hval1, double hval2, int heatmode, int currentid);

/* CR 110929.4 get rid of extraneous args in declaration  */
int swapchains (int swaptries, int swapbetasonly, int currentid, int heatmode);
/* AS - adding this function for swapping heats only */
void swapbetas (int ci, int cj);
/* AS - adding a function to swap heats between processes */
//#ifdef MPI_ENABLED
int swapchains_bwprocesses (/*int *swapper, int *swappee, */int current_id, int step, int swaptries, int swapbetasonly,int chainduration, int burnduration,/* std::ofstream &f1, int swapA, int swapB*/ int heatmode);
//#endif
void printchaininfo (FILE * outto, int heatmode, 
                     double hval1, double hval2, int currentid);

/**** GLOBAL FUNCTIONS IN gtreeprint.c ****/
void gtreeprint (int ci, int li/*, int step , int callsource */ );
void poptreeprint (int ci /*, int step*/);


/**** GLOBAL FUNCTIONS IN update_gtree_common.c ****/
int findperiod (int ci, double t);
int nowedgepop (int ci, struct edge *gtree, double ptime);
void init_treeweight (void);
void free_treeweight (void);
void init_holdgtree (struct genealogy *g, int numgenes);
void free_holdgtree (struct genealogy *g, int numgenes);

void treeweight (int ci, int li);
void initialize_integrate_tree_prob (int ci,
                                     struct genealogy_weights *gweight,
                                     struct probcalc *pcalc);
void integrate_tree_prob (int ci, struct genealogy_weights *gweight,
                          struct genealogy_weights *holdgweight,
                          struct probcalc *pcalc, struct probcalc *holdpcalc,
                          double *holdt);
void copyfraclike (int ci, int li);
void storescalefactors (int ci, int li);
double finishSWupdateA (int ci, int li, int ai, int edge, int downedge,
                        int sisedge, int newsisedge, double u, double *Aterm);
double finishSWupdateAD (int ci, int li, int ai, int edge, int downedge,
                         int sisedge, int newsisedge, double u,
                         double *Aterm);
#if 0
//  CR 110825.1 updateA() not used and no longer compiled into code
double updateA (int ci, int li, int ai, double u, int *count);
#endif   // if 0

double updateAD (int ci, int li, int ai, double u, int *count);
void restorescalefactors (int ci, int li);

/**** GLOBAL FUNCTIONS IN update_gtree.c ****/

void init_updategenealogy (void);
void free_updategenealogy (void);
int updategenealogy (int ci, int li, int *topolchange, int *tmrcachange);

/******* global in update_t_NW **********/
void init_t_NW (void);
void free_t_NW (void);

int changet_NW (int ci, int timeperiod);

/******* global in update_t_RY **********/
void init_t_RY (void);
void free_t_RY (void);

/*   getnewt() duplicate function declaratation.  One in 
 *          update_gtree_common.h kept.
 *double getnewt (int timeperiod, double t_u_prior, double t_d_prior,
 *               double oldt, int whichupdate);
 */

int changet_RY1 (int ci, int p);

/**** GLOBAL FUNCTIONS IN update_mc_params.c ****/

int changeu (int ci, int j, int *k);
int changekappa (int ci);

/**** GLOBAL FUNCTIONS IN utilities.c ****/

void errr (int ci, int li, int i, double val1, double val2);
void nrerror (const char error_text[]);
double mylogcosh (double x);
double mylogsinh (double x);
void setseeds (int seed);
void resetseeds (int seed);     /* use for debugging, provides control over sequence of random numbers */
void unsetseeds ();
double uniform ();
double expo (double c);
double normprob (double mean, double stdev, double val);
double normdev (double mean, double stdev);
int poisson (double param, int condition);
int geometric (double p);
void hpsortmig (struct migstruct *lptr, int n);
void shellhist (struct hlists *hptr, int length);
void indexx (unsigned long n, struct gtreeevent *arr, unsigned long *indx);
void gcf (double *gammcf, double a, double x, double *gln);
void gser (double *gamser, int a, double x, double *gln);
void checktreestring (char *t);
int put_spaces_in_filepaths(char *pathstr);
/* double expint(int n, double x); */
double expint (int n, double x, int *islog);
double uppergamma (int a, double x);
double lowergamma (int a, double x);
int isemptystring(char *s);
void strdelete (char *s, int pos, int len);
void strinsert (char *dest, char *source, int pos);
void strtrunc (char *s, char c);
char *nextwhite (char *c);
int allwhitespace (char *c);
char *nextnonspace (char *textline);
void setlogfact (void);
void ieevent (struct eevent *a);        /* initialize structure */
double bessi (int n, double x);
int **alloc2Dint (int rows, int cols);
double **orig2d_alloc2Ddouble (int rows, int cols);
void orig2d_free2D (void **a, int rows);
extern void eexp (double x, double *m, int *z);
extern void checkmigt (int i, struct edge *gtree);
extern void copymig (struct migstruct *m1, struct migstruct *m2);
extern void checkmig (int i, struct migstruct **mig, int *nmig);
extern int bitran (void);   
extern int randposint (int lessthanval);

double **alt2d_alloc2Ddouble(long m, long n);
double **allocAndInit2Ddouble (int rows, int cols); /* cr 110907.1 */
void alt2d_free2D(double **x/*, long m , long n */);


/* GLOBAL Functions in multi_t_bins */
void setup_multi_t_arrays (int z);
void free_multi_t_arrays (void);
void return_joint_t (double tvals[]);
double joint_t_prob (double *tvals);

/* called by freeanymemory */
void likelihoodDG_init (int ci, int li);
void likelihoodDG_fin ();

/* freemem.c */
void freeanymemory (void);

/**** Global Functions in File: mcffile */
void writemcf (char mcffilename[]);
void readmcf (char mcffilename[]);

/**** Global Functions in File: jointfind */

/* CR 110921.1  Change type declaration of outfile parameter so that this
 * function can pass in the proper parameter type when calling closeopenout()
 */
void findjointpeaks(FILE **outfile,char *outfilename, char *nestfname,int number_of_parameters);

/**** Global Functions in File: readpriorfile*/
void readpriorfile(char priorfilename[],double *popsizepriorvals, double **mpriorvals);

/****  Global functions and structures for simple profiling, functions 
 *          are in utilities.c ****/
/* duplicate function declaration.  One in update_gtree_common.h kept.
 * void free_gtreecommon ();
 */


/* 5/19/2011 JH adding thermodynamic integration */
/**** Global Functions in File: marglike.c */
void initmarginlikecalc();
void summarginlikecalc();
double harmonicmarginlikecalc(int k);
double thermomarginlikecalc(int k);

/**** Global Functions in File: msdirent.c */
/* not needed here,  prototypes are in msdirent.h
#ifdef _MSC_VER
	DIR *opendir(const char *name);
	int closedir(DIR *dir);
	struct dirent *readdir(DIR *dir);
	void rewinddir(DIR *dir);
#endif
*/

int whichiscoldchain(void);
