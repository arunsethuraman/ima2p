/*IMa2p 2009-2015 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, and Arun Sethuraman */
#ifndef _UPDATE_GTREE_H_
#define _UPDATE_GTREE_H_
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS          /* empty */
# define __END_DECLS            /* empty */
#endif

__BEGIN_DECLS
/* DELETE THE FOLLOWING VARIABLES */
// variables for files that include update_gtree.h
/* int holddownA[MAXLINKED]; only accessed by update_gtree_common.c */
/* int medgedrop; only accessed by update_gtree_common.c */
/* int mrootdrop; only accessed by update_gtree.c */
/* int rootmove; */
/* double lmedgedrop; only accessed by update_gtree_common.c */
/* double lmrootdrop; only accessed by update_gtree.c */
/* struct genealogy_weights holdgweight_updategenealogy; only accessed by update_gtree.c */
/* struct genealogy_weights holdallgweight_updategenealogy; only accessed by update_gtree.c */
/* struct probcalc holdallpcalc_updategenealogy; only accessed by update_gtree.c */
/* prototype of functions local to update_gtree.c and update_gtree_covar.c*/
double findjointime (int ci, int slidepop, int sispop, double edgeuptime,
                     double sisuptime);

/* CR 111006.1 10/6/2011 JH   removed oldmigcount, oltlength, newmigcount 
 * and newtlength these had once been used for calculating the migration 
 * weight but this is now done using the  edgemiginfo structures passed to 
 * and from getmprob()
 */
double addmigration (int ci, int li); 

/* There is another addmigration_covar in imamp.h
void addmigration_covar (int ci, int li, int *oldmigcount, double *oldtlength,
                         double *form, double *form_inv, double *backm,
                         double *backm_inv, double migmultiplier);
*/
__END_DECLS
#endif /* _UPDATE_GTREE_H_ */
