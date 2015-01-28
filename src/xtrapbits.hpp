/* $XFree86$ */
/*
 * This include file is designed to be a portable way for systems to define
 * bit field manipulation of arrays of bits.
 */
#ifndef __XTRAPBITS__
#define __XTRAPBITS__ "@(#)xtrapbits.h	1.6 - 90/09/18  "

/*****************************************************************************
Copyright 1987, 1988, 1989, 1990, 1994 by Digital Equipment Corporation, 
Maynard, MA

Permission to use, copy, modify, and distribute this software and its 
documentation for any purpose and without fee is hereby granted, 
provided that the above copyright notice appear in all copies and that
both that copyright notice and this permission notice appear in 
supporting documentation, and that the name of Digital not be
used in advertising or publicity pertaining to distribution of the
software without specific, written prior permission.  

DIGITAL DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING
ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL
DIGITAL BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR
ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS
SOFTWARE.

*****************************************************************************/
/*
 *
 *  CONTRIBUTORS:
 *
 *      Dick Annicchiarico
 *      Robert Chesler
 *      Dan Coutu
 *      Gene Durso
 *      Marc Evans
 *      Alan Jamison
 *      Mark Henry
 *      Ken Miller
 *
 */
typedef unsigned char *UByteP;  /* Pointer to an unsigned byte array */
#define BitsInByte    8L        /* The number of bits in a byte */

#define BitInByte(bit)	        /* Returns the bit mask of a byte */ \
    (1L << (((bit) % BitsInByte)))

#define BitInWord(bit)	        /* Returns the bit mask of a word */ \
    (1L << (((bit) % (BitsInByte * 2L))))

#define BitInLong(bit)	        /* Returns the bit mask of a long */ \
    (1L << (((bit) % (BitsInByte * 4L))))

#define ByteInArray(array,bit)	/* Returns the byte offset to get to a bit */ \
    (((UByteP)(array))[(bit) / BitsInByte])

#define BitIsTrue(array,bit)    /* Test to see if a specific bit is True */ \
    (ByteInArray(array,bit) & BitInByte(bit))

#define BitIsFalse(array,bit)   /* Test to see if a specific bit is False */ \
    (!(BitIsTrue(array,bit)))

#define BitTrue(array,bit)      /* Set a specific bit to be True */ \
    (ByteInArray(array,bit) |= BitInByte(bit))

#define BitFalse(array,bit)     /* Set a specific bit to be False */ \
    (ByteInArray(array,bit) &= ~BitInByte(bit))

#define BitToggle(array,bit)    /* Toggle a specific bit */ \
    (ByteInArray(array,bit) ^= BitInByte(bit))

#define BitCopy(dest,src,bit)   /* Copy a specific bit */ \
    BitIsTrue((src),(bit)) ? BitTrue((dest),(bit)) : BitFalse((dest),(bit))

#define BitValue(array,bit)     /* Return True or False depending on bit */ \
    (BitIsTrue((array),(bit)) ? True : False)

#define BitSet(array,bit,value) /* Set bit to given value in array */ \
    (value) ? BitTrue((array),(bit)) : BitFalse((array),(bit))

#define BitAssign(A,B,size) /* A = B */ \
    memcpy ((A), (B), (size) * sizeof(unsigned char))
/*
    for ( (i) = 0; (i) < (size); (i)++) \
    ((UByteP)(A))[(i)] = ((UByteP)(B))[(i)]
*/

#define BitUnion(A,B,C,size,i) /* A = B union C */ \
    for ( (i) = 0; (i) < (size); (i)++) \
    ((UByteP)(A))[(i)] = ((UByteP)(B))[(i)] | ((UByteP)(C))[(i)]

#define BitIntersection(A,B,C,size,i) /* A = B intersection C */ \
    for ( (i) = 0; (i) < (size); (i)++) \
    ((UByteP)(A))[(i)] = ((UByteP)(B))[(i)] & ((UByteP)(C))[(i)]

#define BitDifference(A,B,C,size,i) /* A = B - C, or B intersection !C */ \
    for ( (i) = 0; (i) < (size); (i)++) \
    ((UByteP)(A))[(i)] = ((UByteP)(B))[(i)] & ~((UByteP)(C))[(i)]

#define BitEmpty(v,A,size,i) /* if A is empty, v = 0  */ \
    (v) = 0; \
    for ( (i) = 0; (i) < (size); (i)++) \
      (v) += ((UByteP)(A))[(i)] | 0

#define BitZero(A,size) /* A = 0 */ \
    memset ((A), 0, (size) * sizeof(unsigned char))
/*
    for ( (i) = 0; (i) < (size); (i)++) \
      ((UByteP)(A))[(i)] = 0
*/

#define BitOne(A,size) /* A */ \
    memset ((A), '\xff', (size) * sizeof(unsigned char))

#define BitSingleton(A,a,size) /* A = {a} */ \
    BitZero((A),(size)); \
    BitTrue((A),(a))

#define BitNumberTrue(v,A,size,i) /* Returns the number of 1's */ \
    (v) = 0; \
    for ((i) = 0; (i) < (size)*8; (i)++) \
      { \
        if ((BitIsTrue((A), (i)))) \
          { \
            (v) = (v) + 1; \
          } \
      }

#define BitNTrues(v,A,size,i) /* Returns the number of 1's */ \
    (v) = 0; \
    for ((i) = 0; (i) < (size)*8; (i)++) \
      { \
        if ((BitIsTrue((A), (i)))) \
          { \
            (v) = (v) + 1; \
          } \
      }

/*
    (v) = 0; \
    for ((i) = 0; (i) < (size); (i)++) \
      { \
        (v) += BITNUMBERTRUE[((UByteP)(A))[(i)]]; \
      } 
*/
/*
    (v) = 0; \
    for ((i) = 0; (i) < (size)*8; (i)++) \
      { \
        if ((BitIsTrue((A), (i)))) \
          { \
            (v) = (v) + 1; \
          } \
      } 
*/

#define Bit1s(v,a,A,size,i) \
    (v) = 0; \
    for ((i) = 0; (i) < (size)*8; (i)++) \
      { \
        (a)[(i)] = -1; \
        if ((BitIsTrue((A), (i)))) \
          { \
            (a)[(v)] = (i); \
            (v) += 1; \
          } \
      }

/* v = 1 if N(A) is less than N(B). 0 otherwise. */
#define BitIsSmaller(v,A,B,size,i,a,b) /* N(A) < N(B) */ \
    BitNumberTrue((a),(A),(size),(i)); \
    BitNumberTrue((b),(B),(size),(i)); \
    (v) = 0; \
    if ((a) < (b)) \
      { \
        (v) = 1; \
      }

/* v = 1 if A is a subset of B. 0 otherwise. */
/* A - B = empty */
#define BitIsSubset(v,A,B,X,size,i) /* A is a subset of B */ \
    BitDifference((X),(A),(B),(size),(i)); \
    BitNTrues((v),(X),(size),(i)); \
    if ((v) == 0) \
      { \
        (v) = 1; \
      } \
    else \
      { \
        (v) = 0; \
      }

/*            
    BitDifference((C),(A),(B),(size),(i)); \
    BitNumberTrue((v),(C),(size),(i)); \
    if ((v) == 0) \
      { \
        (v) = 1; \
      } \
    else \
      { \
        (v) = 0; \
      } 
*/

/* v = 1 if A is a proper subset of B. 0 otherwise. */
#define BitIsProperSubset(v,A,B,X,size,i,a,b) /* A is a subset of B */ \
    BitIsSmaller ((v),(A),(B),(size),(i),(a),(b)); \
    if ((v) == 1) \
      { \
        BitIsSubset ((v),(A),(B),(X),(size),(i)); \
      }

/* 0010 < 0110: 0 0 -> 0, 1 1 -> 0, 0 1 -> 0, 1 0 -> 1 */
/* 1101 < 0110: 1 0 -> 0, 0 1 -> 0, 1 1 -> 0, 0 0 -> 1 */
/* A < B == ~(~A or B) and N(A) < N(B) */
/* 0010 1101 */
/* ---- -or- */
/* 0110 0110 */
/*      1111 */

#define BitPrint(A,size,i) /* print true bit */ \
    for ((i) = 0; (i) < (size)*8; (i)++) \
      { \
        if ((BitIsTrue((A), (i)))) \
          { \
            fprintf (stdout, "%d ", (i)); \
          } \
      }
    //fprintf (stderr, "\n"); 
    //printf ("\n"); 

#define BitBitPrint(A,size,i) /* print true bit */ \
    for ((i) = 0; (i) < (size)*8; (i)++) \
      { \
        printf("%d", (i)%10); \
      } \
    printf ("\n"); \
    for ((i) = 0; (i) < (size)*8; (i)++) \
      { \
        if ((BitIsTrue((A), (i)))) \
          { \
            printf ("1"); \
          } \
        else \
          { \
            printf ("0"); \
          } \
      } \
    printf ("\n");

#define BitSetNew(A,size) /* */ \
    (A) = (UByteP) malloc ((size) * sizeof(unsigned char)); \
    memset ((A), 0, (size) * sizeof(unsigned char))

#define BitSetDelete(A) /* */ \
    free ((A)); \
    (A) = NULL

#define BitSetZero(A,size) /* A = 0 */ \
    memset ((A), 0, (size) * sizeof(unsigned char))

/**  
 * A: a power set
 * size: set element size
 * n: a power set size
 * i: dummy variable 
 */
#define BitPowerSetNew(A,size,n,i) /* */ \
    (A) = (UByteP *) malloc ((n) * sizeof(UByteP)); \
    for ( (i) = 0; (i) < (n); (i)++) \
      { \
        (A)[(i)] = (UByteP) malloc ((size) * sizeof(unsigned char)); \
        memset ((A)[(i)], 0, (size) * sizeof(unsigned char)); \
      }

#define BitPowerSetZero(A,size,n,i) /* */ \
    for ( (i) = 0; (i) < (n); (i)++) \
      { \
        memset ((A)[(i)], 0, (size) * sizeof(unsigned char)); \
      }

/**  
 * A: a power set
 * n: a power set size
 * i: dummy variable 
 */
#define BitPowerSetDelete(A,n,i) /* */ \
    for ( (i) = 0; (i) < (n); (i)++) \
      { \
        free ((A)[(i)]); \
        (A)[(i)] = NULL; \
      } \
    free ((A)); \
    (A) = NULL;

/**  
 * A: a set of power sets
 * size: set element size
 * n: the size of the set of power sets
 * m: a power set size
 * i: dummy variable 
 * j: dummy variable 
 */
#define BitSetPowerSetNew(A,size,n,m,i,j) /* */ \
    (A) = (UByteP **) malloc ((n) * sizeof(UByteP *)); \
    for ( (i) = 0; (i) < (n); (i)++) \
      { \
        (A)[(i)] = (UByteP *) malloc ((m) * sizeof(UByteP)); \
        for ( (j) = 0; (j) < (m); (j)++) \
          { \
            (A)[(i)][(j)] = (UByteP) malloc ((size) * sizeof(unsigned char)); \
            memset ((A)[(i)][(j)], 0, (size) * sizeof(unsigned char)); \
          } \
      }

#define BitSetPowerSetZero(A,size,n,m,i,j) /* */ \
    for ( (i) = 0; (i) < (n); (i)++) \
      { \
        for ( (j) = 0; (j) < (m); (j)++) \
          { \
            memset ((A)[(i)][(j)], 0, (size) * sizeof(unsigned char)); \
          } \
      }

/**  
 * A: a set of power sets
 * n: the size of the set of power sets
 * m: a power set size
 * i: dummy variable 
 * j: dummy variable 
 */
#define BitSetPowerSetDelete(A,n,m,i,j) /* */ \
    for ( (i) = 0; (i) < (n); (i)++) \
      { \
        for ( (j) = 0; (j) < (m); (j)++) \
          { \
            free ((A)[(i)][(j)]); \
            (A)[(i)][(j)] = NULL; \
          } \
        free ((A)[(i)]); \
        (A)[(i)] = NULL; \
      } \
    free ((A)); \
    (A) = NULL;

#define ArrayCopy(A,B,n,i) \
    for ( (i) = 0; (i) < (n); (i)++) \
      { \
        (A)[(i)] = (B)[(i)]; \
      }


#endif /* __XTRAPBITS__ */
