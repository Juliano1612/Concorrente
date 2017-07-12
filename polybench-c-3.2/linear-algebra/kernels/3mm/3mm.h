/**
 * 3mm.h: This file is part of the PolyBench/C 3.2 test suite.
 *
 *
 * Contact: Louis-Noel Pouchet <pouchet@cse.ohio-state.edu>
 * Web address: http://polybench.sourceforge.net
 */
#ifndef _3MM_H
# define _3MM_H

/* Default to STANDARD_DATASET. */
# if !defined(PERFECT_DATASET) && !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET)
#  define STANDARD_DATASET
# endif

/* Do not define anything if the user manually defines the size. */
# if !defined(NI) && !defined(NJ) && !defined(NK)
/* Define the possible dataset sizes. */
#  ifdef MINI_DATASET
#   define NI 4
#   define NJ 4
#   define NK 4
#   define NL 4
#   define NM 4
#  endif

#  ifdef SMALL_DATASET
#   define NI 256
#   define NJ 256
#   define NK 256
#   define NL 256
#   define NM 256
#  endif

# ifdef PERFECT_DATASET
#  define NI 512
#  define NJ 512
#  define NK 512
#  define NL 512
#  define NM 512
# endif

#  ifdef STANDARD_DATASET /* Default if unspecified. */
#   define NI 1024
#   define NJ 1024
#   define NK 1024
#   define NL 1024
#   define NM 1024
#  endif

#  ifdef LARGE_DATASET
#   define NI 2048
#   define NJ 2048
#   define NK 2048
#   define NL 2048
#   define NM 2048
#  endif

#  ifdef EXTRALARGE_DATASET
#   define NI 4000
#   define NJ 4000
#   define NK 4000
#   define NL 4000
#   define NM 4000
#  endif
# endif /* !N */

# define _PB_NI POLYBENCH_LOOP_BOUND(NI,ni)
# define _PB_NJ POLYBENCH_LOOP_BOUND(NJ,nj)
# define _PB_NK POLYBENCH_LOOP_BOUND(NK,nk)
# define _PB_NL POLYBENCH_LOOP_BOUND(NL,nl)
# define _PB_NM POLYBENCH_LOOP_BOUND(NM,nm)

# ifndef DATA_TYPE
#  define DATA_TYPE double
#  define DATA_PRINTF_MODIFIER "%0.2lf "
# endif


#endif /* !_3MM */
