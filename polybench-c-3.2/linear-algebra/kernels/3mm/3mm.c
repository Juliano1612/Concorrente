/**
 * 3mm.c: This file is part of the PolyBench/C 3.2 test suite.
 *
 *
 * Contact: Louis-Noel Pouchet <pouchet@cse.ohio-state.edu>
 * Web address: http://polybench.sourceforge.net
 */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <pthread.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
/* Default data type is double, default size is 4000. */
#include "3mm.h"

#define TIME() struct timespec start, finish; double elapsed;clock_gettime(CLOCK_MONOTONIC, &start);
#define ENDTIME() clock_gettime(CLOCK_MONOTONIC, &finish); elapsed = (finish.tv_sec - start.tv_sec); elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0; printf("%.f\n", elapsed * 1000000000);

int numThreads;

int id;
int ni = NI;
int nj = NJ;
int nk = NK;
int nl = NL;
int nm = NM;

int tamParte = NI;

double **dE;
double **dA;
double **dB;
double **dF;
double **dC;
double **dD;
double **dG;

pthread_barrier_t barrier;


/* Array initialization. */
#include <unistd.h>
static
void init_array(int ni, int nj, int nk, int nl, int nm,
		DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj),
		DATA_TYPE POLYBENCH_2D(C,NJ,NM,nj,nm),
		DATA_TYPE POLYBENCH_2D(D,NM,NL,nm,nl),
		DATA_TYPE POLYBENCH_2D(E,NJ,NM,nj,nm),
		DATA_TYPE POLYBENCH_2D(F,NJ,NM,nj,nm),
		DATA_TYPE POLYBENCH_2D(G,NJ,NM,nj,nm))
{
  int i, j;

	dA = (double **) malloc (sizeof(double) * NI);
	dB = (double **) malloc (sizeof(double) * NI);
	dC = (double **) malloc (sizeof(double) * NI);
	dD = (double **) malloc (sizeof(double) * NI);
	dE = (double **) malloc (sizeof(double) * NI);
	dF = (double **) malloc (sizeof(double) * NI);
	dG = (double **) malloc (sizeof(double) * NI);

  for (i = 0; i < ni; i++){
		dA[i] = A[i];
    for (j = 0; j < nk; j++){
      A[i][j] = ((DATA_TYPE) i*j) / ni;
			dA[i][j] = A[i][j];
		}
	}
  for (i = 0; i < nk; i++){
		dB[i] = B[i];
    for (j = 0; j < nj; j++){
      B[i][j] = ((DATA_TYPE) i*(j+1)) / nj;
			dB[i][j] = B[i][j];
		}
	}
  for (i = 0; i < nj; i++){
		dC[i] = C[i];
    for (j = 0; j < nm; j++){
      C[i][j] = ((DATA_TYPE) i*(j+3)) / nl;
			dC[i][j] = C[i][j];
		}
	}
  for (i = 0; i < nm; i++){
		dD[i] = D[i];
		for (j = 0; j < nl; j++){
			D[i][j] = ((DATA_TYPE) i*(j+2)) / nk;
			dD[i][j] = D[i][j];
		}
	}
	for (i = 0; i < ni; i++){
		dE[i] = E[i];
    for (j = 0; j < nk; j++){
      E[i][j] = ((DATA_TYPE) 0);
			dE[i][j] = E[i][j];
		}
	}
	for (i = 0; i < ni; i++){
		dF[i] = F[i];
    for (j = 0; j < nk; j++){
      F[i][j] = ((DATA_TYPE) 0);
			dF[i][j] = F[i][j];
		}
	}
	for (i = 0; i < ni; i++){
		dG[i] = G[i];
    for (j = 0; j < nk; j++){
      G[i][j] = ((DATA_TYPE) 0);
			dG[i][j] = G[i][j];
		}
	}
}

/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int ni, int nl,
		 DATA_TYPE POLYBENCH_2D(G,NI,NL,ni,nl))
{
  int i, j;

  for (i = 0; i < ni; i++)
    for (j = 0; j < nl; j++) {
	fprintf (stderr, DATA_PRINTF_MODIFIER, G[i][j]);
	if ((i * ni + j) % 20 == 0) fprintf (stderr, "\n");
    }
  fprintf (stderr, "\n");
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
/* Sequential Function */
static
void kernel_3mm(int ni, int nj, int nk, int nl, int nm,
		DATA_TYPE POLYBENCH_2D(E,NI,NJ,ni,nj),
		DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj),
		DATA_TYPE POLYBENCH_2D(F,NJ,NL,nj,nl),
		DATA_TYPE POLYBENCH_2D(C,NJ,NM,nj,nm),
		DATA_TYPE POLYBENCH_2D(D,NM,NL,nm,nl),
		DATA_TYPE POLYBENCH_2D(G,NI,NL,ni,nl))
{
  int i, j, k;

  /* E := A*B */
  for (i = 0; i < _PB_NI; i++){
    for (j = 0; j < _PB_NJ; j++){
			E[i][j] = 0;
			for (k = 0; k < _PB_NK; ++k){
			  E[i][j] += A[i][k] * B[k][j];
			}
    }
	}
  /* F := C*D */
  for (i = 0; i < _PB_NJ; i++){
    for (j = 0; j < _PB_NL; j++){
			F[i][j] = 0;
			for (k = 0; k < _PB_NM; ++k){
			  F[i][j] += C[i][k] * D[k][j];
			}
		}
	}
  /* G := E*F */
  for (i = 0; i < _PB_NI; i++){
    for (j = 0; j < _PB_NL; j++){
			G[i][j] = 0;
			for (k = 0; k < _PB_NJ; ++k){
			  G[i][j] += E[i][k] * F[k][j];
			}
    }
	}

}

/* OpenMP Function */
static
void kernel_3mm_OpenMP(int ni, int nj, int nk, int nl, int nm,
		DATA_TYPE POLYBENCH_2D(E,NI,NJ,ni,nj),
		DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj),
		DATA_TYPE POLYBENCH_2D(F,NJ,NL,nj,nl),
		DATA_TYPE POLYBENCH_2D(C,NJ,NM,nj,nm),
		DATA_TYPE POLYBENCH_2D(D,NM,NL,nm,nl),
		DATA_TYPE POLYBENCH_2D(G,NI,NL,ni,nl))
{

#pragma omp parallel num_threads(numThreads)
{
	//printf("\nThreads OMP %d Var numThreads %d", omp_get_num_threads(), numThreads);
  /* E := A*B */
  #pragma omp for simd
  for (int i = 0; i < _PB_NI; i++){
    for (int j = 0; j < _PB_NJ; j++){
			E[i][j] = 0;
			for (int k = 0; k < _PB_NK; ++k){
			  E[i][j] += A[i][k] * B[k][j];
			}
    }
	}
  /* F := C*D */
  #pragma omp for simd
  for (int i = 0; i < _PB_NJ; i++){
    for (int j = 0; j < _PB_NL; j++){
			F[i][j] = 0;
			for (int k = 0; k < _PB_NM; ++k){
			  F[i][j] += C[i][k] * D[k][j];
			}
		}
	}
  /* G := E*F */
  #pragma omp for simd
  for (int i = 0; i < _PB_NI; i++){
    for (int j = 0; j < _PB_NL; j++){
			G[i][j] = 0;
			for (int k = 0; k < _PB_NJ; ++k){
			  G[i][j] += E[i][k] * F[k][j];
			}
    }
	}
}
}


static
void *kernel_3mm_PThreads(void *id)
{
  int begin = (*(int *) id) * tamParte;
  int end = begin+tamParte;
  int i, j, k;


  /* E := A*B */
	//printf("\nRODEI\n\tBegin %d End %d TamParte %d\n", begin, end, tamParte);
  for (i = begin; i < end; i++){
		//printf("\nI = %d J = %d K = %d", i, j, k);
    for (j = 0; j < _PB_NJ; j++){
			dE[i][j] = 0;
			for (k = 0; k < _PB_NK; ++k){
			  dE[i][j] += dA[i][k] * dB[k][j];
			}
    }
	}
  /* F := C*D */
  for (i = begin; i < end; i++){
    for (j = 0; j < _PB_NL; j++){
			dF[i][j] = 0;
			for (k = 0; k < _PB_NM; ++k){
			  dF[i][j] += dC[i][k] * dD[k][j];
			}
		}
	}
	pthread_barrier_wait(&barrier);
  /* G := E*F */
  for (i = begin; i < end; i++){
    for (j = 0; j < _PB_NL; j++){
			dG[i][j] = 0;
			for (k = 0; k < _PB_NJ; ++k){
			  dG[i][j] += dE[i][k] * dF[k][j];
			}
    }
	}
	return NULL;

}



int main(int argc, char** argv)
{
  /* Retrieve problem size. */

	int op = atoi(argv[2]); // 0 Sequencial, 1 OpenMP, 2 PThreads
	numThreads = atoi(argv[1]);
	//printf("\nnumThreads = %d", numThreads);
	tamParte = tamParte/numThreads;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(E, DATA_TYPE, NI, NJ, ni, nj);
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, NI, NK, ni, nk);
  POLYBENCH_2D_ARRAY_DECL(B, DATA_TYPE, NK, NJ, nk, nj);
  POLYBENCH_2D_ARRAY_DECL(F, DATA_TYPE, NJ, NL, nj, nl);
  POLYBENCH_2D_ARRAY_DECL(C, DATA_TYPE, NJ, NM, nj, nm);
  POLYBENCH_2D_ARRAY_DECL(D, DATA_TYPE, NM, NL, nm, nl);
  POLYBENCH_2D_ARRAY_DECL(G, DATA_TYPE, NI, NL, ni, nl);

  /* Initialize array(s). */
  init_array (ni, nj, nk, nl, nm,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B),
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(D),
				POLYBENCH_ARRAY(E),
				POLYBENCH_ARRAY(F),
				POLYBENCH_ARRAY(G));

  /* Start timer. */
  //polybench_start_instruments;
  /* Run kernel. */
	pthread_t threads[numThreads];
	pthread_barrier_init(&barrier, NULL, numThreads);
	int ids[numThreads];

	TIME()
	switch (op) {
		case 0: //Sequencial
			kernel_3mm (ni, nj, nk, nl, nm,
						POLYBENCH_ARRAY(E),
						POLYBENCH_ARRAY(A),
						POLYBENCH_ARRAY(B),
						POLYBENCH_ARRAY(F),
						POLYBENCH_ARRAY(C),
						POLYBENCH_ARRAY(D),
						POLYBENCH_ARRAY(G));
			break;
		case 2: // OpenMP
			kernel_3mm_OpenMP (ni, nj, nk, nl, nm,
						POLYBENCH_ARRAY(E),
						POLYBENCH_ARRAY(A),
						POLYBENCH_ARRAY(B),
						POLYBENCH_ARRAY(F),
						POLYBENCH_ARRAY(C),
						POLYBENCH_ARRAY(D),
						POLYBENCH_ARRAY(G));
			break;
		case 1: //Pthreads

			//start threads
			for(int i = 0; i < numThreads; i++){
				ids[i] = i;
				pthread_create(&threads[i], NULL, kernel_3mm_PThreads, &ids[i]);
			}
			//join threads
			for (int j = 0; j < numThreads; j++) {
				pthread_join(threads[j], NULL);
			}
			break;
		default:
			printf("\nOpcao invalida! Passe como argumento:\n[0] Sequencial\n[1] OpenMP\n[2] Pthreads");
	}
	ENDTIME()
  /* Stop and print timer. */
  //polybench_stop_instruments;
  //polybench_print_instruments();

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(ni, nl,  POLYBENCH_ARRAY(G)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(E);
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);
  POLYBENCH_FREE_ARRAY(F);
  POLYBENCH_FREE_ARRAY(C);
  POLYBENCH_FREE_ARRAY(D);
  POLYBENCH_FREE_ARRAY(G);

  return 0;
}
