/**
 * 2mm.c: This file is part of the PolyBench/C 3.2 test suite.
 *
 *
 * Contact: Louis-Noel Pouchet <pouchet@cse.ohio-state.edu>
 * Web address: http://polybench.sourceforge.net
 */
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "util.h"

/* Include pthread header */
#include <pthread.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
/* Default data type is double, default size is 4000. */
#include "2mm.h"

/* Global arguments */
int ni, nj, nk, nl;
double g_alpha, g_beta;
double **g_tmp;
double **g_A;
double **g_B;
double **g_C;
double **g_D;
int T = 1;

/* Create barrier */
pthread_barrier_t barrier;

/* Create lock */
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

/* DCE code. Must scan the entire live-out data. Can be used also to check the correctness of the output. */
static void print_array(int ni, int nl, double POLYBENCH_2D(D,NI,NL,ni,nl)) {
	for (int i = 0; i < ni; i++) {
		for (int j = 0; j < nl; j++) {
			fprintf (stderr, DATA_PRINTF_MODIFIER, D[i][j]);
			if ((i * ni + j) % 20 == 0) fprintf (stderr, "\n");
		}
	}
	fprintf (stderr, "\n");
}

/* Array initialization. */
static void init_array(int ni, int nj, int nk, int nl, double *alpha, double *beta, double POLYBENCH_2D(A,NI,NK,ni,nl), double POLYBENCH_2D(B,NK,NJ,nk,nj), double POLYBENCH_2D(C,NL,NJ,nl,nj), double POLYBENCH_2D(D,NI,NL,ni,nl), double POLYBENCH_2D(tmp,NI,NJ,ni,nj)) {
	*alpha = 32412;
	*beta = 2123;
	g_A = (double **) malloc(sizeof(double) * NI);
	g_B = (double **) malloc(sizeof(double) * NI);
	g_C = (double **) malloc(sizeof(double) * NI);
	g_D = (double **) malloc(sizeof(double) * NI);
	g_tmp = (double **) malloc(sizeof(double) * NI);
	for (int i = 0; i < ni; i++) {
		g_A[i] = A[i];
		for (int j = 0; j < nk; j++) {
			A[i][j] = ((double) i*j) / ni;
			g_A[i][j] = A[i][j];
		}
	}
	for (int i = 0; i < nk; i++) {
		g_B[i] = B[i];
		for (int j = 0; j < nj; j++) {
			B[i][j] = ((double) i*(j+1)) / nj;
			g_B[i][j] = B[i][j];
		}
	}
	for (int i = 0; i < nl; i++) {
		g_C[i] = C[i];
		for (int j = 0; j < nj; j++) {
			C[i][j] = ((double) i*(j+3)) / nl;
			g_C[i][j] = C[i][j];
		}
	}
	for (int i = 0; i < ni; i++) {
		g_D[i] = D[i];
		for (int j = 0; j < nl; j++) {
			D[i][j] = ((double) i*(j+2)) / nk;
			g_D[i][j] = D[i][j];
		}
	}
	for (int i = 0; i < ni; i++) {
		g_tmp[i] = tmp[i];
		for (int j = 0; j < nj; j++) {
			g_tmp[i][j] = tmp[i][j];
		}
	}
}


static void *kernel_2mm_pthreads(void *arg) {
	int *id = (int *) arg;
	int stripe = NI / T;
	int init = (*id) * stripe;
	int end = init + stripe;
	for (int i = init; i < end; i++) {
		for (int j = 0; j < NJ; j++) {
			g_tmp[i][j] = 0;
			for (int k = 0; k < NK; ++k) {
				g_tmp[i][j] += g_alpha * g_A[i][k] * g_B[k][j];
			}
		}
	}
	pthread_barrier_wait(&barrier);
	for (int i = init; i < end; i++) {
		for (int j = 0; j < NJ; j++) {
			g_D[i][j] *= g_beta;
			for (int k = 0; k < NJ; ++k) {
				g_D[i][j] += g_tmp[i][k] * g_C[k][j];
			}
		}
	}
}

/* Main computational kernel. The whole function will be timed, including the call and return. */
/* Parallelized using OpenMP */
static void kernel_2mm_openmp(int ni, int nj, int nk, int nl, double alpha, double beta, double POLYBENCH_2D(tmp,NI,NJ,ni,nj), double POLYBENCH_2D(A,NI,NK,ni,nk), double POLYBENCH_2D(B,NK,NJ,nk,nj), double POLYBENCH_2D(C,NL,NJ,nl,nj), double POLYBENCH_2D(D,NI,NL,ni,nl)) {
	/* D := alpha*A*B*C + beta*D */
	/* A[i][j] */
	/* i -> line */
	/* j -> column */
	#pragma omp parallel num_threads(T)
	{
		#pragma omp for simd collapse(2) schedule(guided)
		for (int i = 0; i < ni; i++) {
			for (int j = 0; j < nj; j++) {
				g_tmp[i][j] = 0;
				for (int k = 0; k < nk; ++k) {
					g_tmp[i][j] += alpha * g_A[i][k] * g_B[k][j];
				}
			}
		}
		#pragma omp for simd collapse(2) schedule(guided)
		for (int i = 0; i < ni; i++) {
			for (int j = 0; j < nl; j++) {
				g_D[i][j] *= beta;
				for (int k = 0; k < nj; ++k) {
					g_D[i][j] += g_tmp[i][k] * g_C[k][j];
				}
			}
		}
	}
}

/* Main computational kernel. The whole function will be timed, including the call and return. */
/* Original sequential code */
static void kernel_2mm(int ni, int nj, int nk, int nl, double alpha, double beta, double POLYBENCH_2D(tmp,NI,NJ,ni,nj), double POLYBENCH_2D(A,NI,NK,ni,nk), double POLYBENCH_2D(B,NK,NJ,nk,nj), double POLYBENCH_2D(C,NL,NJ,nl,nj), double POLYBENCH_2D(D,NI,NL,ni,nl)) {
	/* D := alpha*A*B*C + beta*D */
	/* A[i][j] */
	/* i -> line */
	/* j -> column */
	for (int i = 0; i < ni; i++) {
		for (int j = 0; j < nj; j++) {
			g_tmp[i][j] = 0;
			for (int k = 0; k < nk; ++k) {
				g_tmp[i][j] += alpha * g_A[i][k] * g_B[k][j];
			}
		}
	}
	for (int i = 0; i < ni; i++) {
		for (int j = 0; j < nl; j++) {
			g_D[i][j] *= beta;
			for (int k = 0; k < nj; ++k) {
				g_D[i][j] += g_tmp[i][k] * g_C[k][j];
			}
		}
	}
}


int main(int argc, char** argv) {
	/* Retrieve which code to run */
	int prog = atoi(argv[2]);

	/* Retrieve problem size. */
	ni = NI;
	nj = NJ;
	nk = NK;
	nl = NL;

	/* Number of threads */
	T = atoi(argv[1]);

	/* Vector to store threads id */
	int ids[T];

	/* Variable declaration/allocation. */
	POLYBENCH_2D_ARRAY_DECL(tmp,double,NI,NJ,ni,nj);
	POLYBENCH_2D_ARRAY_DECL(A,double,NI,NK,ni,nk);
	POLYBENCH_2D_ARRAY_DECL(B,double,NK,NJ,nk,nj);
	POLYBENCH_2D_ARRAY_DECL(C,double,NL,NJ,nl,nj);
	POLYBENCH_2D_ARRAY_DECL(D,double,NI,NL,ni,nl);

	// Initialize array(s).
	init_array(ni, nj, nk, nl, &g_alpha, &g_beta, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B), POLYBENCH_ARRAY(C), POLYBENCH_ARRAY(D), POLYBENCH_ARRAY(tmp));

	// Start timer. For sequential code.
	polybench_start_instruments;

	switch(prog) {
		case 0: {
			TIME()
			// Run sequential kernel
			kernel_2mm(ni, nj, nk, nl, g_alpha, g_beta, POLYBENCH_ARRAY(tmp), POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B), POLYBENCH_ARRAY(C), POLYBENCH_ARRAY(D));
			ENDTIME()
			break;
		}

		case 1: {
			TIME()
			// Threads declaration
			pthread_t threads[T];

			// Initialize barrier
			pthread_barrier_init(&barrier, NULL, T);

			// Start threads
			start_pthread(thread, threads, parallel_pthread);

			// Wait for threads to finish
			for(int i = 0; i < T; i++) {
				pthread_join(threads[i], NULL);
			}
			ENDTIME()
			break;
		}

		case 2: {
			TIME()
			// Run OpenMP kernel.
			kernel_2mm_openmp(ni, nj, nk, nl, g_alpha, g_beta, POLYBENCH_ARRAY(tmp), POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B), POLYBENCH_ARRAY(C), POLYBENCH_ARRAY(D));
			ENDTIME()
			break;
		}

		default: break;
	}

	/* Prevent dead-code elimination. All live-out data must be printed by the function call in argument. */
	polybench_prevent_dce(print_array(ni, nl, POLYBENCH_ARRAY(D)));

	/* Be clean. */
	POLYBENCH_FREE_ARRAY(tmp);
	POLYBENCH_FREE_ARRAY(A);
	POLYBENCH_FREE_ARRAY(B);
	POLYBENCH_FREE_ARRAY(C);
	POLYBENCH_FREE_ARRAY(D);

	return 0;
}
