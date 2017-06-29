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
#include <mpi.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
/* Default data type is double, default size is 4000. */
#include "3mm.h"

#define TIME() struct timespec start, finish; double elapsed;clock_gettime(CLOCK_MONOTONIC, &start);
#define ENDTIME() clock_gettime(CLOCK_MONOTONIC, &finish); elapsed = (finish.tv_sec - start.tv_sec); elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0; printf("%.f\n", elapsed * 1000000000);


int world_size, world_rank;

int id;
int ni = NI;
int nj = NJ;
int nk = NK;
int nl = NL;
int nm = NM;

int tamParte = NI;

double *dE;
double *dA;
double *dB;
double *dF;
double *dC;
double *dD;
double *dG;

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

	dA = (double *) malloc (sizeof(double) * NI*NI);
	dB = (double *) malloc (sizeof(double) * NI*NI);
	dC = (double *) malloc (sizeof(double) * NI*NI);
	dD = (double *) malloc (sizeof(double) * NI*NI);
	dE = (double *) malloc (sizeof(double) * NI*NI);
	dF = (double *) malloc (sizeof(double) * NI*NI);
	dG = (double *) malloc (sizeof(double) * NI*NI);

	for(i = 0; i < ni; i++){
		for(j = 0; j < ni; j++){
			A[i][j] = ((DATA_TYPE) i*j) / ni;
			dA[(i*ni)+j] = A[i][j];
			
			B[i][j] = ((DATA_TYPE) i*(j+1)) / nj;
			dB[(i*ni)+j] = B[i][j];
			
			C[i][j] = ((DATA_TYPE) i*(j+3)) / nl;
			dC[(i*ni)+j] = C[i][j];
			
			D[i][j] = ((DATA_TYPE) i*(j+2)) / nk;
			dD[(i*ni)+j] = D[i][j];
			
			E[i][j] = ((DATA_TYPE) 0);
			dE[(i*ni)+j] = E[i][j];
			
			F[i][j] = ((DATA_TYPE) 0);
			dF[(i*ni)+j] = F[i][j];
			
			G[i][j] = ((DATA_TYPE) 0);
			dG[(i*ni)+j] = G[i][j];
		}
	}
/*
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
	}*/
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

/*static
void printMatrixLine(double* m){
	int tamLinha  = 0;
	for(int i = 0; i < sizeof(m*sizeof(double)); i++){
		if(tamLinha == NI){
			printf("\n %f", m[i]);
			tamLinha = 0;
		}else{
			printf(" %f", m[i]);
		}
		tamLinha++;

	}
}*/

/* Main computational kernel. The whole function will be timed,
   including the call and return. */
/* Sequential Function */

static
void kernel_3mm_MPI_Second(){
		int begin = world_rank * tamParte;
		int end = begin+tamParte;
		int i, j, k;

		/* G := E*F */
		for (i = begin; i < end; i++){
			for (j = 0; j < _PB_NL; j++){
				dG[(i*NI)+j] = 0;
				for (k = 0; k < _PB_NJ; ++k){
					dG[(i*NI)+j] += dE[(i*NI)+k] * dF[(k*NI)+j];
				}
			}
		}
}


	static
	void kernel_3mm_MPI_First(){

		int begin = world_rank * tamParte;
		int end = begin+tamParte;
		int i, j, k;

  /* E := A*B */
		for (i = begin; i < end; i++){
			for (j = 0; j < _PB_NJ; j++){
				dE[(i*NI)+j] = 0;
				for (k = 0; k < _PB_NK; ++k){
					dE[(i*NI)+j] += dA[(i*NI)+k] * dB[(k*NI)+j];
				}
			}
		}
  /* F := C*D */
		for (i = begin; i < end; i++){
			for (j = 0; j < _PB_NL; j++){
				dF[(i*NI)+j] = 0;
				for (k = 0; k < _PB_NM; ++k){
					dF[(i*NI)+j] += dC[(i*NI)+k] * dD[(k*NI)+j];
				}
			}
		}

		printf("\n%d Concluiu MPI_First", world_rank);
  
	}



	int main(int argc, char** argv)
	{
  /* Retrieve problem size. */



  /* Variable declaration/allocation. */
	
  /* Initialize array(s). */
	
		
		POLYBENCH_2D_ARRAY_DECL(E, DATA_TYPE, NI, NJ, ni, nj);
		POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, NI, NK, ni, nk);
		POLYBENCH_2D_ARRAY_DECL(B, DATA_TYPE, NK, NJ, nk, nj);
		POLYBENCH_2D_ARRAY_DECL(F, DATA_TYPE, NJ, NL, nj, nl);
		POLYBENCH_2D_ARRAY_DECL(C, DATA_TYPE, NJ, NM, nj, nm);
		POLYBENCH_2D_ARRAY_DECL(D, DATA_TYPE, NM, NL, nm, nl);
		POLYBENCH_2D_ARRAY_DECL(G, DATA_TYPE, NI, NL, ni, nl);



  /* Start timer. */
  //polybench_start_instruments;
  /* Run kernel. */

	TIME()
	MPI_Init(NULL, NULL);
			    // Get the number of processes
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
			    // Get the rank of the process
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	tamParte = tamParte/world_size;

	printf("\nrank %d", world_rank);

	if(world_rank == 0){

		printf("\nGerente... rank %d", world_rank);



		init_array (ni, nj, nk, nl, nm,
		POLYBENCH_ARRAY(A),
		POLYBENCH_ARRAY(B),
		POLYBENCH_ARRAY(C),
		POLYBENCH_ARRAY(D),
		POLYBENCH_ARRAY(E),
		POLYBENCH_ARRAY(F),
		POLYBENCH_ARRAY(G));


		for(int i = 1; i < world_size; i++){
			MPI_Send(dA+(i*tamParte), tamParte*NI,MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
			MPI_Send(dB, NI*NI,MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
			MPI_Send(dC+(i*tamParte), tamParte*NI,MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
			MPI_Send(dD, NI*NI,MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
			MPI_Send(dE+(i*tamParte), tamParte*NI,MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
			MPI_Send(dF+(i*tamParte), tamParte*NI,MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
			//MPI_Send(dG+(i*tamParte), tamParte*NI,MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
		}

	}else{

		printf("\nEscravo... rank %d", world_rank);

		dA = (double *) malloc (sizeof(double) * tamParte*NI);
		dB = (double *) malloc (sizeof(double) * NI*NI);
		dC = (double *) malloc (sizeof(double) * tamParte*NI);
		dD = (double *) malloc (sizeof(double) * NI*NI);
		dE = (double *) malloc (sizeof(double) * tamParte*NI);
		dF = (double *) malloc (sizeof(double) * tamParte*NI);
		//dG = (double *) malloc (sizeof(double) * tamParte*NI);


		MPI_Recv(dA+(world_rank*tamParte), tamParte*NI, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(dB, NI*NI, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(dC+(world_rank*tamParte), tamParte*NI, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(dD, NI*NI, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(dE+(world_rank*tamParte), tamParte*NI, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(dF+(world_rank*tamParte), tamParte*NI, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	kernel_3mm_MPI_First();
			    // Finalize the MPI environment.
	MPI_Finalize();

	ENDTIME()

	return 0;
}
