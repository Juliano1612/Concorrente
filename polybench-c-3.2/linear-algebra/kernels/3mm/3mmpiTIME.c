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
double start, end, sumTime = 0, maiorTime  = -INFINITY , menorTime = INFINITY;

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
void malloca_matrix(){

	dA = (double *) malloc (sizeof(double) * NI*NI);
	dB = (double *) malloc (sizeof(double) * NI*NI);
	dC = (double *) malloc (sizeof(double) * NI*NI);
	dD = (double *) malloc (sizeof(double) * NI*NI);
	dE = (double *) malloc (sizeof(double) * NI*NI);
	dF = (double *) malloc (sizeof(double) * NI*NI);
	dG = (double *) malloc (sizeof(double) * NI*NI);

}


static
void init_array(int ni, int nj, int nk, int nl, int nm)
{
	int i, j;


	for(i = 0; i < ni; i++){
		for(j = 0; j < ni; j++){
			//A[i][j] = ((DATA_TYPE) i*j) / ni;
			dA[(i*ni)+j] = ((DATA_TYPE) i*j) / ni;

			//B[i][j] = ((DATA_TYPE) i*(j+1)) / nj;
			dB[(i*ni)+j] = ((DATA_TYPE) i*(j+1)) / nj;

			//C[i][j] = ((DATA_TYPE) i*(j+3)) / nl;
			dC[(i*ni)+j] = ((DATA_TYPE) i*(j+3)) / nl;

			//D[i][j] = ((DATA_TYPE) i*(j+2)) / nk;
			dD[(i*ni)+j] = ((DATA_TYPE) i*(j+2)) / nk;

			//E[i][j] = ((DATA_TYPE) 0);
			dE[(i*ni)+j] = ((DATA_TYPE) 0);

			//F[i][j] = ((DATA_TYPE) 0);
			dF[(i*ni)+j] = ((DATA_TYPE) 0);

			//G[i][j] = ((DATA_TYPE) 0);
			dG[(i*ni)+j] = ((DATA_TYPE) 0);
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

static
void printMatrixLine(double* m){
	int tamLinha  = 0;
	int contaLinha = 0;
	//printf("\nWR%d Linha[%d]", world_rank, contaLinha);
	for(int i = 0; i < NI*NI; i++){
		if(tamLinha == NI){
			contaLinha++;
			tamLinha = 0;
			//printf("\nWR%d Linha[%d] %.2f", world_rank, contaLinha, m[i]);
			printf("\n%.2f", m[i]);

		}else{
			printf(" %f", m[i]);
		}
		tamLinha++;

	}
	printf("\n");
}

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
		//printf("%d Concluiu MPI_Second BEGIN %d END %d\n", world_rank, begin, end);
		//printMatrixLine(dG);
}


	static
	void kernel_3mm_MPI_First(){

		int begin = world_rank * tamParte;
		int end = begin+tamParte;
		int i, j, k;

  /* E := A*B */
		for (i = begin; i < end; i++){
			for (j = 0; j < NI; j++){
				dE[(i*NI)+j] = 0;
				for (k = 0; k < NI; ++k){
					dE[(i*NI)+j] += dA[(i*NI)+k] * dB[(k*NI)+j];
				}
			}
		}
  /* F := C*D */
		for (i = begin; i < end; i++){
			for (j = 0; j < NI; j++){
				dF[(i*NI)+j] = 0;
				for (k = 0; k < NI; ++k){
					dF[(i*NI)+j] += dC[(i*NI)+k] * dD[(k*NI)+j];
				}
			}
		}

		//printf("%d Concluiu MPI_First BEGIN %d END %d\n", world_rank, begin, end);

	}



	int main(int argc, char** argv)
	{

	malloca_matrix();

	MPI_Init(NULL, NULL);
	for(int i = 0; i < 12; i++){
		start = MPI_Wtime();

		MPI_Comm_size(MPI_COMM_WORLD, &world_size);
		MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

		tamParte = ceil(NI/world_size);

		if(world_rank == 0){

			init_array (ni, nj, nk, nl, nm);


			for(int i = 1; i < world_size; i++){
				MPI_Send(dA+((i*tamParte)*NI), tamParte*NI,MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
				MPI_Send(dB, NI*NI,MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
				MPI_Send(dC+((i*tamParte)*NI), tamParte*NI,MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
				MPI_Send(dD, NI*NI,MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
			}

		}else{

			for(int i = 0; i < ni; i++){
				for(int j = 0; j < ni; j++){
					dA[(i*ni)+j] = ((DATA_TYPE) 0);
					dB[(i*ni)+j] = ((DATA_TYPE) 0);
					dC[(i*ni)+j] = ((DATA_TYPE) 0);
					dD[(i*ni)+j] = ((DATA_TYPE) 0);
					dE[(i*ni)+j] = ((DATA_TYPE) 0);
					dF[(i*ni)+j] = ((DATA_TYPE) 0);
					dG[(i*ni)+j] = ((DATA_TYPE) 0);
				}
			}

			MPI_Recv(dA+((world_rank*tamParte)*NI), tamParte*NI, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(dB, NI*NI, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(dC+((world_rank*tamParte)*NI), tamParte*NI, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(dD, NI*NI, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
				kernel_3mm_MPI_First();


		if(world_rank == 0){

			for(int i = 1; i < world_size; i++){
				MPI_Recv(dE+(i*tamParte*NI), tamParte*NI, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(dF+(i*tamParte*NI), tamParte*NI, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}


			for(int i = 1; i < world_size; i++){
				MPI_Send(dE+(i*tamParte*NI), tamParte*NI,MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
				MPI_Send(dF, NI*NI,MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
			}


		}else{

			MPI_Send(dE+((world_rank*NI)*tamParte), tamParte*NI,MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
			MPI_Send(dF+((world_rank*NI)*tamParte), tamParte*NI,MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

			MPI_Recv(dE+((world_rank*NI)*tamParte), tamParte*NI, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(dF, NI*NI, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		}

			kernel_3mm_MPI_Second();

		if(world_rank == 0){

			for(int i = 1; i < world_size; i++){
				MPI_Recv(dG+(i*tamParte*NI), tamParte*NI, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}

		//printMatrixLine(dG);
		}else{
			MPI_Send(dG+(world_rank*tamParte*NI), tamParte*NI,MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
		}

		end = MPI_Wtime();

		if(world_rank == 0){
	      double b, e;
	      for (int i = 1; i < world_size; i++) {
	        MPI_Recv(&b, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	        MPI_Recv(&e, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	        if(b < start) start = b;
	        if(e > end) end = e;
	      }
		  double auxTime = end-start;
		  sumTime += auxTime;
		  if(auxTime < menorTime) menorTime = auxTime;
		  if(auxTime > maiorTime) maiorTime = auxTime;
	      //printf("%f\n", auxTime);
    	}else{
	      MPI_Send(&start, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	      MPI_Send(&end, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	  	}
	  	if(world_rank == 0){
	  		printf("Execucao %d terminada\n", i+1);
	  	}
	}

	if(world_rank == 0){
		sumTime = sumTime - menorTime - maiorTime;
		printf("MEDIA %d computadores %f\n", world_size, sumTime/10);
	}


			    // Finalize the MPI environment.
	MPI_Finalize();




	return 0;
}
