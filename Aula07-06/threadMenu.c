#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <omp.h>

double **A, **B, **C, *V;
int tamanhoMatriz = 100, buscado = 2;
int numThreads = 11;
int tamParte;


void initArrays(){

	A = (double**) malloc(sizeof(double*) * tamanhoMatriz);
	B = (double**) malloc(sizeof(double*) * tamanhoMatriz);
	C = (double**) malloc(sizeof(double*) * tamanhoMatriz);
	for (int i = 0; i < tamanhoMatriz; i++) {
		A[i] = (double*) malloc(sizeof(double) * tamanhoMatriz);
		B[i] = (double*) malloc(sizeof(double) * tamanhoMatriz);
		C[i] = (double*) malloc(sizeof(double) * tamanhoMatriz);
		for (int j = 0; j < tamanhoMatriz; j++) {
			A[i][j] = rand() % 2;
			B[i][j] = rand() % 2;
			C[i][j] = 0;
		}
	}

	V = (double*) malloc(sizeof(double) *tamanhoMatriz);
	for(int i = 0; i < tamanhoMatriz; i++){
		V[i] = rand() % 10;
	}
}

void *mm_PThread(void *id){
  int begin = (*(int *) id) * tamParte;
  int end = begin+tamParte;
  int i, j, k;


  /* C := A*B */
  for (i = begin; i < end; i++){
    	for (j = 0; j < tamanhoMatriz; j++){
			C[i][j] = 0;
			for (k = 0; k < tamanhoMatriz; ++k){
			  C[i][j] += A[i][k] * B[k][j];
			}
    	}
	}
	return NULL;
}


void *vectorSearch_Pthread(void *id){
	int begin = (*(int *) id) *tamParte;
	int end = begin+tamParte;
	int i;


	for(int i = begin; i < end; i++){
		if(V[i] == buscado){
			printf("\nO elemento %d existe no vetor!", buscado);
			break;
		}
	}

	return NULL;
}

void *escolhaTarefaMatrix(void *id){

	pthread_t threads[numThreads-1];
	int ids[numThreads-1];
	for(int i = 1; i < numThreads; i++){
		ids[i] = i;
		pthread_create(&threads[i], NULL, mm_PThread, &ids[i]);
	}

	for(int j = 1; j < numThreads; j++){
		pthread_join(threads[j], NULL);
	}
	printf("\nMatrizes Muliplicadas!");
}

void *escolhaTarefaVector(void *id){
	pthread_t threads[numThreads-1];
	int ids[numThreads-1];
	for(int i = 1; i < numThreads; i++){
		ids[i] = i;
		pthread_create(&threads[i], NULL, vectorSearch_Pthread, &ids[i]);
	}

	for(int j = 1; j < numThreads; j++){
		pthread_join(threads[j], NULL);
	}
	
}



int main(){

	tamParte = tamanhoMatriz/(numThreads-1);

	int option = 0;

	do{
		printf("*********************************\n");
		printf("***Escolha a tarefa desejada*****\n");
		printf("*********************************\n");

		printf("*[1] - Multiplicação de matrizes*\n");
		printf("*[2] - Busca em vetor           *\n");
		scanf("%d", &option);
	}while(option!= 1 || option != 2);


	pthread_t threadMain;
	int id = 0;
	if(option == 1){
		pthread_create(&threadMain, NULL, escolhaTarefaMatrix, &id);
	}else{
		pthread_create(&threadMain, NULL, escolhaTarefaVector, &id);
	}

	pthread_join(threadMain, NULL);


	return 0;
}