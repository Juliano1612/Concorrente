#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <omp.h>
#include <time.h>

#define TIME() struct timespec start, finish; double elapsed; clock_gettime(CLOCK_MONOTONIC, &start);
#define ENDTIME() clock_gettime(CLOCK_MONOTONIC, &finish); elapsed = (finish.tv_sec - start.tv_sec); elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0; printf("%f\n", elapsed);

int tam = 0, threads = 0, tamParte = 0, fator = 2;
pthread_barrier_t barreira;

int **m;

void printMatrix() {
  int i,j;
  for(i = 0; i < tam; i++) {
    for(j = 0; j < tam; j++) {
      printf("%d\t", m[i][j]);
    }
    printf("\n");
  }
}


void multiplicaSequencial() {
  int i, j;

  /*
  A Matriz deve ser varrida linha a linha realizando o relaxamento
  As bordas são ignoradas e será tomado que a posição inicial de processamento
  será a 1,1
  A inicialização da variável j serve para alternar a célula de início do
  processamento, já que em uma linha é começado com Vermelho e na próxima com processamento
  */
  //relaxamento das peças vermelhas
  TIME()
  for(i = 1; i < tam-1; i++) {
    for(j = 2 - (i%2); j < tam-1; j+=2) {
      m[i][j] = ((m[i-1][j] + m[i+1][j] + m[i][j+1] + m[i][j-1])/4)*fator;
    }
  }
  //relaxamento das peças pretas
  for(i = 1; i < tam-1; i++) {
    for(j = (1 + (i%2)); j < tam-1; j+=2) {
      m[i][j] = ((m[i-1][j] + m[i+1][j] + m[i][j+1] + m[i][j-1])/4)*fator;
    }
  }
  ENDTIME()
}

void multiplicaOPENMP() {
  /*
  A Matriz deve ser varrida linha a linha realizando o relaxamento
  As bordas são ignoradas e será tomado que a posição inicial de processamento
  será a 1,1
  A inicialização da variável j serve para alternar a célula de início do
  processamento, já que em uma linha é começado com Vermelho e na próxima com processamento
  */
  //relaxamento das peças vermelhas
  TIME()
  #pragma omp parallel num_threads(threads)
  {

    #pragma omp for simd
    for(int i = 1; i < tam-1; i++) {
      for(int j = 2 - (i%2); j < tam-1; j+=2) {
        m[i][j] = ((m[i-1][j] + m[i+1][j] + m[i][j+1] + m[i][j-1])/4)*fator;
      }
    }
    //relaxamento das peças pretas
    #pragma omp for simd
    for(int i = 1; i < tam-1; i++) {
      for(int j = 1 + (i%2); j < tam-1; j+=2) {
        m[i][j] = ((m[i-1][j] + m[i+1][j] + m[i][j+1] + m[i][j-1])/4)*fator;
      }
    }
  }
  ENDTIME()
}

void *multiplicaParalelo(void *id) {
  int begin = (*(int *) id) * tamParte;
  int end = begin+tamParte;
  /*
  A Matriz deve ser varrida linha a linha realizando o relaxamento
  As bordas são ignoradas e será tomado que a posição inicial de processamento
  será a 1,1
  A inicialização da variável j serve para alternar a célula de início do
  processamento, já que em uma linha é começado com Vermelho e na próxima com processamento
  */
  //relaxamento das peças vermelhas
    for(int i = begin+1; i < end-1; i++) {
      for(int j = 2 - (i%2); j < tam-1; j+=2) {
        m[i][j] = ((m[i-1][j] + m[i+1][j] + m[i][j+1] + m[i][j-1])/4)*fator;
      }
    }
    pthread_barrier_wait(&barreira);
    //relaxamento das peças pretas
    for(int i = begin+1; i < end-1; i++) {
      for(int j = 1 + (i%2); j < tam-1; j+=2) {
        m[i][j] = ((m[i-1][j] + m[i+1][j] + m[i][j+1] + m[i][j-1])/4)*fator;
      }
    }
}


int main(int argc, char** argv){

  srand(time(0));
  int op; //0 openmp 1 pthreads, por default vai sequencial
  int i, j, k;
  tam = atoi(argv[1]);
  threads = atoi(argv[2]);
  op = atoi(argv[3]);
  tamParte = tam/threads;

  m = (int**) malloc(sizeof(int*) * tam);

  for(i = 0; i < tam; i++) {
    m[i] = (int*) malloc(sizeof(int) * tam);
    for(j = 0; j < tam; j++){
      m[i][j] = rand()%10;
    }
  }

  pthread_barrier_init(&barreira, NULL, threads);
  pthread_t arrayThreads[threads];
  int vet[threads];


  switch(op) {
    case 0: {
        multiplicaOPENMP();
      } break;

    case 1: {
        TIME()
        for(i = 0; i < threads; i++){
          vet[i] = i;
          pthread_create(&arrayThreads[i], NULL, multiplicaParalelo, &vet[i]);
        }
        for(j = 0; j < threads; j++){
          pthread_join(arrayThreads[j], NULL);
        }
        ENDTIME()
      }
      break;

    default: {
        multiplicaSequencial();
      }
      break;
  }
  //printMatrix();
  return 0;


}
