#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <omp.h>


int tam = 0, threads = 0, tamParte = 0;

int **m1;
int **m2;
int **m3;


//Precisa código sequencial e com pthreads
//Para o cálculo do speed up


void printMatrix() {
  int i,j;
  for(i = 0; i < tam; i++) {
    for(j = 0; j < tam; j++) {
      printf("%d\t", m3[i][j]);
    }
    printf("\n");
  }
}

void multiplicaSequencial() {
  int i, j ,k;
  for(i = 0; i < tam; i++) {
    for(j = 0; j < tam; j++) {
      for(k = 0; k < tam; k++) {
        m3[i][j] += m1[i][k] * m2[k][j];
      }
    }
  }
}

void multiplicaOPENMP() {
  #pragma omp parallel for num_threads(4)
  for(int i = 0; i < tam; i++) {
    for(int j = 0; j < tam; j++) {
      for(int k = 0; k < tam; k++) {
        m3[i][j] += m1[i][k] * m2[k][j];
      }
    }
  }
}


void *multiplicaParalelo(void *id) {
  int begin = (*(int *) id) * tamParte;
  int end = begin+tamParte;
  int i, j, k;

  for(i = begin; i < end; i++) {
    for(j = 0; j < tam; j++) {
      for(k = 0; k < tam; k++) {
        m3[i][j] += m1[i][k] * m2[k][j];
      }
    }
  }


}

void inicializaMatrix() {
  int i, j;
  for(i = 0; i < tam; i++) {
    for(j = 0; j < tam; j++) {
      m3[i][j] = 0;
    }
  }
}



int main(int argc, char** argv) {
  srand(time(0));
  int op = 3 ; //0 openmp 1 pthreads
  int i, j, k;
  tam = atoi(argv[1]);
  threads = atoi(argv[2]);
  tamParte = tam/threads;

  m1 = (int**) malloc(sizeof(int*) * tam);
  m2 = (int**) malloc(sizeof(int*) * tam);
  m3 = (int**) malloc(sizeof(int*) * tam);

  for(i = 0; i < tam; i++) {
    m1[i] = (int*) malloc(sizeof(int*) * tam);
    m2[i] = (int*) malloc(sizeof(int*) * tam);
    m3[i] = (int*) malloc(sizeof(int*) * tam);
    for(j = 0; j < tam; j++){
      m1[i][j] = rand()%10;
      m2[i][j] = rand()%10;
      m3[i][j] = 0;
    }
  }

  pthread_t arrayThreads[threads];
  int vet[threads];

  switch(op) {
    case 0: {
        multiplicaOPENMP();
      } break;

    case 1: {
        inicializaMatrix();
        for(i = 0; i < threads; i++){
          vet[i] = i;
          pthread_create(&arrayThreads[i], NULL, multiplicaParalelo, &vet[i]);
        }
        for(j = 0; j < threads; j++){
          pthread_join(arrayThreads[j], NULL);
        }
      } break;

      default: {
        multiplicaSequencial();
      } break;
  }
  return 0;
}
