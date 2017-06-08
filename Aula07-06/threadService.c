#include <stdio.h>
#include <stdlib.h>
#include <time.h>
//#include <pthread.h>



struct node
{
    int servico;
    struct node *next;
};

struct node *front = NULL;
struct node *rear = NULL;


void enqueue(int idServico){
    struct node *nptr = malloc(sizeof(struct node));
    nptr->servico = idServico;
    nptr->next = NULL;
    if (rear == NULL)
    {
        front = nptr;
        rear = nptr;
    }
    else
    {
        rear->next = nptr;
        rear = rear->next;
    }
}

int dequeue(){
	int serv;
	if(front == NULL){
		return -1;
	}else{
        struct node *temp;
        temp = front;
        front = front->next;
        //printf("\n\n%d deleted", temp->servico);
        serv = temp->servico;
        free(temp);
        return serv;
    }
}




int tamanhoMatriz = 32, numThreads = 10;
pthread_mutex_t lock;


void mm(){

	srand(time(0));
	int **A, **B, **C;

	A = (int**) malloc(sizeof(int*) * tamanhoMatriz);
	B = (int**) malloc(sizeof(int*) * tamanhoMatriz);
	C = (int**) malloc(sizeof(int*) * tamanhoMatriz);
	for (int i = 0; i < tamanhoMatriz; i++) {
		A[i] = (int*) malloc(sizeof(int) * tamanhoMatriz);
		B[i] = (int*) malloc(sizeof(int) * tamanhoMatriz);
		C[i] = (int*) malloc(sizeof(int) * tamanhoMatriz);
		for (int j = 0; j < tamanhoMatriz; j++) {
			A[i][j] = rand() % 2;
			B[i][j] = rand() % 2;
			C[i][j] = 0;
		}
	}

  /* C := A*B */
  for (int i = 0; i < tamanhoMatriz; i++){
    	for (int j = 0; j < tamanhoMatriz; j++){
			C[i][j] = 0;
			for (int k = 0; k < tamanhoMatriz; ++k){
			  C[i][j] += A[i][k] * B[k][j];
			}
    	}
	}
	printf("\nAcabou multiplicação!");
}


void vectorSearch(){
	srand(time(0));
	int *V;
	int buscado = rand() % 10;

	V = (int*) malloc(sizeof(int) *tamanhoMatriz);
	for(int i = 0; i < tamanhoMatriz; i++){
		V[i] = rand() % 10;
	}

	for(int i = 0; i < tamanhoMatriz; i++){
		if(V[i] == buscado){
			printf("\nO elemento %d existe no vetor!", buscado);
			break;
		}
	}

}

void *executaServicos(){

	pthread_mutex_lock(&lock);
	int idService = dequeue();
	pthread_mutex_unlock(&lock);
	if(dequeue != -1){
		(*func[idService])();
	}

}

void *atribuiServicos(){
	srand(time(0));
	int idServico, count = 10;

	pthread_t threads[numThreads];

	for(int i = 0; i < numThreads; i++){
		pthread_create(&threads[i], NULL, executaServicos, NULL);	
	}

	while(count > 0){
		enqueue(rand() % 2);

		count--;
	}

	for(int j = 0; j < numThreads; j++){
		pthread_join(threads[j], NULL);
	}
}

void (*func[])() = {mm, vectorSearch};

int main(){

	
	pthread_t mainThread;
	pthread_create(&mainThread, NULL, atribuiServicos, NULL);
	pthread_join(mainThread, NULL);


	return 0;
}