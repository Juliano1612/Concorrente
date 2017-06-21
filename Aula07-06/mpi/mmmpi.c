#include <mpi.h>
#include <stdio.h>

int **m1;
int **m2;
int **m3;


void initializeMatrix(int tam){
	m1 = (int**) malloc(sizeof(int*) * tam);
	m2 = (int**) malloc(sizeof(int*) * tam);
  	m3 = (int**) malloc(sizeof(int*) * tam);

  	for(i = 0; i < tam; i++) {
    	m1[i] = (int*) malloc(sizeof(int*) * tam);
    	m2[i] = (int*) malloc(sizeof(int*) * tam);
    	m3[i] = (int*) malloc(sizeof(int*) * tam);
    	for(j = 0; j < tam; j++){
    		m1[i][j] = (i*j)-(i+j);
    		m2[i][j] = (i*j)+(i-j);
	    	m3[i][j] = 0;
    	}
	}
}




int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);




    // Finalize the MPI environment.
    MPI_Finalize();
}
