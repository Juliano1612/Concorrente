#include <stdio.h>
#include <stdlib.h>
//#include <string.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <unistd.h>
#include "3mm.h"

#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"
#define KRST  "\x1b[0m"

#define CURRTIME() time ( &rawtime ); timeinfo = localtime ( &rawtime );

int id;

double *dTrans;
double *dE;
double *dF;
double *dG;

//vetores com os ranks que criarão comunicadores
int rank0[2] = {0,1};
int rank1[2] = {1,2};
int rank2[2] = {3,4};
int rank3[2] = {4,5};
int rank4[2] = {6,7};
int rank5[2] = {7,8};

int rank6[2] = {0,3};
int rank7[2] = {1,4};
int rank8[2] = {2,5};
int rank9[2] = {3,6};
int rank10[2] = {4,7};
int rank11[2] = {5,8};

int rank12[2] = {0,9};
int rank13[2] = {2,10};
int rank14[2] = {8,11};
int rank15[2] = {6,12};

MPI_Comm comm_world, comm0, comm1, comm2, comm3, comm4, comm5, comm6, comm7, comm8, comm9, comm10, comm11, comm12, comm13, comm14, comm15;
MPI_Group group_orig, group0, group1, group2, group3, group4, group5, group6, group7, group8, group9, group10, group11, group12, group13, group14, group15;

int world_size, world_rank;
int worldsizes[4], worldranks[4];
MPI_Comm comms[4];
MPI_Request request[4];
MPI_Status status;
int TAG = 0;

time_t rawtime;
struct tm * timeinfo;


struct tm* getTimeLog(){
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	return timeinfo;
}


void initMatrixTrans(double tamMatrix){
	int desloc = (int) tamMatrix*tamMatrix;
	for(int i = 0; i  < tamMatrix; i++){
		for(int j = 0; j < tamMatrix; j++){
			dTrans[(int)((i*tamMatrix)+j+(desloc*0)+4)] = (double) (i*j) / tamMatrix;//Matriz A
			dTrans[(int)((i*tamMatrix)+j+(desloc*1)+4)] = (double) (i*(j+1)) / tamMatrix;//Matriz B
			dTrans[(int)((i*tamMatrix)+j+(desloc*2)+4)] = (double) (i*(j+3)) / tamMatrix;//Matriz C
			dTrans[(int)((i*tamMatrix)+j+(desloc*3)+4)] = (double) (i*(j+2)) / tamMatrix;//Matriz D
		}
	}
}


// static void printMatrixLine(double* m){
// 	int tamLinha  = 0;
// 	int contaLinha = 0;
// 	//printf("\nWR%d Linha[%d]", world_rank, contaLinha);
// 	for(int i = 0; i < NI*NI; i++){
// 		if(tamLinha == NI){
// 			contaLinha++;
// 			tamLinha = 0;
// 			//printf("\nWR%d Linha[%d] %.2f", world_rank, contaLinha, m[i]);
// 			printf("\n%.2f", m[i]);
//
// 		}else{
// 			printf(" %f", m[i]);
// 		}
// 		tamLinha++;
//
// 	}
// 	printf("\n");
// }



void createGroupsAndComms(){
	MPI_Comm_group(MPI_COMM_WORLD, &group_orig);//cria o grupo original para gerar os outros

	MPI_Group_incl(group_orig, 2, rank0, &group0);//0 1
	MPI_Group_incl(group_orig, 2, rank1, &group1);//1 2
	MPI_Group_incl(group_orig, 2, rank2, &group2);//3 4
	MPI_Group_incl(group_orig, 2, rank3, &group3);//4 5
	MPI_Group_incl(group_orig, 2, rank4, &group4);//6 7
	MPI_Group_incl(group_orig, 2, rank5, &group5);//7 8

	MPI_Group_incl(group_orig, 2, rank6, &group6);//0 3
	MPI_Group_incl(group_orig, 2, rank7, &group7);//1 4
	MPI_Group_incl(group_orig, 2, rank8, &group8);//2 5
	MPI_Group_incl(group_orig, 2, rank9, &group9);//3 6
	MPI_Group_incl(group_orig, 2, rank10, &group10);//4 7
	MPI_Group_incl(group_orig, 2, rank11, &group11);//5 8

	MPI_Group_incl(group_orig, 2, rank12, &group12);//0 09
	MPI_Group_incl(group_orig, 2, rank13, &group13);//2 10
	MPI_Group_incl(group_orig, 2, rank14, &group14);//8 11
	MPI_Group_incl(group_orig, 2, rank15, &group15);//6 12

	MPI_Comm_create(MPI_COMM_WORLD, group0, &comm0);//0 1
	MPI_Comm_create(MPI_COMM_WORLD, group1, &comm1);//1 2
	MPI_Comm_create(MPI_COMM_WORLD, group2, &comm2);//3 4
	MPI_Comm_create(MPI_COMM_WORLD, group3, &comm3);//4 5
	MPI_Comm_create(MPI_COMM_WORLD, group4, &comm4);//6 7
	MPI_Comm_create(MPI_COMM_WORLD, group5, &comm5);//7 8

	MPI_Comm_create(MPI_COMM_WORLD, group6, &comm6);//0 3
	MPI_Comm_create(MPI_COMM_WORLD, group7, &comm7);//1 4
	MPI_Comm_create(MPI_COMM_WORLD, group8, &comm8);//2 5
	MPI_Comm_create(MPI_COMM_WORLD, group9, &comm9);//3 6
	MPI_Comm_create(MPI_COMM_WORLD, group10, &comm10);//4 7
	MPI_Comm_create(MPI_COMM_WORLD, group11, &comm11);//5 8

	MPI_Comm_create(MPI_COMM_WORLD, group12, &comm12);//0 09
	MPI_Comm_create(MPI_COMM_WORLD, group13, &comm13);//2 10
	MPI_Comm_create(MPI_COMM_WORLD, group14, &comm14);//8 11
	MPI_Comm_create(MPI_COMM_WORLD, group15, &comm15);//6 12

}

// void printSizeAndRank(){
// 	printf("WORLD RANK/SIZE: %d/%d \n\t NEW WORLD 1 RANK/SIZE: %d/%d \n\t NEW WORLD 2 RANK/SIZE: %d/%d \n\t NEW WORLD 3 RANK/SIZE: %d/%d \n\t NEW WORLD 4 RANK/SIZE: %d/%d \n",
// 				world_rank, world_size, worldranks[0], worldsizes[0], worldranks[1], worldsizes[1]
// 				, worldranks[2], worldsizes[2], worldranks[3], worldsizes[3]);
//
// }

void initWorldsSizeAndRanks(){

	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	/*Vetor de rank e size -> [0]UP [1]DOWN [2]LEFT [3]RIGHT */
	switch (world_rank) {
		case 0:// 3 comunicadores
			MPI_Comm_rank(comm6, &worldranks[0]);//up
			worldranks[1] = -1;//down
			MPI_Comm_rank(comm12, &worldranks[2]);//left -> processador
			MPI_Comm_rank(comm0, &worldranks[3]);//right

			MPI_Comm_size(comm6, &worldsizes[0]);//up
			worldsizes[1] = -1;//down
			MPI_Comm_size(comm12, &worldsizes[2]);//left -> processador
			MPI_Comm_size(comm0, &worldsizes[3]);//right

			comms[0] = comm6;
			comms[1] = MPI_COMM_NULL;
			comms[2] = comm12;
			comms[3] = comm0;

			break;
		case 1:
			MPI_Comm_rank(comm7, &worldranks[0]);//up
			worldranks[1] = -1;//down
			MPI_Comm_rank(comm0, &worldranks[2]);//left
			MPI_Comm_rank(comm1, &worldranks[3]);//right

			MPI_Comm_size(comm7, &worldsizes[0]);//up
			worldsizes[1] = -1;//down
			MPI_Comm_size(comm0, &worldsizes[2]);//left
			MPI_Comm_size(comm1, &worldsizes[3]);//right

			comms[0] = comm7;
			comms[1] = MPI_COMM_NULL;
			comms[2] = comm0;
			comms[3] = comm1;
			break;
		case 2:
			MPI_Comm_rank(comm8, &worldranks[0]);//up
			worldranks[1] = -1;//down
			MPI_Comm_rank(comm1, &worldranks[2]);//left
			MPI_Comm_rank(comm13, &worldranks[3]);//right  -> processador

			MPI_Comm_size(comm8, &worldsizes[0]);//up
			worldsizes[1] = -1;//down
			MPI_Comm_size(comm1, &worldsizes[2]);//left
			MPI_Comm_size(comm13, &worldsizes[3]);//right -> processador

			comms[0] = comm8;
			comms[1] = MPI_COMM_NULL;
			comms[2] = comm1;
			comms[3] = comm13;

			break;
		case 3:
			MPI_Comm_rank(comm9, &worldranks[0]);//up
			MPI_Comm_rank(comm6, &worldranks[1]);//down
			worldranks[2] = -1;//left
			MPI_Comm_rank(comm2, &worldranks[3]);//right

			MPI_Comm_size(comm9, &worldsizes[0]);//up
			MPI_Comm_size(comm6, &worldsizes[1]);//down
			worldsizes[2] = -1;//left
			MPI_Comm_size(comm2, &worldsizes[3]);//right

			comms[0] = comm9;
			comms[1] = comm6;
			comms[2] = MPI_COMM_NULL;
			comms[3] = comm2;

			break;
		case 4:
			MPI_Comm_rank(comm10, &worldranks[0]);//up
			MPI_Comm_rank(comm7, &worldranks[1]);//down
			MPI_Comm_rank(comm2, &worldranks[2]);//left
			MPI_Comm_rank(comm3, &worldranks[3]);//right

			MPI_Comm_size(comm10, &worldsizes[0]);//up
			MPI_Comm_size(comm7, &worldsizes[0]);//down
			MPI_Comm_size(comm2, &worldsizes[2]);//left
			MPI_Comm_size(comm3, &worldsizes[3]);//right

			comms[0] = comm10;
			comms[1] = comm7;
			comms[2] = comm2;
			comms[3] = comm3;

			break;
		case 5:
			MPI_Comm_rank(comm11, &worldranks[0]);//up
			MPI_Comm_rank(comm8, &worldranks[1]);//down
			MPI_Comm_rank(comm3, &worldranks[2]);//left
			worldranks[3] = -1;//right

			MPI_Comm_size(comm11, &worldsizes[0]);//up
			MPI_Comm_size(comm8, &worldsizes[1]);//down
			MPI_Comm_size(comm3, &worldsizes[2]);//left
			worldsizes[3] = -1;//right

			comms[0] = comm11;
			comms[1] = comm8;
			comms[2] = comm3;
			comms[3] = MPI_COMM_NULL;

			break;
		case 6:
			worldranks[0] = -1;//up
			MPI_Comm_rank(comm9, &worldranks[1]);//down
			MPI_Comm_rank(comm15, &worldranks[2]);//left
			MPI_Comm_rank(comm4, &worldranks[3]);//right

			worldsizes[0]= -1;//up
			MPI_Comm_size(comm9, &worldsizes[1]);//down
			MPI_Comm_size(comm15, &worldsizes[2]) ;//left
			MPI_Comm_size(comm4, &worldsizes[3]);//right

			comms[0] = MPI_COMM_NULL;
			comms[1] = comm9;
			comms[2] = comm15;
			comms[3] = comm4;

			break;
		case 7:
			worldranks[0] = -1;//up
			MPI_Comm_rank(comm10, &worldranks[1]);//down
			MPI_Comm_rank(comm4, &worldranks[2]);//left
			MPI_Comm_rank(comm5, &worldranks[3]);//right

			worldsizes[0] = -1;//up
			MPI_Comm_size(comm10, &worldsizes[1]);//down
			MPI_Comm_size(comm4, &worldsizes[2]);//left
			MPI_Comm_size(comm5, &worldsizes[3]);//right

			comms[0] = MPI_COMM_NULL;
			comms[1] = comm10;
			comms[2] = comm4;
			comms[3] = comm5;

			break;
		case 8:
			worldranks[0] = -1;//up
			MPI_Comm_rank(comm11, &worldranks[1]);//down
			MPI_Comm_rank(comm5, &worldranks[2]);//left
			MPI_Comm_rank(comm14, &worldranks[3]);//right

			worldsizes[0] = -1;//up
			MPI_Comm_size(comm11, &worldsizes[1]);//down
			MPI_Comm_size(comm5, &worldsizes[2]);//left
			MPI_Comm_size(comm14, &worldsizes[3]);//right

			comms[0] = MPI_COMM_NULL;
			comms[1] = comm11;
			comms[2] = comm5;
			comms[3] = comm14;

			break;
		case 9:
			worldranks[0] = -1;//up
			worldranks[1] = -1;//down
			worldranks[2] = -1;//left
			MPI_Comm_rank(comm12, &worldranks[3]);//right

			worldsizes[0] = -1;//up
			worldsizes[1] = -1;//down
			worldsizes[2] = -1;//left
			MPI_Comm_size(comm12, &worldsizes[3]);//right

			comms[0] = MPI_COMM_NULL;
			comms[1] = MPI_COMM_NULL;
			comms[2] = MPI_COMM_NULL;
			comms[3] = comm12;

			break;
		case 10:
			worldranks[0] = -1;//up
			worldranks[1] = -1;//down
			MPI_Comm_rank(comm13, &worldranks[2]);//left
			worldranks[3] = -1;//right

			worldsizes[0] = -1;//up
			worldsizes[1] = -1;//down
			MPI_Comm_size(comm13, &worldsizes[2]);//left
			worldsizes[3] = -1;//right

			comms[0] = MPI_COMM_NULL;
			comms[1] = MPI_COMM_NULL;
			comms[2] = comm13;
			comms[3] = MPI_COMM_NULL;

			break;
		case 11:
			worldranks[0] = -1;//up
			worldranks[1] = -1;//down
			MPI_Comm_rank(comm14, &worldranks[2]);//left
			worldranks[3] = -1;//right

			worldsizes[0] = -1;//up
			worldsizes[1] = -1;//down
			MPI_Comm_size(comm14, &worldsizes[2]);//left
			worldsizes[3] = -1;//right

			comms[0] = MPI_COMM_NULL;
			comms[1] = MPI_COMM_NULL;
			comms[2] = comm14;
			comms[3] = MPI_COMM_NULL;

			break;
		case 12:
			worldranks[0] = -1;//up
			worldranks[1] = -1;//down
			worldranks[2] = -1;//left
			MPI_Comm_rank(comm15, &worldranks[3]);//right

			worldsizes[0] = -1;//up
			worldsizes[1] = -1;//down
			worldsizes[2] = -1;//left
			MPI_Comm_size(comm15, &worldsizes[3]);//right

			comms[0] = MPI_COMM_NULL;
			comms[1] = MPI_COMM_NULL;
			comms[2] = MPI_COMM_NULL;
			comms[3] = comm15;

			break;
	}
	//printSizeAndRank();

}

void kernel_3mm_MPI_Grid(){

	int size = (int) dTrans[2];

	int desloc = size*size;
	for (int i = 0; i < size; i++){
		for (int j = 0; j < size; j++){
			dE[(i*size)+j] = 0;
			for (int k = 0; k < size; ++k){
				dE[(i*size)+j] += dTrans[(i*size)+k+4] * dTrans[desloc+(k*size)+j+4];
			}
		}
	}
	for (int i = 0; i < size; i++){
		for (int j = 0; j < size; j++){
			dF[(i*size)+j] = 0;
			for (int k = 0; k < size; ++k){
				dF[(i*size)+j] += dTrans[(2*desloc)+(i*size)+k+4] * dTrans[(3*desloc)+(k*size)+j+4];
			}
		}
	}
	for (int i = 0; i < size; i++){
		for (int j = 0; j < size; j++){
			dG[(i*size)+j] = 0;
			for (int k = 0; k < size; ++k){
				dG[(i*size)+j] += dE[(i*size)+k] * dF[(k*size)+j];
			}
		}
	}
	for(int i = 0; i < size; i++){
		for(int j = 0; j < size; j++){
			dTrans[(i*size)+j+4] = dG[(i*size)+j];
		}
	}
}


int rand_grid(){
  srand(time(0));
  int pos;
  if(world_rank == 0 || world_rank == 6){
    do{
      pos = rand() % 4;
    }while(comms[pos] == MPI_COMM_NULL || pos == 2);
    return pos;
  }
  else if(world_rank == 2 || world_rank == 8){
    do{
      pos = rand() % 3;
    }while(comms[pos] == MPI_COMM_NULL);
    return pos;
  }
  else{
    do{
      pos = rand() % 4;
    }while(comms[pos] == MPI_COMM_NULL);
    return pos;
  }
}

void loopdispatcher(int req){
	srand(time(0));
	double potencia;
	MPI_Request reqsdispatcher[req];
	for(int i = 0; i < req; i++){
		potencia = (rand() % 7)+1;
		//potencia = 7;
		potencia = pow(2,potencia);
		initMatrixTrans(potencia);
		dTrans[0] = 0.0;//enviando
		dTrans[1] = i+1;//numRequisicao
		dTrans[2] = potencia;//size
		sleep(2);
		CURRTIME()
		printf("[%d:%d:%02d] %sENVIANDO REQUISICAO %d ... MATRIZES DE TAMANHO %.0f%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KYEL , i+1, potencia, KRST);
		MPI_Send(dTrans, 65540, MPI_DOUBLE, 0, TAG, comm12);
		MPI_Irecv(dTrans, 65540, MPI_DOUBLE, 0, TAG, comm12, &reqsdispatcher[i]);
	}
	int contareq = 0, flagdisp;
	while(contareq < req){
		contareq = 0;
		for(int i = 0; i < req; i++){
			MPI_Test(&reqsdispatcher[i], &flagdisp, MPI_STATUS_IGNORE);
			if(flagdisp)
				contareq++;
		}
	}
	dTrans[0] = 2.0;
	MPI_Send(dTrans, 65540, MPI_DOUBLE, 0, TAG, comm12);
	CURRTIME()
	printf("%s[%d:%d:%02d] %sREQUISITOR FINALIZADO%s\n", KRED,timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KYEL, KRST);
}

void loopprocesser(){

	dE = (double *) malloc (sizeof(double) * 128*128);
	dF = (double *) malloc (sizeof(double) * 128*128);
	dG = (double *) malloc (sizeof(double) * 128*128);

	if(world_rank == 10 || world_rank == 11){
		while(1){
			MPI_Recv(dTrans, 65540, MPI_DOUBLE, 0, TAG, comms[2], MPI_STATUS_IGNORE);
			if(dTrans[0] == 0.0){
				CURRTIME()
				printf("[%d:%d:%02d] %sPROCESSADOR %d INICIANDO PROCESSAMENTO DE REQUISICAO %d --- MATRIZES DE TAMANHO %.0f%s\n", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KCYN, world_rank,(int) dTrans[1], dTrans[2],KRST);
				kernel_3mm_MPI_Grid();
				dTrans[0] = 1.0;
				MPI_Send(dTrans, 65540, MPI_DOUBLE, 0, TAG, comms[2]);
			}else{
				CURRTIME()
				printf("%s[%d:%d:%02d] %sPROCESSADOR %d RECEBEU MENSAGEM DE FINALIZACAO... PROCESSADOR FINALIZADO%s\n",KRED, timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KCYN, world_rank, KRST);
				break;
			}
		}
	}else{
		while(1){
			MPI_Recv(dTrans, 65540, MPI_DOUBLE, 0, TAG, comms[3], MPI_STATUS_IGNORE);
			if(dTrans[0] == 0.0){
				CURRTIME()
				printf("[%d:%d:%02d] %sPROCESSADOR %d INICIANDO PROCESSAMENTO DE REQUISICAO %d --- MATRIZES DE TAMANHO %.0f%s\n", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KCYN, world_rank,(int) dTrans[1], dTrans[2], KRST);
				kernel_3mm_MPI_Grid();
				dTrans[0] = 1.0;
				MPI_Send(dTrans, 65540, MPI_DOUBLE, 0, TAG, comms[3]);
			}else{
				CURRTIME()
				printf("%s[%d:%d:%02d] %sPROCESSADOR %d RECEBEU MENSAGEM DE FINALIZACAO... PROCESSADOR FINALIZADO%s\n", KRED, timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KCYN, world_rank, KRST);
				break;
			}
		}
	}
}

void mensagemEnvio(int lado){
	CURRTIME()
	switch (lado) {
		case 0:
			printf("[%d:%d:%02d]%s N0 %d ENVIANDO PARA NO %d%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KGRN, world_rank, world_rank+3, KRST);
			break;
		case 1:
			printf("[%d:%d:%02d]%s N0 %d ENVIANDO PARA NO %d%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KGRN, world_rank, world_rank-3, KRST);
			break;
		case 2:
			switch (world_rank) {
				case 0:
					printf("[%d:%d:%02d]%s N0 %d ENVIANDO PARA NO %d%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KGRN, world_rank, 9, KRST);
					break;
				case 6:
					printf("[%d:%d:%02d]%s N0 %d ENVIANDO PARA NO %d%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KGRN, world_rank, 12, KRST);
					break;
				default:
					printf("[%d:%d:%02d]%s N0 %d ENVIANDO PARA NO %d%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KGRN, world_rank, world_rank-1, KRST);
			}
			break;
		case 3:
			switch (world_rank) {
				case 2:
					printf("[%d:%d:%02d]%s N0 %d ENVIANDO PARA NO %d%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KGRN, world_rank, 10, KRST);
					break;
				case 8:
					printf("[%d:%d:%02d]%s N0 %d ENVIANDO PARA NO %d%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KGRN, world_rank, 11, KRST);
					break;
				default:
					printf("[%d:%d:%02d]%s N0 %d ENVIANDO PARA NO %d%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KGRN, world_rank, world_rank+1, KRST);
			}

	}
}

void mensagemRecebimento(int lado, double req){
	CURRTIME()
	switch (lado) {
		case 0:
			printf("[%d:%d:%02d] %sN0 %d RECEBEU DE NO %d A REQUISICAO %.0f%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KGRN, world_rank, world_rank+3, req, KRST);
			break;
		case 1:
			printf("[%d:%d:%02d] %sN0 %d RECEBEU DE NO %d A REQUISICAO %.0f%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KGRN, world_rank, world_rank-3, req, KRST);
			break;
		case 2:
			switch (world_rank) {
				case 0:
					printf("[%d:%d:%02d] %sN0 %d RECEBEU DE NO %d A REQUISICAO %.0f%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KGRN, world_rank, 9, req, KRST);
					break;
				case 6:
					printf("[%d:%d:%02d] %sN0 %d RECEBEU DE NO %d A REQUISICAO %.0f%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KGRN, world_rank, 12, req, KRST);
					break;
				default:
					printf("[%d:%d:%02d] %sN0 %d RECEBEU DE NO %d A REQUISICAO %.0f%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KGRN, world_rank, world_rank-1, req, KRST);
			}
			break;
		case 3:
			switch (world_rank) {
				case 2:
					printf("[%d:%d:%02d] %sN0 %d RECEBEU DE NO %d A REQUISICAO %.0f%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KGRN, world_rank, 10, req, KRST);
					break;
				case 8:
					printf("[%d:%d:%02d] %sN0 %d RECEBEU DE NO %d A REQUISICAO %.0f%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KGRN, world_rank, 11, req, KRST);
					break;
				default:
					printf("[%d:%d:%02d] %sN0 %d RECEBEU DE NO %d A REQUISICAO %.0f%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KGRN, world_rank, world_rank+1, req, KRST);
			}

	}
}

void mensagemEnvioProcessar(){
	CURRTIME()
	switch (world_rank) {
		case 2:
			printf("[%d:%d:%02d] %sNO %d ENVIANDO PARA PROCESSADOR %d%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KGRN, world_rank, 10, KRST);
			break;
		case 6:
			printf("[%d:%d:%02d] %sNO %d ENVIANDO PARA PROCESSADOR %d%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KGRN, world_rank, 12, KRST);
			break;
		case 8:
			printf("[%d:%d:%02d] %sNO %d ENVIANDO PARA PROCESSADOR %d%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KGRN, world_rank, 11, KRST);
	}
}

void mensagemProcessadorOcupado(int lado){
	CURRTIME()
	switch (world_rank) {
		case 2:
			printf("[%d:%d:%02d] %sPROCESSADOR %d OCUPADO%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KGRN, 10, KRST);
			break;
		case 6:
			printf("[%d:%d:%02d] %sPROCESSADOR %d OCUPADO%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KGRN, 12, KRST);
			break;
		case 8:
			printf("[%d:%d:%02d] %sPROCESSADOR %d OCUPADO%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KGRN, 11, KRST);
	}
	mensagemEnvio(lado);
}

void loopgrid(){
	double test[4][65540];
	int TAG = 0, myproc = 1;

	int indices[4], num_completed, flagsreqs[4];

	for(int i = 0; i < 4; i++){
		request[i] = MPI_REQUEST_NULL;
	}

	while(1){
		sleep(1);
		for(int i = 0; i < 4; i++){
			if(comms[i] != MPI_COMM_NULL){
				if(request[i] == MPI_REQUEST_NULL){
					if(worldranks[i] == 0){
						MPI_Irecv(&test[i], 65540, MPI_DOUBLE, 1, TAG, comms[i], &request[i]);
					}else{
						MPI_Irecv(&test[i], 65540, MPI_DOUBLE, 0, TAG, comms[i], &request[i]);
					}
				}
			}
		}
		num_completed = 0;
		while(num_completed == 0){
			for(int i = 0; i < 4; i++){
				if(comms[i] != MPI_COMM_NULL){
					MPI_Test(&request[i], &flagsreqs[num_completed], &status);
					if(flagsreqs[num_completed]){
						indices[num_completed] = i;
						num_completed++;
					}
				}
			}
		}
		for(int i = 0; i < num_completed; i++){
			double flagCaminho = test[indices[i]][0];
			int pos = indices[i];
			if(flagCaminho == 0.0){
				mensagemRecebimento(indices[i], test[indices[i]][1]);
				mensagemEnvioProcessar();
				if(world_rank == 2 || world_rank == 6 || world_rank == 8){
					//MPI_Request reqProc = MPI_REQUEST_NULL;
					int flagproc = 0;
					switch(world_rank){
						case 2:
						case 8:
							if(indices[i] == 3){//a mensagem veio do processador
								myproc = 1;
								flagproc = 1;
								CURRTIME()
								if(world_rank == 2)
									printf("[%d:%d:%02d] %sPROCESSADOR 10 COMPLETOU A REQUISICAO %d%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KCYN, (int) test[indices[i]][1], KRST);
								else
									printf("[%d:%d:%02d] %sPROCESSADOR 11 COMPLETOU A REQUISICAO %d%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KCYN, (int) test[indices[i]][1], KRST);
							}
							if(myproc && !flagproc){//se o processador está disponível
								MPI_Send(&test[indices[i]], 65540, MPI_DOUBLE, 1, TAG, comms[3]);
								myproc = 0;
							}else if (!flagproc){
								while(pos == indices[i])
									pos = rand_grid();
								mensagemProcessadorOcupado(pos);
								if(worldranks[pos] == 0){
									MPI_Send(&test[indices[i]], 65540, MPI_DOUBLE, 1, TAG, comms[pos]);
								}else{
									MPI_Send(&test[indices[i]], 65540, MPI_DOUBLE, 0, TAG, comms[pos]);
								}
							}
							request[indices[i]] = MPI_REQUEST_NULL;
							break;
						case 6:
							if(indices[i] == 2){//a mensagem veio do processador
								myproc = 1;
								flagproc = 1;
								CURRTIME()
								printf("[%d:%d:%02d] %sPROCESSADOR 12 COMPLETOU A REQUISICAO %d%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KCYN,  (int) test[indices[i]][1], KRST);
							}else if(myproc && !flagproc){//se o processador está disponível
								MPI_Send(&test[indices[i]], 65540, MPI_DOUBLE, 1, TAG, comms[2]);
								myproc = 0;
							}else if(!flagproc){
								while(pos == indices[i])
									pos = rand_grid();
								mensagemProcessadorOcupado(pos);
								if(worldranks[pos] == 0){
									MPI_Send(&test[indices[i]], 65540, MPI_DOUBLE, 1, TAG, comms[pos]);
								}else{
									MPI_Send(&test[indices[i]], 65540, MPI_DOUBLE, 0, TAG, comms[pos]);
								}
							}
							request[indices[i]] = MPI_REQUEST_NULL;
							break;
					}
				}else{
					while(pos == indices[i])
						pos = rand_grid();
					mensagemEnvio(pos);
					if(worldranks[pos] == 0){
						MPI_Send(&test[indices[i]], 65540, MPI_DOUBLE, 1, TAG, comms[pos]);
					}else{
						MPI_Send(&test[indices[i]], 65540, MPI_DOUBLE, 0, TAG, comms[pos]);
					}
				}
			request[indices[i]] = MPI_REQUEST_NULL;
		}else if(flagCaminho == 1.0){//mensagem está voltando
			CURRTIME()
			printf("[%d:%d:%02d] %sNO %d RETORNANDO A RESPOSTA DE %.0f PARA O REQUISITOR%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KWHT, world_rank, test[indices[i]][1], KRST);
			if(world_rank == 2 || world_rank == 6 || world_rank == 8){
				int flagproc = 0;
				switch(world_rank){
					case 2:
					case 8:
						if(indices[i] == 3){//a mensagem veio do processador
							myproc = 1;
							flagproc = 1;
							CURRTIME()
							if(world_rank == 2)
								printf("[%d:%d:%02d] %sPROCESSADOR 10 COMPLETOU A REQUISICAO %d%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KCYN, (int) test[indices[i]][1], KRST);
							else
								printf("[%d:%d:%02d] %sPROCESSADOR 11 COMPLETOU A REQUISICAO %d%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KCYN, (int) test[indices[i]][1], KRST);
						}
						if(world_rank == 2){
							mensagemEnvio(2);
							MPI_Send(&test[indices[i]], 65540, MPI_DOUBLE, 0, TAG, comms[2]);
						}else{
							mensagemEnvio(1);
							MPI_Send(&test[indices[i]], 65540, MPI_DOUBLE, 0, TAG, comms[1]);
						}
						request[indices[i]] = MPI_REQUEST_NULL;
						break;
					case 6:
						if(indices[i] == 2){//a mensagem veio do processador
							myproc = 1;
							flagproc = 1;
							CURRTIME()
							printf("[%d:%d:%02d] %sPROCESSADOR 12 COMPLETOU A REQUISICAO %d%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KCYN, (int) test[indices[i]][1], KRST);
						}
						mensagemEnvio(1);
						MPI_Send(&test[indices[i]], 65540, MPI_DOUBLE, 0, TAG, comms[1]);
						request[indices[i]] = MPI_REQUEST_NULL;
						break;
				}
			}else if(world_rank == 0 || world_rank == 1){
				mensagemEnvio(2);
				if(world_rank == 0){
					MPI_Send(&test[indices[i]], 65540, MPI_DOUBLE, 1, TAG, comms[2]);
				}else{
					MPI_Send(&test[indices[i]], 65540, MPI_DOUBLE, 0, TAG, comms[2]);
				}
			}else{
				mensagemEnvio(1);
				MPI_Send(&test[indices[i]], 65540, MPI_DOUBLE, 0, TAG, comms[1]);
			}
		}else{
			CURRTIME()
			printf("[%d:%d:%02d] %sNO %d RECEBEU MENSAGEM DE FINALIZACAO... REPASSANDO PARA VIZINHOS%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KGRN, world_rank, KRST);
			switch (world_rank) {
				case 0:
				case 1:
				case 2:
					MPI_Send(&test[indices[i]], 65540, MPI_DOUBLE, 1, TAG, comms[0]);
					MPI_Send(&test[indices[i]], 65540, MPI_DOUBLE, 1, TAG, comms[3]);
					break;
				case 3:
				case 4:
				case 5:
					MPI_Send(&test[indices[i]], 65540, MPI_DOUBLE, 1, TAG, comms[0]);
					break;
				case 6:
					MPI_Send(&test[indices[i]], 65540, MPI_DOUBLE, 1, TAG, comms[2]);
					break;
				case 8:
					MPI_Send(&test[indices[i]], 65540, MPI_DOUBLE, 1, TAG, comms[3]);
					break;
			}
			CURRTIME()
			printf("%s[%d:%d:%02d] %sNO %d FINALIZADO%s\n",KRED, timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KGRN, world_rank, KRST);

			return;
		}
		request[indices[i]] = MPI_REQUEST_NULL;
		}
	}

}

	int main(int argc, char** argv)
	{

	MPI_Init(NULL, NULL);
	dTrans = (double *) malloc (sizeof(double) * (128*128*4+4));

	createGroupsAndComms();
	initWorldsSizeAndRanks();

	if(world_rank == 9){
		loopdispatcher(atoi(argv[1]));
	}else if(world_rank >= 0 && world_rank <= 8){
		loopgrid();
	}else{
		loopprocesser();
		free(dE);
		free(dF);
		free(dG);
	}
	free(dTrans);
	MPI_Finalize();
	if(world_rank == 0){
		CURRTIME()
		printf("[%d:%d:%02d] %sExecucao concluida%s\n",timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, KMAG, KRST);
	}
	return 0;
}
