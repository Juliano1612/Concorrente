#include <stdio.h>
#include <mpi.h>




int main (){

	int ranks1[8] = {0,1,2,3,4,5,6,7};
	int ranks2[8] = {8,9,10,11,12,13,14,15};

	MPI_Init(NULL, NULL);


	MPI_Comm comm_world, new_comm1, new_comm2;
	MPI_Group group_orig, group_world1, group_world2;
	int ierr;

	//comm_world = MPI_COMM_WORLD;
	MPI_Comm_group(MPI_COMM_WORLD, &group_orig);

	//MPI_Comm_group(group_orig, &group_world1);
	//MPI_Comm_group(group_orig, &group_world2);

	MPI_Group_incl(group_orig, 8, ranks1, &group_world1);
	MPI_Group_incl(group_orig, 8, ranks2, &group_world2);

	MPI_Comm_create(MPI_COMM_WORLD, group_world1, &new_comm1);
	MPI_Comm_create(MPI_COMM_WORLD, group_world2, &new_comm2);

	int world_rank = 0, world_size = 0;
	int world_rank1 = 0, world_size1 = 0;
	int world_rank2 = 0, world_size2 = 0;

	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);


	if(world_rank < 8 ){

		MPI_Comm_rank(new_comm1, &world_rank1);
		MPI_Comm_size(new_comm1, &world_size1);

		printf("WORLD RANK/SIZE: %d/%d \t NEW WORLD 1 RANK/SIZE: %d/%d \n",
		world_rank, world_size, world_rank1, world_size1);

	}else{

		MPI_Comm_rank(new_comm2, &world_rank2);
		MPI_Comm_size(new_comm2, &world_size2);

		printf("WORLD RANK/SIZE: %d/%d \t NEW WORLD 2 RANK/SIZE: %d/%d\n",
			world_rank, world_size, world_rank2, world_size2);




	}

	
	
	

	/*// Get the rank and size in the original communicator
	int world_rank, world_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int color = world_rank / 4; // Determine color based on row

	// Split the communicator based on the color and use the
	// original rank for ordering
	MPI_Comm row_comm;
	MPI_Comm_split(MPI_COMM_WORLD, color, world_rank, &row_comm);

	int row_rank, row_size;
	MPI_Comm_rank(row_comm, &row_rank);
	MPI_Comm_size(row_comm, &row_size);

	printf("WORLD RANK/SIZE: %d/%d \t ROW RANK/SIZE: %d/%d\n",
		world_rank, world_size, row_rank, row_size);

	MPI_Comm_free(&row_comm);
	*/
	MPI_Finalize();



}