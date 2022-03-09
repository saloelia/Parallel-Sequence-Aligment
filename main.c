#include "utility.h"
#include "calculation.h"

int main(int argc, char *argv[]) {
    float t1 = omp_get_wtime();
    MPI_Status status;
    MPI_Init(&argc, &argv);
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (size != 2) {
       printf("Run the example with two processes only\n");
       MPI_Abort(MPI_COMM_WORLD, __LINE__);
    }

	char* seq1 = NULL;																	//seq1
	double* weights = NULL;																//weights array
	int num_of_sequences = 0;															//number of sequences
	int* sizes = NULL;																	//all sizes of all sequences (seq1, seq2,..)
	char** seqs = NULL;																	//sequences
	Mutant* scores = NULL;																//scores array

	if(rank==0){

	read_file(&weights,&seq1,&seqs,&sizes,&num_of_sequences);						//read file and get arguments
    	scores = (Mutant*)malloc(num_of_sequences*sizeof(Mutant));						//allocate mutant scores array
    	MPI_Send(&num_of_sequences, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);					//send num of seqs to slave
    	MPI_Send(weights, 4, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);							//send Weights
    	MPI_Send(sizes, num_of_sequences+1, MPI_INT, 1, 0, MPI_COMM_WORLD);				//send sizes
    	MPI_Send(seq1,sizes[0],MPI_CHAR,1,0,MPI_COMM_WORLD);							//send seq1


    	for (int i = 0; i < num_of_sequences; ++i) {									//send all sequences
    		MPI_Send(seqs[i],sizes[i+1],MPI_CHAR,1,0,MPI_COMM_WORLD);
		}


	}else{

	MPI_Recv(&num_of_sequences, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);			//recieve num of seqs from master
    	weights = (double*)malloc(4*sizeof(double));									//allocate weigths
    	MPI_Recv(weights, 4, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);				//recieve weights
    	sizes = (int*)malloc((num_of_sequences+1)*sizeof(int));							//allocate sizes
    	MPI_Recv(sizes, num_of_sequences+1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);				//recieve sizes
    	seq1 = (char*)malloc(sizes[0]*sizeof(char));
    	MPI_Recv(seq1, sizes[0], MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
	
    	seqs = (char**)malloc(num_of_sequences*sizeof(char*));
    	for (int i = 0; i < num_of_sequences; ++i) {
       		seqs[i] = (char*)malloc((sizes[i+1])*sizeof(char));
       		MPI_Recv(seqs[i], sizes[i+1], MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
       		
    	}

    }

	get_scores(weights,seq1,seqs,sizes,num_of_sequences,scores,rank,status);
	
	if(rank==0){
		for (int i = 0; i < num_of_sequences; ++i) {
				printf("For Sequence Number %d, the best score is %lf , offset=%d ,MS(%d,%d)\n",i+1,scores[i].score,scores[i].offset,scores[i].mut_start,scores[i].mut_end);
			}
	float t2 = omp_get_wtime();
	printf("Time of calculation %f minutes\n",(t2-t1)/60);		
	free(scores);
	}

    free(seq1);
    free(seqs);
    free(sizes);
    free(weights);

    MPI_Finalize();

    return 0;
}



