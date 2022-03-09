#include "utility.h"


void create_MPI_mutant_type(MPI_Datatype *mutantType){
	Mutant mutant;
	MPI_Datatype type[4] = { MPI_INT , MPI_DOUBLE , MPI_INT , MPI_INT };
    int blocklen[4] = { 1, 1, 1, 1};
    MPI_Aint disp[4];
    //disp[0] = (char *) &mutant.seq -	 (char *) &mutant;
    disp[0] = (char *) &mutant.offset -	 (char *) &mutant;
    disp[1] = (char *) &mutant.score - (char *) &mutant;
    disp[2] = (char *) &mutant.mut_start - (char *) &mutant;
    disp[3] = (char *) &mutant.mut_end - (char *) &mutant;
    MPI_Type_create_struct(4, blocklen, disp, type, mutantType);
    MPI_Type_commit(mutantType);
}

void read_file(double** weights,char** seq1,char*** seqs,int** sizes,int* num_of_sequences){			//read file method
	FILE *file = fopen(FILE_INPUT, "r");																//open file

	if (file==NULL) {
		fprintf(stderr, "OPEN FILE ERROR\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	*(weights) = (double*)malloc(4*sizeof(double));														//allocate weights array
	if(!*weights){
		fprintf(stderr, "Weights Array Allocation Error\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	if (fscanf(file, "%lf%lf%lf%lf", (*weights), (*weights)+1, (*weights)+2, (*weights)+3 ) != 4){		//scan weights to array
		fprintf(stderr, "Weights scan error\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	char c;
	while ((c = fgetc(file)) != '\n');																	//read until \n

	*(sizes) = (int*)malloc(1*sizeof(int));																//allocate sizes array (N+1)
	if(!*sizes){
		fprintf(stderr, "Sizes Array Allocation Error\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	readSequence(file,seq1,*(sizes));																	//read seq1 with readSequence function and add size of seq1 to sizes array[0]


	if (fscanf(file, "%d", num_of_sequences ) != 1){													//scan number of sequences
		fprintf(stderr, "Number of sequences scan error\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	while ((c = fgetc(file)) != '\n');																	//read until \n

	*(seqs) = (char**)malloc(*(num_of_sequences)*sizeof(char*));										//alocate seqs of seqs array
	if(!*seqs){
		fprintf(stderr, "Sizes Array Allocation Error\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	*(sizes) = (int*)realloc(*(sizes),*(num_of_sequences)*sizeof(int));									//realloc on sizes array (add num of seqs)
	if(!*sizes){
		fprintf(stderr, "Sizes Array Allocation Error\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	for (int i = 0; i < *(num_of_sequences); ++i) {														//fill seqs and sizes arrays
		readSequence(file, (*seqs)+i, (*sizes)+i+1);
	}

}

void readSequence(FILE *file,char** seq,int* size_of_seq_i){											//read sequence function
	char c;
	int size = 0;
	*(seq) = NULL;
	while ( ((c = fgetc(file)) != '\n') || (feof(file))){
		size++;
		*(seq) = (char*)realloc(*(seq),size*sizeof(char));
		(*seq)[size-1] = c;
	}
	*(size_of_seq_i) = size;

}

void prepareMutants(Mutant** mutants,int n,int* num_of_mutants_pointer){
	int num_of_mutants = (((n)*(n-1))/ 2 );											//num of mutants in seq is nC2= n*(n-1)/2																	
	//int offsets = end-start;
	
	*(num_of_mutants_pointer) = num_of_mutants;
	//printf("here ok\n");
	*mutants = (Mutant*)malloc(num_of_mutants*sizeof(Mutant));									//allocate mutants array
	if(!(*mutants)){
		fprintf(stderr, "Mutants Array Allocation Error\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	

	fillMutants(*mutants,n);											  		//fill mutant array with mutants

}

void fillMutants(Mutant* mutant_array,int n){
	int index = 0;
	//for (int i = start; i < end; ++i) {
	int counter1 = 0;
	int counter2 = 1;
	while(counter1!=n){
		for (int var = counter2; var < n; ++var) {
				//mutant_array[index].offset = i;
				mutant_array[index].mut_start = counter1;
				mutant_array[index].mut_end = var;
				index++;
			}
			counter1++;
			counter2++;
		}
	//}
}


