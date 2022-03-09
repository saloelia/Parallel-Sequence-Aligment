#include "utility.h"
#include "calculation.h"

void get_scores(double* weights,char* seq1,char** seqs,int* sizes,int num_of_sequences,Mutant* scores,int rank, MPI_Status status){
	MPI_Datatype MutantMPIType;
	create_MPI_mutant_type(&MutantMPIType);

	for (int i = 0; i < num_of_sequences; ++i) {
		int allOffsets = sizes[0] - sizes[i+1] + 1;							//offsets = n-m+1, n is length of seq1, m is length of seq i
		int start = rank * (allOffsets / 2);								//starting index depending on each process
		int end = start + (allOffsets / 2);									//ending index depending on each process

		if(allOffsets%2!=0)													//if number of offsets is odd, rank 1 will do one more offset
			if(rank==1)
				end++;

		Mutant myBestMutant;
		useOpenMp_CUDA(start,end,seq1,seqs[i],weights,&myBestMutant);
		
		if(rank==0){
			Mutant his_best;
			MPI_Recv(&his_best, 1, MutantMPIType, 1, 0, MPI_COMM_WORLD, &status);
			if(myBestMutant.score > his_best.score){
				scores[i] = myBestMutant;
			}else{
				scores[i] = his_best;
			}

		}else{
			MPI_Send(&myBestMutant,1,MutantMPIType,0,0,MPI_COMM_WORLD);
		}
	}
}

void useOpenMp_CUDA(int start,int end,char* seq1,char* seq2,double* weights,Mutant* bestMutant){
	Mutant best_of_each_thread[4];
	
	for (int i = 0; i < 4; ++i) {
		best_of_each_thread[i].score = -5000;
	}

	omp_set_num_threads(4);
	#pragma omp parallel
	{
	Mutant* all_mutants_in_offset = NULL;
	int all_mutants_size = 0;
	prepareMutants(&all_mutants_in_offset,strlen(seq2)-1,&all_mutants_size);
	
	#pragma omp for
	for (int i = start; i < end; ++i) {
		Mutant best_of_offset;
		CUDA(all_mutants_in_offset,all_mutants_size,seq1+i,seq2,strlen(seq2),weights,&best_of_offset);
		best_of_offset.offset = i;
		if(best_of_offset.score > best_of_each_thread[omp_get_thread_num()].score){
			best_of_each_thread[omp_get_thread_num()] = best_of_offset;
		}

	}
	free(all_mutants_in_offset);
	}
	max_in_thread_array(best_of_each_thread,bestMutant);
	
	}


void max_in_thread_array(Mutant* mutants,Mutant* best){
	Mutant max = mutants[0];
	for (int i = 1; i < 4; ++i) {
		if(mutants[i].score>max.score){
			max = mutants[i];
		}
	}
	best->score = max.score;
	best->mut_end = max.mut_end;
	best->mut_start = max.mut_start;
	best->offset = max.offset;
}

