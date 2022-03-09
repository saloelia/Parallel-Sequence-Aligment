/*
 * calculation.h
 *
 *  Created on: 17 Feb 2022
 *      Author: linuxu
 */

#ifndef CALCULATION_H_
#define CALCULATION_H_
void get_scores(double* weights,char* seq1,char** seqs,int* sizes,int num_of_sequences,Mutant* scores,int rank, MPI_Status status);
void useOpenMp_CUDA(int start,int end,char* seq1,char* seq2,double* weights,Mutant* bestMutant);
//void calculate_score(Mutant* mutant,char* seq1,char* seq2,double* weights);
//int conservative_check(char c1,char c2);
//int semi_conservative_check(char c1,char c2);
void max_in_thread_array(Mutant* mutants,Mutant* best);
void CUDA(Mutant* all_mutants_in_offset,int all_mutants_size,char* seq1,char* seq2,int seq2_size,double* weights,Mutant* best_of_offset);
#endif /* CALCULATION_H_ */
