#include <string.h>
#include <stdio.h>
#include <stdlib.h>
typedef struct Mutant
{
    int offset;
    double score;
    int mut_start;
    int mut_end;
}Mutant;
__global__  void start_calculation(Mutant* mutants, int all_mutant_size, char* seq1, char* seq2, double* weights,int seq2_size);
__device__ void calculate_score(Mutant* mutant,char* seq1,char* seq2,double* weights,int seq2_size);
__device__ int conservative_check(char c1,char c2);
__device__ int semi_conservative_check(char c1,char c2);
__host__ void get_max(Mutant* mutants,int all_mutants_size,Mutant* pointer_best);
void CUDA(Mutant* all_mutants_in_offset,int all_mutants_size,char* seq1,char* seq2,int seq2_size,double* weights,Mutant* best_of_offset);

