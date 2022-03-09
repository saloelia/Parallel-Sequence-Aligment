#pragma once

#include <stdio.h>
#include <mpi.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>

#define FILE_INPUT "/home/linuxu/Downloads/finalProject/input.txt"

typedef struct Mutant
{
	int offset;
    double score;
    int mut_start;
    int mut_end;
}Mutant;
void create_MPI_mutant_type(MPI_Datatype *mutantType);
void read_file(double** weights,char** seq1,char*** seqs,int** sizes,int* num_of_sequences);
void readSequence(FILE *file,char** seq,int* size_of_seq_i);
void prepareMutants(Mutant** mutants,int n,int* num_of_mutants_pointer);
void fillMutants(Mutant* mutant_array,int n);
