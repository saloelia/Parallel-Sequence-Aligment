#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <cuda.h>
//#include "calculation.h"
#include "cudaFunctions.h"

__global__  void start_calculation(Mutant* mutants, int all_mutant_size, char* seq1, char* seq2, double* weights,int seq2_size) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    
    if (i < all_mutant_size){
        calculate_score(mutants+i,seq1,seq2,weights,seq2_size);
     }
}

__device__ void calculate_score(Mutant* mutant,char* seq1,char* seq2,double* weights,int seq2_size){
	int stars = 0;
	int colons = 0;
	int points = 0;
	int spaces = 0;
	int index_to_compare = 0;
	for (int i = 0; i < seq2_size; ++i) {
		if(i!=mutant->mut_start && i!=mutant->mut_end){
			if(i>mutant->mut_start && i<mutant->mut_end){
				index_to_compare = i-1;
			}else if(i > mutant->mut_end){
				index_to_compare = i-2;
			}else{
				index_to_compare = i;
			}

			if(seq1[index_to_compare]==seq2[i]){
				stars++;
			}else if(conservative_check(seq1[index_to_compare],seq2[i])){
				colons++;
			}else if(semi_conservative_check(seq1[index_to_compare],seq2[i])){
				points++;
			}else{
				spaces++;
			}
		}

	}
	(*mutant).score = weights[0]*stars - weights[1]*colons - weights[2]*points - weights[3]*spaces;
	
}

__device__ int conservative_check(char c1,char c2){
	char conservative[9][5] = {"NDEQ","MILV","FYW","NEQK","QHRK","HY","STA","NHQK","MILF"};
	for (int i = 0; i < 9; ++i) {
		int match = 0;
		for (int j = 0; j < 5; ++j) {
			if(conservative[i][j]=='\0')
				break;
			
			if(conservative[i][j]==c1 || conservative[i][j]==c2){
				match++;
			}
			if(match==2){
				
				return 1;
			}
		}
	}

	return 0;
}

__device__ int semi_conservative_check(char c1,char c2){
	char semi_conversative[11][7] = {"SAG","SGND","NEQHRK","ATV","STPA","NDEQHK","HFY","CSA","STNK","SNDEQK","FVILM"};
	for (int i = 0; i < 11; ++i) {
		int match = 0;
		for (int j = 0; j < 7; ++j) {
			if(semi_conversative[i][j]=='\0')
				break;
			
			if(semi_conversative[i][j]==c1 || semi_conversative[i][j]==c2){
				match++;
			}
			if(match==2){
				return 1;
			}
		}
	}

	return 0;
}

__host__ void get_max(Mutant* mutants,int all_mutants_size,Mutant* pointer_best){
	Mutant best = mutants[0];
	//printf("here mate?\n");
	for (int i = 1; i < all_mutants_size; ++i){
	   if(mutants[i].score > best.score){
	   	best = mutants[i];
	   }
	
	}
	pointer_best->score = best.score;
	pointer_best->mut_end = best.mut_end;
	pointer_best->mut_start = best.mut_start;
	
}


void CUDA(Mutant* all_mutants_in_offset,int all_mutants_size,char* seq1,char* seq2,int seq2_size,double* weights,Mutant* best){
    // Error code to check return values for CUDA calls
    cudaError_t err = cudaSuccess;

    //size_t size = all_mutants_size * sizeof(int);
    //printf("arraived to cuda\n");
    Mutant *cuda_mutants;
    char* cuda_seq1;
    char* cuda_seq2;
    double* cuda_weights;
    
    // Allocate memory on GPU to copy the data from the host ////////////////////////////////////////////////////////////////////////
    	    
    err = cudaMalloc((void **)&cuda_mutants, all_mutants_size * sizeof(Mutant));
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    
    err = cudaMalloc((void **)&cuda_seq1, strlen(seq1)*sizeof(char));
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
       
    err = cudaMalloc((void **)&cuda_seq2, strlen(seq2)*sizeof(char));
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    
    err = cudaMalloc((void **)&cuda_weights, 4*sizeof(double));
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    
    
    
    // Copy data from host to the GPU memory /////////////////////////////////////////////////////////////////////////////////////////
    
    err = cudaMemcpy(cuda_mutants, all_mutants_in_offset, all_mutants_size*sizeof(Mutant), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to copy mutants from host to device - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    
    err = cudaMemcpy(cuda_seq1, seq1, strlen(seq1)*sizeof(char), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to copy seq1 from host to device - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    
    err = cudaMemcpy(cuda_seq2, seq2, strlen(seq2)*sizeof(char), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to copy seq2 from host to device - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    
    err = cudaMemcpy(cuda_weights, weights, 4*sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to copy weights from host to device - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    
    
 
    // Launch the Kernel ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    int threadsPerBlock = 1024;
    int blocksPerGrid = (all_mutants_size + threadsPerBlock - 1) / threadsPerBlock;
  
    start_calculation<<<blocksPerGrid, threadsPerBlock>>>(cuda_mutants, all_mutants_size, cuda_seq1,cuda_seq2,cuda_weights,seq2_size-1);
    //printf("score %lf\n",cuda_mutants[0].score);
    
    //printf("here tits\n");
    
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to launch vectorAdd kernel -  %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    
    // Copy the  result from GPU to the host memory.////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    err = cudaMemcpy(all_mutants_in_offset, cuda_mutants, all_mutants_size*sizeof(Mutant), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to copy mutants array from device to host -%s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
         
        
    get_max(all_mutants_in_offset,all_mutants_size,best);
    

    // Free allocated memory on GPU /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    if (cudaFree(cuda_mutants) != cudaSuccess) {
        fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    
    if (cudaFree(cuda_seq1) != cudaSuccess) {
        fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    
    if (cudaFree(cuda_seq2) != cudaSuccess) {
        fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    
    if (cudaFree(cuda_weights) != cudaSuccess) {
        fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    
    
   //printf("here maaaaann\n");
    
}

