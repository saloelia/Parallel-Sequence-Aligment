build:
	mpicxx -fopenmp -c main.c -o main.o
	mpicxx -fopenmp -c utility.c -o utility.o
	mpicxx -fopenmp -c calculation.c -o calculation.o
	nvcc -I./inc -c cudaFunctions.cu -o cudaFunctions.o
	mpicxx -fopenmp -o mpiCudaOpemMP  main.o utility.o calculation.o cudaFunctions.o  /usr/local/cuda/lib64/libcudart_static.a -ldl -lrt

clean:
	rm -f *.o ./mpiCudaOpemMP

run:
	mpiexec -np 2 ./mpiCudaOpemMP < input.txt > results.txt

runOn2:
	mpiexec -np 2 -machinefile  mf  -map-by  node  ./mpiCudaOpemMP
