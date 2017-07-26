# Concorrente

Para compilar

$ gcc linear-algebra/kernels/3mm/3mm.c utilities/polybench.c -o obj/3mm -std=gnu99 -Wall -pthread -fopenmp -Iutilities -Ilinear-algebra/kernels/3mm -DSMALL_DATASET -O3 -lm

Para compilar com Papi

$ gcc linear-algebra/kernels/3mm/3mmPapi.c utilities/polybench.c /usr/local/lib/libpapi.a -o obj/3mm -std=gnu99 -Wall -pthread -fopenmp -Iutilities -Ilinear-algebra/kernels/3mm -Iusr/local/include/ -DSMALL_DATASET -O3 -lm

Para compilar com MPI

$ mpicc linear-algebra/kernels/3mm/3mmpiTime.c utilities/polybench.c -o obj/3mmpiTime -std=gnu99 -Wall -Iutilities -Ilinear-algebra/kernels/3mmpiTime -DSMALL_DATASET -O3 -lm && mpirun -n 2 ./obj/3mmpiTime

Para executar diretamente

$ ./obj/3mm numThreads tipoExecucao

Para executar script

$ ./speedup.sh 3mm numThreads parGNUPlot

NumThreads   -> potências de 2
tipoExecução -> 0 Sequencial
                1 PThreads
                2 OpenMP
parGNUPlot   -> .2 (recomendado)


