# Concorrente

Para compilar

$ gcc linear-algebra/kernels/3mm/3mm.c utilities/polybench.c -o obj/3mm -std=gnu99 -Wall -pthread -fopenmp -Iutilities -Ilinear-algebra/kernels/3mm -DSMALL_DATASET -O0 -lm

Para executar diretamente

$ ./obj/3mm numThreads tipoExecucao

Para executar script

$ ./speedup.sh 3mm numThreads parGNUPlot

NumThreads   -> potências de 2
tipoExecução -> 0 Sequencial
                1 PThreads
                2 OpenMP
parGNUPlot   -> .2 (recomendado)
