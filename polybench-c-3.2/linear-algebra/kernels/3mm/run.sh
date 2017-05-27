#! /bin/bash
PROGRAM="3mmPapi"

gcc -I ../../../utilities/ -I /usr/local/include/  ../../../utilities/polybench.c 3mmPapi.c /usr/local/lib/libpapi.a -DLARGE_DATASET -lm -pthread -fopenmp -o 3mmPapi.out

echo "JULIANOOOOOOO"

# SEQUENTIAL
#file="./data/${PROGRAM}_sequencial.dat"
#if [ -f "$file" ]
#then
#	rm $file
#fi
#for i in {1..11}
#do
#   echo "SEQUENTIAL POLYBENCH - EXECUTION $i"
#   if [[ $i == 1 ]]; then
#     ./${PROGRAM}.out 32 0 > /dev/null
#   else
    # ./${PROGRAM}.out $MSIZE 2 >> testfile.data
#     ./${PROGRAM}.out 32 0 >> ./data/${PROGRAM}_sequencial.dat 2>&1
#  fi
#done


# OMP
#for t in 2 4 8 16
#do
#  file="./data/${PROGRAM}_openmp_$t.dat"
#  if [ -f "$file" ]
#  then
#  	rm $file
#  fi
#  for i in {1..11}
#  do
#		 echo "OMP FOR $t THREADS - EXECUTION $i"
#		 if [[ $i == 1 ]]; then
#	       ./${PROGRAM}.out $t 2 > /dev/null
#	     else
#	      # ./${PROGRAM}.out $MSIZE 2 >> testfile.data
#	       ./${PROGRAM}.out $t 2 >> ./data/${PROGRAM}_openmp_$t.dat 2>&1
#	    fi
#  done
#done

# PTHREAD
for t in 16
do
	file="./data/${PROGRAM}_pthread_$t.dat"
	if [ -f "$file" ]
	then
		rm $file
	fi
	for i in {1..11}
	do
		echo "PTHREAD FOR $t THREADS - EXECUTION $i"
		if [[ $i == 1 ]]; then
	      ./${PROGRAM}.out $t 1 > /dev/null
	    else
	     # ./${PROGRAM}.out $MSIZE 2 >> testfile.data
	      ./${PROGRAM}.out $t 1 >> ./data/${PROGRAM}_pthread_$t.dat 2>&1
	   fi
	done
done
