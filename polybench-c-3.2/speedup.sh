#! /bin/bash

declare -a RESULTS #vetor de resultados (tempos)
PROG=$1 #nome do programa
THREADS=$(( $2 + 0 )) #numero de threads inicial
YTIC=$3 #gnuplot - intervalo de y
NOTSEQ=1 #flag para executar sequencial apenas 1x (1 thread)
#NPROC=$(grep -m 1 name /proc/cpuinfo | sed 's/^.*: //' | sed 's/@/\\\\\\@/') #pega nome do processador

rm -rf results/$PROG.* #reseta resultados salvos do arquivo txt .dat

while [[ "$THREADS" -ge "2" ]]; do #executa até num threads ser 2
	for i in $(seq 0 2);	do #escolha de tipo de execucao - S P O
		if [[ "$NOTSEQ" -eq "0" ]]; then #checa se já rodou sequencial
			if [[ "$i" -eq "0" ]]; then #se nao rodou, executa sequencial
				echo -e "\n\x1B[35m$THREADS threads:"
				continue #senao, pula
			fi
		fi
		SUM=0 #sum iniciada em 0
		RESULTS=() #reset o vetor
		for j in $(seq 1 13); do
			ID=$(( 14 - j ))
			if [ "$j" -eq "1" ]; then
				MIN=$(./obj/$PROG $THREADS $i)
				MIN=$(( MIN + 0 ))
				MAX=$(./obj/$PROG $THREADS $i)
				MAX=$(( MAX + 0 ))
			else
				RESULTS[ID]=$(./obj/$PROG $THREADS $i)
				RESULTS[ID]=$(( ${RESULTS[$ID]} + 0 ))
				if [ "${RESULTS[$ID]}" -le "$MIN" ]; then
					MIN=$ID
				fi
				if [ "${RESULTS[$ID]}" -ge "$MAX" ]; then
					MAX=$ID
				fi
			fi
		done

		for j in $(seq 1 12); do
			if [ "${RESULTS[$j]}" -ne "$MIN" ]; then
				if [ "${RESULTS[$j]}" -ne "$MAX" ]; then
					SUM=$(( SUM + ${RESULTS[$j]} + 0 ))
				fi
			fi
		done
		DEC=$(bc -l <<< "($SUM / 10000000000)")
		if [[ "$i" -eq "0" ]]; then
			SEQ=$DEC
			DEC=$(bc -l <<< "$DEC * 10")
			echo -e "\n\n\t\t\x1B[33mSequential:\t\x1B[31mTime: $DEC\n"
			echo -e "\x1B[35m$THREADS threads:"
		else
			if [[ "$i" -eq "1" ]]; then
				SPD=$(bc -l <<< "($SEQ / $DEC)")
				DEC=$(bc -l <<< "$DEC * 10")
				echo -e "\t\x1B[33mPthread:\t\x1B[31mTime: $DEC\t\x1B[36mSpeedup: $SPD"
				printf "%s\t%s\t" "$THREADS" "$SPD" >> results/$PROG.dat
			else
				SPD=$(bc -l <<< "($SEQ / $DEC)")
				DEC=$(bc -l <<< "$DEC * 10")
				echo -e "\t\x1B[33mOpenMP:\t\t\x1B[31mTime: $DEC\t\x1B[36mSpeedup: $SPD"
				printf "%s\t\n" "$SPD" >> results/$PROG.dat
			fi
		fi
	done
	NOTSEQ=0
	THREADS=$(( $THREADS / 2 ))
done
echo -e "\x1B[0m\n"

 gnuplot -c plot.gp $PROG $YTIC "$NPROC" $SIZE
