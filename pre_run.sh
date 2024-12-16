#!/bin/bash
echo "PEP_PATH=/dellfsqd2/ST_OCEAN/USER/yangqi3/02_Project/02_nervous_system/02.cross-species/00.data/01.pep" >> batch_run.sh
echo "RESULT_PATH=/dellfsqd2/ST_OCEAN/USER/yangqi3/02_Project/02_nervous_system/02.cross-species/02.samap/01.blast" >> batch_run.sh

names=(PM CP SE AS MA OM DR PA AM PV GG MM HS)

echo "mkdir -p \${RESULT_PATH}/job_out" >> batch_run.sh
for ((i=0;i<${#names[@]};i++))
do
	for (( j=i+1; j<${#names[@]}; j++ ))
	do
		echo "qsub -cwd -l vf=10g,num_proc=8 -q st.q -o \${RESULT_PATH}/job_out -e \${RESULT_PATH}/job_out -P PARTER run_blast.sh --tr1 \${PEP_PATH}/${names[i]}_longest.pep --t1 prot --n1 ${names[i]} --tr2 \${PEP_PATH}/${names[j]}_longest.pep --t2 prot --n2 ${names[j]} --d d --threads 8" >> batch_run.sh
	done
done
