#!/bin/bash

#parse arguments
die() {
	printf '%s\n' "$1" >&2
	exit 1
}

_help() {
	printf "Usage: map_genes.sh [--tr1] [--t1] [--n1] [--tr2] [--t2] [--n2] [--d] [--threads] \n\t[--tr1]: path to query sequence/transcriptome/proteome 1\n\t[--t1]: is 1 a transcriptome [nucl] or proteome [prot]\n\t[--n1]: two character identifier of 1\n\t[--tr2]: path to transcriptome/proteome 2\n\t[--t2]: is 2 a transcriptome [nucl] or proteome [prot]\n\t[--n2]: two character identifier of 2\n\t[--d]: the direction of blast\n\t\t"s" presents performing the blast in single direction, \n\t\t"d" presents performing the blast in double direction\n\t[--threads]: Number of threads (CPUs) to use in the BLAST search.\n"
}

if [ $# -eq 0 ]; then
    _help
    exit 1
fi

while :; do
    case $1 in
        -h|-\?|--help)
            _help
	    exit
            ;;
        --tr1)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                tr1=$2
                shift
            else
                die 'ERROR: "--tr1" requires a non-empty option argument.'
            fi
            ;;
        --tr2)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                tr2=$2
                shift
            else
                die 'ERROR: "--tr2" requires a non-empty option argument.'
            fi
            ;;
        --n1)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                n1=$2
                shift
            else
                die 'ERROR: "--n1" requires a non-empty option argument.'
            fi
            ;;
        --n2)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                n2=$2
                shift
            else
                die 'ERROR: "--n2" requires a non-empty option argument.'
            fi
            ;;
        --t1)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                t1=$2
                shift
            else
                die 'ERROR: "--t1" requires a non-empty option argument.'
            fi
            ;;
        --t2)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                t2=$2
                shift
            else
                die 'ERROR: "--t2" requires a non-empty option argument.'
            fi
            ;;
	--d)        # 
	    if [ "$2" ]; then
		d=$2
		shift
	    else
		die 'ERROR: "--query" requires a non-empty option argument.'
	    fi
	    ;;
        --threads)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                n_threads=$2
                shift
            fi
            ;;
        --)              # End of all options.
            shift
            break
            ;;
        -?*)
            printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
            ;;
        *)               # Default case: No more options, so break out of the loop.
            break
    esac

    shift
done

if [ -z ${n_threads+x} ]; then
        n_threads=8
fi

n1="${n1:0:2}"
n2="${n2:0:2}"
tr_f1=$(basename ${tr1})
tr_f2=$(basename ${tr2})

if [[ "$d" == "s" ]]
then

mkdir -p "01.blastdb/${n2}/${t2}"
mkdir -p "02.maps/single_direction/${n1}${n2}/${t1}${t2}"

if [ ! -f "01.blastdb/${n2}/${t2}/${tr_f2}.nhr" ] && [ ! -f "01.blastdb/${n2}/${t2}/${tr_f2}.phr" ]
then
/dellfsqd2/ST_OCEAN/USER/yangqi3/01_Software/BLAST/blast-2.9.0/ncbi-blast-2.9.0+/bin/makeblastdb -in "${tr2}" -dbtype $t2 -out "01.blastdb/${n2}/${t2}/${tr_f2}"
fi

if [[ "$t1" == "nucl" && "$t2" == "nucl" ]]
then
echo "Running tblastx from n1 to n2"
/dellfsqd2/ST_OCEAN/USER/yangqi3/01_Software/BLAST/blast-2.9.0/ncbi-blast-2.9.0+/bin/tblastx -query "${tr1}" -db "01.blastdb/${n2}/${t2}/${tr_f2}" -outfmt 6 -out "02.maps/single_direction/${n1}${n2}/${t1}${t2}/${n1}_to_${n2}.txt" -num_threads ${n_threads} -max_hsps 1 -evalue 1e-6
fi

if [[ "$t1" == "nucl" && "$t2" == "prot" ]]
then
echo "Running blastx from n1 to n2"
/dellfsqd2/ST_OCEAN/USER/yangqi3/01_Software/BLAST/blast-2.9.0/ncbi-blast-2.9.0+/bin/blastx -query "${tr1}" -db "01.blastdb/${n2}/${t2}/${tr_f2}" -outfmt 6 -out "02.maps/single_direction/${n1}${n2}/${t1}${t2}/${n1}_to_${n2}.txt" -num_threads ${n_threads} -max_hsps 1 -evalue 1e-6
fi

if [[ "$t1" == "prot" && "$t2" == "nucl" ]]
then
echo "Running tblastn from n1 to n2"
/dellfsqd2/ST_OCEAN/USER/yangqi3/01_Software/BLAST/blast-2.9.0/ncbi-blast-2.9.0+/bin/tblastn -query "${tr1}" -db "01.blastdb/${n2}/${t2}/${tr_f2}" -outfmt 6 -out "02.maps/single_direction/${n1}${n2}/${t1}${t2}/${n1}_to_${n2}.txt" -num_threads ${n_threads} -max_hsps 1 -evalue 1e-6
fi

if [[ "$t1" == "prot" && "$t2" == "prot" ]]
then
echo "Running blastp from n1 to n2"
/dellfsqd2/ST_OCEAN/USER/yangqi3/01_Software/BLAST/blast-2.9.0/ncbi-blast-2.9.0+/bin/blastp -query "${tr1}" -db "01.blastdb/${n2}/${t2}/${tr_f2}" -outfmt 6 -out "02.maps/single_direction/${n1}${n2}/${t1}${t2}/${n1}_to_${n2}.txt" -num_threads ${n_threads} -max_hsps 1 -evalue 1e-6
fi

fi

if [[ "$d" == "d" ]]
then

mkdir -p "01.blastdb/${n1}/${t1}"
mkdir -p "01.blastdb/${n2}/${t2}"
mkdir -p "02.maps/double_directions/${n1}${n2}/${t1}${t2}"

if [ ! -f "01.blastdb/${n1}/${t1}/${tr_f1}.nhr" ] && [ ! -f "01.blastdb/${n1}/${t1}/${tr_f1}.phr" ]
then
/dellfsqd2/ST_OCEAN/USER/yangqi3/01_Software/BLAST/blast-2.9.0/ncbi-blast-2.9.0+/bin/makeblastdb -in "${tr1}" -dbtype $t1 -out "01.blastdb/${n1}/${t1}/${tr_f1}"
fi

if [ ! -f "01.blastdb/${n2}/${t2}/${tr_f2}.nhr" ] && [ ! -f "01.blastdb/${n2}/${t2}/${tr_f2}.phr" ]
then
/dellfsqd2/ST_OCEAN/USER/yangqi3/01_Software/BLAST/blast-2.9.0/ncbi-blast-2.9.0+/bin/makeblastdb -in "${tr2}" -dbtype $t2 -out "01.blastdb/${n2}/${t2}/${tr_f2}"
fi

if [[ "$t1" == "nucl" && "$t2" == "nucl" ]]
then
echo "Running tblastx in both directions"
/dellfsqd2/ST_OCEAN/USER/yangqi3/01_Software/BLAST/blast-2.9.0/ncbi-blast-2.9.0+/bin/tblastx -query "${tr1}" -db "01.blastdb/${n2}/${t2}/${tr_f2}" -outfmt 6 -out "02.maps/double_directions/${n1}${n2}/${t1}${t2}/${n1}_to_${n2}.txt" -num_threads ${n_threads} -max_hsps 1 -evalue 1e-6
/dellfsqd2/ST_OCEAN/USER/yangqi3/01_Software/BLAST/blast-2.9.0/ncbi-blast-2.9.0+/bin/tblastx -query "${tr2}" -db "01.blastdb/${n1}/${t1}/${tr_f1}" -outfmt 6 -out "02.maps/double_directions/${n1}${n2}/${t1}${t2}/${n2}_to_${n1}.txt" -num_threads ${n_threads} -max_hsps 1 -evalue 1e-6
fi

if [[ "$t1" == "nucl" && "$t2" == "prot" ]]
then
echo "Running blastx from 1 to 2 and tblastn from 2 to 1"
/dellfsqd2/ST_OCEAN/USER/yangqi3/01_Software/BLAST/blast-2.9.0/ncbi-blast-2.9.0+/bin/blastx -query "${tr1}" -db "01.blastdb/${n2}/${t2}/${tr_f2}" -outfmt 6 -out "02.maps/double_directions/${n1}${n2}/${t1}${t2}/${n1}_to_${n2}.txt" -num_threads ${n_threads} -max_hsps 1 -evalue 1e-6
/dellfsqd2/ST_OCEAN/USER/yangqi3/01_Software/BLAST/blast-2.9.0/ncbi-blast-2.9.0+/bin/tblastn -query "${tr2}" -db "01.blastdb/${n1}/${t1}/${tr_f1}" -outfmt 6 -out "02.maps/double_directions/${n1}${n2}/${t1}${t2}/${n2}_to_${n1}.txt" -num_threads ${n_threads} -max_hsps 1 -evalue 1e-6
fi

if [[ "$t1" == "prot" && "$t2" == "nucl" ]]
then
echo "Running tblastn from 1 to 2 and blastx from 2 to 1"
/dellfsqd2/ST_OCEAN/USER/yangqi3/01_Software/BLAST/blast-2.9.0/ncbi-blast-2.9.0+/bin/tblastn -query "${tr1}" -db "01.blastdb/${n2}/${t2}/${tr_f2}" -outfmt 6 -out "02.maps/double_directions/${n1}${n2}/${t1}${t2}/${n1}_to_${n2}.txt" -num_threads ${n_threads} -max_hsps 1 -evalue 1e-6
/dellfsqd2/ST_OCEAN/USER/yangqi3/01_Software/BLAST/blast-2.9.0/ncbi-blast-2.9.0+/bin/blastx -query "${tr2}" -db "01.blastdb/${n1}/${t1}/${tr_f1}" -outfmt 6 -out "02.maps/double_directions/${n1}${n2}/${t1}${t2}/${n2}_to_${n1}.txt" -num_threads ${n_threads} -max_hsps 1 -evalue 1e-6
fi

if [[ "$t1" == "prot" && "$t2" == "prot" ]]
then
echo "Running blastp in both directions"
/dellfsqd2/ST_OCEAN/USER/yangqi3/01_Software/BLAST/blast-2.9.0/ncbi-blast-2.9.0+/bin/blastp -query "${tr1}" -db "01.blastdb/${n2}/${t2}/${tr_f2}" -outfmt 6 -out "02.maps/double_directions/${n1}${n2}/${t1}${t2}/${n1}_to_${n2}.txt" -num_threads ${n_threads} -max_hsps 1 -evalue 1e-6
/dellfsqd2/ST_OCEAN/USER/yangqi3/01_Software/BLAST/blast-2.9.0/ncbi-blast-2.9.0+/bin/blastp -query "${tr2}" -db "01.blastdb/${n1}/${t1}/${tr_f1}" -outfmt 6 -out "02.maps/double_directions/${n1}${n2}/${t1}${t2}/${n2}_to_${n1}.txt" -num_threads ${n_threads} -max_hsps 1 -evalue 1e-6
fi

fi
