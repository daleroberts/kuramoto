#!/usr/bin/env bash
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l ncpus=992
#PBS -l mem=80GB
#PBS -j oe
#PBS -r y
#PBS -M dale.roberts@anu.edu.au
#PBS -m abe

NCPUS=${PBS_NCPUS-12}

if [ "`type -t module`" = 'function' ]; then
    module load openmpi/1.8.7
    module load boost/1.57.0
fi

cd ~/src/kuramoto
mkdir -p results

find ~/src/kuramoto -size 0 -print0 | xargs -0 rm

INPUT=params.csv
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }

cp $INPUT ~/src/kuramoto/results
cp *.g6 ~/src/kuramoto/results
cd ~/src/kuramoto/results

# fix line endings
perl -pi -e 's/\r\n|\n|\r/\n/g' $INPUT

# reverse simulations
#mv $INPUT params.old
#head -1 params.old > $INPUT
##cat params.old | sed 1d | tail -r >> $INPUT
#cat params.old | sed 1d | tac >> $INPUT
#rm params.old

OLDIFS=$IFS
IFS=,

sed 1d $INPUT | while read i graphfile seed ngraphs npaths nsteps alpha lambda sigma K max_t outfile
do
    M=`wc -l $graphfile | awk '{print $1}'`

    echo "$i $graphfile seed:$seed ngraphs:$ngraphs npaths:$npaths nsteps:$nsteps alpha:$alpha lambda:$lambda sigma:$sigma K:$K max_t:$max_t outfile:$outfile" | tee doit.status

    [ ! -f $graphfile ] && { echo "$graphfile file not found"; exit 99; }

    if [ ! -f "$outfile.txt" ]; then
    	mpirun -stdin none -np $NCPUS ../kuramoto_mpi $graphfile $M $seed $npaths $nsteps $alpha $lambda $sigma $K $max_t > "$outfile.txt"
    fi
done

IFS=$OLDIFS
