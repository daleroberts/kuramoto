#!/usr/bin/env bash
#PBS -P v10
#PBS -q express
#PBS -l walltime=00:10:00
#PBS -l ncpus=80
#PBS -l mem=20GB
#PBS -j oe
#PBS -M dale.roberts@anu.edu.au
#PBS -m abe

if [ "`type -t module`" = 'function' ]; then
    module load openmpi/1.6.3
    module load gnuplot
fi

cd ~/src/kuramoto

INPUT=params.csv
OLDIFS=$IFS
IFS=,

[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }

perl -pi -e 's/\r\n|\n|\r/\n/g' $INPUT

sed 1d $INPUT | while read i graphfile seed ngraphs npaths nsteps alpha lambda sigma K max_t outfile
do
    M=`wc -l $graphfile | awk '{print $1}'`
    echo "$i $graphfile seed:$seed ngraphs:$M npaths:$npaths nsteps:$nsteps alpha:$alpha lambda:$lambda sigma:$sigma K:$K max_t:$max_t outfile:$outfile" | tee doit.status
    mpirun -bind-to-none -stdin none -np $PBS_NCPUS ./kuramoto_mpi $graphfile $M $seed $npaths $nsteps $alpha $lambda $sigma $K $max_t > "$outfile.txt"
    export outfile=$outfile
    gnuplot -e "outfile='$outfile.png'" -e "infile='$outfile.txt'" plot.gnuplot
done
IFS=$OLDIFS

