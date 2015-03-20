#!/usr/bin/env bash
#
#
INPUT=params.csv
OLDIFS=$IFS
IFS=,

[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }

perl -pi -e 's/\r\n|\n|\r/\n/g' $INPUT

sed 1d $INPUT | while read i graphfile seed ngraphs npaths nsteps alpha lambda sigma K max_t outfile
do
    M=`wc -l $graphfile | awk '{print $1}'`
    echo "$i $graphfile seed:$seed ngraphs:$M npaths:$npaths nsteps:$nsteps alpha:$alpha lambda:$lambda sigma:$sigma K:$K max_t:$max_t outfile:$outfile" | tee doit.status
    mpirun -stdin none -c 12 ./kuramoto_mpi $graphfile $M $seed $npaths $nsteps $alpha $lambda $sigma $K $max_t > "$outfile.txt"
    export outfile=$outfile
    gnuplot -e "outfile='$outfile.png'" -e "infile='$outfile.txt'" plot.gnuplot
done
IFS=$OLDIFS

