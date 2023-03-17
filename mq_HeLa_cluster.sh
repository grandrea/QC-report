#!/bin/bash

#SBATCH --job-name=MaxQuant
#SBATCH --mail-type=NONE
#SBATCH --partition=cpuq
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem 80G

FASTA=/ssu/spssu/fasta/Human_canonical_221113.fasta
JOBDIR=`pwd`
NUMTHREADS=21

source /home/andrea.graziadei/.bashrc

rm template.xml

maxquant -c template.xml

sed -i "s@example.fasta@$FASTA@g" template.xml
sed -i "s/<ibaq>False/<ibaq>True/g" template.xml
sed -i "s/<ibaqLogFit>False/<ibaqLogFit>True/g" template.xml
sed -i "s/<ibaqChargeNormalization>False</<ibaqChargeNormalization>True</g" template.xml
sed -i "s/<numThreads>1/<numThreads>$NUMTHREADS/g" template.xml
sed -i "s/<useDotNetCore>True/<useDotNetCore>False/g" template.xml
sed -i "s/<writeMsScansTable>False/<writeMsScansTable>True/g" template.xml
sed -i "s/<calcPeakProperties>False/<calcPeakProperties>True/g" template.xml

ls *.raw > list.txt
for i in `cat list.txt`
do
	FILEPATH=`echo $JOBDIR/$i`
	sed  "s@file.example.RAW@$FILEPATH@g" template.xml > current_run.xml
	maxquant current_run.xml
	cd combined
	rm -rv proc andromeda ps search ser
	cd ../
	rm -rv `echo $i | sed 's/.raw//g'`
	mv combined combined_$i
done
	

