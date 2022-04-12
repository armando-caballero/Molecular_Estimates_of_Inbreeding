#!/bin/bash
#$ -cwd

rm script_FMOLPIGSLIM.sh.*

if [ $# -ne 3 ]
then
	echo "Usage: $0 <CURRENT> <REPS> <n>" 
	exit 1
fi

#Set arguments
CURRENT=$1
REPS=$2
n=$3

WDIR=$PWD

#Output directory
if [ -d "$WDIR/OUTPUT_CERDOS_$n" ]
then
rm -r $WDIR/OUTPUT_CERDOS_$n
fi

mkdir -p $WDIR/OUTPUT_CERDOS_$n

mkdir -p /state/partition1/armandoZZ$n/$SLURM_JOBID/

###################### TRANSFER TO state/partition1 #########################

cp dataBP.* /state/partition1/armandoZZ$n/$SLURM_JOBID/
cp list_allsnps /state/partition1/armandoZZ$n/$SLURM_JOBID/
cp FMOLPIGSLIM /state/partition1/armandoZZ$n/$SLURM_JOBID/
cp GenealogiaGuadyerbas /state/partition1/armandoZZ$n/$SLURM_JOBID/
cp Animales_genotipados_1206 /state/partition1/armandoZZ$n/$SLURM_JOBID/

touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`
cd /state/partition1/armandoZZ$n/$SLURM_JOBID

###################### run FMOLPIGS #########################

START=$(date +%s)

num=$RANDOM
echo "$num" > seedfile

time ./FMOLPIGSLIM>>out<<@
0
-99
500	NINDNP
18	Lenght genome Morgans (99=free)
18	CHROM (max 50)
0	neutral(1) or not (0)
17	genBP
$CURRENT	current
1	INITIALFREQ05 (1) or not (0)
$REPS	REPS
@

cat Fmolout.dat >> FMOLOUT.DAT

echo "IND      coh   Fped    Fibd    Fvr1BP  Fvr2BP  Fyang1BP Fyang2BP FLH1BP  FLH2BP  FhomBP" > kk1
mv list_individual_Fs_23 kk2
cat kk1 kk2 > list_individual_Fs_23
mv list_individual_Fs_BP_23 kk2
cat kk1 kk2 > list_individual_Fs_BP_23

######################## TRANSFER OF FILES TO DIRECTORY #########################

cp -r /state/partition1/armandoZZ$n/$SLURM_JOBID/dfilename.dat $WDIR/OUTPUT_CERDOS_$n
cp -r /state/partition1/armandoZZ$n/$SLURM_JOBID/data.ped $WDIR/OUTPUT_CERDOS_$n
cp -r /state/partition1/armandoZZ$n/$SLURM_JOBID/data.map $WDIR/OUTPUT_CERDOS_$n
cp -r /state/partition1/armandoZZ$n/$SLURM_JOBID/list* $WDIR/OUTPUT_CERDOS_$n
cp -r /state/partition1/armandoZZ$n/$SLURM_JOBID/data.phe $WDIR/OUTPUT_CERDOS_$n
cp -r /state/partition1/armandoZZ$n/$SLURM_JOBID/FMOLOUT.DAT $WDIR/OUTPUT_CERDOS_$n
cp -r /state/partition1/armandoZZ$n/$SLURM_JOBID/outfile* $WDIR/OUTPUT_CERDOS_$n
cp -r /state/partition1/armandoZZ$n/$SLURM_JOBID/timefile $WDIR/OUTPUT_CERDOS_$n

######################## state/partition1 CLEANING #########################

rm -r /state/partition1/armandoZZ$n/$SLURM_JOBID/
rm $WDIR/$SLURM_JOBID.*
