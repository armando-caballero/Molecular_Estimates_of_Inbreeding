#!/bin/bash
#$ -cwd

rm script_EstimationMolFNeW.sh.*

if [ $# -ne 7 ]  
then
	echo "Usage: $0 <TNIND> <NIND> <EC> <G> <genBP> <REPS> <n>" 
	exit 1
fi

#Set arguments
TNIND=$1
NIND=$2
EC=$3
G=$4
genBP=$5
REPS=$6
n=$7

WDIR=$PWD 
if [ -d "RESULTS.T$TNIND.N$NIND.$n" ]
then
rm -r RESULTS.T$TNIND.N$NIND.$n
fi
mkdir -p $WDIR/RESULTS.T$TNIND.N$NIND.$n

mkdir -p /state/partition1/EstMolNe$n/$SLURM_JOBID/
cp EstimationMolFNeIBD /state/partition1/EstMolNe$n/$SLURM_JOBID/
cp dataBP_W.map /state/partition1/EstMolNe$n/$SLURM_JOBID/dataBP.map
cp dataBP_W.ped /state/partition1/EstMolNe$n/$SLURM_JOBID/dataBP.ped
cp list_allsnps_W /state/partition1/EstMolNe$n/$SLURM_JOBID/list_allsnps
touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`
cd /state/partition1/EstMolNe$n/$SLURM_JOBID

num=$RANDOM
echo "$num" > seedfile

time ./EstimationMolFNeIBD>>out<<@
0
-99
500	NP
$TNIND	TNIND
$NIND	NIND
0.1	var_effects
0.5	ha (positive effect)
2.0	VE
18	L (99: free recombination)
$G	generations
$genBP		genBP
$EC	EC (1), RM (0)
0	neutral (1) for fitness
0	INITIALFREQ05
0.05	MAF
$REPS	replicates
@

echo "gen    ind   QT   Fped   Fibd    Fvr1BP  Fvr2BP  Fyang1BP Fyang2BP FLH1BP  FLH2BP  FhomBP" > kk1
mv list_individual_Fs kk2
cat kk1 kk2 > list_individual_Fs
mv list_individual_FsBP kk2
cat kk1 kk2 > list_individual_FsBP

cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/genfile.dat $WDIR/RESULTS.T$TNIND.N$NIND.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/outfileF.dat $WDIR/RESULTS.T$TNIND.N$NIND.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/outfileID.dat $WDIR/RESULTS.T$TNIND.N$NIND.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/outfileNe.dat $WDIR/RESULTS.T$TNIND.N$NIND.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/outfileR.dat $WDIR/RESULTS.T$TNIND.N$NIND.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/outfileVF.dat $WDIR/RESULTS.T$TNIND.N$NIND.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/outfileFBP.dat $WDIR/RESULTS.T$TNIND.N$NIND.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/outfileIDBP.dat $WDIR/RESULTS.T$TNIND.N$NIND.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/outfileNeBP.dat $WDIR/RESULTS.T$TNIND.N$NIND.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/outfileRBP.dat $WDIR/RESULTS.T$TNIND.N$NIND.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/outfileVFBP.dat $WDIR/RESULTS.T$TNIND.N$NIND.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/outfileFLASTGEN.dat $WDIR/RESULTS.T$TNIND.N$NIND.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/outfileVFLASTGEN.dat $WDIR/RESULTS.T$TNIND.N$NIND.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/frequencies.dat $WDIR/RESULTS.T$TNIND.N$NIND.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/list* $WDIR/RESULTS.T$TNIND.N$NIND.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/GRAPHS $WDIR/RESULTS.T$TNIND.N$NIND.$n
#cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/repfile.dat $WDIR/RESULTS.T$TNIND.N$NIND.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/dfilename* $WDIR/RESULTS.T$TNIND.N$NIND.$n
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/data.map $WDIR/RESULTS.T$TNIND.N$NIND.$n/data$n.map
cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/data.ped $WDIR/RESULTS.T$TNIND.N$NIND.$n/data$n.ped

cp -r /state/partition1/EstMolNe$n/$SLURM_JOBID/out $WDIR/RESULTS.T$TNIND.N$NIND.$n

rm -r /state/partition1/EstMolNe$n/$SLURM_JOBID/
rm $WDIR/$SLURM_JOBID.*