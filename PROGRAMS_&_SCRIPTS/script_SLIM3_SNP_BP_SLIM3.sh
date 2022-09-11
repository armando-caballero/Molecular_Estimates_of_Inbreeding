#!/bin/bash
#$ -cwd

rm script_SLIM3_SNP_BP_SLIM3.sh.*

if [ $# -ne 1 ]
then
	echo "Usage: $0 <n>" 
	exit 1
fi

#Set arguments
n=$1

WDIR=$PWD

#Output directory
if [ -d "$WDIR/OUTPUT_SLIM_$n" ]
then
rm -r $WDIR/OUTPUT_SLIM_$n
fi

mkdir -p $WDIR/OUTPUT_SLIM_$n

mkdir -p /state/partition1/armandoSLIM$n/$SLURM_JOBID/

###################### TRANSFER TO state/partition1 #########################

cp SNP_BP_SLIM3 /state/partition1/armandoSLIM$n/$SLURM_JOBID/
cp slim3_INPUT_N500_s0.1_h0.2_var /state/partition1/armandoSLIM$n/$SLURM_JOBID/
cp transfer_sel_dom.pl /state/partition1/armandoSLIM$n/$SLURM_JOBID/

touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`
cd /state/partition1/armandoSLIM$n/$SLURM_JOBID

###################### run SLIM3 #########################

START=$(date +%s)
module load SLiM/3.3.2
slim -t -m -d "mean_s=-0.1" -d "beta_s=1.0" -d "mean_a=-0.1" -d "beta_a=1.0" -d "rho_value=0.999" -d "ave_hs=0.2" -d "Ne_pop=500" -d "prop_deleterious=0.05" -d "size_genome=2250000000" -d "file_output_slim='./slimout'" -d "file_output_pars='./slim_results_tableEff.txt'" -d "phen_out='0'" ./slim3_INPUT_N500_s0.1_h0.2_var
#slim -t -m -d "mean_s=-0.000001" -d "beta_s=1000000.0" -d "mean_a=-0.1" -d "beta_a=1.0" -d "rho_value=0.0001" -d "ave_hs=0.0" -d "Ne_pop=500" -d "prop_deleterious=0.05" -d "size_genome=2250000000" -d "file_output_slim='./slimout'" -d "file_output_pars='./slim_results_tableEff.txt'" -d "phen_out='0'" ./slim3_INPUT_N500_s0.1_h0.2_var
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "SLIM3 took 	$DIFF seconds" >> timefile

cp -r /state/partition1/armandoSLIM$n/$SLURM_JOBID/slimout $WDIR/OUTPUT_SLIM_$n

###################### SNP_BP #########################

num=$RANDOM
echo "$num" > seedfile

START=$(date +%s)
./SNP_BP_SLIM3<<@
-99
500	N
@
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "SNP_BP took 	$DIFF seconds" >> timefile

######################## TRANSFER OF FILES TO DIRECTORY #########################

cp -r /state/partition1/armandoSLIM$n/$SLURM_JOBID/dataBP.ped $WDIR/OUTPUT_SLIM_$n
cp -r /state/partition1/armandoSLIM$n/$SLURM_JOBID/dataBP.map $WDIR/OUTPUT_SLIM_$n
cp -r /state/partition1/armandoSLIM$n/$SLURM_JOBID/qtBP.phe $WDIR/OUTPUT_SLIM_$n
cp -r /state/partition1/armandoSLIM$n/$SLURM_JOBID/list* $WDIR/OUTPUT_SLIM_$n
cp -r /state/partition1/armandoSLIM$n/$SLURM_JOBID/outfileSLIM $WDIR/OUTPUT_SLIM_$n
cp -r /state/partition1/armandoSLIM$n/$SLURM_JOBID/slim_results_tableEff.txt $WDIR/OUTPUT_SLIM_$n
cp -r /state/partition1/armandoSLIM$n/$SLURM_JOBID/timefile $WDIR/OUTPUT_SLIM_$n

######################## state/partition1 CLEANING #########################

rm -r /state/partition1/armandoSLIM$n/$SLURM_JOBID/
rm $WDIR/$SLURM_JOBID.*
