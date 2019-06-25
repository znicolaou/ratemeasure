#!/bin/bash
#SBATCH -A p30575
#SBATCH -n 50
#SBATCH -p short
#SBATCH -t 04:00:00
#SBATCH --mem=15000
#SBATCH --array=0-99
#SBATCH --output=outs/d0_%a.out
cd /projects/p30575/cantera
source activate cantera_env

ZGN_experiment=dilute0
ZGN_start=10

if [ -f ${ZGN_experiment}/norms${SLURM_ARRAY_TASK_ID}.dat ]; then exit; fi

mkdir -p ${ZGN_experiment}/${SLURM_ARRAY_TASK_ID}

seedmin=$((1+10*SLURM_ARRAY_TASK_ID))
seedmax=$((10+10*SLURM_ARRAY_TASK_ID))
for reac in {37,52}; do
for keep in `seq $ZGN_start 15 160`; do
for remove in `seq $ZGN_start 15 160`; do
for seed in `seq $seedmin $seedmax`; do
echo $reac $keep $remove $seed
cp ${ZGN_experiment}/${ZGN_experiment}ms.dat ${ZGN_experiment}/${SLURM_ARRAY_TASK_ID}/${reac}_${keep}_${remove}_${seed}ms.dat
cp ${ZGN_experiment}/${ZGN_experiment}mi.dat ${ZGN_experiment}/${SLURM_ARRAY_TASK_ID}/${reac}_${keep}_${remove}_${seed}mi.dat

js=`jobs | wc -l`
while [ $js -ge 50 ]; do
sleep 1
js=`jobs | wc -l`
done

if [ ! -f ${ZGN_experiment}/${SLURM_ARRAY_TASK_ID}/${reac}_${keep}_${remove}_${seed}norms.dat ]; then
/projects/p30575/anaconda3/envs/cantera_env/bin/python experiment.py 5e3 ${ZGN_experiment}.dat $reac $remove $keep $seed 0.001 0 ${ZGN_experiment}/${SLURM_ARRAY_TASK_ID}/${reac}_${keep}_${remove}_${seed} &
fi

done
done
done
done
wait

for reac in {37,52}; do
for keep in `seq $ZGN_start 15 160`; do
for remove in `seq $ZGN_start 15 160`; do
for seed in `seq $seedmin $seedmax`; do
cat ${ZGN_experiment}/${SLURM_ARRAY_TASK_ID}/${reac}_${keep}_${remove}_${seed}norms.dat >> ${ZGN_experiment}/norms${SLURM_ARRAY_TASK_ID}.dat
rm ${ZGN_experiment}/${SLURM_ARRAY_TASK_ID}/${reac}_${keep}_${remove}_${seed}norms.dat
rm ${ZGN_experiment}/${SLURM_ARRAY_TASK_ID}/${reac}_${keep}_${remove}_${seed}ms.dat
rm ${ZGN_experiment}/${SLURM_ARRAY_TASK_ID}/${reac}_${keep}_${remove}_${seed}mi.dat
done
done
done
done
rmdir ${ZGN_experiment}/${SLURM_ARRAY_TASK_ID}

source deactivate
