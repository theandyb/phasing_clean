#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=8000MB
#SBATCH --ntasks=1
#SBATCH --time 01:00:00
#SBATCH --job-name=annHet
#SBATCH --array=1-1000
#SBATCH --requeue
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing_clean/output/trio_phase_15/slurm/het_annT-%A_%a.err
#SBATCH -o /net/snowwhite/home/beckandy/research/phasing_clean/output/trio_phase_15/slurm/het_annT-%A_%a.out

input_dir="/net/snowwhite/home/beckandy/research/phasing_clean/output/trio_phase_15/het_loc"
output_dir="${input_dir}/annotated/"

python /net/snowwhite/home/beckandy/research/phasing_clean/code/append_cpg_het.py -c 15 \
  -s ${input_dir}/pair_${SLURM_ARRAY_TASK_ID}_het_loc.txt \
  -o ${output_dir}/pair_${SLURM_ARRAY_TASK_ID}.csv
