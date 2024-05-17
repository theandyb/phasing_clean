#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=8000MB
#SBATCH --ntasks=1
#SBATCH --time 01:00:00
#SBATCH --job-name=sAn
#SBATCH --array=1-602
#SBATCH --requeue
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing_clean/output/trio_phase/slurm/s_ann-%A_%a.err
#SBATCH -o /net/snowwhite/home/beckandy/research/phasing_clean/output/trio_phase/slurm/s_ann-%A_%a.out

input_dir="/net/snowwhite/home/beckandy/research/phasing_clean/output/trio_phase/switch_errors/shapeit"
output_dir="${input_dir}/annotated"

python /net/snowwhite/home/beckandy/research/phasing_clean/code/append_cpg.py -c X \
  -s ${input_dir}/error_${SLURM_ARRAY_TASK_ID}.diff.switch \
  -o ${output_dir}/switch_${SLURM_ARRAY_TASK_ID}.csv
