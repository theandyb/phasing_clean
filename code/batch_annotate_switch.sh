#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=8000MB
#SBATCH --ntasks=1
#SBATCH --time 01:00:00
#SBATCH --job-name=SAn
#SBATCH --array=290
#SBATCH --requeue
#SBATCH -p main
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing_clean/output/trio_phase_1/slurm/ann-%A_%a.err
#SBATCH -o /net/snowwhite/home/beckandy/research/phasing_clean/output/trio_phase_1/slurm/ann-%A_%a.out

chrom=5
input_dir="/net/snowwhite/home/beckandy/research/phasing_clean/output/trio_phase_${chrom}/switch_errors/beagle"
output_dir="${input_dir}/annotated"

python /net/snowwhite/home/beckandy/research/phasing_clean/code/append_cpg.py -c $chrom \
  -s ${input_dir}/error_${SLURM_ARRAY_TASK_ID}.diff.switch \
  -o ${output_dir}/switch_${SLURM_ARRAY_TASK_ID}.csv

input_dir="/net/snowwhite/home/beckandy/research/phasing_clean/output/trio_phase_${chrom}/switch_errors/eagle"
output_dir="${input_dir}/annotated"

python /net/snowwhite/home/beckandy/research/phasing_clean/code/append_cpg.py -c $chrom \
  -s ${input_dir}/error_${SLURM_ARRAY_TASK_ID}.diff.switch \
  -o ${output_dir}/switch_${SLURM_ARRAY_TASK_ID}.csv

input_dir="/net/snowwhite/home/beckandy/research/phasing_clean/output/trio_phase_${chrom}/switch_errors/shapeit"
output_dir="${input_dir}/annotated"

python /net/snowwhite/home/beckandy/research/phasing_clean/code/append_cpg.py -c $chrom \
  -s ${input_dir}/error_${SLURM_ARRAY_TASK_ID}.diff.switch \
  -o ${output_dir}/switch_${SLURM_ARRAY_TASK_ID}.csv
