#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=8000MB
#SBATCH --ntasks=1
#SBATCH --time 01:00:00
#SBATCH --job-name=SAn
#SBATCH --array=1-602
#SBATCH --requeue
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing_clean/output/vote_1/slurm/ann-%A_%a.err
#SBATCH -o /net/snowwhite/home/beckandy/research/phasing_clean/output/vote_1/slurm/ann-%A_%a.out

chrom="1"
input_dir="/net/snowwhite/home/beckandy/research/phasing_clean/output/vote_1/switches"
output_dir="${input_dir}/annotated"

python /net/snowwhite/home/beckandy/research/phasing_clean/code/append_cpg.py -c $chrom \
  -s ${input_dir}/vote_${SLURM_ARRAY_TASK_ID}.diff.switch \
  -o ${output_dir}/switch_${SLURM_ARRAY_TASK_ID}.csv

