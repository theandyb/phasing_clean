#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=8000MB
#SBATCH --ntasks=1
#SBATCH --time 01:00:00
#SBATCH --job-name=SAn
#SBATCH --array=1-13244%200
#SBATCH --requeue
#SBATCH --constraint=avx2
#SBATCH --exclude=r[6320,6319,6333,6331,6323,6334,6335,6321,6302],dl[3616]
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing_clean/output/slurm/ann-%A_%a.err
#SBATCH -o /net/snowwhite/home/beckandy/research/phasing_clean/output/slurm/ann-%A_%a.out

NUM_SAMPLES=602
NUM_CHRS=22
CHRS_PER_TASK=$NUM_SAMPLES

CHR_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) / CHRS_PER_TASK + 1 ))
SAMPLE_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) % CHRS_PER_TASK + 1 ))

CHROMOSOME=$CHR_INDEX

echo "Starting job for Task ID: $SLURM_ARRAY_TASK_ID"
echo "Chromosome: $CHROMOSOME"
echo "Sample Index: $SAMPLE_INDEX"

input_dir="/net/snowwhite/home/beckandy/research/phasing_clean/output/vote_${CHROMOSOME}/switches"
output_dir="${input_dir}/annotated"

python /net/snowwhite/home/beckandy/research/phasing_clean/code/append_cpg.py -c $CHROMOSOME \
  -s ${input_dir}/vote_${SAMPLE_INDEX}.diff.switch \
  -o ${output_dir}/switch_${SAMPLE_INDEX}.csv

