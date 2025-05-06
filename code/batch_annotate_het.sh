#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=8000MB
#SBATCH --ntasks=1
#SBATCH --time 01:00:00
#SBATCH --job-name=anny
#SBATCH --array=0-13243
#SBATCH --exclude=r[6320,6319,6333,6331,6323]
#SBATCH --requeue
#SBATCH -p main
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing_clean/output/slurm/het_annT-%A_%a.err
#SBATCH -o /net/snowwhite/home/beckandy/research/phasing_clean/output/slurm/het_annT-%A_%a.out

# Calculate the chromosome and job index based on the array index
total_jobs_per_chromosome=602
chromosome_index=$(( $SLURM_ARRAY_TASK_ID / total_jobs_per_chromosome ))
job_index_within_chromosome=$(( $SLURM_ARRAY_TASK_ID % total_jobs_per_chromosome ))

chrom=$(( chromosome_index + 1 ))
sample=$(( job_index_within_chromosome + 1 ))

echo "Running job $sample for chromosome $chrom"

input_dir="/net/snowwhite/home/beckandy/research/phasing_clean/output/trio_phase_${chrom}/no_th/het_loc"
#input_dir="/net/snowwhite/home/beckandy/research/phasing_clean/output/switch_errors/het_loc"
output_dir="${input_dir}/annotated/"

python /net/snowwhite/home/beckandy/research/phasing_clean/code/append_cpg_het.py -c ${chrom} \
  -s ${input_dir}/pair_${sample}_het_loc.txt \
  -o ${output_dir}/pair_${sample}.csv
