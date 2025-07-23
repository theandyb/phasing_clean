#!/bin/bash
#
#SBATCH --job-name=ANDY
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --mem=8000MB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --array=0-13243%602
#SBATCH --constraint=avx2
#SBATCH --exclude=r[6320,6319,6333,6331,6323,6313,6332]
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing_clean/output/slurm/sample.%A.%a.err
#SBATCH --output=/net/snowwhite/home/beckandy/research/phasing_clean/output/slurm/sample.%A.%a.out

# Calculate the chromosome and job index based on the array index
total_jobs_per_chromosome=602
chromosome_index=$(( $SLURM_ARRAY_TASK_ID / total_jobs_per_chromosome ))
job_index_within_chromosome=$(( $SLURM_ARRAY_TASK_ID % total_jobs_per_chromosome ))

chrom=$(( chromosome_index + 1 ))
sample=$(( job_index_within_chromosome + 1 ))

echo "Running job $sample for chromosome $chrom"

input_dir="/net/snowwhite/home/beckandy/research/phasing_clean/output/trio_phase_${chrom}/switch_errors/beagle"
output_dir="${input_dir}/annotated"

if [ ! -e "${output_dir}/switch_${sample}.csv" ]; then
  python /net/snowwhite/home/beckandy/research/phasing_clean/code/append_cpg.py -c $chrom \
    -s ${input_dir}/error_${sample}.diff.switch \
    -o ${output_dir}/switch_${sample}.csv
fi


input_dir="/net/snowwhite/home/beckandy/research/phasing_clean/output/trio_phase_${chrom}/switch_errors/eagle"
output_dir="${input_dir}/annotated"

if [ ! -e "${output_dir}/switch_${sample}.csv" ]; then
  python /net/snowwhite/home/beckandy/research/phasing_clean/code/append_cpg.py -c $chrom \
    -s ${input_dir}/error_${sample}.diff.switch \
    -o ${output_dir}/switch_${sample}.csv
fi


input_dir="/net/snowwhite/home/beckandy/research/phasing_clean/output/trio_phase_${chrom}/switch_errors/shapeit"
output_dir="${input_dir}/annotated"

if [ ! -e "${output_dir}/switch_${sample}.csv" ]; then
  python /net/snowwhite/home/beckandy/research/phasing_clean/code/append_cpg.py -c $chrom \
    -s ${input_dir}/error_${sample}.diff.switch \
    -o ${output_dir}/switch_${sample}.csv
fi

