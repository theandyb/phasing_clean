#!/bin/bash
#
#SBATCH --job-name=tootwo
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --array=1-602
#SBATCH --partition=main
#SBATCH --constraint=avx2
#SBATCH --exclude=r[6320,6319,6333,6331]
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing_clean/output/slurm/tripleHet.%A.%a.err
#SBATCH --output=/net/snowwhite/home/beckandy/research/phasing_clean/output/slurm/tripleHet.%A.%a.out

# Code for counting triple hets
chrom=22
base_dir="/net/snowwhite/home/beckandy/research/phasing_clean"
source_vcf="${base_dir}/data/1kgp/chr${chrom}/chr${chrom}_unphased_overlap.bcf"
include_dir="${base_dir}/data/1kgp/chr1/triple_het/trio_lists" # hardcoded to chr1, since we only generate include lists once
out_dir="${base_dir}/data/1kgp/chr${chrom}/triple_het"
include_samples="${include_dir}/sample_${SLURM_ARRAY_TASK_ID}.txt"
target_sample=$(head -n ${SLURM_ARRAY_TASK_ID} /net/snowwhite/home/beckandy/research/phasing_clean/data/1kgp/child_ids.txt  | tail -1)
#working_dir="/net/snowwhite/home/beckandy/scratch"
working_dir="/tmp"

echo "Filtering vcf..."
# filter bcf to trio
bcftools view -S ${include_samples} -Ob ${source_vcf} > ${working_dir}/family_${SLURM_ARRAY_TASK_ID}.bcf

echo "Getting list of triple het positions..."
# get list of triple heterozygous positions
bcftools view -g ^hom ${working_dir}/family_${SLURM_ARRAY_TASK_ID}.bcf | \
bcftools query -f "%CHROM\t%POS\n" > ${out_dir}/sample_${SLURM_ARRAY_TASK_ID}.tsv

# cleanup
rm ${working_dir}/family_${SLURM_ARRAY_TASK_ID}.bcf
