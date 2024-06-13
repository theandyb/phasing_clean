#!/bin/bash
#
#SBATCH --job-name=tripleHet
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --array=1-602
#SBATCH --constraint=avx2
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing_clean/output/trio_phase_1/slurm/tripleHet.%A.%a.err
#SBATCH --output=/net/snowwhite/home/beckandy/research/phasing_clean/output/trio_phase_1/slurm/tripleHet.%A.%a.out

# Code for phasing 602 child samples against a panel with their parents removed
chrom=1
base_dir="/net/snowwhite/home/beckandy/research/phasing_clean"
source_vcf="${base_dir}/data/1kgp/chr${chrom}/chr${chrom}_unphased_overlap.bcf"
include_dir="${base_dir}/data/1kgp/chr${chrom}/triple_het/trio_lists"
out_dir="${base_dir}/data/1kgp/chr${chrom}/triple_het"
include_samples="${include_dir}/sample_${SLURM_ARRAY_TASK_ID}.txt"
target_sample=$(head -n ${SLURM_ARRAY_TASK_ID} /net/snowwhite/home/beckandy/research/phasing_clean/data/1kgp/child_ids.txt  | tail -1)
#working_dir="/net/snowwhite/home/beckandy/scratch"
working_dir="/tmp"

# filter bcf to trio
bcftools view -S ${include_samples} -Ob ${source_vcf} > ${working_dir}/family_${SLURM_ARRAY_TASK_ID}.bcf

# get list of triple heterozygous positions
bcftools view -g ^hom ${working_dir}/family_${SLURM_ARRAY_TASK_ID}.bcf | \
bcftools query -f "%CHROM\t%POS\n" > ${out_dir}/sample_${SLURM_ARRAY_TASK_ID}.tsv

# cleanup
rm ${working_dir}/family_${SLURM_ARRAY_TASK_ID}.bcf
