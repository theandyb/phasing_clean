#!/bin/bash
#
#SBATCH --job-name=compTrio
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --array=0-13243%602
#SBATCH --constraint=avx2
#SBATCH --exclude=r[6320,6319,6333,6331,6323,6332]
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing_clean/output/slurm/compare.%A.%a.err
#SBATCH --output=/net/snowwhite/home/beckandy/research/phasing_clean/output/slurm/compare.%A.%a.out

# Calculate the chromosome and job index based on the array index
total_jobs_per_chromosome=602
chromosome_index=$(( $SLURM_ARRAY_TASK_ID / total_jobs_per_chromosome ))
job_index_within_chromosome=$(( $SLURM_ARRAY_TASK_ID % total_jobs_per_chromosome ))

chrom=$(( chromosome_index + 1 ))
sample=$(( job_index_within_chromosome + 1 ))

echo "Running job $sample for chromosome $chrom"
out_dir="/net/snowwhite/home/beckandy/research/phasing_clean/output/trio_phase_${chrom}"
beagle_vcf="$out_dir/beagle/sample_${sample}.vcf.gz"
eagle_vcf="$out_dir/eagle/sample_${sample}.vcf.gz"
shapeit_vcf="$out_dir/shapeit/sample_${sample}.vcf.gz"
truth_vcf="$out_dir/truth/sample_${sample}.vcf.gz"

#het_sites="$out_dir/het_pos/sample_${SLURM_ARRAY_TASK_ID}.bed"

# Get list of heterozygouse sites
bcftools query -f '%CHROM\t%POS\t[%GT\t]\n' $truth_vcf |\
awk '{if($3 == "1|0" || $3 == "0|1")print($0)}' > $out_dir/het_loc/pair_${sample}_het_loc.txt

# VCFTOOLS DIFF
vcftools --gzvcf $eagle_vcf \
  --gzdiff $truth_vcf \
  --diff-switch-error \
  --out ${out_dir}/switch_errors/eagle/error_${sample}

vcftools --gzvcf $shapeit_vcf \
  --gzdiff $truth_vcf \
  --diff-switch-error \
  --out ${out_dir}/switch_errors/shapeit/error_${sample}

vcftools --gzvcf $beagle_vcf \
  --gzdiff $truth_vcf \
  --diff-switch-error \
  --out  ${out_dir}/switch_errors/beagle/error_${sample}

echo "Done\!"
