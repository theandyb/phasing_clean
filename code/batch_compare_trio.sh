#!/bin/bash
#
#SBATCH --job-name=compare
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --array=1-602
#SBATCH --constraint=avx2
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing_clean/output/trio_phase_15/slurm/check.%A.%a.err
#SBATCH --output=/net/snowwhite/home/beckandy/research/phasing_clean/output/trio_phase_15/slurm/check.%A.%a.out

out_dir="/net/snowwhite/home/beckandy/research/phasing_clean/output/trio_phase_15/"
beagle_vcf="$out_dir/beagle/sample_${SLURM_ARRAY_TASK_ID}.vcf.gz"
eagle_vcf="$out_dir/eagle/sample_${SLURM_ARRAY_TASK_ID}.vcf.gz"
shapeit_vcf="$out_dir/shapeit/sample_${SLURM_ARRAY_TASK_ID}.vcf.gz"
truth_vcf="$out_dir/truth/sample_${SLURM_ARRAY_TASK_ID}.vcf.gz"

#het_sites="$out_dir/het_pos/sample_${SLURM_ARRAY_TASK_ID}.bed"

# Count number of heterozygous sites within each region
## Generate "bed"-like file with sites
# Get list of heterozygouse sites
bcftools query -f '%CHROM\t%POS\t[%GT\t]\n' $truth_vcf |\
awk '{if($3 == "1|0" || $3 == "0|1")print($0)}' > $out_dir/het_loc/pair_${SLURM_ARRAY_TASK_ID}_het_loc.txt

# VCFTOOLS DIFF
vcftools --gzvcf $eagle_vcf \
  --gzdiff $truth_vcf \
  --diff-switch-error \
  --out ${out_dir}switch_errors/eagle/error_${SLURM_ARRAY_TASK_ID}

vcftools --gzvcf $shapeit_vcf \
  --gzdiff $truth_vcf \
  --diff-switch-error \
  --out ${out_dir}switch_errors/shapeit/error_${SLURM_ARRAY_TASK_ID}

vcftools --gzvcf $beagle_vcf \
  --gzdiff $truth_vcf \
  --diff-switch-error \
  --out  ${out_dir}switch_errors/beagle/error_${SLURM_ARRAY_TASK_ID}

echo "Done\!"
