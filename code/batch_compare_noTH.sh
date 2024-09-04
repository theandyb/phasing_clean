#!/bin/bash
#
#SBATCH --job-name=compare
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --array=1-602
#SBATCH --constraint=avx2
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing_clean/output/trio_phase_1/slurm/check_noTH.%A.%a.err
#SBATCH --output=/net/snowwhite/home/beckandy/research/phasing_clean/output/trio_phase_1/slurm/check_noTH.%A.%a.out

chrom=1
base_dir="/net/snowwhite/home/beckandy/research/phasing_clean"
out_dir="${base_dir}/output/trio_phase_${chrom}/no_th"
beagle_vcf="${base_dir}/output/trio_phase_${chrom}/beagle/sample_${SLURM_ARRAY_TASK_ID}.vcf.gz"
eagle_vcf="${base_dir}/output/trio_phase_${chrom}/eagle/sample_${SLURM_ARRAY_TASK_ID}.vcf.gz"
shapeit_vcf="${base_dir}/output/trio_phase_${chrom}/shapeit/sample_${SLURM_ARRAY_TASK_ID}.vcf.gz"
truth_vcf="${base_dir}/output/trio_phase_${chrom}/truth/sample_${SLURM_ARRAY_TASK_ID}.vcf.gz"
triple_het_file="${base_dir}/data/1kgp/chr${chrom}/triple_het/sample_${SLURM_ARRAY_TASK_ID}.tsv"

beagle_new="$out_dir/vcf/beagle_${SLURM_ARRAY_TASK_ID}.vcf.gz"
eagle_new="$out_dir/vcf/eagle_${SLURM_ARRAY_TASK_ID}.vcf.gz"
shapeit_new="$out_dir/vcf/shapeit_${SLURM_ARRAY_TASK_ID}.vcf.gz"
truth_new="$out_dir/vcf/truth_${SLURM_ARRAY_TASK_ID}.vcf.gz"

#het_sites="$out_dir/het_pos/sample_${SLURM_ARRAY_TASK_ID}.bed"

# mask our vcfs
bcftools view -T ^${triple_het_file} -Oz ${beagle_vcf} > $beagle_new
bcftools view -T ^${triple_het_file} -Oz ${eagle_vcf} > $eagle_new
bcftools view -T ^${triple_het_file} -Oz ${shapeit_vcf} > $shapeit_new
bcftools view -T ^${triple_het_file} -Oz ${truth_vcf} > $truth_new

# Get list of heterozygouse sites
bcftools query -f '%CHROM\t%POS\t[%GT\t]\n' $truth_new |\
awk '{if($3 == "1|0" || $3 == "0|1")print($0)}' > $out_dir/het_loc/pair_${SLURM_ARRAY_TASK_ID}_het_loc.txt

# VCFTOOLS DIFF
vcftools --gzvcf $eagle_new \
  --gzdiff $truth_new \
  --diff-switch-error \
  --out ${out_dir}/switch_errors/eagle/error_${SLURM_ARRAY_TASK_ID}

vcftools --gzvcf $shapeit_new \
  --gzdiff $truth_new \
  --diff-switch-error \
  --out ${out_dir}/switch_errors/shapeit/error_${SLURM_ARRAY_TASK_ID}

vcftools --gzvcf $beagle_new \
  --gzdiff $truth_new \
  --diff-switch-error \
  --out  ${out_dir}/switch_errors/beagle/error_${SLURM_ARRAY_TASK_ID}

echo "Done\!"
