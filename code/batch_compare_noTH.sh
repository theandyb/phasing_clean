#!/bin/bash
#
#SBATCH --job-name=ANDY
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --array=0-13243%500
#SBATCH --exclude=r[6320,6319,6333,6331,6323,6313]
#SBATCH --constraint=avx2
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing_clean/output/slurm/check_noTH.%A.%a.err
#SBATCH --output=/net/snowwhite/home/beckandy/research/phasing_clean/output/slurm/check_noTH.%A.%a.out

# Calculate the chromosome and job index based on the array index
total_jobs_per_chromosome=602
chromosome_index=$(( $SLURM_ARRAY_TASK_ID / total_jobs_per_chromosome ))
job_index_within_chromosome=$(( $SLURM_ARRAY_TASK_ID % total_jobs_per_chromosome ))

chrom=$(( chromosome_index + 1 ))
sample=$(( job_index_within_chromosome + 1 ))

echo "Running job $sample for chromosome $chrom"

base_dir="/net/snowwhite/home/beckandy/research/phasing_clean"
out_dir="${base_dir}/output/trio_phase_${chrom}/no_th"
beagle_vcf="${base_dir}/output/trio_phase_${chrom}/no_th/beagle/sample_${sample}.vcf.gz"
eagle_vcf="${base_dir}/output/trio_phase_${chrom}/no_th/eagle/sample_${sample}.vcf.gz"
shapeit_vcf="${base_dir}/output/trio_phase_${chrom}/no_th/shapeit/sample_${sample}.vcf.gz"
truth_vcf="${base_dir}/output/trio_phase_${chrom}/no_th/truth/sample_${sample}.vcf.gz"
triple_het_file="${base_dir}/data/1kgp/chr${chrom}/triple_het/sample_${sample}.tsv"

#index vcf files (if not already indexed)
if [ ! -e "${beagle_vcf}.csi" ]; then
  bcftools index $beagle_vcf
fi

if [ ! -e "${eagle_vcf}.csi" ]; then
  bcftools index $eagle_vcf
fi

if [ ! -e "${shapeit_vcf}.csi" ]; then
  bcftools index $shapeit_vcf
fi

if [ ! -e "${truth_vcf}.csi" ]; then
  bcftools index $truth_vcf
fi

# beagle_new="$out_dir/vcf/beagle_${sample}.vcf.gz"
# eagle_new="$out_dir/vcf/eagle_${sample}.vcf.gz"
# shapeit_new="$out_dir/vcf/shapeit_${sample}.vcf.gz"
# truth_new="$out_dir/vcf/truth_${sample}.vcf.gz"

#het_sites="$out_dir/het_pos/sample_${SLURM_ARRAY_TASK_ID}.bed"

# # mask our vcfs
# if [ ! -e "$truth_new" ]; then
#   bcftools view -T ^${triple_het_file} -Oz ${beagle_vcf} > $beagle_new
#   bcftools view -T ^${triple_het_file} -Oz ${eagle_vcf} > $eagle_new
#   bcftools view -T ^${triple_het_file} -Oz ${shapeit_vcf} > $shapeit_new
#   bcftools view -T ^${triple_het_file} -Oz ${truth_vcf} > $truth_new
# fi

# Get list of heterozygouse sites
if [ ! -e "${out_dir}/het_loc/pair_${sample}_het_loc.txt" ]; then
bcftools query -f '%CHROM\t%POS\t[%GT\t]\n' $truth_vcf |\
awk '{if($3 == "1|0" || $3 == "0|1")print($0)}' > $out_dir/het_loc/pair_${sample}_het_loc.txt
fi

# VCFTOOLS DIFF
if [ ! -e "${out_dir}/switch_errors/eagle/error_${sample}.diff.indv.switch" ]; then
  vcftools --gzvcf $eagle_vcf \
    --gzdiff $truth_vcf \
    --diff-switch-error \
    --out ${out_dir}/switch_errors/eagle/error_${sample}
fi

if [ ! -e "${out_dir}/switch_errors/shapeit/error_${sample}.diff.indv.switch" ]; then
  vcftools --gzvcf $shapeit_vcf \
    --gzdiff $truth_vcf \
    --diff-switch-error \
    --out ${out_dir}/switch_errors/shapeit/error_${sample}
fi

if [ ! -e "${out_dir}/switch_errors/beagle/error_${sample}.diff.indv.switch" ]; then
vcftools --gzvcf $beagle_vcf \
  --gzdiff $truth_vcf \
  --diff-switch-error \
  --out  ${out_dir}/switch_errors/beagle/error_${sample}
fi

echo "Done\!"


