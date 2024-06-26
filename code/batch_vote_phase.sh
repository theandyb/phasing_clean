#!/bin/bash
#
#SBATCH --job-name=vote
#SBATCH --ntasks=1
#SBATCH --time=01:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5GB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --array=1-1000
#SBATCH --constraint=avx2
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing_clean/output/vote_x/slurm/vote_easy.%A_%a.err
#SBATCH --output=/net/snowwhite/home/beckandy/research/phasing_clean/output/vote_x/slurm/vote_easy.%A.%a.out

base_dir="/net/snowwhite/home/beckandy/research/phasing_clean"
in_dir="${base_dir}/output/switch_errors/vcf"
out_dir="${base_dir}/output/vote_x"


if [ ! -f "${out_dir}/vote_${SLURM_ARRAY_TASK_ID}.bcf" ]; then

# Index input bcf
bcftools index $in_dir/shapeit/pair_${SLURM_ARRAY_TASK_ID}.bcf
bcftools index $in_dir/beagle/pair_${SLURM_ARRAY_TASK_ID}.bcf
bcftools index $in_dir/eagle/pair_${SLURM_ARRAY_TASK_ID}.bcf

# Get intersection of sites in the 3 VCFs
bcftools isec -n=3 -c all $in_dir/shapeit/pair_${SLURM_ARRAY_TASK_ID}.bcf $in_dir/beagle/pair_${SLURM_ARRAY_TASK_ID}.bcf $in_dir/eagle/pair_${SLURM_ARRAY_TASK_ID}.bcf >\
  $out_dir/common_pos_${SLURM_ARRAY_TASK_ID}.txt

# Generate vcf.gz files for each algorithm
bcftools view -Oz -T $out_dir/common_pos_${SLURM_ARRAY_TASK_ID}.txt $in_dir/eagle/pair_${SLURM_ARRAY_TASK_ID}.bcf > $out_dir/eagle_${SLURM_ARRAY_TASK_ID}.vcf.gz
bcftools view -Oz -T $out_dir/common_pos_${SLURM_ARRAY_TASK_ID}.txt $in_dir/beagle/pair_${SLURM_ARRAY_TASK_ID}.bcf > $out_dir/beagle_${SLURM_ARRAY_TASK_ID}.vcf.gz
bcftools view -Oz -T $out_dir/common_pos_${SLURM_ARRAY_TASK_ID}.txt $in_dir/shapeit/pair_${SLURM_ARRAY_TASK_ID}.bcf > $out_dir/shapeit_${SLURM_ARRAY_TASK_ID}.vcf.gz

# Vote!
Rscript ${base_dir}/code/vote_phase.R ${SLURM_ARRAY_TASK_ID} $out_dir $out_dir/vote_${SLURM_ARRAY_TASK_ID}.vcf.gz

#whatshap
gunzip $out_dir/vote_${SLURM_ARRAY_TASK_ID}.vcf.gz
bcftools view -Ov $in_dir/pair_${SLURM_ARRAY_TASK_ID}_true.bcf >\
  $out_dir/true_${SLURM_ARRAY_TASK_ID}.vcf

whatshap compare --names truth,vote --tsv-pairwise $out_dir/error_vote_${SLURM_ARRAY_TASK_ID}.tsv $out_dir/true_${SLURM_ARRAY_TASK_ID}.vcf $out_dir/vote_${SLURM_ARRAY_TASK_ID}.vcf

# vcftools
vcftools --vcf $out_dir/vote_${SLURM_ARRAY_TASK_ID}.vcf \
  --diff $out_dir/true_${SLURM_ARRAY_TASK_ID}.vcf \
  --diff-switch-error \
  --out ${out_dir}/switches/vote_${SLURM_ARRAY_TASK_ID}

rm $out_dir/true_${SLURM_ARRAY_TASK_ID}.vcf
rm $out_dir/eagle_${SLURM_ARRAY_TASK_ID}.vcf.gz
rm $out_dir/beagle_${SLURM_ARRAY_TASK_ID}.vcf.gz
rm $out_dir/shapeit_${SLURM_ARRAY_TASK_ID}.vcf.gz
rm $out_dir/common_pos_${SLURM_ARRAY_TASK_ID}.txt

bcftools view -Ob $out_dir/vote_${SLURM_ARRAY_TASK_ID}.vcf > $out_dir/vote_${SLURM_ARRAY_TASK_ID}.bcf
rm $out_dir/vote_${SLURM_ARRAY_TASK_ID}.vcf

else
echo "All good for ${SLURM_ARRAY_TASK_ID}"
fi
