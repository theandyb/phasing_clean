#!/bin/bash
#
#SBATCH --job-name=vote2
#SBATCH --ntasks=1
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5GB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --array=1-602
#SBATCH --constraint=avx2
#SBATCH --exclude=r[6330,6324]
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing_clean/output/vote_2/slurm/vote.%A.%a.err
#SBATCH --output=/net/snowwhite/home/beckandy/research/phasing_clean/output/vote_2/slurm/vote.%A.%a.out

chrom=2

base_dir="/net/snowwhite/home/beckandy/research/phasing_clean"
in_dir="${base_dir}/output/trio_phase_${chrom}"
out_dir="${base_dir}/output/vote_${chrom}"
tmp_dir="/tmp"

if [ ! -f "${out_dir}/error_vote_${SLURM_ARRAY_TASK_ID}.tsv" ]; then

# Index input bcf
bcftools index $in_dir/shapeit/sample_${SLURM_ARRAY_TASK_ID}.vcf.gz
bcftools index $in_dir/beagle/sample_${SLURM_ARRAY_TASK_ID}.vcf.gz
bcftools index $in_dir/eagle/sample_${SLURM_ARRAY_TASK_ID}.vcf.gz

# Get intersection of sites in the 3 VCFs
bcftools isec -n=3 -c all $in_dir/shapeit/sample_${SLURM_ARRAY_TASK_ID}.vcf.gz $in_dir/beagle/sample_${SLURM_ARRAY_TASK_ID}.vcf.gz $in_dir/eagle/sample_${SLURM_ARRAY_TASK_ID}.vcf.gz >\
  $tmp_dir/common_pos_${SLURM_ARRAY_TASK_ID}.txt

# Generate vcf.gz files for each algorithm
bcftools view -Oz -T $tmp_dir/common_pos_${SLURM_ARRAY_TASK_ID}.txt $in_dir/eagle/sample_${SLURM_ARRAY_TASK_ID}.vcf.gz > $tmp_dir/eagle_${SLURM_ARRAY_TASK_ID}.vcf.gz
bcftools view -Oz -T $tmp_dir/common_pos_${SLURM_ARRAY_TASK_ID}.txt $in_dir/beagle/sample_${SLURM_ARRAY_TASK_ID}.vcf.gz > $tmp_dir/beagle_${SLURM_ARRAY_TASK_ID}.vcf.gz
bcftools view -Oz -T $tmp_dir/common_pos_${SLURM_ARRAY_TASK_ID}.txt $in_dir/shapeit/sample_${SLURM_ARRAY_TASK_ID}.vcf.gz > $tmp_dir/shapeit_${SLURM_ARRAY_TASK_ID}.vcf.gz

# Vote!
Rscript ${base_dir}/code/vote_phase_trio.R ${SLURM_ARRAY_TASK_ID} $tmp_dir $tmp_dir/vote_${SLURM_ARRAY_TASK_ID}.vcf.gz

#whatshap
gunzip $tmp_dir/vote_${SLURM_ARRAY_TASK_ID}.vcf.gz
bcftools view -Ov $in_dir/truth/sample_${SLURM_ARRAY_TASK_ID}.vcf.gz >\
  $tmp_dir/true_${SLURM_ARRAY_TASK_ID}.vcf

whatshap compare --names truth,vote --tsv-pairwise $out_dir/error_vote_${SLURM_ARRAY_TASK_ID}.tsv $tmp_dir/true_${SLURM_ARRAY_TASK_ID}.vcf $tmp_dir/vote_${SLURM_ARRAY_TASK_ID}.vcf

# vcftools
vcftools --vcf $tmp_dir/vote_${SLURM_ARRAY_TASK_ID}.vcf \
  --diff $tmp_dir/true_${SLURM_ARRAY_TASK_ID}.vcf \
  --diff-switch-error \
  --out ${out_dir}/switches/vote_${SLURM_ARRAY_TASK_ID}

rm $tmp_dir/true_${SLURM_ARRAY_TASK_ID}.vcf
rm $tmp_dir/eagle_${SLURM_ARRAY_TASK_ID}.vcf.gz
rm $tmp_dir/beagle_${SLURM_ARRAY_TASK_ID}.vcf.gz
rm $tmp_dir/shapeit_${SLURM_ARRAY_TASK_ID}.vcf.gz
rm $tmp_dir/common_pos_${SLURM_ARRAY_TASK_ID}.txt

bcftools view -Ob $tmp_dir/vote_${SLURM_ARRAY_TASK_ID}.vcf > $out_dir/vote_${SLURM_ARRAY_TASK_ID}.bcf
rm $tmp_dir/vote_${SLURM_ARRAY_TASK_ID}.vcf

else
echo "All good for ${SLURM_ARRAY_TASK_ID}"
fi
