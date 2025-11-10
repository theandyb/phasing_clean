#!/bin/bash
#
#SBATCH --job-name=vote2
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5GB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --array=1-13244%200
#SBATCH --constraint=avx2
#SBATCH --exclude=r[6320,6319,6333,6331,6323,6334,6335,6321,6302],dl[3616]
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing_clean/output/slurm/vote.%A.%a.err
#SBATCH --output=/net/snowwhite/home/beckandy/research/phasing_clean/output/slurm/vote.%A.%a.out

NUM_SAMPLES=602
NUM_CHRS=22
CHRS_PER_TASK=$NUM_SAMPLES

CHR_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) / CHRS_PER_TASK + 1 ))
SAMPLE_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) % CHRS_PER_TASK + 1 ))

chrom=$CHR_INDEX

base_dir="/net/snowwhite/home/beckandy/research/phasing_clean"
in_dir="${base_dir}/output/trio_phase_${chrom}"
out_dir="${base_dir}/output/vote_${chrom}"
tmp_dir="/tmp"

if [ ! -f "${out_dir}/switches/vote_${SAMPLE_INDEX}.diff.switch" ]; then

# Index input bcf
bcftools index $in_dir/shapeit/sample_${SAMPLE_INDEX}.vcf.gz
bcftools index $in_dir/beagle/sample_${SAMPLE_INDEX}.vcf.gz
bcftools index $in_dir/eagle/sample_${SAMPLE_INDEX}.vcf.gz

# Get intersection of sites in the 3 VCFs
bcftools isec -n=3 -c all $in_dir/shapeit/sample_${SAMPLE_INDEX}.vcf.gz $in_dir/beagle/sample_${SAMPLE_INDEX}.vcf.gz $in_dir/eagle/sample_${SAMPLE_INDEX}.vcf.gz >\
  $tmp_dir/common_pos_${SAMPLE_INDEX}.txt

# Generate vcf.gz files for each algorithm
bcftools view -Oz -T $tmp_dir/common_pos_${SAMPLE_INDEX}.txt $in_dir/eagle/sample_${SAMPLE_INDEX}.vcf.gz > $tmp_dir/eagle_${SAMPLE_INDEX}.vcf.gz
bcftools view -Oz -T $tmp_dir/common_pos_${SAMPLE_INDEX}.txt $in_dir/beagle/sample_${SAMPLE_INDEX}.vcf.gz > $tmp_dir/beagle_${SAMPLE_INDEX}.vcf.gz
bcftools view -Oz -T $tmp_dir/common_pos_${SAMPLE_INDEX}.txt $in_dir/shapeit/sample_${SAMPLE_INDEX}.vcf.gz > $tmp_dir/shapeit_${SAMPLE_INDEX}.vcf.gz

# Vote!
Rscript ${base_dir}/code/vote_phase_trio.R ${SAMPLE_INDEX} $tmp_dir $tmp_dir/vote_${SAMPLE_INDEX}.vcf.gz

#whatshap
gunzip $tmp_dir/vote_${SAMPLE_INDEX}.vcf.gz
bcftools view -Ov $in_dir/truth/sample_${SAMPLE_INDEX}.vcf.gz >\
  $tmp_dir/true_${SAMPLE_INDEX}.vcf

#whatshap compare --names truth,vote --tsv-pairwise $out_dir/error_vote_${SAMPLE_INDEX}.tsv $tmp_dir/true_${SAMPLE_INDEX}.vcf $tmp_dir/vote_${SAMPLE_INDEX}.vcf

# vcftools
vcftools --vcf $tmp_dir/vote_${SAMPLE_INDEX}.vcf \
  --diff $tmp_dir/true_${SAMPLE_INDEX}.vcf \
  --diff-switch-error \
  --out ${out_dir}/switches/vote_${SAMPLE_INDEX}

rm $tmp_dir/true_${SAMPLE_INDEX}.vcf
rm $tmp_dir/eagle_${SAMPLE_INDEX}.vcf.gz
rm $tmp_dir/beagle_${SAMPLE_INDEX}.vcf.gz
rm $tmp_dir/shapeit_${SAMPLE_INDEX}.vcf.gz
rm $tmp_dir/common_pos_${SAMPLE_INDEX}.txt

bcftools view -Ob $tmp_dir/vote_${SAMPLE_INDEX}.vcf > $out_dir/vote_${SAMPLE_INDEX}.bcf
rm $tmp_dir/vote_${SAMPLE_INDEX}.vcf

else
echo "All good for ${SAMPLE_INDEX}"
fi
