#!/bin/bash
#
#SBATCH --job-name=noTH
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8GB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --array=0-13243
#SBATCH --constraint=avx2
#SBATCH --exclude=r[6320,6319,6333,6331,6323]
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing_clean/output/slurm/sample.%A.%a.err
#SBATCH --output=/net/snowwhite/home/beckandy/research/phasing_clean/output/slurm/sample.%A.%a.out

# Calculate the chromosome and job index based on the array index
total_jobs_per_chromosome=602
chromosome_index=$(( $SLURM_ARRAY_TASK_ID / total_jobs_per_chromosome ))
job_index_within_chromosome=$(( $SLURM_ARRAY_TASK_ID % total_jobs_per_chromosome ))

chrom=$(( chromosome_index + 1 ))
sample=$(( job_index_within_chromosome + 1 ))

echo "Running job $sample for chromosome $chrom"

# Code for phasing 602 child samples against a panel with their parents removed
# also remove triple heterozygous sites

base_dir="/net/snowwhite/home/beckandy/research/phasing_clean"
shapeit_dir="/net/snowwhite/home/beckandy/software/shapeit5/shapeit5" #location where shapeit is compiled
shapeit_map="${shapeit_dir}/resources/maps/b38/chr${chrom}.b38.gmap.gz"
chunk_file="${shapeit_dir}/resources/chunks/b38/4cM/chunks_chr${chrom}.txt"


panel_vcf="${base_dir}/data/1kgp/chr${chrom}/chr${chrom}_phased_overlap_2504.bcf"
phased_vcf="${base_dir}/data/1kgp/chr${chrom}/chr${chrom}_phased_overlap.bcf"
source_vcf="${base_dir}/data/1kgp/chr${chrom}/chr${chrom}_unphased_overlap.bcf"
exclude_dir="${base_dir}/data/1kgp/exclude_samples/"
out_dir="${base_dir}/output/trio_phase_${chrom}/no_th/"
exclude_samples="${exclude_dir}/sample_${sample}.txt"
target_sample=$(head -n ${sample} /net/snowwhite/home/beckandy/research/phasing_clean/data/1kgp/child_ids.txt  | tail -1)
working_dir="/tmp/"
triple_het_file="${base_dir}/data/1kgp/chr${chrom}/triple_het/sample_${sample}.tsv"

echo "Phasing subject ${target_sample} on chromosome ${chrom}"

test_file1="$out_dir/truth/sample_${sample}.vcf.gz"
test_file2="$out_dir/shapeit/sample_${sample}_rare.vcf.gz"

bcftools view -I -S ^${exclude_samples} -T ^${triple_het_file} -Ob $panel_vcf  > ${working_dir}andy_ref_${sample}_chr_${chrom}.bcf #| \
bcftools index ${working_dir}andy_ref_${sample}_chr_${chrom}.bcf

echo "reference panel generated"

bcftools view  -Ob -I -s $target_sample -T ^${triple_het_file} $source_vcf > ${working_dir}andy_target_${sample}_chr_${chrom}.bcf
bcftools index -f ${working_dir}andy_target_${sample}_chr_${chrom}.bcf

# generate "truth" vcf.gz for comparison
bcftools view -I -Ob -s $target_sample $phased_vcf -T ^${triple_het_file} > ${working_dir}andy_truth_${sample}_chr_${chrom}.bcf
bcftools index  ${working_dir}andy_truth_${sample}_chr_${chrom}.bcf


echo "bcf target and reference"

$shapeit_dir/phase_common/bin/phase_common \
    --input ${working_dir}andy_target_${sample}_chr_${chrom}.bcf \
    --map ${shapeit_map} \
    --reference ${working_dir}andy_ref_${sample}_chr_${chrom}.bcf \
    --region "chr${chrom}" \
    --filter-maf 0.001 \
    --thread 20 \
    --output ${working_dir}andy_scaffold_${sample}_chr_${chrom}.bcf


shapeit_dir="/net/snowwhite/home/beckandy/software/shapeit5/shapeit5_old" #location where shapeit is compiled
chunk_file="${shapeit_dir}/resources/chunks/b38/4cM/chunks_chr22.txt"
while read LINE; do
    CHK=$(echo $LINE | awk '{ print $1; }')
    SRG=$(echo $LINE | awk '{ print $3; }')
    IRG=$(echo $LINE | awk '{ print $4; }')
    $shapeit_dir/phase_rare/bin/phase_rare \
      --input ${working_dir}andy_target_${sample}_chr_${chrom}.bcf \
      --scaffold ${working_dir}andy_scaffold_${sample}_chr_${chrom}.bcf \
      --map ${shapeit_map} \
      --input-region ${IRG} \
      --scaffold-region ${SRG} \
      --thread 8 \
      --output ${out_dir}/shapeit/sample_${sample}_chunk_${CHK}.bcf
done < ${chunk_file}

gdb --args $shapeit_dir/phase_rare/bin/phase_rare \
    --input ${working_dir}andy_target_${sample}_chr_${chrom}.bcf \
    --scaffold ${working_dir}andy_scaffold_${sample}_chr_${chrom}.bcf \
    --map ${shapeit_map} \
    --input-region ${IRG} \
    --scaffold-region ${SRG} \
    --thread 1 \
    --output ${out_dir}/shapeit/sample_${sample}_chunk_${CHK}.bcf

