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
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing_clean/output/slurm/sample_noTH.%A.%a.err
#SBATCH --output=/net/snowwhite/home/beckandy/research/phasing_clean/output/slurm/sample_noTH.%A.%a.out

# Calculate the chromosome and job index based on the array index
total_jobs_per_chromosome=602
chromosome_index=$(( $SLURM_ARRAY_TASK_ID / total_jobs_per_chromosome ))
job_index_within_chromosome=$(( $SLURM_ARRAY_TASK_ID % total_jobs_per_chromosome ))

chrom=$(( chromosome_index + 1 ))
sample=$(( job_index_within_chromosome + 1 ))

echo "Running job $sample for chromosome $chrom"

# Code for phasing 602 child samples against a panel with their parents removed
# also remove triple heterozygous sites

eagle_map="/net/snowwhite/home/beckandy/software/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz"
base_dir="/net/snowwhite/home/beckandy/research/phasing_clean"
shapeit_dir="/net/snowwhite/home/beckandy/software/shapeit5/shapeit5" #location where shapeit is compiled
shapeit_map="${shapeit_dir}/resources/maps/b38/chr${chrom}.b38.gmap.gz"
beagle_map="${base_dir}/data/ref_maps/beagle/plink.chr${chrom}.GRCh38.map"
beagle_jar="/net/snowwhite/home/beckandy/bin/beagle.05May22.33a.jar" # location of beagle executive

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
echo "Beagle map: ${beagle_map}"
echo "SHAPEIT map: ${shapeit_map}"

test_file1="$out_dir/truth/sample_${sample}.vcf.gz"
test_file2="$out_dir/shapeit/sample_${sample}.vcf.gz"
test_file3="$out_dir/beagle/sample_${sample}.vcf.gz"
test_file4="$out_dir/eagle/sample_${sample}.vcf.gz"


#if ! test -f $out_dir/truth/sample_${sample}.vcf.gz; then
if [ ! -e "$test_file1" ] || [ ! -e "$test_file2" ] || [ ! -e "$test_file3" ] || [ ! -e "$test_file4" ]; then
  # generate reference panel for phasing
  if ! test -f ${working_dir}andy_ref_${sample}.vcf.gz; then
  bcftools view  -S ^${exclude_samples} -T ^${triple_het_file} -Oz $panel_vcf  > ${working_dir}andy_ref_${sample}_noTH.vcf.gz #| \
  bcftools index ${working_dir}andy_ref_${sample}_noTH.vcf.gz
  fi

  echo "reference panel generated"

  # generate target vcf.gz for phasing
  if ! test -f ${working_dir}andy_target_${sample}.vcf.gz; then
  bcftools view  -Oz -s $target_sample -T ^${triple_het_file} $source_vcf > ${working_dir}andy_target_${sample}_noTH.vcf.gz
  fi
  echo "target vcf generated"

  # generate "truth" vcf.gz for comparison
  if ! test -f ${working_dir}andy_truth_${sample}.vcf.gz; then
  bcftools view  -Oz -s $target_sample $phased_vcf -T ^${triple_het_file} > ${working_dir}andy_truth_${sample}_noTH.vcf.gz
  fi
  echo "truth vcf generated"

  # bcf  version of target and reference
  if ! test -f ${working_dir}andy_ref_${sample}.bcf; then
  bcftools view -Ob ${working_dir}andy_target_${sample}_noTH.vcf.gz > ${working_dir}andy_target_${sample}_noTH.bcf
  bcftools view -Ob ${working_dir}andy_ref_${sample}_noTH.vcf.gz > ${working_dir}andy_ref_${sample}_noTH.bcf
  bcftools index  ${working_dir}andy_target_${sample}_noTH.bcf
  bcftools index ${working_dir}andy_ref_${sample}_noTH.bcf
  fi

  echo "bcf target and reference"

  ## Phasing!
  eagle --vcfTarget ${working_dir}andy_target_${sample}_noTH.bcf \
    --vcfRef ${working_dir}andy_ref_${sample}_noTH.bcf \
    --geneticMapFile=${eagle_map} \
    --vcfOutFormat v \
    --chrom "chr${chrom}" \
    --numThreads 4 \
    --outPrefix=${working_dir}andy_eagle_${sample}_noTH

  $shapeit_dir/phase_common/bin/phase_common \
    --input ${working_dir}andy_target_${sample}_noTH.bcf \
    --map ${shapeit_map} \
    --reference ${working_dir}andy_ref_${sample}_noTH.bcf \
    --region "chr${chrom}" \
    --thread 4 \
    --output ${working_dir}andy_shapeit_${sample}_noTH.bcf

  java -Xmx8g -jar $beagle_jar \
    gt=${working_dir}andy_target_${sample}_noTH.vcf.gz \
    ref=${working_dir}andy_ref_${sample}_noTH.vcf.gz \
    map=${beagle_map} \
    nthreads=4 \
    impute=false \
    out=${working_dir}andy_beagle_${sample}_noTH

  # move results to output directory
  mv ${working_dir}andy_beagle_${sample}_noTH.vcf.gz $out_dir/beagle/sample_${sample}.vcf.gz
  bcftools view -Oz ${working_dir}andy_eagle_${sample}_noTH.vcf > $out_dir/eagle/sample_${sample}.vcf.gz
  bcftools view -Oz ${working_dir}andy_shapeit_${sample}_noTH.bcf > $out_dir/shapeit/sample_${sample}.vcf.gz
  mv ${working_dir}andy_truth_${sample}_noTH.vcf.gz $out_dir/truth/sample_${sample}.vcf.gz

  # clean up!
  rm ${working_dir}andy_shapeit_${sample}_noTH.*
  rm ${working_dir}andy_ref_${sample}_noTH.*
  rm ${working_dir}andy_target_${sample}_noTH.*
  rm ${working_dir}andy_beagle_${sample}_noTH.*
  rm ${working_dir}andy_eagle_${sample}_noTH.*

  echo "complete!"
else
  echo "Already done!"
fi

