#!/bin/bash
#
#SBATCH --job-name=trioPhase
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8GB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --array=1-602
#SBATCH --constraint=avx2
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing_clean/output/trio_phase_1/slurm/sample.%A.%a.err
#SBATCH --output=/net/snowwhite/home/beckandy/research/phasing_clean/output/trio_phase_1/slurm/sample.%A.%a.out

# Code for phasing 602 child samples against a panel with their parents removed
# also remove triple heterozygous sites

chrom=15
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
exclude_samples="${exclude_dir}/sample_${SLURM_ARRAY_TASK_ID}.txt"
target_sample=$(head -n ${SLURM_ARRAY_TASK_ID} /net/snowwhite/home/beckandy/research/phasing_clean/data/1kgp/child_ids.txt  | tail -1)
working_dir="/net/snowwhite/home/beckandy/scratch/"
triple_het_file="${base_dir}/data/1kgp/chr${chrom}/triple_het/sample_${SLURM_ARRAY_TASK_ID}.tsv"

echo "Phasing subject ${target_sample} on chromosome ${chrom}"
echo "Beagle map: ${beagle_map}"
echo "SHAPEIT map: ${shapeit_map}"

if ! test -f $out_dir/truth/sample_${SLURM_ARRAY_TASK_ID}.vcf.gz; then
  # generate reference panel for phasing
  if ! test -f ${working_dir}andy_ref_${SLURM_ARRAY_TASK_ID}.vcf.gz; then
  bcftools view  -S ^${exclude_samples} -T ^${triple_het_file} -Oz $panel_vcf  > ${working_dir}andy_ref_${SLURM_ARRAY_TASK_ID}_noTH.vcf.gz #| \
    #bcftools norm --no-version -m -any | \
    #bcftools norm --no-version -Oz -d none -f /net/snowwhite/home/beckandy/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna> ${working_dir}andy_ref_${SLURM_ARRAY_TASK_ID}.vcf.gz
  bcftools index ${working_dir}andy_ref_${SLURM_ARRAY_TASK_ID}_noTH.vcf.gz
  fi

  echo "reference panel generated"

  # generate target vcf.gz for phasing
  if ! test -f ${working_dir}andy_target_${SLURM_ARRAY_TASK_ID}.vcf.gz; then
  bcftools view  -Oz -s $target_sample -T ^${triple_het_file} $source_vcf > ${working_dir}andy_target_${SLURM_ARRAY_TASK_ID}_noTH.vcf.gz
  fi
  echo "target vcf generated"

  # generate "truth" vcf.gz for comparison
  if ! test -f ${working_dir}andy_truth_${SLURM_ARRAY_TASK_ID}.vcf.gz; then
  bcftools view  -Oz -s $target_sample $phased_vcf -T ^${triple_het_file} > ${working_dir}andy_truth_${SLURM_ARRAY_TASK_ID}_noTH.vcf.gz
  fi
  echo "truth vcf generated"

  # bcf  version of target and reference
  if ! test -f ${working_dir}andy_ref_${SLURM_ARRAY_TASK_ID}.bcf; then
  bcftools view -Ob ${working_dir}andy_target_${SLURM_ARRAY_TASK_ID}_noTH.vcf.gz > ${working_dir}andy_target_${SLURM_ARRAY_TASK_ID}_noTH.bcf
  bcftools view -Ob ${working_dir}andy_ref_${SLURM_ARRAY_TASK_ID}_noTH.vcf.gz > ${working_dir}andy_ref_${SLURM_ARRAY_TASK_ID}_noTH.bcf
  bcftools index  ${working_dir}andy_target_${SLURM_ARRAY_TASK_ID}_noTH.bcf
  bcftools index ${working_dir}andy_ref_${SLURM_ARRAY_TASK_ID}_noTH.bcf
  fi

  echo "bcf target and reference"

  ## Phasing!
  eagle --vcfTarget ${working_dir}andy_target_${SLURM_ARRAY_TASK_ID}_noTH.bcf \
    --vcfRef ${working_dir}andy_ref_${SLURM_ARRAY_TASK_ID}_noTH.bcf \
    --geneticMapFile=${eagle_map} \
    --vcfOutFormat v \
    --chrom "chr${chrom}" \
    --numThreads 4 \
    --outPrefix=${working_dir}andy_eagle_${SLURM_ARRAY_TASK_ID}_noTH

  $shapeit_dir/phase_common/bin/phase_common \
    --input ${working_dir}andy_target_${SLURM_ARRAY_TASK_ID}_noTH.bcf \
    --map ${shapeit_map} \
    --reference ${working_dir}andy_ref_${SLURM_ARRAY_TASK_ID}_noTH.bcf \
    --region "chr${chrom}" \
    --thread 4 \
    --output ${working_dir}andy_shapeit_${SLURM_ARRAY_TASK_ID}_noTH.bcf

  java -Xmx8g -jar $beagle_jar \
    gt=${working_dir}andy_target_${SLURM_ARRAY_TASK_ID}_noTH.vcf.gz \
    ref=${working_dir}andy_ref_${SLURM_ARRAY_TASK_ID}_noTH.vcf.gz \
    map=${beagle_map} \
    nthreads=4 \
    impute=false \
    out=${working_dir}andy_beagle_${SLURM_ARRAY_TASK_ID}_noTH

  # move results to output directory
  mv ${working_dir}andy_beagle_${SLURM_ARRAY_TASK_ID}_noTH.vcf.gz $out_dir/beagle/sample_${SLURM_ARRAY_TASK_ID}.vcf.gz
  bcftools view -Oz ${working_dir}andy_eagle_${SLURM_ARRAY_TASK_ID}_noTH.vcf > $out_dir/eagle/sample_${SLURM_ARRAY_TASK_ID}.vcf.gz
  bcftools view -Oz ${working_dir}andy_shapeit_${SLURM_ARRAY_TASK_ID}_noTH.bcf > $out_dir/shapeit/sample_${SLURM_ARRAY_TASK_ID}.vcf.gz
  mv ${working_dir}andy_truth_${SLURM_ARRAY_TASK_ID}_noTH.vcf.gz $out_dir/truth/sample_${SLURM_ARRAY_TASK_ID}.vcf.gz

  # clean up!
  rm ${working_dir}andy_shapeit_${SLURM_ARRAY_TASK_ID}_noTH.*
  rm ${working_dir}andy_ref_${SLURM_ARRAY_TASK_ID}_noTH.*
  rm ${working_dir}andy_target_${SLURM_ARRAY_TASK_ID}_noTH.*
  rm ${working_dir}andy_beagle_${SLURM_ARRAY_TASK_ID}_noTH.*
  rm ${working_dir}andy_eagle_${SLURM_ARRAY_TASK_ID}_noTH.*

  echo "complete!"
fi

