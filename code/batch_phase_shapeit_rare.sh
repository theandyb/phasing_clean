#!/bin/bash
#
#SBATCH --job-name=pdPhase
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8GB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --array=1-1000
#SBATCH --constraint=avx2
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing_clean/output/switch_errors/slurm/pseudo_dip.%A.%a.err
#SBATCH --output=/net/snowwhite/home/beckandy/research/phasing_clean/output/switch_errors/slurm/pseudo_dip.%A.%a.out

# Code for generating switch error results for pseudo-diploids constructed using male X chromosomes (non-PAR)
# This version uses all samples not used in the generation of the PD as the reference panel (males and females)
# Here we do not do any downsampling to make reference panels of a comparable size to male-only reference panels

base_dir="/net/snowwhite/home/beckandy/research/phasing_clean"
VCF="${base_dir}/data/1kgp/chrX_2504_snps_noPAR_noSing.bcf"
shapeit_dir="/net/snowwhite/home/beckandy/software/shapeit5/shapeit5"
shapeit_map="${shapeit_dir}/resources/maps/b38/chrX.b38.gmap.gz"
chunk_file="${base_dir}/data/glimpse_chunks/chrX_chunks.txt"

out_dir="/net/snowwhite/home/beckandy/research/phasing_clean/output/switch_errors"

# Read line from sample_pairs.csv
line=$(awk -v id=${SLURM_ARRAY_TASK_ID} 'NR==id{ print; exit }' /net/snowwhite/home/beckandy/research/phasing_clean/data/sample_pairs.csv)
arrSub=(${line//,/ })
sub_a=${arrSub[1]}
sub_b=${arrSub[2]}

echo "Our subjects are $sub_a and $sub_b"

#Filter VCF to our subjects and their list of sites, removing positions with not alt in either as well
bcftools view -I -v snps -s $sub_a,$sub_b -Ou $VCF |\
	bcftools view -I -T $out_dir/site_list/sites_${SLURM_ARRAY_TASK_ID}.tsv |\
	bcftools view -I -Oz > /tmp/andy_input_${SLURM_ARRAY_TASK_ID}.vcf.gz

# Generate pseudodiploid from two haplotypes
Rscript /net/snowwhite/home/beckandy/research/phasing_clean/code/diploidifier.R /tmp/andy_input_${SLURM_ARRAY_TASK_ID}.vcf.gz /tmp/andy_input_${SLURM_ARRAY_TASK_ID} ${SLURM_ARRAY_TASK_ID}

bcftools view -I -Ob /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf > /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.bcf
bcftools view -I -Ob /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf > /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.bcf
rm /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf
bcftools index /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.bcf
bcftools index /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.bcf

# Reference VCF
bcftools view -I -Ou $VCF | \
  bcftools view -I -s^$sub_a,$sub_b -Ou | \
  bcftools norm --no-version -Ou -m -any | \
  bcftools norm --no-version -Ob -d none -f /net/snowwhite/home/beckandy/research/phasing_clean/data/ref_GRCh38.fna > /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf
bcftools index /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf


$shapeit_dir/phase_common/bin/phase_common \
  --input /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.bcf \
  --map $shapeit_map \
  --reference /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf \
  --region chrX \
  --filter-maf 0.001 \
  --thread 20 \
  --output /tmp/andy_scaffold_${SLURM_ARRAY_TASK_ID}.bcf


# Get first and last position from input VCF for chunking
first_pos=$(bcftools query -f '%POS\n' /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.bcf | head -1)
last_pos=$(bcftools query -f '%POS\n' /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.bcf | tail -1)
echo "First position: $first_pos"
echo "Last position: $last_pos"

# Get same from scaffolded output
scaf_first_pos=$(bcftools query -f '%POS\n' /tmp/andy_scaffold_${SLURM_ARRAY_TASK_ID}.bcf | head -1)
scaf_last_pos=$(bcftools query -f '%POS\n' /tmp/andy_scaffold_${SLURM_ARRAY_TASK_ID}.bcf | tail -1)
echo "Scaffolded first position: $scaf_first_pos"
echo "Scaffolded last position: $scaf_last_pos"

gdb --args $shapeit_dir/phase_rare/bin/phase_rare \
  --input /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.bcf \
  --scaffold /tmp/andy_scaffold_${SLURM_ARRAY_TASK_ID}.bcf \
  --map $shapeit_map \
  --input-region "chrX:2781513-155700570" \
  --scaffold-region "chrX:2781513-155700570" \
  --thread 20 \
  --output /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.bcf

gdb --args /net/snowwhite/home/beckandy/software/shapeit5/shapeit5_old/phase_rare/bin/phase_rare \
  --input /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.bcf \
  --scaffold /tmp/andy_scaffold_${SLURM_ARRAY_TASK_ID}.bcf \
  --map $shapeit_map \
  --input-region "chrX:2781513-155700570" \
  --scaffold-region "chrX:2781513-155700570" \
  --thread 20 \
  --output /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.bcf

# write phasing results in tsv form
bcftools query -f '%CHROM\t%POS\t[%GT]\n' /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.bcf > $out_dir/phase_results/shapeit_${SLURM_ARRAY_TASK_ID}.txt
# Generate switch error files

bcftools view -Oz /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.bcf > /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.vcf.gz
rm /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.bcf*


vcftools --gzvcf /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.vcf.gz \
  --gzdiff /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf.gz \
  --diff-switch-error \
  --out $out_dir/switch_errors/shapeit/error_${SLURM_ARRAY_TASK_ID}


