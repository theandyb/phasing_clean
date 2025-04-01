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

beagle_map="${base_dir}/data/ref_maps/beagle/plink.chrX.GRCh38.map"
beagle_jar="/net/snowwhite/home/beckandy/bin/beagle.05May22.33a.jar"

eagle_map="/net/snowwhite/home/beckandy/software/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz"

out_dir="/net/snowwhite/home/beckandy/research/phasing_clean/output/switch_errors"
# Read line from sample_pairs.csv
line=$(awk -v id=${SLURM_ARRAY_TASK_ID} 'NR==id{ print; exit }' /net/snowwhite/home/beckandy/research/phasing_clean/data/sample_pairs.csv)
arrSub=(${line//,/ })
sub_a=${arrSub[1]}
sub_b=${arrSub[2]}

echo "Our subjects are $sub_a and $sub_b"

# mkdir -p $out_dir
# mkdir -p $out_dir/site_list
# mkdir -p $out_dir/phase_results
# mkdir -p $out_dir/vcf_n_sites
# mkdir -p $out_dir/het_loc/annotated
# mkdir -p $out_dir/switch_errors/eagle/annotated
# mkdir -p $out_dir/switch_errors/beagle/annotated
# mkdir -p $out_dir/switch_errors/shapeit/annotated
# mkdir -p $out_dir/whatshap/beagle
# mkdir -p $out_dir/whatshap/eagle
# mkdir -p $out_dir/whatshap/shapeit
# mkdir -p $out_dir/vcf/beagle
# mkdir -p $out_dir/vcf/eagle
# mkdir -p $out_dir/vcf/shapeit


#if [ ! -f "${out_dir}/switch_errors/eagle/error_${SLURM_ARRAY_TASK_ID}.diff.switch"  ]; then
if [ ! -f "$out_dir/phase_results/shapeit_${SLURM_ARRAY_TASK_ID}.txt" ]; then

# Get list of sites for our two subjects
bcftools view -v snps -Ou $VCF | \
  bcftools view -s $sub_a,$sub_b -Ou |\
  bcftools query  -f '%CHROM\t%POS\n' | tr -s ' ' |
  awk '!seen[$0]++' > $out_dir/site_list/sites_${SLURM_ARRAY_TASK_ID}.tsv
echo -e "Site list generated"

#Filter VCF to our subjects and their list of sites, removing positions with not alt in either as well
bcftools view -v snps -s $sub_a,$sub_b -Ou $VCF |\
	bcftools view -T $out_dir/site_list/sites_${SLURM_ARRAY_TASK_ID}.tsv |\
	bcftools view -Oz > /tmp/andy_input_${SLURM_ARRAY_TASK_ID}.vcf.gz

# Generate pseudodiploid from two haplotypes
Rscript /net/snowwhite/home/beckandy/research/phasing_clean/code/diploidifier.R /tmp/andy_input_${SLURM_ARRAY_TASK_ID}.vcf.gz /tmp/andy_input_${SLURM_ARRAY_TASK_ID} ${SLURM_ARRAY_TASK_ID}

bcftools view -Oz /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf > /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz
bcftools view -Oz /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf > /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf.gz
rm /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf
bcftools index /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz


# Save number of heterozygotes
n_hets=$(bcftools query -f '[%GT]\n' /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz | awk '{if($1 == "1/0" || $1 == "0/1")count++}END{print(count)}')
echo -e "${SLURM_ARRAY_TASK_ID}\t$n_hets" > $out_dir/vcf_n_sites/pair_${SLURM_ARRAY_TASK_ID}_hets.txt
echo -e "${SLURM_ARRAY_TASK_ID}\t$n_hets"

# Get list of heterozygouse sites
bcftools query -f '%CHROM\t%POS\t[%GT\t]\n' /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz |\
awk '{if($3 == "1/0" || $3 == "0/1")print($0)}' > $out_dir/het_loc/pair_${SLURM_ARRAY_TASK_ID}_het_loc.txt

# Reference VCF
bcftools view -Ou $VCF | \
  bcftools view -s^$sub_a,$sub_b -Ou | \
  bcftools norm --no-version -Ou -m -any | \
  bcftools norm --no-version -Ob -d none -f /net/snowwhite/home/beckandy/research/phasing_clean/data/ref_GRCh38.fna > /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf
bcftools index /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf

n_ref=$(bcftools query -f '%POS\n' /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf | wc -l)
echo -e "${SLURM_ARRAY_TASK_ID}\t$n_ref" > $out_dir/vcf_n_sites/pair_${SLURM_ARRAY_TASK_ID}_ref.txt

### PHASING
eagle --vcfTarget /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz \
  --vcfRef /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf \
  --geneticMapFile=$eagle_map \
  --vcfOutFormat v \
  --chrom chrX \
  --numThreads 4 \
  --outPrefix=/tmp/andy_eagle_${SLURM_ARRAY_TASK_ID}

echo "EAGLE done"

$shapeit_dir/phase_common/bin/phase_common \
  --input /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz \
  --map $shapeit_map \
  --reference /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf \
  --region chrX \
  --thread 4 \
  --output /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.bcf

echo "SHAPEIT done"

# BEAGLE

bcftools view /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf -Ov > /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.vcf

java -Xmx7g -jar $beagle_jar \
  gt=/tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz \
  ref=/tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.vcf \
  map=$beagle_map \
  nthreads=4 \
  impute=false \
  out=/tmp/andy_beagle_${SLURM_ARRAY_TASK_ID}

# write phasing results in tsv form
bcftools query -f '%CHROM\t%POS\t[%GT]\n' /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.bcf > $out_dir/phase_results/shapeit_${SLURM_ARRAY_TASK_ID}.txt
bcftools query -f '%CHROM\t%POS\t[%GT]\n' /tmp/andy_eagle_${SLURM_ARRAY_TASK_ID}.vcf > $out_dir/phase_results/eagle_${SLURM_ARRAY_TASK_ID}.txt
bcftools query -f '%CHROM\t%POS\t[%GT]\n' /tmp/andy_beagle_${SLURM_ARRAY_TASK_ID}.vcf.gz > $out_dir/phase_results/beagle_${SLURM_ARRAY_TASK_ID}.txt

# Generate switch error files

bcftools view -Oz /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.bcf > /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.vcf.gz
rm /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.bcf*

vcftools --vcf /tmp/andy_eagle_${SLURM_ARRAY_TASK_ID}.vcf \
  --gzdiff /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf.gz \
  --diff-switch-error \
  --out $out_dir/switch_errors/eagle/error_${SLURM_ARRAY_TASK_ID}

vcftools --gzvcf /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.vcf.gz \
  --gzdiff /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf.gz \
  --diff-switch-error \
  --out $out_dir/switch_errors/shapeit/error_${SLURM_ARRAY_TASK_ID}

vcftools --gzvcf /tmp/andy_beagle_${SLURM_ARRAY_TASK_ID}.vcf.gz \
  --gzdiff /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf.gz \
  --diff-switch-error \
  --out  $out_dir/switch_errors/beagle/error_${SLURM_ARRAY_TASK_ID}

# Count number of heterozygous sites within each region
# ## Generate "bed"-like file with sites
# bcftools query -i 'GT="1/0" || GT="0/1"' -f '%CHROM\t%POS\t%POS\n' /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz  |\
#   tr -s ' ' |\
#   awk '{print($1"\t"$2-1"\t"$3)}' > /tmp/andy_sites_${SLURM_ARRAY_TASK_ID}.bed

# ## bedtools intersect
# tail -n+2 $out_dir/switch_errors/eagle/error_${SLURM_ARRAY_TASK_ID}.diff.switch |\
#   bedtools intersect -a - -b /tmp/andy_sites_${SLURM_ARRAY_TASK_ID}.bed -c > $out_dir/switch_errors/eagle/hets_per_region_${SLURM_ARRAY_TASK_ID}.bed
#
# tail -n+2 $out_dir/switch_errors/shapeit/error_${SLURM_ARRAY_TASK_ID}.diff.switch |\
#   bedtools intersect -a - -b /tmp/andy_sites_${SLURM_ARRAY_TASK_ID}.bed -c > $out_dir/switch_errors/shapeit/hets_per_region_${SLURM_ARRAY_TASK_ID}.bed
#
# tail -n+2 $out_dir/switch_errors/beagle/error_${SLURM_ARRAY_TASK_ID}.diff.switch |\
#   bedtools intersect -a - -b /tmp/andy_sites_${SLURM_ARRAY_TASK_ID}.bed -c > $out_dir/switch_errors/beagle/hets_per_region_${SLURM_ARRAY_TASK_ID}.bed

# whatshap adds additional summary statistics regarding phase errors
# need gunzipped truth
gunzip /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf.gz
whatshap compare --names truth,eagle --tsv-pairwise $out_dir/whatshap/eagle/eval_${SLURM_ARRAY_TASK_ID}.tsv /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf /tmp/andy_eagle_${SLURM_ARRAY_TASK_ID}.vcf
# need gunzipped shapeit
gunzip /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.vcf.gz
whatshap compare --names truth,shapeit --tsv-pairwise $out_dir/whatshap/shapeit/eval_${SLURM_ARRAY_TASK_ID}.tsv /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf  /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.vcf

gunzip /tmp/andy_beagle_${SLURM_ARRAY_TASK_ID}.vcf.gz
whatshap compare --names truth,beagle --tsv-pairwise $out_dir/whatshap/beagle/eval_${SLURM_ARRAY_TASK_ID}.tsv /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf  /tmp/andy_beagle_${SLURM_ARRAY_TASK_ID}.vcf

# stash truth vcf
bcftools view -Ob /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf > $out_dir/vcf/pair_${SLURM_ARRAY_TASK_ID}_true.bcf

# stash test vcf
bcftools view -Ob /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz > $out_dir/vcf/pair_${SLURM_ARRAY_TASK_ID}_test.bcf

# stash phasing results
bcftools view -Ob /tmp/andy_eagle_${SLURM_ARRAY_TASK_ID}.vcf > $out_dir/vcf/eagle/pair_${SLURM_ARRAY_TASK_ID}.bcf
bcftools view -Ob /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.vcf > $out_dir/vcf/shapeit/pair_${SLURM_ARRAY_TASK_ID}.bcf

sed -i '5s/^/##contig=<ID=chrX>\n/' /tmp/andy_beagle_${SLURM_ARRAY_TASK_ID}.vcf
bcftools view -Ob /tmp/andy_beagle_${SLURM_ARRAY_TASK_ID}.vcf > $out_dir/vcf/beagle/pair_${SLURM_ARRAY_TASK_ID}.bcf

# clean up tmp files
rm /tmp/andy_beagle_${SLURM_ARRAY_TASK_ID}.log
rm /tmp/andy_beagle_${SLURM_ARRAY_TASK_ID}.vcf
rm /tmp/andy_eagle_${SLURM_ARRAY_TASK_ID}.vcf
rm /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test*
rm /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf
rm /tmp/andy_input_${SLURM_ARRAY_TASK_ID}.vcf.gz
rm /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf*
rm /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.vcf*
rm /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.vcf
#rm /tmp/andy_sites_${SLURM_ARRAY_TASK_ID}.bed

else
echo "All good for ${SLURM_ARRAY_TASK_ID}"
fi
