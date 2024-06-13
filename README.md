# Phasing Algorithm Evaluation with Synthetic Diploids

This repository is an accompaniment to the paper ""

## Data included in this repository

-   `data/1kgp/unrelated_subjects.txt`: A list of IDs for the 2,502 unrelated samples from 1kGP phase 3
-   `data/1kgp/unrelated_subj.tsv`: Slightly processed version of the [Tab delimited file list for the 2504 panel](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index)
-   `data/1kgp/subject_info.csv`: SAMPLE_NAME, POPULATION, SUPER, fatherID, motherID, and sex (1 - male, 2 - female)

## Misc Reference Data

Chromosome sizes and windows:

```         
curl https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes > data/hg38.chrom.sizes
cat data/hg38.chrom.sizes | grep 'chrX\s' > data/chrX.sizes
cat data/hg38.chrom.sizes | grep 'chr8\s' > data/chr8.sizes
cat data/hg38.chrom.sizes | grep 'chr15\s' > data/chr15.sizes
cat data/hg38.chrom.sizes | grep 'chr22\s' > data/chr22.sizes
cat data/hg38.chrom.sizes | grep 'chr1\s' > data/chr1.sizes

bedtools makewindows -g data/chrX.sizes -w 10000 > chrX_10k.bed
bedtools makewindows -g data/chr8.sizes -w 10000 > chr8_10k.bed
bedtools makewindows -g data/chr15.sizes -w 10000 > chr15_10k.bed
bedtools makewindows -g data/chr22.sizes -w 10000 > chr22_10k.bed
bedtools makewindows -g data/chr1.sizes -w 10000 > chr1_10k.bed
```

DeCODE sex averaged genetic maps

```
bigWigToBedGraph -chrom=chrX http://hgdownload.soe.ucsc.edu/gbdb/hg38/recombRate/recombAvg.bw stdout > data/decode/chrX_recomb.bed

bigWigToBedGraph -chrom=chr8 http://hgdownload.soe.ucsc.edu/gbdb/hg38/recombRate/recombAvg.bw stdout > data/decode/chr8_recomb.bed

bigWigToBedGraph -chrom=chr15 http://hgdownload.soe.ucsc.edu/gbdb/hg38/recombRate/recombAvg.bw stdout > data/decode/chr15_recomb.bed

bigWigToBedGraph -chrom=chr1 http://hgdownload.soe.ucsc.edu/gbdb/hg38/recombRate/recombAvg.bw stdout > data/decode/chr1_recomb.bed
```

1kgp pilot mask:

```
curl ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/PilotMask/20160622.allChr.pilot_mask.bed | grep "chrX" > data/1kgp/chrX_pilot_mask.bed

curl ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/PilotMask/20160622.allChr.pilot_mask.bed | grep "chr1\s" > data/1kgp/chr1/chr1_pilot_mask.bed

curl ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/PilotMask/20160622.allChr.pilot_mask.bed | grep "chr15" > data/1kgp/chr15/chr15_pilot_mask.bed
```

### Genetic maps for phasing

#### Beagle

For this map, we need to append the string 'chr' to the beginning of each row to match the chromosome naming convention of the VCF and fasta files. Note: these come from HapMap

```
mkdir -p data/ref_maps/beagle
wget https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
unzip plink.GRCh38.map.zip
rm README.txt plink.GRCh38.map.zip

mv plink*.GRCh38.map data/ref_maps/beagle/
sed -i -e 's/^/chr/' data/ref_maps/beagle/plink*.GRCh38.map
```

#### Eagle

The eagle software includes the map in a directory named `tables`. It also is from HapMap, but has been reformatted.

#### SHAPEIT

SHAPEIT5 provides a map in the directory under which it is installed in the subdirectory `resources/maps/b38`


## Downloading and Processing X Chromosome Reference Data

To generate our synthetic diploids, we use start with the [1000 Genomes Project 30x on GRCh38 Phased VCF](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/) for chromosome X:

```         
curl ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz > data/1kgp/chrX.vcf.gz
curl ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz.tbi > data/1kgp/chrX.vcf.gz.tbi
```

We apply a few filtering steps to this VCF:

1.  We subset down to the 2,504 unrelated 1kGP phase 3 samples
2.  We remove the pseudoautosomal region (chrX:10001-2781479 and chrX:155701383-156030895)
3.  We apply the 1kGP pilot genome accessibility mask to remove hard to sequence regions of the chromosome
4.  We only keep snps
5.  We remove excess annotation from the VCF
6.  We also remove multiallelic sites

We download the [1kGP Pilot accessibility map](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/PilotMask/) (note that this bed file is a list of sites that we want to keep in the VCF):

```         
curl ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/PilotMask/20160622.allChr.pilot_mask.bed | grep "chrX" > data/1kgp/chrX_pilot_mask.bed
```

and apply the specified filters using `bcftools` and `bedtools`

```         
bcftools view -Ov --types snps -S data/1kgp/unrelated_subjects.txt data/1kgp/chrX.vcf.gz | \
bcftools view -Ov -t ^chrX:10001-2781479 | \
bcftools view -Ov -t ^chrX:155701382-156030895 | \
bcftools annotate -x ^INFO/AF,INFO/AC,INFO/CM,^FORMAT/GT | \
bedtools intersect -a - -b data/1kgp/chrX_pilot_mask.bed -wa -header | \
bcftools view -Ob > data/1kgp/chrX_2504_snps_noPAR.bcf
```

We also want to remove positions that have more than one entry (multiallelic sites):

```         
bcftools query -f "%POS\n" data/1kgp/chrX_2504_snps_noPAR.bcf | \
awk '{count[$1]++}END{for(key in count){if(count[key]>1){print(key)}}}' | \
sort -n | \
sed -e 's/^/chrX\t/' > data/1kgp/repeat_sites.txt

bcftools view -T ^data/1kgp/repeat_sites.txt -Ob data/1kgp/chrX_2504_snps_noPAR.bcf > data/1kgp/chrX_2504_snps_noPAR2.bcf
rm data/1kgp/chrX_2504_snps_noPAR.bcf
mv data/1kgp/chrX_2504_snps_noPAR2.bcf data/1kgp/chrX_2504_snps_noPAR.bcf
```

And finally, we will remove singletons as well:

```         
bcftools view -Ov data/1kgp/chrX_2504_snps_noPAR.bcf | vcftools --singletons --vcf -
awk '{if($1 == "chrX")print($1"\t"$2)}' out.singletons > data/1kgp/singletons.txt
rm out.log out.singletons

bcftools view -Ob -T ^data/1kgp/singletons.txt data/1kgp/chrX_2504_snps_noPAR.bcf > data/1kgp/chrX_2504_snps_noPAR_noSing.bcf
```

### MAF distribution in X

```
bcftools +fill-tags data/1kgp/chrX_2504_snps_noPAR_noSing.bcf -- -t AF | \
bcftools query -f "%CHROM\t%POS\t%AF\n" | \
  awk '{maf = (1-$3 < $3 ? 1-$3 : $3); print($1"\t"$2"\t"maf)}' > data/1kgp/chrX_freq.tsv
```

### GC Content in Windows on X and 15

```
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz > data/ref_GRCh38.fna.gz
gunzip data/ref_GRCh38.fna.gz
fastaindex -f data/ref_GRCh38.fna -i data/ref_GRCh38.fna.in
fastafetch -f data/ref_GRCh38.fna -i data/ref_GRCh38.fna.in -q "chrX" > data/chrX.fasta
fastafetch -f data/ref_GRCh38.fna -i data/ref_GRCh38.fna.in -q "chr15" > data/chr15.fasta
fastafetch -f data/ref_GRCh38.fna -i data/ref_GRCh38.fna.in -q "chr1" > data/chr1.fasta

bedtools complement -i data/1kgp/chrX_pilot_mask.bed -g data/chrX.sizes | bedtools sort > data/mask_pilot_X_fasta.bed

bedtools complement -i data/1kgp/chr15/chr15_pilot_mask.bed -g data/chr15.sizes | bedtools sort > data/mask_pilot_15_fasta.bed

bedtools complement -i data/1kgp/chr1/chr1_pilot_mask.bed -g data/chr1.sizes | bedtools sort > data/mask_pilot_1_fasta.bed

bedtools maskfasta -fi data/chrX.fasta -bed data/mask_pilot_X_fasta.bed -fo data/chrX_mask.fasta
samtools faidx data/chrX_mask.fasta

bedtools maskfasta -fi data/chr15.fasta -bed data/mask_pilot_15_fasta.bed -fo data/chr15_mask.fasta
samtools faidx data/chr15_mask.fasta

bedtools maskfasta -fi data/chr1.fasta -bed data/mask_pilot_1_fasta.bed -fo data/chr1_mask.fasta
samtools faidx data/chr1_mask.fasta


bedtools makewindows -g data/chrX.sizes -w 1000 | sort -k 1,1V -k2,2n > data/chrX.1kb.sorted.bed
bedtools makewindows -g data/chr15.sizes -w 1000 | sort -k 1,1V -k2,2n > data/chr15.1kb.sorted.bed
bedtools makewindows -g data/chr1.sizes -w 1000 | sort -k 1,1V -k2,2n > data/chr1.1kb.sorted.bed

bedtools nuc -fi data/chrX_mask.fasta -bed data/chrX.1kb.sorted.bed > data/chrX_gc1kb_pilot.bed
bedtools nuc -fi data/chr15_mask.fasta -bed data/chr15.1kb.sorted.bed > data/chr15_gc1kb_pilot.bed
bedtools nuc -fi data/chr1_mask.fasta -bed data/chr1.1kb.sorted.bed > data/chr1_gc1kb_pilot.bed
```

## Downloading and processing Chromosome 15

For our analysis of trios, we need to download both the phased and unphased VCFs from 1kGP:

```
chrom=1
curl "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz" > "data/1kgp/chr${chrom}/chr${chrom}_phased.vcf.gz"

curl "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi" > "data/1kgp/chr${chrom}/chr${chrom}_phased.vcf.gz.tbi"

curl "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chrom}.recalibrated_variants.vcf.gz" > "data/1kgp/chr${chrom}/chr${chrom}_unphased.vcf.gz"

curl "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chrom}.recalibrated_variants.vcf.gz.tbi" > "data/1kgp/chr${chrom}/chr${chrom}_unphased.vcf.gz.tbi"
```

For the phased VCF

1.  Remove excess annotation
2.  Remove non-snps
3.  Remove multiallelic sites
4.  Apply filter=PASS

```         
chrom=1
bcftools annotate -x ^INFO/AF,INFO/AC,^FORMAT/GT "data/1kgp/chr${chrom}/chr${chrom}_phased.vcf.gz" | \
bcftools view -v snps -Ob > "data/1kgp/chr${chrom}/chr${chrom}_phased.bcf"

bcftools query -f "%POS\n" "data/1kgp/chr${chrom}/chr${chrom}_phased.bcf" | \
awk '{count[$1]++}END{for(key in count){if(count[key]>1){print(key)}}}' | \
sort -n | \
sed -e "s/^/chr${chrom}\t/" > "data/1kgp/chr${chrom}/repeat_sites_phased.txt"

bcftools view -T ^"data/1kgp/chr${chrom}/repeat_sites_phased.txt" -Ob "data/1kgp/chr${chrom}/chr${chrom}_phased.bcf" > "data/1kgp/chr${chrom}/chr${chrom}_phased_biallelic.bcf"


bcftools view "data/1kgp/chr${chrom}/chr${chrom}_phased_biallelic.bcf" | \
bedtools intersect -a - -b "data/1kgp/chr${chrom}/chr${chrom}_pilot_mask.bed" -wa -header | \
bcftools view -Ob > "data/1kgp/chr${chrom}/chr${chrom}_phased_mask.bcf"

bcftools view "data/1kgp/chr${chrom}/chr${chrom}_phased_mask.bcf" | vcftools --singletons --out phased --vcf -
awk -v var="$chrom" '{if($1 == "chr"var)print($1"\t"$2)}' phased.singletons > "data/1kgp/chr${chrom}/singletons_phased.txt"
rm phased.log phased.singletons

bcftools view -Ob -T ^"data/1kgp/chr${chrom}/singletons_phased.txt" "data/1kgp/chr${chrom}/chr${chrom}_phased_mask.bcf" > "data/1kgp/chr${chrom}/chr${chrom}_phased_mask_noSing.bcf"

bcftools view -Ob -S data/1kgp/unrelated_subjects.txt "data/1kgp/chr${chrom}/chr${chrom}_phased_mask_noSing.bcf" > "data/1kgp/chr${chrom}/chr${chrom}_phased_mask_noSing_2504.bcf"

bcftools index "data/1kgp/chr${chrom}/chr${chrom}_phased_mask_noSing.bcf"
```

For the (large) unphased VCF we:

1.  Remove excess annotation
2.  Limit to snps
3.  Filter to list of sites in phased bcf

We also remove sites from the phased VCF that are biallelic in the unphased VCF

```
chrom=1
bcftools annotate -x ^INFO/AF,INFO/AC,INFO/DP,^FORMAT/GT "data/1kgp/chr${chrom}/chr${chrom}_unphased.vcf.gz" | \
bcftools view -v snps -Ob > "data/1kgp/chr${chrom}/chr${chrom}_unphased.bcf"

bcftools index "data/1kgp/chr${chrom}/chr${chrom}_unphased.bcf"

bcftools isec -n~11 -c all "data/1kgp/chr${chrom}/chr${chrom}_phased_mask_noSing.bcf" "data/1kgp/chr${chrom}/chr${chrom}_unphased.bcf" > "data/1kgp/chr${chrom}/overlap_sites.txt"


bcftools view -T "data/1kgp/chr${chrom}/overlap_sites.txt" "data/1kgp/chr${chrom}/chr${chrom}_unphased.bcf" | \
bcftools view -m2 -M2 -Ob > "data/1kgp/chr${chrom}/chr${chrom}_unphased_overlap.bcf"


bcftools query -f "%CHROM\t%POS\n" "data/1kgp/chr${chrom}/chr${chrom}_unphased_overlap.bcf" > "data/1kgp/chr${chrom}/overlap_sites_biallelic.txt"


bcftools view -T "data/1kgp/chr${chrom}/overlap_sites_biallelic.txt" -Ob "data/1kgp/chr${chrom}/chr${chrom}_phased_mask_noSing.bcf" > "data/1kgp/chr${chrom}/chr${chrom}_phased_overlap.bcf"

bcftools view -Ob -S data/1kgp/unrelated_subjects.txt "data/1kgp/chr${chrom}/chr${chrom}_phased_overlap.bcf" > "data/1kgp/chr${chrom}/chr${chrom}_phased_overlap_2504.bcf"
```

### Extract additional information from vcf

#### Count triple heterozygous sites per trio

```
chrom=1
mkdir -p "data/1kgp/chr${chrom}/triple_het/trio_lists"

Rscript code/trio_family_lists.R

mkdir -p "output/trio_phase_${chrom}/slurm"
sbatch code/batch_triple_hets.sh


for i in `seq 1 602`; do
wc -l "data/1kgp/chr${chrom}/triple_het/sample_${i}.tsv" >> "output/trio_phase_${chrom}/triple_het.tmp"
done


awk '{print(NR"\t"$1)}' "output/trio_phase_${chrom}/triple_het.tmp" > "output/trio_phase_${chrom}/triple_het.tsv"
#here
```

## Phasing Synthetic Diploids

First, we sample males from the same sub-populations using the script `code/sample_X_pairs.R`. This generates the file `data/sample_pairs.csv`.

We then submit the batch script `code/batch_phase_sd.sh` to generate the synthetic diploids and their unique reference panel, phase using Beagle5, Eagle2, and SHAPEIT5, and compare the inferred phase to the original haplotypes. Note: you will need to change the variables to reflect the paths on your own computer.

Additional steps:

1.Annotate each heterozygous position with CpG status and 3-mer motif (`code/batch_annotate_het.sh`)
2.Generate `phased_size.txt` in the directory `output/switch_errors/het_loc`

```
for i in `seq 1 1000`; do
size=$(awk 'NR == 1 {val = $2}END{print $2 - val}' output/switch_errors/het_loc/pair_${i}_het_loc.txt)
echo $i"\t"$size >> output/switch_errors/het_loc/phased_size.txt
done
```

3.Generate count of heterozygous sites in each synthetic diploid:

```
for i in `seq 1 1000`; do
wc -l output/switch_errors/het_loc/pair_${i}_het_loc.txt >> output/switch_errors/het_pos_count.tmp
done

awk '{print(NR"\t"$1)}' output/switch_errors/het_pos_count.tmp > output/switch_errors/het_pos_count.tsv
```

4.Count number of hets that are at CpG sites for each synthetic diploid

```
for i in `seq 1 1000`; do
awk -F',' '{n_cpg += $3}END{print(n_cpg)}' output/switch_errors/het_loc/annotated/pair_${i}.csv  >> output/switch_errors/het_cpg_count.tmp
done

awk '{print(NR"\t"$1)}' output/switch_errors/het_cpg_count.tmp > output/switch_errors/het_cpg_count.tsv
```

5.Get maf at variant sites

```
bcftools view data/1kgp/chrX_2504_snps_noPAR_noSing.bcf | vcftools --vcf - --freq --out chrX_freq
awk 'NR>1 {split($5,a,":");split($6,b,":");c = (a[2]<b[2])?a[2]:b[2]; print($1"\t"$2"\t"c)}' chrX_freq.frq > data/1kgp/chrX_maf.tsv
```

6. Get proportions of errors in maf categories: `code/maf_distributions.R`

## Phasing Trio Children

For each of the 602 children of the trios, we need to identify which samples need to be excluded from the reference panel when phasing (i.e., either removing the parents, the child themselves, or no one in the case neither the parents nor the child were part of the 2,504 phase 3 samples). The script in `code/trio_exclude_lists.R` generates a file for each child with a list of IDs to exclude from the reference panel for the children.


## Phasing chromosome 15

```
chrom=1
out_dir="/net/snowwhite/home/beckandy/research/phasing_clean/output/trio_phase_${chrom}"
mkdir -p $out_dir/beagle
mkdir $out_dir/eagle
mkdir $out_dir/shapeit
mkdir -p $out_dir/het_loc/annotated
mkdir $out_dir/truth
mkdir -p $out_dir/switch_errors/beagle/annotated
mkdir -p $out_dir/switch_errors/eagle/annotated
mkdir -p $out_dir/switch_errors/shapeit/annotated
```

The batch script `code/batch_phase_trios.sh` submits a job to SLURM for each of the 602 trios and generates 4 phased VCFs: 1 with the "true" phase for the child, and 3 phases infered using only the reference panel. The batch script `code/batch_compare_trio.sh` generates a list of switch error locations. A summary of these switches is generated by the script `code/trio_phase_results.R` and analysis is documented in `analysis/trio_phasing_results.Rmd`. We also extract some additional summary information on the samples (e.g., number of heterozygous sites)

```
for i in `seq 1 602`; do
wc -l output/trio_phase15/het_pos/sample_${i}.bed >> output/trio_phase_15/het_pos_count.tmp 
done

awk '{print(NR"\t"$1)}' output/trio_phase_15/het_pos_count.tmp > output/trio_phase_15/het_pos_count.tsv
```

### Count hets at CpG in each trio

```
chrom=1
for i in `seq 1 602`; do
awk -F',' '{n_cpg += $3}END{print(n_cpg)}' "output/trio_phase_${chrom}/het_loc/annotated/pair_${i}.csv " >> "output/trio_phase_${chrom}/het_cpg_count.tmp"
done

awk '{print(NR"\t"$1)}' output/trio_phase_15/het_cpg_count.tmp > output/trio_phase_15/het_cpg_count.tsv
```

### Get maf at variant sites

```
chrom=22
bcftools +fill-tags "data/1kgp/chr${chrom}/chr${chrom}_phased_overlap.bcf" -- -t AF | \
  bcftools query -f "%CHROM\t%POS\t%AF\n" | \
  awk '{maf = (1-$3 < $3 ? 1-$3 : $3); print($1"\t"$2"\t"maf)}' > "data/1kgp/chr${chrom}/chr${chrom}_freq.tsv"
```
