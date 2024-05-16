# Phasing Algorithm Evaluation with Synthetic Diploids

This repository is an accompaniment to the paper "SUGAR TEA COFFEE FOR FREE"

## Data included in this repository

* `data/1kgp/unrelated_subjects.txt`: A list of IDs for the 2,502 unrelated samples from 1kGP phase 3

## Downloading and Processing Reference Data

To generate our synthetic diploids, we use start with the [1000 Genomes Project 30x on GRCh38 Phased VCF](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/) for chromosome X:

```         
curl ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz > data/1kgp/chrX.vcf.gz
curl ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz.tbi > data/1kgp/chrX.vcf.gz.tbi
```

We apply the [1kGP Pilot accessibility map](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/PilotMask/) to remove hard to sequence regions of the chromosome (note that this bed file is a list of sites that we want to keep in the VCF):

```
curl ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/PilotMask/20160622.allChr.pilot_mask.bed | grep "chrX" > data/1kgp/chrX_pilot_mask.bed
```
