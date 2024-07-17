# CHIP-candidate_mutations
## Build a masked version GRCh38 reference genome
The U2AF1 gene is essential for mRNA splicing and is often mutated in hematopoietic cancers. However, updates in the GRCh38 human reference genome hinder the detection of these mutations by variant calling pipelines as described in [Failure to Detect Mutations in U2AF1 due to Changes in the GRCh38 Reference Sequence](https://doi.org/10.1016/j.jmoldx.2021.10.013).

Based on the instructions detailed in that paper, we built a masked genome for the GRCh38 version of the genome:

1.- [Get GATK Resource Bundle for GRCh38](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)

2.- [Get region to mask from NCBI](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_GRC_exclusions.bed)

3.- Use __bedtools maskfasta__ to mask the fasta file from the GATK Bundle using the .bed file with the regions to be masked

4.- Use __bwa index__ to make the BWA index files for the masked reference fasta file

