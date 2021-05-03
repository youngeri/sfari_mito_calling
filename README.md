# SFARI Mitochondrial Calling Pipeline
This is calling pipeline has been used to call mitochondrial variants on SFARI Genomes. It performs the follow tasks:
* Single sample mitochondrial variant calling from Genome bam/cram files.
* Merging of single sample variant calls (i.e. vcf files) into a combined Hail matrix table. 
* Summary annotation and export as a Hail table for import into Elastic search database used by the SFARI browser. 

All code is based the gnomAD mitochondrial pipeline and scripts and has been adapted to work on Institutional clusters.

## Requirements

* Java 1.8
* Python
* R
* bwa
* GATK
* Picard
* HaplocheckCLI
* Hail








