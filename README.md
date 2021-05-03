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
* GATK (tested on 4.1.7.0)
* Picard (tested on X)
* HaplocheckCLI (custom jar from repository)
* Cromwell (tested on X)
* Hail

### Notes
GATK version used in the Google cloud WDL can be extracted from the respective docker image. In snapshot version 25 this can be extracted from `us.gcr.io/broad-gatk/gatk:4.1.7.0` contained in the path `/root/gatk.jar`  
Similarly Picard jar can be extracted from the image `us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386` contained in the path `/usr/gitc/picard.jar`  
HaplocheckCLI jar can be extracted from the image `us.gcr.io/broad-dsde-methods/haplochecker:haplochecker-0124` contained in the path `/usr/mtdnaserver/haplocheckCLI.jar`. Alternatively this can be built from source from the Leklab forked repository.  


## Single sample variant calling

### Inputs
The main input file is specified in mito_pipeline_inputs.json. In summary this contains the location of the follow input files:
* Input bam/cram file (from genome sequencing)
* GATK jar file (String)
* Picard jar file (String)
* haplocheckCLI jar file (String)
* hg38 mitochondrial reference sequence and associated index files
* hg38 mitochondrial reference sequence (shifted by 8000 bp) and associated index files
* Shift back chain file
* blacklisted sites (BED file)
* Non-control region interval list
* Control region shifted interval list

This needs to be edited to point to location of these files on the local file system.  

The workflow.options input file contains where output files should be copied after successful completion. It also switches on the folliowing convenient features:
* output files are copied to the output directory but are not embedded in the deep directory structure that cromwell/WDL produces.
* Removes intermediate output files that are not part of the final output.
* Caches WDL workflow progress so pipeline can resume from the failure point instead of from scratch. Please refer to advanced configuration for more details.

### Running WDL using Cromwell



## Merging single sample variant calls
Currently in development


## Additional annotation for import into Elastic search database
Currently in development




