# SFARI Mitochondrial Calling Pipeline
This is the calling pipeline that has been used to call mitochondrial variants on SFARI Genomes. It performs the follow tasks:
* Single sample mitochondrial variant calling from Genome bam/cram files.
* Merging of single sample variant calls (i.e. vcf files) into a combined Hail matrix table. 
* Summary annotation and export as a Hail table for import into Elastic search database used by the SFARI browser. 

All code is based the gnomAD mitochondrial pipeline and scripts developed and tested for Google Cloud Platform and has been adapted to work on Institutional clusters.

## Requirements

* Java 1.8
* Python
* R
* bwa
* GATK (tested on 4.1.7.0)
* Picard (tested on 2.18.27)
* HaplocheckCLI (custom jar from [Leklab repository](https://github.com/leklab/haplocheckCLI))
* Cromwell (tested on v56)
* Hail

### Notes
GATK version used in the Google cloud WDL can be extracted from the respective docker image. In snapshot version 25 this can be extracted from `us.gcr.io/broad-gatk/gatk:4.1.7.0` contained in the path `/root/gatk.jar`  


Similarly Picard jar can be extracted from the image `us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386` contained in the path `/usr/gitc/picard.jar`  
  

HaplocheckCLI jar can be extracted from the image `us.gcr.io/broad-dsde-methods/haplochecker:haplochecker-0124` contained in the path `/usr/mtdnaserver/haplocheckCLI.jar`. Alternatively this can be built from source from the [Leklab forked repository](https://github.com/leklab/haplocheckCLI).  


## Single sample variant calling

### Running WDL using Cromwell
The Mitochondrial WDL pipeline is run using cromwell using the following command line.
```
java -Dconfig.file=slurm.conf -jar \
/home/ml2529/shared/tools/jars/cromwell-56.jar run \
MitochondriaPipeline_Snapshot25.wdl \
-i mito_pipeline_inputs.json  \
-o workflow.options

```

The `launch.sh` file is a batch script that can be used to launch in Slurm by the following command  
```
sbatch launch.sh
```

### Inputs
The main input file is specified in `mito_pipeline_inputs.json`. In summary this contains the location of the follow input files:
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

This **needs to be edited** to point to location of these files on the local file system. All input files needed are in the `./wdl_input_files` directory. The GATK, Picard and HaplocheckCLI jar files can be obtained as detailed in Requirements section.

The `workflow.options` input file contains where output files should be copied after successful completion. It also switches on the folliowing convenient features:
* `"use_relative_output_paths": true` - output files are copied to the output directory but are not embedded in the deep directory structure that cromwell/WDL produces.
* `"delete_intermediate_output_files": true` - Removes intermediate output files that are not part of the final output.
* `write_to_cache": true, "read_from_cache": true` - Caches WDL workflow progress so pipeline can resume from the failure point instead of from scratch.  

Please refer to advanced configuration for more details.

### Outputs
The successful run of the WDL based pipeline should produce the following files copied to the output directory:
* Mitochondrial genome (MT) aligned BAM file
* Variant calls (vcf file) from Mutect2 
* VCF file with multi-allelic variants split into bi-allelic variants
* Per base coverage file

It also contains additional useful files/metrics
* Unaligned mitochondrial genome (MT) reads as BAM file extracted from genome BAM
* Input vcf file used for haplocheckCLI
* Mark duplicates metrics
* Coverage metrics
* Theoretical sensitivity metrics
* Mean coverage
* Contamination metrics file (from haplocheckCLI)
* Contamination fraction estimate
* Major haplogroup

### Advanced configuration

#### Scheduling specific configuration
The `slurm.conf` is specific for Slurm scheduling commands. For this to work using a different scheduling backend this will need to be edited accordingly.

#### Resume from failures via call caching
A useful feature that is not configured by default in cromwell based WDL is the ability to resume from a failure point. It is quite difficult to set up due to the lack of documentation. This pipeline uses MySQL to log workflow progress and intermediate files. Please refer to the database section in the `slurm.conf` and edit the below and changing the database `url`, `user` and `password` fields for it to work.

```
database {
  # mysql example
  #driver = "slick.driver.MySQLDriver$" #old way

  profile = "slick.jdbc.MySQLProfile$"


  # see all possible parameters and default values here:
  # http://slick.lightbend.com/doc/3.2.0/api/index.html#slick.jdbc.JdbcBackend$DatabaseFactoryDef@forConfig(String,Config,Driver):Database
  # https://dev.mysql.com/doc/connector-j/8.0/en/connector-j-reference-jdbc-url-format.html

  db {
    driver = "com.mysql.cj.jdbc.Driver"
    url = "jdbc:mysql://chdgenes.org/cromwell?rewriteBatchedStatements=true&useSSL=false"
    user = "cromwell"
    password = "L3kK1ds2018"
    connectionTimeout = 5000
  }

  # For batch inserts the number of inserts to send to the DB at a time
  insert-batch-size = 2000

}
```


## Merging single sample variant calls
The hail scripts used for gnomAD mitochondrial calling pipeline was modified. This is done in two stages:

### 1. Creating a Hail matrix table for coverage
The single sample vcf produced by the above WDL pipeline only shows variant sites. This becomes a challenge when merging vcfs where one sample has a variant and the other does not. This is typically solved by GATK by using gVCF files with the concept of reference blocks. This pipeline uses per base coverage (generated above) to fill in missing information, that is the coverage of a homoplasmic reference site on samples that do not contain a variant at that site.  

This step prepares the reference sample by using the follow script.
`./hail/annotate_coverage.py`

### 2. Merging and filtering single sample vcf files
The following script takes the per base coverage of each sample as a Hail matrix table and other filtering parameters to perform the merging.  
`./hail/combine_vcfs.py`

The final output files are a Hail matrix table and compressed vcf file.  

## Additional annotation for import into Elastic search database
Currently in development. The follow gnomAD hail code will be adapted for the SFARI project.
* `./hail/add_annotations.py`




