version 1.0

import "AlignAndCall.wdl" as AlignAndCall

workflow MitochondriaPipeline {

  meta {
    description: "Takes in an hg38 bam or cram and outputs VCF of SNP/Indel calls on the mitochondria."
    allowNestedInputs: true
  }

  input {
    File wgs_aligned_input_bam_or_cram
    File wgs_aligned_input_bam_or_cram_index
    String contig_name = "chrM"
    Float autosomal_coverage = 30

    String gatk
    String picard

    # Full reference is only requred if starting with a CRAM (BAM doesn't need these files)
    File? ref_fasta
    File? ref_fasta_index
    File? ref_dict

    File mt_dict
    File mt_fasta
    File mt_fasta_index
    File mt_amb
    File mt_ann
    File mt_bwt
    File mt_pac
    File mt_sa

    #Shifted reference is used for calling the control region (edge of mitochondria reference).
    #This solves the problem that BWA doesn't support alignment to circular contigs.
    File mt_shifted_dict
    File mt_shifted_fasta
    File mt_shifted_fasta_index
    File mt_shifted_amb
    File mt_shifted_ann
    File mt_shifted_bwt
    File mt_shifted_pac
    File mt_shifted_sa


  }

  parameter_meta {
    wgs_aligned_input_bam_or_cram: "Full WGS hg38 bam or cram"
    autosomal_coverage: "Median coverage of full input bam"
    out_vcf: "Final VCF of mitochondrial SNPs and INDELs"
    vaf_filter_threshold: "Hard threshold for filtering low VAF sites"
    f_score_beta: "F-Score beta balances the filtering strategy between recall and precision. The relative weight of recall to precision."
    contig_name: "Name of mitochondria contig in reference that wgs_aligned_input_bam_or_cram is aligned to"
  }

  call SubsetBamToChrM {
    input:
      gatk = gatk,
      input_bam = wgs_aligned_input_bam_or_cram,
      input_bai = wgs_aligned_input_bam_or_cram_index,
      contig_name = contig_name,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict
  }

  call RevertSam {
    input:
      picard = picard,
      input_bam = SubsetBamToChrM.output_bam,
  }

  String base_name = basename(SubsetBamToChrM.output_bam, ".bam")

  call AlignAndCall.AlignAndCall as AlignAndCall {
    input:
      unmapped_bam = RevertSam.unmapped_bam,
      autosomal_coverage = autosomal_coverage,
      base_name = base_name,
      picard = picard,
      mt_dict = mt_dict,
      mt_fasta = mt_fasta,
      mt_fasta_index = mt_fasta_index,
      mt_amb = mt_amb,
      mt_ann = mt_ann,
      mt_bwt = mt_bwt,
      mt_pac = mt_pac,
      mt_sa = mt_sa,
      mt_shifted_dict = mt_shifted_dict,
      mt_shifted_fasta = mt_shifted_fasta,
      mt_shifted_fasta_index = mt_shifted_fasta_index,
      mt_shifted_amb = mt_shifted_amb,
      mt_shifted_ann = mt_shifted_ann,
      mt_shifted_bwt = mt_shifted_bwt,
      mt_shifted_pac = mt_shifted_pac,
      mt_shifted_sa = mt_shifted_sa
  }


  output {
    File subset_bam = SubsetBamToChrM.output_bam
    File subset_bai = SubsetBamToChrM.output_bai
    File mt_aligned_bam = AlignAndCall.mt_aligned_bam
    File mt_aligned_bai = AlignAndCall.mt_aligned_bai    
  }
}

task SubsetBamToChrM {
  input {
    String gatk
    File input_bam
    File input_bai
    String contig_name
    String basename = basename(basename(input_bam, ".cram"), ".bam")
    File? ref_fasta
    File? ref_fasta_index
    File? ref_dict

  }

  meta {
    description: "Subsets a whole genome bam to just Mitochondria reads"
  }
  parameter_meta {
    ref_fasta: "Reference is only required for cram input. If it is provided ref_fasta_index and ref_dict are also required."
    input_bam: {
      localization_optional: true
    }
    input_bai: {
      localization_optional: true
    }
  }
  command <<<
    java -Xmx4G -jar ~{gatk} PrintReads \
      ~{"-R " + ref_fasta} \
      -L ~{contig_name} \
      --read-filter MateOnSameContigOrNoMappedMateReadFilter \
      --read-filter MateUnmappedAndUnmappedReadFilter \
      -I ~{input_bam} \
      --read-index ~{input_bai} \
      -O ~{basename}.bam
  >>>

  output {
    File output_bam = "~{basename}.bam"
    File output_bai = "~{basename}.bai"
  }
}

task RevertSam {
  input {
    String picard
    File input_bam
    String basename = basename(input_bam, ".bam")
  }

  meta {
    description: "Removes alignment information while retaining recalibrated base qualities and original alignment tags"
  }
  parameter_meta {
    input_bam: "aligned bam"
  }
  command {
    java -Xmx4G -jar ~{picard} \
    RevertSam \
    INPUT=~{input_bam} \
    OUTPUT_BY_READGROUP=false \
    OUTPUT=~{basename}.bam \
    VALIDATION_STRINGENCY=LENIENT \
    ATTRIBUTE_TO_CLEAR=FT \
    ATTRIBUTE_TO_CLEAR=CO \
    SORT_ORDER=queryname \
    RESTORE_ORIGINAL_QUALITIES=false
  }

  output {
    File unmapped_bam = "~{basename}.bam"
  }
}



