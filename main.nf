#!/usr/bin/env nextflow

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --bams sample.bam [Options]
    
    Inputs Options:
    --input         Input file

    Resource Options:
    --max_cpus      Maximum number of CPUs (int)
                    (default: $params.max_cpus)  
    --max_memory    Maximum memory (memory unit)
                    (default: $params.max_memory)
    --max_time      Maximum time (time unit)
                    (default: $params.max_time)
    See here for more info: https://github.com/lifebit-ai/hla/blob/master/docs/usage.md
    """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

// Define channels from repository files
projectDir = workflow.projectDir

// Define Channels from input
Channel
    .fromPath(params.input)
    .ifEmpty { exit 1, "Cannot find input file : ${params.input}" }
    .splitCsv(skip:1)
    .map {sample_name, file_path -> [ sample_name, file_path ] }
    .set { ch_input }

process samtools {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    set val(sample_name), file(bam_file) from ch_input_1
    file(reference) from ch_reference_1
    
    output:
    file "*.cram"

    script:
    """
    samtools view -T $reference -C -o ${sample_name}_samtools.cram $bam_file
    """
  }

process cramtools {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    set val(sample_name), file(bam_file) from ch_input_2
    file(reference) from ch_reference_2
    
    output:
    file "*.cram"

    script:
    """
    samtools faidx $reference
    java -jar /cramtools/cramtools-2.1.jar cram --lossless-quality-score --capture-all-tags -I $bam_file -R $reference -O ${sample_name}_cramtools.cram
    """
  }

process gzip {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    set val(sample_name), file(bam_file) from ch_input_3
    
    output:
    file "*.gz"

    script:
    """
    gzip -k $bam_file
    mv *.gz ${sample_name}_gzip.gz
    """
  }

process pigz {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    set val(sample_name), file(input_file) from ch_input
    file(run_sh_script) from ch_run_sh_script
    
    output:
    file "*.gz"

    script:
    """
    pigz filename 
    mv *.gz ${sample_name}_pigz.gz
    """
  }

