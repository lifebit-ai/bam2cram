#!/usr/bin/env nextflow

def helpMessage() {
    log.info """
    Usage:
    nextflow run main.nf --input input.csv --reference reference.fasta [Options]
    
    Inputs Options:
    --input         Input csv file with bam paths
    --reference     Reference fasta file

    Resource Options:
    --cpus          Number of CPUs (int)
                    (default: $params.cpus)  
    --max_cpus      Maximum number of CPUs (int)
                    (default: $params.max_cpus)
    --memory        Memory (memory unit)
                    (default: $params.memory)   
    --max_memory    Maximum memory (memory unit)
                    (default: $params.max_memory)
    --time          Time limit (time unit)
                    (default: $params.time)
    --max_time      Maximum time (time unit)
                    (default: $params.max_time)
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
    .splitCsv(header:true)
    .map{ row -> file(row.bam) }
    .set { ch_input }

ch_input.into{ch_input_0;
              ch_input_1;
              ch_input_2;
              ch_input_3;
              ch_input_4;
              ch_input_5;
              ch_input_6;
              ch_input_7;
              ch_input_8;
              ch_input_9;
              ch_input_10;
              ch_input_11}

Channel
    .fromPath(params.reference)
    .ifEmpty { exit 1, "Cannot find input file : ${params.reference}" }
    .set { ch_reference }

ch_reference.into{ch_reference_0;
                  ch_reference_1;
                  ch_reference_2;
                  ch_reference_3;
                  ch_reference_4;
                  ch_reference_5;
                  ch_reference_6;
                  ch_reference_7;
                  ch_reference_8;
                  ch_reference_9;
                  ch_reference_10;
                  ch_reference_11}

process samtools_default_30 {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}/${task.process}/", mode: 'copy'

    input:
    file(bam_file) from ch_input_0
    file(reference) from ch_reference_0
    
    output:
    file "*.cram"

    script:
    """
    samtools view -T $reference -o ${bam_file}.cram -O cram,version=3.0 $bam_file
    """
  }

process samtools_default_31 {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}/${task.process}/", mode: 'copy'

    input:
    file(bam_file) from ch_input_1
    file(reference) from ch_reference_1
    
    output:
    file "*.cram"

    script:
    """
    samtools view --threads $task.cpus -T $reference -o ${bam_file}.cram -O cram,version=3.1 $bam_file
    """
  }

process samtools_normal_30 {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}/${task.process}/", mode: 'copy'

    input:
    file(bam_file) from ch_input_2
    file(reference) from ch_reference_2
    
    output:
    file "*.cram"

    script:
    """
    samtools view --threads $task.cpus -T $reference -o ${bam_file}.cram -O cram,version=3.0 --output-fmt-option seqs_per_slice=10000 $bam_file
    """
  }

process samtools_normal_31 {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}/${task.process}/", mode: 'copy'

    input:
    file(bam_file) from ch_input_3
    file(reference) from ch_reference_3
    
    output:
    file "*.cram"

    script:
    """
    samtools view --threads $task.cpus -T $reference -o ${bam_file}.cram -O cram,version=3.1 --output-fmt-option seqs_per_slice=10000 $bam_file
    """
  }

process samtools_fast_30 {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}/${task.process}/", mode: 'copy'

    input:
    file(bam_file) from ch_input_4
    file(reference) from ch_reference_4
    
    output:
    file "*.cram"

    script:
    """
    samtools view --threads $task.cpus -T $reference -o ${bam_file}.cram -O cram,version=3.0,level=1 --output-fmt-option seqs_per_slice=1000 $bam_file
    """
  }

process samtools_fast_31 {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}/${task.process}/", mode: 'copy'

    input:
    file(bam_file) from ch_input_5
    file(reference) from ch_reference_5
    
    output:
    file "*.cram"

    script:
    """
    samtools view --threads $task.cpus -T $reference -o ${bam_file}.cram -O cram,version=3.1,level=1 --output-fmt-option seqs_per_slice=1000 $bam_file
    """
  }

process samtools_small_30 {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}/${task.process}/", mode: 'copy'

    input:
    file(bam_file) from ch_input_6
    file(reference) from ch_reference_6
    
    output:
    file "*.cram"

    script:
    """
    samtools view --threads $task.cpus -T $reference -o ${bam_file}.cram -O cram,version=3.0,level=6,use_bzip2=1 --output-fmt-option seqs_per_slice=25000 $bam_file
    """
  }

process samtools_small_31 {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}/${task.process}/", mode: 'copy'

    input:
    file(bam_file) from ch_input_7
    file(reference) from ch_reference_7
    
    output:
    file "*.cram"

    script:
    """
    samtools view --threads $task.cpus -T $reference -o ${bam_file}.cram -O cram,version=3.1,level=6,use_bzip2=1,use_fqz=1 --output-fmt-option seqs_per_slice=25000 $bam_file
    """
  }

process samtools_archive_30 {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}/${task.process}/", mode: 'copy'

    input:
    file(bam_file) from ch_input_8
    file(reference) from ch_reference_8
    
    output:
    file "*.cram"

    script:
    """
    samtools view --threads $task.cpus -T $reference -o ${bam_file}.cram -O cram,version=3.0,level=7,use_bzip2=1 --output-fmt-option seqs_per_slice=100000 $bam_file
    """
  }

process samtools_archive_31 {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}/${task.process}/", mode: 'copy'

    input:
    file(bam_file) from ch_input_9
    file(reference) from ch_reference_9
    
    output:
    file "*.cram"

    script:
    """
    samtools view --threads $task.cpus -T $reference -o ${bam_file}.cram -O cram,version=3.1,level=7,use_bzip2=1,use_fqz=1,use_arith=1 --output-fmt-option seqs_per_slice=100000 $bam_file
    """
  }

process samtools_archive_lzma_30 {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}/${task.process}/", mode: 'copy'

    input:
    file(bam_file) from ch_input_10
    file(reference) from ch_reference_10
    
    output:
    file "*.cram"

    script:
    """
    samtools view --threads $task.cpus -T $reference -o ${bam_file}.cram -O cram,version=3.0,level=7,use_bzip2=1,use_lzma=1 --output-fmt-option seqs_per_slice=100000 $bam_file
    """
  }

process samtools_archive_lzma_31 {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}/${task.process}/", mode: 'copy'

    input:
    file(bam_file) from ch_input_11
    file(reference) from ch_reference_11
    
    output:
    file "*.cram"

    script:
    """
    samtools view --threads $task.cpus -T $reference -o ${bam_file}.cram -O cram,version=3.1,level=7,use_bzip2=1,use_fqz=1,use_arith=1,use_lzma=1 --output-fmt-option seqs_per_slice=100000 $bam_file
    """
  }