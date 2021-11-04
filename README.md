# BAM file compression benchmark pipeline

Benchmarks of CRAM 2.1 and 3.0 are using the faster CRAM codecs; primarily deflate and rANS.

Also included is the performance of the proposed CRAM v3.1 standard. This is not yet a ratified GA4GH standard, but these figures give indicative results.

## Rationale

The following Samtools (version 1.14) profiles are tested:

| Profile      | CRAM versions | options                                                             |
|--------------|---------------|---------------------------------------------------------------------|
| default      | 3.0, 3.1      |                                                                     |
| fast         | 3.0, 3.1      | seqs_per_slice=1000, level=1                                        |
| normal       | 3.0, 3.1      | seqs_per_slice=10000                                                |
| small        | 3.0           | seqs_per_slice=25000, level=6,use_bzip2                             |
| small        | 3.1           | seqs_per_slice=25000, level=6,use_bzip2,use_fqz                     |
| archive      | 3.0           | seqs_per_slice=100000,level=7,use_bzip2                             |
| archive      | 3.1           | seqs_per_slice=100000,level=7,use_bzip2,use_fqz,use_arith           |
| archive lzma | 3.0           | seqs_per_slice=100000,level=7,use_bzip2, use_lzma                   |
| archive lzma | 3.1           | seqs_per_slice=100000,level=7,use_bzip2,use_fqz,use_arith, use_lzma |

##
Requirements
This workflow requires at least 2 CPUs and 4GB of memory.
## Usage

    Usage:
    nextflow run main.nf --input input.csv --reference reference.fasta [Options]

    Inputs Options:
    --input         Input csv file with sample_id, bam and bai paths
    --reference     Reference fasta file

    Resource Options:
    --cpus          Number of CPUs (int)
                    (default: 2)  
    --max_cpus      Maximum number of CPUs (int)
                    (default: 2)
    --memory        Memory (memory unit)
                    (default: 4 GB)   
    --max_memory    Maximum memory (memory unit)
                    (default: 4 GB)
    --time          Time limit (time unit)
                    (default: 8h)
    --max_time      Maximum time (time unit)
                    (default: 8h)

