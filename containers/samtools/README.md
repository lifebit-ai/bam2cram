# samtools

This bcftool contaier - `lifebitai/samtools:s3` is especailly designed to working with s3 files.

The warpper written based on example - https://github.com/delagoya/s3samtools

Try it - 

```
docker run -v $PWD:$PWD -w $PWD --rm lifebitai/samtools:s3 samtools view s3://igv.org.test/csi_test/chr10p.bam | head -n 1
```
