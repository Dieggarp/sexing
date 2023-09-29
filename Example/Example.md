This fold contains genomic information of 10 sea lions that can be used as toy examples for our method.

This dataset allows testing the method starting from step 5, 'SNPs Identification.' To use this VCF file, the input command for the 'populations' tool should be:

```{bash,eval=FALSE}
populations --vcf example.Ofla.vcf --max-obs-het 0.7 -R 0.3 -t 8 -O populations --vcf-all
```

For the entire data set consult Peralta et al. (under review) A rapid approach for sex assignment by RAD-seq using a reference genome
