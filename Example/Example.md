This fold contains genomic information of 10 sea lions that can be used as toy examples for our method.

Este set de datos permite probar el metodo a partir del paso 5 "SNPs Identification". Para utilizar este vcf el comando de entrada para populations debe ser:

```{bash,eval=FALSE}
populations --vcf example.Ofla.vcf --max-obs-het 0.7 -R 0.3 -t 8 -O populations --vcf-all
```

For the entire data set consult Peralta et al. (under review) A rapid approach for sex assignment by RAD-seq using a reference genome
