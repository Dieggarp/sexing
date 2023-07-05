# **"When sex scarce: A rapid approach to sexing individuals by RAD-seq using a reference genome"**
**Diego M. Peralta, Juan I. Túnez, Ulises E. Rodríguez Cruz, Santiago G. Ceballos**

2023-06-14



### Abstract

Assigning sex to individuals without previous information is a common objective of molecular ecology. Here, we developed a framework for sexing animals by using two indexes based on the different properties of the mammalian sexual chromosomes. We mapped RAD-seq loci to a reference genome to obtain missingness and coverage depth from chromosomes Y, X and autosomal, which allowed identifying the sex of fur seals from a previous study with previous sex information. Moreover, we sexed 38 sea lions sampled non-invasively, allowing us to discuss our indexes’ reliability at different coverage depths. We believe this approach could extrapolate to any mammal species or taxa with known XY sex chromosome systems and different qualities of the GBS sequencing.

### Objective

Determining the sex of different individuals using GBS and a reference genome, based on the distinct properties of mammalian sex chromosomes X and Y.


### Requirements
The requirements are listed according to their utilization.

* Linux
* BWA (0.7.17-r1188)
* Samtools (1.3.1)
* Stacks (2.60)
* VCFtools (0.1.17)
* R (4.3.0)

Versions may differ from those used here.


### Usage
Here is a detailed step-by-step procedure, beginning with raw sequences and concluding with sex identification.


**Alingment to a Reference Genome**

1.Create a directory, here we call it sexing. Then create the follow subdirectories.
```{bash,eval=FALSE}
sexing/
  files
  bwa
  alignments.bwa
  gstacks.bwa
  populations
```

2.Prepare a tab-delimited population file, popfile.tsv also in the files/ directory. This file describes the populations to which the samples belong.
```{bash,eval=FALSE}
"individual1"TAB"population1"
"individual2"TAB"population2"
```

In our case the file looks like:
```{bash,eval=FALSE}
CB02  1
Q06 1
Q06B	1
PQ03	1
```


**BWA**

3.Create a Reference Genome Index. Download the reference genome file of the most related organism to your species and put it in files directory.
```{bash,eval=FALSE}
bwa index -p bwa/gac files/genome.fasta > bwa/bwa_index.oe
```

|   Single end
```{bash,eval=FALSE}
bwa mem -M -t 16 bwa/gac location_demultiplex/individual1.fq.gz | samtools view -b > ./alignments.bwa/sample1.bam
```

|   Paired end
```{bash,eval=FALSE}
bwa mem -M -t 16 bwa/gac location_demultiplex/individual1.R1.fq.gz location_demultiplex/individual1.R2.fq.gz | samtools view -b > ./alignments.bwa/sample1.bam
```

|   And, next
```{bash,eval=FALSE}
samtools sort -o alignments.bwa/sample1.bam alignments.bwa/sample1.bam
```


3.1.To avoid making the previous steps for each sample, we appealed to a loop. This must be running from the directory where demultiplex files are stored. *nohup* is an option that allows running programs in the background.

|   Single end
```{bash,eval=FALSE}
nohup sh -c 
'while read K;
 do
  bwa mem -M -t 16 ../sexing/bwa/gac $K.1.fq.gz | samtools view -b >  ../sexing/alignments.bwa/$K.bam;
 done < names.list.txt' 
&
```

|   Paired end
```{bash,eval=FALSE}
nohup sh -c 
'while read K;
 do 
  bwa mem -M -t 16 ../sexing/bwa/gac $K.1.fq.gz $K.2.fq.gz | samtools view -b > ../sexing/alignments.bwa/$K.bam;
 done < names.list.txt' 
&
```

|   The file "names.list.txt" lists the names used in the analysis sequences of all the individuals, as follows.
```{bash,eval=FALSE}
individual1
individual2
individual3
```


3.2.After this, we sorted the files for each individual in a single step.
```{bash,eval=FALSE}
nohup sh -c 
'for K in *.bam;
 do 
  samtools sort -o $K $K;
 done' 
&
```


**gstacks**

4.Align reads to reference genome. gstacks module will examine a RAD dataset one locus at a time, looking at all individuals in the metapopulation for that locus.
```{bash,eval=FALSE}
gstacks -I alignments.bwa/ -M files/popfile.tsv -t 24 -O gstacks.bwa/
```


**SNPs Identification**

5.Creates a catalog of SNPs with populations module. The SNPs filtering steps were aimed to retain as many markers associated with sex chromosomes as possible. Thus, considering that Y linked loci are present only in males, we set the minimum proportion of individuals across populations to process a locus (-R) to 0.3, which assumes a minimum of 30 % of males in our data set.
```{bash,eval=FALSE}
populations -P gstacks.bwa --popmap files/popfile.tsv --max-obs-het 0.7 -R 0.3 -t 8 -O populations --vcf
```


**Sex identification**

6.Enter to the "populations" directory and obtain files of missingness and depth of coverage of each individual for chromosomes X, Y and the autosomal. VCFtools --chr flag uses SNPs contained in a specified chromosome. Here we attempted to automate the process. First create a list of the names of chromosomes under interest. To do this, enter to the chosen species genome directory from https://www.ncbi.nlm.nih.gov/genome/ and copy the chromosome identifiers from the column "RefSeq". In our case: Y (NC_045613.1), X (NC_045612.1), and chromosome 10 (NC_045604.1). After that, paste them into a new file we call chromosome.list.txt.
```{bash,eval=FALSE}
nano chromosome.list.txt
```

Our file looks like:
```{bash,eval=FALSE}
NC_045613.1
NC_045612.1
NC_045604.1
```

After saving, make a directory to deposit the new files (we call it "sexing") and make the next loop to create the files of interest in one step.
```{bash,eval=FALSE}
cd sexing

cat ../chromosome.list.txt | while read line;
 do
  vcftools --chr $line --vcf ../populations.snps.vcf --depth --out ${line}

  vcftools --chr $line --vcf ../populations.snps.vcf --missing-indv --out ${line}
 done

rename 's/$/.tsv/' *.i*
```

Simplify the tables leaving only the first and last columns (sample names and data of interest).
```{bash,eval=FALSE}
for file in *.tsv;
 do
  awk '{print $1, $NF }' $file > ${file/.tsv/_2_.tsv};
 done
```

And named the columns with the name of the chromosome identifiers.

```{bash,eval=FALSE}
for f in *.imiss_2_.tsv;
 do
  n=${f%%.imiss_2_.tsv}
  filename=`echo ${f:r}`; sed -i -e "s/F_MISS/F_MISS_$n/" $f;
 done
```

```{bash,eval=FALSE}
for f in *.idepth_2_.tsv;
 do
  n=${f%%.idepth_2_.tsv}
  filename=`echo ${f:r}`; sed -i -e "s/MEAN_DEPTH/M_DEPTH_$n/" $f;
 done
```

6.1.Run the script "sexing.R" inside the "sexing" directory to calculate Index_X and Index_Y and to reveal the sex of each individual.
```{bash,eval=FALSE}
Rscript sexing.R
```

Two new files will appear when it finishes, "final_sexing.csv", containing results and "sexing_plots.pdf", containing the different plots.

![image](https://github.com/Dieggarp/Sexing/assets/88154471/d2a87809-3933-4c8a-ae38-dcbe92db2f90)

