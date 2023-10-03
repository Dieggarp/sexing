# **"A rapid approach for sex assignment by RAD-seq using a reference genome"**
**Diego M. Peralta, Juan I. Túnez, Ulises E. Rodríguez Cruz, Santiago G. Ceballos**

2023-09-28



### Abstract

Sex identification is a common objective in molecular ecology. While many vertebrates display sexual dimorphism, determining the sex can be challenging in certain situations, such as species lacking clear sex-related phenotypic characteristics or in studies using non-invasive methods. In these cases, DNA analyses serve as valuable tools not only for sex determination but also for validating sex assignment based on phenotypic traits. In this study, we developed a bioinformatic framework for sex assignment using genomic data obtained through GBS, and having an available closely related genome assembled at the chromosome level. Our method consists of two ad hoc indexes that rely on the different properties of the mammalian heteromorphic sex chromosomes. For this purpose, we mapped RAD-seq loci to a reference genome and then obtained missingness and coverage depth values for the autosomes and X and Y chromosomes of each individual. Our methodology successfully determined the sex of 165 fur seals that had been phenotypically sexed in a previous study and 40 sea lions sampled in a non-invasive way. Additionally, we evaluated the accuracy of each index in sequences with varying average coverage depths, with Index Y proving greater reliability and robustness in assigning sex to individuals with low-depth coverage. We believe that the approach presented here can be extended to any animal taxa with known heteromorphic XY/ZW sex chromosome systems and that it can tolerate various qualities of GBS sequencing data.

### Objective

Determining the sex of different individuals using GBS and a reference genome, based on the distinct properties of mammalian sex chromosomes X and Y.


### Requirements
The requirements are listed according to their utilization.

* Linux
* BWA (0.7.17-r1188)
* Samtools (1.3.1)
* Stacks (2.64)
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

5. Create a catalog of SNPs with populations module. The SNPs filtering steps were aimed to retain as many markers associated with sex chromosomes as possible. Thus, considering that Y linked loci are present only in males, we set the minimum proportion of individuals across populations to process a locus (-R) to 0.3, which assumes a minimum of 30 % of males in our data set. The flag --vcf-all retrieves all the positions within the RAD-loci, containing fixed and variable sites.
```{bash,eval=FALSE}
populations -P gstacks.bwa --popmap files/popfile.tsv --max-obs-het 0.7 -R 0.3 -t 8 -O populations --vcf-all
```


**Sex identification**

6. Enter to the "populations" directory and obtain files of missingness and depth of coverage of each individual for chromosomes X, Y and the autosomal. VCFtools --chr flag uses loci contained in a specified chromosome, while --not-chr avoids those we chose. Here, we attempted to automate the process. First, create a list of the names of sex chromosomes. To do this, enter to the chosen species genome directory from https://www.ncbi.nlm.nih.gov/genome/ and copy the chromosome identifiers from the column "RefSeq". In our case: Y (NC_045613.1) and X (NC_045612.1). After that, paste them into a new file we call chromosome.list.txt.
```{bash,eval=FALSE}
nano chromosome.list.txt
```

Our file looks like:
```{bash,eval=FALSE}
NC_045613.1
NC_045612.1
```

After saving, make a directory to deposit the new files (we call it "sexing") and make the next loop to create the files of interest in one step.
```{bash,eval=FALSE}
cd sexing

cat ../chromosome.list.txt | while read line;
 do
  vcftools --chr $line --vcf ../populations.all.vcf --depth --out ${line}

  vcftools --chr $line --vcf ../populations.all.vcf --missing-indv --out ${line}
 done

vcftools --not-chr NC_045613.1 --not-chr NC_045612.1 --vcf ../populations.all.vcf --depth --out autosomals

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
Rscript Sexing.R
```

Two new files will appear when it finishes, "final_sexing.csv", containing results and "sexing_plots.pdf", containing the different plots. Next, there are examples of the output files.

First, an example of the table "final_sexing.csv".

![image](https://github.com/Dieggarp/sexing/blob/ce2e44f1648fb7c36b31a94b5044528b97eb42b3/final%20sexig.png)

And second, an example of the figures.

![image](https://github.com/Dieggarp/sexing/blob/0e7e3dc1c820daf3d57f5cedbf192bfd5947484a/PlotArc.jpg) ![image](https://github.com/Dieggarp/sexing/blob/0e7e3dc1c820daf3d57f5cedbf192bfd5947484a/Sup%20Plot%20Arc.jpg)

Figure A corresponds to a dispersion plot of Index Y and X. Figure B corresponds to one of the control figures showing Index X and Y vs Coverage depth. In both cases, red dots include male individuals, and black dots include females.







