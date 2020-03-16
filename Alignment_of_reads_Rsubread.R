library(Rsubread)

##Earlier we put all the sequencing read data (.fastq.gz files) in the data directory. Now we need to find them in order to tell the Rsubread aligner which files to look at. We can search for all .fastq.gz files in the data directory using the list.files command. The pattern argument takes a regular expression. In this case we are using the $ to mean the end of the string, so we will only get files that end in “.fastq.gz”
setwd('C:/Users/Tony Isebe/Desktop/RNASEQANALYSIS/Mouse_mammary')

##load in the fastq files for alignment

fastqfiles <- list.files('C:/Users/Tony Isebe/Desktop/RNASEQANALYSIS/Mouse_mammary', pattern = '.fastq.gz$', full.names = TRUE)

fastqfiles

##Perform the alignment

#First build an index

##Read sequences are stored in compressed (gzipped) FASTQ files. Before the differential expression analysis can proceed, these reads must be aligned to the mouse genome and counted into annotated genes. This can be achieved with functions in the Rsubread package.
##The first step in performing the alignment is to build an index. In order to build an index you need to have the fasta file (.fa), which can be downloaded from the UCSC genome browser. Here we are building the index just for chromosome 1. This may take several minutes to run. Building the full index using the whole mouse genome usually takes about 30 minutes to an hr on a server. We won’t be building the index in the workshop due to time constraints, we have provided the index files for you. The command below assumes the chr1 genome information for mm10 is stored in the “chr1.fa” file.


buildindex(basename = 'chr12_mm10', reference = 'chr12.fa')

##The above command will generate the index files in the working directory. In this example, the prefix for the index files is chr1_mm10. You can see the additional files generated using the dir command, which will list every file in your current working directory.

#ALIGN READS TO CHROMOSOME 1 OF THE REFERENCE GENOME

#Now that we have generated our index, we can align our reads using the align command. There are often numerous mapping parameters that we can specify, but usually the default mapping parameters for the align function are fine. If we had paired end data, we would specify the second read file/s using the readfile2 argument. Our mouse data comprises 100bp single end reads.
#We can specify the output files, or we can let Rsubread choose the output file names for us. The default output file name is the filename with “.subread.BAM” added at the end.
#Now we can align our 12 fastq.gz files using the align command.
#This will align each of the 12 samples one after the other. As we’re only using a subset of 1000 reads per sample, aligning should just take a minute or so for each sample. To run the full samples from this dataset would take several hours per sample. The BAM files are saved in the working directory.
#To see how many parameters you can change try the args function:

#In this example we have kept many of the default settings, which have been optimised to work well under a variety of situations. The default setting for align is that it only keeps reads that uniquely map to the reference genome. For testing differential expression of genes, this is what we want, as the reads are unambigously assigned to one place in the genome, allowing for easier interpretation of the results. Understanding all the different parameters you can change involves doing a lot of reading about the aligner that you are using, and can take a lot of time to understand! Today we won’t be going into the details of the parameters you can change, but you can get more information from looking at the help:
align(index = 'C:/Users/Tony Isebe/Desktop/RNASEQANALYSIS/Mouse_mammary/chr12_mm10', readfile1 = fastqfiles)

##We can get a summary of the proportion of reads that mapped to the reference genome using the propmapped function.

bamfiles2 <- list.files(path = 'C:/Users/Tony Isebe/Desktop/RNASEQANALYSIS/Mouse_mammary', pattern = '.BAM$', full.names = TRUE)
bamfiles

props <- propmapped(files = bamfiles)

props


##QUALITY CONTROL

#Extracting quality scores

qs <- qualityScores(filename = 'SRR1552450.fastq.gz', nreads = 100)

#checking dimensions of the quality scores

dim(qs)

##checking the first few elements of qs with head

head(qs)

##A quality score of 30 corresponds to a 1 in 1000 chance of an incorrect base call. (A quality score of 10 is a 1 in 10 chance of an incorrect base call.) To look at the overall distribution of quality scores across the 100 reads, we can look at a boxplot

boxplot(qs1)

qs1 <- qualityScores(filename = 'SRR1552451.fastq.gz', nreads = 50)

##COUNT OF READS MAPPING TO THE GENOME

##Now that we have figured out where each read comes from in the genome, we need to summarise the information across genes or exons. The alignment produces a set of BAM files, where each file contains the read alignments for each library. In the BAM file, there is a chromosomal location for every read that mapped uniquely. The mapped reads can be counted across mouse genes by using the featureCounts function. featureCounts contains built-in annotation for mouse (mm9, mm10) and human (hg19) genome assemblies (NCBI refseq annotation).
##The code below uses the exon intervals defined in the NCBI refseq annotation of the mm10 genome. Reads that map to exons of genes are added together to obtain the count for each gene, with some care taken with reads that span exon-exon boundaries. featureCounts takes all the BAM files as input, and outputs an object which includes the count matrix. Each sample is a separate column, each row is a gene.

fc <- featureCounts(bamfiles2, annot.inbuilt = 'mm10')


##Checking the slots which are stored in fc

names(fc)

##The statistics of the read mapping can be seen with fc$stats. This reports the numbers of unassigned reads and the reasons why they are not assigned (eg. ambiguity, multi-mapping, secondary alignment, mapping quality, fragment length, chimera, read duplicate, non-junction and so on), in addition to the number of successfully assigned reads for each library. See subread documentation (‘Program output’ section). (We know the real reason why the majority of the reads aren’t mapping - they’re not from chr 1!)

##taking a look at the featurecount stats

fc$stat


fc$counts

##The counts for the samples are stored in fc$counts. Take a look at that.
## Take a look at the dimensions to see the number of genes

dim(fc$counts)
## Take a look at the first 6 lines

head(fc$counts)
##The row names of the fc$counts matrix represent the Entrez gene identifiers for each gene and the column names are the output filenames from calling the align function. The annotation slot shows the annotation information that featureCounts used to summarise reads over genes.

head(fc$annotation)

##FURTHER NOTES

#If you are sequencing your own data, the sequencing facility will almost always provide fastq files.
#For publicly available sequence data from GEO/SRA, the files are usually in the Sequence Read Archive (SRA) format. Prior to read alignment, these files need to be converted into the FASTQ format using the fastq-dump utility from the SRA Toolkit. See http: //www.ncbi.nlm.nih.gov/books/NBK158900 for how to download and use the SRA Toolkit.
#By default, alignment is performed with unique=TRUE. If a read can be aligned to two or more locations, Rsubread will attempt to select the best location using a number of criteria. Only reads that have a unique best location are reported as being aligned. Keeping this default is recommended, as it avoids spurious signal from non-uniquely mapped reads derived from, e.g., repeat regions.
#The Phred offset determines the encoding for the base-calling quality string in the FASTQ file. For the Illumina 1.8 format onwards, this encoding is set at +33. However, older formats may use a +64 encoding. Users should ensure that the correct encoding is specified during alignment. If unsure, one can examine the first several quality strings in the FASTQ file. A good rule of thumb is to check whether lower-case letters are present (+64 encoding) or absent (+33).
#featureCounts requires gene annotation specifying the genomic start and end position of each exon of each gene. Rsubread contains built-in gene annotation for mouse and human. For other species, users will need to read in a data frame in GTF format to define the genes and exons. Users can also specify a custom annotation file in SAF format. See the Rsubread users guide for more information, or try ?featureCounts, which has an example of what an SAF file should like like.

