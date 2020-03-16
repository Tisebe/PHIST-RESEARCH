setwd('C:/Users/Tony Isebe/Desktop/collo_rna_seq/RNAseq_QC_and_DESeq')

library(pasillaBamSubset)
library(pasilla)

if(!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
##else
library(BiocManager)
BiocManager::install('arrayQualityMetrics')

##locating the file


data.file = system.file('extdata/pasilla_gene_counts.tsv', package = 'pasilla', mustWork = TRUE)
data.file

##Reading the file with R

pasillaCountTable <- read.table(data.file, header = TRUE, row.names = 1)
head(pasillaCountTable)

##Here, header=TRUE indicates that the first line contains column names and row.names=1 means that the first column should be used as row names. This leaves us with a data.frame containing integer count values.


##Find a description of the samples

pasillaDesign <- data.frame(
  row.names = colnames(pasillaCountTable),
  condition=c('untreated','untreated','untreated',
              'untreated','treated','treated','treated'),
  libType=c('single-end','single-end','paired-end',
            'paired-end','single-end','paired-end','paired-end')
)

##Here, simply use a chunk of R code to set up this defne. More often, you will read this into R from a spreadsheet table, using the R functions read.table or read.csv; or perhaps even from a relational database table using R' database access facilities.
#To analyse these samples, we will have to account for the fact that we have both single-end and paired-end method.
#To keep things simple at the start, we defer the discussion of this to Section 4 and first demonstrate a simple analysisby using only the paired-end samples.

pasillaDesign

pairedSamples <- pasillaDesign$libType=='paired-end'
countTable <- pasillaCountTable[, pairedSamples]
condition <- pasillaDesign$condition[pairedSamples]


head(countTable)


condition

##load library

library(DESeq)

##instantiate the Count Data Set

cds <- newCountDataSet(countTable, condition)


##Normalization

#determine effective library size

cds <- estimateSizeFactors(cds)
sizeFactors(cds)


##perform the normalization

head(counts(cds, normalized=TRUE))


##estimate the variance

cds <- estimateDispersions(cds)
plotDispEsts(cds)

##the dispersion values which will be used by the subsequent testing are stored in the feature data slot of cds

head(fData(cds)) ##Fixme: Why is the value for FBgn00000014 Inf and not 0 (given that all counts are 0)? You can verify that these values are indeed the maxima from the two value vectors in fitInfo(cds), which we saw on page 5. Advanced users who want to fiddle with the dispersion estimation can change the values in fData(cds) prior to calling the testing function.

##Calling Differential Expression

##Standard comparison between two experimental conditions

##to see whether there is differential expression between conditions\untreated"and\treated",we simply call the function nbinomTest. It performs the tests as described in [1] and returns a data frame with the pvalues and other useful information.

res <- nbinomTest(cds, 'untreated', 'treated')
head(res)

#The interpretation of the columns of data.frame is as follows.
#id                    feature identifier
#baseMean     mean normalised counts, averaged over all samples from both conditions
#baseMeanA   mean normalised counts from condition A
#baseMeanB    mean normalised counts from condition B
#foldChange     fold change from condition A to B
#log2FoldChange   the logarithm (to basis 2) of the fold change
#pval          p value for the statistical signi_cance of this change
#padj       p value adjusted for multiple testing with the Benjamini-Hochberg procedure (see the R function
#p.adjust), which controls false discovery rate (FDR)


##Let us first plot the log2 fold changes against the mean normalised counts, colouring in red those genes that are significant at 10% FDR

plotMA(res)
##It is also instructive to look at the histogram of p values. The enrichment of low p values ##stems from
##the differentially expressed genes, while those not differentially expressed are spread uniformly over ##the range from zero
##to one (except for the p values from genes with very low counts, which take discrete values and so give ##rise to high
##counts for some bins at the right.)

hist(res$pval, breaks = 100, col = 'grey', border = 'slateblue', main = '')

##filtering for significant genes

resSig <- res[res$padj<0.1, ]
resSig

##list the most significant differentially expressed genes

head(resSig[order(resSig$pval), ])

##list the most strongly downregulated of significant genes

head(resSig[order(resSig$foldChange, -resSig$baseMean),])


##or at the most strongly upregulated significant genes

head(resSig[order(-resSig$foldChange, -resSig$baseMean),])


##saving the output to a file

write.csv(res, file = 'My pasilla analysis result table.csv')

##power to detect significantly expressed genes depends on expression strength, for weakly expressed genes, stronger changes need to be made to call significantly expressed genes

##focusing in untreated replicates

ncu <- counts(cds, normalized=TRUE)[, conditions(cds)=='untreated']

##ncu is a matrix with two columns

plotMA(data.frame(baseMean=rownames(ncu),
                  log2FoldChange=log2(ncu[,2]/ncu[,1])),
       col = 'black')
##As one can see in Figure 4, the log fold changes between replicates are stronger for lowly expressed ##genes than for
##highly expressed ones. We ought to conclude that a gene's expression is inuenced by the treatment ##only if the change
##between treated and untreated samples is stronger than what we see between replicates, and hence, ##the dividing line
##between red and black in Figure 2 mimics the shape seen in Figure 4.


##WORKING PARTIALLY WITHOUT REPLICATES

#If you have replicates for one condition and not for the other, then you can still proceed as before., only conditions with replicates will be used in estimation of dispersion

#subset the data to only three objects

cdsUUT <- cds[, 1:3]
pData(cdsUUT)


##now perform the analysis as before

cdsUUT <- estimateSizeFactors(cdsUUT)
cdsUUT <- estimateDispersions(cdsUUT)
resUUT <- nbinomTest(cdsUUT, 'untreated','treated')
##producing the analogous plot as before

plotMA(resUUT)
##Plot shows the same symmetry in up- and down-regulation as before, but a certain asymmetry in the ##boundary
##line for significance. This has an easy explanation: low counts suffer from proportionally stronger shot ##noise than high
##counts, and since there is only one \untreated" sample versus two \treated" ones, a stronger ##downward fold-change is
##required to be called significant than for the upward direction


##WORKING WITHOUT ANY REPLICATES

##Proper replicates are essential to interpret a biological experiment. After all, if one compares two ##conditions and finds
##a difference, how else can one know that this difference is due to the different conditions and would not ##have arisen
##between replicates, as well, just due to experimental or biological noise? Hence, any attempt to work ##without replicates
##will lead to conclusions of very limited reliability.

##seeing how this works, reduce count to just two columns
cds2 <- cds[, c('untreated3','treated3')]

#Now, without any replicates at all, the estimateDispersions function will refuse to proceed unless we #instruct it
#to ignore the condition labels and estimate the variance by treating all samples as if they were replicates #of the same
#condition:


cds2 <- estimateDispersions(cds2, method='blind', sharingMode='fit-only')

##Note the option sharingMode="fit-only". Remember that the default, sharingMode="maximum", takes ##care of
##outliers, i.e., genes with dispersion much larger than the fitted values. Without replicates, we cannot ##catch such
##outliers and so have to disable this functionality.


##finding differential gene expression

res2 <- nbinomTest(cds2, 'untreated','treated')

plotMA(res2) ##Unsurprisingly, we find much fewer hits, as can be seen from the plot

##MULTIFACTOR DESIGN

##Returning to the full pasilla dataset

head(pasillaCountTable)


pasillaDesign

##When creating a count data set with multiple factors, just pass a data frame instead of the condition factor:


cdsFull <- newCountDataSet(pasillaCountTable, pasillaDesign)
cdsFull <- estimateSizeFactors(cdsFull)
cdsFull <- estimateDispersions(cdsFull)

plotDispEsts(cdsFull)
##For inference, we now specify two models by formulas. The full model regresses the genes' expression ##on both the
##library type and the treatment condition, the reduced model regresses them only on the library type. ##For each gene,
##we fit generalized linear models (GLMs) according to the two models, and then compare them in order ##to infer whether
##the additional specification of the treatment improves the fit and hence, whether the treatment has ##significant effect.


fit1 <- fitNbinomGLMs(cdsFull, count ~ libType + condition)
fit0 <- fitNbinomGLMs(cdsFull, count ~ libType)


str(fit1)


##performing the test

pvalsGLM <- nbinomGLMTest(fit1, fit0)
padjGLM <- p.adjust(pvalsGLM, method = 'BH')

head(fit1)

##The first three columns show the fitted coefficients, converted to a logarithm base 2 scale. The log2 ##fold change
##due to the condition is shown in the third column. As indicated by the column name, it is the effect of ##\untreated",
##i.e., the log ratio of \untreated" versus \treated". (This is unfortunately the other way round as before, ##due to the
##peculiarities of contrast coding.) Note that the library type also had noticeable influence on the ##expression, and hence
##would have increased the dispersion estimates (and so reduced power) if we had not fitted an effect for ##it.
##The column deviance is the deviance of the _t. (Comparing the deviances with a _2 likelihood ratio test ##is how
##nbinomGLMTest calculates the p values.) The last column, converged, indicates whether the ##calculation of coefficients
##and deviance has fully converged. (If it is false too often, you can try to change the GLM control ##parameters, as
##explained in the help page to fitNbinomGLMs.) 


##Finally, we show that taking the library type into account was important to have good detection power ##by doing
##the analysis again using the standard workflow, as outlined earlier, and without informing DESeq of the ##library types

cdsFullB <- newCountDataSet(pasillaCountTable, pasillaDesign$condition)
cdsFullB <- estimateSizeFactors(cdsFullB)
cdsFullB <- estimateDispersions(cdsFullB)
resFullB <- nbinomTest(cdsFullB, 'untreated','treated')
plotMA(resFullB)


tab2 <- table(
  'all samples simple'=resFullB$padj<0.1,
  'all samples GLM'=padjGLM<0.1
)


addmargins(tab2)


##INDEPENDENT FILTERING AND MULTIPLE TESTING

##The analyses of the previous sections involve the application of statistical tests, one by one, to each ##row of the data
##set, in order to identify those genes that have evidence for differential expression. The idea of ##independent filtering is
##to filter out those tests from the procedure that have no, or little chance of showing significant evidence, ##without even
##looking at their test statistic. Typically, this results in increased detection power at the same ##experiment-wide type I
##error. Here, we measure experiment-wide type I error in terms of the false discovery rate. A good ##choice for a filtering
##criterion is one that;
#1. is statistically independent from the test statistic under the null hypothesis,
#2. is correlated with the test statistic under the alternative, and
#3. does not notably change the dependence structure {if there is any {between the test statistics of nulls #and
#alternatives.

rs <- rowSums(counts(cdsFull))
theta <- 0.4
use <- (rs <- quantile(rs, probs = theta))
table(use)
use
cdsFilt <- cdsFull[use, ]
##Above, we consider as a filter criterion rs, the overall sum of counts (irrespective of biological ##condition), and remove
##the genes in the lowest 40% quantile (as indicated by the parameter theta). We perform the testing as ##before

fitFilt1 <- fitNbinomGLMs(cdsFilt, count ~ libType + condition)

plot(rank(rs)/length(rs), -log10(pvalsGLM), pch=16, cex=0.45)


##Heat Map of the Count Table

##To explore a count table, it is often instructive to look at it as a heatmap. Below we show how to produce such a
##heatmap from the variance stabilisation transformed data for all 7 samples


cdsFullBlind <- estimateDispersions(cdsFull, method='blind')
vsdFull <- varianceStabilizingTransformation(cdsFullBlind)

library(RColorBrewer)
library(gplots)


select <- order(rowMeans(counts(cdsFull)), decreasing = TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)

heatmap.2(exprs(vsdFull)[select,], col=hmcol, trace = 'none', margin=c(10,6))


##for comparison, let us do the same with the untransformed counts

heatmap.2(counts(cdsFull)[select,],col=hmcol, trace = 'none', margin=c(10,6))


##Heat map of the sample-to-sample distances

##Another use of variance stabilized data is sample clustering. Here, we apply the dist function to the transpose of the transformed count matrix to get sample-to-sample distances.

dists <- dist(t(exprs(vsdFull)))


##a heat map of the distance matrix gives information on similarities and differencs between different samples

mat <- as.matrix(dists)
rownames(mat)=colnames(mat)=with(pData(cdsFullBlind), paste(condition, libType, sep = ':'))
heatmap.2(mat, trace = 'none', col=rev(hmcol), margin=c(13,13))

##The clustering correctly reects our experimental design, i.e., samples are more similar when they have ##the same
##treatment or the same library type. (To avoid potential circularities in this conclusion, it was important to ##re-estimate
##the dispersions with method="blind" in the calculation for cdsFullBlind above, as only then, the variance ##stabilizing
##transformation is not informed about the design, and we can be sure that it is not biased towards a ##result supporting
##the design.)


##PRINCIPAL COMPONENT PLOT OF THE SAMPLES

print(plotPCA(vsdFull, intgroup = c('condition','libType')))

