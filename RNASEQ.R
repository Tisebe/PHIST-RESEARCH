install.packages('Rsubread')
install.packages('limma')
install.packages('tidyverse')



BMI <- data.frame(
  gender =c('Male','Male','Female'),
  height =c(152,171.5,165),
  weight =c(81,93,78),
  age =c(42,38, 26))
print(BMI)

var.1 = c(0,1,2,3)
var.2 <- c('learn', 'R')
c(TRUE,1) -> var.3
var.1

var_x <- 'Hello'
cat('The class of var_x is', class(var_x),'\n')

ls()

v <- c(2,5.5,6)
t <- c(8,3,4)
print(v+t)

v <- c(2,5.5,6)
t <- c(8,3,4)
print(v/t)


v <- c(2,5.5,6)
t <- c(8,3,4)
print(v%%t)

v <- c(2,5.5,6,9)
t <- c(8,2.5,14,9)
print(v>t)

v <- c(2,5.5,6,9)
t <- c(8,2.5,14,9)
print(v<t)

v <- c(2,5.5,6,9)
t <- c(8,2.5,14,9)
print(v==t)


v <- c(2,5.5,6,9)
t <- c(8,2.5,14,9)
print(v!=t)

v  <- c(3,1,TRUE,2+3i)
t <- c(4,1,FALSE,2+3i)
print(v&t)

v  <- c(3,1,TRUE,2+3i)
t <- c(4,1,FALSE,2+3i)
print(v!t)

v <- c(3,0,TRUE,2+2!)
t <- c(1,3, TRUE, 2+3!)
print(v&&t)

x <- 30L
if(is.integer(x)) {
  print('X is an Integer')
}

x <- c('what', 'is','truth')
if('Truth'%in% x){print('Truth is found')} else {print('Truth is not found')}


v <- c('Hello','loop')
cnt <-2
repeat {
  print(v)
  cnt <- cnt+1
if(cnt > 5){
  break
}
}


a <- "Hello"
b <- 'How'
c <- "are you"

print(paste(a,b,c))
print(paste(a,b,c, sep = "-"))
print(paste(a,b,c, sep = "", collapse = ""))

print(seq(5,9, by=0.4))

t <- c('sun','mon','tue','wed','thurs','fri','sat')
u <- t[c(2,3,6)]
print(u)

M <- matrix(c(3:14), nrow = 4, byrow = TRUE)
print(M)

vector1 <- c(5,9,3)
vector2 <- c(10,11,12,13,14,15)

result <- array(c(vector1,vector2),dim = c(3,3,2))
print(result)
vector1 <- c(5,9,3)
vector2 <- c(10,11,12,13,14,15)

column.names <- c('COL1','COL2','COL3')
row.names <- c('ROW1','ROW2','ROW3')
matrix.names <- c('Matrix1','Matrix2')

result <- array(c(vector1,vector2),dim = c(3,3,2),dimnames = list(row.names,column.names,matrix.names))
print(result)


data <- c('East','West','East','North','East','North','West','East','West','West','West','East','North')
print(data)

print(is.factor(data))

factor_data <- factor(data)
print(factor_data)
print(is.factor(factor_data))

height <- c(132,151,162,139,166,147,122)
weight <- c(48,49,66,53,67,52,40)
gender <- c('male','male','female','female','male','female','male')

input_data <- data.frame(height,weight,gender)
print(input_data)

print(is.factor(input_data$gender))
print(input_data$gender)


new_order_data <- factor(factor_data,levels = c('East','West','North'))
print(new_order_data)


emp.data <- data.frame(
  emp_id = c(1:5),
  emp_name = c('Rick','Dan','Michelle','Ryan','Gary'),
  salary = c(623.3,515.2,611.0,729.0,843.25),
  
  start_date = as.Date(c('2012-01-01','2013-09-23','2014-11-15','2014-05-11','2015-03-27')),
  stringsAsFactors = FALSE
)
print(emp.data)

emp.data <- data.frame(
  emp_id = c(1:5),
  emp_name = c('Rick','Dan','Michelle','Ryan','Gary'),
  salary = c(623.3,515.2,611.0,729.0,843.25),
  
  start_date = as.Date(c('2012-01-01','2013-09-23','2014-11-15','2014-05-11','2015-03-27')),
  stringsAsFactors = FALSE
)
print(summary(emp.data))

result <- data.frame(emp.data$emp_name,emp.data$salary)
print(result)

result <- emp.data[1:2,]
print(result)

city <- c('Tampa','Seattle','Hartford','Denver')
state <- c('FL','WA','CT','CO')
zipcode <- c(33602,98104,06161,80294)

addresses <- cbind(city,state,zipcode)
cat('#### The First data frame\n')

print(addresses)

new.address <- data.frame(
  city = c('Lowry','Charlotte'),
  state = c('CO','FL'),
  zipcode = c('80230','33949'),
  stringsAsFactors = FALSE
)
cat('### The Second data frame\n')
print(new.address)

all.addresses <- rbind(addresses,new.address)
cat('### The combined data frame\n')
print(all.addresses)
library(MASS)
merged.Pima <- merge(x=Pima.te, y= Pima.tr,
                     by.x = c('bp','bmi'),
                     by.y = c('bp','bmi')
                     )
print(merged.Pima)
nrow(merged.Pima)

library(MASS)
print(ships)
molten.ships <- melt(ships,id = c('type','year'))
print(molten.ships)

molten.ships <- molten(ships, id = c("type","year"))
print(molten.ships)



print(getwd())
setwd('/web/com')
print(getwd())

data <- read.csv('input.csv')
print(data)
print(is.data.frame(data))
print(ncol(data))
print(nrow(data))

data <- read.csv('input.csv')
sal <- max(data$salary)

data <- read.csv('input.csv')
sal <- max(data$salary)
retval <- subset(data, salary == max(salary))
print(retval)

retval <- subset(data, as.Date(start_date) > as.Date('2014-01-01'))
print(retval)

data <- read.csv('input.csv')
retval <- subset(data, as.Date(start_date) > as.Date('2014-01-01'))
write.csv(retval,'output.csv', row.names = FALSE)
newdata <- read.csv('output.csv')
print(newdata)
 install.packages('xlsx')
any(grepl('xlsx',installed.packages()))
library("xlsx")

fun(libname, pkgname)

install.packages('xlsxjars')
library("xlsx")
install.packages("xlsx")
library("xlsx")
data <- read.xls("input.xlsx", sheetIndex = 1)
print(data)
install.packages("XML")

library('XML')
library('methods')
result <- xmlParse(file = "input.xml")
print(result)

rootnode <- xmlRoot(result)
rootsize <- xmlSize(rootnode)

print(rootsize)
xmldataframe <- xmlToDataFrame("input.xml")
print(xmldataframe)

install.packages('rjson')


library('rjson')
result <- fromJSON(file = 'input.json')
print(result)
json_data_frame <- as.data.frame(result)
print(json_data_frame)

install.packages('RCurl')
install.packages('XML')
install.packages('plyr')

url <- 'http://www.geos.ed.ac.uk/~weather/jcmb_ws/'
links <- getHTMLLinks(links, 'JCMB_2015')


install.packages('RMySQL')
mysqlconnection = dbConnect(MySQL(), user = 'root', password = 'sakila', host = 'localhost')

mysqlconnection = dbConnect(MySQL(), user = 'root', password = '', dbname = 'sakila',
  
pie(x,labels,radius,main,col,clockwise)                                                   host = 'localhost')
x <- c(21,62,10,53)
labels <- c('London','New York', 'Singapore','Mumbai')
png(file = 'city.jpg')
pie (x,labels)
dev.off()

x <- c(21, 62, 10, 53)
labels <- c("London", "New York", "Singapore", "Mumbai")

# Give the chart file a name.
png(file = "city.jpg")

# Plot the chart.
pie(x,labels)

# Save the file.
dev.off()

H <- c(7,12,28,3,41)
png(file = 'barchart.png')
barplot(H)
dev.off()

for f in `cat files`; do STAR --genomeDir ../STAR/ENSEMBL.homo_sapiens.release-75 \
--readFilesIn fastq/$f\_1.fastq fastq/$f\_2.fastq \
--runThreadN 12 --outFileNamePrefix aligned/$f.; done

library("airway")
indir <- system.file("extdata", package="airway", mustWork=TRUE)

indir <- system.file("extdata", package="airway", mustWork=TRUE)
list.files(indir)

csvfile <- file.path(indir, "sample_table.csv")
sampleTable <- read.csv(csvfile, row.names = 1)
sampleTable

filenames <- file.path(indir, paste0(sampleTable$Run, "_subset.bam"))
file.exists(filenames)

bamfiles <- BamFileList(filenames, yieldSize=2000000)library("Rsamtools")



library("Rsamtools")
bamfiles <- BamFileList(filenames, yieldSize=2000000)

seqinfo(bamfiles[1])

library("GenomicFeatures")
gtffile <- file.path(indir,"Homo_sapiens.GRCh37.75_subset.gtf")
txdb <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character())
txdb

ebg <- exonsBy(txdb, by="gene")
ebg

library("GenomicAlignments")
library("BiocParallel")

register(SerialParam())
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )


se
dim(se)
assayNames(se)
head(assay(se), 3)
colSums(assay(se))
rowRanges(se)

str(metadata(rowRanges(se)))
colData(se)
colData(se) <- DataFrame(sampleTable)
colData(se)

se$cell
se$dex
library("magrittr")
se$dex %<>% relevel("untrt")
se$dex
se$dex <- relevel(se$dex, "untrt")

vignette('airway')

data("airway")
se <- airway

se$dex %<>% relevel("untrt")
se$dex
round( colSums(assay(se)) / 1e6, 1 )
colData(se)
library("DESeq2")
dds <- DESeqDataSet(se, design = ~ cell + dex)

countdata <- assay(se)
head(countdata, 3)
coldata <- colData(se)
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                                         design = ~ cell + dex)


nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

lambda <- 10^seq(from = -1, to = 2, length = 1000)
cts <- matrix(rpois(1000*100, lambda), ncol = 100)
library("vsn")
meanSdPlot(cts, ranks = FALSE)

install.packages('hexbin')

log.cts.one <- log2(cts + 1)
meanSdPlot(log.cts.one, ranks = FALSE)
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)


colData(vsd)
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)
library("dplyr")
library("ggplot2")


dds <- estimateSizeFactors(dds)
df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))


colnames(df)[1:2] <- c("x", "y")

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

sampleDists <- dist(t(assay(vsd)))
sampleDists
install.packages("pheatmap")
library("RColorBrewer")

library("pheatmap")
library("RColorBrewer")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))

samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$dex, dds$cell, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

plotPCA(vsd, intgroup = c("dex", "cell"))
pcaData <- plotPCA(vsd, intgroup = c( "dex", "cell"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = dex, shape = cell)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()

mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = dex, shape = cell)) +
  geom_point(size = 3) + coord_fixed()

mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = dex, shape = cell)) +
  geom_point(size = 3) + coord_fixed()

dds <- DESeq(dds)
res <- results(dds)
res
res <- results(dds, contrast=c("dex","trt","untrt"))
mcols(res, use.names = TRUE)

summary(res)
res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)
results(dds, contrast = c("cell", "N061011", "N61311"))
sum(res$pvalue < 0.05, na.rm=TRUE)
sum(res$padj < 0.1, na.rm=TRUE)
resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])

topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("dex"))

library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("dex","cell"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = cell)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)


ggplot(geneCounts, aes(x = dex, y = count, color = cell, group = cell)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()


library("apeglm")
resultsNames(dds)

res <- lfcShrink(dds, coef="dex_trt_vs_untrt", type="apeglm")
plotMA(res, ylim = c(-5, 5))

res.noshr <- results(dds, name="dex_trt_vs_untrt")
plotMA(res.noshr, ylim = c(-5, 5))

plotMA(res, ylim = c(-5,5))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})


hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)

mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("cell","dex")])
pheatmap(mat, annotation_col = anno)

qs <- c(0, quantile(resLFC1$baseMean[resLFC1$baseMean > 0], 0:6/6))
bins <- cut(resLFC1$baseMean, qs)
levels(bins) <- paste0("~", round(signif((qs[-1] + qs[-length(qs)])/2, 2)))
fractionSig <- tapply(resLFC1$pvalue, bins, function(p)
  mean(p < .05, na.rm = TRUE))
barplot(fractionSig, xlab = "mean normalized count",
        ylab = "fraction of small p values")

library("AnnotationDbi")
library("org.Hs.eg.db")

columns(org.Hs.eg.db)

res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

resOrdered <- res[order(res$pvalue),]
head(resOrdered)

resOrderedDF <- as.data.frame(resOrdered)[1:100, ]
write.csv(resOrderedDF, file = "results.csv")
library("ReportingTools")
htmlRep <- HTMLReport(shortName="report", title="My report",
                      reportDirectory="./report")
publish(resOrderedDF, htmlRep)
url <- finish(htmlRep)
browseURL(url)

resGR <- results(dds, name="dex_trt_vs_untrt", format="GRanges")
resGR$log2FoldChange <- res$log2FoldChange
resGR

resGR$symbol <- mapIds(org.Hs.eg.db, names(resGR), "SYMBOL", "ENSEMBL")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("devtools")

library("Gviz")

window <- resGR[topGene] + 1e6
strand(window) <- "*"
resGRsub <- resGR[resGR %over% window]
naOrDup <- is.na(resGRsub$symbol) | duplicated(resGRsub$symbol)
resGRsub$group <- ifelse(naOrDup, names(resGRsub), resGRsub$symbol)
status <- factor(ifelse(resGRsub$padj < 0.1 & !is.na(resGRsub$padj),
                        "sig", "notsig"))

options(ucscChromosomeNames = FALSE)
g <- GenomeAxisTrack()
a <- AnnotationTrack(resGRsub, name = "gene ranges", feature = status)
d <- DataTrack(resGRsub, data = "log2FoldChange", baseline = 0,
               type = "h", name = "log2 fold change", strand = "+")
plotTracks(list(g, d, a), groupAnnotation = "group",
           notsig = "grey", sig = "hotpink")



library("sva")
dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ dex, colData(dds))
mod0 <- model.matrix(~   1, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv = 2)

svseq$sv

par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:2) {
  stripchart(svseq$sv[, i] ~ dds$cell, vertical = TRUE, main = paste0("SV", i))
  abline(h = 0)
}

ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV1 + SV2 + dex



library("RUVSeq")
set <- newSeqExpressionSet(counts(dds))
idx  <- rowSums(counts(set) > 5) >= 2
set  <- set[idx, ]
set <- betweenLaneNormalization(set, which="upper")
not.sig <- rownames(res)[which(res$pvalue > .1)]
empirical <- rownames(set)[ rownames(set) %in% not.sig ]
set <- RUVg(set, empirical, k=2)
pData(set)

par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:2) {
  stripchart(pData(set)[, i] ~ dds$cell, vertical = TRUE, main = paste0("W", i))
  abline(h = 0)
}

ddsruv <- dds
ddsruv$W1 <- set$W_1
ddsruv$W2 <- set$W_2
design(ddsruv) <- ~ W1 + W2 + dex


library("fission")
data("fission")
ddsTC <- DESeqDataSet(fission, ~ strain + minute + strain:minute)

ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ strain + minute)
resTC <- results(ddsTC)
resTC$symbol <- mcols(ddsTC)$symbol
head(resTC[order(resTC$padj),], 4)
fiss <- plotCounts(ddsTC, which.min(resTC$padj), 
                   intgroup = c("minute","strain"), returnData = TRUE)
fiss$minute <- as.numeric(as.character(fiss$minute))
ggplot(fiss,
       aes(x = minute, y = count, color = strain, group = strain)) + 
  geom_point() + stat_summary(fun.y=mean, geom="line") +
  scale_y_log10()
resultsNames(ddsTC)
res30 <- results(ddsTC, name="strainmut.minute30", test="Wald")
res30[which.min(resTC$padj),]
betas <- coef(ddsTC)
colnames(betas)

topGenes <- head(order(resTC$padj),20)
mat <- betas[topGenes, -c(1,2)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE)


devtools::session_info()











