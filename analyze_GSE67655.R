#############################################################
# Analyzes the GSE67655 microarray experiment.
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67665
# 
#############################################################

#############################################################
# Constants
#############################################################

# How fold changes should be depicted: use -1 if earlier stages first; 1 otherwise.
pairContrastCoefficient <- -1
source ("./constants.R")

#############################################################
# Set up dependencies.
#############################################################


# To ensure proper library load, first setRepositories(), then select CRAN, BioC, BioC extra, CRAN (extras)
# Function to ensure packages always load on native script run:

# package.load <- function(pkg){
#   new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
#   if (length(new.pkg))
#     install.packages(new.pkg, dependencies = TRUE)
#   sapply(pkg, require, character.only = TRUE)
# }
# packages<- c("oligo","makecdfenv","maanova","qvalue")
# package.load(packages)

library(oligo)
library(makecdfenv)
library(maanova)
library(qvalue)
# Libraries below only necessary on the first run to create the custom cdf package
# library(pdInfoBuilder)
# library(ff)
# library(doMC)
# registerDoMC(3)



#############################################################
# Get CEL files from Gene Expression Omnibus
#############################################################

# download.file(url='http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE67665&format=file',
#               destfile='GSE67665_RAW.tar', method='curl')
# # Create RawFiles directory, untar and drop in
# dir.create("./RawFiles")
# untar("GSE67665_RAW.tar",exdir="./RawFiles")

# Copy RawFiles out, do a bunch of steps to make them match designFile downstream

files<- list.files(path="./RawFiles",pattern="*.gz",full.names=TRUE)
file.names<- list.files(path="./RawFiles",pattern="*.gz",full.names=FALSE)
file.newnames<-sapply(strsplit(file.names, split='_', fixed=TRUE), function(x) (x[1]))
file.newnames<-paste(file.newnames,".CEL.gz",sep="")
file.copy(files,".")
file.rename(from=list.files(pattern="*.gz"),to=file.newnames)

#############################################################
# Get necessary brainarray cdf to use with oligo
#############################################################


download.file("http://mbni.org/customcdf/20.0.0/enst.download/zebgene10st_Dr_ENST_20.0.0.zip","tmp.zip")
unzip("tmp.zip")
dir()
z <- cdf2table("zebgene10st_Dr_ENST.cdf")
seed <- new("GenericPDInfoPkgSeed",table=z,author="Kyle Scot Shank",email="kshank@mdibl.org", species = "Danio Rerio", pkgName = "pd.zebgene10st.dr.ENST.20")
makePdInfoPackage(seed)
install.packages("pd.zebgene10st.dr.ENST.20/",repos=NULL,type="source")

#############################################################
# Create Design Info
#############################################################

design<- data.frame(Array = sapply(strsplit(file.names, split='_', fixed=TRUE), function(x) (x[1])),
                     Group = c(rep("CTRL",3),rep("04HR",3),rep("DY01",3),rep("DY03",3),rep("DY07",3),rep("DY14",3),rep("DY90",3)),
                     Dye = rep(1,length(file.names)),
                     Sample=seq(1:length(file.names)))

write.table(design,"file.dat")

(nSampleGroups <- length(unique(design$Group)))

#############################################################
# Process raw data.
#############################################################
# Read all CEL files in cwd in directory, label them, read them in via read.celfiles()
celFiles <- list.celfiles(".",full.names=TRUE,listGzipped=TRUE)
celData <- read.celfiles(celFiles,pkgname="pd.zebgene10st.dr.ENST.20")

# # CDF (Chip Description File) from UMichigan BrainArray Lab.
# celData@cdfName <- "Zebrafish_Dr_ENST"
# Above is unnecessary as oligo package initiated download of pd.zebgene.1.0.st file

# Process CEL data using RMA, and quantile normalization across samples.
celDataRMA <- rma(celData, normalize=TRUE)
# class(celDataRMA)
# [1] "ExpressionSet"
# attr(,"package")
# [1] "Biobase"

# Export data from the ExpressionSet into a matrix.
celDataExn <- exprs(celDataRMA)

# Convert the matrix into a data.frame.
celDataExnDF <- data.frame(ProbeID=row.names(celDataExn), celDataExn)
(names(celDataExnDF) <- c(gsub("\\.CEL(\\.gz)?", "", names(celDataExnDF), ignore.case=TRUE)))
# head(celDataExnDF)

# Save the normalized exn data to file.
write.table(celDataExnDF, normalizedDataFile, sep=dataDelimiter, row.names=FALSE, quote=FALSE)

# Sample count == column count minus the probe col.
(nSamples <- ncol(celDataExnDF) - 1)


#############################################################
# Check design.
#############################################################
# Ensure design file sample count is the same as the number read.
stopifnot(all.equal(nSamples, nrow(design)))

# Ensure that samples line up between the design and the read values.
celDataCols <- colnames(celData)
for(ix in 1:nSamples){
  celVal <- unlist(strsplit(celDataCols[ix], "[.]"))[1]
  stopifnot(celVal == design$Array[ix])
}


#############################################################
# Create a boxplot of log2 'probe set' data.
#############################################################
png(boxplotFile)
boxplot(celDataExnDF[-1], main="Quant. Norm. Log2 Probe Set Intensities", 
        names=c(1:nSamples), xlab="Array", ylab="Intensity", col="red")
dev.off()


#############################################################
# Scatterplots of normalized data.
#############################################################
# Use all sample groups directly but remove ProbeID column.
celDataExnDataOnlyDF <- celDataExnDF[,c(1:nSamples+1)]
png("scatterplot_norm.png", width=plotWidth, height=plotHeight)
plot(celDataExnDataOnlyDF, main="Quant. Norm. Control, 4 Hours, 1, 3, 7, 14, and 90 dpa")
dev.off()


#############################################################
# MA plots.
#############################################################
# Calc avg exn values for each gene in each sample in each sample group.
# Creates vars: (matrix w/ single column) but really a 1-D array: 1 row w/ column for exn value of each gene
maCtrl <- (celDataExn[, 1] + celDataExn[, 2] + celDataExn[, 3])/3
ma04HR <- (celDataExn[, 4] + celDataExn[, 5] + celDataExn[, 6])/3
maDY01 <- (celDataExn[, 7] + celDataExn[, 8] + celDataExn[, 9])/3
maDY03 <- (celDataExn[, 10] + celDataExn[, 11] + celDataExn[, 12])/3
maDY07 <- (celDataExn[, 13] + celDataExn[, 14] + celDataExn[, 15])/3
maDY14 <- (celDataExn[, 16] + celDataExn[, 17] + celDataExn[, 18])/3
maDY90 <- (celDataExn[, 19] + celDataExn[, 20] + celDataExn[, 21])/3


# Plots:
# X: Average values for each gene in the entire array.
# Y: Fold Change.

# 04HR v ctrl
X <- (ma04HR + maCtrl)/2
Y <- (ma04HR - maCtrl)
png("MAplot_04hr_vs_ctrl.png")
smoothScatter(X, Y, main="04HR vs Control", xlab="(04HR + CTRL)/2", ylab="(04HR - CTRL)", sub=Sys.time())
dev.off()


# DY01 v ctrl
X <- (maDY01 + maCtrl)/2
Y <- (maDY01 - maCtrl)
png("MAplot_DY01_vs_ctrl.png")
smoothScatter(X, Y, main="DY01 vs Control", xlab="(DY01 + CTRL)/2", ylab="(DY01 - CTRL)", sub=Sys.time())
dev.off()

# DY03 v ctrl
X <- (maDY03 + maCtrl)/2
Y <- (maDY03 - maCtrl)
png("MAplot_DY03_vs_ctrl.png")
smoothScatter(X, Y, main="DY03 vs Control", xlab="(DY03 + CTRL)/2", ylab="(DY03 - CTRL)", sub=Sys.time())
dev.off()

# DY07 v ctrl
X <- (maDY07 + maCtrl)/2
Y <- (maDY07 - maCtrl)
png("MAplot_DY07_vs_ctrl.png")
smoothScatter(X, Y, main="DY07 vs Control", xlab="(DY07 + CTRL)/2", ylab="(DY07 - CTRL)", sub=Sys.time())

# DY14 v ctrl
X <- (maDY14 + maCtrl)/2
Y <- (maDY14 - maCtrl)
png("MAplot_DY14_vs_ctrl.png")
smoothScatter(X, Y, main="DY14 vs Control", xlab="(DY14 + CTRL)/2", ylab="(DY14 - CTRL)", sub=Sys.time())
dev.off()

# DY90 v ctrl
X <- (maDY90 + maCtrl)/2
Y <- (maDY90 - maCtrl)
png("MAplot_DY90_vs_ctrl.png")
smoothScatter(X, Y, main="DY90 vs Control", xlab="(DY90 + CTRL)/2", ylab="(DY190- CTRL)", sub=Sys.time())
dev.off()


#############################################################
# Read in data and design file so it can be read in by R/maanova
#############################################################
# design.dat maps samples to sample groups.
# uses std settings for Affymetric files: 
#   probeid is in first col
#   intensity values start at col 2 and go on to the right
#   don't need to log transform values (already done)
#   only need one-color data
maData <- read.madata(normalizedDataFile, design, probeid=1, intensity=2, 
                      arrayType="oneColor", log.tran=FALSE, spot=FALSE)


#############################################################
# Fit ANOVA model
#############################################################
# Using the design, applies formula for every gene.
#   Formula = '~' + sample-group column header in the design file.
fitMaanova <- fitmaanova(maData, formula=~Group, method="REML", verbose=TRUE, subCol=FALSE)
# Ignore warnings: "coercing argument of type 'double' to logical"


#############################################################
# Do all possible pairwise comparisons.
#############################################################
# Get the contrast matrix of all possible pairwise comparisons.
contrastMatrix <- pairContrastCoefficient * PairContrast(nSampleGroups)

# Analyze intensities using the analysis of variance model (defined above).
# Computes all possible comparisons b/w all our sample groups
#   Output has all the stats we need: fold changes and q-values
#   'perm' param sets no. of permutations, more creates more accurate pvalues, but takes longer.
maTest <- matest(maData, fitMaanova, Contrast=contrastMatrix, term=sampleGroupColName, 
                 n.perm=1000, MME.method="REML", shuffle.method="sample", test.type="ttest")
# Ignore warnings: In any(parsed.formula$random) : coercing argument of type 'double' to logical

maTest <- adjPval(maTest, method='adaptive')

summary <- summarytable(maTest)
# Ignore warnings; In summarytable.ttest(matestobj, smethod, test, whichTest, threshold,:  No whichTest infomation.

qvalue <- qvalue(maTest$Fs$Pvalperm)

# A check to ensure order is the same (should == TRUE)
table(qvalue$pval == maTest$Fs$Pvalperm)

### qsummary(qvalue)
### Function not working: possibly no longer supported? Below works
summary(qvalue)

allResults <- data.frame(summary, qvalue$qval)

# Name the qvalue columns.
(resultsNamedColCount <- grep("^X1$", names(allResults)) - 1)
(names(allResults) <- c(
  names(allResults)[1:resultsNamedColCount],
  gsub("Fold\\.change", qValueColName, grep("Fold\\.change", names(allResults), value=TRUE))
  ))


#############################################################
# Annotate results.
#############################################################

##### Create all_paireiwse_results.txt file:

# Funky stuff to get the first column of ProbeIDs.
annotation <- data.frame(rownames(allResults))
(names(annotation) <- c(probeColName))
allResultsAnno <- merge(allResults, annotation, by.x=0, by.y=probeColName, all.x=TRUE)
colnames(allResultsAnno)[1] <- probeColName
names(allResultsAnno)
# Fix ProbeID Names
allResultsAnno$ProbeID<-strtrim(allResultsAnno$ProbeID,nchar(allResultsAnno$ProbeID)-3)

# Error Check for Annotations
stopifnot(all.equal(nrow(annotation), nrow(celDataExnDF)))

# Merge with Ensembl Biomart Output
ensembl_output <- read.csv(file="mart_export.txt",head=TRUE,stringsAsFactors=FALSE)
allResultsAnno<- merge(allResultsAnno, ensembl_output, by.x="ProbeID", by.y="Ensembl.Transcript.ID", all.x=TRUE)

write.table(allResultsAnno, pairwiseContrastsFile, sep=dataDelimiter, quote=FALSE, row.names=FALSE)

#############################################################
# Output Annotated Result 
#############################################################

raw.csv<-read.csv(file="summarytable.csv",head=T,stringsAsFactors=FALSE)
colnames(raw.csv)[1]<-"ProbeID"
raw.csv$ProbeID<-strtrim(raw.csv$ProbeID,nchar(raw.csv$ProbeID)-3)
update.csv<-merge(raw.csv, ensembl_output, by.x="ProbeID", by.y="Ensembl.Transcript.ID", all.x=TRUE)
write.csv(update.csv,file="summarytable.csv",row.names=FALSE)

# Output Each Pairwise Group

X04HR.Ctrl <- subset(update.csv,select=c(1,2:4,65:72))
X04HR.DY01 <- subset(update.csv,select=c(1,5:7,65:72))
X04HR.DY03 <- subset(update.csv,select=c(1,8:10,65:72))
X04HR.DY07 <- subset(update.csv,select=c(1,11:13,65:72))
X04HR.DY14 <- subset(update.csv,select=c(1,14:16,65:72))
X04HR.DY90 <- subset(update.csv,select=c(1,17:19,65:72))

CTRL.DY01 <- subset(update.csv,select=c(1,20:22,65:72))
CTRL.DY03 <- subset(update.csv,select=c(1,23:25,65:72))
CTRL.DY07 <- subset(update.csv,select=c(1,26:28,65:72))
CTRL.DY14 <- subset(update.csv,select=c(1,29:31,65:72))
CTRL.DY90 <- subset(update.csv,select=c(1,32:34,65:72))

DY01.DY03 <- subset(update.csv,select=c(1,35:37,65:72))
DY01.DY07 <- subset(update.csv,select=c(1,38:40,65:72))
DY01.DY14 <- subset(update.csv,select=c(1,41:43,65:72))
DY01.DY90 <- subset(update.csv,select=c(1,44:46,65:72))

DY03.DY07 <- subset(update.csv,select=c(1,47:49,65:72))
DY03.DY14 <- subset(update.csv,select=c(1,50:52,65:72))
DY03.DY90 <- subset(update.csv,select=c(1,53:55,65:72))

DY07.DY14 <- subset(update.csv,select=c(1,56:58,65:72))
DY07.DY90 <- subset(update.csv,select=c(1,59:61,65:72))

DY14.DY90 <- subset(update.csv,select=c(1,62:64,65:72))