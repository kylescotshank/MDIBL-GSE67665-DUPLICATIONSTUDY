#############################################################
# Analyzes the GSE67655 microarray experiment.
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67665
# 
# @Author: Kyle Shank - Bioinformatics Training Specialist, 
#                       MDI Biological Laboratory
#############################################################

#############################################################
# Constants
#############################################################

# How fold changes should be depicted: use -1 if earlier stages first; 1 otherwise.
pairContrastCoefficient <- -1
# constants.R is for re-use throughout RegenDB
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

download.file(url='http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE67665&format=file',
              destfile='GSE67665_RAW.tar', method='curl')
# Create RawFiles directory, untar and drop in
dir.create("./RawFiles")
untar("GSE67665_RAW.tar",exdir="./RawFiles")

# Copy RawFiles out, do a bunch of steps to make them match designFile downstream

files<- list.files(path="./RawFiles",pattern="*.gz",full.names=TRUE)
file.names<- list.files(path="./RawFiles",pattern="*.gz",full.names=FALSE)
file.newnames<-sapply(strsplit(file.names, split='_', fixed=TRUE), function(x) (x[1]))
file.newnames<-paste(file.newnames,".CEL.gz",sep="")
file.copy(files,".")
file.rename(from=list.files(pattern="*.gz"),to=file.newnames)

#############################################################
# Get necessary brainarray cdf to use with oligo, create local
# package source to use.
#
# only need to run commented-out code the first time through
#############################################################


# download.file("http://mbni.org/customcdf/20.0.0/enst.download/zebgene10st_Dr_ENST_20.0.0.zip","tmp.zip")
# unzip("tmp.zip")
# dir()
# z <- cdf2table("zebgene10st_Dr_ENST.cdf")
# seed <- new("GenericPDInfoPkgSeed",table=z,author="me",email="me@theinternet.com", species = "Danio Rerio", pkgName = "pd.zebgene10st.dr.ENST.20")
# makePdInfoPackage(seed)
install.packages("pd.zebgene10st.dr.ENST.20/",repos=NULL,type="source")

#############################################################
# Create Design Info
#############################################################

design<- data.frame(Array = sapply(strsplit(file.names, split='_', fixed=TRUE), function(x) (x[1])),
                     Group = c(rep("DPI00_CTRL",3),rep("DPI00_FOURHPI",3),rep("DPI01",3),rep("DPI03",3),rep("DPI07",3),rep("DPI14",3),rep("DPI90",3)),
                     Dye = rep(1,length(file.names)),
                     Sample=seq(1:length(file.names)))

write.table(design,"design.dat",sep="\t",row.names=FALSE)

# Check number of sample groups - should be 7
(nSampleGroups <- length(unique(design$Group)))

#############################################################
# Process raw data.
#############################################################
# Read all CEL files in cwd in directory, label them, read them in via read.celfiles()
# pckname is to ensure that latest brainarray version is used
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
maCTRL <- (celDataExn[, 1] + celDataExn[, 2] + celDataExn[, 3])/3
ma4HPI <- (celDataExn[, 4] + celDataExn[, 5] + celDataExn[, 6])/3
ma1DPI <- (celDataExn[, 7] + celDataExn[, 8] + celDataExn[, 9])/3
ma3DPI <- (celDataExn[, 10] + celDataExn[, 11] + celDataExn[, 12])/3
ma7DPI <- (celDataExn[, 13] + celDataExn[, 14] + celDataExn[, 15])/3
ma14DPI <- (celDataExn[, 16] + celDataExn[, 17] + celDataExn[, 18])/3
ma90DPI <- (celDataExn[, 19] + celDataExn[, 20] + celDataExn[, 21])/3


# Plots:
# X: Average values for each gene in the entire array.
# Y: Fold Change.

# 4HPI v ctrl
X <- (ma4HPI + maCTRL)/2
Y <- (ma4HPI - maCTRL)
png("MAplot_4HPI_vs_CTRL.png")
smoothScatter(X, Y, main="4HPI vs Control", xlab="(4HPI + CTRL)/2", ylab="(4HPI - CTRL)", sub=Sys.time())
dev.off()


# 1DPI v CTRL
X <- (ma1DPI + maCTRL)/2
Y <- (ma1DPI - maCTRL)
png("MAplot_1DPI_vs_CTRL.png")
smoothScatter(X, Y, main="1DPI vs Control", xlab="(1DPI + CTRL)/2", ylab="(1DPI - CTRL)", sub=Sys.time())
dev.off()

# 3DPI v CTRL
X <- (ma3DPI + maCTRL)/2
Y <- (ma3DPI - maCTRL)
png("MAplot_3DPI_vs_CTRL.png")
smoothScatter(X, Y, main="3DPI vs Control", xlab="(3DPI + CTRL)/2", ylab="(3DPI - CTRL)", sub=Sys.time())
dev.off()

# 7DPI v CTRL
X <- (ma7DPI + maCTRL)/2
Y <- (ma7DPI - maCTRL)
png("MAplot_7DPI_vs_CTRL.png")
smoothScatter(X, Y, main="7DPI vs Control", xlab="(7DPI + CTRL)/2", ylab="(7DPI - CTRL)", sub=Sys.time())

# 14DPI v CTRL
X <- (ma14DPI + maCTRL)/2
Y <- (ma14DPI - maCTRL)
png("MAplot_14DPI_vs_CTRL.png")
smoothScatter(X, Y, main="DY14 vs Control", xlab="(14DPI + CTRL)/2", ylab="(14DPI - CTRL)", sub=Sys.time())
dev.off()

# 90DPI v CTRL
X <- (ma90DPI + maCTRL)/2
Y <- (ma90DPI - maCTRL)
png("MAplot_90DPI_vs_CTRL.png")
smoothScatter(X, Y, main="90DPI vs Control", xlab="(90DPI + CTRL)/2", ylab="(90DPI- CTRL)", sub=Sys.time())
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
colnames(allResultsAnno)[128:135]<-c("EnsemblGeneID","AssociatedGeneName","Description","GeneType","ChromsomeName","GeneStartBP","GeneEndBP","Strand")

write.table(allResultsAnno, pairwiseContrastsFile, sep=dataDelimiter, quote=FALSE, row.names=FALSE)

#############################################################
# Output Annotated Result 
#############################################################

update.csv<-allResultsAnno
write.csv(update.csv,file="summarytable.csv",row.names=FALSE)

# Output Each Pairwise Group

DPI00_CTRL.DPI00_FOURHPI <- subset(update.csv,select=c(1,2:6,107,128:135))
write.csv(DPI00_CTRL.DPI00_FOURHPI,file="CTRL_v_4HPI.csv",row.names=FALSE)
DPI00_CTRL.DPI01 <- subset(update.csv,select=c(1,7:11,108,128:135))
write.csv(DPI00_CTRL.DPI01,file="CTRL_v_1DPI.csv",row.names=FALSE)
DPI00_CTRL.DPI03 <- subset(update.csv,select=c(1,12:16,109,128:135))
write.csv(DPI00_CTRL.DPI03,file="CTRL_v_3DPI.csv",row.names=FALSE)
DPI00_CTRL.DPI07 <- subset(update.csv,select=c(1,17:21,110,128:135))
write.csv(DPI00_CTRL.DPI07,file="CTRL_v_7DPI.csv",row.names=FALSE)
DPI00_CTRL.DPI14 <- subset(update.csv,select=c(1,22:26,111,128:135))
write.csv(DPI00_CTRL.DPI14,file="CTRL_v_14DPI.csv",row.names=FALSE)
DPI00_CTRL.DPI90 <- subset(update.csv,select=c(1,27:31,112,128:135))
write.csv(DPI00_CTRL.DPI90,file="CTRL_v_90DPI.csv",row.names=FALSE)

DPI00_FOURHPI.DPI01 <- subset(update.csv,select=c(1,32:36,113,128:135))
write.csv(DPI00_FOURHPI.DPI01 ,file="4HPI_v_1DPI.csv",row.names=FALSE)
DPI00_FOURHPI.DPI03 <- subset(update.csv,select=c(1,37:41,114,128:135))
write.csv(DPI00_FOURHPI.DPI03,file="4HPI_v_3DPI.csv",row.names=FALSE)
DPI00_FOURHPI.DPI07 <- subset(update.csv,select=c(1,42:46,115,128:135))
write.csv(DPI00_FOURHPI.DPI07,file="4HPI_v_7DPI.csv",row.names=FALSE)
DPI00_FOURHPI.DPI14 <- subset(update.csv,select=c(1,47:51,116,128:135))
write.csv(DPI00_FOURHPI.DPI14,file="4HPI_v_14DPI.csv",row.names=FALSE)
DPI00_FOURHPI.DPI90 <- subset(update.csv,select=c(1,52:56,117,128:135))
write.csv(DPI00_FOURHPI.DPI90,file="4HPI_v_90DPI.csv",row.names=FALSE)

DPI01.DPI03 <- subset(update.csv,select=c(1,57:61,118,128:135))
write.csv(DPI01.DPI03,file="1DPI_v_3DPI.csv",row.names=FALSE)
DPI01.DPI07 <- subset(update.csv,select=c(1,62:66,119,128:135))
write.csv(DPI01.DPI07,file="1DPI_v_7DPI.csv",row.names=FALSE)
DPI01.DPI14 <- subset(update.csv,select=c(1,67:71,120,128:135))
write.csv(DPI01.DPI14,file="1DPI_v_14DPI.csv",row.names=FALSE)
DPI01.DPI90 <- subset(update.csv,select=c(1,72:76,121,128:135))
write.csv(DPI01.DPI90,file="1DPI_v_90DPI.csv",row.names=FALSE)

DPI03.DPI07 <- subset(update.csv,select=c(1,77:81,122,128:135))
write.csv(DPI03.DPI07,file="3DPI_v_7DPI.csv",row.names=FALSE)
DPI03.DPI14 <- subset(update.csv,select=c(1,82:86,123,128:135))
write.csv(DPI03.DPI14,file="3DPI_v_14DPI.csv",row.names=FALSE)
DPI03.DPI90 <- subset(update.csv,select=c(1,87:91,124,128:135))
write.csv(DPI03.DPI90,file="3DPI_v_90DPI.csv",row.names=FALSE)

DPI07.14DPI <- subset(update.csv,select=c(1,92:96,125,128:135))
write.csv(DPI07.14DPI,file="7DPI_v_14DPI.csv",row.names=FALSE)
DPI07.90DPI <- subset(update.csv,select=c(1,97:101,126,128:135))
write.csv(DPI07.90DPI,file="7DPI_v_90DPI.csv",row.names=FALSE)

DPI14.DPI90 <- subset(update.csv,select=c(1,102:106,127,128:135))
write.csv(DPI14.DPI90,file="14DPI_v_90DPI.csv",row.names=FALSE)
# 
# 
# #############################################################
# # Fetch .sdrf and .idf files (if necessary)
# #############################################################
# # download.file("http://www.ebi.ac.uk/arrayexpress/files/E-GEOD-67665/E-GEOD-67665.sdrf.txt","E-GEOD-67665.sdrf.txt")
# # download.file("http://www.ebi.ac.uk/arrayexpress/files/E-GEOD-67665/E-GEOD-67665.idf.txt", "E-GEOD-67665.idf.txt")


datanorm<-read.table(file="data_norm.txt",header=T)
