#############################################################
# Constants for use in all experiment analysis R scripts.
# 
# Include by:
#  source("../bin/constants.R")
#
# @author Michael C Rosenstein
# @version $Id$
#############################################################

designFile <- "design.dat"
sampleGroupColName <- "Group"

boxplotFile <- "boxplot_log2_norm.png"
normalizedDataFile <- "data_norm.txt"
pairwiseContrastsFile <- "all_pairwise_results.txt"
dataDelimiter <- "\t"
plotWidth <- 960
plotHeight <- 960

# Column names for the pairwiseConstrastsFile.
probeColName <- "ProbeID"
qValueColName <- "Qval"
