#!/usr/bin/env Rscript
#SBATCH --job-name='LOLA'
#SBATCH --output=/groups/gerlich/experiments/Experiments_005700/005752/LOLA/LOLA_CTCF/LOLA.out
#SBATCH -c 12
#SBATCH --mem 80G
#SBATCH --qos=short

library("LOLA")
library("simpleCache")
library("GenomicRanges")

# Lola Core

regionDB = loadRegionDB("/groups/gerlich/experiments/Experiments_004500/004506/LOLACore/hg19")
bedDir = "/groups/gerlich/experiments/Experiments_005700/005752/LOLA/LOLA_CTCF/"
outFolder = "/groups/gerlich/experiments/Experiments_005700/005752/LOLA/LOLA_CTCF/LOLA/"

binVector = list()
# load in bedfiles
binVector["TopTADs"] = readBed(paste(bedDir, "ctcf_top10.bed", sep=""))

# load in universe
allBins = readBed(paste(bedDir,"universe_ctcf.bed", sep=""))
# run analysis
userSets = GRangesList(binVector)
locResults = runLOLA(userSets, allBins, regionDB, cores=12)
# output results
writeCombinedEnrichment(locResults, outFolder=outFolder)


# LOLA EXT

regionDB = loadRegionDB("/groups/gerlich/experiments/Experiments_004500/004506/LOLAExt/hg19")
bedDir = "/groups/gerlich/experiments/Experiments_005700/005752/LOLA/LOLA_CTCF/"
outFolderEXT = "/groups/gerlich/experiments/Experiments_005700/005752/LOLA/LOLA_CTCF/LOLA_ext/"

binVector = list()
# load in bedfiles
binVector["TopTADs"] = readBed(paste(bedDir, "ctcf_top10.bed", sep=""))
#binVector["BottomTADs"] = readBed(paste(bedDir, "BottomTADsTransAmount.bed", sep=""))

# load in universe
allBins = readBed(paste(bedDir,"universe_ctcf.bed", sep=""))
# run analysis
userSets = GRangesList(binVector)
locResults = runLOLA(userSets, allBins, regionDB, cores=12)
# output results
writeCombinedEnrichment(locResults, outFolder=outFolderEXT)

