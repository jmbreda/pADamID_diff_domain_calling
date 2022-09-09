# Load dependencies
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(knitr)

library("optparse")
library("tools")
library(HMMt)

setwd('~/Workspace/jb20220118_Stefano_pADamID')
# Prepare output
output.dir <- "ts200406_domain_calling_DRB_Top1"
dir.create(output.dir, showWarnings = FALSE)

# Input tracks
input.dir <- "ts200326_differential_analysis_regions/bigwig"

diff_mean_DRB <- import(file.path(input.dir, "diff_mean_DRB.bw"))
diff_mean_Top1 <- import(file.path(input.dir, "diff_mean_Top1.bw"))
diff_smooth_mean_DRB <- import(file.path(input.dir, "diff_smooth_mean_DRB.bw"))
diff_smooth_mean_Top1 <- import(file.path(input.dir, "diff_smooth_mean_Top1.bw"))

# Input example normalized bins
input.dir <- "results/normalized/bin-10kb/"
bins_example <- read_tsv(file.path(input.dir,
                                   "pADamID-RPE_170_CT_r1_Lmnb2-10kb.norm.txt.gz"),
                         col_names = c("seqnames", "start", "end", "score"))
bins_example$start <- bins_example$start + 1


opts_chunk$set(fig.width = 10, fig.height = 4, 
               dev=c('png', 'pdf'), fig.path = file.path(output.dir, "figures/")) 
pdf.options(useDingbats = FALSE)




## HMM calling
# Write data files
HMM.dir <- file.path(output.dir, "HMM")
dir.create(HMM.dir, showWarnings = F)

WriteGRToTable <- function(gr, example, dir, name) {
  
  # Convert to tibble
  tib <- as_tibble(gr) %>%
    dplyr::select(seqnames, start, end, score)
  
  # Match with example
  idx <- match(paste(example$seqnames, example$start),
               paste(tib$seqnames, tib$start))
  
  # Change score of example
  example$score <- tib$score[idx]
  example$start <- example$start - 1
  
  # Write file
  f_name <- file.path(dir, paste0(name, ".norm.txt.gz"))
  #gz1 <- gzfile(f_name, "w")
  #write_tsv(example, gz1, col_names = F)
  #close(gz1)
  gz1 <- gzfile(f_name, "w")
  write.table(example, gz1, row.names=F, col.names=F, sep="\t", quote = F)
  close(gz1)
  
  # Return file name
  f_name
  
}

diff_mean_DRB.file <- WriteGRToTable(diff_mean_DRB, bins_example, HMM.dir, "diff_mean_DRB")
diff_mean_Top1.file <- WriteGRToTable(diff_mean_Top1, bins_example, HMM.dir, "diff_mean_Top1")
diff_smooth_mean_DRB.file <- WriteGRToTable(diff_smooth_mean_DRB, bins_example, HMM.dir, "diff_smooth_mean_DRB")
diff_smooth_mean_Top1.file <- WriteGRToTable(diff_smooth_mean_Top1, bins_example, HMM.dir, "diff_smooth_mean_Top1")

# Run HMM script
RunHMM <- function(f_name, HMM.dir) {
  # Prepare command
  command <- paste("Rscript",
                   "bin/HMM_calling/HMM.R",
                   "-n", f_name,
                   "-o", HMM.dir)
  print(command)
  system(command)
  
}

#RunHMM(diff_mean_DRB.file, HMM.dir)
#RunHMM(diff_mean_Top1.file, HMM.dir)
#RunHMM(diff_smooth_mean_DRB.file, HMM.dir)
#RunHMM(diff_smooth_mean_Top1.file, HMM.dir)




# Run script Elzo
sign_regions.dir <- file.path(output.dir, "significant_regions")
dir.create(sign_regions.dir, showWarnings = F)

domainogram.dir <- file.path(output.dir, "domainograms")
dir.create(domainogram.dir, showWarnings = F)

RunDomainoGram <- function(f_name, name) {
  # Prepare command
  command <- paste("Rscript",
                   "ts200406_domain_calling_DRB_Top1/diff_damid_domainogram.R",
                   f_name,
                   "FALSE",
                   name)
  
  print(command)
  system(command)
  
}

#RunDomainoGram(diff_mean_DRB.file, "diff_mean_DRB")


RunSignRegions <- function(f_name, sign_regions.dir, name,
                           invert = "FALSE") {
  
  # Prepare command
  command <- paste("Rscript",
                   "ts200406_domain_calling_DRB_Top1/significant_regions.R",
                   f_name,
                   invert,
                   file.path(sign_regions.dir,
                             paste0(name, "_sign_regions.bed")))
  
  print(command)
  system(command)
  
}

RunSignRegions(diff_mean_DRB.file, sign_regions.dir, "diff_mean_DRB")
RunSignRegions(diff_smooth_mean_DRB.file, sign_regions.dir, "diff_smooth_mean_DRB")
