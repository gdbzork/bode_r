#!/usr/bin/Rscript

# Load up all the tRNA data, generate summaries

# Load some needed data and functions
source("~/r/lib/trna_types.R")
source("~/r/lib/countlib.R")

# Load the raw counts, generate some lists
filenames <- commandArgs(TRUE)
raw.data <- loadAll(filenames)
species <- levels(raw.data$species)
tissue <- levels(raw.data$tissue)
filename_tags <- filename2annot(filenames)
library_tags <- annot2tags(filename_tags)

# Create by-species and by-tissue-and-species tables, and normalize them
for (spec in species) {
  for (tiss in tissue) {
    tmp_tbl = raw2counts(raw.data,library_tags,species=spec,tissue=tiss)
    if (sum(c(0,0) == dim(tmp_tbl)) == 0) {
      tag = sprintf("%s_%s",spec,tiss)
      tagn = sprintf("%s_%s_norm",spec,tiss)
      assign(tag,tmp_tbl)
      assign(tagn,makeNormal(tmp_tbl))
    }
  }
  tmp_tbl = raw2counts(raw.data,library_tags,species=spec)
  if (sum(c(0,0) == dim(tmp_tbl)) == 0) {
    tag = sprintf("%s_all",spec)
    tagn = sprintf("%s_all_norm",spec)
    assign(tag,tmp_tbl)
    assign(tagn,makeNormal(tmp_tbl))
  }
}

isotype_counts <- tabulateTypes(species,class='isotype')
anticodon_counts <- tabulateTypes(species,class='anticodon')

anticodon_levels <- mergeSets(list(cfa_liver_norm,hsa_liver_norm,mdo_liver_norm,mml_liver_norm,mmu_liver_norm,rno_liver_norm),c('cfa','hsa','mdo','mml','mmu','rno'),anticodons.core,fld='anticodon')
isotype_levels <- anticodons2isotypes(anticodon_levels)

isotype_norm <- makeNormal(isotype_levels)
anticodon_norm <- makeNormal(anticodon_levels)
