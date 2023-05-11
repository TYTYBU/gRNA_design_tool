setwd("~/Library/CloudStorage/OneDrive-Personal/Tian_Workfile/BWH/BE_design")

library(openxlsx)
library(tidyverse)
source("./gRNA_design_tool/hepG2_design_functions_v2.R")


# required columns for df.main: target_id, chr, hg19_pos
df.raw = read.xlsx("0922_Lipid_Topmed_22gene_08sentinels_prox_finemapping_eQTL.xlsx")
df.main = data.frame(target_id = df.raw$prioritizedVariantID, chr = paste0("chr", df.raw$hg38Chromosome), hg19_pos = df.raw$hg19Position)



gRNA_len = 20
gRNA_ct = 5
detect_range = 50

# subset phase info enough for the input data
hepG2 = readRDS("HepG2_phase_hg19.rds")
ls.phase = get_hepG2_phase(df.main, hepG2, detect_range)

# output gRNAs
df.out = get_hepG2_gRNA_CRISPRi(df.main)


## positive control gRNAs
gRNA_ct = 12

# load gene names
df.raw = read.xlsx("0922_IGVFcholesterol_pilotscreenloci.xlsx")
gene_names = df.raw$id

# construct data frame for pos ctrl gRNAs
df.pos_ctrl = read_tsv("CRISPRi_pos_ctrl-sgrna-designs.txt")
df.pos_ctrl = get_hepG2_pos_ctrl_gRNA_CRISPRi(df.pos_ctrl)


## negative control gRNAs
df.neg_ctrl = read.xlsx("110422_18locus_BEgRNAs.xlsx", sheet = "ABE_negcontrols")



wb <- createWorkbook()
addWorksheet(wb, "CRISPRi_gRNAs")
writeData(wb, "CRISPRi_gRNAs", df.out)
addWorksheet(wb, "pos_control")
writeData(wb, "pos_control", df.pos_ctrl)
addWorksheet(wb, "neg_control")
writeData(wb, "neg_control", df.neg_ctrl)

saveWorkbook(wb, file = "CRISPRi_gRNA_controls_011623.xlsx", overwrite = TRUE)
















