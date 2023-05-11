setwd("~/Library/CloudStorage/OneDrive-Personal/Tian_Workfile/BWH/BE_design")
library(openxlsx)
source("./gRNA_design_tool/hepG2_design_functions_v2.R")

# import gene list
df.raw = read.xlsx("1022_LDLcodingvariantBE_genes.xlsx", sheet = 1)
ind = which(df.raw$`Screening.priority.tier.(1=highest)` == 1)
gene_list = unique(df.raw$Gene[ind])

# load reference hg19 genome and Hepg2 phase info
hepG2 = readRDS("HepG2_phase_hg19.rds")

# get cds and exon coordinates
df.mane_hg19 = get_mane_hg19(gene_list)
df.mane_hg19_cds = get_mane_hg19_cds(df.mane_hg19)
df.mane_hg19_exon = get_mane_hg19_exon(df.mane_hg19)

# get hepG2 phase info for cds plus detect range, and exon splice sites plus detect range
ls.phase = get_hepG2_phase_cds_tiling(df.mane_hg19_cds, df.mane_hg19_exon, hepG2)

# extract variants from vcf
library(VariantAnnotation)
df.variants = c()
for (gene in gene_list){
  vcf = readVcf(paste0("./vcf/", gene, "_variants.vcf.gz"))
  variants = names(vcf)
  temp = stringr::str_split(variants, pattern = "_", simplify = T)
  colnames(temp) = c("chr", "hg38_pos", "ref", "alt")
  temp = data.frame(gene, temp)
  df.variants = rbind(df.variants, temp)
}

#filter for missense variant
ind = which(nchar(df.variants$ref) == 1 & nchar(df.variants$alt) == 1)
df.variants = df.variants[ind,]

# liftover from hg38 coordinates to hg19
pos.hg38 = GRanges(seqnames = df.variants$chr, ranges = as.character(df.variants$hg38_pos))
ch = import.chain(system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain"))
pos.hg19 = unlist(liftOver(pos.hg38, ch))
df.variants$hg19_pos = start(pos.hg19)

# get gRNAs for UKB variant codons
df.out = get_hegG2_gRNA_reporter_cds_UKBB_variant(df.mane_hg19_cds, df.variants, ls.phase)

# parse output
df.UKB_variant_out = parse_UKB_variant_out(df.out, col_to_keep = 1:4, gRNA_prefix = "gUKB", homozygous_only = T)

wb <- createWorkbook()
addWorksheet(wb, "UKB_variant_codon")
writeData(wb, "UKB_variant_codon", df.UKB_variant_out)
# write.xlsx(df.UKB_variant_out, file = "UKB_variant_codon_gRNA_reporter_121622.xlsx", overwrite = TRUE)


# get splicing sites targets
df.main_all = get_splice_pos_v2(df.mane_hg19_exon)

# only target acceptor A on sense strand, and donor T on antisense strand
ind = grep(pattern = "_acceptorA_sense|_donorT_antisense", df.main_all$target_id)
df.main_all = df.main_all[ind,]

# get gRNAs for splicing sites
df.out_all = get_hegG2_gRNA_reporter_splice_sites_wrapper(df.mane_hg19_exon, df.main_all, ls.phase, relative_target_pos_range=5:7, detect_range=45)

# parse output
ind = match(df.out_all$target_id, table = df.main_all$target_id)
df.out2 = cbind(df.main_all[ind,], df.out_all[,-(1:3)])
df.ctrl_out = parse_splice_site_out(df.out2, col_to_keep = 1:7, gRNA_prefix = "gCtrl", homozygous_only = T)

addWorksheet(wb, "splice_site_ctrl")
writeData(wb, "splice_site_ctrl", df.ctrl_out)
saveWorkbook(wb, file = "UKB_variant_codon_gRNA_reporter_121622.xlsx", overwrite = TRUE)




