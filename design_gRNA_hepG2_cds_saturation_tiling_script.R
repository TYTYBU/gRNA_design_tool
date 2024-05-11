setwd("~/Library/CloudStorage/OneDrive-Personal/Tian_Workfile/BWH/BE_design")

library(openxlsx)
source("./gRNA_design_tool/hepG2_design_functions_v2.R")


## for 1/2 tiled screening
# set parameters
relative_target_pos_range = 6
tiling_step = 2
cds_start_flank = 6
cds_end_flank = 0

# import gene list
df.raw = read.xlsx("1022_LDLcodingvariantBE_genes.xlsx", sheet = 1)
ind = which(df.raw$`Screening.priority.tier.(1=highest)` == 1)
gene_list = unique(df.raw$Gene[ind])

# load reference hg19 genome and Hepg2 phase info
hepG2 = readRDS("HepG2_phase_hg19.rds")

# get cds and exon coordinates
df.mane_hg19 = get_mane_hg19(gene_list)
df.mane_hg19_cds = get_mane_hg19_cds(df.mane_hg19, detect_range = (50+max(cds_start_flank, cds_end_flank)))
df.mane_hg19_exon = get_mane_hg19_exon(df.mane_hg19)

# get hepG2 phase info for cds plus detect range, and exon splice sites plus detect range
ls.phase = get_hepG2_phase_cds_tiling(df.mane_hg19_cds, df.mane_hg19_exon, hepG2,
                                      detect_range = (50+max(cds_start_flank, cds_end_flank)))

# get cds tiling targets
df.main_all = get_cds_tiling_pos(df.mane_hg19_cds, 
                                 tiling_step = tiling_step, 
                                 cds_start_flank = cds_start_flank, 
                                 cds_end_flank = cds_end_flank)

# get gRNAs for cds tiling
df.out_all = get_hegG2_gRNA_reporter_cds_tiling_wrapper(df.mane_hg19_cds, df.main_all, ls.phase, 
                                                        relative_target_pos_range=relative_target_pos_range,
                                                        detect_range = 50)
# parse output
df.out2 = cbind(df.main_all[,1:8], df.out_all[,-(1:3)])
df.tile_out = parse_tiling_out(df.out2, col_to_keep = 1:9, gRNA_prefix = "gTile", homozygous_only = T)
write.xlsx(df.tile_out, file = "cds_saturation_tiling_gRNA_reporter_120722.xlsx", overwrite = TRUE)


# wb <- createWorkbook()
# addWorksheet(wb, "cds_tiling")
# writeData(wb, "cds_tiling", df.tile_out)
# 
# # get splicing sites targets
# df.main_all = get_splice_pos_v2(df.mane_hg19_exon)
# 
# # get gRNAs for splicing sites
# df.out_all = get_hegG2_gRNA_reporter_splice_sites_wrapper(df.mane_hg19_exon, df.main_all, ls.phase, relative_target_pos_range=4, detect_range=45)
# 
# # parse output
# ind = match(df.out_all$target_id, table = df.main_all$target_id)
# df.out2 = cbind(df.main_all[ind,], df.out_all[,-(1:3)])
# df.ctrl_out = parse_tiling_out(df.out2, col_to_keep = 1:7, gRNA_prefix = "gCtrl", homozygous_only = T)
# 
# addWorksheet(wb, "splicing_ctrl")
# writeData(wb, "splicing_ctrl", df.ctrl_out)
# 
# saveWorkbook(wb, file = "cds_saturation_tiling_gRNA_reporter.xlsx", overwrite = TRUE)
