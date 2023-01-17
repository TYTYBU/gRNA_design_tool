setwd("C:/Users/ElninoYu/OneDrive/Tian_Workfile/BWH/BE_design")
setwd("~/Library/CloudStorage/OneDrive-Personal/Tian_Workfile/BWH/BE_design")

library(openxlsx)
library(tidyverse)
source("hepG2_design_functions_v2.R")



# load reference hg19 genome and Hepg2 phase info
hg19 = readRDS("hg19_string.rds")
hepG2 = readRDS("HepG2_phase_hg19.rds")
base_complement = list("A" = "T", "T" = "A", "C" = "G", "G" = "C")

# # check all phase variety
# phase_tb = c()
# for (chr in names(hepG2)){
#   phase_tb = c(phase_tb, hepG2[[chr]]$Phase)
# }
# phase_tb = table(phase_tb)
# barplot(phase_tb)




## for GWAS variants
# prepare input data frame
# set parameters
relative_target_pos_range = 4:8
gRNA_len = 20
reporter_extend_len = 6
detect_range = 50

# required columns for df.main: target_id, chr, hg19_pos
df.raw = read.xlsx("0922_Lipid_Topmed_22gene_08sentinels_prox_finemapping_eQTL.xlsx")
df.main = data.frame(target_id = df.raw$prioritizedVariantID, chr = paste0("chr", df.raw$hg38Chromosome), hg19_pos = df.raw$hg19Position)

# subset phase info enough for the input data
ls.phase = get_hepG2_phase(df.main, hepG2, detect_range)
# output gRNAs and reporters
df.out = get_hegG2_gRNA_reporter(df.main, hg19, detect_range, relative_target_pos_range, gRNA_len, reporter_extend_len)

# temp = get_hegG2_gRNA_reporter_v2(df.main, hg19, detect_range, relative_target_pos_range, gRNA_len, reporter_extend_len)
# df.out = temp[[1]]
# ls.detect = temp[[2]]

# parse output
df.variant_out = parse_out(df.out, gRNA_prefix = "gVariant")



## for GWAS positive control genes
df.raw = read.xlsx("0922_IGVFcholesterol_pilotscreenloci.xlsx")
gene_list = unique(df.raw$id)

# load mane v1.0 track from UCSC genome browser
df.mane = read.table("mane_track_UCSC_v1.0.txt", sep = "\t", header = F) 

# get splice site coordinates and populate df.main
df.main = get_splice_pos(gene_list, df.mane)
# subset phase info enough for the input data
ls.phase = get_hepG2_phase(df.main, hepG2, detect_range)
# output gRNAs and reporters
df.out = get_hegG2_gRNA_reporter(df.main, hg19, detect_range, relative_target_pos_range, gRNA_len, reporter_extend_len)
# parse output
df.control_out = parse_out(df.out, gRNA_prefix = "gControl")

write.xlsx(rbind(df.variant_out, df.control_out), "GWAS_gRNA_reporter.xlsx")











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
# df.out2 = cbind(df.main_all, df.out_all[,-(1:3)])
# df.ctrl_out = parse_tiling_out(df.out2, col_to_keep = 1:7, gRNA_prefix = "gCtrl", homozygous_only = T)
# 
# addWorksheet(wb, "splicing_ctrl")
# writeData(wb, "splicing_ctrl", df.ctrl_out)
# 
# saveWorkbook(wb, file = "cds_saturation_tiling_gRNA_reporter.xlsx", overwrite = TRUE)











  
  
  














# draw a plot showing gRNA and reporter position relative to a target position


relative_target_pos_range = 4:8
gRNA_len = 20
reporter_extend_len = 6
detect_range = 50

strands = c("+", "-") 


blk_height = 1
blk_level = 1
blk_margin = 1

x_max = gRNA_len-min(relative_target_pos_range)+reporter_extend_len
y_max = length(relative_target_pos_range)*(blk_height*2+blk_margin) + blk_height
plot(NA, xlim=c(-x_max,x_max), ylim=c(-y_max,y_max))

# draw genome
if ("+" %in% strands){
  rect(xleft = -detect_range:(detect_range-1), ybottom = 0, xright = -(detect_range-1):detect_range, ytop = blk_height, col = add.alpha("grey", 0.5))
  rect(xleft = 0, ybottom = 0, xright = 1, ytop = blk_height, col = add.alpha("red", 1)) # draw target site
}

if ("-" %in% strands){
  rect(xleft = -detect_range:(detect_range-1), ybottom = 0, xright = -(detect_range-1):detect_range, ytop = -blk_height, col = add.alpha("grey", 0.5))
  rect(xleft = 0, ybottom = 0, xright = 1, ytop = -blk_height, col = add.alpha("red", 1)) # draw target site
}


# draw gRNA and reporter pair
for (i in 1:length(relative_target_pos_range)){
  blk_level = i
  
  if ("+" %in% strands){
    relative_target_pos = relative_target_pos_range[i]
    rect_left_pos = -(relative_target_pos-1)
    rect_left_range = rect_left_pos:(rect_left_pos+gRNA_len-1)
    rect_bottom_pos = (blk_height*2+blk_margin)*blk_level
    rect(xleft = rect_left_range, ybottom = rect_bottom_pos, xright = rect_left_range+1, ytop = rect_bottom_pos+blk_height, col = add.alpha("blue", 0.5))
    rect(xleft = 0, ybottom = rect_bottom_pos, xright = 1, ytop = rect_bottom_pos+blk_height, col = add.alpha("red", 1)) # draw target site
    
    rect_left_pos = -(relative_target_pos-1+reporter_extend_len)
    rect_left_range = rect_left_pos:(rect_left_pos+gRNA_len-1+reporter_extend_len*2)
    rect_bottom_pos = (blk_height*2+blk_margin)*blk_level-blk_height
    rect(xleft = rect_left_range, ybottom = rect_bottom_pos, xright = rect_left_range+1, ytop = rect_bottom_pos+blk_height, col = add.alpha("orange", 0.5))
    rect(xleft = 0, ybottom = rect_bottom_pos, xright = 1, ytop = rect_bottom_pos+blk_height, col = add.alpha("red", 1)) # draw target site
  }
  
  if ("-" %in% strands){
    relative_target_pos = relative_target_pos_range[i]
    rect_left_pos = -(gRNA_len-relative_target_pos)
    rect_left_range = rect_left_pos:(rect_left_pos+gRNA_len-1)
    rect_bottom_pos = (blk_height*2+blk_margin)*blk_level
    rect(xleft = rect_left_range, ybottom = -rect_bottom_pos, xright = rect_left_range+1, ytop = -(rect_bottom_pos+blk_height), col = add.alpha("blue", 0.5))
    rect(xleft = 0, ybottom = -rect_bottom_pos, xright = 1, ytop = -(rect_bottom_pos+blk_height), col = add.alpha("red", 1)) # draw target site
    
    rect_left_pos = -(gRNA_len-relative_target_pos+reporter_extend_len)
    rect_left_range = rect_left_pos:(rect_left_pos+gRNA_len-1+reporter_extend_len*2)
    rect_bottom_pos = (blk_height*2+blk_margin)*blk_level-blk_height
    rect(xleft = rect_left_range, ybottom = -rect_bottom_pos, xright = rect_left_range+1, ytop = -(rect_bottom_pos+blk_height), col = add.alpha("orange", 0.5))
    rect(xleft = 0, ybottom = -rect_bottom_pos, xright = 1, ytop = -(rect_bottom_pos+blk_height), col = add.alpha("red", 1)) # draw target site
  }
  

}




text(x = 1:10, y = 1, labels = c("X", "A", "B"))









polygon(c(0.3, 0.4, 0.5), c(0.05, 0.4, 0.05), lwd=4)
polygon (c(0, 0.1, 0.5, 0.5, 0.1), c(0.4, 0.5, 0.5, 0.3, 0.3), col="green")

