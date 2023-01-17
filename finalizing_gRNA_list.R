setwd("~/Library/CloudStorage/OneDrive-Personal/Tian_Workfile/BWH/BE_design")
library(openxlsx)
library(tidyverse)
library(ggplot2)

# import gene list
df.raw = read.xlsx("1022_LDLcodingvariantBE_genes.xlsx", sheet = 1)
ind = which(df.raw$`Screening.priority.tier.(1=highest)` == 1)
gene_list = unique(df.raw$Gene[ind])

gene_list_1 = gene_list[1:8]
gene_list_2 = gene_list[c(9:13,15)]

cols = c("unique_id", "target_id", "gene", "chr", "gene_strand", "gRNA_strand" , 
         "relative_target_pos", "target_base", "gRNA_seq", "reporter_seq")

df1 = read.xlsx("cds_saturation_tiling_gRNA_reporter_120722.xlsx")
df1$gene = stringr::str_split(df1$target_id, "_", simplify = T)[,1]
df1 = df1[df1$gene %in% gene_list_1, cols]

df2 = read.xlsx("UKB_variant_codon_gRNA_reporter_121622.xlsx", sheet = 1)
df2$gene = stringr::str_split(df2$target_id, "_", simplify = T)[,1]
df2 = df2[df2$gene %in% gene_list_2, cols]

df3 = read.xlsx("UKB_variant_codon_gRNA_reporter_121622.xlsx", sheet = 2)
df3$gene = stringr::str_split(df3$target_id, "_", simplify = T)[,1]
df3 = df3[df3$gene %in% gene_list_2, cols]
df3$unique_id = gsub(pattern = "gCtrl", replacement = "gPosCtrl", df3$unique_id)

df4 = read.xlsx("110422_18locus_BEgRNAs.xlsx", sheet = "ABE_negcontrols")
df4[,cols] = NA
df4$unique_id = df4$Full.gRNA_Rep.name
df4$gRNA_seq = df4$gRNA
df4$reporter_seq = df4$`Reporter.32-nt`
df4$relative_target_pos = df4$Target.base.position.in.gRNA
df4 = df4[,cols]

df.all = rbind(df1, df2, df3, df4)

# remove all gRNAs that contain TTTT
ind = grep(pattern = "TTTT", x = df.all$gRNA_seq)
df.all = df.all[-ind,]

# Change gRNA_seq to be called 20nt_gRNA_seq
# Create a new column beside it called 19-20nt_gRNA_seq
# Whenever the initial base in the 20nt_gRNA_seq is a G, remove it to create a 19-nt gRNA.
df.all$`20nt_gRNA_seq` = df.all$gRNA_seq
df.all$`19-20nt_gRNA_seq` = sub(pattern = "^G", replacement = "", x = df.all$gRNA_seq)

# concatenate full oilgo sequence
# order: U6 stub, gRNA_seq (19-20 nt), hairpin_term, reporter_seq, 4nt barcode, r2seq RC
U6_stub = toupper("tggaaaggacgaaacaccg")
hairpin_term = "GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTT"
r2seq_RC = "AGATCGGAAGAGCACACG"
temp = read.xlsx("110422_18locus_BEgRNAs.xlsx", sheet = "1122_18LDLlocus_ABE_gRNAlib")
barcodes = temp$`4nt.barcode`

df.out = tibble(unique_id = df.all$unique_id,
                "20nt_gRNA_seq" = df.all$`20nt_gRNA_seq`,
                U6_stub, 
                "19-20nt_gRNA_seq" = df.all$`19-20nt_gRNA_seq`,
                hairpin_term,
                reporter_seq = df.all$reporter_seq,
                "4nt_barcode" = NA,
                r2seq_RC)

barcodes = c(rep(barcodes, floor(nrow(df.out)/length(barcodes))), barcodes[1:(nrow(df.out) %% length(barcodes))])
df.out$`4nt_barcode` = barcodes
df.out$Full_oligo = paste0(df.out$U6_stub, df.out$`19-20nt_gRNA_seq`, df.out$hairpin_term, df.out$reporter_seq,
                           df.out$`4nt_barcode`, df.out$r2seq_RC)


# check for duplicated gRNAs
temp = table(df.out$`19-20nt_gRNA_seq`)
ind = which(as.numeric(temp)>1)
duplicated_seqs = names(temp)[ind]
temp = c()
for (gRNA_seq in duplicated_seqs){
  ind = which(df.out$`19-20nt_gRNA_seq` == gRNA_seq)
  temp2 = data.frame(gRNA_seq, invloved_gRNA_ids = paste0(df.out$unique_id[ind], collapse = ", "))
  temp = rbind(temp, temp2)
}
write.xlsx(temp, "dulplicated_gRNAs.xlsx")

unique_gRNA_seqs = unique(df.out$`19-20nt_gRNA_seq`)
ind = match(unique_gRNA_seqs, table = df.out$`19-20nt_gRNA_seq`)
df.out2 = df.out[ind,]
df.out2$duplication = NA
ind = match(duplicated_seqs, table = df.out2$`19-20nt_gRNA_seq`)
df.out2$duplication[ind] = temp$invloved_gRNA_ids
write.xlsx(df.out2, file = "final_gRNA_oligos_122222.xlsx")






## draw gRNA count graph
# import gene list
df.raw = read.xlsx("1022_LDLcodingvariantBE_genes.xlsx", sheet = 1)
ind = which(df.raw$`Screening.priority.tier.(1=highest)` == 1)
df.raw = df.raw[ind,c("Gene", "Saturation.tiling.gRNAs", "Variant.focused.gRNAs", "Splice.control.gRNAs")]

df = c()
for (total_genes in 10:15){
  for (tiling_genes in 1:total_genes){
    total_gRNAs = sum(df.raw$Saturation.tiling.gRNAs[1:tiling_genes])
    variant_genes = total_genes - tiling_genes
    if (variant_genes>0){
      total_gRNAs = total_gRNAs + sum(df.raw$Variant.focused.gRNAs[(tiling_genes+1):total_genes]) + sum(df.raw$Splice.control.gRNAs[(tiling_genes+1):total_genes])
    }
    
    temp = data.frame(total_genes, tiling_genes, total_gRNAs)
    df = rbind(df, temp)
  }
}

df$total_genes = as.character(df$total_genes)
gRNA_ct_thres = 13000
x_axis_labels = c(df.raw$Gene[1], paste0("+", df.raw$Gene[2:nrow(df.raw)]))
plt = df %>% ggplot(aes(x=tiling_genes, y=total_gRNAs, group=total_genes, color=total_genes)) + geom_line() +
  geom_point(shape=21, size=1) + geom_hline(yintercept=gRNA_ct_thres, linetype="dashed", color = "red") + 
  scale_x_continuous(breaks = 1:15, labels = x_axis_labels) + 
  scale_y_continuous(breaks = sort(c(seq(5000, 30000, length.out=5), gRNA_ct_thres))) + 
  labs(x = "Tiling genes", y = "gRNA count", color = "Total genes") + theme_light() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
ggsave(filename = "tiling_vs_variant_focused.png", plot = plt, device = "png", width = 8, height = 4)  




