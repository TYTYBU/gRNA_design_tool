get_mane_v2 <- function(gene_list, hg19=F, current_url = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.2.ensembl_genomic.gtf.gz"){
  library(rtracklayer)
  library(liftOver)
  
  # import mane gtf
  df.mane_hg38 = import(current_url)
  
  if (!all(gene_list %in% df.mane_hg38$gene_name)){
    ind = which(!(gene_list %in% df.mane_hg38$gene_name))
    simpleError(message = paste0("Some genes are not found in MANE: ", paste(gene_list[ind], collapse = ",")))
  }
  df.mane_hg38 = df.mane_hg38[(df.mane_hg38$gene_name %in% gene_list),]
  
  if (hg19){
    # liftover from hg38 coordinates to hg19 coordinates
    ch = import.chain(system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain"))
    df.mane_hg19 = unlist(liftOver(df.mane_hg38, ch))
    
    return(df.mane_hg19)
  } else {
    return(df.mane_hg38)
  }
}

get_mane_hg19 <- function(gene_list, current_url = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.2.ensembl_genomic.gtf.gz"){
  library(rtracklayer)
  library(liftOver)
  
  # import mane gtf
  df.mane_hg38 = import(current_url)
  
  if (!all(gene_list %in% df.mane_hg38$gene_name)){
    ind = which(!(gene_list %in% df.mane_hg38$gene_name))
    simpleError(message = paste0("Some genes are not found in MANE: ", paste(gene_list[ind], collapse = ",")))
  }
  df.mane_hg38 = df.mane_hg38[(df.mane_hg38$gene_name %in% gene_list),]
  
  # liftover from hg38 coordinates to hg19 coordinates
  ch = import.chain(system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain"))
  df.mane_hg19 = unlist(liftOver(df.mane_hg38, ch))
  
  return(df.mane_hg19)
}

get_mane_hg19_cds <- function(df.mane_hg19, detect_range=50){
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  df.mane_hg19_cds = df.mane_hg19[df.mane_hg19$type == "CDS",]
  
  # import hg19 genome
  genome <- BSgenome.Hsapiens.UCSC.hg19
  
  # get cds sequences on + strand, plus detect_range
  temp_cds = df.mane_hg19_cds
  strand(temp_cds) = "+"
  start(temp_cds) = start(temp_cds) - detect_range
  end(temp_cds) = end(temp_cds) + detect_range
  df.mane_hg19_cds$cds_extend_start = start(temp_cds)
  df.mane_hg19_cds$cds_extend_end = end(temp_cds)
  df.mane_hg19_cds$cds_extend_seq = as.character(getSeq(genome, temp_cds))
  
  return(df.mane_hg19_cds)
}

get_mane_hg19_exon <- function(df.mane_hg19, detect_range=50){
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  df.mane_hg19_exon = df.mane_hg19[df.mane_hg19$type == "exon",]
  
  # import hg19 genome
  genome <- BSgenome.Hsapiens.UCSC.hg19
  
  # get exon_start, exon_end sequences, plus detect_range
  temp_exon = df.mane_hg19_exon
  strand(temp_exon) = "+"
  exon_start = start(temp_exon)
  start(temp_exon) = exon_start - detect_range
  end(temp_exon) = exon_start + detect_range
  df.mane_hg19_exon$exon_start_extend_start = start(temp_exon)
  df.mane_hg19_exon$exon_start_extend_end = end(temp_exon)
  df.mane_hg19_exon$exon_start_extend_seq = as.character(getSeq(genome, temp_exon))
  
  temp_exon = df.mane_hg19_exon
  strand(temp_exon) = "+"
  exon_end = end(temp_exon)
  start(temp_exon) = exon_end - detect_range
  end(temp_exon) = exon_end + detect_range
  df.mane_hg19_exon$exon_end_extend_start = start(temp_exon)
  df.mane_hg19_exon$exon_end_extend_end = end(temp_exon)
  df.mane_hg19_exon$exon_end_extend_seq = as.character(getSeq(genome, temp_exon))
  
  # get max exon number for each gene
  df.mane_hg19_exon$total_exon_number = 1
  for (gene_name in unique(df.mane_hg19_exon$gene_name)){
    ind = which(df.mane_hg19_exon$gene_name == gene_name)
    df.mane_hg19_exon$total_exon_number[ind] = max(df.mane_hg19_exon$exon_number[ind])
  }
  
  return(df.mane_hg19_exon)
}

get_hepG2_phase <- function(df.main, hepG2, detect_range){
  # subset phase info enough for the input data
  ls.pos_range = list()
  pb = txtProgressBar(min = 0, max = nrow(df.main), initial = 0, style = 3)
  
  for (i in 1:nrow(df.main)){
    chr = df.main$chr[i]
    pos = df.main$hg19_pos[i]
    ls.pos_range[[chr]] = c(ls.pos_range[[chr]], (pos-detect_range):(pos+detect_range))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  ls.phase = list()
  for (chr in names(ls.pos_range)){
    ls.pos_range[[chr]] = unique(ls.pos_range[[chr]])
    ind = which(hepG2[[chr]]$Pos %in% ls.pos_range[[chr]])
    ls.phase[[chr]] = hepG2[[chr]][ind,]
  }
  
  return(ls.phase)
}

get_hepG2_phase_cds_tiling <- function(df.mane_hg19_cds, df.mane_hg19_exon=NULL, hepG2, detect_range=50){
  df.cds_meta = mcols(df.mane_hg19_cds)
  ls.pos_range = list()
  for (i in 1:length(df.mane_hg19_cds)){
    chr = as.character(seqnames(df.mane_hg19_cds)[i])
    cds_start = start(df.mane_hg19_cds)[i]
    cds_end = end(df.mane_hg19_cds)[i]
    ls.pos_range[[chr]] = c(ls.pos_range[[chr]], (cds_start-detect_range):(cds_end+detect_range))
  }
  
  if (!is.null(df.mane_hg19_exon)){
    for (i in 1:length(df.mane_hg19_exon)){
      chr = as.character(seqnames(df.mane_hg19_exon)[i])
      exon_start = start(df.mane_hg19_exon)[i]
      exon_end = end(df.mane_hg19_exon)[i]
      ls.pos_range[[chr]] = c(ls.pos_range[[chr]], (exon_start-detect_range):(exon_start+detect_range), (exon_end-detect_range):(exon_end+detect_range))
    }
  }
  
  ls.phase = list()
  for (chr in names(ls.pos_range)){
    ls.pos_range[[chr]] = unique(ls.pos_range[[chr]])
    ind = which(hepG2[[chr]]$Pos %in% ls.pos_range[[chr]])
    ls.phase[[chr]] = hepG2[[chr]][ind,]
  }
  
  return(ls.phase)
}

get_hegG2_gRNA_reporter <- function(df.main, hg19, detect_range, relative_target_pos_range, gRNA_len, reporter_extend_len){
  base_complement = list("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  
  # extract reference seq around target site
  df.main$ref_seq = NA
  pb = txtProgressBar(min = 0, max = nrow(df.main), initial = 0, style = 3)
  
  for (i in 1:nrow(df.main)){
    chr = df.main$chr[i]
    pos = df.main$hg19_pos[i]
    df.main$n_phase[i] = sum((pos-detect_range):(pos+detect_range) %in% ls.phase[[chr]]$Pos)
    df.main$ref_seq[i] = substr(hg19[[chr]], start = pos-detect_range, stop = pos+detect_range)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # output gRNAs and reporters
  df.out = c()
  pb = txtProgressBar(min = 0, max = nrow(df.main), initial = 0, style = 3)
  
  for (i in 1:nrow(df.main)){
    # derive haplotype 1 and haplotype 2
    chr = df.main$chr[i]
    pos = df.main$hg19_pos[i]
    df.detect = data.frame(hg19_pos = (pos-detect_range):(pos+detect_range), 
                           hg19_base = unlist(str_split(df.main$ref_seq[i], "")))
    
    ind = match(df.detect$hg19_pos, table = ls.phase[[chr]]$Pos)
    temp = ls.phase[[chr]][ind, 3:5]
    colnames(temp) = paste0("HepG2_", colnames(temp))
    df.detect = cbind(df.detect, temp)
    
    df.detect$h1 = df.detect$hg19_base
    df.detect$h2 = df.detect$hg19_base
    for (j in grep(pattern = "\\d+\\|\\d+", x = df.detect$HepG2_Phase)){ # consider only high-quality phase info, e.g. 0|1
      # when HepG2_Ref is more than 1 base, remove bases from following positions to compensate
      hepG2_ref_len = nchar(df.detect$HepG2_Ref[j])
      if (hepG2_ref_len > 1 & j < nrow(df.detect)){
        add_max = hepG2_ref_len-1
        if (add_max > nrow(df.detect)-j){
          add_max = nrow(df.detect)-j
        }
        df.detect$h1[j+1:add_max] = ""
        df.detect$h2[j+1:add_max] = ""
      }
      
      phase = as.numeric(unlist(str_split(df.detect$HepG2_Phase[j], pattern = "\\|")))
      hepG2_genotypes = c(df.detect$HepG2_Ref[j], unlist(str_split(df.detect$HepG2_Alt[j], pattern = ",")))
      df.detect$h1[j] = hepG2_genotypes[phase[1]+1]
      df.detect$h2[j] = hepG2_genotypes[phase[2]+1]
    }
    
    # recreate df.detect to ensure haplotype 1 and 2 have single base at each position
    for (h in c("h1", "h2")){
      h_corrected = paste0(h, "_corrected")
      df.detect[[h_corrected]] = NA
      
      temp = unlist(str_split(df.detect[[h]][1:detect_range], pattern = "")) # 50 bases before target position
      if (length(temp) < detect_range){
        temp = c(rep("", detect_range-length(temp)), temp)
      }
      df.detect[[h_corrected]][1:detect_range] = temp[(length(temp)-detect_range+1):length(temp)]
      
      temp = unlist(str_split(df.detect[[h]][(detect_range+1):(detect_range*2+1)], pattern = "")) # 51 bases on and after target position
      if (length(temp) < (detect_range+1)){
        temp = c(temp, rep("", detect_range+1-length(temp)))
      }
      df.detect[[h_corrected]][(detect_range+1):(detect_range*2+1)] = temp[1:(detect_range+1)]
    }
    
    # output gRNAs and reporters
    temp = data.frame(target_id = df.main$target_id[i], chr = df.main$chr[i], target_hg19_pos = df.main$hg19_pos[i], relative_target_pos = relative_target_pos_range)
    for (h in c("h1", "h2")){
      h_corrected = paste0(h, "_corrected")
      h_seq = df.detect[[h_corrected]]
      
      h_target_base = paste0(h, "_target_base")
      temp[[h_target_base]] = h_seq[detect_range+1]
      
      h_gRNA_strand = paste0(h, "_gRNA_strand")
      temp[[h_gRNA_strand]] = "+"
      
      h_gRNA_seq = paste0(h, "_gRNA_seq")
      temp[[h_gRNA_seq]] = NA
      
      h_reporter_seq = paste0(h, "_reporter_seq")
      temp[[h_reporter_seq]] = NA
      
      if ("gRNA_strand" %in% colnames(df.main)){
        if (df.main$gRNA_strand[i] == "+"){
          final_seq = h_seq
        } else {
          final_seq = rev(unlist(base_complement[h_seq]))
          temp[[h_gRNA_strand]] = "-"
        }
      } else {
        if (h_seq[detect_range+1] %in% c("A", "C")){
          final_seq = h_seq
          
        } else {
          final_seq = rev(unlist(base_complement[h_seq]))
          temp[[h_gRNA_strand]] = "-"
        }
      }
      
      for (m in 1:nrow(temp)){
        gRNA_start_ind = (detect_range+1)-(temp$relative_target_pos[m]-1)
        gRNA_seq = paste0(final_seq[gRNA_start_ind:(gRNA_start_ind+gRNA_len-1)], collapse = "")
        temp[[h_gRNA_seq]][m] = gRNA_seq
        reporter_seq = paste0(final_seq[(gRNA_start_ind-reporter_extend_len):(gRNA_start_ind+gRNA_len-1+reporter_extend_len)], collapse = "")
        temp[[h_reporter_seq]][m] = reporter_seq
      }
    }
    
    df.out = rbind(df.out, temp)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  df.out$is_gRNA_same = F
  ind = which(df.out$h1_gRNA_seq == df.out$h2_gRNA_seq)
  df.out$is_gRNA_same[ind] = T
  df.out$is_reporter_same = F
  ind2 = which(df.out$h1_reporter_seq == df.out$h2_reporter_seq)
  df.out$is_reporter_same[ind2] = T
  
  return(df.out)
}

get_hegG2_gRNA_reporter_cds_tiling <- function(df.main, df.detect_cds, relative_target_pos_range=c(4), 
                                               gRNA_len=20, reporter_extend_len=6, detect_range=50){
  base_complement = list("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  
  df.out = c()
  # output gRNAs and reporters
  for (k in 1:nrow(df.main)){
    # derive haplotype 1 and haplotype 2
    pos = df.main$hg19_pos[k]
    ind = which(df.detect_cds$hg19_pos >= (pos-detect_range) & df.detect_cds$hg19_pos <= (pos+detect_range))
    df.detect = df.detect_cds[ind,]
    df.detect$relative_to_target = df.detect$hg19_pos - pos
    
    df.detect$h1 = df.detect$hg19_base
    df.detect$h2 = df.detect$hg19_base
    for (j in grep(pattern = "\\d+\\|\\d+", x = df.detect$HepG2_Phase)){ # consider only high-quality phase info, e.g. 0|1
      # when HepG2_Ref is more than 1 base, remove bases from following positions to compensate
      hepG2_ref_len = nchar(df.detect$HepG2_Ref[j])
      if (hepG2_ref_len > 1 & j < nrow(df.detect)){
        add_max = hepG2_ref_len-1
        if (add_max > nrow(df.detect)-j){
          add_max = nrow(df.detect)-j
        }
        df.detect$h1[j+1:add_max] = ""
        df.detect$h2[j+1:add_max] = ""
      }
      
      phase = as.numeric(unlist(stringr::str_split(df.detect$HepG2_Phase[j], pattern = "\\|")))
      hepG2_genotypes = c(df.detect$HepG2_Ref[j], unlist(stringr::str_split(df.detect$HepG2_Alt[j], pattern = ",")))
      df.detect$h1[j] = hepG2_genotypes[phase[1]+1]
      df.detect$h2[j] = hepG2_genotypes[phase[2]+1]
    }
    
    # recreate df.detect to ensure haplotype 1 and 2 have single base at each position
    for (h in c("h1", "h2")){
      h_corrected = paste0(h, "_corrected")
      df.detect[[h_corrected]] = NA
      
      temp = unlist(stringr::str_split(df.detect[[h]][1:detect_range], pattern = "")) # 50 bases before target position
      if (length(temp) < detect_range){
        temp = c(rep("", detect_range-length(temp)), temp)
      }
      df.detect[[h_corrected]][1:detect_range] = temp[(length(temp)-detect_range+1):length(temp)]
      
      temp = unlist(stringr::str_split(df.detect[[h]][(detect_range+1):(detect_range*2+1)], pattern = "")) # 51 bases on and after target position
      if (length(temp) < (detect_range+1)){
        temp = c(temp, rep("", detect_range+1-length(temp)))
      }
      df.detect[[h_corrected]][(detect_range+1):(detect_range*2+1)] = temp[1:(detect_range+1)]
    }
    
    # output gRNAs and reporters
    temp = data.frame(target_id = df.main$target_id[k], chr = df.main$chr[k], target_hg19_pos = df.main$hg19_pos[k], relative_target_pos = relative_target_pos_range)
    for (h in c("h1", "h2")){
      h_corrected = paste0(h, "_corrected")
      h_seq = df.detect[[h_corrected]]
      
      h_target_base = paste0(h, "_target_base")
      temp[[h_target_base]] = h_seq[detect_range+1]
      
      h_gRNA_strand = paste0(h, "_gRNA_strand")
      temp[[h_gRNA_strand]] = "+"
      
      h_gRNA_seq = paste0(h, "_gRNA_seq")
      temp[[h_gRNA_seq]] = NA
      
      h_reporter_seq = paste0(h, "_reporter_seq")
      temp[[h_reporter_seq]] = NA
      
      if ("gRNA_strand" %in% colnames(df.main)){
        if (df.main$gRNA_strand[k] == "+"){
          final_seq = h_seq
        } else {
          final_seq = rev(unlist(base_complement[h_seq]))
          temp[[h_gRNA_strand]] = "-"
        }
      } else {
        if (h_seq[detect_range+1] %in% c("A", "C")){
          final_seq = h_seq
          
        } else {
          final_seq = rev(unlist(base_complement[h_seq]))
          temp[[h_gRNA_strand]] = "-"
        }
      }
      
      for (m in 1:nrow(temp)){
        gRNA_start_ind = (detect_range+1)-(temp$relative_target_pos[m]-1)
        gRNA_seq = paste0(final_seq[gRNA_start_ind:(gRNA_start_ind+gRNA_len-1)], collapse = "")
        temp[[h_gRNA_seq]][m] = gRNA_seq
        reporter_seq = paste0(final_seq[(gRNA_start_ind-reporter_extend_len):(gRNA_start_ind+gRNA_len-1+reporter_extend_len)], collapse = "")
        temp[[h_reporter_seq]][m] = reporter_seq
      }
    }
    
    df.out = rbind(df.out, temp)
  }
  return(df.out)
}

get_hegG2_gRNA_reporter_cds_tiling_wrapper <- function(df.mane_hg19_cds, df.main_all, ls.phase, relative_target_pos_range=c(4), 
                                                       gRNA_len=20, reporter_extend_len=6, detect_range=50){
  # output gRNAs and reporters
  df.out_all = c()
  pb = txtProgressBar(min = 0, max = length(df.mane_hg19_cds), initial = 0, style = 3)
  df.cds_meta = mcols(df.mane_hg19_cds)
  
  for (i in 1:length(df.mane_hg19_cds)){
    # for each cds:
    # construct a df.detect_cds to hold hepG2 phase info
    chr = as.character(seqnames(df.mane_hg19_cds)[i])
    df.detect_cds = data.frame(hg19_pos = df.cds_meta$cds_extend_start[i]:df.cds_meta$cds_extend_end[i], 
                               hg19_base = unlist(stringr::str_split(df.cds_meta$cds_extend_seq[i], "")))
    ind = match(df.detect_cds$hg19_pos, table = ls.phase[[chr]]$Pos)
    temp = ls.phase[[chr]][ind, 3:5]
    colnames(temp) = paste0("HepG2_", colnames(temp))
    df.detect_cds = cbind(df.detect_cds, temp)  
    
    # for each cds:
    # construct a df.main for target_id, target chr/pos/strand/base
    ind = which(df.main_all$mane_hg19_cds_ind == i)
    df.main = df.main_all[ind,]
    
    df.out = get_hegG2_gRNA_reporter_cds_tiling(df.main, df.detect_cds, relative_target_pos_range, gRNA_len, reporter_extend_len, detect_range)
    df.out_all = rbind(df.out_all, df.out)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  
  df.out_all$is_gRNA_same = F
  ind = which(df.out_all$h1_gRNA_seq == df.out_all$h2_gRNA_seq)
  df.out_all$is_gRNA_same[ind] = T
  df.out_all$is_reporter_same = F
  ind2 = which(df.out_all$h1_reporter_seq == df.out_all$h2_reporter_seq)
  df.out_all$is_reporter_same[ind2] = T
  
  return(df.out_all)
}

get_hegG2_gRNA_reporter_cds_UKBB_variant <- function(df.mane_hg19_cds, df.variants, ls.phase, 
                                                     gRNA_len=20, reporter_extend_len=6, detect_range=50){
  base_complement = list("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  
  # filter for coding variants
  df.variants$mane_hg19_cds_ind = NA
  df.mane_hg19_cds$HepG2_ct = 0
  df.mane_hg19_cds$UKBB_ct = 0
  for (i in 1:length(df.mane_hg19_cds)){
    chr = as.character(seqnames(df.mane_hg19_cds)[i])
    cds_start = start(df.mane_hg19_cds)[i]
    cds_end = end(df.mane_hg19_cds)[i]
    ind = which(df.variants$chr == chr & df.variants$hg19_pos >= cds_start & df.variants$hg19_pos <= cds_end)
    if (length(ind)>0){
      df.variants$mane_hg19_cds_ind[ind] = i
    }
    ind = which(ls.phase[[chr]]$Pos >= cds_start & ls.phase[[chr]]$Pos <= cds_end)
    df.mane_hg19_cds$HepG2_ct[i] = length(ind)
    ind = which(df.variants$chr == chr & df.variants$hg19_pos >= cds_start & df.variants$hg19_pos <= cds_end)
    df.mane_hg19_cds$UKBB_ct[i] = length(ind)
  }
  df.variants = df.variants[!is.na(df.variants$mane_hg19_cds_ind),]
  
  ind_UKBB = which(df.mane_hg19_cds$UKBB_ct > 0)
  pb = txtProgressBar(min = 0, max = max(ind_UKBB), initial = 0, style = 3)
  df.out = c()
  for (i in ind_UKBB){
    # for each cds, construct a df.detect_cds
    gene = df.mane_hg19_cds$gene_name[i]
    chr = as.character(seqnames(df.mane_hg19_cds)[i])
    chr_num = sub(pattern = "chr", replacement = "", chr)
    strand = as.character(strand(df.mane_hg19_cds)[i])
    cds_start = start(df.mane_hg19_cds)[i]
    cds_end = end(df.mane_hg19_cds)[i]
    cds_phase = df.mane_hg19_cds$phase[i]
    df.detect_cds = data.frame(hg19_pos = df.mane_hg19_cds$cds_extend_start[i]:df.mane_hg19_cds$cds_extend_end[i], 
                               hg19_base = unlist(stringr::str_split(df.mane_hg19_cds$cds_extend_seq[i], "")))
    df.detect_cds$has_UKBB_variant = F
    ind = which(df.variants$chr == chr & df.variants$hg19_pos >= cds_start & df.variants$hg19_pos <= cds_end)
    ind2 = match(df.variants$hg19_pos[ind], table = df.detect_cds$hg19_pos)
    df.detect_cds$has_UKBB_variant[ind2] = T
    
    df.detect_cds$is_target_pos_strand = F
    df.detect_cds$is_target_neg_strand = F
    
    # attach hepG2 phase info
    ind = match(df.detect_cds$hg19_pos, table = ls.phase[[chr]]$Pos)
    temp = ls.phase[[chr]][ind, 3:5]
    colnames(temp) = paste0("HepG2_", colnames(temp))
    df.detect_cds = cbind(df.detect_cds, temp)
    
    df.detect_cds$h1 = df.detect_cds$hg19_base
    df.detect_cds$h2 = df.detect_cds$hg19_base
    for (j in grep(pattern = "\\d+\\|\\d+", x = df.detect_cds$HepG2_Phase)){ # consider only high-quality phase info, e.g. 0|1
      # indels are not considered in the current analysis
      hepG2_ref_len = nchar(df.detect_cds$HepG2_Ref[j])
      hepG2_alt_len = nchar(df.detect_cds$HepG2_Alt[j])
      if (hepG2_ref_len == 1 & hepG2_alt_len == 1){
        phase = as.numeric(unlist(stringr::str_split(df.detect_cds$HepG2_Phase[j], pattern = "\\|")))
        hepG2_genotypes = c(df.detect_cds$HepG2_Ref[j], unlist(stringr::str_split(df.detect_cds$HepG2_Alt[j], pattern = ",")))
        df.detect_cds$h1[j] = hepG2_genotypes[phase[1]+1]
        df.detect_cds$h2[j] = hepG2_genotypes[phase[2]+1]
      }
    }
    
    # attach codon pos
    df.detect_cds$codon_position = NA
    if (cds_phase == 0){
      codon_position = c(1,2,3)
    } else if (cds_phase == 1){
      codon_position = c(3,1,2)
    } else if (cds_phase == 2){
      codon_position = c(2,3,1)
    }
    
    ind = match(x = cds_start:cds_end, table = df.detect_cds$hg19_pos)
    if (strand == "-"){
      ind = rev(ind)
    } 
    df.detect_cds$codon_position[ind] = codon_position
    
    # get codon on positive strand and negative strand
    # get target position (inital nucleotide in codon) on both strands
    df.detect_cds$h1_pos_strand_codon = NA
    df.detect_cds$h2_pos_strand_codon = NA
    for (n in which(!is.na(df.detect_cds$codon_position))){
      if (strand == "+"){ # when gene on positive strand
        if (df.detect_cds$codon_position[n] == 1){
          df.detect_cds$h1_pos_strand_codon[n] = paste0(df.detect_cds$h1[n:(n+2)], collapse = "")
          df.detect_cds$h2_pos_strand_codon[n] = paste0(df.detect_cds$h2[n:(n+2)], collapse = "")
          if (any(df.detect_cds$has_UKBB_variant[n:(n+2)])){
            df.detect_cds$is_target_pos_strand[n] = T
            df.detect_cds$is_target_neg_strand[n+2] = T
          }
        } else if (df.detect_cds$codon_position[n] == 2){
          df.detect_cds$h1_pos_strand_codon[n] = paste0(df.detect_cds$h1[(n-1):(n+1)], collapse = "")
          df.detect_cds$h2_pos_strand_codon[n] = paste0(df.detect_cds$h2[(n-1):(n+1)], collapse = "")
          if (any(df.detect_cds$has_UKBB_variant[(n-1):(n+1)])){
            df.detect_cds$is_target_pos_strand[n-1] = T
            df.detect_cds$is_target_neg_strand[n+1] = T
          }
        } else if (df.detect_cds$codon_position[n] == 3){
          df.detect_cds$h1_pos_strand_codon[n] = paste0(df.detect_cds$h1[(n-2):n], collapse = "")
          df.detect_cds$h2_pos_strand_codon[n] = paste0(df.detect_cds$h2[(n-2):n], collapse = "")
          if (any(df.detect_cds$has_UKBB_variant[(n-2):n])){
            df.detect_cds$is_target_pos_strand[n-2] = T
            df.detect_cds$is_target_neg_strand[n] = T
          }
        }
      } else { # when gene on negative strand
        if (df.detect_cds$codon_position[n] == 1){
          df.detect_cds$h1_pos_strand_codon[n] = paste0(df.detect_cds$h1[(n-2):n], collapse = "")
          df.detect_cds$h2_pos_strand_codon[n] = paste0(df.detect_cds$h2[(n-2):n], collapse = "")
          if (any(df.detect_cds$has_UKBB_variant[(n-2):n])){
            df.detect_cds$is_target_pos_strand[n-2] = T
            df.detect_cds$is_target_neg_strand[n] = T
          }
        } else if (df.detect_cds$codon_position[n] == 2){
          df.detect_cds$h1_pos_strand_codon[n] = paste0(df.detect_cds$h1[(n-1):(n+1)], collapse = "")
          df.detect_cds$h2_pos_strand_codon[n] = paste0(df.detect_cds$h2[(n-1):(n+1)], collapse = "")
          if (any(df.detect_cds$has_UKBB_variant[(n-1):(n+1)])){
            df.detect_cds$is_target_pos_strand[n-1] = T
            df.detect_cds$is_target_neg_strand[n+1] = T
          }
        } else if (df.detect_cds$codon_position[n] == 3){
          df.detect_cds$h1_pos_strand_codon[n] = paste0(df.detect_cds$h1[n:(n+2)], collapse = "")
          df.detect_cds$h2_pos_strand_codon[n] = paste0(df.detect_cds$h2[n:(n+2)], collapse = "")
          if (any(df.detect_cds$has_UKBB_variant[n:(n+2)])){
            df.detect_cds$is_target_pos_strand[n] = T
            df.detect_cds$is_target_neg_strand[n+2] = T
          }
        }
      }
    }
    
    df.detect_cds$h1_neg_strand_codon = sapply(df.detect_cds$h1_pos_strand_codon, get_revcom_seq, USE.NAMES = F)
    df.detect_cds$h2_neg_strand_codon = sapply(df.detect_cds$h2_pos_strand_codon, get_revcom_seq, USE.NAMES = F)
    
    # determine codon types
    df.detect_cds$h1_pos_strand_codon_type = gsub(pattern = "T|C|G", replacement = "B", df.detect_cds$h1_pos_strand_codon)
    df.detect_cds$h2_pos_strand_codon_type = gsub(pattern = "T|C|G", replacement = "B", df.detect_cds$h2_pos_strand_codon)
    df.detect_cds$h1_neg_strand_codon_type = gsub(pattern = "T|C|G", replacement = "B", df.detect_cds$h1_neg_strand_codon)
    df.detect_cds$h2_neg_strand_codon_type = gsub(pattern = "T|C|G", replacement = "B", df.detect_cds$h2_neg_strand_codon)
    
    # output gRNAs and reporters
    for (gRNA_strand in c("pos", "neg")){
      for (k in which(df.detect_cds[[paste0("is_target_", gRNA_strand, "_strand")]])){
        pos = df.detect_cds$hg19_pos[k]
        ind = which(df.detect_cds$hg19_pos >= (pos-detect_range) & df.detect_cds$hg19_pos <= (pos+detect_range))
        df.detect = df.detect_cds[ind,]
        df.detect$relative_to_target = df.detect$hg19_pos - pos
        target_id = paste(gene, chr_num, pos, sep = "_")
      
        temp = data.frame(target_id, chr, target_hg19_pos = pos, gene_strand = strand)
        for (h in c("h1", "h2")){
          h_codon_type = df.detect_cds[[paste(h, gRNA_strand, "strand_codon_type", sep = "_")]][k]
          h_pos_strand_codon = df.detect_cds[[paste0(h, "_pos_strand_codon")]][k]
          h_pos_strand_codon_type = df.detect_cds[[paste0(h, "_pos_strand_codon_type")]][k]
          h_neg_strand_codon = df.detect_cds[[paste0(h, "_neg_strand_codon")]][k]
          h_neg_strand_codon_type = df.detect_cds[[paste0(h, "_neg_strand_codon_type")]][k]
          
          relative_target_pos_range = rep(NA,3)
          if (h_codon_type %in% c("AAA", "ABA", "AAB", "BAB")){
            relative_target_pos_range = 4:6
          } else if (h_codon_type %in% c("BAA", "BBA")){
            relative_target_pos_range = 3:5
          } else if (h_codon_type %in% c("ABB")){
            relative_target_pos_range = 5:7
          } 
          
          h_seq = df.detect[[h]]
          temp[[paste0(h, "_target_base")]] = h_seq[detect_range+1]
          
          h_gRNA_strand = paste0(h, "_gRNA_strand")
          temp[[h_gRNA_strand]] = "+"
          
          final_seq = h_seq
          if (gRNA_strand == "neg"){
            final_seq = rev(unlist(base_complement[h_seq]))
            temp[[h_gRNA_strand]] = "-"
          }
          
          temp[[paste0(h, "_pos_strand_codon")]] = h_pos_strand_codon
          temp[[paste0(h, "_pos_strand_codon_type")]] = h_pos_strand_codon_type
          temp[[paste0(h, "_neg_strand_codon")]] = h_neg_strand_codon
          temp[[paste0(h, "_neg_strand_codon_type")]] = h_neg_strand_codon_type
          
          h_relative_target_pos = paste0(h, "_relative_target_pos")
          temp2 = data.frame(relative_target_pos_range)
          colnames(temp2) = h_relative_target_pos
          temp = cbind(temp, temp2)
          
          h_gRNA_seq = paste0(h, "_gRNA_seq")
          temp[[h_gRNA_seq]] = NA
          
          h_reporter_seq = paste0(h, "_reporter_seq")
          temp[[h_reporter_seq]] = NA
          
          if (!all(is.na(relative_target_pos_range))){
            for (m in 1:nrow(temp)){
              gRNA_start_ind = (detect_range+1)-(temp[[h_relative_target_pos]][m]-1)
              gRNA_seq = paste0(final_seq[gRNA_start_ind:(gRNA_start_ind+gRNA_len-1)], collapse = "")
              temp[[h_gRNA_seq]][m] = gRNA_seq
              reporter_seq = paste0(final_seq[(gRNA_start_ind-reporter_extend_len):(gRNA_start_ind+gRNA_len-1+reporter_extend_len)], collapse = "")
              temp[[h_reporter_seq]][m] = reporter_seq
            }
          }
        }
        
        df.out = rbind(df.out, temp)
      }
    }
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  df.out$is_gRNA_same = F
  ind = which(df.out$h1_gRNA_seq == df.out$h2_gRNA_seq)
  df.out$is_gRNA_same[ind] = T
  df.out$is_reporter_same = F
  ind2 = which(df.out$h1_reporter_seq == df.out$h2_reporter_seq)
  df.out$is_reporter_same[ind2] = T
  
  return(df.out)
}

get_hepG2_gRNA_CRISPRi <- function(df.main, detect_range=50, gRNA_len=20, gRNA_ct=5){
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg19)
  base_complement = list("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  
  # import hg19 genome
  genome <- BSgenome.Hsapiens.UCSC.hg19
  
  # get sequences in detect_range around variant on + strand
  temp = GRanges(seqnames = df.main$chr, 
                 ranges = IRanges(start = df.main$hg19_pos - detect_range, end = df.main$hg19_pos + detect_range),
                 strand = "+")
  df.main$ref_seq = getSeq(genome, temp)
  
  # output gRNAs and reporters
  df.out = c()
  pb = txtProgressBar(min = 0, max = nrow(df.main), initial = 0, style = 3)
  
  for (i in 1:nrow(df.main)){
    # derive haplotype 1 and haplotype 2
    chr = df.main$chr[i]
    pos = df.main$hg19_pos[i]
    ref_seq = df.main$ref_seq[i]
    df.detect = data.frame(hg19_pos = (pos-detect_range):(pos+detect_range), 
                           hg19_base = unlist(str_split(as.character(ref_seq), "")))
    
    ind = match(df.detect$hg19_pos, table = ls.phase[[chr]]$Pos)
    temp = ls.phase[[chr]][ind, 3:5]
    colnames(temp) = paste0("HepG2_", colnames(temp))
    df.detect = cbind(df.detect, temp)
    
    df.detect$h1 = df.detect$hg19_base
    df.detect$h2 = df.detect$hg19_base
    for (j in grep(pattern = "\\d+\\|\\d+", x = df.detect$HepG2_Phase)){ # consider only high-quality phase info, e.g. 0|1
      # when HepG2_Ref is more than 1 base, remove bases from following positions to compensate
      hepG2_ref_len = nchar(df.detect$HepG2_Ref[j])
      if (hepG2_ref_len > 1 & j < nrow(df.detect)){
        add_max = hepG2_ref_len-1
        if (add_max > nrow(df.detect)-j){
          add_max = nrow(df.detect)-j
        }
        df.detect$h1[j+1:add_max] = ""
        df.detect$h2[j+1:add_max] = ""
      }
      
      phase = as.numeric(unlist(str_split(df.detect$HepG2_Phase[j], pattern = "\\|")))
      hepG2_genotypes = c(df.detect$HepG2_Ref[j], unlist(str_split(df.detect$HepG2_Alt[j], pattern = ",")))
      df.detect$h1[j] = hepG2_genotypes[phase[1]+1]
      df.detect$h2[j] = hepG2_genotypes[phase[2]+1]
    }
    
    # recreate df.detect to ensure haplotype 1 and 2 have single base at each position
    for (h in c("h1", "h2")){
      h_corrected = paste0(h, "_corrected")
      df.detect[[h_corrected]] = NA
      
      temp = unlist(str_split(df.detect[[h]][1:detect_range], pattern = "")) # 50 bases before target position
      if (length(temp) < detect_range){
        temp = c(rep("", detect_range-length(temp)), temp)
      }
      df.detect[[h_corrected]][1:detect_range] = temp[(length(temp)-detect_range+1):length(temp)]
      
      temp = unlist(str_split(df.detect[[h]][(detect_range+1):(detect_range*2+1)], pattern = "")) # 51 bases on and after target position
      if (length(temp) < (detect_range+1)){
        temp = c(temp, rep("", detect_range+1-length(temp)))
      }
      df.detect[[h_corrected]][(detect_range+1):(detect_range*2+1)] = temp[1:(detect_range+1)]
    }
    
    # Design 5 gRNAs as close as possible to each variant
    # 23-nt gRNA sequence is NNNNNNNNNNNNNNNNNNNNNGG where last 3 nt is the PAM
    # can be in either strand
    # The 20-nt gRNA and 3-nt PAM is homozygous in the HepG2 genome
    # The 20-nt gRNA does NOT have a TTTT sequence in it
    
    df.candidate = c()
    df.detect$h1_corrected_shifted = c(df.detect$h1_corrected[-1],NA)
    temp = c()
    
    # for gRNAs on the + strand
    for (k in which(df.detect$h1_corrected == "G" & df.detect$h1_corrected_shifted == "G")){
      ind_start = (k-2)-gRNA_len+1
      ind_end = k-2
      if (ind_start > 0){
        if (all(df.detect$h1_corrected[ind_start:ind_end] == df.detect$h2_corrected[ind_start:ind_end])){
          gRNA_seq = paste0(df.detect$h1_corrected[ind_start:ind_end], collapse = "")
          if (!grepl(pattern = "TTTT", gRNA_seq)){
            temp = data.frame(target_id = df.main$target_id[i], chr = df.main$chr[i], target_hg19_pos = df.main$hg19_pos[i], 
                              gRNA_strand = "+", gRNA_start_pos = df.detect$hg19_pos[ind_start], target_distance_to_gRNA_start = abs((detect_range+1)-ind_start+1),
                              gRNA_seq)
            df.candidate = rbind(df.candidate, temp)
          }
        }
      }
    }
    
    # for gRNA on the - strand
    for (k in which(df.detect$h1_corrected == "C" & df.detect$h1_corrected_shifted == "C")){
      ind_start = (k+3)+gRNA_len-1
      ind_end = k+3
      if (ind_start <= nrow(df.detect)){
        if (all(df.detect$h1_corrected[ind_start:ind_end] == df.detect$h2_corrected[ind_start:ind_end])){
          gRNA_seq = paste0(base_complement[df.detect$h1_corrected[ind_start:ind_end]], collapse = "")
          if (!grepl(pattern = "TTTT", gRNA_seq)){
            temp = data.frame(target_id = df.main$target_id[i], chr = df.main$chr[i], target_hg19_pos = df.main$hg19_pos[i], 
                              gRNA_strand = "-", gRNA_start_pos = df.detect$hg19_pos[ind_start], target_distance_to_gRNA_start = abs((detect_range+1)-ind_start+1),
                              gRNA_seq)
            df.candidate = rbind(df.candidate, temp)
          }
        }
      }
    }
    
    if (!is.null(df.candidate)){
      df.candidate = df.candidate[order(df.candidate$target_distance_to_gRNA_start),] # order candidate gRNAs by distance to target
      # select top 5 gRNAs (if available)
      if (gRNA_ct <= nrow(df.candidate)){
        df.candidate = df.candidate[1:gRNA_ct,]
      }
      df.candidate = data.frame(gRNA_id = paste0(df.candidate$target_id, "_", 1:nrow(df.candidate)), df.candidate)
      df.out = rbind(df.out, df.candidate)
    }
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  return(df.out)
}

get_hepG2_pos_ctrl_gRNA_CRISPRi <- function(df.pos_ctrl, gRNA_ct = 12){
  df.pos_ctrl = df.pos_ctrl %>% mutate(gRNA_id = NA, gene = Input, chr = NA, gRNA_strand = `Strand of sgRNA`, cut_pos = `sgRNA 'Cut' Position`, 
                                       gRNA_PAM_start = NA, gRNA_PAM_end = NA, gRNA_seq = `sgRNA Sequence`) 
  ind = match(df.pos_ctrl$gene, table = gene_names)
  df.pos_ctrl$chr = df.raw$HG38.chrom[ind]
  
  # sgRNA on + strand, cut-off site at 18
  # sgRNA on - strand, cut-off site at 4
  ind = which(df.pos_ctrl$gRNA_strand == "+")
  df.pos_ctrl$gRNA_PAM_start[ind] = df.pos_ctrl$cut_pos[ind] - 17
  df.pos_ctrl$gRNA_PAM_end[ind] = df.pos_ctrl$cut_pos[ind] + 5
  ind = which(df.pos_ctrl$gRNA_strand == "-")
  df.pos_ctrl$gRNA_PAM_start[ind] = df.pos_ctrl$cut_pos[ind] - 6
  df.pos_ctrl$gRNA_PAM_end[ind] = df.pos_ctrl$cut_pos[ind] + 19
  
  # liftover from hg38 coordinates to hg19 coordinates
  temp = GRanges(seqnames = df.pos_ctrl$chr, 
                 ranges = IRanges(start = df.pos_ctrl$gRNA_PAM_start, end = df.pos_ctrl$gRNA_PAM_end),
                 strand = "+")
  ch = import.chain(system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain"))
  temp_hg19 = unlist(liftOver(temp, ch))
  df.pos_ctrl$gRNA_PAM_start_hg19 = start(temp_hg19)
  df.pos_ctrl$gRNA_PAM_end_hg19 = end(temp_hg19)
  
  # remove gRNAs where the gRNA seq or PAM contain heterozygous hegG2 sites
  ind_heterozygous = c()
  df.pos_ctrl$n_phase = 0
  for (i in 1:nrow(df.pos_ctrl)){
    chr = df.pos_ctrl$chr[i]
    ind2 = which(ls.phase[[chr]]$Pos >= df.pos_ctrl$gRNA_PAM_start_hg19[i] & ls.phase[[chr]]$Pos <= df.pos_ctrl$gRNA_PAM_end_hg19[i])
    df.pos_ctrl$n_phase[i] = length(ind2)
    if (length(ind2)>0){
      # if(any(ls.phase[[chr]]$Phase[ind2] == "0|1") | any(ls.phase[[chr]]$Phase[ind2] == "1|0") | any(ls.phase[[chr]]$Phase[ind2] == "0/1") | any(ls.phase[[chr]]$Phase[ind2] == "1/0")){
      #   ind_heterozygous = c(ind_heterozygous, i)
      # }
      ind_heterozygous = c(ind_heterozygous, i)
    }
  }
  
  df.pos_ctrl = df.pos_ctrl[-ind_heterozygous,]
  
  # create gRNA ids
  # select top 12 gRNA from each gene
  df.pos_ctrl2 = c()
  for (gene in gene_names){
    ind = which(df.pos_ctrl$gene == gene)
    if (length(ind)>0){
      if (length(ind) > gRNA_ct){
        ind = ind[1:gRNA_ct]
      }
      temp = df.pos_ctrl[ind,]
      temp$gRNA_id = paste0("gPosCtrl_", gene, "_", 1:length(ind))
    }
    
    df.pos_ctrl2 = rbind(df.pos_ctrl2, temp)
  }
  
  df.pos_ctrl2 = df.pos_ctrl2 %>% dplyr::select(gRNA_id, gene, chr, gRNA_strand, cut_pos, gRNA_PAM_start, gRNA_PAM_end, gRNA_PAM_start_hg19, gRNA_PAM_end_hg19, gRNA_seq)
  return(df.pos_ctrl2)
}

get_splice_pos_v2 <- function(df.mane_hg19_exon){
  df.main = c()
  for (i in 1:length(df.mane_hg19_exon)){
    gene_name = df.mane_hg19_exon$gene_name[i]
    exon_number = df.mane_hg19_exon$exon_number[i]
    total_exon_number = df.mane_hg19_exon$total_exon_number[i]
    chr = as.character(seqnames(df.mane_hg19_exon)[i])
    gene_strand = as.character(strand(df.mane_hg19_exon)[i])
    exon_start = start(df.mane_hg19_exon)[i]
    exon_end = end(df.mane_hg19_exon)[i]
    
    if (gene_strand == "+"){
      gRNA_sense = c("sense", "antisense")
      gRNA_strand = c("+", "-")
      # first exon doesn't have acceptor (upstream) site
      if (exon_number > 1){
        temp = data.frame(target_id = paste0(gene_name, "_", exon_number, "_acceptorA_", gRNA_sense), 
                          chr, gene_strand, gRNA_sense, gRNA_strand, hg19_pos = exon_start-2)
        df.main = rbind(df.main, temp)
        temp = data.frame(target_id = paste0(gene_name, "_", exon_number, "_acceptorG_", gRNA_sense), 
                          chr, gene_strand, gRNA_sense, gRNA_strand, hg19_pos = exon_start-1)
        df.main = rbind(df.main, temp)
      }
      # last exon doesn't have donor (downstream) site
      if (exon_number < total_exon_number){
        temp = data.frame(target_id = paste0(gene_name, "_", exon_number, "_donorG_", gRNA_sense), 
                          chr, gene_strand, gRNA_sense, gRNA_strand, hg19_pos = exon_end+1)
        df.main = rbind(df.main, temp)
        temp = data.frame(target_id = paste0(gene_name, "_", exon_number, "_donorT_", gRNA_sense), 
                          chr, gene_strand, gRNA_sense, gRNA_strand, hg19_pos = exon_end+2)
        df.main = rbind(df.main, temp)
      }
    } else {
      gRNA_sense = c("sense", "antisense")
      gRNA_strand = c("-", "+")
      # first exon doesn't have acceptor (upstream) site
      if (exon_number > 1){
        temp = data.frame(target_id = paste0(gene_name, "_", exon_number, "_acceptorA_", gRNA_sense), 
                          chr, gene_strand, gRNA_sense, gRNA_strand, hg19_pos = exon_end+2)
        df.main = rbind(df.main, temp)
        temp = data.frame(target_id = paste0(gene_name, "_", exon_number, "_acceptorG_", gRNA_sense), 
                          chr, gene_strand, gRNA_sense, gRNA_strand, hg19_pos = exon_end+1)
        df.main = rbind(df.main, temp)
      }
      # last exon doesn't have donor (downstream) site
      if (exon_number < total_exon_number){
        temp = data.frame(target_id = paste0(gene_name, "_", exon_number, "_donorG_", gRNA_sense), 
                          chr, gene_strand, gRNA_sense, gRNA_strand, hg19_pos = exon_start-1)
        df.main = rbind(df.main, temp)
        temp = data.frame(target_id = paste0(gene_name, "_", exon_number, "_donorT_", gRNA_sense), 
                          chr, gene_strand, gRNA_sense, gRNA_strand, hg19_pos = exon_start-2)
        df.main = rbind(df.main, temp)
      }
    }
  }
  
  return(df.main)
}

get_hegG2_gRNA_reporter_splice_sites_wrapper <- function(df.mane_hg19_exon, df.main_all, ls.phase, relative_target_pos_range=c(4), 
                                                         gRNA_len=20, reporter_extend_len=6, detect_range=50){
  # output gRNAs and reporters
  df.out_all = c()
  pb = txtProgressBar(min = 0, max = length(df.mane_hg19_cds), initial = 0, style = 3)
  df.exon_meta = mcols(df.mane_hg19_exon)
  
  for (i in 1:length(df.mane_hg19_exon)){
    # for each exon:
    # construct a df.detect_exon to hold hepG2 phase info
    df.detect_exon_start = data.frame(hg19_pos = df.exon_meta$exon_start_extend_start[i]:df.exon_meta$exon_start_extend_end[i], 
                                      hg19_base = unlist(stringr::str_split(df.exon_meta$exon_start_extend_seq[i], "")))
    df.detect_exon_end = data.frame(hg19_pos = df.exon_meta$exon_end_extend_start[i]:df.exon_meta$exon_end_extend_end[i], 
                                    hg19_base = unlist(stringr::str_split(df.exon_meta$exon_end_extend_seq[i], "")))
    df.detect_exon = rbind(df.detect_exon_start, df.detect_exon_end)
    
    
    exon_strand = as.character(strand(df.mane_hg19_exon)[i])
    df.detect_exon_start$relative_to_exon_start = df.detect_exon_start$hg19_pos - start(df.mane_hg19_exon)[i]
    df.detect_exon_end$relative_to_exon_end = df.detect_exon_end$hg19_pos - end(df.mane_hg19_exon)[i]
    
    
    chr = as.character(seqnames(df.mane_hg19_exon)[i])
    ind = match(df.detect_exon$hg19_pos, table = ls.phase[[chr]]$Pos)
    temp = ls.phase[[chr]][ind, 3:5]
    colnames(temp) = paste0("HepG2_", colnames(temp))
    df.detect_exon = cbind(df.detect_exon, temp)  
    
    # for each exon:
    # construct a df.main for target_id, target chr/pos/strand/base
    exon_start_extend_start = df.exon_meta$exon_start_extend_start[i]
    exon_end_extend_end = df.exon_meta$exon_end_extend_end[i]
    ind = which(df.main_all$chr == chr & df.main_all$hg19_pos >= exon_start_extend_start & df.main_all$hg19_pos <= exon_end_extend_end)
    df.main = df.main_all[ind,]
    
    df.out = get_hegG2_gRNA_reporter_cds_tiling(df.main, df.detect_exon, relative_target_pos_range, gRNA_len, reporter_extend_len, detect_range)
    df.out_all = rbind(df.out_all, df.out)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  
  df.out_all$is_gRNA_same = F
  ind = which(df.out_all$h1_gRNA_seq == df.out_all$h2_gRNA_seq)
  df.out_all$is_gRNA_same[ind] = T
  df.out_all$is_reporter_same = F
  ind2 = which(df.out_all$h1_reporter_seq == df.out_all$h2_reporter_seq)
  df.out_all$is_reporter_same[ind2] = T
  
  return(df.out_all)
}

get_splice_pos <- function(gene_list, df.mane){
  df.main = c()
  for (gene_name in gene_list){
    ind = which(df.mane$V19 == gene_name)[1] # sometime there's two entries for the same gene, mane select and mane plus clinical
    strand = df.mane$V6[ind]
    
    temp = unlist(str_split(df.mane$V12[ind], pattern = ","))
    exon_start = df.mane$V2[ind] + as.numeric(temp[-length(temp)])
    temp = unlist(str_split(df.mane$V11[ind], pattern = ","))
    exon_end = exon_start + as.numeric(temp[-length(temp)])
    
    exon_number = df.mane$V10[ind]
    # first exon doesn't have acceptor (upstream) site
    temp = data.frame(target_id = paste0(gene_name, "_", 2:exon_number, "_acceptorA"), 
                      chr = df.mane$V1[ind], strand, hg38_pos = exon_start[2:exon_number] - 1)
    df.main = rbind(df.main, temp)
    temp = data.frame(target_id = paste0(gene_name, "_", 2:exon_number, "_acceptorG"), 
                      chr = df.mane$V1[ind], strand, hg38_pos = exon_start[2:exon_number])
    df.main = rbind(df.main, temp)
    # last exon doesn't have donor (downstream) site
    temp = data.frame(target_id = paste0(gene_name, "_", 1:(exon_number-1), "_donorG"), 
                      chr = df.mane$V1[ind], strand, hg38_pos = exon_end[1:(exon_number-1)] + 1)
    df.main = rbind(df.main, temp)
    temp = data.frame(target_id = paste0(gene_name, "_", 1:(exon_number-1), "_donorT"), 
                      chr = df.mane$V1[ind], strand, hg38_pos = exon_end[1:(exon_number-1)] + 2)
    df.main = rbind(df.main, temp)
  }
  
  # liftover from hg38 coordinates to hg19
  pos.hg38 = GRanges(seqnames = df.main$chr, ranges = as.character(df.main$hg38_pos))
  ch = import.chain(system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain"))
  pos.hg19 = unlist(liftOver(pos.hg38, ch))
  df.main$hg19_pos = start(pos.hg19)
  
  df.main$hg19_Ref = NA
  for (i in 1:nrow(df.main)){
    chr = df.main$chr[i]
    pos = df.main$hg19_pos[i]
    df.main$hg19_Ref[i] = substr(hg19[[chr]], start = pos, stop = pos)
  }
  
  return(df.main)
}



get_cds_tiling_pos <- function(df.mane_hg19_cds, tiling_step=1, cds_start_flank=0, cds_end_flank=0){
  df.cds_meta = mcols(df.mane_hg19_cds)
  
  df.main = c()
  for (i in 1:length(df.mane_hg19_cds)){
    chr = as.character(seqnames(df.mane_hg19_cds)[i])
    gene_name = df.cds_meta$gene_name[i]
    gene_strand = as.character(strand(df.mane_hg19_cds)[i])
    cds_start = start(df.mane_hg19_cds)[i]
    cds_end = end(df.mane_hg19_cds)[i]
    exon_number = df.cds_meta$exon_number[i]
    
    # for each cds:
    # construct a df.main for target_id, target chr/pos/strand/base
    if (gene_strand == "+"){
      gene_strand_factor = 1
      cds_start_actual = cds_start
      cds_end_actual = cds_end
      
      gRNA_sense = "sense"
      gRNA_sense_factor = 1
      gRNA_strand = "+"
      pos_start = cds_start_actual - cds_start_flank*gene_strand_factor*gRNA_sense_factor
      pos_end = cds_end_actual + cds_end_flank*gene_strand_factor*gRNA_sense_factor
      pos_range = seq(pos_start, pos_end, tiling_step*gene_strand_factor*gRNA_sense_factor)
      pos_range = unique(c(pos_range, pos_end))
      temp = data.frame(target_id = paste(gene_name, exon_number, pos_range, gRNA_sense, sep = "_"), 
                        chr, gene_strand, gRNA_sense, gRNA_strand, hg19_pos = pos_range,
                        relative_to_cds_start = (pos_range-cds_start_actual)*gene_strand_factor, 
                        relative_to_cds_end = (pos_range-cds_end_actual)*gene_strand_factor,
                        mane_hg19_cds_ind = i)
      df.main = rbind(df.main, temp)
      
      gRNA_sense = "antisense"
      gRNA_sense_factor = -1
      gRNA_strand = "-"
      pos_start = cds_end_actual - cds_start_flank*gene_strand_factor*gRNA_sense_factor
      pos_end = cds_start_actual + cds_end_flank*gene_strand_factor*gRNA_sense_factor
      pos_range = seq(pos_start, pos_end, tiling_step*gene_strand_factor*gRNA_sense_factor)
      pos_range = unique(c(pos_range, pos_end))
      temp = data.frame(target_id = paste(gene_name, exon_number, pos_range, gRNA_sense, sep = "_"), 
                        chr, gene_strand, gRNA_sense, gRNA_strand, hg19_pos = pos_range, 
                        relative_to_cds_start = (pos_range-cds_start_actual)*gene_strand_factor, 
                        relative_to_cds_end = (pos_range-cds_end_actual)*gene_strand_factor,
                        mane_hg19_cds_ind = i)
      df.main = rbind(df.main, temp)
    } else {
      gene_strand_factor = -1
      cds_start_actual = cds_end
      cds_end_actual = cds_start
      
      gRNA_sense = "sense"
      gRNA_sense_factor = 1
      gRNA_strand = "-"
      pos_start = cds_start_actual - cds_start_flank*gene_strand_factor*gRNA_sense_factor
      pos_end = cds_end_actual + cds_end_flank*gene_strand_factor*gRNA_sense_factor
      pos_range = seq(pos_start, pos_end, tiling_step*gene_strand_factor*gRNA_sense_factor)
      pos_range = unique(c(pos_range, pos_end))
      temp = data.frame(target_id = paste(gene_name, exon_number, pos_range, gRNA_sense, sep = "_"), 
                        chr, gene_strand, gRNA_sense, gRNA_strand, hg19_pos = pos_range, 
                        relative_to_cds_start = (pos_range-cds_start_actual)*gene_strand_factor, 
                        relative_to_cds_end = (pos_range-cds_end_actual)*gene_strand_factor,
                        mane_hg19_cds_ind = i)
      df.main = rbind(df.main, temp)
      
      gRNA_sense = "antisense"
      gRNA_sense_factor = -1
      gRNA_strand = "+"
      pos_start = cds_end_actual - cds_start_flank*gene_strand_factor*gRNA_sense_factor
      pos_end = cds_start_actual + cds_end_flank*gene_strand_factor*gRNA_sense_factor
      pos_range = seq(pos_start, pos_end, tiling_step*gene_strand_factor*gRNA_sense_factor)
      pos_range = unique(c(pos_range, pos_end))
      temp = data.frame(target_id = paste(gene_name, exon_number, pos_range, gRNA_sense, sep = "_"), 
                        chr, gene_strand, gRNA_sense, gRNA_strand, hg19_pos = pos_range, 
                        relative_to_cds_start = (pos_range-cds_start_actual)*gene_strand_factor, 
                        relative_to_cds_end = (pos_range-cds_end_actual)*gene_strand_factor,
                        mane_hg19_cds_ind = i)
      df.main = rbind(df.main, temp)
    }
  }
  
  return(df.main)
}

# get_exon_tiling_pos <- function(gene_list, df.mane, exon_start_flank = 2, exon_end_flank = 2, tiling_step = 1){
#   df.main = c()
#   pb = txtProgressBar(min = 0, max = length(gene_list), initial = 0, style = 3)
#   j=0
#   
#   for (gene_name in gene_list){
#     ind = which(df.mane$V19 == gene_name)[1] # sometime there's two entries for the same gene, mane select and mane plus clinical
#     gene_strand = df.mane$V6[ind]
#     
#     temp = unlist(str_split(df.mane$V12[ind], pattern = ","))
#     exon_start = df.mane$V2[ind] + as.numeric(temp[-length(temp)])
#     temp = unlist(str_split(df.mane$V11[ind], pattern = ","))
#     exon_end = exon_start + as.numeric(temp[-length(temp)])
#     exon_number = df.mane$V10[ind]
#     
#     for (i in 1:exon_number){
#       # mane exon range: actual_exon_start+1 to actual_exon_end
#       if (gene_strand == "+"){
#         gene_strand_factor = 1
#         exon_start_actual = exon_start[i]+1
#         exon_end_actual = exon_end[i]
#         
#         gRNA_sense = "sense"
#         gRNA_sense_factor = 1
#         gRNA_strand = "+"
#         pos_start = exon_start_actual - exon_start_flank*gene_strand_factor*gRNA_sense_factor
#         pos_end = exon_end_actual + exon_end_flank*gene_strand_factor*gRNA_sense_factor
#         pos_range = seq(pos_start, pos_end, tiling_step*gene_strand_factor*gRNA_sense_factor)
#         pos_range = unique(c(pos_range, pos_end))
#         temp = data.frame(target_id = paste(gene_name, i, pos_range, gRNA_sense, sep = "_"), 
#                           chr = df.mane$V1[ind], gene_strand, gRNA_sense, gRNA_strand, hg38_pos = pos_range, 
#                           relative_to_exon_start = (pos_range-exon_start_actual)*gene_strand_factor, 
#                           relative_to_exon_end = (pos_range-exon_end_actual)*gene_strand_factor)
#         df.main = rbind(df.main, temp)
#         
#         gRNA_sense = "antisense"
#         gRNA_sense_factor = -1
#         gRNA_strand = "-"
#         pos_start = exon_end_actual - exon_start_flank*gene_strand_factor*gRNA_sense_factor
#         pos_end = exon_start_actual + exon_end_flank*gene_strand_factor*gRNA_sense_factor
#         pos_range = seq(pos_start, pos_end, tiling_step*gene_strand_factor*gRNA_sense_factor)
#         pos_range = unique(c(pos_range, pos_end))
#         temp = data.frame(target_id = paste(gene_name, i, pos_range, gRNA_sense, sep = "_"), 
#                           chr = df.mane$V1[ind], gene_strand, gRNA_sense, gRNA_strand, hg38_pos = pos_range, 
#                           relative_to_exon_start = (pos_range-exon_start_actual)*gene_strand_factor, 
#                           relative_to_exon_end = (pos_range-exon_end_actual)*gene_strand_factor)
#         df.main = rbind(df.main, temp)
#       } else {
#         gene_strand_factor = -1
#         exon_start_actual = exon_end[i]
#         exon_end_actual = exon_start[i]+1
#         
#         gRNA_sense = "sense"
#         gRNA_sense_factor = 1
#         gRNA_strand = "-"
#         pos_start = exon_start_actual - exon_start_flank*gene_strand_factor*gRNA_sense_factor
#         pos_end = exon_end_actual + exon_end_flank*gene_strand_factor*gRNA_sense_factor
#         pos_range = seq(pos_start, pos_end, tiling_step*gene_strand_factor*gRNA_sense_factor)
#         pos_range = unique(c(pos_range, pos_end))
#         temp = data.frame(target_id = paste(gene_name, i, pos_range, gRNA_sense, sep = "_"), 
#                           chr = df.mane$V1[ind], gene_strand, gRNA_sense, gRNA_strand, hg38_pos = pos_range, 
#                           relative_to_exon_start = (pos_range-exon_start_actual)*gene_strand_factor, 
#                           relative_to_exon_end = (pos_range-exon_end_actual)*gene_strand_factor)
#         df.main = rbind(df.main, temp)
#         
#         gRNA_sense = "antisense"
#         gRNA_sense_factor = -1
#         gRNA_strand = "+"
#         pos_start = exon_end_actual - exon_start_flank*gene_strand_factor*gRNA_sense_factor
#         pos_end = exon_start_actual + exon_end_flank*gene_strand_factor*gRNA_sense_factor
#         pos_range = seq(pos_start, pos_end, tiling_step*gene_strand_factor*gRNA_sense_factor)
#         pos_range = unique(c(pos_range, pos_end))
#         temp = data.frame(target_id = paste(gene_name, i, pos_range, gRNA_sense, sep = "_"), 
#                           chr = df.mane$V1[ind], gene_strand, gRNA_sense, gRNA_strand, hg38_pos = pos_range, 
#                           relative_to_exon_start = (pos_range-exon_start_actual)*gene_strand_factor, 
#                           relative_to_exon_end = (pos_range-exon_end_actual)*gene_strand_factor)
#         df.main = rbind(df.main, temp)
#       }
#     }
#     
#     j=j+1
#     setTxtProgressBar(pb, j)
#   }
#   close(pb)
#   
#   pos.hg38 = GRanges(seqnames = df.main$chr, ranges = as.character(df.main$hg38_pos))
#   ch = import.chain(system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain"))
#   pos.hg19 = unlist(liftOver(pos.hg38, ch))
#   df.main$hg19_pos = start(pos.hg19)
#   
#   df.main$hg19_Ref = NA
#   pb = txtProgressBar(min = 0, max = nrow(df.main), initial = 0, style = 3)
#   
#   for (i in 1:nrow(df.main)){
#     chr = df.main$chr[i]
#     pos = df.main$hg19_pos[i]
#     df.main$hg19_Ref[i] = substr(hg19[[chr]], start = pos, stop = pos)
#     setTxtProgressBar(pb, i)
#   }
#   close(pb)
#   
#   return(df.main)
# }

parse_out <- function(df.out, gRNA_prefix="gRNA", collapse_haplotypes="gRNA"){
  # each line contains one gRNA-reporter pair
  # unique_id, chr, target_hg19_pos, target_relative_pos, editor, haplotype, gRNA_seq, reporter_seq
  strand_name = c("pos", "neg")
  names(strand_name) = c("+", "-")
  
  df.out2 = c()
  pb = txtProgressBar(min = 0, max = nrow(df.main), initial = 0, style = 3)
  
  for (i in 1:nrow(df.out)){
    temp = data.frame(unique_id = NA, df.out[i,1:4], target_base=df.out$h1_target_base[i], haplotype=1, 
                      editor=ifelse(df.out$h1_target_base[i] %in% c("A", "T"), yes = "ABE", no = "CBE"), gRNA_strand=df.out$h1_gRNA_strand[i], 
                      gRNA_seq=df.out$h1_gRNA_seq[i], reporter_seq = df.out$h1_reporter_seq[i])
    
    if ((collapse_haplotypes == "gRNA" & df.out$is_gRNA_same[i] == F) | (collapse_haplotypes == "reporter" & df.out$is_reporter_same[i] == F)){
      temp2 = data.frame(unique_id = NA, df.out[i,1:4], target_base=df.out$h2_target_base[i], haplotype=2, 
                         editor=ifelse(df.out$h2_target_base[i] %in% c("A", "T"), yes = "ABE", no = "CBE"), gRNA_strand=df.out$h2_gRNA_strand[i], 
                         gRNA_seq=df.out$h2_gRNA_seq[i], reporter_seq = df.out$h2_reporter_seq[i])
      temp = rbind(temp, temp2)
    } else if (!collapse_haplotypes %in% c("gRNA", "reporter")){
      simpleError(message = "collapse_haplotypes can only 'gRNA' or 'reporter'!")
    }
    
    temp$unique_id = paste0(gRNA_prefix, "__", temp$target_id, "__", temp$editor, "_", strand_name[temp$gRNA_strand], "_p", temp$relative_target_pos, "_h", temp$haplotype)
    
    df.out2 = rbind(df.out2, temp)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  return(df.out2)
}

parse_tiling_out <- function(df.out, col_to_keep, gRNA_prefix="gRNA", collapse_haplotypes="gRNA", homozygous_only=T){
  # each line contains one gRNA-reporter pair
  # unique_id, chr, target_hg19_pos, target_relative_pos, editor, haplotype, gRNA_seq, reporter_seq
  strand_name = c("pos", "neg")
  names(strand_name) = c("+", "-")
  
  if (!collapse_haplotypes %in% c("gRNA", "reporter")){
    simpleError(message = "collapse_haplotypes can only 'gRNA' or 'reporter'!")
  }
  
  ind = which((collapse_haplotypes == "gRNA" & df.out$is_gRNA_same == T) | (collapse_haplotypes == "reporter" & df.out$is_reporter_same == T))
  temp = data.frame(unique_id = NA, df.out[ind, col_to_keep], target_base=df.out$h1_target_base[ind], haplotype=1, 
                    editor=ifelse(df.out$h1_target_base[ind] %in% c("A", "T"), yes = "ABE", no = "CBE"), 
                    # gRNA_strand=df.out$h1_gRNA_strand[ind],
                    gRNA_seq=df.out$h1_gRNA_seq[ind], reporter_seq = df.out$h1_reporter_seq[ind])
  
  if (!homozygous_only){
    ind = which((collapse_haplotypes == "gRNA" & df.out$is_gRNA_same == F) | (collapse_haplotypes == "reporter" & df.out$is_reporter_same == F))
    temp2 = data.frame(unique_id = NA, df.out[ind, col_to_keep], target_base=df.out$h1_target_base[ind], haplotype=1, 
                       editor=ifelse(df.out$h1_target_base[ind] %in% c("A", "T"), yes = "ABE", no = "CBE"), 
                       # gRNA_strand=df.out$h1_gRNA_strand[ind], 
                       gRNA_seq=df.out$h1_gRNA_seq[ind], reporter_seq = df.out$h1_reporter_seq[ind])
    
    ind = which((collapse_haplotypes == "gRNA" & df.out$is_gRNA_same == F) | (collapse_haplotypes == "reporter" & df.out$is_reporter_same == F))
    temp3 = data.frame(unique_id = NA, df.out[ind, col_to_keep], target_base=df.out$h2_target_base[ind], haplotype=2, 
                       editor=ifelse(df.out$h2_target_base[ind] %in% c("A", "T"), yes = "ABE", no = "CBE"), 
                       # gRNA_strand=df.out$h2_gRNA_strand[ind], 
                       gRNA_seq=df.out$h2_gRNA_seq[ind], reporter_seq = df.out$h2_reporter_seq[ind])
    temp = rbind(temp, temp2, temp3)
  }
  
  temp$unique_id = paste0(gRNA_prefix, "__", temp$target_id, "__", temp$editor, "_", strand_name[temp$gRNA_strand], "_p", temp$relative_target_pos, "_h", temp$haplotype)
  
  return(temp)
}

parse_splice_site_out <- function(df.out, col_to_keep, gRNA_prefix="gControl", collapse_haplotypes="gRNA", homozygous_only=T){
  # each line contains one gRNA-reporter pair
  # unique_id, chr, target_hg19_pos, target_relative_pos, editor, haplotype, gRNA_seq, reporter_seq
  strand_name = c("pos", "neg")
  names(strand_name) = c("+", "-")
  
  if (!collapse_haplotypes %in% c("gRNA", "reporter")){
    simpleError(message = "collapse_haplotypes can only 'gRNA' or 'reporter'!")
  }
  
  ind = which((collapse_haplotypes == "gRNA" & df.out$is_gRNA_same == T) | (collapse_haplotypes == "reporter" & df.out$is_reporter_same == T))
  temp = data.frame(unique_id = NA, df.out[ind, col_to_keep], target_base=df.out$h1_target_base[ind], haplotype=1, 
                    # gRNA_strand=df.out$h1_gRNA_strand[ind],
                    gRNA_seq=df.out$h1_gRNA_seq[ind], reporter_seq = df.out$h1_reporter_seq[ind])
  
  if (!homozygous_only){
    ind = which((collapse_haplotypes == "gRNA" & df.out$is_gRNA_same == F) | (collapse_haplotypes == "reporter" & df.out$is_reporter_same == F))
    temp2 = data.frame(unique_id = NA, df.out[ind, col_to_keep], target_base=df.out$h1_target_base[ind], haplotype=1, 
                       # gRNA_strand=df.out$h1_gRNA_strand[ind], 
                       gRNA_seq=df.out$h1_gRNA_seq[ind], reporter_seq = df.out$h1_reporter_seq[ind])
    
    ind = which((collapse_haplotypes == "gRNA" & df.out$is_gRNA_same == F) | (collapse_haplotypes == "reporter" & df.out$is_reporter_same == F))
    temp3 = data.frame(unique_id = NA, df.out[ind, col_to_keep], target_base=df.out$h2_target_base[ind], haplotype=2, 
                       # gRNA_strand=df.out$h2_gRNA_strand[ind], 
                       gRNA_seq=df.out$h2_gRNA_seq[ind], reporter_seq = df.out$h2_reporter_seq[ind])
    temp = rbind(temp, temp2, temp3)
  }
  
  temp$unique_id = paste0(gRNA_prefix, "__", temp$target_id, "__", strand_name[temp$gRNA_strand], "_p", temp$relative_target_pos, "_h", temp$haplotype)
  
  return(temp)
}

parse_UKB_variant_out <- function(df.out, col_to_keep, gRNA_prefix="gRNA", collapse_haplotypes="gRNA", homozygous_only=T){
  # each line contains one gRNA-reporter pair
  # unique_id, chr, target_hg19_pos, target_relative_pos, editor, haplotype, gRNA_seq, reporter_seq
  strand_name = c("pos", "neg")
  names(strand_name) = c("+", "-")
  
  if (!collapse_haplotypes %in% c("gRNA", "reporter")){
    simpleError(message = "collapse_haplotypes can only 'gRNA' or 'reporter'!")
  }
  
  ind = which((collapse_haplotypes == "gRNA" & df.out$is_gRNA_same == T) | (collapse_haplotypes == "reporter" & df.out$is_reporter_same == T))
  temp = data.frame(unique_id = NA, df.out[ind, col_to_keep], target_base=df.out$h1_target_base[ind], 
                    gRNA_strand=df.out$h1_gRNA_strand[ind], 
                    pos_strand_codon = df.out$h1_pos_strand_codon[ind], pos_strand_codon_type = df.out$h1_pos_strand_codon_type[ind],
                    neg_strand_codon = df.out$h1_neg_strand_codon[ind], neg_strand_codon_type = df.out$h1_neg_strand_codon_type[ind],
                    relative_target_pos = df.out$h1_relative_target_pos[ind], haplotype=1, 
                    gRNA_seq=df.out$h1_gRNA_seq[ind], reporter_seq = df.out$h1_reporter_seq[ind])
  
  if (!homozygous_only){
    ind = which((collapse_haplotypes == "gRNA" & df.out$is_gRNA_same == F) | (collapse_haplotypes == "reporter" & df.out$is_reporter_same == F))
    temp2 = data.frame(unique_id = NA, df.out[ind, col_to_keep], target_base=df.out$h1_target_base[ind], 
                       gRNA_strand=df.out$h1_gRNA_strand[ind], 
                       pos_strand_codon = df.out$h1_pos_strand_codon[ind], pos_strand_codon_type = df.out$h1_pos_strand_codon_type[ind],
                       neg_strand_codon = df.out$h1_neg_strand_codon[ind], neg_strand_codon_type = df.out$h1_neg_strand_codon_type[ind],
                       relative_target_pos = df.out$h1_relative_target_pos[ind], haplotype=1, 
                       gRNA_seq=df.out$h1_gRNA_seq[ind], reporter_seq = df.out$h1_reporter_seq[ind])
    
    ind = which((collapse_haplotypes == "gRNA" & df.out$is_gRNA_same == F) | (collapse_haplotypes == "reporter" & df.out$is_reporter_same == F))
    temp3 = data.frame(unique_id = NA, df.out[ind, col_to_keep], target_base=df.out$h2_target_base[ind],  
                       gRNA_strand=df.out$h2_gRNA_strand[ind], codon = df.out$h2_codon[ind], 
                       pos_strand_codon = df.out$h2_pos_strand_codon[ind], pos_strand_codon_type = df.out$h2_pos_strand_codon_type[ind],
                       neg_strand_codon = df.out$h2_neg_strand_codon[ind], neg_strand_codon_type = df.out$h2_neg_strand_codon_type[ind],
                       relative_target_pos = df.out$h2_relative_target_pos[ind], haplotype=2, 
                       gRNA_seq=df.out$h2_gRNA_seq[ind], reporter_seq = df.out$h2_reporter_seq[ind])
    temp = rbind(temp, temp2, temp3)
  }
  
  temp$unique_id = paste0(gRNA_prefix, "__", temp$target_id, "__", strand_name[temp$gRNA_strand], "_p", temp$relative_target_pos, "_h", temp$haplotype)
  
  return(temp)
}

get_hegG2_gRNA_reporter_v2 <- function(df.main, hg19, detect_range, relative_target_pos_range, gRNA_len, reporter_extend_len){
  base_complement = list("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  
  # extract reference seq around target site
  df.main$ref_seq = NA
  for (i in 1:nrow(df.main)){
    chr = df.main$chr[i]
    pos = df.main$hg19_pos[i]
    df.main$n_phase[i] = sum((pos-detect_range):(pos+detect_range) %in% ls.phase[[chr]]$Pos)
    df.main$ref_seq[i] = substr(hg19[[chr]], start = pos-detect_range, stop = pos+detect_range)
  }
  
  # output gRNAs and reporters
  df.out = c()
  ls.detect = list()
  pb = txtProgressBar(min = 0, max = nrow(df.main), initial = 0, style = 3)
  
  for (i in 1:nrow(df.main)){
    # derive haplotype 1 and haplotype 2
    chr = df.main$chr[i]
    pos = df.main$hg19_pos[i]
    df.detect = data.frame(hg19_pos = (pos-detect_range):(pos+detect_range), 
                           hg19_base = unlist(str_split(df.main$ref_seq[i], "")))
    
    ind = match(df.detect$hg19_pos, table = ls.phase[[chr]]$Pos)
    temp = ls.phase[[chr]][ind, 3:5]
    colnames(temp) = paste0("HepG2_", colnames(temp))
    df.detect = cbind(df.detect, temp)
    
    df.detect$h1 = df.detect$hg19_base
    df.detect$h2 = df.detect$hg19_base
    for (j in grep(pattern = "\\d+\\|\\d+", x = df.detect$HepG2_Phase)){ # consider only high-quality phase info, e.g. 0|1
      # when HepG2_Ref is more than 1 base, remove bases from following positions to compensate
      hepG2_ref_len = nchar(df.detect$HepG2_Ref[j])
      if (hepG2_ref_len > 1){
        df.detect$h1[j+1:(hepG2_ref_len-1)] = ""
        df.detect$h2[j+1:(hepG2_ref_len-1)] = ""
      }
      
      phase = as.numeric(unlist(str_split(df.detect$HepG2_Phase[j], pattern = "\\|")))
      hepG2_genotypes = c(df.detect$HepG2_Ref[j], unlist(str_split(df.detect$HepG2_Alt[j], pattern = ",")))
      df.detect$h1[j] = hepG2_genotypes[phase[1]+1]
      df.detect$h2[j] = hepG2_genotypes[phase[2]+1]
    }
    
    # recreate df.detect to ensure haplotype 1 and 2 have single base at each position
    for (h in c("h1", "h2")){
      h_corrected = paste0(h, "_corrected")
      df.detect[[h_corrected]] = NA
      
      temp = unlist(str_split(df.detect[[h]][1:detect_range], pattern = "")) # 50 bases before target position
      if (length(temp) < detect_range){
        temp = c(rep("", detect_range-length(temp)), temp)
      }
      df.detect[[h_corrected]][1:detect_range] = temp[(length(temp)-detect_range+1):length(temp)]
      
      temp = unlist(str_split(df.detect[[h]][(detect_range+1):(detect_range*2+1)], pattern = "")) # 51 bases on and after target position
      if (length(temp) < (detect_range+1)){
        temp = c(temp, rep("", detect_range+1-length(temp)))
      }
      df.detect[[h_corrected]][(detect_range+1):(detect_range*2+1)] = temp[1:(detect_range+1)]
    }
    
    # save df.detect
    ls.detect[[df.main$target_id[i]]] = df.detect
    
    # output gRNAs and reporters
    temp = data.frame(target_id = df.main$target_id[i], chr = df.main$chr[i], target_hg19_pos = df.main$hg19_pos[i], relative_target_pos = relative_target_pos_range)
    for (h in c("h1", "h2")){
      h_corrected = paste0(h, "_corrected")
      h_seq = df.detect[[h_corrected]]
      
      h_target_base = paste0(h, "_target_base")
      temp[[h_target_base]] = h_seq[detect_range+1]
      
      h_gRNA_strand = paste0(h, "_gRNA_strand")
      temp[[h_gRNA_strand]] = "+"
      
      h_gRNA_seq = paste0(h, "_gRNA_seq")
      temp[[h_gRNA_seq]] = NA
      
      h_reporter_seq = paste0(h, "_reporter_seq")
      temp[[h_reporter_seq]] = NA
      
      if (h_seq[detect_range+1] %in% c("A", "C")){
        final_seq = h_seq
        
      } else {
        final_seq = rev(unlist(base_complement[h_seq]))
        temp[[h_gRNA_strand]] = "-"
      }
      
      for (m in 1:nrow(temp)){
        gRNA_start_ind = (detect_range+1)-(temp$relative_target_pos[m]-1)
        gRNA_seq = paste0(final_seq[gRNA_start_ind:(gRNA_start_ind+gRNA_len-1)], collapse = "")
        temp[[h_gRNA_seq]][m] = gRNA_seq
        reporter_seq = paste0(final_seq[(gRNA_start_ind-reporter_extend_len):(gRNA_start_ind+gRNA_len-1+reporter_extend_len)], collapse = "")
        temp[[h_reporter_seq]][m] = reporter_seq
      }
    }
    
    df.out = rbind(df.out, temp)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  df.out$is_gRNA_same = F
  ind = which(df.out$h1_gRNA_seq == df.out$h2_gRNA_seq)
  df.out$is_gRNA_same[ind] = T
  df.out$is_reporter_same = F
  ind2 = which(df.out$h1_reporter_seq == df.out$h2_reporter_seq)
  df.out$is_reporter_same[ind2] = T
  
  return(list(df.out, ls.detect))
}

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

get_revcom_seq <- function(seq){
  if (is.na(seq)){
    seq_rc = NA
  } else {
    base_complement = list("A" = "T", "T" = "A", "C" = "G", "G" = "C")
    seq = unlist(stringr::str_split(seq, pattern = ""))
    seq_rc = paste0(rev(unlist(base_complement[seq])), collapse = "")
  }
  return(seq_rc)
}