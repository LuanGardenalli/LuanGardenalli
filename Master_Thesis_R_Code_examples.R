#Pre processing data:
pisaRex.preprocess <- function(data.peptides,
                               data.proteins,
                               cell.line,
                               type,
                               dimension){
  library(tidyverse)
  library(readxl)
  library(limma)
  library(FactoMineR)
  library(factoextra)
  library(missMDA)
  library(EnhancedVolcano)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(patchwork)
  
  ##################################################
  if (dimension == "REX") {
    data <- data.peptides
    data <- data %>% 
      mutate(peptideRowID = row_number()) %>% 
      filter(Contaminant=="FALSE") %>% 
      dplyr::select(peptideRowID,`Annotated Sequence` : `Modifications in Master Proteins`, 
                    contains("Abundance:") & 
                      contains(cell.line) & 
                      contains(type) & 
                      contains(dimension)) %>% 
      drop_na(starts_with("Abundance:")) %>% 
      mutate(across(starts_with("Abundance:"), ~.x / sum(.x)))
    
    #Non-Cys peptide normalization:
    cys.peptides <- data %>% 
      filter(grepl('C', `Annotated Sequence`)) %>% 
      filter(grepl("Carbamidomethyl", Modifications))
    non.cys.peptides <- data %>% 
      filter(!grepl('C', `Annotated Sequence`))
    non.cys.accessions <- non.cys.peptides %>%
      dplyr::select(`Master Protein Accessions`) %>% 
      pull()
    non.cys.id <- non.cys.peptides %>% 
      mutate(non.cys.id = row_number())
    cys.id <- cys.peptides %>% 
      mutate(cys.id = row_number()) %>% 
      inner_join(non.cys.id %>% dplyr::select(non.cys.id, `Master Protein Accessions`), relationship = "many-to-many")
    denominators <- cys.id %>%
      dplyr::select(cys.id, non.cys.id) %>% 
      left_join(non.cys.id, by = "non.cys.id") %>%
      dplyr::select(cys.id , starts_with("Abundance:")) %>% 
      group_by(cys.id) %>% 
      summarise_all(~ sum(.))
    data <- cys.id %>% 
      dplyr::select(!c(non.cys.id, cys.id)) %>% 
      distinct() %>% 
      mutate(across(starts_with("Abundance:"), ~ .x / denominators$.x))
    data$`Master Protein Accessions` <- sub("-.*", "", data$`Master Protein Accessions`)
    data <- data %>% mutate(`Master Protein Accessions` = word(`Master Protein Accessions`, 1, sep = ";"))
    
  } else if (dimension == "PISA") {
    denominators <- data.proteins %>% 
      filter(`Protein FDR Confidence: Combined`=="High" & Contaminant=="FALSE") %>% 
      dplyr::select(`Accession`, `Gene Symbol`,
                    contains("Abundance:") & 
                      contains(cell.line) & 
                      contains(type) & 
                      contains("EXP")) %>% 
      drop_na() %>% 
      mutate(across(starts_with("Abundance:"), ~.x / sum(.x))) %>% 
      mutate(mean.control = rowMeans(pick(3:5), na.rm = TRUE),
             mean.sample = rowMeans(pick(6:8), na.rm = TRUE)) %>% 
      dplyr::select(Accession, `Gene Symbol`, mean.control, mean.sample)
    
    data <- data.proteins %>% 
      filter(`Protein FDR Confidence: Combined`=="High" & Contaminant=="FALSE") %>% 
      dplyr::select(`Accession`, `Gene Symbol`,
                    contains("Abundance:") & 
                      contains(cell.line) & 
                      contains(type) & 
                      contains(dimension)) %>% 
      drop_na() %>%
      mutate(across(starts_with("Abundance:"), ~.x / sum(.x)))
    
    data <- merge(data, denominators, .by = c("Accession", "Gene Symbol")) %>% 
      mutate(across(contains("Control") & contains("Abundance:"), ~ .x / mean.control),
             across(contains("Sample") & contains("Abundance:"), ~ .x / mean.sample)) %>% 
      dplyr::select(!mean.control & !mean.sample)
    
  } else if ( dimension == "EXP") {
    data <- data.proteins
    data <- data %>% 
      filter(`Protein FDR Confidence: Combined`=="High" & Contaminant=="FALSE") %>% 
      dplyr::select(`Accession`, `Gene Symbol`,
                    contains("Abundance:") & 
                      contains(cell.line) & 
                      contains(type) & 
                      contains(dimension)) %>% 
      drop_na() %>% 
      mutate(across(starts_with("Abundance:"), ~.x / sum(.x)))
  }
  
  #CV filter and log2 transformation:
  Control <- names(data %>% dplyr::select(contains("Control")))
  Sample <- names(data %>% dplyr::select(contains("Sample")))
  data <- data %>% 
    mutate(`Control CV (%)` = (apply(.[Control], 1, sd, na.rm=TRUE) / rowMeans(.[Control], na.rm=TRUE)) * 100, 
           `Sample CV (%)` = (apply(.[Sample], 1, sd, na.rm=TRUE) / rowMeans(.[Sample], na.rm=TRUE)) * 100) %>% 
    filter(`Control CV (%)` < 30 & `Sample CV (%)` < 30)
  data[c(Control, Sample)][data[c(Control, Sample)] == 0] <- NA
  data[c(Control, Sample)] <- log2(data[c(Control, Sample)])
  return(data)
} 

#Statistical analysis and QC plots:
pisaRex.pairwiseReport <- function(data.peptides,
                                   data.proteins,
                                   cell.line,
                                   type,
                                   dimension){
  library(tidyverse)
  library(readxl)
  library(limma)
  library(FactoMineR)
  library(factoextra)
  library(missMDA)
  library(EnhancedVolcano)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(patchwork)
  
  data <- pisaRex.preprocess(data.peptides, data.proteins, cell.line, type, dimension)
  
  Control <- names(data %>% 
                     dplyr::select(contains("Control")) %>% 
                     dplyr::select(!`Control CV (%)`))
  
  Sample <- names(data %>% 
                    dplyr::select(contains("Sample")) %>% 
                    dplyr::select(!`Sample CV (%)`)) 
  
  #Differential Expression with Limma
  samples = grep('Abundance:', colnames(data), value=T)
  treatment <- c(rep("Control", length(Control)), rep("Sample", length(Sample)))
  replicates <- c(rep(c(1:length(Sample)), 2))
  sample.table <- tibble("Sample" = samples,
                         "Treatment" = treatment,
                         "Replicates" = replicates)
  sample.table <- column_to_rownames(sample.table, var = "Sample")
  design <- sapply(unique(sample.table$Treatment), function(x) sample.table$Treatment %in% x)
  design = cbind(design)+0
  if (dimension == "REX") {
    contrast = makeContrasts(
      Control - Sample,
      levels = design)
  } else {
    contrast = makeContrasts(
      Sample - Control,
      levels = design)
  }
  fit <- lmFit(data[c(Control, Sample)], design)
  fit = contrasts.fit(fit, contrast)
  fit = eBayes(fit)
  if (dimension == "REX") {
    DEGs <- topTable(fit, 
                     number = Inf, 
                     genelist = data$`Master Protein Accessions`, 
                     adjust.method = "BH", 
                     sort.by = "none", 
                     p.value = 1, 
                     lfc = 0, 
                     confint = F) %>% 
      mutate(`ID` = word(ID, 1, sep = ";"))
  } else {
    DEGs <- topTable(fit, 
                     number = Inf, 
                     genelist = data$`Gene Symbol`, 
                     adjust.method = "BH", 
                     sort.by = "none", 
                     p.value = 1, 
                     lfc = 0, 
                     confint = F) %>% 
      mutate(`ID` = word(ID, 1, sep = ";"))
  }
  DEGs <- DEGs %>% 
    mutate(FC.rank = rank(-abs(logFC))) %>% 
    mutate(P.value.rank = rank(adj.P.Val)) %>% 
    mutate(Rank_sum = P.value.rank + FC.rank) %>% 
    mutate(Rank = rank(Rank_sum)) %>% 
    dplyr::select(!FC.rank & !P.value.rank & !Rank_sum)
  
  if (dimension == "REX") {
    DEGs <- DEGs %>% mutate(peptideRowID = data$peptideRowID)
    DEGs$ID <- sub("-.*", "", DEGs$ID)
    gene.map <- as_tibble(bitr(DEGs$ID, 
                               fromType = "UNIPROT",  
                               toType = c("SYMBOL"),    
                               OrgDb = org.Hs.eg.db))
    colnames(gene.map) <- c("ID", "Symbol")
    DEGs <- merge(DEGs, gene.map, by.x = "ID", by.y = "ID", all.x = TRUE)
    
    rowID <- data %>% 
      dplyr::select(peptideRowID, `Positions in Master Proteins`, `Annotated Sequence`)
    
    DEGs <- merge(DEGs, rowID, by.x = "peptideRowID", by.y = "peptideRowID", all.x = F)
  } else {
    DEGs <- DEGs %>% 
      dplyr::rename(Symbol = ID) %>% 
      mutate(proteinRowID = row_number())
  }
    return(DEGs %>% 
             drop_na())
}


#Cysteine Heatmap:
{
  log2ratios <- function(data.peptides,
                         data.proteins,
                         cell.line,
                         type,
                         proteins) {
    limma  <- pisaRex.pairwiseReport.accession(data.peptides,
                                               data.proteins,
                                               cell.line,
                                               type,
                                               dimension="REX", 
                                               output = "genes") %>% 
      filter(Accession %in% proteins)
    
    replicates <- pisaRex.preprocess(data.peptides,
                                     data.proteins,
                                     cell.line,
                                     type,
                                     dimension="REX") %>% 
      filter(`Master Protein Accessions` %in% proteins) %>% 
      add_column(FDR = limma$adj.P.Val, P = limma$P.Value) %>% 
      mutate(start_pos = as.numeric(str_extract(`Positions in Master Proteins`, "(?<=\\[)\\d+(?=-)")),
             clean_seq = str_extract(`Annotated Sequence`, "(?<=\\.)[A-Z]+(?=\\.)"),
             `Cys Position` = map2_chr(clean_seq, start_pos, ~ {
               
               c_rel_pos <- str_locate_all(.x, "C")[[1]][, "start"]
               if (length(c_rel_pos) == 0) {
                 return(NA_character_)
               }
               
               abs_pos <- .y + c_rel_pos - 1
               paste0("Cys", paste(abs_pos, collapse = "/"))})) %>%
      select(-start_pos, -clean_seq) %>%
      mutate(`Cys Position` = paste0(`Master Protein Accessions`,": ", `Cys Position`))
    
    log2ratios <- replicates %>% 
      dplyr::select(`Master Protein Accessions`, `Cys Position`, FDR, P, contains("Abundance:")) %>% 
      arrange(FDR) %>% 
      distinct(`Cys Position`, .keep_all = T) %>% 
      rowwise() %>% 
      mutate(control_mean = rowMeans(pick(contains("Control")), na.rm = TRUE)) %>% 
      mutate(across(contains("Sample"), ~ control_mean - .x)) %>% 
      dplyr::select(!control_mean & !contains("Control")) %>% 
      mutate(`Cys Position` = str_remove(`Cys Position`, "^.*?: "))
    
    return(log2ratios)
  }
  
  #Input:
  {
    cell.line <- "THP1"
    TA <- log2ratios(data.peptides,
                     data.proteins,
                     cell.line,
                     type="Alpha",
                     proteins)
    
    TB <- log2ratios(data.peptides,
                     data.proteins,
                     cell.line,
                     type="Beta",
                     proteins)
    
    TG <- log2ratios(data.peptides,
                     data.proteins,
                     cell.line,
                     type="Gamma",
                     proteins)
    
    cell.line <- "HL60"
    HA <- log2ratios(data.peptides,
                     data.proteins,
                     cell.line,
                     type="Alpha",
                     proteins)
    
    HB <- log2ratios(data.peptides,
                     data.proteins,
                     cell.line,
                     type="Beta",
                     proteins)
    
    HG <- log2ratios(data.peptides,
                     data.proteins,
                     cell.line,
                     type="Gamma",
                     proteins)
    
    heatmap.df.T <- full_join(TA, 
                              TB, 
                              by = join_by(`Master Protein Accessions`, `Cys Position`), 
                              keep = F,
                              suffix = c(".alpha", ".beta")) %>% 
      full_join(TG, 
                by = join_by(`Master Protein Accessions`, `Cys Position`), 
                keep = F) %>% 
      rename(FDR.gamma = FDR,
             P.gamma = P)
    
    heatmap.df.H <- full_join(HA, 
                              HB, 
                              by = join_by(`Master Protein Accessions`, `Cys Position`), 
                              keep = F,
                              suffix = c(".alpha", ".beta")) %>% 
      full_join(HG, 
                by = join_by(`Master Protein Accessions`, `Cys Position`), 
                keep = F) %>% 
      rename(FDR.gamma = FDR,
             P.gamma = P)
    
    heatmap.df <- full_join(heatmap.df.T,
                            heatmap.df.H,
                            by = join_by(`Master Protein Accessions`, `Cys Position`),
                            keep = F,
                            suffix = c(".THP1", ".HL"))
    heatmap.df[9,2] <- " Cys189"
    
    FDRs <- heatmap.df %>% 
      select(contains("FDR"))
    heatmap.df <- heatmap.df %>% select(!contains("FDR"))
    
    Ps <- heatmap.df %>% 
      select(contains("P."))
    heatmap.df <- heatmap.df %>% select(!contains("P."))
    
    replicate.long <- heatmap.df %>%
      select(!`Master Protein Accessions`) %>% 
      pivot_longer(-`Cys Position`,names_to = "sample", values_to = "log2fc") %>% 
      mutate(group = str_remove(sample, "^.*?, Sample, "))
    
    colnames(FDRs) <- unique(replicate.long$group)
    FDRs$`Cys Position` <- heatmap.df$`Cys Position`
    
    colnames(Ps) <- unique(replicate.long$group)
    Ps$`Cys Position` <- heatmap.df$`Cys Position`
    
    fdr.long <- FDRs %>% 
      pivot_longer(-`Cys Position`,names_to = "group", values_to = "fdr")
    
    p.long <- Ps %>% 
      pivot_longer(-`Cys Position`,names_to = "group", values_to = "p")
    
    fdr.matrix <- replicate.long %>%
      left_join(fdr.long, by = c("Cys Position", "group")) %>%
      select(`Cys Position`, sample, fdr) %>%
      pivot_wider(names_from = sample, values_from = fdr) %>%
      column_to_rownames("Cys Position") %>%
      as.matrix()
    
    p.matrix <- replicate.long %>%
      left_join(p.long, by = c("Cys Position", "group")) %>%
      select(`Cys Position`, sample, p) %>%
      pivot_wider(names_from = sample, values_from = p) %>%
      column_to_rownames("Cys Position") %>%
      as.matrix()
    
    
    heatmap.df[is.na(heatmap.df)] <- 0
    
    heatmap.matrix <- heatmap.df %>%
      select(!`Master Protein Accessions`) %>% 
      column_to_rownames("Cys Position") %>%
      as.matrix()
    
    gene.map <- as_tibble(bitr(heatmap.df$`Master Protein Accessions`, 
                               fromType = "UNIPROT",  
                               toType = c("SYMBOL"),    
                               OrgDb = org.Hs.eg.db))
    heatmap.df <- left_join(heatmap.df, gene.map, join_by(`Master Protein Accessions` == UNIPROT)) %>% 
      distinct(`Cys Position`, .keep_all = T)
  }
  
  
  #HM:
  {
    names <- c("", "THP1, Alpha", "",
               "", "THP1, Beta",  "",
               "", "THP1, Gamma", "",
               
               "", "HL60, Alpha", "",
               "", "HL60, Beta",  "",
               "", "HL60, Gamma", "")
    
    my_colors <- colorRamp2(c(-1, 0, 2), 
                            c("#2166AC", "#F7F7F7","#B2182B"))
    
    rows <- c(4,6,43,33,3,7,23,29,10,12,39,38,5,1,36,31,28,8,9,41,32,40,21,11,19,30,18,46,44)
    highlight <- heatmap.df %>% 
      select(SYMBOL, `Cys Position`) %>% 
      rowid_to_column("row") %>% 
      mutate(color = if_else(row %in% rows, "red", "black"))
    
    
    
    ha = rowAnnotation(foo = anno_empty(border = FALSE, 
                                        width = max_text_width(proteins)))
    HM <- Heatmap(heatmap.matrix,
                  name = "Log2\nRatio",
                  col = my_colors,
                  cluster_rows = T,
                  cluster_columns = F,
                  
                  
                  row_split = heatmap.df$SYMBOL,       
                  row_title = "Oxidized Cysteine Positions",              
                  show_row_names = T,
                  row_names_side = "left",
                  row_names_gp = gpar(fontsize = 8, col = highlight$color),
                  right_annotation = ha,
                  show_row_dend = F,
                  
                  show_column_names = TRUE,
                  column_labels = names,
                  column_names_rot = 35,
                  column_names_gp = gpar(fontsize = 8),
                  show_heatmap_legend = F,
                  
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    if (!is.na(fdr.matrix[i, j]) && fdr.matrix[i, j] < 0.05 && fdr.matrix[i, j] > 0.001) {
                      grid.text("*", x, y, gp = gpar(fontsize = 10, col = "black"))
                    } else if (!is.na(fdr.matrix[i, j]) && fdr.matrix[i, j] < 0.001) {
                      grid.text("**", x, y, gp = gpar(fontsize = 10, col = "black"))
                    } else if (!is.na(p.matrix[i, j]) && p.matrix[i, j] < 0.05) {
                      grid.text("#", x, y, gp = gpar(fontsize = 7, col = "black"))
                    }
                  })
    
    
    
    draw(HM)
    names <- names(row_order(HM))
    for(i in 1:9) {
      decorate_annotation("foo", slice = i, {
        grid.rect(x = 0, width = unit(1, "mm"), gp = gpar(fill = "black", col = NA), just = "right")
        grid.text(paste(names[[i]], collapse = "\n"), x = unit(1, "mm"), just = "left",
                  gp = gpar(
                    fontsize = 8, 
                    fontfamily = "sans", 
                    fontface = "bold"))})
    }
  }
}




#DIABLO Analysis:
{
  format_for_diablo <- function(df) {
    df %>%
      pivot_longer(cols = -Accession, names_to = "Sample", values_to = "Intensity") %>%
      pivot_wider(names_from = Accession, values_from = Intensity) %>%
      column_to_rownames("Sample") %>%
      as.matrix()
}
  
  diablo.input <- function(data.peptides,
                           data.proteins,
                           cell.line,
                           type) {
    
    names <- c("Accession", 
               paste0("Control1_",type), 
               paste0("Control2_",type), 
               paste0("Sample1_",type), 
               paste0("Sample2_",type))
    names_rex <- c("Accession", 
                   paste0("Sample1_",type), 
                   paste0("Sample2_",type), 
                   paste0("Control1_",type), 
                   paste0("Control2_",type))
    
    exp <- pisaRex.preprocess(data.peptides,
                              data.proteins,
                              cell.line,
                              type,
                              dimension="EXP")
    
    
    exp <- exp %>% dplyr::select(Accession, starts_with("Abundance"))
    exp <- exp %>% select(1:3, 5:6)
    colnames(exp) <- names
    
    pisa <- pisaRex.preprocess(data.peptides,
                               data.proteins,
                               cell.line,
                               type,
                               dimension="PISA")
    pisa <- pisa %>% dplyr::select(Accession, starts_with("Abundance"))
    colnames(pisa) <- names
    
    
    rex <- pisaRex.preprocess(data.peptides,
                              data.proteins,
                              cell.line,
                              type,
                              dimension="REX")
    
    rex <- rex %>% dplyr::select(peptideRowID, starts_with("Abundance"))
    rex <- rex %>% select(1:3, 5:6)
    colnames(rex) <- names_rex
    rex <- rex %>% select(Accession,
                          paste0("Control1_",type), 
                          paste0("Control2_",type), 
                          paste0("Sample1_",type), 
                          paste0("Sample2_",type))
    return(list("exp"=exp, "pisa"=pisa, "rex"=rex))
}
  
  
  #Run DIABLO on all IFN types:
  cell.line <-  "THP1"
  {
    #Input:
    type <- "Alpha"
    alpha <- diablo.input(data.peptides,
                          data.proteins,
                          cell.line,
                          type)
    type <- "Beta"
    beta <- diablo.input(data.peptides,
                         data.proteins,
                         cell.line,
                         type)
    type <- "Gamma"
    gamma <- diablo.input(data.peptides,
                          data.proteins,
                          cell.line,
                          type)
    
    X_exp_alpha <- format_for_diablo(alpha$exp)
    X_pisa_alpha <- format_for_diablo(alpha$pisa)
    X_rex_alpha <- format_for_diablo(alpha$rex)
    
    X_exp_beta <- format_for_diablo(beta$exp)
    X_pisa_beta <- format_for_diablo(beta$pisa)
    X_rex_beta <- format_for_diablo(beta$rex)
    
    X_exp_gamma <- format_for_diablo(gamma$exp)
    X_pisa_gamma <- format_for_diablo(gamma$pisa)
    X_rex_gamma <- format_for_diablo(gamma$rex)
    
    
    common_exp_features <- Reduce(intersect, list(
      colnames(X_exp_alpha), 
      colnames(X_exp_beta), 
      colnames(X_exp_gamma)
    ))
    X_exp_alpha_aligned <- X_exp_alpha[, common_exp_features]
    X_exp_beta_aligned  <- X_exp_beta[, common_exp_features]
    X_exp_gamma_aligned <- X_exp_gamma[, common_exp_features]
    
    X_exp_combined <- rbind(X_exp_alpha_aligned, X_exp_beta_aligned, X_exp_gamma_aligned)
    
    
    common_pisa_features <- Reduce(intersect, list(
      colnames(X_pisa_alpha), 
      colnames(X_pisa_beta), 
      colnames(X_pisa_gamma)
    ))
    X_pisa_alpha_aligned <- X_pisa_alpha[, common_pisa_features]
    X_pisa_beta_aligned  <- X_pisa_beta[, common_pisa_features]
    X_pisa_gamma_aligned <- X_pisa_gamma[, common_pisa_features]
    
    X_pisa_combined <- rbind(X_pisa_alpha_aligned, X_pisa_beta_aligned, X_pisa_gamma_aligned)
    
    
    common_rex_features <- Reduce(intersect, list(
      colnames(X_rex_alpha), 
      colnames(X_rex_beta), 
      colnames(X_rex_gamma)
    ))
    X_rex_alpha_aligned <- X_rex_alpha[, common_rex_features]
    X_rex_beta_aligned  <- X_rex_beta[, common_rex_features]
    X_rex_gamma_aligned <- X_rex_gamma[, common_rex_features]
    
    X_rex_combined <- rbind(X_rex_alpha_aligned, X_rex_beta_aligned, X_rex_gamma_aligned)
    
    X_list <- list(
      expression = X_exp_combined,
      solubility = X_pisa_combined,
      redox = X_rex_combined)
    
    Y_condition <- factor(c(
      "Control", "Control", "Alpha", "Alpha", 
      "Control", "Control", "Beta", "Beta",    
      "Control", "Control", "Gamma", "Gamma"  
    ))
    
    paired_design <- data.frame(
      sample_pair = factor(c(
        "Pair_A1", "Pair_A2", "Pair_A1", "Pair_A2",  
        "Pair_B1", "Pair_B2", "Pair_B1", "Pair_B2",  
        "Pair_C1", "Pair_C2", "Pair_C1", "Pair_C2" )))
    
    X_list_within <- list(
      expression = withinVariation(X = X_list$expression, design = paired_design),
      solubility = withinVariation(X = X_list$solubility, design = paired_design),
      redox      = withinVariation(X = X_list$redox,      design = paired_design)
    )
    
    design_matrix <- matrix(0.1, ncol = length(X_list), nrow = length(X_list), 
                            dimnames = list(names(X_list), names(X_list)))
    diag(design_matrix) <- 0
    
    keepX_grid <- c(10, 20, 30, 50, 75)
    
    tune_parameters <- list(
      solubility = keepX_grid,
      redox = keepX_grid,
      expression = keepX_grid)
    
    library(BiocParallel)
    parallel_setup <- SnowParam(workers = 6)
    
    diablo_tuning <- tune.block.splsda(
      X = X_list_within, 
      Y = Y_condition, 
      ncomp = 3, 
      test.keepX = tune_parameters, 
      design = design_matrix, 
      validation = 'loo', 
      dist = "centroids.dist", 
      BPPARAM = parallel_setup,
      progressBar = TRUE)
    
    optimal_keepX <- diablo_tuning$choice.keepX
    plot(diablo_tuning)
    
    final_multiclass_diablo <- block.splsda(
      X = X_list_within, 
      Y = Y_condition, 
      ncomp = 3, 
      keepX = optimal_keepX, 
      design = design_matrix)
    
}
}