#!/usr/bin/env Rscript
rm(list = ls())
if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse", repos = "https://cloud.r-project.org/")
}
library(optparse)
# ---- Define command line options ----
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Input MAGeCK.count.txt",default = "E:/2025_UdeM_phd/Haley_lab/CRISPR-Cas9/data/251010_SP1_Jurkat_mageck.count.txt"),
  make_option(c("-o", "--output"), type = "character", help = "Create all output files in the specified output directory",default = paste0(Sys.Date(),"_HaleyLab_CRISPRa")),
  make_option(c("-p", "--output_prefix"), type = "character", help = "The prefix of the output file", default = as.character(Sys.Date())),
  make_option(c("-l", "--sample_label"), type = "character", help = "The SP and cell line of samples, with format SP_cellline, e.g., SP1_K562", default = "SP_cell"),
  make_option(c("-c", "--contrast"), type = "character", help = "Contrast array between samples at different time points, separated by comma, e.g., T17-T0,T10-T0,T24_DAR-T19_DMSO", default = "T17-T0,T10-T0")
              
  
)

# ---- Parse arguments ----
opt_parser <- OptionParser(option_list = option_list,
                           usage = "usage: Rscript CRISPRa_workflow_T0med3mad_normbyNTC.R \n
                           -i SP1_K562_mageck.count.txt\n
                           -o SP1_K562_v1\n
                           -p SP1_K562\n
                           -l SP1_K562\n
                           -c T17-T0,T10-T0",
                           description = "\nHi sweet, this workflow was developed for Genome-wide CRISPR activation screening data based on R-4.4.3.
                           \nPlease let Xiaozhen know if you have any questions!
                           \nImportant notes:
                           \n1. The input should be raw counts produced by MAGeCK count with columns separated by tab;
                           \n2. Two replicates for each time point are required;
                           \n3. The column names for input should be: sgRNA, Gene, SPn_Cellline_Tn_repn. For example: SP1_K562_T0_rep1, or SP1_H23_T24_DAR_rep1;
                           \n4. Each paramater should be specified.\n
                           ")
# print_help(opt_parser)
opt <- parse_args(opt_parser)
outDir <- paste0("voom_output/",opt$output,'/')
prefix <- paste0(opt$output_prefix,"_")
sampletitle<-opt$sample_label
contrast_array <- unlist(strsplit(opt$contrast, ","))
if (file.exists("voom_output")){
  cat("filefold voom_output exists\n")
} else {
  dir.create("voom_output")
  cat("creat filefold voom_output\n")
}
if (file.exists(outDir)){
  cat(paste0("filefold ",outDir," exists\n"))
} else {
  dir.create(outDir)
  cat(paste0("creat filefold ",outDir,"\n"))
}

#load the required libraries
pkg_check <- function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE,quietly = TRUE))
}

# List of required packages
packages1 <- c(
  "data.table",
  "ggplot2",
  "stringr",
  "dplyr",
  "cowplot",
  "ggpubr",
  "reshape2",
  "ggpointdensity",
  "viridis"
)

# Install and load package via BiocManager
pkg_check_BiocManager <- function(pkg) {
  # Check and install BiocManager if not available
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org/")
  }
  
  # Check and install target package if not available
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
  
  # Load package quietly
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE, quietly = TRUE)
  )
}

# List of required packages
packages2 <- c(
  "edgeR",
  "limma",
  "SummarizedExperiment",
  "EnhancedVolcano"
  )

# Check and load all
lapply(packages1, pkg_check)
lapply(packages2, pkg_check_BiocManager)
cat("----------------------------------------------------------------------------------------\n")
cat("Start analysis...\n")
suppressMessages(suppressWarnings({
  # import sgRNA raw counts output from MAGeCK count
  dat_dropout_sample_all<-read.csv(opt$input,sep = '\t')
  rownames(dat_dropout_sample_all)<-dat_dropout_sample_all$sgRNA
  dat_dropout_sample<-dat_dropout_sample_all[,grep('Gene|sgRNA',colnames(dat_dropout_sample_all),invert = T)]
  
  # Create sample metadata (must match number of columns in `counts`)
  sample_info <-   data.frame(
    row.names = colnames(dat_dropout_sample),
    Time = factor(gsub(".*_(T\\d+.*)_(rep.*)", "\\1", colnames(dat_dropout_sample))),
    Replicate = gsub(".*_(T\\d+.*)_(rep.*)", "\\2", colnames(dat_dropout_sample))
  )
  print(sample_info)
  cat("\n")
  
  # step 1 filtering
  se <- SummarizedExperiment(dat_dropout_sample, colData=sample_info)
  is.ref <- colData(se)[["Time"]] %in% "T0"
  ref.ab <- aveLogCPM(assay(se, 1)[,is.ref])
  med.ab <- median(ref.ab)
  mad.ab <- mad(ref.ab, center=med.ab)
  threshold <- med.ab - 3*mad.ab
  filtered <- ref.ab >= threshold
  se <- se[filtered,]
  
  # hist plot
  p<-ggplot(data.frame(Ref.aveLogCPM=ref.ab),aes(Ref.aveLogCPM))+
    geom_histogram(binwidth = 0.5, alpha = 0.7, fill = "grey", color = "black") +
    theme_cowplot()+
    ggtitle("Distribution of Ref.aveLogCPM",paste0(sampletitle," n=",length(ref.ab)))+
    geom_vline(xintercept = c(threshold,med.ab), color = c("red","blue"), linetype = "dashed", size = 1) +
    scale_x_continuous(breaks = seq(-4,4,2))+
    theme(strip.text = element_text(size = 10),
          strip.background = element_blank(),
          plot.title = element_text(size=12))
    
  ggsave(paste0(outDir,prefix,"Ref.aveLogCPM_dist.pdf"),
         width = 4,height = 3, limitsize= F, p)
  
  
  # Step 2: Normalization with NTCs
  use.for.norm <- grepl("NTC_", rownames(se))
  Y <- log2(assays(se)[[1]]+1)
  factors <- apply(Y[use.for.norm,], 2, median)
  factors <- factors-median(factors)
  norm_factors<-(2^factors)*10^6-1
  print(norm_factors)
  cat("\n")
  
  log_counts<-reshape2::melt(Y,
                             value.name = "log_counts",
                             measure.vars = colnames(Y),
                             variable.name = "Time")
  dat<-log_counts
  p1<-ggplot(dat,aes(log2(log_counts+1),color=Time))+
    geom_density(linewidth = 0.3) +
    theme_cowplot()+
    scale_color_manual(values = rainbow(ncol(Y)))+
    ggtitle("Distribution of log_counts",paste0(sampletitle," n=",nrow(Y)))+
     scale_x_continuous(breaks = seq(0,12,4))+
    theme(strip.text = element_text(size = 10),
          strip.background = element_blank(),
          plot.title = element_text(size=12))
  
  
  # Step 3: differential abundance analysis
  # Creating design matrix for differential analysis
  design <- model.matrix(~ 0 + Time, data = sample_info)
  colnames(design) <- levels(sample_info$Time)
  
  # Make pairwise contrasts
  # Create contrast matrix using the string
  contrast_matrix <- makeContrasts(contrasts = contrast_array, levels = design)
  colnames(contrast_matrix)<-gsub("-","_vs_",colnames(contrast_matrix))
  print(contrast_matrix)
  cat("\n")
  
  # voom transformation with weights
  pdf(paste0(outDir,prefix,"sgRNA_differential_abundance_limma_mean-variance.pdf"),width = 4*2, height = 4)
  if(length(levels(sample_info$Time))>=3){
    v <- voomWithQualityWeights(assays(se)[[1]],
                              design=design,
                              plot=TRUE,
                              block=colData(se)$Replicate,
                              lib.size = norm_factors)
  }else{
    v <- voomWithQualityWeights(assays(se)[[1]],
                                design=design,
                                plot=TRUE,
                                lib.size = norm_factors)
  }
  dev.off()
  # Access normalized logCPM values
  logCPM <- as.data.frame(v$E)
  dat<-reshape2::melt(logCPM,
                      value.name = "log_norm_counts",
                      measure.vars = colnames(logCPM),
                      variable.name = "Time")
  p3<-ggplot(dat,aes(log_norm_counts,color=Time))+
    geom_density(linewidth = 0.3) +
    theme_cowplot()+
    scale_color_manual(values = rainbow(ncol(logCPM)))+
    ggtitle("Distribution of log_norm_counts",paste0(sampletitle," n=",nrow(logCPM)))+
    # geom_vline(xintercept = c(threshold,med.ab), color = c("red","blue"), linetype = "dashed", size = 1) +
    scale_x_continuous(breaks = seq(0,12,4))+
    theme(strip.text = element_text(size = 10),
          strip.background = element_blank(),
          plot.title = element_text(size=12))
  
  p<-plot_grid(p1,p3,nrow = 1)
  ggsave(paste0(outDir,prefix,"log_norm_counts_dist.pdf"),
         width = 4*2,height = 3,limitsize= F, p)
  
  
  # Fit linear model
  fit <- lmFit(v, design)
  fit <- eBayes(fit, robust=TRUE)
  # Create gene sets: group sgRNAs by promoter name
  dat_dropout_sample_all$promoter<-gsub("(.*)_(.*)_(.*)$","\\1_\\2",dat_dropout_sample_all$sgRNA)
  df_sg<-dat_dropout_sample_all[grep('^NTC$',dat_dropout_sample_all$Gene,invert = T),c("sgRNA","promoter")]
  df_sg<-df_sg %>%
    dplyr::filter(sgRNA %in% rownames(Y))
  promoter_sets <- split(df_sg$sgRNA, df_sg$promoter)
  
  
  #correlation between replicates
  library(ggpubr)
  dat<-as.data.frame(logCPM)
  CPM_corr<-list()
  con_list<-unique(gsub("_rep\\d","",rownames(sample_info)))
  for (i in con_list) {
    con_sel<-i
    CPM_corr[[i]]<-ggplot(dat,aes_string(x=paste0(con_sel,"_rep1"),y=paste0(con_sel,"_rep2")))+
      geom_point(size=0.001)+
      geom_pointdensity(size=0.001) +
      scale_color_viridis("Density")+
      ggtitle("log2(normalized counts) corr",paste0("\n",sampletitle," ",con_sel," n=",nrow(dat)))+
      theme_cowplot()+
      stat_cor(method = "pearson",label.x = -1, label.y = 13)+
      geom_abline(slope = 1,intercept = 0,alpha=.5,linetype="dashed")+
      coord_cartesian(ylim = c(-2,15),xlim = c(-2,15))
  }
  CPM_corrall<-plot_grid(plotlist = CPM_corr,nrow = 1)
  ggsave(paste0(outDir,prefix,"sgRNA_differential_abundance_limma_log2normedC_corr.pdf"),
         width = 4.5*length(levels(sample_info$Time)),height = 4,limitsize= F, CPM_corrall)
  
  # contrast
  my_compa<-colnames(contrast_matrix)
  my_compa_result<-list()
  my_fry_results<-list()
  my_fry_results_all<-list()
  my_compa_result_plot<-list()
  my_fry_results_plot<-list()
  for (k in 1:length(my_compa)) {
    # Make pairwise contrasts
    contrast_matrix_sel<-as.matrix(contrast_matrix[,k])
    colnames(contrast_matrix_sel)<-my_compa[k]
    
    # Apply contrasts
    fit2 <- contrasts.fit(fit, contrast_matrix_sel)
    
    # Apply empirical Bayes
    fit2 <- eBayes(fit2, robust = TRUE)
    my_compa_result[[k]] <- topTableF(fit2, number = Inf, adjust.method = "BH")
    write.table(cbind(sgRNA=rownames(my_compa_result[[k]])
                      ,my_compa_result[[k]]),paste0(outDir, prefix,my_compa[k],"_sgRNA_differential_abundance_limma.txt"),sep = '\t',row.names = F)
    
    # Volcano plot
    dat<-my_compa_result[[k]] %>%
      mutate(log10_pvalue=-log10(P.Value),
             type=case_when(grepl('^NTC_',rownames(my_compa_result[[k]])) ~ "NTC",
                            .default = "Target")) %>%
      mutate(type=factor(type, levels = c("Target","NTC"),ordered = T)) %>%
      arrange(type)
    
    my_compa_result_plot[[k]]<-ggplot(dat, aes_string(x = my_compa[k], y = "log10_pvalue")) +
      # Add significance threshold lines
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", alpha = 0.8) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50", alpha = 0.8) +
      # Add points colored by guide type
      geom_point(aes(color = type), alpha = 0.7, size = 0.5) +
      # Custom colors: black for NTC, grey for others
      scale_color_manual(
        name = "Guide Type",
        values = c("NTC" = "black", "Target" = "grey50"),
        guide = guide_legend(override.aes = list(size = 3, alpha = 1))
      ) +
      
      # Axis labels and limits
      labs(
        x = paste0("Log2 Fold Change (", my_compa[k], ")"),
        y = expression(-Log[10]~italic(P)),
        title = paste0(sampletitle, " ", my_compa[k])
      ) +
      annotate("text", x = 3, y = max(-log10(my_compa_result[[k]]$P.Value)+0.5, na.rm = TRUE), 
               label = paste0("up: ", nrow(my_compa_result[[k]][which(my_compa_result[[k]][,my_compa[k]]>=1 & my_compa_result[[k]]$P.Value < 0.05),])), 
               hjust = 0) +
      annotate("text", x = -5, y = max(-log10(my_compa_result[[k]]$P.Value)+0.5, na.rm = TRUE), 
               label = paste0("down: ", nrow(my_compa_result[[k]][which(my_compa_result[[k]][,my_compa[k]]<=(-1) & my_compa_result[[k]]$P.Value < 0.05),])), 
               hjust = 0)+
      # Set axis limits
      coord_cartesian(xlim = c(-10, 10)) +
        # Styling
      theme_cowplot() +
      theme(
        plot.title = element_text(size = 12),
        legend.position = "top",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10)
      )
    
    #compute the median of LFC based on all of guides
    my_compa_result[[k]]$promoter<-gsub("(.*)_(.*)_(.*)$","\\1_\\2",rownames(my_compa_result[[k]]))
    my_compa_result[[k]]$Gene<-gsub("(.*)_(.*)_(.*)$","\\1",rownames(my_compa_result[[k]]))
    my_compa_result_seq_median<-my_compa_result[[k]] %>%
      group_by(Gene,promoter) %>%
      summarise(LFC_median=median(.data[[my_compa[k]]], na.rm = TRUE),.groups = 'drop')
  
    # Perform gene-set testing using fry
    my_fry_results[[k]] <- fry(
      y = v,
      index = promoter_sets,
      design = design,
      contrast = contrast_matrix[ , my_compa[k], drop = FALSE],
      adjust.method="BH"
    )
    my_fry_results_all[[k]]<-merge(my_compa_result_seq_median,my_fry_results[[k]],by.x="promoter",by.y="row.names")
    write.table(my_fry_results_all[[k]],paste0(outDir, prefix,my_compa[k],"_promoter_level_fry_results.txt"),sep = '\t',row.names = F)
    
    dat<-my_fry_results_all[[k]]
    my_fry_results_plot[[k]]<-EnhancedVolcano(dat , x="LFC_median", y="PValue" , lab = dat$promoter, subtitle = NULL, 
                                              title =paste0(sampletitle," ",my_compa[k]) , 
                                              ylab = bquote(~-Log[10]~italic(P)) ,
                                              FCcutoff = log2(1.5) , pCutoff = 0.05,labSize=1
    ) + 
      annotate("text", x = 2, y = max(-log10(dat$PValue)+0.5, na.rm = TRUE), 
               label = paste0("up: ", nrow(dat[which(dat$LFC_median>=log2(1.5) & dat$PValue < 0.05),])), 
               hjust = 0) +
      annotate("text", x = -2, y = max(-log10(dat$PValue)+0.5, na.rm = TRUE), 
               label = paste0("down: ", nrow(dat[which(dat$LFC_median<=(-log2(1.5)) & dat$PValue < 0.05),])), 
               hjust = 0)+
      coord_cartesian(xlim = c(-6,6))+
      theme(plot.title = element_text(size=12))+
      theme_cowplot() +
      theme(
        plot.title = element_text(size = 12),
        legend.position = "top",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10)
      )
    
  }
  my_compa_result_plotall<-plot_grid(plotlist = my_compa_result_plot,nrow = 1)
  ggsave(paste0(outDir,prefix,"sgRNA_differential_abundance_limma_volcano.pdf"),
         width = 7*length(my_compa),height = 5.6,limitsize= F, my_compa_result_plotall)
  
  my_fry_results_plotall<-plot_grid(plotlist = my_fry_results_plot,nrow = 1)
  ggsave(paste0(outDir,prefix,"promoter_differential_abundance_limma_volcano.pdf"),
         width = 6*length(my_compa),height = 5.6,limitsize= F, my_fry_results_plotall)
  
  cat(rep("🐱", nrow(sample_info)), "\n")
  cat(paste0("The analysis results have been saved in ",outDir), "\n")
  cat("Meow! 😺\n")
  cat("Well done! Your are excellent!\n")

}))

