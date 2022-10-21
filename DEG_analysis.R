library(DESeq2)
library(tidyverse)
library(ggplot2)
library(caret)
library(dplyr)    
library(pheatmap)
library(RColorBrewer)

      counts_data <- read.csv("feature_counts_all.txt", header = TRUE, row.names = 1, sep = ";")
      
      
      
      # read in sample info
      condition <- c( "Case", "Case", "Control", "Control", "Control")
      
      
      # Create METAdata frame
      colData <- data.frame(condition)
      
      rownames(colData) <- colnames(counts)
      
      # making sure the row names in colData matches to column names in counts_data
      all(colnames(counts_data) %in% rownames(colData))
      
      # are they in the same order?
      all(colnames(counts_data) == rownames(colData))
      
      # Step 2: construct a DESeqDataSet object ----------
      dds <- DESeqDataSetFromMatrix(countData = counts_data,
                                    colData = colData,
                                    design = ~ condition)
      
      
      # pre-filtering: removing rows with low gene counts
      # keeping rows that have at least 10 reads total
      
      keep <- rowSums(counts(dds)) >= 10
      dds <- dds[keep,]
      
      
      # set the factor level
      dds$condition <- relevel(dds$condition, ref = "Control")
      
      
      
      # Step 3: Run DESeq ----------------------
      dds <- DESeq(dds)
      res0.05 <- results(dds, alpha = 0.05)
      
      #Make dataframe
      data <- data.frame(res0.05)
      
      # Set a boolean column for significance
      data$significant <- ifelse(data$padj < 0.05 & data$log2FoldChange > 0 , "Significant-up", 
                                 ifelse(data$padj < 0.05 & data$log2FoldChange < 0,"Significant-down", "Non-sig"))
      
      
      # Create MA plot and main has the groups---------------------
      
      ma <- ggplot(data, aes(baseMean, log2FoldChange, colour=significant, title= "g" )) +
        geom_point(size=1) + 
        scale_x_log10() +
        #geom_text(label=data$label) +
        ggtitle(id) +
        geom_hline(yintercept = 0, colour="tomato1", size=0) +
        labs(x="mean of normalized counts", y="log fold change") +
        scale_colour_manual(name="q-value", values=c("Significant-up"="blue","Significant-down"="red", "Non-sig"="grey50")) +
        theme_bw()
      
      ggsave(paste (id,"MA.jpeg", sep = "_"), plot = ma)
      
      # Create vOLCANO plot and main has the groups--------------------
      
      
      vo <- ggplot(data, aes(y=-log10(padj), x= log2FoldChange, colour=significant)) +
        geom_point(size=1) + 
        ggtitle(id) +
        #geom_text(label=data$label, nudge_x = 2, nudge_y = 2, size= 2) +
        geom_hline(yintercept = 1.3, colour="grey", size=0, linetype = "dashed") +
        geom_vline(xintercept = 1.1, colour="grey", size=1, linetype = "dashed") +
        geom_vline(xintercept = -1.1, colour="grey", size=1,linetype = "dashed") +
        labs(y="-log padj", x="log fold change") +
        scale_colour_manual(name="q-value", values=c("Significant-up"="blue","Significant-down"="red", "Non-sig"="grey50")) +
        theme_bw()
      
      ggsave(paste (id,"volcano.jpeg", sep = "_"), plot = vo)
      
      write.csv(data,"DEG_result.csv",row.names = T)
      
      # Plot heatmap-------------------------------
      
      dds <- estimateSizeFactors(dds)
      
      res_sig <- subset(data, padj < 0.05)
      
      # Extract the normalized counts
      normalized_counts <- counts(dds, normalized=TRUE)
      # Transform the normalized counts 
      vsd_normalized_counts <- vst(normalized_counts, blind = TRUE)
      
      
      
      # Plot the PCA of PCA
      P <- plotPCA(vsd_normalized_counts)
      ggsave(paste (id,"pca.jpeg", sep = "_"), plot = P)
      
      # Subset normalized counts to significant genes
      sig_norm_counts <- normalized_counts[rownames(res_sig), ]
      
      # Choose heatmap color palette
      heat_colors <- brewer.pal( name = "RdBu", n=6)
      
      
      f <- pheatmap(sig_norm_counts, 
                    color = heat_colors, 
                    cluster_rows = T, 
                    show_rownames = F,
                    annotation_col = dplyr::select(colData,condition), 
                    scale = "row",
                    angle_col= 0)
      ggsave(paste (id,"heatmap.jpeg", sep = "_"), plot = f)
      
      
      
    },
    # Specifying error message
    error=function(w){         
      print(paste(i," is missed",sep=""))
    }
    
  )  
}
