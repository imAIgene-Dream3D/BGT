###########################################################################
#
# Peter Brazda
# 23 February 2024
# Helper functions for SORTseq pre-processing
#
###########################################################################



suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(readxl))
suppressMessages(library(ggplot2))
suppressMessages(library(platetools))
suppressMessages(library(scales))
suppressMessages(library(plyr))
suppressMessages(library(png))
suppressMessages(library(gridExtra))
suppressMessages(library(platetools))
suppressMessages(library(data.table))
suppressMessages(library(viridis))
suppressMessages(library(readxl))

Prepare_SORTseq_QC_report <- function(plateID, matrices, primer, gencode) {
  # Load and prepare count data
  counts_file <- paste0(matrices, plateID, ".transcripts.txt")
  counts <- fread(counts_file, header = TRUE)
  counts <- separate(counts, col = "GENEID", into = c("ensembl_gene", "chr"), sep = "__")
  
  # Merge with gencode and filter
  DT <- merge(gencode, counts, by = "ensembl_gene")
  DT <- subset(DT, !ensembl_gene %in% c("ENSG00000213380",
                                        "ENSG00000180155",
                                        "ENSG00000214338",
                                        "ENSG00000158427",
                                        "ENSG00000239605"))
  setnames(DT, c("gene_name", "chr.x"), c("gene", "chr"))
  DT <- subset(DT, gene_type %in% c("ERCC_spikein", "protein_coding"))
  
  # Identify and remove cluster names
  remove_cluster_name <- DT[grepl("cluster", DT[["chr.y"]]), ]
  DT <- DT[!DT$gene %in% remove_cluster_name$gene, ]
  
  # Organize reads and remove unwanted genes
  reads1 <- DT[, c(2, 3, 10:ncol(DT))]
  ERCC <- reads1[reads1$chr %like% "ERCC", ]
  chrM <- subset(gencode, chr == "chrM")
  chrM <- reads1[reads1$gene %in% intersect(reads1$gene, chrM$gene_name), ]
  
  # Primer and plate distribution preparation
  dTV_numbers <- as.vector(t(primer))
  plate <- data.frame(dTV_numbers)
  plate$plate_order <- 1:384
  
  # Filter reads for mRNA, excluding chrM and ERCC genes
  reads <- subset(reads1, !gene %in% union(chrM$gene, ERCC$gene))
  
  # Set column names and transform data frames
  colnames_transform <- function(df) {
    colnames(df) <- c("gene", "chr", plate$plate_order)
    df %>% remove_rownames() %>% column_to_rownames(var = "gene") %>% select(-chr)
  }
  
  reads <- colnames_transform(reads)
  chrM <- colnames_transform(chrM)
  ERCC <- colnames_transform(ERCC)
  
  # Data transformation for mRNA, mitRNA, and ERCC
  transform_data <- function(df) {
    df_t <- as.data.frame(t(df))
    rownames(df_t) <- as.integer(rownames(df_t))
    df_t <- cbind(rownames(df_t), df_t, row.names = NULL)
    colnames(df_t)[1] <- "wells"
    df_t <- merge(plate, df_t, by.x = "dTV_numbers", by.y = "wells")
    df_t$dTV_numbers <- NULL
    row.names(df_t) <- make.names(df_t[,1], unique = TRUE)
    df_t[,1] <- NULL
    row.names(df_t) <- sub("X", "", row.names(df_t))
    as.data.frame(t(df_t))
  }
  
  mRNA <- transform_data(reads)
  chrM <- transform_data(chrM)
  ERCC <- transform_data(ERCC)
  
  # Compute total reads
  compute_counts <- function(df, name) {
    counts <- data.frame(colSums(df))
    colnames(counts) <- name
    setDT(counts, keep.rownames = TRUE)[]
  }
  
  mRNA_count <- compute_counts(mRNA, "mRNA_reads")
  mit_count <- compute_counts(chrM, "mit_reads")
  ERCC_count <- compute_counts(ERCC, "ERCC_reads")
  
  # Merge all read counts
  all_reads <- merge(mRNA_count, mit_count, by = "rn", all = TRUE)
  all_reads <- merge(all_reads, ERCC_count, by = "rn", all = TRUE)
  colnames(all_reads) <- c("wells", "mRNA_reads", "mit_reads", "ERCC_reads")
  
  # Prepare for plotting
  all_reads_dt <- as.data.table(all_reads)
  plot_reads <- melt(all_reads_dt, id.vars = "wells")
  plot_reads$wells <- as.numeric(as.character(plot_reads$wells))
  

  ggplot(plot_reads, aes(fill=variable, y=value, x=wells)) + 
    geom_bar(stat="identity") + theme_classic() + scale_x_continuous(breaks = c(1,50,100,150,200,250,300,350,384))+theme_classic() + ggtitle("Total read counts")
  ggsave("01_Total read counts_full.png", plot = last_plot(), width = 10, height = 7)
  
  
  ggplot(plot_reads, aes(fill=variable, y=value, x=wells)) +
    geom_bar(stat="identity", position="fill") + scale_x_continuous(breaks = c(1,50,100,150,200,250,300,350,384)) + 
    scale_y_continuous(labels = scales::percent) + theme_classic() + ylab("percentage") +theme_classic() + ggtitle("Proportional read counts")
  ggsave("02_Proportional read counts.png", plot = last_plot(), device = "png", width = 10, height = 7)
  
  
  ### percentage of molecules
  all_reads$wells <- as.numeric(gsub("CelSeq.dTV", "", all_reads$wells))
  all_reads$wells <- as.numeric(all_reads$wells)   # make the wells numeric
  all_reads$total <- rowSums(all_reads[,-1])    # sum of all rows per well 
  all_reads <- all_reads[order(all_reads$wells)]   # order the wells from 1 to 384
  
  # calculate percentages
  percentage <- all_reads
  percentage$mRNA <- percentage$mRNA_reads/percentage$total*100
  percentage$mit <- percentage$mit_reads/percentage$total*100
  percentage$ERCC <- percentage$ERCC_reads/percentage$total*100
  #percentage[,2:5] <- NULL
  percentage[,2:4] <- NULL # keep total column
  
  #mRNA
  wells <- num_to_well(percentage$wells, plate = 384)
  vals <- percentage$mRNA
  df384_mRNA <- data.frame(wells, vals)
  
  raw_map(data = df384_mRNA$vals,
          well = df384_mRNA$wells,
          plate = 384) +
    ggtitle("mRNA counts (%)") +
    scale_fill_gradientn(colors = viridis_pal(option = "C")(100)) + theme_classic()
  ggsave("03_Percentage of mRNA molecules.png", plot = last_plot(), device = "png", width = 5.5, height = 3.5) 
  
  
  #mitochondrial
  wells <- num_to_well(1:384, plate = 384)
  vals <- percentage$mit
  df384_mit <- data.frame(wells, vals)
  
  raw_map(data = df384_mit$vals,
          well = df384_mit$wells,
          plate = 384) +
    ggtitle("mitoRNA counts (%)") +
    scale_fill_gradientn(colors = viridis_pal(option = "C")(100)) + theme_classic()
  ggsave("04_Percentage of mito molecules.png", plot = last_plot(), device = "png", width = 5.5, height = 3.5) 
  
  
  #ERCC
  wells <- num_to_well(1:384, plate = 384)
  vals <- percentage$ERCC
  df384_ERCC <- data.frame(wells, vals)
  
  raw_map(data = df384_ERCC$vals,
          well = df384_ERCC$wells,
          plate = 384) +
    ggtitle("ERCC counts (%)") +
    scale_fill_gradientn(colors = viridis_pal(option = "C")(100)) + theme_classic()
  ggsave("05_Percentage of ERCC molecules.png", plot = last_plot(), device = "png", width = 5.5, height = 3.5) 
  
  
  ##### number of molecules
  #mRNA
  wells <- num_to_well(1:384, plate = 384)
  vals <- all_reads$mRNA_reads
  df384_mRNA <- data.frame(wells, vals)
  
  raw_map(data = df384_mRNA$vals,
          well = df384_mRNA$wells,
          plate = 384) +
    ggtitle("mRNA counts (#)") +
    scale_fill_gradientn(colors = viridis_pal(option = "C")(100)) + theme_classic()
  ggsave("06_Number of mRNA molecules.png", plot = last_plot(), device = "png", width = 5.5, height = 3.5) 
  
  
  #mitochondrial
  wells <- num_to_well(1:384, plate = 384)
  vals <- all_reads$mit_reads
  df384_mit <- data.frame(wells, vals)
  
  raw_map(data = df384_mit$vals,
          well = df384_mit$wells,
          plate = 384) +
    ggtitle("mitoRNA counts (#)") +
    scale_fill_gradientn(colors = viridis_pal(option = "C")(100)) + theme_classic()
  ggsave("07_Number of mito molecules.png", plot = last_plot(), device = "png", width = 5.5, height = 3.5) 
  
  
  #ERCC
  wells <- num_to_well(1:384, plate = 384)
  vals <- all_reads$ERCC_reads
  df384_ERCC <- data.frame(wells, vals)
  
  raw_map(data = df384_ERCC$vals,
          well = df384_ERCC$wells,
          plate = 384) +
    ggtitle("ERCC counts (#)") +
    scale_fill_gradientn(colors = viridis_pal(option = "C")(100)) + theme_classic()
  ggsave("08_Number of ERCC molecules.png", plot = last_plot(), device = "png", width = 5.5, height = 3.5) 
  
  
  
  ###### gene coverage
  cover <- mRNA
  cover[cover>1] <- 1   # every count more than 1 change to 1
  genes <- data.frame(colSums(cover))    # sum of genes per cell, number of genes per cell
  setDT(genes, keep.rownames = TRUE)[]
  colnames(genes) <- c("wells","gene_coverage")
  genes$wells <- as.numeric(gsub("CelSeq.dTV", "", genes$wells))
  genes <- genes[order(genes$wells)]
  
  #plate
  wells <- num_to_well(1:384, plate = 384)
  vals <- genes$gene_coverage
  df384_genes <- data.frame(wells, vals)
  
  raw_map(data = df384_genes$vals,
          well = df384_genes$wells,
          plate = 384) +
    ggtitle("Gene coverage (#)") +
    scale_fill_gradientn(colors = viridis_pal(option = "C")(100)) + theme_classic()
  ggsave("09_Gene coverage.png", plot = last_plot(), device = "png", width = 5.5, height = 3.5)
  
  
  ##### gene complexity
  ggplot(data=genes, aes(x=genes$gene_coverage)) + 
    geom_histogram(breaks=seq(0,5000, by=50), col = NA, fill = "#E87D72", aes()) +
    # scale_fill_gradientn(colours = c("red","pink", "lightblue", "#6E9BF8"), values = c(0,0.01,0.1,1)) + 
    labs(title="Gene complexity", x="Gene complexity", y="Frequency") + theme_classic()
  ggsave("10_Gene_complexity_hist.png", plot = last_plot(), device = "png", width = 10, height = 7)
  
  
  #### ERC complexity
  ERCC_cover <- ERCC
  ERCC_cover[ERCC_cover>1] <- 1
  
  coverage <- data.frame(colSums(ERCC_cover))
  setDT(coverage, keep.rownames = TRUE)[]
  colnames(coverage) <- c("wells", "gene_coverage")
  coverage$wells <- as.numeric(gsub("CelSeq.dTV", "", coverage$wells))
  coverage <- coverage[order(coverage$wells)]
  
  ggplot(data=coverage, aes(x=coverage$gene_coverage)) + 
    geom_histogram(breaks=seq(0,60, by=1), col = NA, alpha = 1, fill = "#6E9BF8") +
    labs(title="ERCC complexity", x="Gene coverage", y="Frequency") + theme_classic()
  ggsave("11_ERCC complexity_hist.png", plot = last_plot(), device = "png", width = 10, height = 7)
  
  
  ######### ERCC (and mit) filtering 
  ERCC15 <- percentage
  ggplot(ERCC15, aes(x = reorder(wells, -ERCC), y = ERCC)) + geom_bar(stat = "identity", col = NA, fill = "#6E9BF8", alpha = 1) + scale_x_discrete(name="Wells") + theme(axis.text.x = element_text(size=5, angle = 90, hjust = 1)) + scale_y_continuous(name = "ERCC(%)")+theme_classic()
  ggsave("12_ERCC percentage_rank.png", plot = last_plot(), device = "png", width = 10, height = 7)
  
  ggplot(ERCC_count, aes(x = reorder(rn, -ERCC_reads), y = ERCC_reads)) + geom_bar(stat = "identity", col = NA, fill = "#6E9BF8", alpha = 1) + scale_x_discrete(name="Wells") + theme(axis.text.x = element_text(size=5, angle = 90, hjust = 1)) + scale_y_continuous(name = "ERCC(#)")+theme_classic()
  ggsave("13_ERCC count_rank.png", plot = last_plot(), device = "png", width = 10, height = 7)
  
  
  
  ggplot(percentage, aes(x = reorder(wells, -mit), y = mit)) + geom_bar(stat = "identity", col = NA, fill = "#53B74C", alpha = 1) + scale_x_discrete(name = "wells") + scale_y_continuous(name = "mito(%)") + theme(axis.text.x = element_text(size=5, angle = 90, hjust = 1))+theme_classic()
  ggsave("14_mito percentage_rank.png", plot = last_plot(), device = "png", width = 10, height = 7) 
  
  ggplot(mit_count, aes(x = reorder(rn, -mit_reads), y = mit_reads)) + geom_bar(stat = "identity", col = NA, fill = "#53B74C", alpha = 1) + scale_x_discrete(name="Wells") + theme(axis.text.x = element_text(size=5, angle = 90, hjust = 1)) + scale_y_continuous(name = "mito(#)")+theme_classic()
  ggsave("15_mito count_rank.png", plot = last_plot(), device = "png", width = 10, height = 7)
  
  
  
  ggplot(percentage, aes(x = reorder(wells, -mRNA), y = mRNA)) + geom_bar(stat = "identity", col = NA, fill = "#E87D72", alpha = 1) + scale_x_discrete(name = "wells") + scale_y_continuous(name = "mRNA(%)") + theme(axis.text.x = element_text(size=5, angle = 90, hjust = 1))+theme_classic()
  ggsave("16_mRNA percentage_rank.png", plot = last_plot(), device = "png", width = 10, height = 7)
  
  
  ggplot(mRNA_count, aes(x = reorder(rn, -mRNA_reads), y = mRNA_reads)) + geom_bar(stat = "identity", col = NA, fill = "#E87D72", alpha = 1) + scale_x_discrete(name="Wells") + theme(axis.text.x = element_text(size=5, angle = 90, hjust = 1)) + scale_y_continuous(name = "mRNA(#)")+theme_classic()
  ggsave("17_mRNA count_rank.png", plot = last_plot(), device = "png", width = 10, height = 7)
  
  
}

Filter_SORTseq_tables <- function(plate, input_location, output_location_nME, output_location_nE, mito_pct_max = 30, ercc_content_max = 25) {
  file <- paste(plate, ".transcripts.txt", sep = "")
  output_file_nME <- paste(plate, "_noME.tab", sep = "")
  output_file_nE <- paste(plate, "_noE.tab", sep = "")
  
  DT <- fread(paste(input_location, file, sep = ""), header = TRUE)
  DT <- separate(DT, col = "GENEID", into = c("ensembl_gene", "chr"), sep = "__")
  DT <- merge(gencode, DT, by = "ensembl_gene")
  DT <- subset(DT, gene_type %in% c("ERCC_spikein", "protein_coding"))
  removeclustername <- DT[grepl("cluster", DT[["chr.y"]]), ]
  DT <- DT[!DT$gene_name %in% removeclustername$gene_name, ] # '-cluster' names
  DT <- subset(DT, !ensembl_gene %in% c("ENSG00000213380", "ENSG00000180155", "ENSG00000214338", "ENSG00000158427", "ENSG00000239605"))
  
  DT$sum <- rowSums(DT[, 10:ncol(DT)])
  
  mat <- data.matrix(DT[, 10:ncol(DT)])
  rownames(mat) <- DT$gene_name
  
  # Get mitochondrial genes and ERCC spikeins
  mito.genes <- subset(gencode, chr == "chrM")
  ERCC <- subset(gencode, chr == "ERCC")
  mito.ID <- which(rownames(mat) %in% mito.genes$gene_name)
  ERCC.ID <- which(rownames(mat) %in% ERCC$gene_name)
  
  # Filter cells based on mitochondrial content and ERCC content
  cell.high.mito.content <- which(colSums(mat[mito.ID, ]) / colSums(mat[-ERCC.ID, ]) * 100 > mito_pct_max)
  cell.high.ERCC.content <- which(colSums(mat[ERCC.ID, ]) / colSums(mat) * 100 > ercc_content_max)
  out <- union(cell.high.mito.content, cell.high.ERCC.content)
  mat.filtered <- mat[, -out]
  
  DTX <- as.data.table(mat.filtered, keep.rownames = "gene_name")
  DTX$sum <- NULL
  
  `%notin%` <- Negate(`%in%`)
  DTY <- DTX[!(gene_name %in% ERCC$gene_name)]
  DTM <- DTY[!(gene_name %in% mito.genes$gene_name)]
  
  # Export filtered data
  fwrite(DTM, file = paste(output_location_nME, output_file_nME, sep = ""), sep = "\t")
  fwrite(DTY, file = paste(output_location_nE, output_file_nE, sep = ""), sep = "\t")
}
