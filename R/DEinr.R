DEinr <- function(DEdirectory, isPairedEnd = TRUE, pval.threshold = 0.05, log2FC.threshold = 1) {
  DEfiles <- list.files(DEdirectory)
  bamfiles <- DEfiles[grep("bam", DEfiles)]
  gtffile <- DEfiles[grep("gtf", DEfiles)]

  #==========================================================
  #Get the count table
  #==========================================================
  if (isPairedEnd){
    R1bam <- bamfiles[grep("_R1_", bamfiles)]
    R2bam <- bamfiles[grep("_R2_", bamfiles)]

    R1bam <- R1bam[order(R1bam)]
    R2bam <- R2bam[order(R2bam)]

    samplenames <- unlist(strsplit(R1bam, "_R1"))[seq(1, 2*length(R1bam), 2)]

    R1R2 <- list()
    for (i in seq_along(samplenames)){
      Rsamtools::mergeBam(file.path(DEdirectory, c(R1bam[i], R2bam[i])),
                                                    destination = file.path(DEdirectory, paste0(samplenames[i], ".bam"))
                                                    , overwrite = TRUE)
      R1R2[[samplenames[i]]] <- file.path(DEdirectory, paste0(samplenames[i], ".bam"))
    }

    readcounts <- lapply(R1R2, Rsubread::featureCounts,
                         annot.ext = file.path(DEdirectory, gtffile),
                         isGTFAnnotationFile = TRUE,
                         nthreads = 2)

    readcounts <- lapply(readcounts, function(x) {return(x$counts)})
    readcounts <- do.call(cbind, readcounts)
  } else {

    samplenames <- unlist(strsplit(bamfiles, '\\.bam$'))

    R0 <- as.list(file.path(DEdirectory, bamfiles))
    names(R0) <- samplenames

    readcounts <- lapply(R0, Rsubread::featureCounts,
                         annot.ext = file.path(DEdirectory, gtffile),
                         isGTFAnnotationFile = TRUE,
                         nthreads = 2)

    readcounts <- lapply(readcounts, function(x) {return(x$counts)})
    readcounts <- do.call(cbind, readcounts)
  }

  #==========================================================
  #Perform the analysis with edgeR
  #==========================================================
  ctrlgroup <- c(grep("control", colnames(readcounts), ignore.case = TRUE),
                 grep("ctrl", colnames(readcounts), ignore.case = TRUE)) # BAM files pointing to control samples should include control or ctrl in their name (case-insensitive)
  treatgroup <- seq_len(ncol(readcounts))[!(seq_len(ncol(readcounts)) %in% ctrlgroup)]
  group <- rep("ctrl", ncol(readcounts))
  group[treatgroup] <- "treatment"
  group <- as.factor(group)
  y <- edgeR::DGEList(counts=readcounts,group=group)
  y <-  edgeR::calcNormFactors(y)
  design <- stats::model.matrix(~group)
  y <- edgeR::estimateDisp(y,design)

  fit <- edgeR::glmQLFit(y,design) #Quasi likelihood test
  qlf <- edgeR::glmQLFTest(fit,coef=2)
  DEtable <- topTags(qlf, n = nrow(readcounts), p.value = 1)
  DEtable.edgeR <- DEtable$table

  #==========================================================
  #Perform the analysis with DESeq2
  #==========================================================
  if (isPairedEnd){
    coldata <- data.frame(condition = group,
                          type = "paired-end",
                          row.names = colnames(readcounts))
  } else {
    coldata <- data.frame(condition = group,
                          type = "single-read",
                          row.names = colnames(readcounts))
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = readcounts,
                                        colData = coldata,
                                        design = ~ condition)
  featureData <- data.frame(gene=rownames(readcounts))
  S4Vectors::mcols(dds) <- S4Vectors::DataFrame(mcols(dds), featureData)
  dds <- DESeq2::DESeq(dds)
  DEtable.DESeq2 <- DESeq2::results(dds, pAdjustMethod = "fdr")
  DEtable.DESeq2 <- DEtable.DESeq2[!is.na(DEtable.DESeq2$padj),]

  #==========================================================
  #Compute TPM
  #==========================================================
  gtfinfo <- rtracklayer::import(file.path(DEdirectory, gtffile))
  gtfinfo <- gtfinfo[gtfinfo$type == "exon",]
  gtfinfo <- data.frame(gene = gtfinfo$gene_id,
                        width = width(gtfinfo))
  gtfinfo$gene <- as.factor(gtfinfo$gene)
  genelengths <- tapply(gtfinfo$width, gtfinfo$gene, sum)
  genelengths <- genelengths[rownames(readcounts)]

  genelengths <- as.vector(genelengths)
  tpmcounts <- rpkm(y, gene.length = genelengths)
  tpmcounts <- t( t(tpmcounts) / colSums(tpmcounts) ) * 1e6

  #==========================================================
  #Select DE genes
  #==========================================================

  DE.edgeR <- rownames(DEtable.edgeR[(DEtable.edgeR$FDR <= pval.threshold) &
                              (abs(DEtable.edgeR$logFC) >= log2FC.threshold), ])

  DE.DESeq2 <- rownames(DEtable.edgeR[(DEtable.DESeq2$padj <= pval.threshold) &
                              (abs(DEtable.DESeq2$log2FoldChange) >= log2FC.threshold), ])

  DEgenes <- c(DE.edgeR, DE.DESeq2)
  DEgenes <- DEgenes[!duplicated(DEgenes)]

  DE.table <- cbind(as.matrix(tpmcounts[DEgenes,]),
                    as.matrix(DEtable.edgeR[DEgenes, c("logFC", "FDR")]),
                    as.matrix(DEtable.DESeq2[DEgenes, c("log2FoldChange", "padj")]))

  colnames(DE.table) <- gsub("[.]bam", "", colnames(DE.table))
  colnames(DE.table)[(ncol(DE.table)-3):ncol(DE.table)] <- c("Log2FC_edgeR", "FDR_edgeR", "Log2FC_DESeq2", "FDR_DESeq2")

  DE.table <- data.frame(Gene_name = rownames(DE.table),
                         DE.table)
  DE.table <- DE.table[!is.na(DE.table[,ncol(DE.table)]),]

  DE.table <- DE.table[order(apply(DE.table[, c("Log2FC_edgeR", "Log2FC_DESeq2")], 1, mean)),]

  return(DE.table)
}
