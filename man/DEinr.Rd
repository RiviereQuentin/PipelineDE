\name{DEinr}
\alias{DEinr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform Differentially Expression analysis in R
}
\description{
Identify differentially expressed genes based directly on sequence alignment maps (in bam format). The function performs the analysis with edgeR and DESeq2.
}
\usage{
DEinr(DEdirectory, isPairedEnd = TRUE, pval.threshold = 0.05, log2FC.threshold = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
%%     ~~Describe \code{x} here~~
  \item{DEdirectory}{character indicating the path to the directory containing the bam files and the gtf file of the reference genome.\code{The bam files corresponding to the control condition must contain in their names "ctrl" or "control" (case-insensitive). The gtf file should be the file used for the read alignments, i.e. the chromosome names in the gtf file must be the same than in the bam files.}}
  \item{isPairedEnd}{boolean indicating whether the sequencing was made in paired-end mode or not. Default = TRUE.\code{If paired-end read data are considered, it is necessary to put in the DEdirectory a bam file for each pair of the read. The names of the bam files must contain "_R1_" for the first mates, "_R2_" for the second ones.}}
  \item{pval.threshold}{numeric indicating the adjusted p-value threshold (false-discovery rate) to detect differentiall expressd genes. Default = 0.05.}
  \item{log2FC.threshold}{numeric indicating the log2 fold-change threshold (in absolute value) to detect differentiall expressd genes. Default = 1.}
}
\details{
With edgeR, Genewise Negative Binomial Generalized Linear Models with Quasi-likelihood Tests are applied; with DESeq2, Negative Binomial GLM fitting and Wald significance tests are performed.}
\value{
data frame with the first column giving the names of the differentially expressed genes (\code(Gene_name)) and the following columns indicating the expression values in each sample (in TPM) and the log2 fold-change and false-discovery rates computed with egeR and DESeq2 (\code(Log2FC_edgeR, FDR_edgeR, Log2FC_DESeq2, FDR_DESeq2)). Note: the columns related to the expression values are named according to the sample names (determined based on the names of the bam files.)
}
\references{
  Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140

  McCarthy DJ, Chen Y and Smyth GK (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. Nucleic Acids Research 40,
  4288-4297

  Chen Y, Lun ATL, Smyth GK (2016). From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline.
  F1000Research 5, 1438

  Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)
}
\author{
Quentin Riviere
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# Run a differential expression analysis based on a toy dataset obtained in Desmodesmus quadricauda.
# In this example, the sequencing was performed in paired-end mode.
# Default parameters are considered: FDR threshold set to 0.05 and |log2 FC| to 1.

DE.res <- DEinr(DEdirectory = system.file("extdata", package = "PipelineDE"))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
