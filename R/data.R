#' Estimated Absolute Counts for Yeast Genes
#'
#' This contains estimates for absolute copy numbers per cell of yeast genes.
#' Cells were grown for 24 hours in either Edinburgh minimal media (EMM) or
#' EMM lacking a nitrogen source, resulting in a quiescent state.
#' The copy numbers of 49 calibration mRNAs were estimated using a NanoString
#' assay with 7 external controls, interpolating the log2-transformed
#' calibration mRNA values with the values from the 7 external controls.
#' These were then used to normalize total RNA-Seq data (without any poly-A
#' enrichment or rRNA depletion step), resulting in the genome-wide estimates
#' reported here.
#'
#' @details An important note: not all of the IDs in this table are also found in
#' the \code{\link{yeastAnnos}} \code{data.frame}, but all of the genes that
#' are filtered in a standard sleuth analysis have the same IDs in both tables.
#'
#' @format A \code{data.frame} object with 7289 rows and 7 columns
#' \describe{
#'   \item{PombaseID}{The Systematic ID of the yeast gene on Pombase,
#'     the S. pombe database, found at \url{https://www.pombase.org}}
#'   \item{GeneName}{The common name for the yeast gene}
#'   \item{Ctr1_CPC}{The estimated copy numbers per cell in control
#'     sample #1. \code{NA} values indicate that the gene was not detected}
#'   \item{Ctr2_CPC}{The estimated copy numbers per cell in control
#'     sample #2. \code{NA} values indicate that the gene was not detected}
#'   \item{noN1_CPC}{The estimated copy numbers per cell in nitrogen-starved
#'     sample #1. \code{NA} values indicate that the gene was not detected}
#'   \item{noN2_CPC}{The estimated copy numbers per cell in nitrogen-starved
#'     sample #2. \code{NA} values indicate that the gene was not detected}
#'   \item{fold_change}{Calculated as the ratio between the average nitrogen-starved
#'     copy numbers per cell and the average control copy numbers per cell}
#' }
#' @source The absolute count estimates were taken from columns A, B, K, N,
#'   Q, and T from supplemental table S2 of Marguerat et al. Cell 2012, found
#'   at \url{https://doi.org/10.1016/j.cell.2012.09.019}.
#' @source The direct URL for the supplemental table is 
#'   \url{https://ars.els-cdn.com/content/image/1-s2.0-S0092867412011269-mmc1.xlsx}.
#' @seealso \link{denom} 
"absCounts"

#' Denominator for sleuth-ALR
#'
#' This contains the transcript name of the yeast gene with the most
#' consistent absolute counts, as estimated by a NanoString assay
#' conducted by Marguerat et al. This gene can be treated as
#' a "validated reference gene" to normalize the RNA-Seq data generated
#' in this study using Compositional Normalization.
#'
#' @format A \code{character} vector with one entry containing
#'   the transcript name of the yeast gene.
#'
#' @source The absolute count data was taken from supplemental table S2 of
#'   Marguerat et al. Cell 2012, found at
#'   \url{https://doi.org/10.1016/j.cell.2012.09.019}.
#' @source Direct URL for the raw data is
#'   \url{https://ars.els-cdn.com/content/image/1-s2.0-S0092867412011269-mmc1.xlsx}.
#' @seealso \link{absCounts} 
"denom"

#' Schizosaccharomyces pombe Annotations
#'
#' This contains annotation information for all fission yeast coding and non-coding
#' genes, as described in Ensembl Fungi Genomes Release 37. This can serve as the
#' \code{target_mapping} \code{data.frame} when doing a \code{sleuth} analysis.
#'
#' @details An important note: not all of the IDs in this table are also found in
#' the \code{\link{absCounts}} \code{data.frame}, but all of the genes that
#' are filtered in a standard sleuth analysis have the same IDs in both tables.
#'
#' @format A \code{data.frame} with 7269 rows and 6 columns:
#' \describe{
#'   \item{target_id}{The Ensembl transcript ID for the mRNA}
#'   \item{gene_id}{The Systematic ID used by PomBase, found at
#'     \url{https://www.pombase.org}}
#'   \item{gene_symbol}{The standard gene name used by PomBase}
#'   \item{gene_biotype}{The biotype for the gene. Possible entries are
#'     "ncRNA", "protein_coding", "pseudogene", "RNase_MRP_RNA",
#'     "RNase_P_RNA", "rRNA", "snoRNA", "snRNA", "SRP_RNA", or "tRNA"}
#'   \item{transcript_biotype}{The biotype for the transcript. Same possible
#'     entries as gene_biotype}
#'   \item{full_name}{The full name for the gene product}
#' }
#' @source This information was taken from the headers of the FASTA files containing
#'  coding and non-coding RNAs from Ensembl Fungi Genomes Release 37, found at
#'  \url{http://oct2017-fungi.ensembl.org/}
#' @source The direct URL for the coding cDNA FASTA file is
#'  \url{ftp://ftp.ensemblgenomes.org/pub/release-37/fungi/fasta/schizosaccharomyces_pombe/cdna/Schizosaccharomyces_pombe.ASM294v2.cdna.all.fa.gz}
#' @source The direct URL for the non-coding cDNA FASTA file is
#'  \url{ftp://ftp.ensemblgenomes.org/pub/release-37/fungi/fasta/schizosaccharomyces_pombe/ncrna/Schizosaccharomyces_pombe.ASM294v2.ncrna.fa.gz}
"yeastAnnos"

#' Marguerat et al Fission Yeast Sample to Covariates Table
#'
#' This \code{data.frame} contains the metadata for the poly-A-enriched
#' RNA-Seq data from Marguerat et al. This can be used as the sample-to-covariates
#' table for sleuth analyses.
#'
#' @format A \code{data.frame} with 4 rows and 3 columns:
#' \describe{
#'   \item{sample}{The sample name in ArrayExpress. "pA" indicates poly-A
#'     enriched. These names match the sample names used in Supplemental
#'     Table S2 of Marguerat et al.}
#'   \item{accession}{The accession number on ArrayExpress. This matches the
#'     file name for the raw FASTQ data and the kallisto results directories.}
#'   \item{condition}{The experimental condition, either "control" (Edinburgh
#'     minimal media, EMM) or "noN" (EMM without a nitrogen source)}
#' }
#' @source The raw data and metadata for the Marguerat study is found at
#'   ArrayExpress, \url{https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1154/}
#' @source The metadata is taken directly from the SDRF file from ArrayExpress,
#'   \url{https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-1154/E-MTAB-1154.sdrf.txt}.
"yeastS2C"
