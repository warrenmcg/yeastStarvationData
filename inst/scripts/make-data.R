## Step 1A: Download the Yeast Annotations and create the kallisto index
main_url <- "ftp://ftp.ensemblgenomes.org/pub/release-37/fungi/fasta/schizosaccharomyces_pombe"
cdna_url <- file.path(main_url, 'cdna/Schizosaccharomyces_pombe.ASM294v2.cdna.all.fa.gz')
ncrna_url <- file.path(main_url, 'ncrna/Schizosaccharomyces_pombe.ASM294v2.ncrna.fa.gz')
download.file(cdna_url, 'ASM294v2.cdna.all.fa.gz')
download.file(ncrna_url, 'ASM294v2.ncrna.fa.gz')

cdnas <- Biostrings::readDNAStringSet('ASM294v2.cdna.all.fa.gz')
ncrnas <- Biostrings::readDNAStringSet('ASM294v2.ncrna.fa.gz')
file.remove('ASM294v2.cdna.all.fa.gz')
file.remove('ASM294v2.ncrna.fa.gz')

all_rnas <- c(cdnas, ncrnas)
Biostrings::writeXStringSet(all_rnas, 'ASM294v2.pombase.all.fa.gz', compress = TRUE)

if (Sys.which('kallisto') == "") {
  stop("kallisto is not part of your system's $PATH. ",
       "Please make sure it is installed and on your $PATH.")
} else {
  v_string <- system2("kallisto", "version", stdout = TRUE)
  message(paste("Using", v_string))
}
args <- c('index', '-i ASM294v2.pombase.all.kidx', 'ASM294v2.pombase.all.fa.gz')
system2('kallisto', args)

## Step 1B: Download the raw FASTQ files from ArrayExpress

data_dir <- '../extdata'
samples <- paste0('ERR13590', 6:9)
url <- 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR135'
sample_files <- paste0(samples, '.fastq.gz')
url_files <- file.path(url, samples, sample_files)

out_files <- file.path(data_dir, sample_files)

for (i in 1:length(samples)) {
  download.file(url_files[i], out_files[i])
  ## Step 2: Process the FASTQ files using kallisto
  ## parameters for kallisto run:
  ## kallisto quant -i {YEAST_KALLISTO_INDEX} -b 100 -l 200 -s 20 --single \
  ##   --fr-stranded -o {OUTPUT_FILE} {INPUT_FILE}
  args <- c('quant', '-i ASM294v2.pombase.all.kidx', '-b 100', '-l 200',
            '-s 20', '--single', '--fr-stranded',
            paste('-o', file.path(data_dir, samples[i])),
            shQuote(out_files[i]))
  system2('kallisto', args)
  file.remove(out_files[i])
}

## Remove the FASTA file and kallisto index, which are no longer needed
file.remove('ASM294v2.pombase.all.kidx')
file.remove('ASM294v2.pombase.all.fa.gz')

## Step 3: Download the absolute count data from Lovell et al 2015

url <- 'https://ars.els-cdn.com/content/image/1-s2.0-S0092867412011269-mmc1.xlsx'
download.file(url, '../extdata/mmc1.xlsx')
table <- openxlsx::read.xlsx('../extdata/mmc1.xlsx', sheet = 3, startRow = 2)
file.remove('../extdata/mmc1.xlsx')

table <- dplyr::select(table, 1, 2, 11, 14, 17, 20)
table$fold_change <- log2(rowMeans(as.matrix(table[,5:6])) / rowMeans(as.matrix(table[,3:4])))

# Replace old IDs with current IDs
old_ids <- c("SPBC713.13", "SPAC823.02", "SPAC1556.06.1",
             "SPNCRNA.98", "SPSNRNA.03", "SPSNORNA.40")
new_ids <- c("SPBC713.14c", "SPAC823.17", "SPAC1556.06",
             "ENSRNA049676787", "ENSRNA049676282", "ENSRNA049677584")
table$Systematic.name[match(old_ids, table$Systematic.name)] <- new_ids

# Determine which gene is most consistent, indicating it is a valid reference gene
cov <- matrixStats::rowSds(as.matrix(table[,3:6])) / rowMeans(as.matrix(table[,3:6]))
names(cov) <- table$Systematic.name
denom_gene <- names(which(cov == min(cov[which(cov>0)], na.rm = T)))

absCounts <- table
colnames(absCounts) <- c("PombaseID", "GeneName", "Ctr1_CPC", "Ctr2_CPC",
  "noN2_1_CPC", "noN2_2_CPC", "fold_change")
save(absCounts, file = '../../data/absCounts.rda')
denom <- annos$target_id[which(annos$gene_id == denom_gene)]
save(denom, file = '../../data/denom.rda')

## Step 4: Collect the metadata and yeast annotations necessary for sleuth
## Metadata
url <- "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-1154/E-MTAB-1154.sdrf.txt"

table <- read.table(url, sep = "\t", header = T, stringsAsFactors = F)
all_samples <- dplyr::select(table, sample = Source.Name,
                          accession = Comment.ENA_RUN.,
                          condition = FactorValue.MEDIA.)
yeastS2C <- all_samples[grepl("pA$", all_samples$sample), ]
yeastS2C$condition <- as.factor(yeastS2C$condition)
levels(yeastS2C$condition) <- c('control', 'noN2')

save(yeastS2C, file = '../../data/yeastS2C.rda')

## Yeast Annotations
splits <- strsplit(names(all_rnas), split = " ", fixed = T)
data <- lapply(splits, function(x) {
   target_id <- x[1]
   gene_id <- strsplit(x[4], ":", fixed = T)[[1]][2]
   gene_symbol <- strsplit(x[7], ":", fixed = T)[[1]][2]
   gene_biotype <- strsplit(x[5], ":", fixed = T)[[1]][2]
   transcript_biotype <- strsplit(x[6], ":", fixed = T)[[1]][2]
   if (length(x) < 8)
     description <- NA
   else {
     description <- paste(x[8:length(x)], collapse = " ")
     description <- gsub("description:", "", description)
   }
   c(target_id = target_id, gene_id = gene_id, gene_symbol = gene_symbol,
     gene_biotype = gene_biotype, transcript_biotype = transcript_biotype,
     full_name = description)
})

annos <- do.call(rbind, data)
yeastAnnos <- as.data.frame(annos, stringsAsFactors = FALSE)

save(yeastAnnos, file = '../../data/yeastAnnos.rda')
