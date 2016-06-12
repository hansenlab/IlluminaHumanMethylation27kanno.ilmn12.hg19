processUCSCsnp <- function(snpfile) {
    require(GenomicRanges)
    cat("Reading file\n")
    df <- read.delim(gzfile(snpfile), header = FALSE,
                     stringsAsFactors = FALSE)
    names(df) <- c("chr", "start", "end", "name", "strand",
                   "refNCBI", "class", "alleleFreqs")
    print(table(df$chr))
    cat("Only keeping chrs 1-22, X, Y\n")
    df <- df[df$chr %in% paste0("chr", c(1:22, "X", "Y")),]
    print(table(df$class))
    cat("Only keeping class 'single'\n")
    df <- df[df$class == "single",]
    cat("Computing MAF\n")
    df$alleleFreqs <- sub(",$", "", df$alleleFreqs)
    sp <- strsplit(df$alleleFreqs, ",")
    minFreq <- sapply(sp, function(xx) min(as.numeric(xx)))
    cat("Instantiating object\n")
    grSNP <- GRanges(seqnames = df$chr, strand = df$strand, ranges = IRanges(start = df$start + 1, end = df$end),
                     MAF = minFreq, ref = df$refNCBI)
    names(grSNP) <- df$name
    grSNP
}

grSnp146CommonSingle <- processUCSCsnp("extdata/snp146Common_small.txt.gz")
save(grSnp146CommonSingle, file = "extdata/grSnp146CommonSingle.rda")

grSnp144CommonSingle <- processUCSCsnp("extdata/snp144Common_small.txt.gz")
save(grSnp144CommonSingle, file = "extdata/grSnp144CommonSingle.rda")

grSnp142CommonSingle <- processUCSCsnp("extdata/snp142Common_small.txt.gz")
save(grSnp142CommonSingle, file = "extdata/grSnp142CommonSingle.rda")

grSnp141CommonSingle <- processUCSCsnp("extdata/snp141Common_small.txt.gz")
save(grSnp141CommonSingle, file = "extdata/grSnp141CommonSingle.rda")

grSnp138CommonSingle <- processUCSCsnp("extdata/snp138Common_small.txt.gz")
save(grSnp138CommonSingle, file = "extdata/grSnp138CommonSingle.rda")