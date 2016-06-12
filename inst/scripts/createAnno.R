library(minfi)
manifestFile <- "../../../files_27k/humanmethylation27_270596_v1-2.csv"
maniTmp <- minfi:::read.manifest.27k(manifestFile)
anno <- maniTmp$manifest
manifestList <- maniTmp$manifestList

## Manifest package
IlluminaHumanMethylation27kmanifest <- do.call(IlluminaMethylationManifest,
                                                list(TypeI = manifestList$TypeI,
                                                     TypeII = manifestList$TypeII,
                                                     TypeControl = manifestList$TypeControl,
                                                     TypeSnpI = manifestList$TypeSnpI,
                                                     TypeSnpII = manifestList$TypeSnpII,
                                                     annotation = "IlluminaHumanMethylation27k"))
## Annotation package
anno$IlmnID <- NULL
nam <- names(anno)
names(nam) <- nam
nam[c("AddressA_ID", "AddressB_ID", "AlleleA_ProbeSeq", "AlleleB_ProbeSeq",
"Next_Base", "Color_Channel", "IlmnStrand")] <-  c("AddressA", "AddressB",
                                                                         "ProbeSeqA", "ProbeSeqB",
                                                                         "NextBase", "Color", "strand")

names(nam) <- NULL
names(anno) <- nam
rownames(anno) <- anno$Name
anno <- anno[getManifestInfo(IlluminaHumanMethylation27kmanifest, type = "locusNames"),]
anno$Type <- rep("I", nrow(anno))


Locations <- anno[, c("Chr", "MAPINFO")]
names(Locations) <- c("chr", "pos")
Locations$pos <- as.integer(Locations$pos)
Locations$chr <- paste("chr", Locations$chr, sep = "")
Locations$strand <- ifelse(anno$strand == "F", "+", "-")
table(Locations$chr, exclude = NULL)
rownames(Locations) <- anno$Name
Locations <- as(Locations, "DataFrame")

Manifest <- anno[, c("Name", "AddressA", "AddressB",
                     "ProbeSeqA", "ProbeSeqB", "Type", "NextBase", "Color")]
Manifest <- as(Manifest, "DataFrame")

usedColumns <- c(names(Manifest), names(SNPs.Illumina), 
                 c("CHR", "MAPINFO", "Strand",
                   "Chromosome_36", "Coordinate_36", "Genome_Build"),
                 c("UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island"))
usedColumns <- c(names(Manifest), 
                 c("CHR", "MAPINFO", "Strand"))
Other <- anno[, setdiff(names(anno), usedColumns)]
nam <- names(Other)
Other <- as(Other, "DataFrame")

## We now use an exisitng grSnp object containing a GRanges of relevant SNPs.
## This is created in a separate script

##
## SNP overlap
##

map <- cbind(Locations, Manifest)
map <- GRanges(seqnames = map$chr, ranges = IRanges(start = map$pos, width = 1),
               Strand = map$strand, Type = map$Type)
map <- minfi:::getProbePositionsDetailed(map)
names(map) <- rownames(Locations)

## dbSNP
load("extdata/grSnp146CommonSingle.rda")
SNPs.146CommonSingle <- minfi:::.doSnpOverlap(map, grSnp146CommonSingle)
load("extdata/grSnp144CommonSingle.rda")
SNPs.144CommonSingle <- minfi:::.doSnpOverlap(map, grSnp144CommonSingle)
load("extdata/grSnp142CommonSingle.rda")
SNPs.142CommonSingle <- minfi:::.doSnpOverlap(map, grSnp142CommonSingle)
load("extdata/grSnp141CommonSingle.rda")
SNPs.141CommonSingle <- minfi:::.doSnpOverlap(map, grSnp141CommonSingle)
load("extdata/grSnp138CommonSingle.rda")
SNPs.138CommonSingle <- minfi:::.doSnpOverlap(map, grSnp138CommonSingle)


annoStr <- c(array = "IlluminaHumanMethylation27k",
             annotation = "ilmn12",
             genomeBuild = "hg19")
defaults <- c("Locations", "Manifest",
              "SNPs.146CommonSingle", 
              "Islands.UCSC", "Other")
defaults <- c("Locations", "Manifest",
              "SNPs.146CommonSingle", "Other")

annoObj <-
    IlluminaMethylationAnnotation(list(Locations = Locations,
                                       Manifest = Manifest,
                                       #Islands.UCSC = Islands.UCSC,
                                       Other = Other,
                                       #SNPs.Illumina = SNPs.Illumina,
                                       SNPs.146CommonSingle = SNPs.146CommonSingle,
                                       SNPs.144CommonSingle = SNPs.144CommonSingle,
                                       SNPs.142CommonSingle = SNPs.142CommonSingle,
                                       SNPs.141CommonSingle = SNPs.141CommonSingle,
                                       SNPs.138CommonSingle = SNPs.138CommonSingle
                                       ),
                                  annotation = annoStr, defaults = defaults)
validObject(annoObj)

annoName <- sprintf("%sanno.%s.%s", annoStr["array"], annoStr["annotation"],
                    annoStr["genomeBuild"])
cat("creating object:", annoName, "\n")
assign(annoName, annoObj)
save(list = annoName,
     file = file.path("../../data",paste(annoName, "rda", sep = ".")), compress = "xz")

