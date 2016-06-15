library(minfi)
manifestFile <- "../../../IlluminaHumanMethylation27k_files/data/HumanMethylation27_270596_v.1.2.bpm"
if(!file.exists(manifestFile) || !file.exists("extdata")) {
    cat("Missing files, quitting\n")
    q(save = "no")
}

maniTmp <- minfi:::read.manifest.27k(manifestFile)
anno <- maniTmp$manifest
manifestList <- maniTmp$manifestList

# Adding colors to Type I SNP probes as discovered experimentally by Tim Triche:
manifestList$TypeSnpI[match(
        c("rs1019916", "rs10457834", "rs1416770", "rs1941955", "rs2125573", 
        "rs2235751", "rs2521373", "rs264581", "rs2804694", "rs2959823", 
        "rs5931272", "rs6546473", "rs739259", "rs798149", "rs845016", 
        "rs866884"), manifestList$TypeSnpI$Name), "Color"] <- 
        c("Red", "Red", "Red", "Red", "Red", "Grn", "Grn", "Red", "Red", 
        "Red", "Red", "Grn", "Red", "Grn", "Grn", "Red")

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


Locations <- anno[, c("Chr", "MapInfo")]
names(Locations) <- c("chr", "pos")
Locations$pos <- as.integer(Locations$pos)
Locations$chr <- paste("chr", Locations$chr, sep = "")
Locations$strand <- ifelse(anno$strand == "TOP", "+", "-")
table(Locations$chr, exclude = NULL)
rownames(Locations) <- anno$Name
Locations <- as(Locations, "DataFrame")

Manifest <- anno[, c("Name", "AddressA", "AddressB",
                     "ProbeSeqA", "ProbeSeqB", "Type", "NextBase", "Color")]
Manifest <- as(Manifest, "DataFrame")


usedColumns <- c(names(Manifest), 
                 c("Chr", "MapInfo", "Strand"))
Other <- anno[, setdiff(names(anno), usedColumns)]
nam <- names(Other)
Other <- as(Other, "DataFrame")
Other[,"X"] <- NULL
Other[,"X.1"] <- NULL
## We now use an exisitng grSnp object containing a GRanges of relevant SNPs.
## This is created in a separate script

##
## SNP overlap
##

map <- cbind(Locations, Manifest)
map <- GRanges(seqnames = map$chr, ranges = IRanges(start = map$pos, width = 1),
               Strand = map$strand, Type = map$Type)
map <- minfi:::.getProbePositionsDetailed(map)
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
load("extdata/grSnp137CommonSingle.rda")
SNPs.137CommonSingle <- minfi:::.doSnpOverlap(map, grSnp137CommonSingle)
load("extdata/grSnp135CommonSingle.rda")
SNPs.135CommonSingle <- minfi:::.doSnpOverlap(map, grSnp135CommonSingle)
load("extdata/grSnp132CommonSingle.rda")
SNPs.132CommonSingle <- minfi:::.doSnpOverlap(map, grSnp132CommonSingle)


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
                                       Other = Other,
                                       SNPs.146CommonSingle = SNPs.146CommonSingle,
                                       SNPs.144CommonSingle = SNPs.144CommonSingle,
                                       SNPs.142CommonSingle = SNPs.142CommonSingle,
                                       SNPs.141CommonSingle = SNPs.141CommonSingle,
                                       SNPs.138CommonSingle = SNPs.138CommonSingle,
                                       SNPs.137CommonSingle = SNPs.137CommonSingle,
                                       SNPs.135CommonSingle = SNPs.135CommonSingle,
                                       SNPs.132CommonSingle = SNPs.132CommonSingle
                                       ),
                                  annotation = annoStr, defaults = defaults)
validObject(annoObj)

annoName <- sprintf("%sanno.%s.%s", annoStr["array"], annoStr["annotation"],
                    annoStr["genomeBuild"])
cat("creating object:", annoName, "\n")
assign(annoName, annoObj)
save(list = annoName,
     file = file.path("../../data",paste(annoName, "rda", sep = ".")), compress = "xz")

sessionInfo()
q(save = "no")

