if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

list.of.packages <- c( "optparse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

list.of.bioconductorpackages <-c("panelcn.mops")
new.packages <- list.of.bioconductorpackages[!(list.of.bioconductorpackages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

library(panelcn.mops)
#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("-i", "--inputbam"), type="character", default=NULL, 
              help="input bam file", metavar="character"),
  make_option(c("-c", "--controlbam"), type="character", default=NULL, 
              help="control bam file", metavar="character"),
  make_option(c("-b", "--bed"), type="character", default=NULL, 
              help="gene bed interval file", metavar="character"),
  make_option(c("-o", "--out"), type="character",  
              help="output file name", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$inputbam)){
  print_help(opt_parser)
  stop("input bam missing.n", call.=FALSE)
} else if (is.null(opt$controlbam)){
  print_help(opt_parser)
  stop("control bam missing.n", call.=FALSE)
} else if (is.null(opt$bed)){
  print_help(opt_parser)
  stop("gene bed file missing.n", call.=FALSE)
} 
########################################
Bam_HD753_50ng = "~/harddisk/SSH_server_Data/snakemake_pipeline/snakemake_DNA_seq_pipeline/Output/DNA-HD753-50ng_S8_001/recal/DNA-HD753-50ng_S8.bam"
Bam_control_NA12878 = "~/harddisk/SSH_server_Data/NA12878_control/NA12878.exome.bam"
targetbed= "~/harddisk/SSH_server_Data/variantAnnotation/genesHG38.bed"  ###target bed file

inputBam <- opt$inputbam
controlBam <- opt$controlbam
bed <- opt$bed
#######################################

countWindows <- getWindows(bed)
countWindows_seqName <- c("1", "2",  "3","4" , "5" , "6" , "7" , "8", "9" , "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
                          "20", "21", "22" ,"X"  ,"Y") 
countWindows<- countWindows[countWindows$chromosome %in% countWindows_seqName,]
BAMFiles<- c(inputBam,controlBam)

Bam_HD753test <- countBamListInGRanges(countWindows = countWindows,
                                       bam.files = BAMFiles)

selectedGenes <- countWindows$gene
XandCB <- Bam_HD753test

result <- panelcn.mops(XandCB,testi = 1,
                       classes=c("CN0", "CN1", "CN2", "CN3", "CN4", "CN5", "CN6",
                                 "CN7","CN8","CN16","CN32","CN64","CN128"),
                       I = c(0.025, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 8, 16, 32, 64),norm=1,normType = "poisson",sizeFactor = "mean")
resultlist <- list(result)

sampleNames <- colnames(elementMetadata(Bam_HD753test))
resulttable <- createResultTable(resultlist = resultlist, XandCB = XandCB,
                                 countWindows = countWindows,
                                 selectedGenes = selectedGenes,
                                 sampleNames = sampleNames)


CNV_table <- resulttable[[1]]
CNV_table <- CNV_table[CNV_table$CN!= "CN2",]
CNV_table$cn_alteration <- "amplification"
CNV_table[CNV_table$CN == "CN0" | CNV_table$CN == "CN1",]$cn_alteration <- "deletion"
outputCNV <-CNV_table[c("Gene","cn_alteration")]
outputCNV <- unique(outputCNV)


if (is.null(opt$inputbam)){
  outputName <- paste(tools::file_path_sans_ext(basename(inputBam)), "_cnv",".csv", sep="")
  write.csv(outputCNV,outputName,row.names = F,col.names = T)
} else {
  dir.create(dirname(opt$out))
  write.csv(outputCNV,opt$out,row.names = F,col.names = T)
}

##################

