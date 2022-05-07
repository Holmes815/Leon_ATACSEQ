#Cargar las librerías
library(BSgenome.Hsapiens.UCSC.hg19)
library(Rsubread)
#Creando genoma de referencia para alinear las lecturas
mainChromosomes <- paste0("chr", c(1:21, "X", "Y","M"))
mainChrSeq <- lapply(mainChromosomes,
                     function(x)BSgenome.Hsapiens.UCSC.hg19[[x]])
names(mainChrSeq) <- mainChromosomes
mainChrSeqSet <- DNAStringSet(mainChrSeq)
writeXStringSet(mainChrSeqSet,
                "BSgenome.Hsapiens.UCSC.hg19.mainChrs.fa")
#Creando indice Rsubread previo a alineación
buildindex("BSgenome.Hsapiens.UCSC.hg19.mainChrs",
           "BSgenome.Hsapiens.UCSC.hg19.mainChrs.fa",
            indexSplit = TRUE,
            memory = 10000)

#Creamos grupo de lecuras para cada grupo experimental
#MCF7E2
MCF7E2read1<-"C:/Users/Kokia/OneDrive/Documentos/Bioquímica clases maes/Leon-ATAC/ATAC-SEQ/18134-118-02_S103_L002_R1_001.fastq"
MCF7E2read2<-"C:/Users/Kokia/OneDrive/Documentos/Bioquímica clases maes/Leon-ATAC/ATAC-SEQ/18134-118-02_S103_L002_R2_001.fastq"
#MCF7Tam
MCF7Tamread1<- "C:/Users/Kokia/OneDrive/Documentos/Bioquímica clases maes/Leon-ATAC/ATAC-SEQ/18134-118-03_S0_L001_R1_001.fastq"
MCF7Tamread2<-"C:/Users/Kokia/OneDrive/Documentos/Bioquímica clases maes/Leon-ATAC/ATAC-SEQ/18134-118-03_S0_L001_R2_001.fastq"
#MCF7Veh
MCF7Vehread1<-"C:/Users/Kokia/OneDrive/Documentos/Bioquímica clases maes/Leon-ATAC/ATAC-SEQ/18134-118-01_S102_L002_R1_001.fastq"
MCF7Vehread2<-"C:/Users/Kokia/OneDrive/Documentos/Bioquímica clases maes/Leon-ATAC/ATAC-SEQ/18134-118-01_S102_L002_R2_001.fastq"
#MNE2
MNE2read1<-"C:/Users/Kokia/OneDrive/Documentos/Bioquímica clases maes/Leon-ATAC/ATAC-SEQ/18134-118-05_S105_L002_R1_001.fastq"
MNE2read2<-"C:/Users/Kokia/OneDrive/Documentos/Bioquímica clases maes/Leon-ATAC/ATAC-SEQ/18134-118-05_S105_L002_R2_001.fastq"
#MNTam
MNTamread1<-"C:/Users/Kokia/OneDrive/Documentos/Bioquímica clases maes/Leon-ATAC/ATAC-SEQ/18134-118-06_S106_L002_R1_001.fastq"
MNTamread2<-"C:/Users/Kokia/OneDrive/Documentos/Bioquímica clases maes/Leon-ATAC/ATAC-SEQ/18134-118-06_S106_L002_R2_001.fastq"

#Alineamos nuestros fastq con el indice previamente creado
#MCF7E2
align("BSgenome.Hsapiens.UCSC.hg19.mainChrs",
      readfile1=MCF7E2read1,readfile2=MCF7E2read2,
      output_file = "ATAC_50K_2.bam",
      nthreads=2,type=1,
      unique=TRUE,maxFragLength = 2000)
#MCF7Tam
align("BSgenome.Hsapiens.UCSC.hg19.mainChrs",
      readfile1=MCF7Tamread1,readfile2=MCF7Tamread2,
      output_file = "ATAC_50K_2.bam",
      nthreads=2,type=1,
      unique=TRUE,maxFragLength = 2000)
#MCF7Veh
align("BSgenome.Hsapiens.UCSC.hg19.mainChrs",
      readfile1=MCF7Vehread1,readfile2=MCF7Vehread2,
      output_file = "ATAC_50K_2.bam",
      nthreads=2,type=1,
      unique=TRUE,maxFragLength = 2000)
#MNE2
align("BSgenome.Hsapiens.UCSC.hg19.mainChrs",
      readfile1=MNE2read1,readfile2=MNE2read2,
      output_file = "ATAC_50K_2.bam",
      nthreads=2,type=1,
      unique=TRUE,maxFragLength = 2000)
#MNTam
align("BSgenome.Hsapiens.UCSC.hg19.mainChrs",
      readfile1=MNTamread1,readfile2=MNTamread2,
      output_file = "ATAC_50K_2.bam",
      nthreads=2,type=1,
      unique=TRUE,maxFragLength = 2000)
#Una vez que hemos creado los bam para cada grupo procedemos a ahcer sorting e
#indexarlos para su uso en otras librerias o programas como IGV
library(Rsamtools)
#MCF7E2
outBAM1<-"ATAC_MCF7_E2.bam"
sortedBAM1 <- file.path(dirname(outBAM1),
                       paste0("Sorted_",basename(outBAM1))
)
sortBam(outBAM1,gsub("\\.bam","",basename(sortedBAM1)))
indexBam(sortedBAM1)
#MCF7Tam
outBAM2<-"ATAC_MCF7_Tam.bam"
sortedBAM2 <- file.path(dirname(outBAM2),
                       paste0("Sorted_",basename(outBAM2))
)
sortBam(outBAM2,gsub("\\.bam","",basename(sortedBAM2)))
indexBam(sortedBAM2)
#MCF7Veh
outBAM3<-"ATAC_MCF7_Veh.bam"
sortedBAM3 <- file.path(dirname(outBAM3),
                        paste0("Sorted_",basename(outBAM3))
)
sortBam(outBAM3,gsub("\\.bam","",basename(sortedBAM3)))
indexBam(sortedBAM3)
#MNE2
outBAM4<-"ATAC_MN_E2.bam"
sortedBAM4 <- file.path(dirname(outBAM4),
                        paste0("Sorted_",basename(outBAM4))
)
sortBam(outBAM4,gsub("\\.bam","",basename(sortedBAM4)))
indexBam(sortedBAM4)
#MNTam
outBAM5<- "ATAC_MN_Tam.bam"
sortedBAM5 <- file.path(dirname(outBAM5),
                        paste0("Sorted_",basename(outBAM5))
)
sortBam(outBAM5,gsub("\\.bam","",basename(sortedBAM5)))
indexBam(sortedBAM5)
#Reasignamos identidades
Sort_MCF7_E2 <- sortedBAM1
Sort_MCF7_Tam <- sortedBAM2
Sort_MCF7_Veh <- sortedBAM3
Sort_MN_E2 <- sortedBAM4
Sort_MN_Tam <- sortedBAM5
#Checar cantidad de lecturas mapeadas a cromosomas correspondientes
mappedreads1<- idxstatsBam(Sort_MCF7_E2)
mappedreads2<- idxstatsBam(Sort_MCF7_Tam)
mappedreads3<- idxstatsBam(Sort_MCF7_Veh)
mappedreads4<- idxstatsBam(Sort_MN_E2)
mappedreads5<- idxstatsBam(Sort_MN_Tam)
#Ploteamos para conocer la distribucion
library(ggplot2)
ggplot(mappedreads1, aes(seqnames, mapped, fill = seqnames)) + geom_bar(stat = "identity") +
   coord_flip() + ggtitle("MCF7 E2")

ggplot(mappedreads2, aes(seqnames, mapped, fill = seqnames)) + geom_bar(stat = "identity") +
   coord_flip() + ggtitle("MCF7 Tam")

ggplot(mappedreads3, aes(seqnames, mapped, fill = seqnames)) + geom_bar(stat = "identity") +
   coord_flip() + ggtitle("MCF7 Veh")

ggplot(mappedreads4, aes(seqnames, mapped, fill = seqnames)) + geom_bar(stat = "identity") +
   coord_flip() + ggtitle("MN E2")

ggplot(mappedreads5, aes(seqnames, mapped, fill = seqnames)) + geom_bar(stat = "identity") +
   coord_flip() + ggtitle("MN Tam")
#Procesamiento post alineamiento
library(GenomicAlignments)
flags = scanBamFlag(isProperPair = TRUE)
myParam = ScanBamParam(flag = flags, what = c("qname", "mapq", "isize"))

#Creamos el subset de lecturas unicas y pareadas en nuestros datos
atacReadsMCF7E2 <- readGAlignmentPairs(sortedBAM1, param = myParam)
atacReadsMCF7Tam <- readGAlignmentPairs(sortedBAM2, param = myParam)
atacReadsMCF7Veh <- readGAlignmentPairs(sortedBAM3, param = myParam)
atacReadsMNE2 <- readGAlignmentPairs(sortedBAM4, param = myParam)
atacReadsMNTam <- readGAlignmentPairs(sortedBAM5, param = myParam)

#Definimos el tamaños de nuestros insetos en los distintos grupos
#MCF7E2
read1_MCF7E2<-first(atacReadsMCF7E2)
read2_MCF7E2<- second(atacReadsMCF7E2)
atacReadsMCF7E2_read1 <- first(atacReadsMCF7E2)
insertSizes_MCF7E2 <- abs(elementMetadata(atacReadsMCF7E2_read1)$isize)

#MCF7Tam
read1_MCF7Tam<-first(atacReadsMCF7Tam)
read2_MCF7Tam<- second(atacReadsMCF7Tam)
atacReadsMCF7Tam_read1 <- first(atacReadsMCF7Tam)
insertSizes_MCF7Tam <- abs(elementMetadata(atacReadsMCF7Tam_read1)$isize)

#MCF7Veh
read1_MCF7Veh<-first(atacReadsMCF7Veh)
read2_MCF7Veh<- second(atacReadsMCF7Veh)
atacReadsMCF7Veh_read1 <- first(atacReadsMCF7Veh)
insertSizes_MCF7Veh <- abs(elementMetadata(atacReadsMCF7Veh_read1)$isize)

#MNE2
read1_MNE2<-first(atacReadsMNE2)
read2_MNE2<- second(atacReadsMNE2)
atacReadsMNE2_read1 <- first(atacReadsMNE2)
insertSizes_MNE2 <- abs(elementMetadata(atacReadsMNE2_read1)$isize)

#MNTam
read1_MNTam<-first(atacReadsMNTam)
read2_MNTam<- second(atacReadsMNTam)
atacReadsMNTam_read1 <- first(atacReadsMNTam)
insertSizes_MNTam <- abs(elementMetadata(atacReadsMNTam_read1)$isize)

#El paso posterior es crear una serie de variables que contengan información de
#interes como podría ser zonas libres de nucelosomas, mono ocupados y diocupados.

#MCF7 E2
atacReads_NucFree_MCF7E2 <- atacReadsMCF7E2[insertSizes_MCF7E2 < 100, ]
atacReads_MonoNuc_MCF7E2 <- atacReadsMCF7E2[insertSizes_MCF7E2 > 180 & insertSizes_MCF7E2 < 240, ]
atacReads_diNuc_MCF7E2 <- atacReadsMCF7E2[insertSizes_MCF7E2 > 315 & insertSizes_MCF7E2 < 437, ]

#MCF7 Tam
atacReads_NucFree_MCF7Tam <- atacReadsMCF7Tam[insertSizes_MCF7Tam < 100, ]
atacReads_MonoNuc_MCF7Tam <- atacReadsMCF7Tam[insertSizes_MCF7Tam > 180 & insertSizes_MCF7Tam < 240, ]
atacReads_diNuc_MCF7Tam <- atacReadsMCF7Tam[insertSizes_MCF7Tam > 315 & insertSizes_MCF7Tam < 437, ]

#MCF7 Veh
atacReads_NucFree_MCF7Veh <- atacReadsMCF7Veh[insertSizes_MCF7Veh < 100, ]
atacReads_MonoNuc_MCF7Veh <- atacReadsMCF7Veh[insertSizes_MCF7Veh > 180 & insertSizes_MCF7Veh < 240, ]
atacReads_diNuc_MCF7Veh <- atacReadsMCF7Veh[insertSizes_MCF7Veh > 315 & insertSizes_MCF7Veh < 437, ]

#MN E2
atacReads_NucFree_MNE2 <- atacReadsMNE2[insertSizes_MNE2 < 100, ]
atacReads_MonoNuc_MNE2 <- atacReadsMNE2[insertSizes_MNE2 > 180 & insertSizes_MNE2 < 240, ]
atacReads_diNuc_MNE2 <- atacReadsMNE2[insertSizes_MNE2 > 315 & insertSizes_MNE2 < 437, ]

#MN Tam
atacReads_NucFree_MNTam <- atacReadsMNTam[insertSizes_MNTam < 100, ]
atacReads_MonoNuc_MNTam <- atacReadsMNTam[insertSizes_MNTam > 180 & insertSizes_MNTam < 240, ]
atacReads_diNuc_MNTam <- atacReadsMNTam[insertSizes_MNTam > 315 & insertSizes_MNTam < 437, ]

#Exportamos los archivos a format BAM para poder realizar visualizaciones en otros softwares
#como IGV

#MCF7 E2
nucFreeRegionBam_MCF7E2 <- gsub("\\.bam", "_nucFreeRegions\\.bam", sortedBAM1)
monoNucBam_MCF7E2 <- gsub("\\.bam", "_monoNuc\\.bam", sortedBAM1)
diNucBam_MCF7E2 <- gsub("\\.bam", "_diNuc\\.bam", sortedBAM1)
library(rtracklayer)
export(atacReads_NucFree_MCF7E2, nucFreeRegionBam_MCF7E2, format = "bam")
export(atacReads_MonoNuc_MCF7E2, monoNucBam_MCF7E2, format = "bam")
export(atacReads_diNuc_MCF7E2, diNucBam_MCF7E2, format = "bam")

#MCF7 Tam
nucFreeRegionBam_MCF7Tam <- gsub("\\.bam", "_nucFreeRegions\\.bam", sortedBAM2)
monoNucBam_MCF7Tam <- gsub("\\.bam", "_monoNuc\\.bam", sortedBAM2)
diNucBam_MCF7Tam <- gsub("\\.bam", "_diNuc\\.bam", sortedBAM2)
library(rtracklayer)
export(atacReads_NucFree_MCF7Tam, nucFreeRegionBam_MCF7Tam, format = "bam")
export(atacReads_MonoNuc_MCF7Tam, monoNucBam_MCF7Tam, format = "bam")
export(atacReads_diNuc_MCF7Tam, diNucBam_MCF7Tam, format = "bam")

#MCF7 Veh
nucFreeRegionBam_MCF7Veh <- gsub("\\.bam", "_nucFreeRegions\\.bam", sortedBAM3)
monoNucBam_MCF7Veh <- gsub("\\.bam", "_monoNuc\\.bam", sortedBAM3)
diNucBam_MCF7Veh <- gsub("\\.bam", "_diNuc\\.bam", sortedBAM3)
library(rtracklayer)
export(atacReads_NucFree_MCF7Veh, nucFreeRegionBam_MCF7Veh, format = "bam")
export(atacReads_MonoNuc_MCF7Veh, monoNucBam_MCF7Veh, format = "bam")
export(atacReads_diNuc_MCF7Veh, diNucBam_MCF7Veh, format = "bam")

#MN E2
nucFreeRegionBam_MNE2 <- gsub("\\.bam", "_nucFreeRegions\\.bam", sortedBAM4)
monoNucBam_MNE2 <- gsub("\\.bam", "_monoNuc\\.bam", sortedBAM4)
diNucBam_MNE2 <- gsub("\\.bam", "_diNuc\\.bam", sortedBAM4)
library(rtracklayer)
export(atacReads_NucFree_MNE2, nucFreeRegionBam_MNE2, format = "bam")
export(atacReads_MonoNuc_MNE2, monoNucBam_MNE2, format = "bam")
export(atacReads_diNuc_MNE2, diNucBam_MNE2, format = "bam")

#MN Tam
nucFreeRegionBam_MNTam <- gsub("\\.bam", "_nucFreeRegions\\.bam", sortedBAM5)
monoNucBam_MNTam <- gsub("\\.bam", "_monoNuc\\.bam", sortedBAM5)
diNucBam_MNTam <- gsub("\\.bam", "_diNuc\\.bam", sortedBAM5)
library(rtracklayer)
export(atacReads_NucFree_MNTam, nucFreeRegionBam_MNTam, format = "bam")
export(atacReads_MonoNuc_MNTam, monoNucBam_MNTam, format = "bam")
export(atacReads_diNuc_MNTam, diNucBam_MNTam, format = "bam")

#Creamos el parametro atac fragments para visualizar en los BigWig
atacFragments_MCF7E2 <- granges(atacReadsMCF7E2)
atacFragments_MCF7Tam <- granges(atacReadsMCF7Tam)
atacFragments_MCF7Veh <- granges(atacReadsMCF7Veh)
atacFragments_MNE2 <- granges(atacReadsMNE2)
atacFragments_MNTam <- granges(atacReadsMNTam)

#Creamos los archivos Bigwig para visualizar en IGV
#MCF7 E2
openRegionRPMBigWig_MCF7E2 <- gsub("\\.bam", "_openRegionRPM\\.bw", sortedBAM1)
myCoverage_MCF7E2 <- coverage(atacFragments_MCF7E2, weight = (10^6/length(atacFragments_MCF7E2)))
export.bw(myCoverage_MCF7E2, openRegionRPMBigWig_MCF7E2)

#MCF7 Tam
openRegionRPMBigWig_MCF7Tam <- gsub("\\.bam", "_openRegionRPM\\.bw", sortedBAM2)
myCoverage_MCF7Tam <- coverage(atacFragments_MCF7Tam, weight = (10^6/length(atacFragments_MCF7Tam)))
export.bw(myCoverage_MCF7Tam, openRegionRPMBigWig_MCF7Tam)

#MCF7 Veh
openRegionRPMBigWig_MCF7Veh <- gsub("\\.bam", "_openRegionRPM\\.bw", sortedBAM3)
myCoverage_MCF7Veh <- coverage(atacFragments_MCF7Veh, weight = (10^6/length(atacFragments_MCF7Veh)))
export.bw(myCoverage_MCF7Veh, openRegionRPMBigWig_MCF7Veh)

#MNE2
openRegionRPMBigWig_MNE2 <- gsub("\\.bam", "_openRegionRPM\\.bw", sortedBAM4)
myCoverage_MNE2 <- coverage(atacFragments_MNE2, weight = (10^6/length(atacFragments_MNE2)))
export.bw(myCoverage_MNE2, openRegionRPMBigWig_MNE2)

#MNTam
openRegionRPMBigWig_MNTam <- gsub("\\.bam", "_openRegionRPM\\.bw", sortedBAM5)
myCoverage_MNTam <- coverage(atacFragments_MNTam, weight = (10^6/length(atacFragments_MNTam)))
export.bw(myCoverage_MNTam, openRegionRPMBigWig_MNTam)

#La parte siguiente es el peakcalling, uno de los peakcallers más populares es MACS2
#por lo que el código siguiente se debe introducir en macs2 en una terminal linux
#sustituyendo el nombre del archvo y el directorio de salida por el que les convenga
macspeaks <- dir("./Diferencial", pattern = "*.narrowPeak", full.names = TRUE)

#La elección entre narrowpeaks y broadpeaks depende de que es lo que su busca,
#el formato narrowpeaks es ideal para encontrar TF o marcas de histonas muy especificas
#por ejemplo

#Calidad y limpieza de peaks
BiocManager::install(c("ChIPQC","rtracklayer", "DT","dplyr","tidiyr"))
library(ChIPQC)
library(rtracklayer)
library(DT)
library(dplyr)
library(tidyr)
#Importamos el bed de las regiones a excluir
blkList<- import.bed("./consensusBlacklist.bed")
#MCF7E2
openRegionPeaksMCF7E2<- "./Sorted_ATAC_MCF7_E2_peaks.narrowPeak"
qcResMCF7E2 <- ChIPQCsample("C:/Users/Kokia/OneDrive/Documentos/Bioquímica clases maes/Leon-ATAC/ATAC-SEQ/Sorted_ATAC_MCF7_E2.bam",
                            peaks = openRegionPeaksMCF7E2, annotation = "hg19", blacklist = blkList, verboseT = FALSE)
#Removemos las secuencias blaklisted de nuestros archivos MCF7E2
MacsCallsMCF7_e2 <- granges(qcResMCF7E2[seqnames(qcResMCF7E2) %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7",
                                                                     "chr8","chr9","chr10","chr11","chr12","chr13","chr14",
                                                                     "chr15","chr16","chr17","chr18","chr19","chr20","chr21",
                                                                     "chrX","chrY","chrM")])
data.frame(Blacklisted = sum(MacsCallsMCF7_e2 %over% blkList), Not_Blacklisted = sum(!MacsCallsMCF7_e2 %over% blkList))
MacsCallsMCF7_e2<- MacsCallsMCF7_e2[!MacsCallsMCF7_e2 %over% blkList]

#MCF7 Tam

openRegionPeaksMCF7Tam<- "./Sorted_ATAC_MCF7_Tam_peaks.narrowPeak"
qcResMCF7Tam <- ChIPQCsample("C:/Users/Kokia/OneDrive/Documentos/Bioquímica clases maes/Leon-ATAC/ATAC-SEQ/Sorted_ATAC_MCF7_Tam.bam",
                            peaks = openRegionPeaksMCF7Tam, annotation = "hg19", blacklist = blkList, verboseT = FALSE)
#Removemos las secuencias blaklisted de nuestros archivos MCF7E2
MacsCallsMCF7_tam <- granges(qcResMCF7Tam[seqnames(qcResMCF7Tam) %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7",
                                                                     "chr8","chr9","chr10","chr11","chr12","chr13","chr14",
                                                                     "chr15","chr16","chr17","chr18","chr19","chr20","chr21",
                                                                     "chrX","chrY","chrM")])

data.frame(Blacklisted = sum(MacsCallsMCF7_tam %over% blkList), Not_Blacklisted = sum(!MacsCallsMCF7_tam %over% blkList))
MacsCallsMCF7_tam<- MacsCallsMCF7_tam[!MacsCallsMCF7_tam %over% blkList]

#MCF7 Veh
openRegionPeaksMCF7Veh<- "./Sorted_ATAC_MCF7_Veh_peaks.narrowPeak"
qcResMCF7Veh <- ChIPQCsample("C:/Users/Kokia/OneDrive/Documentos/Bioquímica clases maes/Leon-ATAC/ATAC-SEQ/Sorted_ATAC_MCF7_Veh.bam",
                            peaks = openRegionPeaksMCF7Veh, annotation = "hg19", blacklist = blkList, verboseT = FALSE)
#Removemos las secuencias blaklisted de nuestros archivos MCF7E2
MacsCallsMCF7_Veh <- granges(qcResMCF7Veh[seqnames(qcResMCF7Veh) %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7",
                                                                     "chr8","chr9","chr10","chr11","chr12","chr13","chr14",
                                                                     "chr15","chr16","chr17","chr18","chr19","chr20","chr21",
                                                                     "chrX","chrY","chrM")])
data.frame(Blacklisted = sum(MacsCallsMCF7_Veh %over% blkList), Not_Blacklisted = sum(!MacsCallsMCF7_Veh %over% blkList))
MacsCallsMCF7_Veh<- MacsCallsMCF7_Veh[!MacsCallsMCF7_Veh %over% blkList]

#MN E2

openRegionPeaksMNE2<- "./Sorted_ATAC_MN_E2_peaks.narrowPeak"
qcResMNE2 <- ChIPQCsample("C:/Users/Kokia/OneDrive/Documentos/Bioquímica clases maes/Leon-ATAC/ATAC-SEQ/Sorted_ATAC_MN_E2.bam",
                             peaks = openRegionPeaksMNE2, annotation = "hg19", blacklist = blkList, verboseT = FALSE)
#Removemos las secuencias blaklisted de nuestros archivos MN E2
MacsCallsMNE2<- granges(qcResMNE2[seqnames(qcResMNE2) %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7",
                                                                        "chr8","chr9","chr10","chr11","chr12","chr13","chr14",
                                                                        "chr15","chr16","chr17","chr18","chr19","chr20","chr21",
                                                                        "chrX","chrY","chrM")])
data.frame(Blacklisted = sum(MacsCallsMNE2 %over% blkList), Not_Blacklisted = sum(!MacsCallsMNE2 %over% blkList))
MacsCallsMNE2<- MacsCallsMNE2[!MacsCallsMNE2 %over% blkList]

#MN Tam

openRegionPeaksMN_Tam<- "./Sorted_ATAC_MN_Tam_peaks.narrowPeak"
qcResMN_Tam <- ChIPQCsample("C:/Users/Kokia/OneDrive/Documentos/Bioquímica clases maes/Leon-ATAC/ATAC-SEQ/Sorted_ATAC_MN_Tam.bam",
                          peaks = openRegionPeaksMN_Tam, annotation = "hg19", blacklist = blkList, verboseT = FALSE)
#Removemos las secuencias blaklisted de nuestros archivos MN Tam
MacsCallsMN_Tam<- granges(qcResMN_Tam[seqnames(qcResMN_Tam) %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7",
                                                             "chr8","chr9","chr10","chr11","chr12","chr13","chr14",
                                                             "chr15","chr16","chr17","chr18","chr19","chr20","chr21",
                                                             "chrX","chrY","chrM")])
data.frame(Blacklisted = sum(MacsCallsMN_Tam %over% blkList), Not_Blacklisted = sum(!MacsCallsMN_Tam %over% blkList))
MacsCallsMN_Tam<- MacsCallsMN_Tam[!MacsCallsMN_Tam %over% blkList]

#Anotación de los picos encontrados
#MCF7 E2
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
MacsCallsMCF7_e2_Anno<- annotatePeak(MacsCallsMCF7_e2, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)
plotAnnoPie(MacsCallsMCF7_e2_Anno)
#Análisis funcional de los picos
library(rGREAT)
great_Job_MCF7_E2 <- submitGreatJob(MacsCallsMCF7_e2, species = "hg19")
great_ResultTable_MCF7_E2 = getEnrichmentTables(great_Job_MCF7_E2, category = "GO")
great_ResultTable_MCF7_E2[["GO Biological Process"]][1:30, ]

#MCF7 Tam
MacsCallsMCF7_Tam_Anno<- annotatePeak(MacsCallsMCF7_tam, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)
plotAnnoPie(MacsCallsMCF7_Tam_Anno)
#Análisis funcional de los picos
library(rGREAT)
great_Job_MCF7_Tam <- submitGreatJob(MacsCallsMCF7_tam, species = "hg19")
great_ResultTable_MCF7_Tam = getEnrichmentTables(great_Job_MCF7_Tam, category = "GO")
great_ResultTable_MCF7_Tam[["GO Biological Process"]][1:40, ]

#MCF7 Veh
MacsCallsMCF7_Veh_Anno<- annotatePeak(MacsCallsMCF7_Veh, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)
plotAnnoPie(MacsCallsMCF7_Veh_Anno)
#Análisis funcional de los picos
library(rGREAT)
great_Job_MCF7_Veh <- submitGreatJob(MacsCallsMCF7_Veh, species = "hg19")
great_ResultTable_MCF7_Veh = getEnrichmentTables(great_Job_MCF7_Veh, category = "GO")
great_ResultTable_MCF7_Veh[["GO Biological Process"]][1:40, ]

#MN E2
MacsCallsMNE2_Anno<- annotatePeak(MacsCallsMNE2, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)
plotAnnoPie(MacsCallsMNE2_Anno)
#Análisis funcional de los picos
library(rGREAT)
great_Job_MNE2 <- submitGreatJob(MacsCallsMNE2, species = "hg19")
great_ResultTable_MNE2 = getEnrichmentTables(great_Job_MNE2, category = "GO")
great_ResultTable_MNE2[["GO Biological Process"]][1:40, ]

#MN Tam
MacsCallsMN_Tam_Anno<- annotatePeak(MacsCallsMN_Tam, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)
plotAnnoPie(MacsCallsMN_Tam_Anno)
#Análisis funcional de los picos
library(rGREAT)
great_Job_MN_Tam <- submitGreatJob(MacsCallsMN_Tam, species = "hg19")
great_ResultTable_MN_Tam = getEnrichmentTables(great_Job_MN_Tam, category = "GO")
great_ResultTable_MN_Tam[["GO Biological Process"]][1:40, ]


#Análisis diferencial de datos
#Hacemos un set que incuya todos los grupo de replicas y datos de ATAC Seq para
#buscar sets de picos no redundantes

peaks <- dir("./Diferencial", pattern = "*.narrowPeak", full.names = TRUE)
myPeaks <- lapply(peaks, ChIPQC:::GetGRanges, simple = TRUE)
myPeaks
allPeaksSet_nR<- reduce(unlist(GRangesList(myPeaks)))
overlap <- list()
for (i in 1:length(myPeaks)) {
   overlap[[i]] <- allPeaksSet_nR %over% myPeaks[[i]]
}
overlapMatrix <- do.call(cbind, overlap)
colnames(overlapMatrix) <- basename(peaks)
mcols(allPeaksSet_nR) <- overlapMatrix
#Hacemos blacklisting de los datos
nrToCount <- allPeaksSet_nR[!allPeaksSet_nR %over% blkList & !seqnames(allPeaksSet_nR) %in%
                               "chrM"]
nrToCount
#Conteo de los picos diferenciales
library(Rsubread)
occurrences <- rowSums(as.data.frame(elementMetadata(nrToCount)))
nrToCount <- nrToCount[occurrences >= 2, ]
nrToCount
library(GenomicAlignments)
bamsToCount <- dir("./Diferencial", full.names = TRUE,
                   pattern = "*.\\.bam$")
consensusToCount <- soGGi:::runConsensusRegions(GRangesList(myPeaks), "none")
myCounts <- summarizeOverlaps(nrToCount, bamsToCount, singleEnd = FALSE)
myCounts
#Iniciamos el análisis diferencial
library(DESeq2)
colnames(myCounts) <- c("MCF7_E2_1","MCF7_E2_2", "MCF7_Tam_1","MCF7_Tam_2","MN_E2_1", "MN_E2_2","MN_Tam_1", "MN_Tam_2")
library(DESeq2)
Group <- factor(c("MCF7_E2","MCF7_E2", "MCF7_Tam","MCF7_Tam","MN_E2", "MN_E2","MN_Tam", "MN_Tam"))
metaData <- data.frame(Group, row.names = colnames(myCounts))
atacDDS <- DESeqDataSetFromMatrix(assay(myCounts), metaData, ~Group, rowRanges = rowRanges(myCounts))
atacDDS <- DESeq(atacDDS)

#Ya que creamos nuestra variable atacDDS podemos proceder a hacer las comparaciones entre los grupos
#
#MCF7 E2 vs MN E2
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(writexl)
MCF7_E2_vs_MN_E2 <- results(atacDDS, c("Group", "MCF7_E2", "MN_E2"), format = "GRanges")
MCF7_E2_vs_MN_E2 <- MCF7_E2_vs_MN_E2[order(MCF7_E2_vs_MN_E2$pvalue)]
MCF7_E2_vs_MN_E2
anno_MCF7vsMN_E2 <- annotatePeak(MCF7_E2_vs_MN_E2, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,verbose = FALSE)
MCF7vsMN_E2_ATAC <- as.data.frame(anno_MCF7vsMN_E2)
MCF7vsMN_E2_ATAC[1, ]
write_xlsx(MCF7vsMN_E2_ATAC, "./Diferencial/MCF7_MN_E2.xlsx")
library(clusterProfiler)
MCF7vsMN_E2_GO <- enrichGO(MCF7vsMN_E2_ATAC$geneId, OrgDb = "org.Hs.eg.db", ont = "BP", maxGSSize = 5000)
MCF7vsMN_E2_GO_DF <- as.data.frame(MCF7vsMN_E2_GO)
write_xlsx(MCF7vsMN_E2_GO_DF, "./Diferencial/MCF7_MN_E2_GO.xlsx")

#MCF7 Tam vs MN Tam
MCF7_Tam_vs_MN_Tam <- results(atacDDS, c("Group", "MCF7_Tam", "MN_Tam"), format = "GRanges")
MCF7_Tam_vs_MN_Tam <- MCF7_Tam_vs_MN_Tam[order(MCF7_Tam_vs_MN_Tam$pvalue)]
MCF7_Tam_vs_MN_Tam
anno_MCF7vsMN_Tam <- annotatePeak(MCF7_Tam_vs_MN_Tam, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,verbose = FALSE)
MCF7vsMN_Tam_ATAC <- as.data.frame(anno_MCF7vsMN_Tam)
MCF7vsMN_Tam_ATAC[1, ]
write_xlsx(MCF7vsMN_Tam_ATAC, "./Diferencial/MCF7_MN_Tam.xlsx")
library(clusterProfiler)
MCF7vsMN_Tam_GO <- enrichGO(MCF7vsMN_Tam_ATAC$geneId, OrgDb = "org.Hs.eg.db", ont = "BP", maxGSSize = 500)
MCF7vsMN_Tam_GO_DF <- as.data.frame(MCF7vsMN_Tam_GO)
write_xlsx(MCF7vsMN_Tam_GO_DF, "./Diferencial/MCF7_MN_Tam_GO_2.xlsx")
