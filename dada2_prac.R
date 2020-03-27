library(dada2)

path <-"110_16S/"
list.files(path)

#import files
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#QC
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

#Filter+trim
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,225), trimLeft = c(17,21),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)
head(out)

#Learn error rate
errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)

plotErrors(errF, nominalQ=TRUE)

#sample inference
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)
dadaFs[[1]]

# merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=FALSE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#remove chimers
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=FALSE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)


#track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "tax/silva_nr_v132_train_set.fa.gz", multithread=FALSE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#phyloseq handoff