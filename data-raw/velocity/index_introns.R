library(BUSpaRse)
library(AnnotationHub)
library(BSgenome.Hsapiens.UCSC.hg38)

# query AnnotationHub for human Ensembl annotation
ah <- AnnotationHub()
query(ah, pattern = c("Ensembl", "94", "Homo sapiens", "EnsDb"))
edb <- ah[["AH64923"]]

# for 10x v1/v2 L=98
# for 10x v3    L=91
# make sure annotation and  genome use same version (e.g. GRCh38)
out_path <- 'velocity/EnsDb.Hsapiens.v94.91nt'
dir.create(out_path, recursive = TRUE)

get_velocity_files(edb,
                   L = 91,
                   Genome = BSgenome.Hsapiens.UCSC.hg38,
                   out_path = out_path,
                   isoform_action = "separate")

# Intron index
# cd data-raw/velocity/EnsDb.Hsapiens.v94
# kallisto index -i ./cDNA_introns.idx ./cDNA_introns.fa
