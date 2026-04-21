# ============================================================
# 00_prepare_alignment.R
# Read FASTA, align, export files, save objects
#
# Terminal: cd into R/ folder, then:
#   Rscript 00_prepare_alignment.R
#
# Outputs:
#   out/rdata/alignment_objects.RData
#   FASTA/gulf_fish_COI_aligned.fasta
#   MrBayes/gulf_fish_COI_aligned.nex
# ============================================================

# install.packages(c("ape", "phangorn"))
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("msa")
# BiocManager::install("Biostrings")

library(ape)
library(phangorn)
library(msa)

dir.create("out/rdata", recursive = TRUE, showWarnings = FALSE)
dir.create("out/trees", recursive = TRUE, showWarnings = FALSE)


cat("\n--- STEP 1: Reading FASTA ---\n")

dna_strings <- Biostrings::readDNAStringSet(
  "FASTA/gulf_fish_COI_combined.fasta"
)

cat("  Sequences read:", length(dna_strings), "\n")
cat("  Tip labels:\n")
print(names(dna_strings))


cat("\n--- STEP 2: Running ClustalW alignment ---\n")

msa_result <- msa(dna_strings, method = "ClustalW", type = "dna")

aln_length <- ncol(as.matrix(msa_result))
cat("  Alignment length:", aln_length, "bp\n")
cat("  Number of taxa:  ", nrow(as.matrix(msa_result)), "\n")


cat("\n--- STEP 3: Converting alignment formats ---\n")

aln_dnabin <- msaConvert(msa_result, type = "ape::DNAbin")
aln_phydat <- msaConvert(msa_result, type = "phangorn::phyDat")

cat("  DNAbin dimensions: ", nrow(aln_dnabin), "taxa x",
    ncol(aln_dnabin), "sites\n")


family_colors <- c(
  "Red_Drum"          = "#2166ac",
  "Spec_Trout"        = "#2166ac",
  "Black_Drum"        = "#2166ac",
  "Atlantic_Croaker"  = "#2166ac",
  "Southern_Flounder" = "#d95f02",
  "Sheepshead"        = "#1b9e77",
  "Flathead_Mullet"   = "#888888",
  "Bull_Shark"        = "#e41a1c",
  "Rainbow_Trout"     = "#FA8072",
  "Brown_Trout"       = "#FA8072"
)

tip_color <- function(labels) {
  species <- sub("_[0-9]+$", "", labels)
  unname(family_colors[species])
}


cat("\n--- STEP 4: Saving out/rdata/alignment_objects.RData ---\n")

save(aln_dnabin, aln_phydat, family_colors, tip_color,
     file = "out/rdata/alignment_objects.RData")

cat("  Saved: out/rdata/alignment_objects.RData\n")


cat("\n--- STEP 5: Writing alignment files ---\n")

write.dna(
  aln_dnabin,
  file   = "FASTA/gulf_fish_COI_aligned.fasta",
  format = "fasta",
  colsep = ""
)
cat("  Saved: FASTA/gulf_fish_COI_aligned.fasta\n")

aln_mat  <- as.character(aln_dnabin)
aln_list <- setNames(
  lapply(seq_len(nrow(aln_mat)), function(i) aln_mat[i, ]),
  rownames(aln_mat)
)

write.nexus.data(
  aln_list,
  file        = "MrBayes/gulf_fish_COI_aligned.nex",
  format      = "dna",
  interleaved = FALSE
)
cat("  Saved: MrBayes/gulf_fish_COI_aligned.nex\n")


cat("\n============================================================\n")
cat("ALIGNMENT SUMMARY\n")
cat("============================================================\n")
cat("  Input sequences : ", length(dna_strings), "\n")
cat("  Aligned length  : ", aln_length, "bp\n")
cat("  Taxa            : ", nrow(aln_dnabin), "\n")
cat("  Sites (aligned) : ", ncol(aln_dnabin), "\n")

aln_chars     <- as.character(aln_dnabin)
variable_sites <- sum(apply(aln_chars, 2, function(col) {
  length(unique(col[col != "-" & col != "n"])) > 1
}))
cat("  Variable sites  : ", variable_sites, "\n")
cat("============================================================\n")
