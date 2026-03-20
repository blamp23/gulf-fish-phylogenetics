# ============================================================
# 00_prepare_alignment.R
# Shared setup: read FASTA, align, export files, save objects
# Author: [Your Name]  |  Texas A&M University
#
# Terminal: cd into R/ folder, then:
#   Rscript 00_prepare_alignment.R
#
# Outputs:
#   out/rdata/alignment_objects.RData  <- loaded by all three tree scripts
#   FASTA/gulf_fish_COI_aligned.fasta
#   MrBayes/gulf_fish_COI_aligned.nex  <- used by 03_tree_Bayesian.R
# ============================================================

# -- Install packages (run once, then leave commented out) ---
# install.packages(c("ape", "phangorn"))
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("msa")
# BiocManager::install("Biostrings")
# ------------------------------------------------------------


# ============================================================
# STEP 1: Load libraries
# WHY: Each library handles a different piece of the pipeline.
#   ape       — reading/writing sequences, DNAbin format, plotting
#   phangorn  — phyDat format (needed by ML and model testing)
#   msa       — MUSCLE multiple sequence alignment (Bioconductor)
#   Biostrings — DNAStringSet format (input to msa)
# ============================================================
library(ape)
library(phangorn)
library(msa)

# Create output directories if they don't exist
dir.create("out/rdata", recursive = TRUE, showWarnings = FALSE)
dir.create("out/trees", recursive = TRUE, showWarnings = FALSE)


# ============================================================
# STEP 2: Read unaligned sequences
# WHY: readDNAStringSet() reads directly into Biostrings format,
# which is what msa() expects. The FASTA already has clean tip
# labels (e.g., "Red_Drum_1"), so no renaming is needed.
# Avoiding the read.dna() -> as.character() -> DNAStringSet()
# round-trip prevents encoding artifacts for ambiguous bases.
# ============================================================
cat("\n--- STEP 2: Reading FASTA ---\n")

dna_strings <- Biostrings::readDNAStringSet(
  "FASTA/gulf_fish_COI_combined.fasta"
)

cat("  Sequences read:", length(dna_strings), "\n")
cat("  Tip labels:\n")
print(names(dna_strings))


# ============================================================
# STEP 3: Multiple sequence alignment with MUSCLE
# WHY: Before counting differences between sequences, every
# position must correspond to the same homologous site in every
# other sequence (i.e., they must be positionally homologous).
# MUSCLE is fast and accurate for mitochondrial protein-coding
# genes like COI.
# ============================================================
cat("\n--- STEP 3: Running MUSCLE alignment (~30 s) ---\n")

msa_result <- msa(dna_strings, method = "Muscle")

aln_length <- ncol(as.matrix(msa_result))
cat("  Alignment length:", aln_length, "bp\n")
cat("  Number of taxa:  ", nrow(as.matrix(msa_result)), "\n")


# ============================================================
# STEP 4: Convert alignment to ape and phangorn formats
# WHY: Different tree-building functions require different
# internal representations of the same data.
#   aln_dnabin — used by ape::dist.dna() for NJ distances
#   aln_phydat — used by phangorn::pml() and modelTest() for ML
# msaConvert() handles the translation cleanly.
# ============================================================
cat("\n--- STEP 4: Converting alignment formats ---\n")

aln_dnabin <- msaConvert(msa_result, type = "ape::DNAbin")
aln_phydat <- msaConvert(msa_result, type = "phangorn::phyDat")

cat("  DNAbin dimensions: ", nrow(aln_dnabin), "taxa x",
    ncol(aln_dnabin), "sites\n")


# ============================================================
# STEP 5: Define shared color palette and tip-color helper
# WHY: Storing these in out/rdata/alignment_objects.RData means all three
# tree scripts load the identical palette — one edit here
# propagates everywhere.
#
# Color assignments by family:
#   Sciaenidae      (drums/croakers/trout) — blue   #2166ac
#   Paralichthyidae (flounders)            — orange #d95f02
#   Sparidae        (sheepshead/porgies)   — green  #1b9e77
#   Mugilidae       (mullets)              — gray   #888888
#   Carcharhinidae  (sharks — outgroup)    — red    #e41a1c
# ============================================================
family_colors <- c(
  "Red_Drum"          = "#2166ac",
  "Spec_Trout"        = "#2166ac",
  "Black_Drum"        = "#2166ac",
  "Atlantic_Croaker"  = "#2166ac",
  "Southern_Flounder" = "#d95f02",
  "Sheepshead"        = "#1b9e77",
  "Flathead_Mullet"   = "#888888",
  "Bull_Shark"        = "#e41a1c"
)

# tip_color(): given a vector of tip labels like "Red_Drum_1",
# strips the trailing _1/_2/_3 suffix and returns the hex color
# for that species' family.
tip_color <- function(labels) {
  species <- sub("_[0-9]+$", "", labels)
  unname(family_colors[species])
}


# ============================================================
# STEP 6: Save alignment objects to .RData
# WHY: save() serializes R objects to a binary file.
# load() in any downstream script restores them exactly,
# including the tip_color function and family_colors vector.
# This avoids re-running the ~30 s alignment every time.
# ============================================================
cat("\n--- STEP 6: Saving out/rdata/alignment_objects.RData ---\n")

save(aln_dnabin, aln_phydat, family_colors, tip_color,
     file = "out/rdata/alignment_objects.RData")

cat("  Saved: out/rdata/alignment_objects.RData\n")
cat("  Objects stored: aln_dnabin, aln_phydat,",
    "family_colors, tip_color\n")


# ============================================================
# STEP 7: Export aligned sequences to FASTA and NEXUS
# WHY: Exporting the alignment lets you inspect it in SeqView
# or AliView, and provides the NEXUS file that MrBayes reads.
# ============================================================
cat("\n--- STEP 7: Writing alignment files ---\n")

# -- Aligned FASTA (for inspection / archiving) --------------
write.dna(
  aln_dnabin,
  file   = "FASTA/gulf_fish_COI_aligned.fasta",
  format = "fasta",
  colsep = ""
)
cat("  Saved: FASTA/gulf_fish_COI_aligned.fasta\n")

# -- NEXUS for MrBayes ---------------------------------------
# WHY: MrBayes reads NEXUS format. write.nexus.data() writes the
# DATA block; 03_tree_Bayesian.R appends the MRBAYES block.
# interleaved = FALSE writes each sequence as one continuous
# line, which avoids rare parsing issues in some MrBayes versions.
aln_mat  <- as.character(aln_dnabin)   # character matrix ntax x nchar
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


# ============================================================
# STEP 8: Print alignment summary
# ============================================================
cat("\n============================================================\n")
cat("ALIGNMENT SUMMARY\n")
cat("============================================================\n")
cat("  Input sequences : ", length(dna_strings), "\n")
cat("  Aligned length  : ", aln_length, "bp\n")
cat("  Taxa            : ", nrow(aln_dnabin), "\n")
cat("  Sites (aligned) : ", ncol(aln_dnabin), "\n")

# Count variable sites (columns where not all bases are identical)
aln_chars     <- as.character(aln_dnabin)
variable_sites <- sum(apply(aln_chars, 2, function(col) {
  length(unique(col[col != "-" & col != "n"])) > 1
}))
cat("  Variable sites  : ", variable_sites, "\n")

cat("\nNext steps:\n")
cat("  NJ tree   : Rscript 01_tree_NJ.R\n")
cat("  ML tree   : Rscript 02_tree_ML.R\n")
cat("  Bayesian  : Rscript 03_tree_Bayesian.R\n")
cat("============================================================\n")
