# ============================================================
# 02_tree_ML.R
# Maximum Likelihood phylogeny of Gulf of Mexico COI sequences
# Author: [Your Name]  |  Texas A&M University
#
# Terminal: Rscript 02_tree_ML.R
# Requires: out/rdata/alignment_objects.RData (run 00_prepare_alignment.R first)
# Note: 1000 bootstrap replicates will take 10-20 minutes
#
# Outputs:
#   out/rdata/trees_ML.RData  — rooted ML tree with bootstrap node labels
#   out/trees/tree_ML.pdf     — single-panel figure, 12x9, bootstrap % on nodes
# ============================================================


# ============================================================
# STEP 1: Load libraries and alignment objects
# ============================================================
library(ape)
library(phangorn)

dir.create("out/rdata", recursive = TRUE, showWarnings = FALSE)
dir.create("out/trees", recursive = TRUE, showWarnings = FALSE)

if (!file.exists("out/rdata/alignment_objects.RData")) {
  stop("out/rdata/alignment_objects.RData not found.\n",
       "Run 00_prepare_alignment.R first.")
}
load("out/rdata/alignment_objects.RData")
cat("Loaded: aln_dnabin, aln_phydat, family_colors, tip_color\n")


# ============================================================
# STEP 2: Build NJ starting tree (K80)
# WHY: ML optimization needs a starting topology to begin its
# search. A reasonable starting tree (rather than a random one)
# reduces the risk of getting trapped in a poor local optimum
# and speeds up convergence. We compute NJ quickly here — this
# is *not* the final tree, just a launch point for ML.
# ============================================================
cat("\n--- STEP 2: Building NJ starting tree ---\n")

dist_k80      <- dist.dna(aln_dnabin, model = "K80")
tree_nj_start <- nj(dist_k80)
outgrp_shark  <- grep("Bull_Shark", tree_nj_start$tip.label, value = TRUE)
tree_nj_start <- root(tree_nj_start, outgroup = outgrp_shark,
                       resolve.root = TRUE)

cat("  NJ starting tree built.\n")


# ============================================================
# STEP 3: Model selection with modelTest()
#
# WHAT THIS IS DOING:
#   modelTest() fits dozens of substitution models to the data
#   and ranks them by AIC (Akaike Information Criterion).
#   Lower AIC = better balance of model fit and complexity.
#
# WHY bother if we are using GTR+Gamma+I anyway?
#   Printing the best AIC model is good scientific practice —
#   it tells you whether your chosen GTR+Gamma+I is justified
#   or whether a simpler model (e.g., HKY+Gamma) fits just as
#   well. If modelTest picks something much simpler, note it
#   in your Methods section.
# ============================================================
cat("\n--- STEP 3: Model selection (modelTest) ---\n")
cat("  Testing substitution models — may take a minute...\n")

mt         <- modelTest(aln_phydat, tree = tree_nj_start,
                         multicore = FALSE)
best_model <- mt$Model[which.min(mt$AIC)]

cat("  Best-fit model by AIC:", best_model, "\n")
cat("  (Using GTR+Gamma+I for ML — conservative, publication-standard)\n")
cat("  Top 5 models by AIC:\n")
print(head(mt[order(mt$AIC), c("Model", "AIC", "AICw")], 5))


# ============================================================
# STEP 4: Fit Maximum Likelihood tree (GTR+Gamma+I)
#
# WHAT ML IS DOING:
#   ML finds the tree topology T and branch lengths b that
#   maximize the probability of observing the alignment A given
#   the model M: argmax P(A | T, b, M).
#   Unlike NJ, ML evaluates the entire alignment simultaneously
#   rather than compressing it into pairwise distances first.
#
# Parameter choices:
#   model = "GTR"      — General Time Reversible: 6 free
#                        substitution rates (most flexible standard model)
#   optGamma = TRUE    — estimates the Gamma shape parameter α,
#                        which allows substitution rate to vary
#                        across sites (fast-evolving 3rd codon
#                        positions vs. slow 1st/2nd positions)
#   optInv = TRUE      — estimates the proportion of invariant
#                        sites (sites that never change); helps
#                        with COI, which has conserved regions
#   rearrangement = "stochastic" — combines random NNI moves with
#                        nearest-neighbor interchange; more thorough
#                        than plain NNI at the cost of more time
# ============================================================
cat("\n--- STEP 4: Fitting ML tree (GTR+Gamma+I, stochastic) ---\n")
cat("  This may take a few minutes...\n")

fit_start <- pml(tree_nj_start, data = aln_phydat)
fit_ml    <- optim.pml(
  fit_start,
  model         = "GTR",
  optGamma      = TRUE,
  optInv        = TRUE,
  rearrangement = "stochastic",
  control       = pml.control(trace = 0)
)

cat("  Log-likelihood:", round(fit_ml$logLik, 2), "\n")
cat("  Gamma shape (α):", round(fit_ml$shape, 4), "\n")
cat("  Prop. invariant:", round(fit_ml$inv, 4), "\n")
cat("  ML tree fitted.\n")


# ============================================================
# STEP 5: Bootstrap (1000 replicates)
#
# WHAT BOOTSTRAP IS DOING:
#   Each replicate randomly resamples columns from the alignment
#   with replacement (same total length), rebuilds the ML tree,
#   and records which internal nodes appear. The bootstrap value
#   at a node = the percentage of replicates in which that exact
#   split was recovered.
#
# WHY 1000?
#   100 replicates is the traditional minimum; 1000 gives more
#   stable estimates, especially for nodes with intermediate
#   support (50-80%). For a class project, 100 is acceptable
#   if time is short.
#
# Interpreting bootstrap values:
#   >= 95%  — strong support, confidently report
#   70-94%  — moderate support, worth reporting with caveat
#   50-69%  — weak support, report cautiously
#   < 50%   — not supported; often collapsed in published trees
# ============================================================
cat("\n--- STEP 5: Bootstrap (1000 replicates) ---\n")
cat("  Estimated time: 10-20 minutes. Go get coffee.\n")

bs_ml <- bootstrap.pml(
  fit_ml,
  bs      = 1000,
  optNni  = TRUE,
  control = pml.control(trace = 0)
)

# Attach bootstrap values as node labels; type = "n" = no plot yet
tree_ml_bs <- plotBS(fit_ml$tree, bs_ml, type = "n")

cat("  Bootstrap done.\n")


# ============================================================
# STEP 6: Root on Bull Shark
# WHY: See 01_tree_NJ.R STEP 3 for full rationale.
# ============================================================
cat("\n--- STEP 6: Rooting ML tree on Bull Shark ---\n")

outgrp_ml      <- grep("Bull_Shark", tree_ml_bs$tip.label, value = TRUE)
tree_ml_rooted <- root(tree_ml_bs, outgroup = outgrp_ml,
                        resolve.root = TRUE)

cat("  Done.\n")


# ============================================================
# STEP 7: Save tree object
# ============================================================
cat("\n--- STEP 7: Saving out/rdata/trees_ML.RData ---\n")

save(tree_ml_rooted, fit_ml, bs_ml, file = "out/rdata/trees_ML.RData")

cat("  Saved: out/rdata/trees_ML.RData\n")


# ============================================================
# STEP 8: Plot and save PDF
# ============================================================
cat("\n--- STEP 8: Saving out/trees/tree_ML.pdf ---\n")

pdf("out/trees/tree_ML.pdf", width = 12, height = 9)

layout(
  matrix(c(1,
            2), nrow = 2, byrow = TRUE),
  heights = c(7.5, 1.5)
)

# ── Panel 1: ML tree with bootstrap values ──────────────────
par(mar = c(1, 1, 3, 1))
plot.phylo(
  tree_ml_rooted,
  main       = "Maximum Likelihood (GTR+\u0393+I)\nOutgroup: Bull Shark | Node labels = bootstrap %",
  tip.color  = tip_color(tree_ml_rooted$tip.label),
  cex        = 0.78,
  edge.width = 1.6,
  no.margin  = FALSE
)



bs_labels <- round(as.numeric(tree_ml_rooted$node.label) * 100)
bs_labels[is.na(bs_labels) | bs_labels < 1] <- ""



nodelabels(
  bs_labels,
  adj   = c(1.1, -0.3),
  cex   = 0.55,
  frame = "none",
  col   = "gray25"
)
# WHY: Branch lengths in ML trees represent substitutions per site
# under GTR+Gamma+I. The scale bar lets the reader compare
# evolutionary distances between clades.
add.scale.bar(cex = 0.7, col = "gray30")

# ── Panel 2: Shared legend strip ────────────────────────────
par(mar = c(0, 1, 0, 1))
plot.new()
legend(
  "center",
  legend  = c("Sciaenidae (drums / croakers / trout)",
               "Paralichthyidae (flounders)",
               "Sparidae (sheepshead / porgies)",
               "Mugilidae (mullets)",
               "Carcharhinidae (sharks — outgroup)", 
              "Salmonidae (Rainbow / Brown Trout)"),
  col     = c("#2166ac", "#d95f02", "#1b9e77", "#888888", "#e41a1c",  "#FA8072"),
  pch     = 15,
  cex     = 0.75,
  horiz   = FALSE,
  ncol    = 3,
  bty     = "n",
  pt.cex  = 1.2
)

dev.off()
cat("  Saved: out/trees/tree_ML.pdf\n")


# ============================================================
# STEP 9: Interpretation notes
# ============================================================
cat("\n")
cat("============================================================\n")
cat("ML INTERPRETATION NOTES\n")
cat("============================================================\n\n")

cat("How to read bootstrap values:\n")
cat("  Each number at an internal node (0-100) is the percentage\n")
cat("  of 1000 bootstrap replicates that recovered that exact\n")
cat("  split. Higher = more reliable.\n")
cat("  Rule of thumb: report nodes >= 70% as 'supported'.\n\n")

cat("Comparing ML to NJ (from 01_tree_NJ.R):\n")
cat("  Nodes with high ML bootstrap support (>= 70%) that also\n")
cat("  appear in both NJ trees are the most trustworthy results.\n")
cat("  If a strongly-supported ML node is absent from NJ, check\n")
cat("  whether the NJ outgroup placement might be distorting it.\n\n")

cat("Sciaenidae monophyly:\n")
cat("  COI evidence supports Sciaenidae as monophyletic if all\n")
cat("  blue tips (Red Drum, Black Drum, Speckled Trout, Atlantic\n")
cat("  Croaker) cluster together with bootstrap >= 70%.\n\n")

cat("Model parameters:\n")
cat("  Log-likelihood:", round(fit_ml$logLik, 2), "\n")
cat("  Gamma shape α :", round(fit_ml$shape,  4),
    " (< 1 = many slow sites + few fast sites)\n")
cat("  Prop. invariant:", round(fit_ml$inv,   4),
    " (fraction of sites assumed to never change)\n")

cat("============================================================\n")
