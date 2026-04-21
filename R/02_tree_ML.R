# ============================================================
# 02_tree_ML.R
# Maximum Likelihood phylogeny of Gulf of Mexico COI sequences
#
# Terminal: Rscript 02_tree_ML.R
# Requires: out/rdata/alignment_objects.RData (run 00_prepare_alignment.R first)
# Note: 1000 bootstrap replicates will take 10-20 minutes
#
# Outputs:
#   out/rdata/trees_ML.RData
#   out/trees/tree_ML.pdf
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


cat("\n--- STEP 1: Building NJ starting tree ---\n")

dist_k80      <- dist.dna(aln_dnabin, model = "K80")
tree_nj_start <- nj(dist_k80)
outgrp_shark  <- grep("Bull_Shark", tree_nj_start$tip.label, value = TRUE)
tree_nj_start <- root(tree_nj_start, outgroup = outgrp_shark,
                       resolve.root = TRUE)

cat("  NJ starting tree built.\n")


cat("\n--- STEP 2: Model selection (modelTest) ---\n")
cat("  Testing substitution models...\n")

mt         <- modelTest(aln_phydat, tree = tree_nj_start, multicore = FALSE)
best_model <- mt$Model[which.min(mt$AIC)]

cat("  Best-fit model by AIC:", best_model, "\n")
cat("  Top 5 models by AIC:\n")
print(head(mt[order(mt$AIC), c("Model", "AIC", "AICw")], 5))


cat("\n--- STEP 3: Fitting ML tree (GTR+Gamma+I, stochastic) ---\n")
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


cat("\n--- STEP 4: Bootstrap (1000 replicates) ---\n")
cat("  Estimated time: 10-20 minutes.\n")

bs_ml <- bootstrap.pml(
  fit_ml,
  bs      = 1000,
  optNni  = TRUE,
  control = pml.control(trace = 0)
)

tree_ml_bs <- plotBS(fit_ml$tree, bs_ml, type = "n")

cat("  Bootstrap done.\n")


cat("\n--- STEP 5: Rooting ML tree on Bull Shark ---\n")

outgrp_ml      <- grep("Bull_Shark", tree_ml_bs$tip.label, value = TRUE)
tree_ml_rooted <- root(tree_ml_bs, outgroup = outgrp_ml,
                        resolve.root = TRUE)

cat("  Done.\n")


cat("\n--- STEP 6: Saving out/rdata/trees_ML.RData ---\n")

save(tree_ml_rooted, fit_ml, bs_ml, file = "out/rdata/trees_ML.RData")

cat("  Saved: out/rdata/trees_ML.RData\n")


cat("\n--- STEP 7: Saving out/trees/tree_ML.pdf ---\n")

pdf("out/trees/tree_ML.pdf", width = 12, height = 9)

layout(
  matrix(c(1, 2), nrow = 2, byrow = TRUE),
  heights = c(7.5, 1.5)
)

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
add.scale.bar(cex = 0.7, col = "gray30")

par(mar = c(0, 1, 0, 1))
plot.new()
legend(
  "center",
  legend = c("Sciaenidae (drums / croakers / trout)",
              "Paralichthyidae (flounders)",
              "Sparidae (sheepshead / porgies)",
              "Mugilidae (mullets)",
              "Carcharhinidae (sharks — outgroup)",
              "Salmonidae (Rainbow / Brown Trout)"),
  col    = c("#2166ac", "#d95f02", "#1b9e77", "#888888", "#e41a1c", "#FA8072"),
  pch    = 15,
  cex    = 0.75,
  horiz  = FALSE,
  ncol   = 3,
  bty    = "n",
  pt.cex = 1.2
)

dev.off()
cat("  Saved: out/trees/tree_ML.pdf\n")
