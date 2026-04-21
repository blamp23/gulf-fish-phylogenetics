# ============================================================
# 01_tree_NJ.R
# Neighbor-Joining phylogeny of Gulf of Mexico COI sequences
#
# Terminal: Rscript 01_tree_NJ.R
# Requires: out/rdata/alignment_objects.RData (run 00_prepare_alignment.R first)
#
# Outputs:
#   out/rdata/trees_NJ.RData
#   out/trees/tree_NJ_bullshark.pdf
#   out/trees/tree_NJ_mullet.pdf
# ============================================================

library(ape)

dir.create("out/rdata", recursive = TRUE, showWarnings = FALSE)
dir.create("out/trees", recursive = TRUE, showWarnings = FALSE)

if (!file.exists("out/rdata/alignment_objects.RData")) {
  stop("out/rdata/alignment_objects.RData not found.\n",
       "Run 00_prepare_alignment.R first.")
}
load("out/rdata/alignment_objects.RData")
cat("Loaded: aln_dnabin, aln_phydat, family_colors, tip_color\n")


cat("\n--- STEP 1: Computing K80 distance matrix ---\n")

dist_k80 <- dist.dna(aln_dnabin, model = "K80")

cat("  Distance matrix:", attr(dist_k80, "Size"), "x",
    attr(dist_k80, "Size"), "\n")
cat("  Mean pairwise distance:",
    round(mean(dist_k80, na.rm = TRUE), 4), "\n")


cat("\n--- STEP 2: NJ tree (Bull Shark outgroup) ---\n")

tree_nj_raw   <- nj(dist_k80)
outgrp_shark  <- grep("Bull_Shark", tree_nj_raw$tip.label, value = TRUE)
tree_nj_shark <- root(tree_nj_raw, outgroup = outgrp_shark,
                      resolve.root = TRUE)

cat("  Tips:", length(tree_nj_shark$tip.label), "\n")


cat("\n--- STEP 3: NJ tree (Flathead Mullet outgroup) ---\n")

tree_nj_mullet_raw    <- nj(dist_k80)
outgrp_mullet         <- grep("Flathead_Mullet",
                               tree_nj_mullet_raw$tip.label, value = TRUE)
tree_nj_mullet_rooted <- root(tree_nj_mullet_raw, outgroup = outgrp_mullet,
                               resolve.root = TRUE)

cat("  Done.\n")


cat("\n--- STEP 4: Bootstrap (1000 replicates) ---\n")

boot_nj_shark <- boot.phylo(
  tree_nj_shark,
  aln_dnabin,
  function(x) root(
    nj(dist.dna(x, model = "K80")),
    outgroup = grep("Bull_Shark", rownames(x), value = TRUE),
    resolve.root = TRUE
  ),
  B = 1000,
  quiet = TRUE
)

boot_nj_mullet <- boot.phylo(
  tree_nj_mullet_rooted,
  aln_dnabin,
  function(x) root(
    nj(dist.dna(x, model = "K80")),
    outgroup = grep("Flathead_Mullet", rownames(x), value = TRUE),
    resolve.root = TRUE
  ),
  B = 1000,
  quiet = TRUE
)

cat("  Done.\n")


cat("\n--- STEP 5: Saving out/rdata/trees_NJ.RData ---\n")

save(tree_nj_shark, tree_nj_mullet_rooted, boot_nj_shark, boot_nj_mullet,
     dist_k80, file = "out/rdata/trees_NJ.RData")

cat("  Saved: out/rdata/trees_NJ.RData\n")


.legend_strip <- function() {
  par(mar = c(0, 1, 0, 1))
  plot.new()
  legend(
    "center",
    legend = c("Sciaenidae (drums / croakers / trout)",
               "Paralichthyidae (flounders)",
               "Sparidae (sheepshead / porgies)",
               "Mugilidae (mullets)",
               "Carcharhinidae (sharks — outgroup)",
               "Salmonidae (Rainbow, Brown Trout)"),
    col    = c("#2166ac", "#d95f02", "#1b9e77", "#888888", "#e41a1c", "#FA8072"),
    pch    = 15,
    cex    = 0.75,
    horiz  = FALSE,
    ncol   = 3,
    bty    = "n",
    pt.cex = 1.2
  )
}

cat("\n--- STEP 6a: Saving out/trees/tree_NJ_bullshark.pdf ---\n")

pdf("out/trees/tree_NJ_bullshark.pdf", width = 9, height = 9)
layout(matrix(c(1, 2), nrow = 2), heights = c(7.5, 1.5))

par(mar = c(1, 1, 3, 1))
plot.phylo(
  tree_nj_shark,
  main       = "Neighbor-Joining (K80)\nOutgroup: Bull Shark | Node labels = bootstrap %",
  tip.color  = tip_color(tree_nj_shark$tip.label),
  cex        = 0.78,
  edge.width = 1.6,
  no.margin  = FALSE
)
nodelabels(boot_nj_shark / 10,
           adj   = c(1.1, -0.3),
           cex   = 0.55,
           frame = "none",
           col   = "gray25")
add.scale.bar(cex = 0.7, col = "gray30")
.legend_strip()
dev.off()
cat("  Saved: out/trees/tree_NJ_bullshark.pdf\n")

cat("\n--- STEP 6b: Saving out/trees/tree_NJ_mullet.pdf ---\n")

pdf("out/trees/tree_NJ_mullet.pdf", width = 9, height = 9)
layout(matrix(c(1, 2), nrow = 2), heights = c(7.5, 1.5))

par(mar = c(1, 1, 3, 1))
plot.phylo(
  tree_nj_mullet_rooted,
  main       = "Neighbor-Joining (K80)\nOutgroup: Flathead Mullet | Node labels = bootstrap %",
  tip.color  = tip_color(tree_nj_mullet_rooted$tip.label),
  cex        = 0.78,
  edge.width = 1.6,
  no.margin  = FALSE
)
nodelabels(boot_nj_mullet / 10,
           adj   = c(1.1, -0.3),
           cex   = 0.55,
           frame = "none",
           col   = "gray25")
add.scale.bar(cex = 0.7, col = "gray30")
.legend_strip()
dev.off()
cat("  Saved: out/trees/tree_NJ_mullet.pdf\n")
