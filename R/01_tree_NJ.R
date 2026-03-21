# ============================================================
# 01_tree_NJ.R
# Neighbor-Joining phylogeny of Gulf of Mexico COI sequences
# Author: [Your Name]  |  Texas A&M University
#
# Terminal: Rscript 01_tree_NJ.R
# Requires: out/rdata/alignment_objects.RData (run 00_prepare_alignment.R first)
#
# Outputs:
#   out/rdata/trees_NJ.RData                   — tree + bootstrap objects
#   out/trees/tree_NJ_bullshark.pdf            — NJ tree, Bull Shark outgroup
#   out/trees/tree_NJ_mullet.pdf               — NJ tree, Flathead Mullet outgroup
#   out/trees/tree_nj_outgroup_*.png (x8)      — per-outgroup PNGs for website widget
# ============================================================


# ============================================================
# STEP 1: Load libraries and alignment objects
# WHY: alignment_objects.RData contains aln_dnabin, aln_phydat,
# family_colors, and tip_color — everything built by the shared
# setup script. Loading it here avoids re-running the ~30 s
# alignment.
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


# ============================================================
# STEP 2: Compute pairwise K80 distance matrix
#
# WHAT THIS IS DOING:
#   dist.dna() calculates the evolutionary distance between every
#   pair of sequences. model = "K80" applies the Kimura
#   2-parameter correction.
#
# WHY K80 for COI?
#   Raw p-distance (fraction of differing sites) underestimates
#   the true number of substitutions because the same site can
#   change more than once (multiple hits). K80 corrects for this.
#   It also accounts for transition/transversion bias — transitions
#   (A<->G, C<->T) occur ~2-3x more often than transversions in
#   mitochondrial DNA. K80 handles this with just two parameters,
#   striking a good balance between accuracy and simplicity.
# ============================================================
cat("\n--- STEP 2: Computing K80 distance matrix ---\n")

dist_k80 <- dist.dna(aln_dnabin, model = "K80")

cat("  Distance matrix:", attr(dist_k80, "Size"), "x",
    attr(dist_k80, "Size"), "\n")
cat("  Mean pairwise distance:",
    round(mean(dist_k80, na.rm = TRUE), 4), "\n")


# ============================================================
# STEP 3: Build NJ tree — Bull Shark outgroup
#
# WHAT NJ IS DOING:
#   Neighbor-Joining (Saitou & Nei 1987) builds a tree by
#   iteratively joining the pair of taxa whose connection
#   minimizes the total branch length of the tree (the
#   "minimum evolution" criterion). It is fast — O(n^2) — and
#   gives a good first approximation of the true topology.
#
# WHY root on Bull Shark?
#   Bull Shark (Carcharhinus leucas, Carcharhinidae) is a
#   cartilaginous fish (Chondrichthyes), while all other taxa
#   are ray-finned fishes (Actinopterygii). The two lineages
#   diverged ~420 Mya — far deeper than any split within
#   the ingroup. Placing the root here gives the tree a
#   meaningful evolutionary direction (ancestor -> descendant).
# ============================================================
cat("\n--- STEP 3: NJ tree (Bull Shark outgroup) ---\n")

tree_nj_raw      <- nj(dist_k80)
tree_nj_unrooted <- tree_nj_raw   # unrooted copy — used in STEP 8 widget export
outgrp_shark     <- grep("Bull_Shark", tree_nj_raw$tip.label, value = TRUE)
tree_nj_shark    <- root(tree_nj_raw, outgroup = outgrp_shark,
                          resolve.root = TRUE)

cat("  Tips:", length(tree_nj_shark$tip.label), "\n")
cat("  Done.\n")


# ============================================================
# STEP 4: Build NJ tree — Flathead Mullet alternative outgroup
#
# WHY a second outgroup?
#   Outgroup choice determines where the root falls, which
#   changes which ingroup clades appear "sister" to each other.
#   Flathead Mullet (Mugil cephalus, Mugilidae) is a ray-finned
#   fish like the ingroup but belongs to a different order
#   (Mugiliformes vs. Perciformes / Pleuronectiformes).
#
#   Comparing this tree with the Bull Shark-rooted tree shows
#   which ingroup groupings are *stable* (appear in both) vs.
#   which depend on outgroup placement (appear in only one).
#   Stable groupings are more trustworthy to report.
#
# NOTE: We reuse the same dist_k80 matrix — only the rooting changes.
# ============================================================
cat("\n--- STEP 4: NJ tree (Flathead Mullet outgroup) ---\n")

tree_nj_mullet_raw    <- nj(dist_k80)
outgrp_mullet         <- grep("Flathead_Mullet",
                               tree_nj_mullet_raw$tip.label, value = TRUE)
tree_nj_mullet_rooted <- root(tree_nj_mullet_raw, outgroup = outgrp_mullet,
                               resolve.root = TRUE)

cat("  Done.\n")


# Bootstrap NJ — Bull Shark outgroup
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

# Bootstrap NJ — Mullet outgroup
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


# ============================================================
# STEP 5: Save tree objects
# WHY: Storing the rooted trees lets you reload them later for
# comparison with ML/Bayesian results without re-running NJ.
# ============================================================
cat("\n--- STEP 5: Saving out/rdata/trees_NJ.RData ---\n")

save(tree_nj_shark, tree_nj_mullet_rooted, boot_nj_shark, boot_nj_mullet, dist_k80,
     tree_nj_unrooted,
     file = "out/rdata/trees_NJ.RData")

cat("  Saved: out/rdata/trees_NJ.RData\n")


# ============================================================
# STEP 6: Plot and save PDFs
# WHY: One PDF per outgroup keeps figures focused and makes it
# easier to place them side-by-side in a report without rescaling.
# layout() gives each tree its own panel plus a shared legend strip.
# ============================================================

# Reusable legend strip — identical across both PDFs
.legend_strip <- function() {
  par(mar = c(0, 1, 0, 1))
  plot.new()
  legend(
    "center",
    legend  = c("Sciaenidae (drums / croakers / trout)",
                "Paralichthyidae (flounders)",
                "Sparidae (sheepshead / porgies)",
                "Mugilidae (mullets)",
                "Carcharhinidae (sharks — outgroup)", 
                "Salmonidae (Rainbow, Brown Trout)"
                ),
    col     = c("#2166ac", "#d95f02", "#1b9e77", "#888888", "#e41a1c", "#FA8072" ),
    pch     = 15,
    cex     = 0.75,      # slightly smaller text
    horiz   = FALSE,     # stack vertically instead of horizontal
    ncol    = 3,         # but arrange in 3 columns so it's compact
    bty     = "n",
    pt.cex  = 1.2
  )
}
# ── PDF 1: Bull Shark outgroup ───────────────────────────────
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
# WHY bootstrap on NJ? boot.phylo() resamples alignment columns,
# rebuilds the NJ tree each time, and counts how often each node
# appears. Values >= 70% indicate that node is reliably recovered
# under this distance method and dataset.
nodelabels(boot_nj_shark/10,
           adj   = c(1.1, -0.3),
           cex   = 0.55,
           frame = "none",
           col   = "gray25")
# WHY scale bar? Branch lengths in NJ trees represent substitutions
# per site under the K80 model. The scale bar lets the reader gauge
# how much sequence divergence separates any two taxa.
add.scale.bar(cex = 0.7, col = "gray30")

.legend_strip()
dev.off()
cat("  Saved: out/trees/tree_NJ_bullshark.pdf\n")

# ── PDF 2: Flathead Mullet outgroup ──────────────────────────
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
# WHY bootstrap on NJ? See note above — same logic applies here.
# Comparing bootstrap values between the two outgroup trees shows
# which nodes are stable regardless of how the tree is rooted.
nodelabels(boot_nj_mullet/10,
           adj   = c(1.1, -0.3),
           cex   = 0.55,
           frame = "none",
           col   = "gray25")
# WHY scale bar? Scale bar shows substitutions per site under K80 —
# the unit of branch length in this tree.
add.scale.bar(cex = 0.7, col = "gray30")

.legend_strip()
dev.off()
cat("  Saved: out/trees/tree_NJ_mullet.pdf\n")


# ============================================================
# STEP 7: Interpretation notes
# ============================================================
cat("\n")
cat("============================================================\n")
cat("NJ INTERPRETATION NOTES\n")
cat("============================================================\n\n")

cat("What NJ is and is not:\n")
cat("  NJ is a fast heuristic — it finds a tree that approximates\n")
cat("  minimum total branch length. It does NOT evaluate an\n")
cat("  explicit model of evolution, so it gives no node-support\n")
cat("  values. Use it as a quick sanity check and a starting\n")
cat("  topology for ML (02_tree_ML.R).\n\n")

cat("Outgroup comparison:\n")
cat("  Groupings that appear in BOTH panels are robust to outgroup\n")
cat("  placement and are the safest to report. Groupings present\n")
cat("  in only one panel should be treated cautiously and checked\n")
cat("  against the ML bootstrap values.\n\n")

cat("Expected pattern for Sciaenidae:\n")
cat("  Red Drum, Black Drum, Speckled Trout, and Atlantic Croaker\n")
cat("  (all blue) should cluster together as a monophyletic group\n")
cat("  in both trees if COI resolves family-level relationships.\n")
cat("  If any blue tip falls outside the blue cluster, note it —\n")
cat("  it may indicate a misidentified sequence or true paraphyly.\n\n")

cat("Long Bull Shark branch:\n")
cat("  Expect a noticeably long branch leading to Bull Shark in\n")
cat("  Panel 1. This reflects ~420 Myr of independent evolution\n")
cat("  (Chondrichthyes vs. Actinopterygii) and is biologically\n")
cat("  expected, not an error.\n")

cat("============================================================\n")


# ============================================================
# STEP 8: Export per-outgroup trees for interactive website widget
#
# WHY: The explore-phylogeny.html page includes an interactive
# widget that lets students pick any species as outgroup and see
# how the tree changes. Rather than running R in the browser,
# we pre-render one PNG per outgroup choice here and serve them
# as static images. This keeps the website fast and dependency-free
# while still delivering a dynamic, educational experience.
#
# Each PNG uses cladogram style (use.edge.length=FALSE) so branch
# lengths don't dominate — the goal is topology comparison, not
# measuring divergence. No bootstrap labels are added so the images
# stay clean; bootstrap information is discussed in the page text.
# ============================================================
cat("\n--- STEP 8: Exporting per-outgroup PNGs for website widget ---\n")

# Legend strip sized for square PNGs: ncol=2 to fit entries compactly
.legend_strip_png <- function() {
  par(mar = c(0, 1, 0, 1))
  plot.new()
  legend(
    "center",
    legend  = c("Sciaenidae (drums / croakers / trout)",
                "Paralichthyidae (flounders)",
                "Sparidae (sheepshead / porgies)",
                "Mugilidae (mullets)",
                "Carcharhinidae (sharks — outgroup)",
                "Salmonidae (trout — outgroup)"),
    col     = c("#2166ac", "#d95f02", "#1b9e77", "#888888", "#e41a1c", "#FA8072"),
    pch     = 15,
    cex     = 0.7,
    horiz   = FALSE,
    ncol    = 2,
    bty     = "n",
    pt.cex  = 1.2
  )
}

# Map: outgroup tip label -> output filename
outgrp_map <- list(
  list(tip = "Bull_Shark_1",        file = "tree_nj_outgroup_bullshark.png"),
  list(tip = "Flathead_Mullet_1",   file = "tree_nj_outgroup_mullet.png"),
  list(tip = "Southern_Flounder_1", file = "tree_nj_outgroup_flounder.png"),
  list(tip = "Sheepshead_1",        file = "tree_nj_outgroup_sheepshead.png"),
  list(tip = "Red_Drum_1",          file = "tree_nj_outgroup_reddrum.png"),
  list(tip = "Black_Drum_1",        file = "tree_nj_outgroup_blackdrum.png"),
  list(tip = "Spec_Trout_1",        file = "tree_nj_outgroup_trout.png"),
  list(tip = "Atlantic_Croaker_1",  file = "tree_nj_outgroup_croaker.png"),
  list(tip = "Brown_Trout_1",        file = "tree_nj_outgroup_browntrout.png"),
  list(tip = "Rainbow_Trout_1",        file = "tree_nj_outgroup_rainbowtrout.png")
)

for (og in outgrp_map) {

  # Guard: skip gracefully if the tip label doesn't exist in this dataset
  if (!og$tip %in% tree_nj_unrooted$tip.label) {
    cat("  WARNING: tip '", og$tip, "' not in tree — skipping.\n", sep = "")
    next
  }

  # Re-root the unrooted NJ tree on the single _1 replicate of this species.
  # WHY _1 only: one representative tip gives a clean, unambiguous root
  # placement. Using all three replicates can create a trifurcation at
  # the root that obscures the topology the widget is meant to illustrate.
  tree_og <- root(tree_nj_unrooted, outgroup = og$tip, resolve.root = TRUE)

  out_path <- paste0("out/trees/", og$file)
  png(out_path, width = 900, height = 900, res = 120)
  layout(matrix(c(1, 2), nrow = 2), heights = c(8.5, 1.5))

  # WHY use.edge.length=FALSE: cladogram style makes topology comparison
  # easier — students focus on which species cluster together rather than
  # on branch lengths, which vary across outgroup choices in ways that
  # distract from the topological lesson. No title: the website widget
  # provides the label dynamically.
  par(mar = c(1, 1, 1, 1))
  plot.phylo(
    tree_og,
    use.edge.length = FALSE,
    tip.color       = tip_color(tree_og$tip.label),
    cex             = 0.72,
    edge.width      = 1.4,
    no.margin       = FALSE
  )

  .legend_strip_png()
  dev.off()

  cat("  Saved: out/trees/", og$file, "\n", sep = "")
}

cat("\nSTEP 8 complete — 8 PNGs written to out/trees/\n")
cat("Copy to images/ folder to activate the interactive widget.\n")
