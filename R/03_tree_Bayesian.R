# ============================================================
# 03_tree_Bayesian.R
# Bayesian phylogeny via MrBayes — Gulf of Mexico COI sequences
# Author: [Your Name]  |  Texas A&M University
#
# WORKFLOW:
#   Step 1 (R):        Rscript 03_tree_Bayesian.R
#                        <- checks the .nex file, appends MrBayes
#                           block if needed, then exits with instructions
#   Step 2 (Terminal): cd R/MrBayes
#                      mb gulf_fish_COI_aligned.nex
#                        <- runs MrBayes (~5-10 min)
#   Step 3 (R):        Rscript 03_tree_Bayesian.R
#                        <- detects output, reads consensus tree,
#                           saves RData + PDF
#
# Requires: out/rdata/alignment_objects.RData (run 00_prepare_alignment.R first)
#
# MrBayes output files live in: R/MrBayes/
#   gulf_fish_COI_aligned.nex.con.tre  <- consensus tree (what we read)
#   gulf_fish_COI_aligned.nex.t        <- sampled trees
#   gulf_fish_COI_aligned.nex.p        <- parameter log
#
# Outputs (after Step 3):
#   out/rdata/trees_Bayesian.RData
#   out/trees/tree_Bayesian.pdf  — 12x9, posterior probabilities on nodes
# ============================================================


# ============================================================
# STEP 1: Load libraries and alignment objects
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
# STEP 2: Verify the NEXUS file exists
# WHY: 00_prepare_alignment.R writes the NEXUS data block to
# MrBayes/gulf_fish_COI_aligned.nex. If that file is missing,
# something went wrong in the setup script.
# ============================================================
nex_path <- "MrBayes/gulf_fish_COI_aligned.nex"
con_path <- "MrBayes/gulf_fish_COI_aligned.nex.con.tre"

if (!file.exists(nex_path)) {
  stop(
    "NEXUS file not found: ", nex_path, "\n",
    "Run 00_prepare_alignment.R first to generate it."
  )
}
cat("\nFound NEXUS file:", nex_path, "\n")


# ============================================================
# STEP 3: Append MrBayes block to the NEXUS file (if needed)
#
# WHAT THIS IS DOING:
#   MrBayes reads analysis settings from a 'begin mrbayes;'
#   block appended to the NEXUS file. We check whether it is
#   already present before appending to avoid duplicating it
#   on reruns.
#
# MrBayes block settings explained:
#   set autoclose=yes      — don't pause for user input at the end
#   outgroup Bull_Shark_1  — sets the root for displayed trees
#   lset nst=6 rates=gamma — GTR (6 substitution types) + Gamma
#                            rate variation; matches our ML model
#   mcmc ngen=500000       — run 500,000 MCMC generations
#        samplefreq=500    — save a tree every 500 generations
#                            (= 1000 sampled trees total)
#        printfreq=1000    — print progress every 1000 generations
#   sumt burnin=250        — summarize trees, discarding first 250
#                            samples as burn-in (~25% of run)
#   sump burnin=250        — summarize parameters, same burn-in
# ============================================================
cat("\n--- STEP 3: Checking MrBayes block ---\n")

nex_lines    <- readLines(nex_path)
has_mb_block <- any(grepl("begin mrbayes", nex_lines, ignore.case = TRUE))

if (has_mb_block) {
  cat("  MrBayes block already present — skipping append.\n")
} else {
  cat("  Appending MrBayes block to", nex_path, "\n")

  mb_block <- c(
    "",
    "begin mrbayes;",
    "  set autoclose=yes;",
    "  outgroup Bull_Shark_1;",
    "  lset nst=6 rates=gamma;",
    "  mcmc ngen=500000 samplefreq=500 printfreq=1000;",
    "  sumt burnin=250;",
    "  sump burnin=250;",
    "end;"
  )

  write(mb_block, file = nex_path, append = TRUE)
  cat("  MrBayes block appended.\n")
}


# ============================================================
# STEP 4: Check whether MrBayes has already been run
#
# WHY: This script handles two separate invocations cleanly:
#   (a) First run — MrBayes hasn't been run yet. Print clear
#       instructions and exit gracefully so the user knows
#       exactly what to do next.
#   (b) Second run — MrBayes output exists. Read the consensus
#       tree and produce the figure.
# ============================================================
cat("\n--- STEP 4: Checking for MrBayes output ---\n")

if (!file.exists(con_path)) {

  # ── Case (a): MrBayes not yet run ─────────────────────────
  cat("\n")
  cat("============================================================\n")
  cat("MrBayes output not found. Run MrBayes now:\n")
  cat("============================================================\n\n")
  cat("  In a terminal:\n")
  cat("    cd R/MrBayes\n")
  cat("    mb gulf_fish_COI_aligned.nex\n\n")
  cat("  Expected runtime: ~5-10 minutes on a laptop.\n\n")
  cat("  MrBayes will create these output files in R/MrBayes/:\n")
  cat("    gulf_fish_COI_aligned.nex.con.tre  <- consensus tree\n")
  cat("    gulf_fish_COI_aligned.nex.t        <- sampled trees\n")
  cat("    gulf_fish_COI_aligned.nex.p        <- parameter log\n\n")
  cat("  After MrBayes finishes, re-run this script:\n")
  cat("    cd R/\n")
  cat("    Rscript 03_tree_Bayesian.R\n")
  cat("============================================================\n")

  # Exit gracefully — nothing more to do until MrBayes runs
  quit(save = "no", status = 0)
}

# ── Case (b): MrBayes output exists ───────────────────────
cat("  Found:", con_path, "\n")
cat("  Proceeding to read and plot consensus tree.\n")


# ============================================================
# STEP 5: Read MrBayes consensus tree
#
# WHAT THIS IS DOING:
#   read.nexus() parses the .con.tre file that MrBayes wrote.
#   This is the 50% majority-rule consensus tree — a node
#   appears only if it was found in > 50% of post-burn-in
#   sampled trees. Node labels are posterior probabilities (PP).
#
# Interpreting posterior probabilities:
#   PP is the probability that the clade is true given the data
#   and model. Unlike bootstrap (a resampling frequency), PP is
#   a direct probability statement.
#   >= 0.95 — strong support (commonly used threshold)
#   0.90-0.94 — moderate-strong
#   < 0.90  — weak; treat cautiously
#
# WHY: Bayesian PP tends to be higher than ML bootstrap for the
# same node — it is not directly comparable to bootstrap %.
# ============================================================
cat("\n--- STEP 5: Reading MrBayes consensus tree ---\n")

tree_bayes_raw <- read.nexus(con_path)
cat("  Tree read successfully.\n")
cat("  Tips:", length(tree_bayes_raw$tip.label), "\n")


# ============================================================
# STEP 6: Root on Bull Shark
# ============================================================
cat("\n--- STEP 6: Rooting on Bull Shark ---\n")

outgrp_bayes    <- grep("Bull_Shark", tree_bayes_raw$tip.label, value = TRUE)

if (length(outgrp_bayes) == 0) {
  warning("Bull_Shark tips not found in consensus tree. ",
          "Check tip label format in the .con.tre file.")
  tree_bayes_rooted <- tree_bayes_raw
} else {
  tree_bayes_rooted <- root(tree_bayes_raw, outgroup = outgrp_bayes,
                             resolve.root = TRUE)
  cat("  Rooted on:", paste(outgrp_bayes, collapse = ", "), "\n")
}


# ============================================================
# STEP 7: Save tree object
# ============================================================
cat("\n--- STEP 7: Saving out/rdata/trees_Bayesian.RData ---\n")

save(tree_bayes_rooted, file = "out/rdata/trees_Bayesian.RData")

cat("  Saved: out/rdata/trees_Bayesian.RData\n")


# ============================================================
# STEP 8: Plot and save PDF
#
# NOTE ON NODE LABELS:
#   MrBayes .con.tre files store posterior probabilities in the
#   node labels. Depending on the MrBayes version, the label
#   may be a raw probability (e.g., "0.98") or may be absent for
#   the root node. The nodelabels() call below plots whatever
#   is stored; empty labels will simply not appear.
# ============================================================
cat("\n--- STEP 8: Saving out/trees/tree_Bayesian.pdf ---\n")

pdf("out/trees/tree_Bayesian.pdf", width = 12, height = 9)

layout(
  matrix(c(1,
            2), nrow = 2, byrow = TRUE),
  heights = c(7.5, 1.5)
)

# ── Panel 1: Bayesian consensus tree ────────────────────────
par(mar = c(1, 1, 3, 1))
plot.phylo(
  tree_bayes_rooted,
  main       = "Bayesian Inference (GTR+\u0393, MrBayes)\nOutgroup: Bull Shark | Node labels = posterior probability",
  tip.color  = tip_color(tree_bayes_rooted$tip.label),
  cex        = 0.78,
  edge.width = 1.6,
  no.margin  = FALSE
)

# Extract posterior probabilities from raw MrBayes annotations
# WHY: read.nexus() strips the [&prob=...] annotations that MrBayes
# embeds in the Newick string. We extract the prob(percent) values
# directly from the raw file using regex, then match them to nodes
# in the order they appear in the tree string.

tree_line <- raw_lines[grep("tree con_50_majrule", raw_lines)]

# Extract all prob(percent)="XX" values in the order they appear
# These correspond to nodes in the tree, left to right
pp_raw <- regmatches(
  tree_line,
  gregexpr('prob\\(percent\\)="([0-9]+)"', tree_line)
)[[1]]

pp_vals <- as.numeric(gsub('prob\\(percent\\)="([0-9]+)"', "\\1", pp_raw))

cat("  Extracted", length(pp_vals), "posterior probability values\n")
cat("  Range:", min(pp_vals), "to", max(pp_vals), "\n")

# The number of internal nodes in the rooted tree
n_nodes <- tree_bayes_rooted$Nnode
cat("  Internal nodes in rooted tree:", n_nodes, "\n")

# Match PP values to nodes — take the last n_nodes values
# WHY: MrBayes annotates every node including tips; internal node
# annotations come after tip annotations in the Newick string.
# The final n_nodes values correspond to internal nodes.
if (length(pp_vals) >= n_nodes) {
  pp_nodes <- tail(pp_vals, n_nodes)
  pp_labels <- ifelse(pp_nodes >= 50, as.character(pp_nodes), "")
  
  nodelabels(pp_labels,
             adj   = c(1.1, -0.3),
             cex   = 0.55,
             frame = "none",
             col   = "gray25")
  cat("  Node labels applied.\n")
} else {
  cat("  Warning: could not match PP values to nodes — plotting without labels.\n")
}

# WHY: Branch lengths in Bayesian trees represent expected
# substitutions per site. Scale bar allows comparison of
# evolutionary distances across clades.
add.scale.bar(
  x   = 0,      # left edge of plot area
  y   = .5,     # below the lowest tip (tip 1 is at y=1)
  cex = 0.7,
  col = "gray30"
)

# ── Panel 2: Shared legend strip ────────────────────────────
par(mar = c(0, 1, 0, 1))
plot.new()
legend(
  "center",
  legend  = c("Sciaenidae (drums / croakers / trout)",
               "Paralichthyidae (flounders)",
               "Sparidae (sheepshead / porgies)",
               "Mugilidae (mullets)",
               "Carcharhinidae (sharks — outgroup)"),
  col     = c("#2166ac", "#d95f02", "#1b9e77", "#888888", "#e41a1c"),
  pch     = 15,
  cex     = 0.75,
  horiz   = FALSE,
  ncol    = 3,
  bty     = "n",
  pt.cex  = 1.2
)

dev.off()
cat("  Saved: out/trees/tree_Bayesian.pdf\n")


# ============================================================
# STEP 9: Interpretation notes
# ============================================================
cat("\n")
cat("============================================================\n")
cat("BAYESIAN INTERPRETATION NOTES\n")
cat("============================================================\n\n")

cat("What Bayesian inference is doing:\n")
cat("  MrBayes uses Markov chain Monte Carlo (MCMC) to sample\n")
cat("  trees in proportion to their posterior probability —\n")
cat("  P(tree | data, model). After discarding the burn-in,\n")
cat("  the remaining sampled trees form the posterior distribution.\n")
cat("  The consensus tree summarizes which clades appeared most\n")
cat("  often across those samples.\n\n")

cat("Posterior probability vs. ML bootstrap:\n")
cat("  PP (0-1) and bootstrap % (0-100) both measure node\n")
cat("  support, but they are NOT directly equivalent:\n")
cat("    - Bootstrap: how often a node appears when data are\n")
cat("      resampled (frequentist)\n")
cat("    - PP: the probability the clade is real given the\n")
cat("      data and model (Bayesian)\n")
cat("  PP tends to be numerically higher than bootstrap for the\n")
cat("  same node. A node with PP = 0.95 and bootstrap = 75%\n")
cat("  would both be considered 'supported'.\n\n")

cat("MCMC diagnostics to check in R/MrBayes/:\n")
cat("  Open the .p file in Tracer (free, from BEAST website)\n")
cat("  and verify that ESS > 200 for all parameters and that\n")
cat("  the likelihood trace looks like a 'fuzzy caterpillar'\n")
cat("  (well-mixed chains). If ESS is low, increase ngen.\n")

cat("============================================================\n")
