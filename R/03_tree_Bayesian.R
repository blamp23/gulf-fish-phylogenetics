# ============================================================
# 03_tree_Bayesian.R
# Bayesian phylogeny via MrBayes — Gulf of Mexico COI sequences
#
# WORKFLOW:
#   Step 1 (R):        Rscript 03_tree_Bayesian.R
#                        <- checks the .nex file, appends MrBayes
#                           block if needed, then exits with instructions
#   Step 2 (Terminal): cd R/MrBayes
#                      mb gulf_fish_COI_aligned.nex
#   Step 3 (R):        Rscript 03_tree_Bayesian.R
#                        <- reads consensus tree, saves RData + PDF
#
# Requires: out/rdata/alignment_objects.RData (run 00_prepare_alignment.R first)
#
# Outputs (after Step 3):
#   out/rdata/trees_Bayesian.RData
#   out/trees/tree_Bayesian.pdf
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


nex_path <- "MrBayes/gulf_fish_COI_aligned.nex"
con_path <- "MrBayes/gulf_fish_COI_aligned.nex.con.tre"

if (!file.exists(nex_path)) {
  stop("NEXUS file not found: ", nex_path, "\n",
       "Run 00_prepare_alignment.R first to generate it.")
}
cat("\nFound NEXUS file:", nex_path, "\n")


cat("\n--- STEP 1: Checking MrBayes block ---\n")

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


cat("\n--- STEP 2: Checking for MrBayes output ---\n")

if (!file.exists(con_path)) {
  cat("\n")
  cat("============================================================\n")
  cat("MrBayes output not found. Run MrBayes now:\n")
  cat("============================================================\n\n")
  cat("  In a terminal:\n")
  cat("    cd R/MrBayes\n")
  cat("    mb gulf_fish_COI_aligned.nex\n\n")
  cat("  After MrBayes finishes, re-run this script:\n")
  cat("    Rscript 03_tree_Bayesian.R\n")
  cat("============================================================\n")
  quit(save = "no", status = 0)
}

cat("  Found:", con_path, "\n")


cat("\n--- STEP 3: Reading MrBayes consensus tree ---\n")

tree_bayes_raw <- read.nexus(con_path)
cat("  Tips:", length(tree_bayes_raw$tip.label), "\n")


cat("\n--- STEP 4: Rooting on Bull Shark ---\n")

outgrp_bayes <- grep("Bull_Shark", tree_bayes_raw$tip.label, value = TRUE)

if (length(outgrp_bayes) == 0) {
  warning("Bull_Shark tips not found in consensus tree. ",
          "Check tip label format in the .con.tre file.")
  tree_bayes_rooted <- tree_bayes_raw
} else {
  tree_bayes_rooted <- root(tree_bayes_raw, outgroup = outgrp_bayes,
                             resolve.root = TRUE)
  cat("  Rooted on:", paste(outgrp_bayes, collapse = ", "), "\n")
}


cat("\n--- STEP 5: Saving out/rdata/trees_Bayesian.RData ---\n")

save(tree_bayes_rooted, file = "out/rdata/trees_Bayesian.RData")

cat("  Saved: out/rdata/trees_Bayesian.RData\n")


cat("\n--- STEP 6: Saving out/trees/tree_Bayesian.pdf ---\n")

pdf("out/trees/tree_Bayesian.pdf", width = 12, height = 9)

layout(
  matrix(c(1, 2), nrow = 2, byrow = TRUE),
  heights = c(7.5, 1.5)
)

par(mar = c(1, 1, 3, 1))
plot.phylo(
  tree_bayes_rooted,
  main       = "Bayesian Inference (GTR+\u0393, MrBayes)\nOutgroup: Bull Shark | Node labels = posterior probability",
  tip.color  = tip_color(tree_bayes_rooted$tip.label),
  cex        = 0.78,
  edge.width = 1.6,
  no.margin  = FALSE
)

raw_lines <- readLines(con_path)
tree_line <- raw_lines[grep("tree con_50_majrule", raw_lines)]

pp_raw <- regmatches(
  tree_line,
  gregexpr('prob\\(percent\\)="([0-9]+)"', tree_line)
)[[1]]

pp_vals <- as.numeric(gsub('prob\\(percent\\)="([0-9]+)"', "\\1", pp_raw))

cat("  Extracted", length(pp_vals), "posterior probability values\n")

n_nodes <- tree_bayes_rooted$Nnode

if (length(pp_vals) >= n_nodes) {
  pp_nodes  <- tail(pp_vals, n_nodes)
  pp_labels <- ifelse(pp_nodes >= 50, as.character(pp_nodes), "")

  nodelabels(pp_labels,
             adj   = c(1.1, -0.3),
             cex   = 0.55,
             frame = "none",
             col   = "gray25")
} else {
  cat("  Warning: could not match PP values to nodes — plotting without labels.\n")
}

add.scale.bar(x = 0, y = .5, cex = 0.7, col = "gray30")

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

dev.off()
cat("  Saved: out/trees/tree_Bayesian.pdf\n")
