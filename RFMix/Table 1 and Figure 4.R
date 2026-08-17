# ============================================================
# Gene Ontology enrichment analysis
# Honey bee ancestry-associated genomic regions
#
# Ancestries: A, C, M, O
#
# Candidate genes were identified from genomic regions
# meeting the 2-SD regional SNP/ancestry criterion.
#
# GO enrichment:
#   Ontology   = Biological Process (BP)
#   Algorithm  = weight01
#   Statistic  = Fisher's exact test
#   nodeSize   = 10
#
# Outputs:
#   1. GO_BP_A_results.tsv
#   2. GO_BP_C_results.tsv
#   3. GO_BP_M_results.tsv
#   4. GO_BP_O_results.tsv
#   5. Table1_GO_enrichment.csv
#   6. Figure4A_GO_BP_A.tiff
#
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(topGO)
library(dplyr)
library(stringr)
library(ggplot2)


# ------------------------------------------------------------
# 2. Load GO mapping
# ------------------------------------------------------------

load("GOmapping_bee.RData")

# GOmapping.bee contains the honey bee gene-to-GO mapping
geneUniverse <- names(GOmapping.bee)

cat("GO universe:", length(geneUniverse), "genes\n")


# ------------------------------------------------------------
# 3. Candidate-gene files
# ------------------------------------------------------------

candidate.files <- c(
  A = "candidate_A_genes.txt",
  C = "candidate_C_genes.txt",
  M = "candidate_M_genes.txt",
  O = "candidate_O_genes.txt"
)


# ------------------------------------------------------------
# 4. Output directory
# ------------------------------------------------------------

dir.create(
  "GO_results",
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 5. Run GO enrichment for A, C, M, and O
# ------------------------------------------------------------

go_results <- list()

for (anc in names(candidate.files)) {
  
  cat("\n========================================\n")
  cat("Processing ancestry:", anc, "\n")
  cat("========================================\n")
  
  
  # ----------------------------------------------------------
  # Read candidate genes
  # ----------------------------------------------------------
  
  candidate_genes <- readLines(
    candidate.files[anc]
  )
  
  candidate_genes <- unique(
    na.omit(candidate_genes)
  )
  
  cat(
    "Candidate genes before GO filtering:",
    length(candidate_genes),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # Match candidate genes to GO universe
  # ----------------------------------------------------------
  
  candidate_genes_GO <- intersect(
    candidate_genes,
    geneUniverse
  )
  
  cat(
    "GO-matched genes:",
    length(candidate_genes_GO),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # Create gene list
  # ----------------------------------------------------------
  
  geneList <- factor(
    as.integer(
      geneUniverse %in% candidate_genes_GO
    )
  )
  
  names(geneList) <- geneUniverse
  
  
  # ----------------------------------------------------------
  # Create topGO object
  # ----------------------------------------------------------
  
  GOdata <- new(
    "topGOdata",
    ontology = "BP",
    allGenes = geneList,
    annot = annFUN.gene2GO,
    gene2GO = GOmapping.bee,
    nodeSize = 10
  )
  
  
  # ----------------------------------------------------------
  # Fisher enrichment test
  # ----------------------------------------------------------
  
  result <- runTest(
    GOdata,
    algorithm = "weight01",
    statistic = "fisher"
  )
  
  
  # ----------------------------------------------------------
  # Generate GO results
  # ----------------------------------------------------------
  
  res <- GenTable(
    GOdata,
    weightFisher = result,
    topNodes = 100
  )
  
  
  # ----------------------------------------------------------
  # Calculate FDR, logP, and gene ratio
  # ----------------------------------------------------------
  
  res <- res %>%
    mutate(
      weightFisher = as.numeric(
        gsub("< ", "", weightFisher)
      ),
      FDR = p.adjust(
        weightFisher,
        method = "BH"
      ),
      logP = -log10(weightFisher),
      GeneRatio = Significant / Annotated
    )
  
  
  # ----------------------------------------------------------
  # Save complete GO result table
  # ----------------------------------------------------------
  
  output_file <- file.path(
    "GO_results",
    paste0("GO_BP_", anc, "_results.tsv")
  )
  
  write.table(
    res,
    output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  cat(
    "Saved:",
    output_file,
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # Print significant GO terms
  # ----------------------------------------------------------
  
  significant <- res %>%
    filter(FDR < 0.05) %>%
    arrange(FDR)
  
  cat(
    "GO terms with FDR < 0.05:",
    nrow(significant),
    "\n"
  )
  
  if (nrow(significant) > 0) {
    
    print(
      significant %>%
        dplyr::select(
          GO.ID,
          Term,
          Annotated,
          Significant,
          Expected,
          weightFisher,
          FDR,
          GeneRatio
        )
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Store results
  # ----------------------------------------------------------
  
  go_results[[anc]] <- list(
    results = res,
    candidate_genes = candidate_genes,
    candidate_genes_GO = candidate_genes_GO,
    GOdata = GOdata,
    result = result
  )
}


# ============================================================
# 6. Create Table 1
# ============================================================

# Table 1 contains significant terms for A, C, M, and O.
#
# If A has more significant terms than desired for the
# manuscript table, the complete A results remain available
# in GO_BP_A_results.tsv.


table1 <- bind_rows(
  lapply(
    names(go_results),
    function(anc) {
      
      res <- go_results[[anc]]$results
      
      candidate_n <- length(
        go_results[[anc]]$candidate_genes
      )
      
      go_n <- length(
        go_results[[anc]]$candidate_genes_GO
      )
      
      res %>%
        filter(FDR < 0.05) %>%
        arrange(FDR) %>%
        mutate(
          Ancestry = anc,
          Candidate_genes = candidate_n,
          GO_matched_genes = go_n
        ) %>%
        dplyr::select(
          Ancestry,
          Candidate_genes,
          GO_matched_genes,
          GO.ID,
          Term,
          Annotated,
          Significant,
          GeneRatio,
          weightFisher,
          FDR
        )
    }
  )
)


# ------------------------------------------------------------
# Rename columns for manuscript table
# ------------------------------------------------------------

table1 <- table1 %>%
  rename(
    `Ancestry` = Ancestry,
    `Candidate genes` = Candidate_genes,
    `GO-matched genes` = GO_matched_genes,
    `GO ID` = GO.ID,
    `GO term` = Term,
    `Annotated` = Annotated,
    `Significant` = Significant,
    `Gene ratio` = GeneRatio,
    `P` = weightFisher,
    `FDR` = FDR
  )


# ------------------------------------------------------------
# Round numerical values
# ------------------------------------------------------------

table1 <- table1 %>%
  mutate(
    `Gene ratio` = round(`Gene ratio`, 3),
    P = signif(P, 3),
    FDR = signif(FDR, 3)
  )


# ------------------------------------------------------------
# Save Table 1
# ------------------------------------------------------------

write.csv(
  table1,
  "GO_results/Table1_GO_enrichment.csv",
  row.names = FALSE
)

cat(
  "\nSaved: GO_results/Table1_GO_enrichment.csv\n"
)


# ============================================================
# 7. Figure 4A — A-lineage GO enrichment
# ============================================================

A_results_final <- go_results$A$results


# ------------------------------------------------------------
# Select significant A terms
# ------------------------------------------------------------

plot_A <- A_results_final %>%
  filter(FDR < 0.05) %>%
  arrange(FDR)


# ------------------------------------------------------------
# Order terms by significance
# ------------------------------------------------------------

plot_A <- plot_A %>%
  mutate(
    Term = factor(
      Term,
      levels = rev(Term)
    )
  )


# ------------------------------------------------------------
# Create Figure 4A
# ------------------------------------------------------------

Figure4A <- ggplot(
  plot_A,
  aes(
    x = GeneRatio,
    y = Term
  )
) +
  geom_point(
    aes(
      size = Significant,
      color = logP
    )
  ) +
  scale_color_viridis_c(
    option = "plasma"
  ) +
  theme_bw(
    base_size = 14
  ) +
  labs(
    x = "Gene Ratio",
    y = NULL,
    color = expression(
      -log[10](italic(P))
    ),
    size = "Genes"
  ) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(
      size = 11,
      color = "black"
    ),
    axis.text.x = element_text(
      size = 12,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 13,
      color = "black"
    )
  )


# ------------------------------------------------------------
# Display figure
# ------------------------------------------------------------

Figure4A


# ------------------------------------------------------------
# Save Figure 4A
# ------------------------------------------------------------

ggsave(
  filename = "GO_results/Figure4A_GO_BP_A.tiff",
  plot = Figure4A,
  device = "tiff",
  width = 8,
  height = 7,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

cat(
  "Saved: GO_results/Figure4A_GO_BP_A.tiff\n"
)


# ============================================================
# 8. Print final summary
# ============================================================

cat("\n========================================\n")
cat("GO ANALYSIS COMPLETED\n")
cat("========================================\n")

for (anc in names(go_results)) {
  
  candidate_n <- length(
    go_results[[anc]]$candidate_genes
  )
  
  go_n <- length(
    go_results[[anc]]$candidate_genes_GO
  )
  
  n_sig <- sum(
    go_results[[anc]]$results$FDR < 0.05,
    na.rm = TRUE
  )
  
  cat(
    anc,
    ":",
    candidate_n,
    "candidate genes;",
    go_n,
    "GO-matched genes;",
    n_sig,
    "significant GO terms\n"
  )
}