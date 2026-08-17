
library(topGO)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)

# ----------------------------------------------------------
# GO mapping
# ----------------------------------------------------------

load("/Users/mulusewfikere/Downloads/GOmapping_bee.RData")

geneUniverse <- names(GOmapping.bee)

# ----------------------------------------------------------
# Peak overlap files
# ----------------------------------------------------------

peak.files <- c(
  C = "peak_C_gene_overlap_chr.txt",
  M = "peak_M_gene_overlap_chr.txt",
  O = "peak_O_gene_overlap_chr.txt"
)

# ----------------------------------------------------------
# Output storage
# ----------------------------------------------------------

GO_results <- list()
plot_list <- list()

# ----------------------------------------------------------
# Run GO analysis for C, M, O
# ----------------------------------------------------------

for (anc in names(peak.files)) {
  
  cat("\n========================================\n")
  cat("Processing ancestry:", anc, "\n")
  cat("========================================\n")
  
  # --------------------------------------------------------
  # Read overlap file
  # --------------------------------------------------------
  
  annot <- read.delim(
    peak.files[anc],
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = ""
  )
  
  # --------------------------------------------------------
  # Extract gene IDs
  # --------------------------------------------------------
  
  candidate_genes <- str_extract(
    annot$V12,
    'gene_id "[^"]+"'
  )
  
  candidate_genes <- str_remove(
    candidate_genes,
    'gene_id "'
  )
  
  candidate_genes <- str_remove(
    candidate_genes,
    '"'
  )
  
  candidate_genes <- unique(
    na.omit(candidate_genes)
  )
  
  cat(
    "Candidate genes before GO filtering:",
    length(candidate_genes),
    "\n"
  )
  
  # --------------------------------------------------------
  # Keep only genes in GO universe
  # --------------------------------------------------------
  
  candidate_genes <- intersect(
    candidate_genes,
    geneUniverse
  )
  
  cat(
    "Matched GO genes:",
    length(candidate_genes),
    "\n"
  )
  
  # --------------------------------------------------------
  # Skip if no genes
  # --------------------------------------------------------
  
  if (length(candidate_genes) == 0) {
    
    cat(
      "No genes matched the GO universe. Skipping",
      anc,
      "\n"
    )
    
    next
  }
  
  # --------------------------------------------------------
  # Gene list
  # --------------------------------------------------------
  
  geneList <- factor(
    as.integer(
      geneUniverse %in% candidate_genes
    )
  )
  
  names(geneList) <- geneUniverse
  
  # --------------------------------------------------------
  # topGO object
  # --------------------------------------------------------
  
  GOdata <- new(
    "topGOdata",
    ontology = "BP",
    allGenes = geneList,
    annot = annFUN.gene2GO,
    gene2GO = GOmapping.bee,
    nodeSize = 10
  )
  
  # --------------------------------------------------------
  # Fisher enrichment
  # --------------------------------------------------------
  
  result <- runTest(
    GOdata,
    algorithm = "weight01",
    statistic = "fisher"
  )
  
  # --------------------------------------------------------
  # GO results
  # --------------------------------------------------------
  
  res <- GenTable(
    GOdata,
    weightFisher = result,
    topNodes = 100
  )
  
  # --------------------------------------------------------
  # Statistics
  # --------------------------------------------------------
  
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
    ) %>%
    arrange(weightFisher)
  
  # --------------------------------------------------------
  # Save complete results
  # --------------------------------------------------------
  
  write.table(
    res,
    paste0(
      "GO_BP_",
      anc,
      "_results.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  GO_results[[anc]] <- res
  
  # --------------------------------------------------------
  # Print summary
  # --------------------------------------------------------
  
  cat(
    "GO terms tested:",
    nrow(res),
    "\n"
  )
  
  cat(
    "Minimum nominal P:",
    min(res$weightFisher, na.rm = TRUE),
    "\n"
  )
  
  cat(
    "Minimum FDR:",
    min(res$FDR, na.rm = TRUE),
    "\n"
  )
  
  cat(
    "Terms with FDR < 0.05:",
    sum(res$FDR < 0.05, na.rm = TRUE),
    "\n"
  )
  
  # --------------------------------------------------------
  # Top 20 terms for plotting
  #
  # If significant terms exist, plot those.
  # Otherwise plot top 20 by nominal P.
  # --------------------------------------------------------
  
  significant_df <- res %>%
    filter(FDR < 0.05)
  
  if (nrow(significant_df) > 0) {
    
    plot_df <- significant_df %>%
      arrange(weightFisher) %>%
      slice_head(n = 20)
    
    plot_type <- "FDR < 0.05"
    
  } else {
    
    plot_df <- res %>%
      arrange(weightFisher) %>%
      slice_head(n = 20)
    
    plot_type <- "Top 20 nominal P-values"
  }
  
  # --------------------------------------------------------
  # Make term order
  # --------------------------------------------------------
  
  plot_df <- plot_df %>%
    mutate(
      Term = factor(
        Term,
        levels = rev(Term)
      )
    )
  
  # --------------------------------------------------------
  # Dot plot
  # --------------------------------------------------------
  
  p <- ggplot(
    plot_df,
    aes(
      x = GeneRatio,
      y = Term
    )
  ) +
    geom_point(
      aes(
        size = Significant,
        color = -log10(weightFisher)
      )
    ) +
    scale_color_viridis_c(
      option = "plasma"
    ) +
    labs(
      x = "Gene Ratio",
      y = NULL,
      color = expression(
        -log[10](italic(P))
      ),
      size = "Genes"
    ) +
    theme_bw(
      base_size = 14
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
  
  plot_list[[anc]] <- p
  
  # --------------------------------------------------------
  # Save individual figure
  # --------------------------------------------------------
  
  ggsave(
    filename = paste0(
      "GO_BP_dotplot_",
      anc,
      ".tiff"
    ),
    plot = p,
    device = "tiff",
    width = 8,
    height = 7,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )
  
  cat(
    "Saved GO_BP_dotplot_",
    anc,
    ".tiff\n"
  )
  
  # --------------------------------------------------------
  # Print top 20
  # --------------------------------------------------------
  
  cat(
    "\nTop GO terms for",
    anc,
    ":\n"
  )
  
  print(
    plot_df[
      c(
        "GO.ID",
        "Term",
        "Annotated",
        "Significant",
        "weightFisher",
        "FDR",
        "GeneRatio"
      )
    ]
  )
}

# ----------------------------------------------------------
# Combined C/M/O figure
# ----------------------------------------------------------

combined_plot_CMO <-
  plot_list$C +
  plot_list$M +
  plot_list$O +
  plot_layout(
    ncol = 1
  )

ggsave(
  "GO_BP_dotplots_CMO.tiff",
  combined_plot_CMO,
  device = "tiff",
  width = 10,
  height = 12,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

cat("\n========================================\n")
cat("Analysis completed.\n")
cat("========================================\n")