# ============================================================
# Figure 2: Latitude, Ancestry Proportion, and Haplotype Length
# ============================================================

# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(grid)


# ------------------------------------------------------------
# 2. File paths
# ------------------------------------------------------------
# Run this script from the project root directory.
#
# Expected project structure:
#
# project/
# ├── data/
# │   ├── metadata_USonly_2APR25.csv
# │   └── ibd/
# │       └── tracts_with_pos_df.csv
# ├── output/
# └── scripts/
#
# Update these paths if your repository uses a different
# directory structure.

metadata_file <- "data/metadata_USonly_2APR25.csv"
tracts_file   <- "data/ibd/tracts_with_pos_df.csv"
output_dir    <- "output"


# Create output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# ------------------------------------------------------------
# 3. Define lineage colors
# ------------------------------------------------------------

lineage_colors <- c(
  A = "#4daf4a",
  C = "#377eb8",
  M = "#e41a1c",
  O = "gold"
)


# ------------------------------------------------------------
# 4. Read metadata
# ------------------------------------------------------------

meta <- read.csv(metadata_file, header = TRUE)

meta_sub <- meta %>%
  select(
    vcf_id,
    individual_id,
    state,
    managed,
    lat,
    lon,
    mito
  ) %>%
  filter(complete.cases(.)) %>%
  rename(sample = vcf_id)


# ------------------------------------------------------------
# 5. A. Global ancestry proportion vs. latitude
# ------------------------------------------------------------
# `mean_rfmix_ACMO` is assumed to have already been created
# from the RFMix analysis and contain:
# sample, A, C, M, O

mean_rfmix_ACMO <- mean_rfmix_ACMO %>%
  select(sample, A, C, M, O)

mean_rfmix_lat <- mean_rfmix_ACMO %>%
  left_join(meta_sub, by = "sample")


# Convert ancestry proportions to long format
mean_rfmix_long <- mean_rfmix_lat %>%
  pivot_longer(
    cols = c(A, C, M, O),
    names_to = "ancestry",
    values_to = "proportion"
  )


# ------------------------------------------------------------
# 6. Regression statistics for ancestry proportion
# ------------------------------------------------------------

ancestry_results <- mean_rfmix_long %>%
  group_by(ancestry) %>%
  summarise(
    model = list(lm(proportion ~ lat, data = cur_data())),
    .groups = "drop"
  ) %>%
  mutate(
    slope = sapply(model, function(x) coef(x)[["lat"]]),
    R2 = sapply(model, function(x) summary(x)$r.squared),
    r = sapply(model, function(x) {
      cor(
        model.frame(x)$proportion,
        model.frame(x)$lat,
        use = "complete.obs"
      )
    }),
    p_value = sapply(model, function(x) {
      coef(summary(x))["lat", "Pr(>|t|)"]
    }),
    f2 = R2 / (1 - R2)
  ) %>%
  select(-model)

print(ancestry_results)


# ------------------------------------------------------------
# 7. Plot A: Ancestry proportion vs. latitude
# ------------------------------------------------------------

plot_A <- ggplot(
  mean_rfmix_long,
  aes(x = lat, y = proportion, color = ancestry)
) +
  geom_point(
    alpha = 0.7,
    size = 3
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    linewidth = 1
  ) +
  facet_wrap(
    ~ ancestry,
    scales = "free_y"
  ) +
  scale_color_manual(
    values = lineage_colors
  ) +
  labs(
    title = "A",
    x = NULL,
    y = "Ancestral Proportion"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    legend.key.width = unit(1, "cm"),
    legend.title = element_text(
      size = 14,
      color = "black"
    ),
    legend.text = element_text(
      size = 14,
      color = "black"
    ),
    strip.text = element_text(
      size = 16,
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold",
      size = 16,
      color = "black"
    ),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(
      size = 14,
      color = "black"
    )
  )


# ------------------------------------------------------------
# 8. B. Haplotype length vs. latitude
# ------------------------------------------------------------

tracts <- read.csv(
  tracts_file,
  header = TRUE
)


# Summarise mean tract length within each individual and
# ancestry lineage
tracts_summary <- tracts %>%
  group_by(Individual, Ancestry) %>%
  summarise(
    mean_length = mean(Length, na.rm = TRUE),
    .groups = "drop"
  )


# Add metadata
ancestry_long_lat <- tracts_summary %>%
  left_join(
    meta_sub %>%
      rename(Individual = individual_id),
    by = "Individual"
  )


# Keep the four focal lineages
haplotype_data <- ancestry_long_lat %>%
  filter(
    Ancestry %in% c("A", "C", "M", "O")
  )


# ------------------------------------------------------------
# 9. Plot B: Mean haplotype length vs. latitude
# ------------------------------------------------------------

plot_B <- ggplot(
  haplotype_data,
  aes(
    x = lat,
    y = mean_length,
    color = Ancestry
  )
) +
  geom_point(
    alpha = 0.7,
    size = 3
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    linewidth = 1
  ) +
  facet_wrap(
    ~ Ancestry,
    scales = "free_y"
  ) +
  scale_color_manual(
    values = lineage_colors
  ) +
  labs(
    title = "B",
    x = "Latitude",
    y = "Average Ancestry Haplotype Length"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(
      size = 16,
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold",
      size = 16,
      color = "black"
    )
  )


# ------------------------------------------------------------
# 10. Combine panels
# ------------------------------------------------------------

figure_2 <- plot_A / plot_B


# ------------------------------------------------------------
# 11. Save figure
# ------------------------------------------------------------

ggsave(
  filename = file.path(output_dir, "figure_2_AB.pdf"),
  plot = figure_2,
  width = 9,
  height = 11,
  units = "in"
)

ggsave(
  filename = file.path(output_dir, "figure_2_AB.png"),
  plot = figure_2,
  width = 9,
  height = 11,
  units = "in",
  dpi = 600
)


# ------------------------------------------------------------
# 12. Display figure
# ------------------------------------------------------------

figure_2