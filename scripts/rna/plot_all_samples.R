#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (!length(args) %in% c(2, 3)) {
  stop("Usage: Rscript plot_all_samples.R <summary.tsv> <output_dir> [gene]")
}

input_file <- args[1]
output_dir <- args[2]
target_gene <- if (length(args) == 3) args[3] else NULL

required <- c(
  "sample_id", "gene_id", "transcript_id",
  "transcript_start", "transcript_end",
  "total_junctions", "coverage"
)

d <- read_tsv(
  input_file,
  show_col_types = FALSE,
  col_types = cols(
    .default = col_guess(),
    transcript_start = col_double(),
    transcript_end = col_double(),
    total_junctions = col_double(),
    coverage = col_double()
  )
)

missing <- setdiff(required, names(d))
if (length(missing)) {
  stop("Missing columns: ", paste(missing, collapse = ", "))
}

# Optional gene selection; do not rely on MSTRG identifiers.
if (!is.null(target_gene)) {
  by_name <- if ("gene_name" %in% names(d)) {
    !is.na(d$gene_name) & d$gene_name == target_gene
  } else {
    rep(FALSE, nrow(d))
  }

  d <- d[by_name | d$gene_id == target_gene, ]
}

# Include reconstructed structures even when coverage is zero.
# Zero coverage is not evidence of transcript expression.
d <- d %>%
  filter(
    !is.na(sample_id),
    !is.na(transcript_id),
    is.finite(coverage),
    is.finite(transcript_start),
    is.finite(transcript_end)
  ) %>%
  mutate(
    genomic_span = transcript_end - transcript_start + 1
  )

if (!nrow(d)) stop("No transcript models remain after filtering.")

if (anyDuplicated(d[c("sample_id", "transcript_id")])) {
  stop("Duplicate transcript IDs within a sample.")
}

if (any(!is.finite(d$genomic_span)) ||
    any(d$genomic_span <= 0) ||
    any(!is.finite(d$total_junctions))) {
  stop("Invalid transcript coordinates or junction counts.")
}

samples <- sort(unique(d$sample_id))
d$sample_id <- factor(d$sample_id, levels = samples)

base_colours <- c(
  "#24476D", "#159A99", "#8C6AA6",
  "#D49B24", "#A95B60", "#637D88"
)

palette <- if (length(samples) <= length(base_colours)) {
  setNames(base_colours[seq_along(samples)], samples)
} else {
  setNames(grDevices::hcl.colors(length(samples), "Dark 3"), samples)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

format_bp <- function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

theme_report <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    panel.spacing = grid::unit(1.2, "cm"),
    plot.margin = margin(12, 18, 12, 12)
  )

save_figure <- function(p, filename, width = 12, height = 7) {
  ggsave(
    file.path(output_dir, paste0(filename, ".png")),
    plot = p, width = width, height = height,
    dpi = 300, bg = "white"
  )

  ggsave(
    file.path(output_dir, paste0(filename, ".pdf")),
    plot = p, width = width, height = height,
    device = cairo_pdf
  )
}

# Figure 1: Genomic-span distribution across samples.
p_span <- ggplot(
  d, aes(x = genomic_span, fill = sample_id)
) +
  geom_histogram(bins = 50, colour = "white", linewidth = 0.2) +
  facet_wrap(
    ~sample_id,
    ncol = min(2, length(samples)),
    scales = "free_y"
  ) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(labels = format_bp) +
  labs(
    title = "Transcript genomic-span distribution",
    x = "Genomic span (bp)",
    y = "Number of transcript models"
  ) +
  theme_report

save_figure(p_span, "01_genomic_span")

# Figure 2: Splice-junction distribution across samples.
p_junctions <- ggplot(
  d, aes(x = factor(total_junctions), fill = sample_id)
) +
  geom_bar(width = 0.7) +
  facet_wrap(
    ~sample_id,
    ncol = min(2, length(samples)),
    scales = "free_y"
  ) +
  scale_fill_manual(values = palette) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(
    title = "Splice-junction distribution",
    x = "Number of splice junctions",
    y = "Number of transcript models"
  ) +
  theme_report

save_figure(p_junctions, "02_junction_counts")

# Figure 3: Transcript start and end positions.
positions <- d %>%
  select(sample_id, transcript_start, transcript_end) %>%
  pivot_longer(
    cols = c(transcript_start, transcript_end),
    names_to = "position_type",
    values_to = "position"
  ) %>%
  mutate(
    position_type = factor(
      position_type,
      levels = c("transcript_start", "transcript_end"),
      labels = c("Start", "End")
    )
  )

p_boundaries <- ggplot(
  positions, aes(x = position, fill = sample_id)
) +
  geom_histogram(bins = 60, colour = "white", linewidth = 0.2) +
  facet_grid(
    position_type ~ sample_id,
    scales = "free_y"
  ) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(labels = format_bp) +
  labs(
    title = "Transcript start and end positions",
    x = "Genomic coordinate (bp)",
    y = "Number of transcript models"
  ) +
  theme_report +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

save_figure(
  p_boundaries, "03_transcript_boundaries",
  width = max(12, 3.8 * length(samples)), height = 8
)

# Figure 4: Genomic span versus splice-junction count.
# Individual transcript models shown because many groups are sparse.

p_span_junctions <- ggplot(
  d,
  aes(
    x = total_junctions,
    y = genomic_span,
    colour = sample_id
  )
) +
  geom_point(
    position = position_jitter(
      width = 0.16,
      height = 0,
      seed = 42
    ),
    size = 2.6,
    alpha = 0.8
  ) +
  facet_wrap(~sample_id, nrow = 1) +
  scale_colour_manual(values = palette) +
  scale_x_continuous(
    breaks = sort(unique(d$total_junctions)),
    guide = guide_axis(n.dodge = 2)
  ) +
  scale_y_continuous(labels = format_bp) +
  labs(
    title = "Genomic span by splice-junction count",
    x = "Number of splice junctions",
    y = "Genomic span (bp)"
  ) +
  theme_report +
  theme(legend.position = "none")

save_figure(p_span_junctions, "04_span_vs_junctions")


cat("Multi-sample figures generated successfully.\n")
cat("Samples:", paste(samples, collapse = ", "), "\n")
print(d %>% count(sample_id, name = "transcript_models"))
