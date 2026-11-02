#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (!length(args) %in% c(3, 4)) {
  stop(paste(
    "Usage: Rscript plot_single_sample.R",
    "<summary.tsv> <sample_id> <output_dir> [MANE_transcript_id]"
  ))
}

input_file <- args[1]
sample_name <- args[2]
output_dir <- args[3]

mane_id <- if (length(args) == 4) {
  sub("\\.[0-9]+$", "", args[4])
} else {
  NA_character_
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required <- c(
  "sample_id", "transcript_id", "spliced_length",
  "genomic_span", "total_junctions", "coverage",
  "transcript_start", "transcript_end"
)

data <- read_tsv(input_file, show_col_types = FALSE)

missing <- setdiff(required, names(data))
if (length(missing)) {
  stop("Missing columns: ", paste(missing, collapse = ", "))
}

data <- data %>%
  filter(sample_id == sample_name) %>%
  mutate(
    transcript_class = case_when(
      !is.na(mane_id) & transcript_id == mane_id ~ "MANE Select",
      grepl("^MSTRG\\.|^Novel\\.", transcript_id) ~
        "Other transcript models",
      grepl("^ENST", transcript_id) ~ "Reference-identified transcripts",
      TRUE ~ "Other"
    ),
    transcript_class = factor(
      transcript_class,
      levels = c(
        "Reference-identified transcripts",
        "Other transcript models",
        "MANE Select",
        "Other"
      )
    )
  )

if (!nrow(data)) {
  stop("No transcript models found for ", sample_name)
}

if (!any(is.finite(data$coverage) & data$coverage > 0)) {
  message(
    "No positive-coverage transcript models for ", sample_name,
    ". Figures describe reconstructed structures, not confirmed expression."
  )
}

if (anyDuplicated(data$transcript_id)) {
  stop("Duplicate transcript IDs within sample: ", sample_name)
}

colours <- c(
  "Reference-identified transcripts" = "#24476D",
  "Other transcript models" = "#159A99",
  "MANE Select" = "#D49B24",
  "Other" = "#888888"
)

theme_report <- function() {
  theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(
        face = "bold", size = 16,
        margin = margin(b = 8)
      ),
      plot.subtitle = element_text(
        colour = "#596575", size = 10
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.title = element_text(colour = "#263546"),
      axis.text = element_text(colour = "#465465"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box.just = "center",
      legend.title = element_blank(),
      plot.margin = margin(15, 20, 15, 15)
    )
}

save_plot <- function(p, filename, width = 11, height = 7) {
  ggsave(
    filename = file.path(output_dir, filename),
    plot = p,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )

  ggsave(
    filename = file.path(
      output_dir,
      sub("\\.png$", ".pdf", filename)
    ),
    plot = p,
    width = width,
    height = height,
    device = cairo_pdf
  )
}

# 1. Ranked spliced transcript lengths

ranked <- data %>%
  arrange(spliced_length, transcript_id) %>%
  mutate(
    transcript_label = factor(
      transcript_id,
      levels = transcript_id
    )
  )

p_lengths <- ggplot(
  ranked,
  aes(x = spliced_length, y = transcript_label)
) +
  geom_segment(
    aes(x = 0, xend = spliced_length, yend = transcript_label),
    colour = "#DCE3E9",
    linewidth = 0.6
  ) +
  geom_point(
    aes(colour = transcript_class),
    size = 3
  ) +
  scale_colour_manual(values = colours, drop = TRUE) +
  scale_x_continuous(
    labels = function(x) format(x, big.mark = ",", trim = TRUE),
    expand = expansion(mult = c(0, 0.06))
  ) +
  labs(
    title = "Transcript length profiles",
    x = "Spliced transcript length (bp)",
    y = NULL
  ) +
  theme_report()

save_plot(
  p_lengths,
  "transcript_length.png",
  height = max(8, nrow(ranked) * 0.38)
)

# 2. Junction-count frequency

junction_counts <- data %>%
  count(total_junctions, name = "n")

p_junctions <- ggplot(
  junction_counts,
  aes(x = total_junctions, y = n)
) +
  geom_segment(
    aes(xend = total_junctions, yend = 0),
    linewidth = 1.1,
    colour = "#BACBD8"
  ) +
  geom_point(
    colour = "#24476D",
    size = 4
  ) +
  scale_x_continuous(
    breaks = sort(unique(junction_counts$total_junctions)),
    guide = guide_axis(n.dodge = 3)
  ) +
  scale_y_continuous(
    breaks = function(limits) pretty(limits, n = 6),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = "Splice-junction complexity",
    x = "Splice junctions per transcript",
    y = "Number of transcript models"
  ) +
  theme_report()

if (!is.na(mane_id) && any(data$transcript_id == mane_id)) {
  mane_junctions <- data$total_junctions[
    match(mane_id, data$transcript_id)
  ]

  p_junctions <- p_junctions +
    geom_vline(
      data = data.frame(ref_junctions = mane_junctions),
      aes(xintercept = ref_junctions, linetype = "MANE Select"),
      inherit.aes = FALSE,
      colour = "#D49B24",
      linewidth = 0.8,
      show.legend = TRUE
    ) +
    scale_linetype_manual(
      name = NULL,
      values = c("MANE Select" = "dashed")
    ) +
    guides(
      linetype = guide_legend(
        override.aes = list(colour = "#D49B24")
      )
    )
}

save_plot(p_junctions, "junction_counts.png")

# 3. Genomic interval plot

intervals <- data %>%
  arrange(transcript_start, transcript_end, transcript_id) %>%
  mutate(
    transcript_label = factor(
      transcript_id,
      levels = rev(transcript_id)
    )
  )

p_boundaries <- ggplot(
  intervals,
  aes(
    y = transcript_label,
    colour = transcript_class
  )
) +
  geom_segment(
    aes(
      x = transcript_start,
      xend = transcript_end,
      yend = transcript_label
    ),
    linewidth = 1.1,
    lineend = "round"
  ) +
  geom_point(
    aes(x = transcript_start),
    size = 1.5
  ) +
  geom_point(
    aes(x = transcript_end),
    size = 1.5
  ) +
  scale_colour_manual(values = colours, drop = TRUE) +
  scale_x_continuous(
    labels = function(x) format(x, big.mark = ",", trim = TRUE)
  ) +
  labs(
    title = "Genomic extent of transcript models",
    x = "Genomic coordinate (bp)",
    y = NULL
  ) +
  theme_report()

save_plot(
  p_boundaries,
  "transcript_boundaries.png",
  height = max(8, nrow(intervals) * 0.38)
)

# 4. Transcript structural complexity

p_complexity <- ggplot(
  data,
  aes(
    x = total_junctions,
    y = spliced_length,
    colour = transcript_class,
    size = coverage
  )
) +
  geom_point(alpha = 0.8) +
  scale_colour_manual(values = colours, drop = TRUE) +
  scale_size_continuous(range = c(2, 8)) +
  scale_y_continuous(
    labels = function(x) format(x, big.mark = ",", trim = TRUE)
  ) +
  labs(
    title = "Transcript structural complexity",
    x = "Number of splice junctions",
    y = "Spliced transcript length (bp)",
    colour = "Transcript class",
    size = "Coverage"
  ) +
  guides(
    colour = guide_legend(order = 1, nrow = 1),
    size = guide_legend(order = 2, nrow = 1)
  ) +
  theme_report() +
  theme(
    legend.box = "vertical",
    legend.box.just = "center",
    legend.spacing.y = grid::unit(10, "pt"),
    legend.title = element_text(face = "bold")
  )

save_plot(
  p_complexity,
  "transcript_complexity.png"
)

# Original bubble plot retained; approved Figure 1 is generated separately.

cat("Sample:", sample_name, "\n")
cat("Transcript models:", nrow(data), "\n")
cat("MANE ID:", mane_id, "\n")
cat("Figures generated: 4\n")
cat("Output directory:", output_dir, "\n")
