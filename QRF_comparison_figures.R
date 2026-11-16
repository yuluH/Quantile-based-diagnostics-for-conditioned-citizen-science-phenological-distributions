# ============================================================
# Cross-species QRF evaluation figures
# ============================================================
# This script is independent of QRF_phenology_analysis.R.
# MPIW-ICP, pinball loss and MAE figures use manually compiled
# summary statistics from the manuscript analyses. These values
# are not recomputed by this script.
# The optional DOY figure reads Hc.csv, Zm.csv, Fl.csv, Mf.csv
# from the current working directory.
# Run from the project root in R or RStudio.

required_packages <- c("ggplot2", "dplyr", "patchwork")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install the required packages before running: ",
       paste(missing_packages, collapse = ", "))
}
invisible(lapply(required_packages, library, character.only = TRUE))

# Final comparison figures (manually compiled cross-species summaries).

# ======================================================
# Figure: MPIW-ICP trajectory under RCV and SCV
# Nature-style formatting
# ======================================================

base_font <- "sans"

theme_nature <- theme_bw(base_size = 18, base_family = base_font) +
  theme(
    text = element_text(family = base_font),

    panel.grid.minor = element_blank(),

    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(size = 20, face = "bold"),

    axis.title.x = element_text(size = 21),
    axis.title.y = element_text(size = 21),
    axis.text.x  = element_text(size = 18),
    axis.text.y  = element_text(size = 18),

    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 17),
    legend.text  = element_text(size = 15),
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 10,
      l = 50
    )
  )

df_traj <- data.frame(
  Species = rep(
    c("H. chinensis", "Z. moltrechti", "F. limnocharis", "M. fissipes"),
    each = 6
  ),
  Model = rep(
    c("10 CSCIs", "VIP > 1.0 CSCIs"),
    times = 4,
    each = 3
  ),
  Stage = rep(c("In-sample", "RCV", "SCV"), times = 8),
  MPIW = c(
    7.268, 16.73, 18.45,
    15.729, 19.19, 20.27,
    13.702, 30.93, 35.98,
    34.164, 38.55, 42.94,
    7.988, 13.42, 15.78,
    15.970, 17.56, 16.64,
    9.970, 18.66, 22.31,
    20.728, 22.59, 24.78
  ),
  ICP = c(
    1.000, 0.762, 0.458,
    0.959, 0.713, 0.557,
    1.000, 0.718, 0.628,
    0.835, 0.540, 0.496,
    1.000, 0.676, 0.573,
    0.850, 0.570, 0.493,
    1.000, 0.693, 0.585,
    0.907, 0.587, 0.539
  )
)

df_traj <- df_traj %>%
  mutate(
    Model = recode(
      Model,
      "10 CSCIs" = "bold(a)",
      "VIP > 1.0 CSCIs" = "bold(b)"
    ),
    Model = factor(Model, levels = c("bold(a)", "bold(b)")),
    Species_label = case_when(
      Species == "H. chinensis" ~ "italic(H.~chinensis)",
      Species == "Z. moltrechti" ~ "italic(Z.~moltrechti)",
      Species == "F. limnocharis" ~ "italic(F.~limnocharis)",
      Species == "M. fissipes" ~ "italic(M.~fissipes)"
    )
  )

df_lines <- df_traj %>%
  filter(Stage %in% c("RCV", "SCV")) %>%
  left_join(
    df_traj %>%
      filter(Stage == "In-sample") %>%
      select(Species, Model, MPIW_in = MPIW, ICP_in = ICP),
    by = c("Species", "Model")
  ) %>%
  mutate(
    Validation = factor(Stage, levels = c("RCV", "SCV"))
  )

df_traj$Stage <- factor(
  df_traj$Stage,
  levels = c("In-sample", "RCV", "SCV")
)

p_traj <- ggplot() +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    color = "gray45",
    linewidth = 0.8
  ) +
  geom_segment(
    data = df_lines,
    aes(
      x = MPIW_in,
      y = ICP_in,
      xend = MPIW,
      yend = ICP,
      color = Species_label,
      linetype = Validation
    ),
    linewidth = 1.2,
    alpha = 0.85
  ) +
  geom_point(
    data = df_traj,
    aes(
      x = MPIW,
      y = ICP,
      color = Species_label,
      shape = Stage
    ),
    size = 4.5,
    stroke = 1.2
  ) +
  facet_wrap(
    ~ Model,
    scales = "free_x",
    ncol = 2,
    labeller = label_parsed
  ) +
  scale_color_discrete(
    labels = function(x) parse(text = x)
  ) +
  scale_linetype_manual(
    values = c("RCV" = "dashed", "SCV" = "dotdash")
  ) +
  scale_shape_manual(
    values = c("In-sample" = 16, "RCV" = 1, "SCV" = 17)
  ) +
  scale_y_continuous(
    limits = c(0.4, 1.05),
    breaks = seq(0.4, 1.0, by = 0.1)
  ) +
  labs(
    x = "Mean Prediction Interval Width (MPIW, days)",
    y = "Interval Coverage Probability (ICP)",
    color = "Species",
    linetype = "Validation",
    shape = "Evaluation"
  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 1),
    shape = guide_legend(nrow = 1)
  ) +
  theme_nature

print(p_traj)

ggsave(
  filename = "Figure_MPIW_ICP_Trajectory_RCV_SCV.png",
  plot = p_traj,
  width = 14,
  height = 8,
  dpi = 500
)

ggsave("Figure_MPIW_ICP_Trajectory_RCV_SCV.pdf", plot = p_traj, width = 14, height = 8)


# ======================================================
# Figure: Pinball Loss change with error bars
# Nature-style formatting
# ======================================================

base_font <- "sans"

theme_nature <- theme_bw(base_size = 18, base_family = base_font) +
  theme(
    text = element_text(family = base_font),

    panel.grid.minor = element_blank(),

    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(face = "bold", size = 18),

    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x  = element_text(size = 17),
    axis.text.y  = element_text(size = 17),

    legend.position = "bottom",
    legend.title = element_text(size = 17),
    legend.text  = element_text(size = 15),

    plot.margin = margin(t = 10, r = 10, b = 10, l = 15)
  )

df_pinball <- data.frame(
  Species = rep(
    c("H. chinensis", "Z. moltrechti", "F. limnocharis", "M. fissipes"),
    each = 18
  ),
  Model = rep(
    c("10 CSCIs", "VIP > 1.0 CSCIs"),
    times = 4,
    each = 9
  ),
  Stage = rep(
    c(rep("In-sample", 3), rep("RCV", 3), rep("SCV", 3)),
    times = 8
  ),
  Quantile = rep(c("25th", "50th", "75th"), times = 24),
  Value = c(
    # H. chinensis, 10 CSCIs
    0.838, 0.289, 0.979,
    3.077, 3.118, 2.641,
    5.076, 6.694, 7.409,

    # H. chinensis, VIP > 1.0 CSCIs
    2.097, 2.503, 2.135,
    3.580, 4.627, 3.180,
    5.644, 8.111, 7.917,

    # Z. moltrechti, 10 CSCIs
    1.606, 0.406, 1.820,
    5.929, 6.203, 6.053,
    7.394, 8.337, 6.783,

    # Z. moltrechti, VIP > 1.0 CSCIs
    4.837, 5.298, 5.489,
    9.086, 11.778, 8.639,
    11.121, 11.691, 9.956,

    # F. limnocharis, 10 CSCIs
    1.058, 0.618, 0.939,
    2.626, 3.416, 3.545,
    3.010, 4.185, 3.526,

    # F. limnocharis, VIP > 1.0 CSCIs
    2.701, 3.639, 3.581,
    3.859, 4.920, 5.215,
    4.370, 5.602, 5.182,

    # M. fissipes, 10 CSCIs
    1.253, 0.233, 1.239,
    4.104, 5.273, 4.474,
    5.455, 7.009, 5.904,

    # M. fissipes, VIP > 1.0 CSCIs
    3.682, 3.610, 4.093,
    6.028, 7.590, 8.345,
    7.863, 9.063, 8.735
  ),
  SD = c(
    # H. chinensis, 10 CSCIs
    NA, NA, NA,
    1.202, 0.730, 0.844,
    3.467, 5.486, 7.420,

    # H. chinensis, VIP > 1.0 CSCIs
    NA, NA, NA,
    1.290, 1.538, 0.864,
    3.637, 6.562, 8.172,

    # Z. moltrechti, 10 CSCIs
    NA, NA, NA,
    2.795, 0.887, 1.142,
    4.110, 2.836, 2.793,

    # Z. moltrechti, VIP > 1.0 CSCIs
    NA, NA, NA,
    1.914, 1.601, 1.674,
    3.807, 4.308, 3.357,

    # F. limnocharis, 10 CSCIs
    NA, NA, NA,
    0.722, 1.155, 0.995,
    0.557, 1.024, 0.494,

    # F. limnocharis, VIP > 1.0 CSCIs
    NA, NA, NA,
    1.189, 2.008, 2.004,
    1.055, 0.776, 0.841,

    # M. fissipes, 10 CSCIs
    NA, NA, NA,
    2.215, 2.342, 1.779,
    1.792, 2.592, 2.629,

    # M. fissipes, VIP > 1.0 CSCIs
    NA, NA, NA,
    2.405, 3.129, 2.951,
    1.766, 2.971, 3.561
  )
)

df_pinball <- df_pinball %>%
  mutate(
    Species = factor(
      Species,
      levels = c(
        "H. chinensis",
        "Z. moltrechti",
        "F. limnocharis",
        "M. fissipes"
      )
    ),
    Model = recode(
      Model,
      "10 CSCIs" = "bold(a)",
      "VIP > 1.0 CSCIs" = "bold(b)"
    ),
    Model = factor(
      Model,
      levels = c("bold(a)", "bold(b)")
    ),
    Stage = factor(
      Stage,
      levels = c("In-sample", "RCV", "SCV")
    ),
    Quantile = factor(
      Quantile,
      levels = c("25th", "50th", "75th")
    ),
    Species_label = case_when(
      Species == "H. chinensis" ~ "italic(H.~chinensis)",
      Species == "Z. moltrechti" ~ "italic(Z.~moltrechti)",
      Species == "F. limnocharis" ~ "italic(F.~limnocharis)",
      Species == "M. fissipes" ~ "italic(M.~fissipes)"
    ),
    Quantile_num = as.numeric(Quantile),
    Stage_offset = case_when(
      Stage == "In-sample" ~ -0.26,
      Stage == "RCV" ~ 0,
      Stage == "SCV" ~ 0.26
    ),
    x_pos = Quantile_num + Stage_offset
  )

p_pinball <- ggplot(df_pinball) +
  geom_col(
    aes(x = x_pos, y = Value, fill = Stage),
    width = 0.22,
    color = "black",
    linewidth = 0.25,
    alpha = 0.88
  ) +
  geom_errorbar(
    data = df_pinball %>% filter(!is.na(SD)),
    aes(
      x = x_pos,
      ymin = pmax(Value - SD, 0),
      ymax = Value + SD
    ),
    width = 0.08,
    linewidth = 0.45,
    color = "black"
  ) +
  scale_x_continuous(
    breaks = c(1, 2, 3),
    labels = c("25th", "50th", "75th")
  ) +
  facet_grid(
    Species_label ~ Model,
    scales = "free_y",
    labeller = label_parsed
  ) +
  labs(
    x = "Quantile",
    y = "Pinball Loss",
    fill = "Evaluation"
  ) +
  theme_nature

print(p_pinball)

ggsave(
  filename = "Figure_Pinball_Loss_Change_with_ErrorBars.png",
  plot = p_pinball,
  width = 11,
  height = 11.5,
  dpi = 500
)

ggsave("Figure_Pinball_Loss_Change_with_ErrorBars.pdf", plot = p_pinball, width = 11, height = 11.5)


# ======================================================
# Figure: MAE change with error bars
# Nature-style formatting
# ======================================================

base_font <- "sans"

theme_nature <- theme_bw(base_size = 18, base_family = base_font) +
  theme(
    text = element_text(family = base_font),

    panel.grid.minor = element_blank(),

    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(face = "bold", size = 18),

    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x  = element_text(size = 17),
    axis.text.y  = element_text(size = 17),

    legend.position = "bottom",
    legend.title = element_text(size = 17),
    legend.text  = element_text(size = 15),

    plot.margin = margin(t = 10, r = 10, b = 10, l = 15)
  )

df_mae <- data.frame(
  Species = rep(
    c("H. chinensis", "Z. moltrechti", "F. limnocharis", "M. fissipes"),
    each = 18
  ),
  Model = rep(
    c("10 CSCIs", "VIP > 1.0 CSCIs"),
    times = 4,
    each = 9
  ),
  Stage = rep(
    c(rep("In-sample", 3), rep("RCV", 3), rep("SCV", 3)),
    times = 8
  ),
  Quantile = rep(c("25th", "50th", "75th"), times = 24),
  Value = c(
    # H. chinensis, 10 CSCIs
    3.35, 0.58, 3.92,
    10.24, 6.24, 9.56,
    16.44, 13.39, 17.76,

    # H. chinensis, VIP > 1.0 CSCIs
    8.37, 5.01, 7.96,
    11.91, 9.25, 11.21,
    18.37, 16.22, 18.88,

    # Z. moltrechti, 10 CSCIs
    6.42, 0.81, 7.28,
    19.06, 12.41, 20.37,
    24.58, 16.67, 21.77,

    # Z. moltrechti, VIP > 1.0 CSCIs
    17.80, 10.60, 19.93,
    26.38, 23.56, 28.34,
    33.78, 23.38, 29.85,

    # F. limnocharis, 10 CSCIs
    4.23, 1.24, 3.76,
    8.64, 6.83, 10.41,
    9.51, 8.37, 11.45,

    # F. limnocharis, VIP > 1.0 CSCIs
    9.89, 7.28, 10.66,
    12.25, 9.84, 14.67,
    13.18, 11.20, 14.25,

    # M. fissipes, 10 CSCIs
    5.01, 0.47, 4.96,
    12.96, 10.55, 13.53,
    15.56, 14.02, 18.31,

    # M. fissipes, VIP > 1.0 CSCIs
    13.60, 7.22, 12.32,
    19.07, 15.18, 20.97,
    22.31, 18.13, 23.28
  ),
  SD = c(
    # H. chinensis, 10 CSCIs
    NA, NA, NA,
    3.28, 1.46, 3.40,
    13.87, 10.97, 8.45,

    # H. chinensis, VIP > 1.0 CSCIs
    NA, NA, NA,
    4.51, 3.08, 3.74,
    14.58, 13.12, 9.34,

    # Z. moltrechti, 10 CSCIs
    NA, NA, NA,
    8.39, 1.77, 3.97,
    11.37, 5.67, 8.88,

    # Z. moltrechti, VIP > 1.0 CSCIs
    NA, NA, NA,
    3.28, 3.20, 5.57,
    12.14, 8.62, 9.18,

    # F. limnocharis, 10 CSCIs
    NA, NA, NA,
    2.56, 2.31, 1.44,
    1.80, 2.05, 1.63,

    # F. limnocharis, VIP > 1.0 CSCIs
    NA, NA, NA,
    1.189, 2.008, 2.004,
    1.91, 1.55, 2.10,

    # M. fissipes, 10 CSCIs
    NA, NA, NA,
    5.67, 4.68, 4.51,
    1.84, 5.18, 7.73,

    # M. fissipes, VIP > 1.0 CSCIs
    NA, NA, NA,
    6.10, 6.26, 4.87,
    5.13, 5.94, 5.79
  )
)

df_mae <- df_mae %>%
  mutate(
    Species = factor(
      Species,
      levels = c(
        "F. limnocharis",
        "H. chinensis",
        "M. fissipes",
        "Z. moltrechti"
      )
    ),
    Model = recode(
      Model,
      "10 CSCIs" = "bold(a)",
      "VIP > 1.0 CSCIs" = "bold(b)"
    ),
    Model = factor(
      Model,
      levels = c("bold(a)", "bold(b)")
    ),
    Stage = factor(
      Stage,
      levels = c("In-sample", "RCV", "SCV")
    ),
    Quantile = factor(
      Quantile,
      levels = c("25th", "50th", "75th")
    ),
    Species_label = case_when(
      Species == "F. limnocharis" ~ "italic(F.~limnocharis)",
      Species == "H. chinensis" ~ "italic(H.~chinensis)",
      Species == "M. fissipes" ~ "italic(M.~fissipes)",
      Species == "Z. moltrechti" ~ "italic(Z.~moltrechti)"
    ),
    Quantile_num = as.numeric(Quantile),
    Stage_offset = case_when(
      Stage == "In-sample" ~ -0.26,
      Stage == "RCV" ~ 0,
      Stage == "SCV" ~ 0.26
    ),
    x_pos = Quantile_num + Stage_offset
  )

p_mae <- ggplot(df_mae) +
  geom_col(
    aes(x = x_pos, y = Value, fill = Stage),
    width = 0.22,
    color = "black",
    linewidth = 0.25,
    alpha = 0.88
  ) +
  geom_errorbar(
    data = df_mae %>% filter(!is.na(SD)),
    aes(
      x = x_pos,
      ymin = pmax(Value - SD, 0),
      ymax = Value + SD
    ),
    width = 0.08,
    linewidth = 0.45,
    color = "black"
  ) +
  scale_x_continuous(
    breaks = c(1, 2, 3),
    labels = c("25th", "50th", "75th")
  ) +
  facet_grid(
    Species_label ~ Model,
    scales = "free_y",
    labeller = label_parsed
  ) +
  labs(
    x = "Quantile",
    y = "Mean Absolute Error (days)",
    fill = "Evaluation"
  ) +
  theme_nature

print(p_mae)

ggsave(
  filename = "Figure_MAE_Change_with_ErrorBars.png",
  plot = p_mae,
  width = 11,
  height = 11.5,
  dpi = 500
)

ggsave("Figure_MAE_Change_with_ErrorBars.pdf", plot = p_mae, width = 11, height = 11.5)

# ======================================================
# Figure: DOY distribution of four species
# Bin width = 5 days
# Output: PDF + 500 dpi PNG
# ======================================================


# Install and load 'patchwork' if not already installed

# ---------- Nature-style theme ----------
base_font <- "sans"

theme_nature <- theme_bw(base_size = 18, base_family = base_font) +
  theme(
    text = element_text(family = base_font),

    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),

    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),

    plot.title = element_text(
      size = 18,
      face = "italic",
      hjust = 0
    ),

    plot.tag = element_text(
      size = 24,
      face = "bold"
    ),

    panel.grid.minor = element_blank()
  )

# ======================================================
# Read four species datasets
# ======================================================

species_files <- c("Hc.csv", "Zm.csv", "Fl.csv", "Mf.csv")
missing_species_files <- species_files[!file.exists(species_files)]
if (length(missing_species_files)) {
  message("Skipping DOY distribution figure. Missing files: ",
          paste(missing_species_files, collapse = ", "))
} else {

species_info <- data.frame(
  file = c(
    "Hc.csv",
    "Zm.csv",
    "Fl.csv",
    "Mf.csv"
  ),
  species = c(
    "H. chinensis",
    "Z. moltrechti",
    "F. limnocharis",
    "M. fissipes"
  )
)

# ======================================================
# Read data and calculate DOY
# ======================================================

read_species_data <- function(file, species_name){

  df <- read.csv(
    file,
    stringsAsFactors = FALSE
  )

  df %>%
    mutate(
      time_observed_at = as.Date(time_observed_at),
      DOY = as.numeric(format(time_observed_at, "%j")),
      Species = species_name
    ) %>%
    filter(!is.na(DOY))
}

df_all <- bind_rows(
  mapply(
    read_species_data,
    species_info$file,
    species_info$species,
    SIMPLIFY = FALSE
  )
)

# ======================================================
# Single panel function
# ======================================================

make_doy_plot <- function(data, species_name){

  ggplot(
    data %>% filter(Species == species_name),
    aes(x = DOY)
  ) +
    geom_histogram(
      binwidth = 5,
      boundary = 0,
      fill = "grey55",
      color = "white",
      linewidth = 0.3
    ) +
    scale_x_continuous(
      limits = c(0,365),
      breaks = seq(0,365,50)
    ) +
    labs(
      title = species_name,
      x = "Day of year (DOY)",
      y = "Frequency"
    ) +
    theme_nature
}

# ======================================================
# Create panels
# ======================================================

p_a <- make_doy_plot(df_all, "H. chinensis")
p_b <- make_doy_plot(df_all, "Z. moltrechti")
p_c <- make_doy_plot(df_all, "F. limnocharis")
p_d <- make_doy_plot(df_all, "M. fissipes")

# ======================================================
# Combine vertically
# ======================================================

p_final <- (
  p_a /
  p_b /
  p_c /
  p_d
) +
  plot_annotation(
    tag_levels = "a"
  ) &
  theme(
    plot.tag = element_text(
      family = base_font,
      face = "bold",
      size = 24
    )
  )

print(p_final)

# ======================================================
# Export
# ======================================================

ggsave(
  "Figure_DOY_Distribution_4Species.png",
  plot = p_final,
  width = 8,
  height = 14,
  dpi = 500
)

ggsave(
  "Figure_DOY_Distribution_4Species.pdf",
  plot = p_final,
  width = 8,
  height = 14
)

cat("Figure exported successfully\n")
cat("1. Figure_DOY_Distribution_4Species.png (500 dpi)\n")

}
