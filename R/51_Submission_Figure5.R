# ============================================================================
# 51_Submission_Figure5.R — Regenerates submission Figure 5 (bank-level hedging
# channel, two panels) from stored pipeline outputs, with the exposure label
# harmonized to "baseline sector FX-credit exposure" (revision round 7).
# Panel (a): hedging gradient, full sample vs COVID-excluded
#            (Rev_DesignC_Gradient_Robustness.csv: gradient_baseline / _exCOVID)
# Panel (b): split-sample USD-vs-PYG differential, hedged vs unhedged
#            (Rev_DesignC_Split_Robustness.csv: split_baseline)
# 90% bands use bank-clustered SEs (one sensitivity convention; see Table 8).
# Output: ../output/submission/Figure_9.png (manuscript Figure 5)
# ============================================================================

suppressMessages({ library(ggplot2); library(gridExtra) })

grad  <- read.csv("../output/revision/csv/Rev_DesignC_Gradient_Robustness.csv")
split <- read.csv("../output/revision/csv/Rev_DesignC_Split_Robustness.csv")

z90 <- qnorm(0.95)

ga <- subset(grad, variant %in% c("gradient_baseline", "gradient_exCOVID") &
                    param == "shk_usd_fxsh")
ga$series <- ifelse(ga$variant == "gradient_baseline",
                    "Gradient (full sample)", "Gradient (COVID excluded)")
ga$lo <- ga$b - z90 * ga$se_bank
ga$hi <- ga$b + z90 * ga$se_bank

pa <- ggplot(ga, aes(h, b, colour = series, fill = series)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
  scale_colour_manual(values = c("Gradient (full sample)" = "#2C5F8A",
                                 "Gradient (COVID excluded)" = "#C0392B")) +
  scale_fill_manual(values = c("Gradient (full sample)" = "#2C5F8A",
                               "Gradient (COVID excluded)" = "#C0392B")) +
  labs(title = "(a) Hedging gradient",
       subtitle = "Shock x USD x sector FX share (per 1 SD of baseline sector FX-credit exposure)",
       x = "Horizon (months)", y = "pp", colour = NULL, fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey40", size = 10))

sb <- subset(split, variant == "split_baseline")
sb$series <- ifelse(sb$group == "hedged",
                    "Hedged (dollar-earning) sectors", "Unhedged sectors")
sb$lo <- sb$b - z90 * sb$se_bank
sb$hi <- sb$b + z90 * sb$se_bank

pb <- ggplot(sb, aes(h, b, colour = series, fill = series)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
  scale_colour_manual(values = c("Hedged (dollar-earning) sectors" = "#1E8449",
                                 "Unhedged sectors" = "#5D6D7E")) +
  scale_fill_manual(values = c("Hedged (dollar-earning) sectors" = "#1E8449",
                               "Unhedged sectors" = "#5D6D7E")) +
  labs(title = "(b) Split sample: USD-vs-PYG differential by borrower type",
       subtitle = "Shock x USD estimated separately on hedged and unhedged sectors",
       x = "Horizon (months)", y = "pp", colour = NULL, fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey40", size = 10))

g <- gridExtra::arrangeGrob(pa, pb, ncol = 2)
ggsave("../output/submission/Figure_9.png", g, width = 15, height = 6.2,
       dpi = 150, bg = "white")
cat("Figure_9.png regenerated from stored CSVs.\n")
