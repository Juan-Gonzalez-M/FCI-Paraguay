# ============================================================================
# 66_Submission_Figure5_Predetermined.R — Round-15 colleague item: the
# state-dependence exhibits must use the PREDETERMINED terms-of-trade state.
#
# PROBLEM.  Main-text Table 9 and Figure 5 display the interaction with
# CONTEMPORANEOUS ToT growth, while the text itself calls the predetermined
# state "the cleaner specification".  A contemporaneous state can reflect the
# same shock or information that moves the FCI, so the exhibits headline the
# least defensible timing while the prose endorses the other one.
#
# FIX.  Re-estimate the interaction LP with the predetermined state (12-month
# trailing MA of d_ToT, lagged one month -- known before t; identical to
# 58_Small_Robustness_Items.R:83) and rebuild the marginal-effect figure on it.
# The estimator, controls and NW(h+1) inference replicate script 58 Part 2;
# the marginal-effect grid and delta-method standard errors replicate
# 21_FCI_Commodity_Puzzle.R:631-666.
#
# NOTE ON THE FIGURE'S SECOND AXIS.  sd(ToT_MA12) is NOT sd(d_ToT), so the
# "1 SD = x pp" conversion in the Figure 5 caption changes.  Both SDs are
# printed below; the caption must be updated from this script's output.
#
# Outputs
#   output/submission/Figure_5.png                    (predetermined; main text)
#   output/submission/Figure_A_ToT_Contemporaneous.png (appendix sensitivity)
#   output/revision/csv/Rev_Marginal_Effects_Predetermined.csv
#   output/revision/csv/Rev_ToT_Interaction_BothStates.csv
#
# ADDITIVE apart from the deliberate overwrite of output/submission/Figure_5.png
# (the previous file is snapshotted in output/archive/pre_round15/).
# ============================================================================

t0 <- Sys.time()
source("revision_helpers.R")
suppressMessages({library(zoo); library(ggplot2)})
cat("=== 66: predetermined-ToT state dependence (Table 9 / Figure 5) ===\n")

SUB_DIR <- file.path(micro_paths$root, "output", "submission")
Z90     <- qnorm(0.95)   # 90% delta-method bands, matching the existing figure

d <- load_ext_data()
setorder(d, fecha)
d[, ToT_MA12 := shift(frollmean(d_ToT, 12, align = "right"), 1)]

XVAR <- "FCI_exCredit_AVG"   # main-text interaction index

# ---------------------------------------------------------------------------
# Interaction LP returning the full HAC covariance block needed for the
# delta-method marginal effects (script 58 keeps only the interaction row).
# ---------------------------------------------------------------------------
int_lp_cov <- function(data, xvar, state, h) {
  dh <- as.data.frame(data)
  dh$y_fwd  <- dplyr::lead(dh$Cred_Real_Total, h)
  dh$y_lag1 <- dplyr::lag(dh$Cred_Real_Total, 1)
  dh$y_lag2 <- dplyr::lag(dh$Cred_Real_Total, 2)
  dh$x_lag1 <- dplyr::lag(dh[[xvar]], 1)
  dh$FCIxS  <- dh[[xvar]] * dh[[state]]
  rhs <- c(xvar, state, "FCIxS", "y_lag1", "y_lag2", "x_lag1", EXOG_MACRO)
  rd  <- na.omit(dh[, c("y_fwd", rhs)])
  if (nrow(rd) < 50) return(NULL)
  m  <- lm(as.formula(paste("y_fwd ~", paste(rhs, collapse = "+"))), data = rd)
  V  <- sandwich::NeweyWest(m, lag = h + 1, prewhite = FALSE)
  ct <- lmtest::coeftest(m, vcov = V)
  list(
    tab = data.table(
      state = state, horizon = h, n = nrow(rd),
      b_fci    = ct[xvar, 1],   se_fci    = ct[xvar, 2],   p_fci    = ct[xvar, 4],
      b_state  = ct[state, 1],  se_state  = ct[state, 2],  p_state  = ct[state, 4],
      b_inter  = ct["FCIxS", 1], se_inter = ct["FCIxS", 2], p_inter = ct["FCIxS", 4]),
    beta1 = unname(coef(m)[xvar]),
    beta3 = unname(coef(m)["FCIxS"]),
    v1    = V[xvar, xvar],
    v3    = V["FCIxS", "FCIxS"],
    c13   = V[xvar, "FCIxS"])
}

# ---------------------------------------------------------------------------
# Marginal effects over a +/- 2 SD grid of the state (script 21 convention)
# ---------------------------------------------------------------------------
marginal_grid <- function(fits, state_sd, label) {
  out <- data.table()
  grid <- seq(-2 * state_sd, 2 * state_sd, length.out = 41)
  for (hh in names(fits)) {
    f <- fits[[hh]]
    if (is.null(f)) next
    me    <- f$beta1 + f$beta3 * grid
    se_me <- sqrt(f$v1 + grid^2 * f$v3 + 2 * grid * f$c13)
    out <- rbind(out, data.table(
      state = label, horizon = as.integer(hh),
      state_value = grid, state_sd_units = grid / state_sd,
      marginal_effect = me, se = se_me,
      ci_lower = me - Z90 * se_me, ci_upper = me + Z90 * se_me))
  }
  out
}

make_plot <- function(me, state_sd, x_lab, sec_lab) {
  me <- copy(me)
  me[, horizon_label := factor(paste0("h = ", horizon),
                               levels = c("h = 6", "h = 12", "h = 18"))]
  ggplot(me, aes(x = state_sd_units, y = marginal_effect)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "steelblue") +
    geom_line(color = "steelblue", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray70") +
    facet_wrap(~horizon_label, ncol = 3) +
    # Secondary ticks are pinned to the primary SD ticks; with free-standing
    # breaks the facet edges collide (e.g. "30-30") at this aspect ratio.
    scale_x_continuous(name = x_lab, breaks = -2:2,
                       sec.axis = sec_axis(~ . * state_sd, name = sec_lab,
                                           breaks = round((-2:2) * state_sd))) +
    labs(y = "Marginal FCI effect on real credit growth (pp per index unit)") +
    theme_minimal(base_size = 11) +
    theme(strip.text = element_text(face = "bold"))
}

# ---------------------------------------------------------------------------
# Estimate both states
# ---------------------------------------------------------------------------
HS <- c(6, 12, 18)
fits_pre <- setNames(lapply(HS, function(h) int_lp_cov(d, XVAR, "ToT_MA12", h)), HS)
fits_con <- setNames(lapply(HS, function(h) int_lp_cov(d, XVAR, "d_ToT",    h)), HS)

tab <- rbind(rbindlist(lapply(fits_pre, `[[`, "tab")),
             rbindlist(lapply(fits_con, `[[`, "tab")))
tab[, state_label := fifelse(state == "ToT_MA12",
                             "predetermined (MA12 d_ToT, lagged)",
                             "contemporaneous d_ToT")]
write_rev_csv(tab, "Rev_ToT_Interaction_BothStates.csv")

cat("\n--- Interaction coefficients, both states ---\n")
print(tab[, .(state_label, horizon, b_inter = round(b_inter, 4),
              se_inter = round(se_inter, 4), p_inter = round(p_inter, 4), n)])

# SDs used for the grids and for the figure's secondary axis
sd_pre <- sd(d$ToT_MA12, na.rm = TRUE)
sd_con <- sd(d$d_ToT,    na.rm = TRUE)
cat(sprintf("\nState SDs:  predetermined MA12 = %.4f pp   contemporaneous = %.4f pp\n",
            sd_pre, sd_con))
cat(sprintf("CAPTION INPUT: 1 SD of the predetermined state = %.1f pp of ToT growth\n",
            sd_pre))

me_pre <- marginal_grid(fits_pre, sd_pre, "predetermined")
me_con <- marginal_grid(fits_con, sd_con, "contemporaneous")
write_rev_csv(rbind(me_pre, me_con), "Rev_Marginal_Effects_Predetermined.csv")

cat("\nMarginal effects at key state values (h = 12, predetermined):\n")
m12 <- me_pre[horizon == 12]
for (sv in c(-2, -1, 0, 1, 2)) {
  r <- m12[which.min(abs(state_sd_units - sv))]
  cat(sprintf("  state = %+d SD: ME = %+.2f pp  [%.2f, %.2f]\n",
              sv, r$marginal_effect, r$ci_lower, r$ci_upper))
}

# ---------------------------------------------------------------------------
# Figures. No embedded title/subtitle: the caption carries them in the
# manuscript (round-14 convention, cf. 62_Submission_Figure1.R).
# ---------------------------------------------------------------------------
p_pre <- make_plot(me_pre, sd_pre,
                   "Predetermined terms-of-trade state (SD units)",
                   "12-month trailing ToT growth, lagged (pp)")
ggsave(file.path(SUB_DIR, "Figure_5.png"), p_pre, width = 14, height = 5, dpi = 300)
cat("  [png]  output/submission/Figure_5.png (predetermined state)\n")

p_con <- make_plot(me_con, sd_con,
                   "Contemporaneous terms-of-trade growth (SD units)",
                   "ToT growth (pp)")
ggsave(file.path(SUB_DIR, "Figure_A_ToT_Contemporaneous.png"), p_con,
       width = 14, height = 5, dpi = 300)
cat("  [png]  output/submission/Figure_A_ToT_Contemporaneous.png (appendix sensitivity)\n")

cat(sprintf("\nDone in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
