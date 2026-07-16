# ============================================================================
# 53_Currency_Migration_Test.R
# Direct test of the currency-"migration" claim (round-8, colleague item S1).
#
# The manuscript inferred a shift in the sensitive credit margin from local-
# to foreign-currency credit from significance patterns and cross-sample
# z-tests (MN p=0.094, USD p=0.984). Both colleagues: comparing significance
# is not a direct test. This script estimates the FCI x post-IT x USD triple
# interaction in one stacked local projection:
#
#   y_{c,t+h} = cur-specific { a + b*FCI_t + g*D_t + d*(FCI_t x D_t)
#               + lags + macro controls } ,  c in {MN, USD}, D = post-IT
#
# with all coefficients currency-specific (equivalent to equation-by-equation
# OLS) and Driscoll-Kraay standard errors (sandwich::vcovPL, cluster = month,
# lag = h+1) so that the cross-currency covariance at the same month and the
# LP-induced serial correlation are both accounted for. The migration test is
#   H0: d_USD - d_MN = 0.
#
# Index: full-sample composite FCI_exCredit_AVG (the "consistent index" of
# the Appendix E.2/E.3 memo rows). Credit series constructed exactly as in
# script 15 (real MN via IPC deflation; USD-equivalent book).
# Output: output/revision/csv/Rev_Currency_Migration_Test.csv
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(sandwich); library(lmtest); library(readxl)
})

POSTIT <- as.Date("2011-05-01")

d <- read.csv("../output/csv/New_External_Variables.csv")
d$fecha <- as.Date(d$fecha)

raw <- read_excel("../data/FCI_data_1.xlsx", sheet = "Main_variables")
mac <- read_excel("../data/FCI_data_1.xlsx", sheet = "Datos_macro")
fecha_col  <- names(raw)[1]; fecha_colm <- names(mac)[1]
cred <- raw %>% rename(fecha = !!sym(fecha_col)) %>%
  mutate(fecha = as.Date(fecha)) %>% arrange(fecha) %>%
  left_join(mac %>% rename(fecha = !!sym(fecha_colm)) %>%
              mutate(fecha = as.Date(fecha)) %>% dplyr::select(fecha, IPC),
            by = "fecha") %>%
  mutate(
    Cred_USD = (Creditos_Sector_privado_USD_equivalente /
                  lag(Creditos_Sector_privado_USD_equivalente, 12) - 1) * 100,
    Cred_Real_MN = ((Creditos_Sector_privado_MN / IPC) /
                      lag(Creditos_Sector_privado_MN / IPC, 12) - 1) * 100) %>%
  dplyr::select(fecha, Cred_USD, Cred_Real_MN)

d <- d %>% left_join(cred, by = "fecha") %>% arrange(fecha) %>%
  mutate(D = as.numeric(fecha >= POSTIT),
         IMAEP_yoy_L1 = lag(IMAEP_yoy, 1), IPC_yoy_L1 = lag(IPC_yoy, 1))

cat("=== 53: FCI x post-IT x USD migration test ===\n\n")

run_h <- function(h) {
  base <- d %>% mutate(
    yUSD_fwd = lead(Cred_USD, h),      yMN_fwd = lead(Cred_Real_MN, h),
    yUSD_l1  = lag(Cred_USD, 1),       yMN_l1  = lag(Cred_Real_MN, 1),
    yUSD_l2  = lag(Cred_USD, 2),       yMN_l2  = lag(Cred_Real_MN, 2),
    x_lag1   = lag(FCI_exCredit_AVG, 1))
  long <- bind_rows(
    base %>% transmute(fecha, cur = "USD", y_fwd = yUSD_fwd, y_l1 = yUSD_l1,
                       y_l2 = yUSD_l2, FCI = FCI_exCredit_AVG, D, x_lag1,
                       IMAEP_yoy, IMAEP_yoy_L1, IPC_yoy, IPC_yoy_L1),
    base %>% transmute(fecha, cur = "MN", y_fwd = yMN_fwd, y_l1 = yMN_l1,
                       y_l2 = yMN_l2, FCI = FCI_exCredit_AVG, D, x_lag1,
                       IMAEP_yoy, IMAEP_yoy_L1, IPC_yoy, IPC_yoy_L1)) %>%
    na.omit()
  long$cur <- factor(long$cur, levels = c("MN", "USD"))
  m <- lm(y_fwd ~ 0 + cur + cur:(FCI * D + y_l1 + y_l2 + x_lag1 +
            IMAEP_yoy + IMAEP_yoy_L1 + IPC_yoy + IPC_yoy_L1), data = long)
  V <- sandwich::vcovPL(m, cluster = long$cur, order.by = long$fecha,
                        lag = h + 1, adjust = TRUE)
  cn <- names(coef(m))
  pick <- function(pat) cn[grepl(pat, cn, fixed = TRUE)]
  b  <- coef(m)
  ct <- lmtest::coeftest(m, vcov = V)
  # per-currency FCI effects pre and post, and the interaction shifts
  k_fci_mn  <- "curMN:FCI";  k_fci_usd <- "curUSD:FCI"
  k_int_mn  <- pick("curMN:FCI:D");  k_int_usd <- pick("curUSD:FCI:D")
  stopifnot(length(k_int_mn) == 1, length(k_int_usd) == 1)
  # triple difference d_USD - d_MN
  L <- rep(0, length(b)); names(L) <- cn
  L[k_int_usd] <- 1; L[k_int_mn] <- -1
  est <- sum(L * b); se <- sqrt(as.numeric(t(L) %*% V %*% L))
  z <- est / se; p <- 2 * pnorm(-abs(z))
  data.frame(
    horizon = h, n_stacked = nrow(long), n_months = length(unique(long$fecha)),
    b_MN_pre = b[k_fci_mn], b_MN_shift = b[k_int_mn],
    p_MN_shift = ct[k_int_mn, 4],
    b_USD_pre = b[k_fci_usd], b_USD_shift = b[k_int_usd],
    p_USD_shift = ct[k_int_usd, 4],
    tripleDiff = est, se = se, z = z, p = p)
}

out <- do.call(rbind, lapply(c(6, 12, 18), run_h))
print(out %>% mutate(across(where(is.numeric), ~round(., 4))), row.names = FALSE)
write.csv(out, "../output/revision/csv/Rev_Currency_Migration_Test.csv", row.names = FALSE)
cat("\nSaved: Rev_Currency_Migration_Test.csv\n")
h12 <- out[out$horizon == 12, ]
cat(sprintf("\nMigration verdict (h=12): triple diff = %.2f (p=%.3f) -> %s\n",
    h12$tripleDiff, h12$p,
    ifelse(h12$p < 0.10, "direct evidence of a currency-composition shift",
           "NOT directly established; use 'point estimates suggest a relative shift' wording")))

# ---------------------------------------------------------------------------
# Part 2 (round 9): validation via the USD-minus-MN DIFFERENCE regression.
# Colleague concern: Driscoll-Kraay with only two currency units invites
# questions. The clean equivalent is a single time-series LP on the monthly
# difference g^USD - g^MN:
#   (g^USD - g^MN)_{t+h} = a + b*FCI_t + g*D_t + delta*(FCI_t x D_t)
#                          + lags of the difference + FCI lag + macro controls
# with Newey-West(h+1) and a moving-block bootstrap over months (block 18,
# recentered percentile-t). delta targets the same triple-difference contrast.
# Part 3: the four implied slopes (MN pre/post, USD pre/post) with 95% CIs
# from the stacked model's DK covariance, for the redesigned Table E.3b.
# ---------------------------------------------------------------------------
cat("\n--- Part 2: USD-minus-MN difference regression ---\n")
set.seed(20260714)
dd <- d %>% mutate(gdiff = Cred_USD - Cred_Real_MN)

run_diff <- function(h, B = 4999, block = 18) {
  dh <- dd %>% mutate(
    y_fwd = lead(gdiff, h),
    y_l1 = lag(gdiff, 1), y_l2 = lag(gdiff, 2),
    x_lag1 = lag(FCI_exCredit_AVG, 1),
    FCIxD = FCI_exCredit_AVG * D)
  rhs <- c("FCI_exCredit_AVG", "D", "FCIxD", "y_l1", "y_l2", "x_lag1",
           "IMAEP_yoy", "IMAEP_yoy_L1", "IPC_yoy", "IPC_yoy_L1")
  rd <- dh %>% dplyr::select(y_fwd, all_of(rhs)) %>% na.omit()
  m <- lm(as.formula(paste("y_fwd ~", paste(rhs, collapse = "+"))), data = rd)
  ct <- lmtest::coeftest(m, vcov = sandwich::NeweyWest(m, lag = h + 1, prewhite = FALSE))
  b_obs <- ct["FCIxD", 1]; t_obs <- ct["FCIxD", 3]
  # moving-block bootstrap over months (recentered percentile-t)
  n <- nrow(rd); nb <- ceiling(n / block)
  tb <- rep(NA_real_, B)
  for (bb in seq_len(B)) {
    starts <- sample(seq_len(n - block + 1), nb, replace = TRUE)
    idx <- unlist(lapply(starts, function(s) s:(s + block - 1)))[1:n]
    rb <- rd[idx, ]
    mb <- tryCatch(lm(as.formula(paste("y_fwd ~", paste(rhs, collapse = "+"))), data = rb),
                   error = function(e) NULL)
    if (is.null(mb) || !"FCIxD" %in% names(coef(mb))) next
    cb <- tryCatch(lmtest::coeftest(mb, vcov = sandwich::NeweyWest(mb, lag = h + 1,
                                                                   prewhite = FALSE)),
                   error = function(e) NULL)
    if (!is.null(cb)) tb[bb] <- (cb["FCIxD", 1] - b_obs) / cb["FCIxD", 2]
  }
  p_mb <- mean(abs(tb) >= abs(t_obs), na.rm = TRUE)
  data.frame(horizon = h, n = n, delta = b_obs, se_nw = ct["FCIxD", 2],
             p_nw = ct["FCIxD", 4], p_movingblock = p_mb,
             b_FCI_pre = ct["FCI_exCredit_AVG", 1],
             p_FCI_pre = ct["FCI_exCredit_AVG", 4])
}
diffres <- do.call(rbind, lapply(c(6, 12, 18), run_diff))
print(diffres %>% mutate(across(where(is.numeric), ~round(., 4))), row.names = FALSE)
write.csv(diffres, "../output/revision/csv/Rev_Currency_Migration_Diff.csv", row.names = FALSE)

cat("\n--- Part 3: four implied slopes from the stacked model (DK 95% CIs) ---\n")
slopes_h <- function(h) {
  base <- d %>% mutate(
    yUSD_fwd = lead(Cred_USD, h),      yMN_fwd = lead(Cred_Real_MN, h),
    yUSD_l1  = lag(Cred_USD, 1),       yMN_l1  = lag(Cred_Real_MN, 1),
    yUSD_l2  = lag(Cred_USD, 2),       yMN_l2  = lag(Cred_Real_MN, 2),
    x_lag1   = lag(FCI_exCredit_AVG, 1))
  long <- bind_rows(
    base %>% transmute(fecha, cur = "USD", y_fwd = yUSD_fwd, y_l1 = yUSD_l1,
                       y_l2 = yUSD_l2, FCI = FCI_exCredit_AVG, D, x_lag1,
                       IMAEP_yoy, IMAEP_yoy_L1, IPC_yoy, IPC_yoy_L1),
    base %>% transmute(fecha, cur = "MN", y_fwd = yMN_fwd, y_l1 = yMN_l1,
                       y_l2 = yMN_l2, FCI = FCI_exCredit_AVG, D, x_lag1,
                       IMAEP_yoy, IMAEP_yoy_L1, IPC_yoy, IPC_yoy_L1)) %>% na.omit()
  long$cur <- factor(long$cur, levels = c("MN", "USD"))
  m <- lm(y_fwd ~ 0 + cur + cur:(FCI * D + y_l1 + y_l2 + x_lag1 +
            IMAEP_yoy + IMAEP_yoy_L1 + IPC_yoy + IPC_yoy_L1), data = long)
  V <- sandwich::vcovPL(m, cluster = long$cur, order.by = long$fecha,
                        lag = h + 1, adjust = TRUE)
  b <- coef(m); cn <- names(b)
  lin <- function(w) {
    L <- rep(0, length(b)); names(L) <- cn; L[names(w)] <- w
    est <- sum(L * b); se <- sqrt(as.numeric(t(L) %*% V %*% L))
    c(est = est, lo95 = est - 1.96 * se, hi95 = est + 1.96 * se)
  }
  kmn <- "curMN:FCI"; kusd <- "curUSD:FCI"
  kmnD <- grep("curMN:FCI:D", cn, value = TRUE); kusdD <- grep("curUSD:FCI:D", cn, value = TRUE)
  rows <- rbind(
    MN_pre   = lin(setNames(1, kmn)),
    MN_post  = lin(setNames(c(1, 1), c(kmn, kmnD))),
    USD_pre  = lin(setNames(1, kusd)),
    USD_post = lin(setNames(c(1, 1), c(kusd, kusdD))),
    MN_shift  = lin(setNames(1, kmnD)),
    USD_shift = lin(setNames(1, kusdD)),
    triple_diff = lin(setNames(c(1, -1), c(kusdD, kmnD))))
  data.frame(horizon = h, quantity = rownames(rows), rows, row.names = NULL)
}
slopes <- do.call(rbind, lapply(c(6, 12, 18), slopes_h))
print(slopes %>% mutate(across(where(is.numeric), ~round(., 3))), row.names = FALSE)
write.csv(slopes, "../output/revision/csv/Rev_Currency_Migration_Slopes.csv", row.names = FALSE)
cat("\nSaved: Rev_Currency_Migration_Diff.csv, Rev_Currency_Migration_Slopes.csv\n")
v12 <- diffres[diffres$horizon == 12, ]
cat(sprintf("\nVALIDATION (h=12): difference-regression delta = %.2f (NW p=%.4f, moving-block p=%.4f) vs stacked -14.60 (p=0.001) -> %s\n",
    v12$delta, v12$p_nw, v12$p_movingblock,
    ifelse(v12$p_nw < 0.10 & sign(v12$delta) < 0, "AGREES", "DISAGREES - weaken manuscript language")))
