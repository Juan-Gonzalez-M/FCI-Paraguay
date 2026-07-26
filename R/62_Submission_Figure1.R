# ============================================================================
# 62_Submission_Figure1.R — Round-14 colleague item B2.
#
# PROBLEM.  Submission Figure 1 (\label{fig:fci}) pointed at Figure_2.png,
# which is the FOUR-METHOD COMPARISON chart (legend DFM/PCA/VAR/ZSCORE,
# internal title "Comprehensive FCI: Method Comparison").  The caption and
# Section 3's five-episode paragraph describe a different exhibit: the
# composite FCI annotated with the five stress episodes.  That figure did not
# exist in the pipeline.  This script produces it.
#
# Also per the same item: no embedded chart title or subtitle — journals
# typeset the caption, and an internal title duplicates it.
#
# Source: output/csv/FCI_Complete_Results.csv (FCI_COMP_AVG, the four-method
# composite described in Section 3).  No estimation; plotting only.
# Output: ../output/submission/Figure_1.png (manuscript Figure 1)
# ============================================================================

suppressMessages({ library(ggplot2); library(dplyr) })

d <- read.csv("../output/csv/FCI_Complete_Results.csv")
d$fecha <- as.Date(d$fecha)
d <- d %>% dplyr::select(fecha, FCI = FCI_COMP_AVG) %>% filter(!is.na(FCI))
stopifnot(nrow(d) == 360)

# Five stress episodes named in Section 3 (identified without ex post tuning)
ep <- data.frame(
  label = c("Domestic banking\ncrisis", "Regional\ncontagion",
            "Global financial\ncrisis", "Commodity bust +\ndepreciation",
            "COVID-19"),
  start = as.Date(c("1998-01-01", "2002-01-01", "2008-09-01",
                    "2015-01-01", "2020-03-01")),
  end   = as.Date(c("1999-12-31", "2003-12-31", "2009-12-31",
                    "2016-12-31", "2021-12-31")))
ep$mid <- ep$start + (ep$end - ep$start) / 2

ymax <- max(d$FCI, na.rm = TRUE)
ymin <- min(d$FCI, na.rm = TRUE)
ep$y  <- ymax + 0.42

p <- ggplot(d, aes(fecha, FCI)) +
  geom_rect(data = ep, inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "grey70", alpha = 0.28) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 0.55, colour = "#1F3B57") +
  geom_text(data = ep, inherit.aes = FALSE,
            aes(x = mid, y = y, label = label),
            size = 2.9, lineheight = 0.92, colour = "grey25") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y",
               expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(breaks = seq(-3, 3, 1)) +
  coord_cartesian(ylim = c(ymin - 0.15, ymax + 0.85)) +
  labs(x = NULL, y = "FCI (index units)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_blank(), plot.subtitle = element_blank())

ggsave("../output/submission/Figure_1.png", p, width = 10, height = 4.6,
       dpi = 200, bg = "white")
cat(sprintf("Figure_1.png written: composite FCI %s to %s, %d obs, 5 episodes.\n",
            min(d$fecha), max(d$fecha), nrow(d)))
