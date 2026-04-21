# Analyzer for disturbance simulations

library(tidyverse)
library(patchwork)
library(plotrix)

data <- read.csv(file.path("Data", "3-dist", "dist_results_s1.csv"))

ggplot(
    data = data,
    mapping = aes(x = separation_distance, y = Relative_concentration)
) +
    geom_hline(yintercept = 0.5, linetype = "longdash", color="darkgrey") +

    stat_summary(fun.data = "mean_se", geom="errorbar", width=0.35, mapping = aes(color=factor(disturb_freq))) +

    geom_point(stat = "summary", fun="mean", width=0.0005, mapping = aes(color=factor(disturb_freq), shape=factor(disturb_freq))) +
    geom_path(stat = "summary", fun="mean", mapping = aes(color=factor(disturb_freq))) +

    labs(
        x = "D",
        y = "Rc",
        color = quote(lambda), shape = quote(lambda)
    ) +

    coord_cartesian(xlim=c(1,9), ylim=c(0, 1)) +
    scale_x_continuous(breaks=c(1,3,5,7,9)) +
    theme_bw()
    ggsave(file.path("Analysis", "Figures", "dist_graph_s1.svg"), width = 6, height = 4)
