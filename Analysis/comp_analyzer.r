# Analyzer for competition simulations

library(tidyverse)
library(patchwork)
library(plotrix)

data <- read.csv(file.path("Data", "2-comp", "comp_results.csv"))


ggplot(
    data = data,
    mapping = aes(x = separation_distance, y = Relative_concentration)
) +
    geom_hline(yintercept = 0.5, linetype = "longdash", color="darkgrey") +
    geom_point(stat = "summary", fun="mean", width=0.0005) +
    geom_path(stat = "summary", fun="mean") +
    stat_summary(fun.data = "mean_se", geom="errorbar", width=0.00025) +
    labs(
        x = "Propagule Formation Rate",
        y = "Rc",
    ) +
    ylim(0, 1) + theme_bw()
    ggsave(file.path("Analysis", "Figures", "comp_graph.svg"), width = 8, height = 4)
