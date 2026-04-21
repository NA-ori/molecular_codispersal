# Analyzer for competition simulations

library(tidyverse)
library(patchwork)
library(plotrix)

data <- read.csv(file.path("Data", "2-comp", "results_s1.csv"))

# s1 results

ggplot(
    data = data,
    mapping = aes(x = separation_distance, y = Relative_concentration)
) +
    geom_hline(yintercept = 0.5, linetype = "longdash", color="darkgrey") +
    geom_point(stat = "summary", fun="mean", color="#74036a", width=0.00025) +
    geom_path(stat = "summary", fun="mean", color="#74036a") +
    stat_summary(fun.data = "mean_se", geom="errorbar", width=0.25) +
    labs(
        x = "D",
        y = "Rc",
    ) +
    coord_cartesian(xlim=c(1,9), ylim=c(0, 1)) + theme_bw()
    ggsave(file.path("Analysis", "Figures", "comp_graph_s1.svg"), width = 6, height = 4)


# s2 results

data <- read.csv(file.path("Data", "2-comp", "results_s2.csv"))

ggplot(
    data = data,
    mapping = aes(x = separation_distance, y = Relative_concentration)
) +
    geom_hline(yintercept = 0.5, linetype = "longdash", color="darkgrey") +
    geom_point(stat = "summary", fun="mean", color="#74036a", width=0.00025) +
    geom_path(stat = "summary", fun="mean", color="#74036a") +
    stat_summary(fun.data = "mean_se", geom="errorbar", width=0.25) +
    labs(
        x = "D",
        y = "Rc",
    ) +
    coord_cartesian(xlim=c(1,9), ylim=c(0, 1)) + theme_bw()
    ggsave(file.path("Analysis", "Figures", "comp_graph_s2.svg"), width = 6, height = 4)



########## Merged results #########

s1_d <- read.csv(file.path("Data", "2-comp", "results_s1.csv"))
s2_d <- read.csv(file.path("Data", "2-comp", "results_s2.csv"))
data <- rbind(s1_d, s2_d)

ggplot(
    data = data,
    mapping = aes(x = separation_distance, y = Relative_concentration)
) +
    geom_hline(yintercept = 0.5, linetype = "longdash", color="darkgrey") +

    geom_point(stat = "summary", fun="mean", width=0.00025, mapping = aes(color=factor(sites))) +
    geom_path(stat = "summary", fun="mean", mapping = aes(color=factor(sites))) +

    stat_summary(fun.data = "mean_se", geom="errorbar", width=0.25, mapping = aes(color=factor(sites))) +

    labs(
        x = "D",
        y = "Rc",
        color="Rings",
    ) +
    #facet_wrap(~sites) +
    coord_cartesian(xlim=c(1,9), ylim=c(0, 1)) +
    scale_x_continuous(breaks=c(1,3,5,7,9)) +
    theme_bw()
    ggsave(file.path("Analysis", "Figures", "comp_graphs_merged.svg"), width = 6, height = 4)
