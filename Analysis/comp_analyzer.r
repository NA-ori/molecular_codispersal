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
    geom_hline(yintercept = 0.5, linetype = "longdash", color = "grey50") +

    stat_summary(fun.data = "mean_se", geom="errorbar", width=0.35, color="#74036a") +

    geom_point(stat = "summary", fun="mean", width=1, color="#74036a") +
    geom_path(stat = "summary", fun="mean", color="#74036a") +

    #guides(shape="none") +
    labs(
        x = "D",
        y = bquote(A[rc]),
    ) +
    #facet_wrap(~sites) +
    coord_cartesian(xlim = c(1, 9), ylim = c(0, 1)) +
    scale_x_continuous(breaks = c(1, 3, 5, 7, 9)) +
    scale_color_manual(values = c("#A40062")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
    ggsave(file.path("Analysis", "Figures", "comp_graph_s2.svg"), width = 4, height = 3)



########## Merged results #########

s1_d <- read.csv(file.path("Data", "2-comp", "results_s1.csv"))
s2_d <- read.csv(file.path("Data", "2-comp", "results_s2.csv"))
data <- rbind(s1_d, s2_d)

ggplot(
    data = data,
    mapping = aes(x = separation_distance, y = Relative_concentration)
) +
    geom_hline(yintercept = 0.5, linetype = "longdash", color = "grey50") +

    stat_summary(fun.data = "mean_se", geom="errorbar", width=0.35, mapping = aes(color=factor(sites))) +

    geom_point(stat = "summary", fun="mean", width=1, mapping = aes(color=factor(sites), shape=factor(sites))) +
    geom_path(stat = "summary", fun="mean", mapping = aes(color=factor(sites))) +

    #guides(shape="none") +
    labs(
        x = "D",
        y = bquote(A[rc]),
        color = "Rings", shape = "Rings"
    ) +
    #facet_wrap(~sites) +
    coord_cartesian(xlim = c(1, 9), ylim = c(0, 1)) +
    scale_x_continuous(breaks = c(1, 3, 5, 7, 9)) +
    scale_color_manual(values = c("#5BCFFB", "#A40062")) +
    theme_bw() +
    theme(legend.position = c(0.9, 0.2), panel.grid.minor = element_blank(),
    legend.background = element_rect(colour="grey50", linewidth=0.3))
    ggsave(file.path("Analysis", "Figures", "comp_graphs_merged.svg"), width = 4, height = 3)


# Facultative Mutualism plot >:D

data <- read.csv(file.path("Data", "fac_results.csv"))
old_data <- read.csv(file.path("Data", "redo.csv"))
data <- rbind(data, old_data)
data <- subset(data, separation_distance <= 5)

ggplot(
    data = data,
    mapping = aes(x = separation_distance, y = Relative_concentration)
) +
    geom_hline(yintercept = 0.5, linetype = "longdash", color = "grey50") +

    stat_summary(fun.data = "mean_se", geom="errorbar", width=0.35, mapping=aes(color=factor(Facultative_modifier))) +

    geom_point(stat = "summary", fun="mean", width=1, mapping=aes(color=factor(Facultative_modifier), shape=factor(Facultative_modifier))) +
    geom_path(stat = "summary", fun="mean", mapping=aes(color=factor(Facultative_modifier))) +

    #guides(shape="none") +
    labs(
        x = "D",
        y = bquote(A[rc]),
        color = "M", shape = "M"
    ) +
    coord_cartesian(xlim = c(1, 9), ylim = c(0, 1)) +
    scale_x_continuous(breaks = c(1, 3, 5, 7, 9)) +
    scale_color_manual(values = c("#D62800", "#FF9B56", "#1AB3FF", "#FFD700", "#A40062")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
    ggsave(file.path("Analysis", "Figures", "fac_graph.svg"), width = 6, height = 3)

# Larger loops plot >:(

data <- read.csv(file.path("Data", "largeloops_results.csv"))
old_data <- read.csv(file.path("Data", "redo.csv"))
data <- rbind(data, old_data)

ggplot(
    data = data,
    mapping = aes(x = separation_distance, y = Relative_concentration)
) +
    geom_hline(yintercept = 0.5, linetype = "longdash", color = "grey50") +

    stat_summary(fun.data = "mean_se", geom="errorbar", width=0.35, mapping=aes(color=factor(Max_cycles))) +

    geom_point(stat = "summary", fun="mean", width=1, mapping=aes(color=factor(Max_cycles), shape=factor(Max_cycles))) +
    geom_path(stat = "summary", fun="mean", mapping=aes(color=factor(Max_cycles))) +

    #guides(shape="none") +
    labs(
        x = "D",
        y = bquote(A[rc]),
        color = "Subcycles", shape = "Subcycles"
    ) +
    coord_cartesian(xlim = c(1, 9), ylim = c(0, 1)) +
    scale_x_continuous(breaks = c(1, 3, 5, 7, 9)) +
    scale_color_manual(values = c("#1AB3FF", "#A40062", "#D62800")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
    ggsave(file.path("Analysis", "Figures", "largeloops_graph.svg"), width = 6, height = 3)


# Check redo

data <- read.csv(file.path("Data", "redo.csv"))

ggplot(
    data = data,
    mapping = aes(x = separation_distance, y = Relative_concentration)
) +
    geom_hline(yintercept = 0.5, linetype = "longdash", color = "grey50") +

    stat_summary(fun.data = "mean_se", geom="errorbar", width=0.35, color="#74036a") +

    geom_point(stat = "summary", fun="mean", width=1, color="#74036a") +
    geom_path(stat = "summary", fun="mean", color="#74036a") +

    #guides(shape="none") +
    labs(
        x = "D",
        y = bquote(A[rc]),
    ) +
    #facet_wrap(~sites) +
    coord_cartesian(xlim = c(1, 9), ylim = c(0, 1)) +
    scale_x_continuous(breaks = c(1, 3, 5, 7, 9)) +
    scale_color_manual(values = c("#A40062")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
    ggsave(file.path("Analysis", "Figures", "redo_check.svg"), width = 4, height = 3)
