# Analyzer for disturbance simulations

library(tidyverse)
library(patchwork)
library(plotrix)

data <- read.csv(file.path("Data", "3-dist", "dist_results_s1.csv"))

# s1 results

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
        y = bquote(A[rc]),
        color = quote(lambda), shape = quote(lambda)
    ) +

    coord_cartesian(xlim=c(1, 9), ylim=c(0, 1)) +
    scale_x_continuous(breaks = c(1, 3, 5, 7, 9)) +
    scale_color_manual(values = c("#A40062", "#5BCFFB", "#D462A6", "#D62800")) +
    theme_bw() +
    theme(legend.position = c(0.9, 0.2), panel.grid.minor = element_blank(),
    legend.background = element_rect(colour="grey50", linewidth=0.3))
    ggsave(file.path("Analysis", "Figures", "dist_graph_s1.svg"), width = 6, height = 4)


# s2 results

data <- read.csv(file.path("Data", "3-dist", "dist_results_s2.csv"))

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
        y = bquote(A[rc]),
        color = quote(lambda), shape = quote(lambda)
    ) +

    coord_cartesian(xlim = c(1, 9), ylim = c(0, 1)) +
    scale_x_continuous(breaks = c(1, 3, 5, 7, 9)) +
    scale_color_manual(values = c("#A40062", "#5BCFFB", "#D462A6", "#D62800")) +
    theme_bw() +
    theme(legend.position = c(0.9, 0.2), panel.grid.minor = element_blank(),
    legend.background = element_rect(colour="grey50", linewidth=0.3))
    ggsave(file.path("Analysis", "Figures", "dist_graph_s2.svg"), width = 6, height = 4)


# Merged results

s1_d <- read.csv(file.path("Data", "3-dist", "dist_results_s1.csv"))
s2_d <- read.csv(file.path("Data", "3-dist", "dist_results_s2.csv"))
data <- rbind(s1_d, s2_d)

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
        y = bquote(A[rc]),
        color = quote(lambda), shape = quote(lambda)
    ) +

    coord_cartesian(xlim = c(1, 9), ylim = c(0, 1)) +
    scale_x_continuous(breaks = c(1, 3, 5, 7, 9)) +
    scale_color_manual(values = c("#A40062", "#5BCFFB", "#D462A6", "#D62800")) +
    scale_shape_manual(values = c(19, 17, 15, 18)) +
    theme_bw() +
    theme(legend.position = c(0.93, 0.22), panel.grid.minor = element_blank(),
    legend.background = element_rect(colour="grey50", linewidth=0.3)) +
    facet_wrap(~sites)
    ggsave(file.path("Analysis", "Figures", "dist_graph_merged.svg"), width = 8, height = 4)
