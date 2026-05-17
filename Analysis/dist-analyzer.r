# Analyzer for disturbance simulations

library(tidyverse)
library(patchwork)
library(plotrix)

data <- read.csv(file.path("Data", "3-dist", "dist_results_s1_complete.csv"))

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

data <- read.csv(file.path("Data", "3-dist", "dist_results_s2_complete.csv"))

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

s1_d <- read.csv(file.path("Data", "3-dist", "dist_results_s1_complete.csv"))
s2_d <- read.csv(file.path("Data", "3-dist", "dist_results_s2_complete.csv"))
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


# Extinction results

ex_data_s1 <- read.csv(file.path("Data", "3-dist", "dist_results_s1_complete.csv"))
ex_data_s2 <- read.csv(file.path("Data", "3-dist", "dist_results_s2_complete.csv"))
ex_data <- rbind(ex_data_s1, ex_data_s2)

# Process extinction data so it is ggplot-usable

ex_data$extinction_total <- 0
ex_data$extinction_adsorbed <- 0
for (row in seq(nrow(ex_data))) {
    # work out adsorbed extinction rates
    if (ex_data$ad_A_extinct[row] == 1) {
        ex_data$extinction_adsorbed[row] = 'adsorbed A extinct'
    } else if (ex_data$ad_NA_extinct [row]== 1) {
        ex_data$extinction_adsorbed[row] = 'adsorbed NA extinct'
    } else if (ex_data$ad_none_extinct[row] == 1) {
        ex_data$extinction_adsorbed[row] = 'neither extinct'
    } else if (ex_data$ad_all_extinct[row] == 1) {
        ex_data$extinction_adsorbed[row] = 'both extinct'
    }
    
    # work out total extinction rates
    if (ex_data$tot_A_extinct[row] == 1) {
        ex_data$extinction_total[row] = 'total A extinct'
    } else if (ex_data$tot_NA_extinct [row]== 1) {
        ex_data$extinction_total[row] = 'total NA extinct'
    } else if (ex_data$tot_none_extinct[row] == 1) {
        ex_data$extinction_total[row] = 'neither extinct'
    } else if (ex_data$tot_all_extinct[row] == 1) {
        ex_data$extinction_total[row] = 'both extinct'
    }

}

# Make different plots for rings and total/adsorbed

# adsorbed

#ggplot(
#    data = subset(ex_data, sites==1),
#    mapping = aes(x = separation_distance, fill = extinction_adsorbed)) + 
#    geom_bar(position = "fill") +
#    facet_wrap(~factor(disturb_freq)) +
#    labs(
#        x = "D",
#        y = "Proportion of Runs",
#        color = "Result"
#    ) +
#    scale_x_continuous(breaks = c(1, 3, 5, 7, 9)) +
#    theme_bw() +
#    theme(panel.grid.minor = element_blank())
#    ggsave(file.path("Analysis", "Figures", "extinction_pies_s1_adsorbed.svg"), width = 8, height = 4)

#ggplot(
#    data = subset(ex_data, sites==2),
 #   mapping = aes(x = separation_distance, fill = extinction_adsorbed)) + 
 #   geom_bar(position = "fill") +
 #   facet_wrap(~factor(disturb_freq)) +
 #   labs(
 #       x = "D",
 #       y = "Proportion of Runs",
 #       color = "Result"
 #   ) +
 #   scale_x_continuous(breaks = c(1, 3, 5, 7, 9)) +
 #   theme_bw() +
 #   theme(panel.grid.minor = element_blank())
 #   ggsave(file.path("Analysis", "Figures", "extinction_pies_s2_adsorbed.svg"), width = 8, height = 4)


# total

ggplot(
    data = subset(ex_data, sites==1),
    mapping = aes(x = separation_distance, fill = extinction_total)) + 
    geom_bar(position = "fill") +
    facet_wrap(~factor(disturb_freq)) +
    labs(
        x = "D",
        y = "Proportion of Runs",
        fill = "Result"
    ) +
    scale_x_continuous(breaks = c(1, 3, 5, 7, 9)) +
    scale_fill_manual(values = c("#A40062", "#5BCFFB", "#D462A6", "darkgrey")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
    ggsave(file.path("Analysis", "Figures", "extinction_pies_s1_total.svg"), width = 8, height = 4)


ggplot(
    data = subset(ex_data, sites==2),
    mapping = aes(x = separation_distance, fill = extinction_total)) + 
    geom_bar(position = "fill") +
    facet_wrap(~factor(disturb_freq)) +
    labs(
        x = "D",
        y = "Proportion of Runs",
        fill = "Result"
    ) +
    scale_x_continuous(breaks = c(1, 3, 5, 7, 9)) +
    scale_fill_manual(values = c("#D62800", "#FF9B56", "#D462A6", "#A40062")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
    ggsave(file.path("Analysis", "Figures", "extinction_pies_s2_total.svg"), width = 8, height = 4)
