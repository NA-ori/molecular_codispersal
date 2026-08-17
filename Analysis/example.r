# Plot some example data from a simulation

library(tidyverse)
library(patchwork)
library(plotrix)
library(svglite)

data <- read.csv(file.path("Data/raw/2-comp/curves_s2/complete_curve_99276178.csv"))
data <- subset(data, pixel==0)
data <- data[c("t", "sp_r1_1_ad", "sp_r1_2_ad", "sp_r1_3_ad", "sp_r2_1_ad", "sp_r2_2_ad", "sp_r2_3_ad")]
data <- data %>%
    pivot_longer(cols= c("sp_r1_1_ad", "sp_r1_2_ad", "sp_r1_3_ad", "sp_r2_1_ad", "sp_r2_2_ad", "sp_r2_3_ad"),
    names_to="species", values_to="concentration")

ggplot(
    data = data,
    mapping = aes(x = t, y = concentration)
) +
    geom_point(stat = "summary", fun="mean", width=1, mapping=aes(color=factor(species), shape=factor(species))) +
    geom_path(stat = "summary", fun="mean", mapping=aes(color=factor(species))) +
    #guides(shape="none") +
    labs(
        x = "Time",
        y = "Count",
        color = "Species", shape = "Species"
    ) +
    coord_cartesian(xlim = c(0, 1000), ylim = c(0, 1000)) +
    scale_x_continuous(breaks = c(0, 250, 500, 750, 1000)) +
    scale_color_manual(values = c("#D62800", "#FF9B56", "#D462A6", "#A40062", "#5BCFFB", "#F5ABB9")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
    ggsave(file.path("Analysis", "Figures", "exampe_graph.svg"), width = 6, height = 3)


# Make color palette for example figure
cols <- colorRampPalette(c("#DC3220", "#005AB5"))
colorlist <- cols(100)

svglite("gradient.svg", width=4, height=4)
plot(rep(1,100),col=(colorlist), pch=15,cex=2)
dev.off()
