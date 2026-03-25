# Analyzer for well-mixed simulations

library(tidyverse)
library(patchwork)
library(plotrix)

data <- read.csv(file.path("Data", "1-well-mixed", "wellm_results.csv"))

# Something to check standard error bars
serrors <- aggregate(Relative_concentration ~ prop_break_rate, data = data,
        FUN = function(x) c(mean = mean(x), se = std.error(x)))
# Join this to the rest of the data


ggplot(
    data = data,
    mapping = aes(x = prop_break_rate, y = Relative_concentration)
) +
    geom_hline(yintercept = 0.5, linetype = "longdash", color="darkgrey") +
    geom_point(stat = "summary", fun="mean", color="#0005D5", width=0.0005) +
    #geom_point(color="#0005D5") +
    geom_path(stat = "summary", fun="mean", color="#0005D5") +
    stat_summary(fun.data = "mean_se", geom="errorbar", width=0.00025) +
    #geom_errorbar(stat="summary",
    #fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},
    #fun.ymax=function(x) {mean(x)+sd(x)/sqrt(length(x))},
    #width = 0.01, linewidth = 0.5) +
    labs(
        x = "Propagule Formation Rate",
        y = "Rc",
    ) +
    xlim(0,0.05) + ylim(0, 0.5) + theme_bw()
    ggsave(file.path("Analysis", "Figures", "wellm_graph.svg"), width = 4, height = 4)
