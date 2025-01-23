# Code for template plots of CSF proteomic module analysis, LM11A-31 Phase 2a follow-up manuscript
# Figure 3C - Placebo vs Drug comparisons of Module changes
# Written by Hayley R.C. Shanks. Last update: Jan 23, 2025
################################################################################
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(gglot2, ggpubr, ggbeeswarm, ggsignif)
csf_module_template_plot <- function(df, var, to_compare, c, w, h, flname, ylabel, plt_title, font_sz){
  # df (data frame): data frame with the variables to be plotted (requires a column, "Group"
  # with grouping variable of interest)
  # var (str): name of the variable within df to be plotted (continuous)
  # c (numeric): cex value for ggbeeswarm, indicating spread of the dots. between [0,1]
  # w (numeric): width of final plot in units
  # h (numeric): height of final plot in units
  # flname (str): output filename, including path (if different from current wd)
  # ylabel (str): y axis label
  # plt_title (str): title of plot
  # font_sz (numeric): font size for final figure
  # compute y axis range so the y limits can be adjusted to ensure geom_signif is not cut off by plot title
  ymin <- min(df[[var]])
  ymax <- max(df[[var]])
  increment <- (ymax-ymin)*0.1
  p <- ggplot(df, aes(fill=Group, y=.data[[var]], x=Group)) +
    geom_hline(yintercept =0, linetype = "dashed", linewidth = 0.25) +
    geom_beeswarm(aes(group=Group, color = Group), cex = c, size =0.75, show.legend =  FALSE) + 
    geom_violin(aes(fill=Group, y=.data[[var]], x=Group, color = Group), alpha = 0.5, scale = "count", show.legend =  FALSE) +
    geom_boxplot(width=0.2, fill = NA, show.legend =  FALSE, outlier.colour = NA, notch = TRUE) + 
    # add the between group significance. will default to 2 sided
    geom_signif(comparisons = to_compare,  na.rm = TRUE, test = "wilcox.test", 
                textsize=1.85, map_signif_level = function(p) sprintf("p = %.2g", p)) +
    xlab("") +  ylab(ylabel) + ggtitle(plt_title) +  theme_classic() + 
    ylim(ymin-increment, ymax+increment) + 
    theme(text = element_text(size = font_sz, color ="black")) +  
    theme(axis.text.y = element_text(size=font_sz, color ="black")) +  
    theme(axis.text.x = element_text(size=font_sz, color ="black")) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) 
  ggsave(flname, plot=p, width = w, height = h, dpi = 300)
}
