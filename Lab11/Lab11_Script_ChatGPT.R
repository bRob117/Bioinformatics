# Bryce Robinson
# BSC4434C Bioinformatics
# Lab 11 BEAST2 Tree Attempt
#
# NOTE:
# ChatGPT was used to help generate and troubleshoot this script.
# My original plan was to use FigTree, but I ran into problems getting it
# to run correctly with my Java setup. Because of that, I recreated the
# Figure 21 visualization in R using the BEAST output file.
#
# The main goal of this script is to reproduce the same information that
# FigTree displays:
# - the topology of the MCC tree
# - posterior support values at internal nodes
# - 95% HPD intervals for node ages
# - a time axis measured as time before present

# Set working directory so R can find the tree file and save the image
setwd("~/Documents/FGCU/BSC4434C/GitHub/Lab11")

# Load packages
# ggplot2 = draws the final plot
# ggtree   = calculates the tree layout (where nodes/tips go on the page)
# treeio   = reads BEAST tree files with metadata
# ape      = tree manipulation functions like ladderize()
# tibble   = helps store tree data in table form
# dplyr    = used for joins, filters, and data manipulation
library(ggplot2)
library(ggtree)
library(treeio)
library(ape)
library(tibble)
library(dplyr)

# -------------------------------------------------------------------
# 1. Read in the BEAST tree
# -------------------------------------------------------------------
# I use read.beast() instead of read.tree() because BEAST output includes
# extra information attached to nodes, such as posterior probabilities
# and HPD intervals. A simpler tree reader would mostly just read the
# branching structure and labels.
beast_tree <- read.beast("Primates.MCC.trees")

# -------------------------------------------------------------------
# 2. Reorder the nodes to match the FigTree display
# -------------------------------------------------------------------
# ladderize() changes the visual order of branches, but it does NOT change
# the actual phylogeny or relationships among taxa.
# I used this because the instructions for Figure 21 say to order nodes.
beast_tree@phylo <- ape::ladderize(beast_tree@phylo, right = FALSE)

# -------------------------------------------------------------------
# 3. Extract the BEAST metadata into a table
# -------------------------------------------------------------------
# as_tibble() turns the tree object into a table where each row is a node or tip.
# This table includes metadata like:
# - node number
# - parent node
# - posterior support
# - HPD intervals for node height
#
# Important idea:
# This table has the statistical/biological info, but it does NOT yet have
# all the plotting coordinates needed to draw the tree in the exact layout.
td <- as_tibble(beast_tree)

# -------------------------------------------------------------------
# 4. Pull apart the 95% HPD intervals
# -------------------------------------------------------------------
# In this file, the HPD interval for node age is stored as a list with
# two numbers for each row:
#   [lower bound, upper bound]
#
# I split those into separate columns because it is much easier to plot
# a horizontal bar when the start and end are in separate numeric columns.
#
# hpd_low  = younger end of the node-age interval
# hpd_high = older end of the node-age interval
td$hpd_low <- sapply(td$height_0.95_HPD, function(z) {
  if (is.numeric(z) && length(z) == 2) z[1] else NA_real_
})

td$hpd_high <- sapply(td$height_0.95_HPD, function(z) {
  if (is.numeric(z) && length(z) == 2) z[2] else NA_real_
})

# -------------------------------------------------------------------
# 5. Get plotting coordinates from ggtree
# -------------------------------------------------------------------
# ggtree() calculates where each node and tip should appear on the plot.
# The object p$data contains the plotting layout, especially x and y positions.
#
# Important difference:
# - td has metadata like posterior and HPD
# - p$data has plotting coordinates like x and y
#
# I need both later, so I save the layout table here.
p <- ggtree(beast_tree, linewidth = 0.4)
layout_df <- p$data

# -------------------------------------------------------------------
# 6. Join metadata and plotting coordinates together
# -------------------------------------------------------------------
# This is one of the most important parts of the whole script.
#
# Why join them?
# - td alone is not enough because it does not contain plotting positions
# - layout_df alone is not enough because it does not contain all BEAST metadata
#
# So I join by node number so each node has:
# - where it goes on the plot (x, y)
# - its posterior value
# - its HPD interval
plot_df <- left_join(
  layout_df,
  td %>% select(node, parent, posterior, hpd_low, hpd_high),
  by = "node",
  suffix = c("", ".meta")
)

# Sometimes the join creates two parent columns.
# If that happens, I keep one clean version.
if ("parent.meta" %in% names(plot_df)) {
  plot_df$parent <- plot_df$parent.meta
  plot_df$parent.meta <- NULL
}

# -------------------------------------------------------------------
# 7. Convert x positions to "time before present"
# -------------------------------------------------------------------
# In the raw layout from ggtree, x is just the plotted horizontal position.
# The tips are at the far right side of the plot.
#
# I want the present to be 0, like in the FigTree display.
# So I subtract the maximum x value from every x coordinate.
#
# Visual meaning of this step:
# - all tips end at 0
# - older nodes move left into negative values
# - the x-axis now behaves like time before present
max_x <- max(plot_df$x, na.rm = TRUE)
plot_df$x_plot <- plot_df$x - max_x

# -------------------------------------------------------------------
# 8. Flip the vertical order
# -------------------------------------------------------------------
# The same tree can be drawn with taxa in different top-to-bottom orders
# without changing the actual relationships.
#
# I reversed the y positions so the tree looks more like the tutorial figure.
max_y <- max(plot_df$y, na.rm = TRUE)
plot_df$y_plot <- max_y - plot_df$y + 1

# -------------------------------------------------------------------
# 9. Build the branch segments manually
# -------------------------------------------------------------------
# I originally hoped ggtree would do all the drawing directly, but I ended up
# building the branch segments myself because I wanted more control over the
# final appearance and over how the extra metadata lined up with the tree.
#
# The idea is:
# - every node has a parent
# - if I know the node position and the parent position,
#   I can draw the branch between them
#
# First, make a table of parent coordinates
parent_coords <- plot_df %>%
  select(node, x_plot, y_plot) %>%
  rename(parent = node, x_parent = x_plot, y_parent = y_plot)

# Join each node to its parent's coordinates
branch_df <- plot_df %>%
  left_join(parent_coords, by = "parent") %>%
  filter(!is.na(parent))

# Horizontal segments:
# These run from the parent x-position to the child x-position
# at the child's y-position.
#
# In a phylogenetic tree, the horizontal part mainly shows branch length / time.
hseg <- branch_df %>%
  transmute(
    x = x_parent,
    xend = x_plot,
    y = y_plot,
    yend = y_plot
  )

# Vertical segments:
# These connect branches up and down at the parent split.
#
# In a phylogenetic tree, the vertical part is just the connector showing
# the branching event, not extra time.
vseg <- branch_df %>%
  transmute(
    x = x_parent,
    xend = x_parent,
    y = y_parent,
    yend = y_plot
  )

# -------------------------------------------------------------------
# 10. Prepare HPD bars for internal nodes
# -------------------------------------------------------------------
# The blue bars in FigTree represent the 95% HPD interval for node age.
#
# 95% HPD means the range of node ages that contains 95% of the highest
# posterior density from the Bayesian analysis. It is similar in spirit to
# an uncertainty interval, but specifically comes from the posterior distribution.
#
# These only make sense for internal nodes because internal nodes represent
# divergence events with estimated ages. The tips are present-day taxa, so
# they do not get the same kind of node-age bar here.
#
# Why negative values?
# Since the x-axis was converted to "time before present", older ages should
# appear farther left. That means older values must be more negative.
#
# So:
# - the left end of the bar is -hpd_high   (older bound = farther left)
# - the right end of the bar is -hpd_low   (younger bound = closer to 0)
bar_df <- plot_df %>%
  filter(!isTip, !is.na(hpd_low), !is.na(hpd_high)) %>%
  mutate(
    bar_xmin = -hpd_high,
    bar_xmax = -hpd_low
  )

# Internal node labels:
# posterior values tell how strongly that node/clade is supported
# across the posterior sample of trees
node_df <- plot_df %>%
  filter(!isTip, !is.na(posterior))

# Tip labels:
# these are the species names at the ends of the tree
tip_df <- plot_df %>%
  filter(isTip)

# -------------------------------------------------------------------
# 11. Draw the final figure
# -------------------------------------------------------------------
# I build the plot in layers:
# 1. horizontal branches
# 2. vertical branch connectors
# 3. HPD bars
# 4. posterior labels
# 5. tip labels
#
# This is basically reproducing what FigTree shows, but manually in ggplot2.
fig21 <- ggplot() +
  
  # Horizontal branch segments
  geom_segment(
    data = hseg,
    aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.4,
    color = "black"
  ) +
  
  # Vertical connector segments
  geom_segment(
    data = vseg,
    aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.4,
    color = "black"
  ) +
  
  # Blue HPD bars for node-age uncertainty
  geom_segment(
    data = bar_df,
    aes(x = bar_xmin, xend = bar_xmax, y = y_plot, yend = y_plot),
    linewidth = 0.4,
    color = "blue",
    alpha = 0.7
  ) +
  
  # Posterior values at internal nodes
  # posterior is the support for that node/clade in the posterior sample
  geom_text(
    data = node_df,
    aes(
      x = x_plot,
      y = y_plot,
      label = format(round(as.numeric(posterior), 4), nsmall = 0)
    ),
    size = 2.5,
    hjust = 0.5,
    vjust = -0.3
  ) +
  
  # Species names at the tips
  geom_text(
    data = tip_df,
    aes(x = 0.6, y = y_plot, label = label),
    hjust = 0,
    size = 3.2
  ) +
  
  # X-axis as time before present
  scale_x_continuous(
    name = "Time before present",
    breaks = seq(-100, 0, by = 25),
    limits = c(min(bar_df$bar_xmin, na.rm = TRUE) - 5, 5)
  ) +
  
  # Remove y-axis clutter so the figure looks more like FigTree
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text.x = element_text(size = 10),
    plot.margin = margin(10, 60, 10, 10)
  ) +
  
  # Allow text to extend slightly outside the plot panel
  coord_cartesian(clip = "off")

# Show plot in R
print(fig21)

# Save image file
ggsave("Primates_MCC_Figure21_match.png", fig21, width = 11, height = 7.5, dpi = 300)

# -------------------------------------------------------------------
# Final summary of the logic
# -------------------------------------------------------------------
# What this script is doing overall:
#
# 1. Read the BEAST tree with all of its metadata.
# 2. Reorder the visual layout to look more like FigTree.
# 3. Extract posterior and HPD information from the tree.
# 4. Get x/y plotting coordinates from ggtree.
# 5. Join metadata and coordinates together by node.
# 6. Shift the x-axis so the tips are at 0 and the tree reads as time before present.
# 7. Rebuild the branch segments manually from node-parent relationships.
# 8. Plot the HPD bars and posterior values on the correct internal nodes.
#
# Main reason it got more complicated than expected:
# the metadata and plotting layout were stored in different places, so I had
# to combine them before I could reproduce the FigTree-style display.
#
# If asked what the plot means:
# - The branching pattern shows the relationships among the primates.
# - The node labels show posterior support for each internal clade.
# - The blue bars show the 95% HPD interval for the estimated age of each node.
# - The x-axis shows time before present.