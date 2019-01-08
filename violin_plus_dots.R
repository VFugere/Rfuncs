## GOAL:
## re-create a figure similar to Fig. 2 in Wilson et al. (2018), 
## Nature 554: 183-188. Available from:
## https://www.nature.com/articles/nature25479#s1
##
## combines a boxplot (or violin) with the raw data, by splitting each
## category location in two (box on left, raw data on right)


# initial set-up ----------------------------------------------------------

## load source code
source("flatviolin.R")

## set plotting theme
theme_set(theme_bw())

## import data
iris <- iris
skim(iris)


# half violin plot with raw data ------------------------------------------

## create a violin plot of Sepal.Length per species
## using the custom function geom_flat_violin()

ggplot(data = subset(fm, weeks == 3), 
       mapping = aes(x = leaf.deployment, y = prop.decomp, fill = leaf.deployment)) + 
  geom_flat_violin(scale = "count", trim = FALSE) + 
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "pointrange", position = position_nudge(0.05)) + 
  geom_dotplot(binaxis = "y", dotsize = 0.5, stackdir = "down", binwidth = 0.1, 
               position = position_nudge(-0.025)) + 
  theme(legend.position = "none") + 
  labs(x = "Leaf deployment", y = "Proportion decomposed")

