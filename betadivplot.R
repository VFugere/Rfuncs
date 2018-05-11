mybetadivplot <- function (x, cex=1) 
{
  g <- scores(x, choices = 1:2)
  plot(g$sites, asp = 1, type = "n", ylab='Dimension 2', xlab='Dimension 1', cex.lab=0.7, cex.axis=0.7)
  points(g$sites, pch = c(21,22,21)[site$region], bg=c('white','white','black')[site$region])
  points(g$centroids, pch = 16, cex = 1.5, col = "red")
  for (i in levels(x$group)) {
    ch <- chull(g$sites[x$group == i, ])
    ch <- c(ch, ch[1])
    lines(x$vectors[x$group == i, 1:2][ch, ], col = "black", lty = "dashed")
  }
}