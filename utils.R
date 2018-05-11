multipan <- function(r,c){
  par(mfrow=c(r,c))
  par(mar=c(4,5,1,1))
  print('dimension feuille en excluant les marges = 15.23cm X 22.85cm ou 6in X 9in')
}

outplot <- function(x,p=F,id=F){
  plot(x=x,y=1:length(x),pch=1,col='grey',xlab='value',ylab='index')
  abline(v=mean(x),lwd=3)
  abline(v=(mean(x)-sd(x)),lty=2,lwd=2)
  abline(v=(mean(x)+sd(x)),lty=2,lwd=2)
  abline(v=(mean(x)-2*sd(x)),lty=2)
  abline(v=(mean(x)+2*sd(x)),lty=2)
  if (p)
  {
    sx = scale(x)
    sx = abs(sx)
    i = which(sx > 2)
    r = data.frame('row' = i,'value' = x[i])
    print(r)
    rm(sx,i)
  }
  if (id)
    identify(x=x,y=1:length(x))
}

spmat <- function (x, c=T) {
  if (c)
  pairs(x, lower.panel = panel.smooth, upper.panel = panel.cor, diag.panel = panel.hist)
  else
  pairs(x, lower.panel = panel.smooth, diag.panel = panel.hist)
}

fitline <- function(model,ci=T,glm=F,t=1,c='black',w=1){
  lmTEMP <- model
  if (glm)
  LMTpred <- predict(lmTEMP, se.fit=T, type='response')
  else
  LMTpred <- predict(lmTEMP, se=T)
  ovTEMP <- order(lmTEMP$model[,1])
  lines(lmTEMP$model[,1][ovTEMP], LMTpred$fit[ovTEMP],lty=t,col=c,lwd=w)
  if (ci){
  lines(lmTEMP$model[,1][ovTEMP], LMTpred$fit[ovTEMP]+2*LMTpred$se[ovTEMP],lty=2,col=c)
  lines(lmTEMP$model[,1][ovTEMP], LMTpred$fit[ovTEMP]-2*LMTpred$se[ovTEMP],lty=2,col=c)
  rm(lmTEMP,LMTpred,ovTEMP)
  }
}

panel.cor <- function(x, y, method="pearson", digits=3, cex.cor=1.2, no.col=FALSE)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, method=method)
  ra <- cor.test(x, y, method=method)$p.value
  txt <- round(r, digits)
  prefix <- ""
  if(ra <= 0.1) prefix <- "."
  if(ra <= 0.05) prefix <- "*"
  if(ra <= 0.01) prefix <- "**"
  if(ra <= 0.001) prefix <- "***"
  if(no.col)
  {
    color <- 1
    if(r < 0) { if(ra <= 0.001) sig <- 4 else sig <- 3 }
    else { if(ra <= 0.001) sig <- 2 else sig <- 1 }
  }
  else
  {
    sig <- 1
    if(ra <= 0.001) sig <- 2
    color <- 2
    if(r < 0) color <- 4
  }
  txt <- paste(txt, prefix, sep="\n")
  text(0.5, 0.5, txt, cex = cex.cor, font=sig, col=color)
}

panel.hist <- function(x, no.col=FALSE, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  his <- hist(x, plot=FALSE)
  breaks <- his$breaks; nB <- length(breaks)
  y <- his$counts
  y <- y/max(y)
  if(no.col) rect(breaks[-nB], 0, breaks[-1], y, col="gray", ...)
  else rect(breaks[-nB], 0, breaks[-1], y, col="dark grey", ...)
}

'%!in%' <- function(x,y)!('%in%'(x,y))

fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
 
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
 
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
 
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
 
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
 
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
 
  x1 <- switch(pos,
    topleft     =x[1] + sw, 
    left        =x[1] + sw,
    bottomleft  =x[1] + sw,
    top         =(x[1] + x[2])/2,
    center      =(x[1] + x[2])/2,
    bottom      =(x[1] + x[2])/2,
    topright    =x[2] - sw,
    right       =x[2] - sw,
    bottomright =x[2] - sw)
 
  y1 <- switch(pos,
    topleft     =y[2] - sh,
    top         =y[2] - sh,
    topright    =y[2] - sh,
    left        =(y[1] + y[2])/2,
    center      =(y[1] + y[2])/2,
    right       =(y[1] + y[2])/2,
    bottomleft  =y[1] + sh,
    bottom      =y[1] + sh,
    bottomright =y[1] + sh)
 
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
 
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}

detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}