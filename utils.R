make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))
                                                
'%!in%' <- function(x,y)!('%in%'(x,y))
                                                
se <- function(x){sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))}
                                                log10p <- function(x) log10(1+x)
plot_chull<-function(xcoord, ycoord, lcolor, llty){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor, lty = llty)
}  
                                                
poly <- function(x, upper, lower, fill){
  polygon(x=c(x, rev(x), x[1]), y=c(upper, rev(lower), upper[1]),border=NA,col=fill)
}
                                                
mean.nona <- function(x){mean(x, na.rm=T)}                                                
