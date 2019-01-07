make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))
'%!in%' <- function(x,y)!('%in%'(x,y))
se <- function(x){sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))}