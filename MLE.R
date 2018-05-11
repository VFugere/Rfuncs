#Useful likelihood functions

nll.power <- function(pars,xobs,yobs) {
  #Unpack pars
  a <- pars[1]
  b <- pars[2]
  s <- pars[3]
  #Calculate predicted values
  yfit <- a*(xobs^b)
  #Calculate and return NLL
  nll <- -sum(dnorm(x=yobs, mean=yfit, sd=s, log=T))
  return(nll)
}

nll.lm <- function(pars,xobs,yobs) {
  #Unpack pars
  b0 <- pars[1]
  b1 <- pars[2]
  s <- pars[3]
  #Calculate predicted values
  yfit <- b0+b1*xobs
  #Calculate and return NLL
  nll <- -sum(dnorm(x=yobs, mean=yfit, sd=s, log=T))
  return(nll)
}