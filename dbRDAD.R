# Appendix to:
# Legendre, P. (2014) Interpreting the replacement and richness difference
# components of beta diversity. Global Ecology and Biogeography, 23, xxx–xxx.
# Appendix S4
#
# R function to compute the dbRDA F-test of significance between response data
# represented by a Euclidean or non-Euclidean dissimilarity matrix and a matrix
# of explanatory variables, following McArdle and Anderson (2001).
dbRDA.D <- function(D, X, nperm=999, option=3, compute.eig=FALSE, coord=FALSE,
                    rda.coord=2, positive.RDA.values=FALSE)
  #
  # Description --
  #
  # Compute the dbRDA F-test of significance between response data represented by
  # a Euclidean or non-Euclidean dissimilarity matrix and a matrix of explanatory
  # variables, using the method of McArdle and Anderson (2001).
  #
  # Usage --
  #
  # dbRDA.D(D, X, nperm=999, option=3, compute.eig=FALSE, coord=FALSE,
  # rda.coord=2, positive.RDA.values=FALSE)
#
# Arguments --
#
# D : Distance matrix representing the response data. D may be non-Euclidean.
# X : Matrix of explanatory variables for the RDA, class 'data.frame' or
# 'matrix'. Factors must be recoded as dummy variables or Helmert contrasts.
# nperm : Number of permutations for the test of significance.
#
# option=1 : Original McArdle-Anderson (2001) equation 4. Slow, not recommended.
# option=2 : McArdle-Anderson equation, simplified.
# option=3 : Least-squares after orthogonalizing X.
# SSY = sum(diag(G)), where G is the Gower-centred distance matrix,
# SSYhat = sum(diag(H %*% G %*% H)), where H is the projector matrix.
#
# Option=1 -- The original F statistic of McArdle and Anderson (2001), eq. 4:
# F = SSYhat / sum(diag(I.minus.H %*% G %*% I.minus.H))
# Degrees of freedom are added to this equation in the output list.
# Option=2 -- Simplified equation:
# F = SSYhat/(SSY-SSYhat)
# Option=3 -- Orthogonalize matrix X by PCA before computing H. No inversion.
# Compute SSYhat as above, then F = SSYhat/(SSY-SSYhat)
# Opt. 2 and 3 are equivalent; they require half the computing time of option 1.
#
# compute.eig=TRUE : the eigenvalues and eigenvectors of D are computed.
# => Do NOT use with very large matrices (slow).
# coord=TRUE : compute the principal coordinates corresponding to the
# positive eigenvalues of D. Requires that compute.eig=TRUE.
# rda.coord : Number of RDA ordination coordinates to compute, for example 2.
# positive.RDA.values=TRUE : store only positive RDA eigenvalues in output list.
# =FALSE: store all RDA eigenvalues in output list.
#
#
# Details --
#
# Compute the dbRDA F-test of significance. The response is represented by a
# Euclidean or non-Euclidean dissimilarity matrix; X is a matrix of explanatory
# variables, as in regular RDA.
#
# The F-statistic is obtained without prior computation of the eigenvalues and
# eigenvectors of the dissimilarity matrix, hence no correction has to be made
# to eliminate the negative eigenvalues. Three computation methods are
# available, all derived from McArdle and Anderson (2001).
#
# The eigenvalues and eigenvectors of D are computed if compute.eig=TRUE.
# If coord=TRUE, the principal coordinates corresponding to the positive
# eigenvalues of D are computed.
#
# The function may fail to produce a meaningful RDA test of significance and
# ordination axes if D is extremely non-Euclidean. This is the case with some
# forms of genomic distances.
#
# Value --
#
# F : F-statistic.
# Rsquare : R-square and adjusted R-square statistics.
# P.perm : Permutational p-value of RDA R-square (test based on F).
# SS.total : Trace of matrix G, equal to the total sum of squares of Y and the
# sum of the eigenvalues of D.
# PCoA.values : Eigenvalues (if they are computed, i.e. if compute.eig=TRUE).
# PCoA.vectors : Principal coordinates for the positive eigenvalues of D.
# RDA.values : RDA eigenvalues.
# RDA.rel.values : RDA relative eigenvalues.
# RDA.cum.values : RDA cumulative relative eigenvalues.
# RDA.coord : Ordination coordinates of objects on selected RDA axes.
#
# References --
#
# Legendre, P. (2014) Interpreting the replacement and richness difference
# components of beta diversity. Global Ecology and Biogeography, 23, xxx–xxx.
#
# Legendre, P. & Legendre, L. (2012) Numerical ecology, 3rd English edition.
# Elsevier Science BV, Amsterdam.
#
# McArdle, B.H. & Anderson, M.J. (2001) Fitting multivariate models to
# community data: a comment on distance-based redundancy analysis.
# Ecology, 82, 290–297.
#
# Example -- Six sites from the mite data available in the vegan package.
#
# library(vegan)
# Load function dbRDA.D()
# data(mite)
# data(mite.env)
# sel = c(14,24,31,41,49,64)
# mite.BC = vegdist(mite[sel,], "bray") # Two negative eigenvalues
# res = dbRDA.D(mite.BC, mite.env[sel,1:2], nperm=999, compute.eig=TRUE)
# plot(res$RDA.coord)
# text(res$RDA.coord, labels=rownames(mite.env[sel6,]), pos=3)
#
# License: GPL-2
# Author:: Pierre Legendre, March 2013
{
  D <- as.matrix(D)
  X <- as.matrix(X)
  n <- nrow(D)
  epsilon <- .Machine$double.eps
  #
  # Gower centring, matrix formula. Legendre & Legendre (2012), equation 9.42
  One <- matrix(1,n,n)
  mat <- diag(n) - One/n
  G <- -0.5 * mat %*% (D^2) %*% mat
  SSY <- sum(diag(G))
  # LCBD <- diag(G)
  #
  # Principal coordinate analysis after eigenvalue decomposition of D
  if(compute.eig) {
    eig <- eigen(G, symmetric=TRUE)
    values <- eig$values # All eigenvalues
    vectors <- eig$vectors # All eigenvectors, scaled to lengths 1
    if(coord) {
      select <- which(values > epsilon)
      princ.coord <- vectors[,select] %*% diag(sqrt(values[select]))
    } else { princ.coord <- NA }
  } else {
    values <- princ.coord <- NA
  }
  #
  # Compute projector matrix H ("hat" matrix in the statistical literature)
  X.c <- scale(X, center=TRUE, scale=FALSE) # Centre matrix X
  m <- qr(X.c, tol=1e-6)$rank # m = rank of X.c
  cat("Rank of X centred =",m,"\n")
  if(m==1) {
    H <- (X.c[,1] %*% t(X.c[,1]))/((t(X.c[,1]) %*% X.c[,1])[1,1])
  } else {
    if(option<3) {
      # if(det(t(X.c)%*%X.c)<epsilon) stop ('Collinearity detected in X')
      if(m < ncol(X.c)) stop ('Collinearity detected in X')
      H <- X.c %*% solve(t(X.c) %*% X.c) %*% t(X.c)
      #
      # option=3: compute projector H from orthogonalized X; no inversion
    } else {
      X.eig <- eigen(cov(X.c))
      k <- length(which(X.eig$values > epsilon))
      X.ortho <- X.c %*% X.eig$vectors[,1:k] # F matrix of PCA
      XprX <- t(X.ortho) %*% X.ortho
      H <- X.ortho %*% diag(diag(XprX)^(-1)) %*% t(X.ortho)
    }
  }
  #
  # Compute the F statistic: McArdle & Anderson (2001), equation 4 modified
  HGH <- H %*% G %*% H
  SSYhat <- sum(diag(HGH))
  #
  if(option==1) {
    I.minus.H <- diag(n) - H
    den1 <- sum(diag(I.minus.H %*% G %*% I.minus.H))
    F <- SSYhat/den1 # F statistic without the degrees of freedom
    Rsquare <- F/(F+1)
  } else {
    F <- SSYhat/(SSY-SSYhat) # F statistic without the degrees of freedom
    Rsquare <- SSYhat/SSY # or equivalent: Rsquare <- F/(F+1)
  }
  RsqAdj <- 1-((1-Rsquare)*(n-1)/(n-1-m))
  #
  # Permutation test of F
  if(nperm > 0) {
    nGE=1
    for(i in 1:nperm) {
      order <- sample(n)
      Gperm <- G[order, order]
      H.Gperm.H <- H %*% Gperm %*% H
      SSYhat.perm <- sum(diag(H.Gperm.H))
      #
      if(option==1) {
        den <- sum(diag(I.minus.H %*% Gperm %*% I.minus.H))
        F.perm <- SSYhat.perm/den
      } else {
        F.perm <- SSYhat.perm/(SSY-SSYhat.perm)
      }
      if(F.perm >= F) nGE=nGE+1
    }
    P.perm <- nGE/(nperm+1)
  } else { P.perm <- NA }
  #
  # Compute RDA ordination coordinates
  if(rda.coord > 0) {
    HGH.eig <- eigen(HGH, symmetric=TRUE)
    # kk <- length(which(HGH.eig$values > epsilon))
    RDA.values <- HGH.eig$values
    rel.eig <- RDA.values/SSY
    cum.eig <- cumsum(rel.eig)
    kk <- length(which(rel.eig > epsilon))
    if(positive.RDA.values) {
      RDA.values <- RDA.values[1:kk]
      rel.eig <- rel.eig[1:kk]
      cum.eig <- cum.eig[1:kk]
    }
    k <- min(rda.coord, kk)
    if(k >= 2) {
      RDA.coord <-sweep(HGH.eig$vectors[,1:k],2,sqrt(RDA.values[1:k]),FUN="*")
    } else {
      RDA.coord <- NA
      cat("k =",k," -- Fewer than two RDA eigenvalues > 0\n")
    }
  } else { RDA.values <- rel.eig <- cum.eig <- RDA.coord <- NA }
  #
  list(F=F*(n-m-1)/m, Rsquare=c(Rsquare,RsqAdj), P.perm=P.perm, SS.total=SSY,
       PCoA.values=values, PCoA.vectors=princ.coord, RDA.values=RDA.values/(n-1),
       RDA.rel.values=rel.eig, RDA.cum.values=cum.eig, RDA.coord=RDA.coord)
}