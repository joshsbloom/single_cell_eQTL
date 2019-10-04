# parallelized standardise function
standardise2 <- function (x, center = TRUE, scale = TRUE) {
  if ( center & scale ) {
    y <- t(x) - Rfast::colmeans(x,parallel=T)
    y <- y / sqrt( Rfast::rowsums(y^2) ) * sqrt( (dim(x)[1] - 1) )
	y <- t(y)
  } else if ( center & !scale ) {
    m <- Rfast::colmeans(x,parallel=T)
    y <- eachrow(x, m, oper ="-" )
  } else if ( !center & scale ) {
    s <- Rfast::colVars(x, std = TRUE,parallel=T)
    y <- eachrow(x, s, oper = "/")
  }
  y
} 

chroms=c('I', 'II', 'III', 'IV', 'V', 'X')

acomb1 <- function(...) abind(..., along=1)
