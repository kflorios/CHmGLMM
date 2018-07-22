#' Generates the Data
#'
#' It acts on the primitives of the estimation problem and returns the Data
#'@param id     individual number
#'@param times  one-dimensional array of time points
#'@param n      number of individuals
#'@param X      fixed effects design matrix. Usually includes an intercept and time.
#'@param Z      random effects design matrix. Usually includes an intercept and time
#'@param betas  the fixed effects parameters for the intercept and the slope.
#'@param b      the random effects taken from a multivariate normal with mean 0 and variance matrix D
#'@param Q      the number of items. Set this to four.
#'@return Returns the generated Data as data.frame
#'@export
generateData2 <- function(id,times,n,X,Z,betas,b,Q) {
  
  ncz <- ncol(Z)
  Data <- data.frame(id = id, time = rep(times, n))
  for (q in seq_len(Q)) {
    eta <- c(X %*% betas) + rowSums(Z * b[id, c(q*ncz - 1, q*ncz)])
    ###eta <- c(X %*% betas) + rowSums(Z * b[id, c(q*ncz)])
    Data[[paste("y", q, sep = "")]] <- rbinom(length(eta), 1, plogis(eta))
  }
  
  
  
  return(Data)
}