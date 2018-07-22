#' A Demo Example of the CHmGLMM package
#'
#' It sets up an estimation problem, and performs the estimation with methods AVE, DWAVE, WAVE, CH-EXP and CH-ECDF
#'@return Returns the demo example output. It runs as a script with no input.
#'@export
demoExample2 <- function() {
  
  library(MASS)
  
  Scenario <- 1
  M <- 1
  
  n <- 100  # sample size
  Q <- 4 # number of longitudinal binary outcomes
  P <- choose(Q,2) # number of pairs of items
  times <<- seq(0, 5, length.out = 6)
  times <<- times - (0+times[length(times)])/2    # center around zero
  id <- rep(seq_len(n), each = length(times))
  betas <- c(-1, 0.5)
  D <- bdiag(rep(list(cbind(c(1, 0.2), c(0.2, 0.25))), Q))
  D[D == 0] <- c(0.1)
  D[3,1]=0.15
  D[5,1]=0.15
  D[7,1]=0.15
  D[5,3]=0.15
  D[7,3]=0.15
  D[7,5]=0.15
  
  D[1,3]=0.15
  D[1,5]=0.15
  D[1,7]=0.15
  D[3,5]=0.15
  D[3,7]=0.15
  D[5,7]=0.15
  
  
  X <- cbind(1, rep(times, n))
  Z <- cbind(1, rep(times, n))
  ncz <- ncol(Z)
  
  m <- 1
  set.seed(1000 + m)
  
  b <- mvrnorm(n, rep(0, Q*2), D)
  
  Data <- generateData2(id,times,n,X,Z,betas,b,Q)
  
  DataRaw <<- Data
  
  modelFit <- estimateModelFit2(Data,Q,n)
  
  betasN=c(-1,0.5,1,0.25,0.2,  -1,0.5,1,0.25,0.2, -1,0.5,1,0.25,0.2, -1,0.5,1,0.25,0.2)
  
  modelFitAndTarget <- as.data.frame(cbind("True" = betasN, modelFit))
  
  outD <- modelFitAndTarget 
  
  require(xtable)
  
  # LaTeX output
  #outD <- cbind("Param" = rep(paste("$\\theta_", 1:5, "$", sep = ""), len = nrow(outD)), outD)
  
  mychar <- character(5)
  mychar[1] <- "$\\beta_0$"
  mychar[2] <- "$\\beta_1$"
  mychar[3] <- "$\\sigma_{11}$"
  mychar[4] <- "$\\sigma_{22}$"
  mychar[5] <- "$\\sigma_{21}$"
  
  outD <- cbind("Param" = rep(paste(mychar[1:5],sep = " "),
                              len = nrow(outD)), outD)
  
  outD <- cbind(" " = c(sapply(1:Q, function (i) c(paste("Outcome", i), 
                                                   rep("", 5 - 1)))), outD)
  
  
   cap <- paste("Scenario ", Scenario, ", based on ", M, " datasets: ", 
                "Using ",length(times)," time points and ",Q," items: ",
                "Columns 'Val(ave)', 'Val(wave)', 'Val(dwave)', 'Val(ch-exp)' and 'Val(ch-ecdf)' give the average of the parameters,",
                "and columns SE(ave), SE(wave), SE(dwave), SE(ch-exp) and SE(ch-ecdf) give the average of the standard errors ",
                "over the simulated datasets based on ",
                "the simple, weighted average and diagonal weighted average estimators as well as the two clic heuristic estimators, respectively. ", sep = "")
                                                   
  print(xtable(outD, caption = cap,digits=c(3)), math.style.negative = TRUE, include.rownames = FALSE, 
        sanitize.text.function = function(x) x)
  
  
  return(modelFitAndTarget)
  
}