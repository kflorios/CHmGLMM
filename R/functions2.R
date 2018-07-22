gauher <- function (n) {
  EPS <- 3e-14
  PIM4 <- 0.751125544464943
  MAXIT <- 10
  m <- trunc((n + 1) / 2)
  x <- w <- rep(-1, n)
  for (i in 1:m) {
    if (i == 1) {
      z <- sqrt(2 * n + 1) - 1.85575 * (2 * n + 1)^(-0.16667)
    }
    else if (i == 2) {
      z <- z - 1.14 * n^0.426 / z
    }
    else if (i == 3) {
      z <- 1.86 * z - 0.86 * x[1]
    }
    else if (i == 4) {
      z <- 1.91 * z - 0.91 * x[2]
    }
    else {
      z <- 2 * z - x[i - 2]
    }
    for (its in 1:MAXIT) {
      p1 <- PIM4
      p2 <- 0
      for (j in 1:n) {
        p3 <- p2
        p2 <- p1
        p1 <- z * sqrt(2/j) * p2 - sqrt((j - 1)/j) * p3
      }
      pp <- sqrt(2 * n) * p2
      z1 <- z
      z <- z1 - p1/pp
      if (abs(z - z1) <= EPS) 
        break
    }
    x[i] <- z
    x[n + 1 - i] <- -z
    w[i] <- 2 / (pp * pp)
    w[n + 1 - i] <- w[i]
  }
  list(x = x, w = w)
}

dmvnorm <- function (x, mu, Sigma, log = FALSE) {
  if (!is.matrix(x))
    x <- rbind(x)
  p <- length(mu)
  ed <- eigen(Sigma, symmetric = TRUE)
  ev <- ed$values
  evec <- ed$vectors
  if (!all(ev >= -1e-06 * abs(ev[1]))) 
    stop("'Sigma' is not positive definite")
  ss <- x - rep(mu, each = nrow(x))
  inv.Sigma <- evec %*% (t(evec) / ev)
  quad <- 0.5 * rowSums((ss %*% inv.Sigma) * ss)
  fact <- - 0.5 * (p * log(2 * pi) + determinant(Sigma, logarithm = TRUE)$modulus)
  if (log)
    as.vector(fact - quad)
  else
    as.vector(exp(fact - quad))
}

logLik.bin <- function (thetas, id, y, X, Z, GHk = 5, extraParam) {
  #thetas <- relist(thetas, lis.thetas)
  #betas <- thetas$betas
  #ncz <- ncol(Z)
  #D <- matrix(0, ncz, ncz) 
  #D[lower.tri(D, TRUE)] <- thetas$D
  #D <- D + t(D)
  #diag(D) <- diag(D) / 2
  #
  betas<-thetas[1:4]
  Sigma2YiRI<-thetas[5]
  Sigma2YiRS<-thetas[6]
  Sigma2YiRIRS<-thetas[7]
  Sigma2YjRI<-thetas[8]
  Sigma2YjRS<-thetas[9]
  Sigma2YjRIRS<-thetas[10]
  
  Dold <- matrix(0,ncol=4,nrow=4)
  Dnew <- matrix(0,ncol=4,nrow=4)
  
  Dold <- extraParam   # extraParam is VarCorr(Models[[i]])$id
  #re-align columns, rows
  #old: glmer, 1,2,3,4
  #new: Florios, 1,3,2,4
  
  Dnew[1,1] = Dold[1,1]
  Dnew[2,2] = Dold[3,3]
  Dnew[3,3] = Dold[2,2]
  Dnew[4,4] = Dold[4,4]
  
  Dnew[2,1] = Dold[3,1]
  Dnew[3,2] = Dold[2,3]
  Dnew[4,3] = Dold[4,2]
  Dnew[3,1] = Dold[2,1]
  Dnew[4,2] = Dold[4,3]
  Dnew[4,1] = Dold[4,1]
  
  #now: fill in upper triangular: D new is symmetric
  for (ii in 1:4) {
    for (jj in ii:4) {
      Dnew[ii,jj] <- Dnew[jj,ii]
    }
  }
  
  ###D<- cbind(c(Sigma2YiRI,Sigma2YiRIRS,0,0),c(Sigma2YiRIRS,Sigma2YiRS,0,0),c(0,0,Sigma2YjRI,Sigma2YjRIRS),c(0,0,Sigma2YjRIRS,Sigma2YjRS))
  #instead: put full D matrix, as obtained from glmer() run
  
  #now plug-in the thetas, in order to perform numerical derivatives wrt theta (scores, hessians)
  Dnew[1,1] = Sigma2YiRI
  Dnew[2,2] = Sigma2YiRS
  Dnew[2,1] = Sigma2YiRIRS
  Dnew[1,2] = Dnew[2,1]
  
  Dnew[3,3] = Sigma2YjRI
  Dnew[4,4] = Sigma2YjRS
  Dnew[4,3] = Sigma2YjRIRS
  Dnew[3,4] = Dnew[4,3]
  
  #simplify notation: D=Dnew, and proceed
  
  D <-nearPD(Dnew)
  ncz <- ncol(Z)
  GH <- gauher(GHk)
  b <- as.matrix(expand.grid(rep(list(GH$x), ncz)))
  dimnames(b) <- NULL
  k <- nrow(b)
  wGH <- as.matrix(expand.grid(rep(list(GH$w), ncz)))
  wGH <- 2^(ncz/2) * apply(wGH, 1, prod) * exp(rowSums(b * b))
  b <- sqrt(2) * b
  Ztb <- Z %*% t(b)
  #
  mu.y <- plogis(as.vector(X %*% betas) + Ztb)
  logBinom <- dbinom(y, 1, mu.y, TRUE)
  log.p.yb <- rowsum(logBinom, id)
  log.p.b <- dmvnorm(b, rep(0, ncol(Z)), D, TRUE)
  p.yb <- exp(log.p.yb + rep(log.p.b, each = nrow(log.p.yb)))
  p.y <- c(p.yb %*% wGH)
  #-sum(log(p.y), na.rm = TRUE)  #logLik as min, original
  sum(log(p.y), na.rm = TRUE)  #logLik as max, CLIC heuristic
}

score.bin <- function (thetas, id, y, X, Z, GHk = 5, extraParam) {
  fd (thetas, logLik.bin, id = id, y = y, X = X, Z = Z, GHk = GHk, extraParam)
}

cd <- function (x, f, ..., eps = 0.001) {
  n <- length(x)
  res <- numeric(n)
  ex <- pmax(abs(x), 1)
  for (i in 1:n) {
    x1 <- x2 <- x
    x1[i] <- x[i] + eps * ex[i]
    x2[i] <- x[i] - eps * ex[i]
    diff.f <- c(f(x1, ...) - f(x2, ...))
    diff.x <- x1[i] - x2[i]
    res[i] <- diff.f/diff.x
  }
  res
}

fd <- function (x, f, ..., eps = 1e-04) {
  f0 <- f(x, ...)
  n <- length(x)
  res <- numeric(n)
  ex <- pmax(abs(x), 1)
  for (i in 1:n) {
    x1 <- x
    x1[i] <- x[i] + eps * ex[i]
    diff.f <- c(f(x1, ...) - f0)
    diff.x <- x1[i] - x[i]
    res[i] <- diff.f/diff.x
  }
  res
}

hess.bin <- function (thetas, id, y, X, Z, GHk = 5, extraParam) {
  #cdd(thetas, logLik.bin, id = id, y = y, X = X, Z = Z, GHk = GHk)
  res <- matrix(nrow=length(thetas),ncol=length(thetas))
  
  #for (i in 1:length(thetas)) {
  #for (j in 1:length(thetas)) {
  thetasIJ <- thetas
  res <- cdd(thetasIJ, logLik.bin, id = id, y = y, X = X, Z = Z, GHk = 5, extraParam)
  #}
  #}
  
  res
}


cdd <- function (x, f, ..., eps = 0.0005) {
  n <- length(x)
  res <- matrix(nrow=n,ncol=n)
  ex <- pmax(abs(x), 1)
  celem=0
  for (i in 1:n) {
    #for (j in 1:n){  hess: symmetric, economy: compute upper triangular part only
    for (j in i:n){
      celem=celem+1
      cat("Computing Element No...:", celem, "of hessian\n")
      x1 <- x2 <- x3 <- x4 <- x
      if (i != j)  { 
        # mixed 2nd order derivatives
        # where i--1 first move
        # where j--2 second move
        for (k in 1:n) {
          if ( i== k) {
            x1[k] = x[k] + eps * ex[i]
            x2[k] = x[k] + eps * ex[i]
            x3[k] = x[k] - eps * ex[i]
            x4[k] = x[k] - eps * ex[i]
          }
          if ( j== k) {
            x1[k] = x[k] + eps * ex[i]
            x2[k] = x[k] - eps * ex[i]
            x3[k] = x[k] + eps * ex[i]
            x4[k] = x[k] - eps * ex[i]
          }
          #x1 [1] = x[1] + eps * ex[i]
          #x1 [2] = x[2] + eps * ex[i]
          #x1 [3] = x[3]
          #x2 [1] = x[1] + eps * ex[i]
          #x2 [2] = x[2] - eps * ex[i]
          #x2 [3] = x[3]          
          #x3 [1] = x[1] - eps * ex[i]
          #x3 [2] = x[2] + eps * ex[i]
          #x3 [3] = x[3]          
          #x4 [1] = x[1] - eps * ex[i]
          #x4 [2] = x[2] - eps * ex[i]
          #x4 [3] = x[3]          
          
          diff.f <- c(f(x1, ...) -f(x2, ...) -f(x3, ...) +f(x4, ...))
          diff.x <- 2 * eps * ex[i]
          diff.y <- 2 * eps * ex[i]
          res[i,j] <- diff.f / (diff.x * diff.y)   # (1 -1 -1 + 1 / h^2)      
        }
        
      }
      if (i == j) { 
        # pure 2nd order derivatives
        # where i--1 first move
        # where j--2 second move
        for (k in 1:n) {
          if ( i== k) {
            x1[k] = x[k] + eps * ex[i]
            x2[k] = x[k]
            x3[k] = x[k] - eps * ex[i]
          }
        }
        #x1 [1] = x[1] + eps * ex[i]
        #x1 [2] = x[2] 
        #x1 [3] = x[3]          
        #x2 [1] = x[1] 
        #x2 [2] = x[2] 
        #x2 [3] = x[3]          
        #x3 [1] = x[1] - eps * ex[i]
        #x3 [2] = x[2] 
        #x3 [3] = x[3]                    
        
        diff.f <- c(f(x1, ...) -2*f(x2, ...) +f(x3, ...) )          
        diff.x <- eps * ex[i]
        res[i,i] <- diff.f / (diff.x * diff.x)   # (1 -2 + 1 / h^2)                
      }
    }
  }
  
  #now compute lower triangular also
  for (i in 1:n) {
    for (j in 1:i) {
      res[i,j] <- res[j,i]
    }
  }
  res # return Hessian in res matrix
  #as.vector(res) # return Hessian in res matrix as a vector by Rows
}


nearPD <- function (M, eig.tol = 1e-06, conv.tol = 1e-07, posd.tol = 1e-08, maxits = 100) {
  # based on function nearcor() submitted to R-help by Jens Oehlschlagel on 2007-07-13, and
  # function posdefify() from package `sfsmisc'
  if (!(is.numeric(M) && is.matrix(M) && identical(M, t(M))))
    stop("Input matrix M must be square and symmetric.\n")
  inorm <- function (x) max(rowSums(abs(x)))
  n <- ncol(M)
  U <- matrix(0, n, n)
  X <- M
  iter <- 0
  converged <- FALSE
  while (iter < maxits && !converged) {
    Y <- X
    T <- Y - U
    e <- eigen(Y, symmetric = TRUE)
    Q <- e$vectors
    d <- e$values
    D <- if (length(d) > 1) diag(d) else as.matrix(d)
    p <- (d > eig.tol * d[1])
    QQ <- Q[, p, drop = FALSE]
    X <- QQ %*% D[p, p, drop = FALSE] %*% t(QQ)
    U <- X - T
    X <- (X + t(X)) / 2
    conv <- inorm(Y - X) / inorm(Y)
    iter <- iter + 1
    converged <- conv <= conv.tol
  }
  X <- (X + t(X)) / 2
  e <- eigen(X, symmetric = TRUE)
  d <- e$values
  Eps <- posd.tol * abs(d[1])
  if (d[n] < Eps) {
    d[d < Eps] <- Eps
    Q <- e$vectors
    o.diag <- diag(X)
    X <- Q %*% (d * t(Q))
    D <- sqrt(pmax(Eps, o.diag) / diag(X))
    X[] <- D * X * rep(D, each = n)
  }
  (X + t(X)) / 2
}

bdiag <- function (...) {
  mlist<-list(...)
  ## handle case in which list of matrices is given
  if (length(mlist) == 1)
    mlist <- unlist(mlist, rec = FALSE)
  csdim <- rbind(c(0,0), apply(sapply(mlist,dim), 1, cumsum ))
  ret <- array(0,dim=csdim[length(mlist)+1,])
  add1 <- matrix(rep(1:0,2),nc=2)
  for(i in seq(along=mlist)){
    indx<-apply(csdim[i:(i+1),]+add1,2,function(x)x[1]:x[2])
    ## non-square matrix
    if(is.null(dim(indx)))ret[indx[[1]],indx[[2]]]<-mlist[[i]]
    ## square matrix
    else ret[indx[,1],indx[,2]]<-mlist[[i]]
  }
  ret
}

#' Computes the AVE, DWAVE, WAVE, CH-EXP and CH-ECDF methods for multivariate GLMMs (m-GLMMs)
#'
#' It acts on the Models structure for all pairs of items and returns the estimates for the m-GLMM parameters with 5 different methods
#'@param Models A list which contains the lme4 model objects taken from the pairwise separate estimations (list of size Q*(Q-1)/2)
#'@param ModelsOne A list which contains the lme4 model objects taken from the univariate separate estimations (list of size Q)
#'@param Data a data.frame with the data. 1st column id, 2nd column time, remaining Q columns are the y 0/1 values (Q items)
#'@param GHk Number of Gauss-Hermite points per dimension of integration
#'@param n   number of individuals
#'@param Q   the number of items. Set this to four.
#'@param extraParam  a helper list which depends on pairwise estimates
#'@param extraParamOne a helper list which depends on univariate estimates
#'@param m   an integer which is useful for subsequent runs
#'@return The estimated parameters of the model with methods AVE, DWAVE, WAVE, CH-EXP and CH-ECDF
#'@export
aveThetas2 <- function (Models,ModelsOne, Data, GHk = 5, n, Q, extraParam, extraParamOne, m=1) {
  # Compute K matrix
  #environment(logLik.bin) <- environment(score.bin) <- environment()
  #environment(logLik.bin.BIG) <- environment(score.bin.BIG) <- environment()
  #environment(logLik.bin.BIG) <- environment(score.bin.BIG) <- environment(estimateModelFit2)
  #environment(logLik.bin.BIG) <- environment(score.bin.BIG) <- environment(estimateModelFit2) <- environment()
  environment(logLik.bin.BIG) <- environment(hess.bin.BIG.Vect) <- environment(score.bin.BIG) <- 
    environment(logLik.bin)   <- environment(hess.bin.Vect)    <- environment(score.bin) <- 
    environment(estimateModelFit2) <- environment()
  
  
  DataRaw <- Data

  Sigma2YiRIone <- vector("list", Q)
  Sigma2YiRSone <- vector("list", Q)
  Sigma2YiRIRSone <- vector("list", Q)  
  for (ii in 1:Q) {
    Sigma2YiRIone[[ii]]   <- VarCorr(ModelsOne[[ii]])$id[1,1]
    Sigma2YiRSone[[ii]]   <- VarCorr(ModelsOne[[ii]])$id[2,2]
    Sigma2YiRIRSone[[ii]] <- VarCorr(ModelsOne[[ii]])$id[1,2]      
  }  
  
  Klis <- vector("list", n)
  for (i in 1:n) {
    cat("Individual No...:", i, "for scores computation\n")
    pairs <- combn(Q, 2)
    P <- ncol(pairs)  
    Scores    <- vector("list", P)
    Sigma2YiRI   <- vector("list", P)
    Sigma2YiRS   <- vector("list", P)    
    Sigma2YiRIRS <- vector("list", P)
    Sigma2YjRI <- vector("list", P)    
    Sigma2YjRS <- vector("list", P)
    Sigma2YjRIRS <- vector("list", P)

    
    
    for (p in 1:P) {
      cp <-p
      Sigma2YiRI[[cp]]   <-   VarCorr(Models[[p]])$id[1,1]
      Sigma2YiRS[[cp]]   <-   VarCorr(Models[[p]])$id[3,3]
      Sigma2YiRIRS[[cp]] <-   VarCorr(Models[[p]])$id[3,1]
      Sigma2YjRI[[cp]]   <-   VarCorr(Models[[p]])$id[2,2]
      Sigma2YjRS[[cp]]   <-   VarCorr(Models[[p]])$id[4,4]
      Sigma2YjRIRS[[cp]] <-   VarCorr(Models[[p]])$id[4,2]  
      prs <- pairs[, p]
      yi <- Data[[paste("y", prs[1], sep = "")]]
      yj <- Data[[paste("y", prs[2], sep = "")]]
      DD <- do.call(rbind, list(Data[1:3], Data[1:3])) 
      DD$outcome <- gl(2, nrow(Data))
      DD$y <- c(yi, yj)
      DD.i <- DD[(ind.i <- DD$id == i), ]
      D <- nearPD(VarCorr(Models[[p]])$id)
      indices <- 1:Q
      indic   <- indices[c(-prs[1],-prs[2])]        
      #Scores[[p]] <- score.bin(fixef(Models[[p]]), id = DD.i$id, y = DD.i$y,       
      #                         X = model.matrix(Models[[p]])[ind.i, ], 
      #                         Z = model.matrix(Models[[p]])[ind.i, seq_len(2*ncz)], 
      #                         GHk = GHk)            
      ###Scores[[p]] <- score.bin(c(fixef(Models[[p]]),SigmaYi[[p]],SigmaYj[[p]],rhoYiYj[[p]]), id = DD.i$id, y = DD.i$y, 
      ###                         X = model.matrix(Models[[p]])[ind.i, ], 
      ###                         Z = cbind(model.matrix(Models[[p]])[ind.i, seq_len(2*ncz)]), 
      ###                         GHk = GHk)            
      
      Scores[[p]] <- #score.bin(c(fixef(Models[[p]]),Sigma2YiRI[[p]],Sigma2YiRS[[p]],Sigma2YiRIRS[[p]],Sigma2YjRI[[p]],Sigma2YjRS[[p]],Sigma2YjRIRS[[p]]),
        #score.bin.BIG(c(fixef(Models[[p]]),Sigma2YiRI[[p]],Sigma2YiRS[[p]],Sigma2YiRIRS[[p]],Sigma2YjRI[[p]],Sigma2YjRS[[p]],Sigma2YjRIRS[[p]]),
        ##score.bin.BIG(c(fixef(Models[[p]]),Sigma2YiRI[[p]],Sigma2YiRS[[p]],Sigma2YiRIRS[[p]],Sigma2YjRI[[p]],Sigma2YjRS[[p]],Sigma2YjRIRS[[p]]),                                   
        score.bin.BIG(c(fixef(Models[[p]]),Sigma2YiRI[[p]],Sigma2YiRS[[p]],Sigma2YiRIRS[[p]],Sigma2YjRI[[p]],Sigma2YjRS[[p]],Sigma2YjRIRS[[p]],
                        fixef(ModelsOne[[indic[1]]]),Sigma2YiRIone[[indic[1]]],Sigma2YiRSone[[indic[1]]],Sigma2YiRIRSone[[indic[1]]],
                        fixef(ModelsOne[[indic[2]]]),Sigma2YiRIone[[indic[2]]],Sigma2YiRSone[[indic[2]]],Sigma2YiRIRSone[[indic[2]]]),                                                                      
                      id = DD.i$id, y = DD.i$y, 
                      X = model.matrix(Models[[p]])[ind.i, ], 
                      Z = cbind(model.matrix(Models[[p]])[ind.i, ]), 
                      GHk = GHk,
                      extraParam = extraParam[[p]],
                      Data = DD.i,
                      r = p,
                      i = i)            
      
      #thetas <- c(fixef(Models[[p]]),SigmaYi[[p]],SigmaYj[[p]],rhoYiYj[[p]])
      #ii<- prs[1]
      #jj<- prs[2]
      #ic[ii]  <- thetas[1]
      #ic[jj]  <- thetas[2]
      #sl[ii]  <- thetas[3]
      #sl[jj]  <- thetas[4]
      #sig[ii] <- thetas[5]
      #sig[jj] <- thetas[6]
      #rh[p]   <- thetas[7]
      
      
    }
    #kkk <- P * length( c( fixef(Models[[1]]) , SigmaYi[[1]], SigmaYj[[1]], rhoYiYj[[1]] )) 
    kkk <- P * length( c( fixef(Models[[1]]) , Sigma2YiRI[[1]], Sigma2YiRS[[1]], Sigma2YiRIRS[[1]], Sigma2YjRI[[1]], Sigma2YjRS[[1]], Sigma2YjRIRS[[1]],
                          rep(c(fixef(ModelsOne[[1]]), Sigma2YiRIone[[1]], Sigma2YiRSone[[1]], Sigma2YiRIRSone[[1]]),2))) 
    K <- matrix(0, kkk, kkk)
    ee <- expand.grid(1:P, 1:P)
    ss <- sapply(Scores, length)
    ss2 <- cumsum(ss)
    ss1 <- c(1, ss2[-P] + 1)
    for (ii in 1:nrow(ee)) {
      k <- ee$Var1[ii]
      j <- ee$Var2[ii]
      row.ind <- seq(ss1[k], ss2[k])
      col.ind <- seq(ss1[j], ss2[j])
      K[row.ind, col.ind] <- Scores[[k]] %o% Scores[[j]]
      K[row.ind, col.ind] <- Scores[[k]] %o% Scores[[j]]
    }
    Klis[[i]] <- K
  }
  K <- Reduce("+", Klis)
  # Extract J matrix
  
  ###J <- lapply(Models, vcov)
  ##J<-matrix(0,18,18)
  Hessians<- vector("list", P) 
  for (p in 1:P) {
    cat("Hessian No: ...",p,"\n")
    prs <- pairs[, p]
    yi <- Data[[paste("y", prs[1], sep = "")]]
    yj <- Data[[paste("y", prs[2], sep = "")]]
    DD <- do.call(rbind, list(Data[1:3], Data[1:3])) 
    DD$outcome <- gl(2, nrow(Data))
    DD$y <- c(yi, yj)
    #DD.i <- DD[(ind.i <- DD$id == i), ]  #take all observations now, from all i=1:n
    D <- nearPD(VarCorr(Models[[p]])$id)
    indices <- 1:Q
    indic   <- indices[c(-prs[1],-prs[2])]        
    
    #Hessians[[p]] <- hess.bin(c( fixef(Models[[p]]),Sigma2YiRI[[p]],Sigma2YiRS[[p]],Sigma2YiRIRS[[p]],Sigma2YjRI[[p]],Sigma2YjRS[[p]],Sigma2YjRIRS[[p]] ),
    #                          id = DD$id, y = DD$y, 
    #                          X = model.matrix(Models[[p]])[, ], 
    #                          Z = cbind(model.matrix(Models[[p]])[, seq_len(2*ncz)]), 
    #                          GHk = GHk,
    #                          extraParam = extraParam[[p]])  
    Hessians[[p]] <- #hess.bin.Vect(c( fixef(Models[[p]]),Sigma2YiRI[[p]],Sigma2YiRS[[p]],Sigma2YiRIRS[[p]],Sigma2YjRI[[p]],Sigma2YjRS[[p]],Sigma2YjRIRS[[p]] ),
      hess.bin.BIG.Vect(c( fixef(Models[[p]]),Sigma2YiRI[[p]],Sigma2YiRS[[p]],Sigma2YiRIRS[[p]],Sigma2YjRI[[p]],Sigma2YjRS[[p]],Sigma2YjRIRS[[p]],
                           fixef(ModelsOne[[indic[1]]]),Sigma2YiRIone[[indic[1]]],Sigma2YiRSone[[indic[1]]],Sigma2YiRIRSone[[indic[1]]],
                           fixef(ModelsOne[[indic[2]]]),Sigma2YiRIone[[indic[2]]],Sigma2YiRSone[[indic[2]]],Sigma2YiRIRSone[[indic[2]]]),               
                        id = DD$id, y = DD$y,
                        X = model.matrix(Models[[p]])[, ],
                        Z = cbind(model.matrix(Models[[p]])[, ]),
                        GHk = GHk,
                        extraParam = extraParam[[p]],
                        Data = DD,
                        r = p,
                        i = n+1)  
    
    #thetas <- c(fixef(Models[[p]]),SigmaYi[[p]],SigmaYj[[p]],rhoYiYj[[p]])
    #ii<- prs[1]
    #jj<- prs[2]
    #ic[ii]  <- thetas[1]
    #ic[jj]  <- thetas[2]
    #sl[ii]  <- thetas[3]
    #sl[jj]  <- thetas[4]
    #sig[ii] <- thetas[5]
    #sig[jj] <- thetas[6]
    #rh[p]   <- thetas[7]
    
  }
  
  J <- do.call(bdiag, Hessians)
  
  
  ## Now J is redefined as the block matrix 108x108, or (6x18) x (6x18)
  
  ###J <- bdiag(lapply(J, function (x) solve(matrix(x@x, x@Dim[1], x@Dim[2]))))
  
  # Compute A matrix
  nbetas <- length( c(fixef(Models[[1]]) , Sigma2YiRI[[1]], Sigma2YiRS[[1]], Sigma2YiRIRS[[1]], Sigma2YjRI[[1]], Sigma2YjRS[[1]], Sigma2YjRIRS[[1]],
                      rep(c(fixef(ModelsOne[[1]]), Sigma2YiRIone[[1]], Sigma2YiRSone[[1]], Sigma2YiRIRSone[[1]]),2) ) )
  #thetas <- c( sapply(Models, fixef) )
  thetas <- numeric(0)
  for (ii in 1:20*P) {
    thetas[ii]=0
  }
  for (ii in 1:P) {
    thetas[(ii-1)*20+1]=fixef(Models[[ii]])[1]
    thetas[(ii-1)*20+2]=fixef(Models[[ii]])[2]     
    thetas[(ii-1)*20+3]=fixef(Models[[ii]])[3]     
    thetas[(ii-1)*20+4]=fixef(Models[[ii]])[4]             
    thetas[(ii-1)*20+5]=Sigma2YiRI[[ii]]   
    thetas[(ii-1)*20+6]=Sigma2YiRS[[ii]]   
    thetas[(ii-1)*20+7]=Sigma2YiRIRS[[ii]]   
    thetas[(ii-1)*20+8]=Sigma2YjRI[[ii]]         
    thetas[(ii-1)*20+9]=Sigma2YjRS[[ii]]         
    thetas[(ii-1)*20+10]=Sigma2YjRIRS[[ii]]               
    p<-ii
    prs <- pairs[, p]    
    indices <- 1:Q
    indic   <- indices[c(-prs[1],-prs[2])]        
    thetas[(ii-1)*20+11]=fixef(ModelsOne[[indic[1]]])[1]
    thetas[(ii-1)*20+12]=fixef(ModelsOne[[indic[1]]])[2]
    thetas[(ii-1)*20+13]=Sigma2YiRIone[[indic[1]]]
    thetas[(ii-1)*20+14]=Sigma2YiRSone[[indic[1]]]
    thetas[(ii-1)*20+15]=Sigma2YiRIRSone[[indic[1]]]
    
    thetas[(ii-1)*20+16]=fixef(ModelsOne[[indic[2]]])[1]
    thetas[(ii-1)*20+17]=fixef(ModelsOne[[indic[2]]])[2]
    thetas[(ii-1)*20+18]=Sigma2YiRIone[[indic[2]]]
    thetas[(ii-1)*20+19]=Sigma2YiRSone[[indic[2]]]
    thetas[(ii-1)*20+20]=Sigma2YiRIRSone[[indic[2]]]
  }
  
  Kone <-    vector("list", P) 
  HelpMat <- vector("list", P) 
  CLIC    <- vector("list", P) 
  wtCLIC     <- vector("list", P) 
  wtCLIC_B    <- vector("list", P) 
  CLIC_stdz    <- vector("list", P) 
  
  for (p in 1:6) {
    Kone[[p]] <- K[((p-1)*nbetas+1):((p-1)*nbetas+nbetas),((p-1)*nbetas+1):((p-1)*nbetas+nbetas)]
  }
  
  for (p in 1:6) {
    HelpMat[[p]] <- Kone[[p]] %*% solve(Hessians[[p]])
  }
  
  for (p in 1:6) {
    #CLIC[[p]] <- logLik(Models[[p]]) + sum(diag(HelpMat[[p]]))   #sum(diag(A)) is trace(A) in R language for A square matrix
    
    prs <- pairs[, p]
    yi <- Data[[paste("y", prs[1], sep = "")]]
    yj <- Data[[paste("y", prs[2], sep = "")]]
    DD$y <- c(yi, yj)
    #DD.i <- DD[(ind.i <- DD$id == i), ]  #take all observations now, from all i=1:n
    D <- nearPD(VarCorr(Models[[p]])$id)
    indices <- 1:Q
    indic   <- indices[c(-prs[1],-prs[2])]        
    
    CLIC[[p]] <- logLik.bin.BIG(c( fixef(Models[[p]]),Sigma2YiRI[[p]],Sigma2YiRS[[p]],Sigma2YiRIRS[[p]],Sigma2YjRI[[p]],Sigma2YjRS[[p]],Sigma2YjRIRS[[p]],
                                   fixef(ModelsOne[[indic[1]]]),Sigma2YiRIone[[indic[1]]],Sigma2YiRSone[[indic[1]]],Sigma2YiRIRSone[[indic[1]]],               
                                   fixef(ModelsOne[[indic[2]]]),Sigma2YiRIone[[indic[2]]],Sigma2YiRSone[[indic[2]]],Sigma2YiRIRSone[[indic[2]]]),               
                                id = DD$id, y = DD$y,
                                X = model.matrix(Models[[p]])[, ],
                                Z = cbind(model.matrix(Models[[p]])[, ]),
                                GHk = GHk,
                                extraParam = extraParam[[p]],
                                Data = DD,
                                r = p,
                                i = n+1) + sum(diag(HelpMat[[p]]))   #sum(diag(A)) is trace(A) in R language for A square matrix
  }
  
  #logLik_Models <- c(logLik(Models[[1]]),logLik(Models[[2]]),logLik(Models[[3]]),
  #                   logLik(Models[[4]]),logLik(Models[[5]]),logLik(Models[[6]]))
  
  trace_Models <- c(sum(diag(HelpMat[[1]])),sum(diag(HelpMat[[2]])),sum(diag(HelpMat[[3]])),
                    sum(diag(HelpMat[[4]])),sum(diag(HelpMat[[5]])),sum(diag(HelpMat[[6]])))
  
  CLIC_Models <- c(CLIC[[1]],CLIC[[2]],CLIC[[3]],CLIC[[4]],CLIC[[5]],CLIC[[6]])
  
  ## 7.9.2015: KJF. Standardize the CLIC values to mean 0 and sd 1
  
  mu_CLIC <-mean(CLIC_Models)
  sigma_CLIC <- sd(CLIC_Models)
  
  CLIC_Models_stdz <- (CLIC_Models - mu_CLIC) / sigma_CLIC
  
  logLik_Models <- CLIC_Models - trace_Models  #simplify, not call again logLik.bin
  
  write.table(logLik_Models,paste(m,"_logLik_Models.txt",sep=""),row.names=F,col.names=F)
  write.table(trace_Models,paste(m,"_trace_Models.txt",sep=""),row.names=F,col.names=F)
  write.table(CLIC_Models,paste(m,"_CLIC_Models.txt",sep=""),row.names=F,col.names=F)
  
  #weights <- exp(CLIC)
  #A2 <- weights/sum(weights) 
  
  #7.9.2015. KJF. standardization in CLIC values for weights
  for (p in 1:6) {
    CLIC_stdz[[p]] <- CLIC_Models_stdz[p]
  }
  for (p in 1:6) {
    #wtCLIC[[p]]    <- exp(CLIC[[p]])    #VV idea
    ##wtCLIC[[p]]    <- exp(mpfr(CLIC[[p]],80))    #KF idea
    wtCLIC[[p]]    <- exp(CLIC_stdz[[p]])    #KF idea, standardization 7.9.2015
    #wtCLIC[[p]]    <- as.numeric(exp(as.brob(CLIC[[p]])))     #KF idea
    ###wtCLIC[[p]] <- - 1 / CLIC[[p]]    #KF idea
  }
  
  # 4.10.2015, KJF add.
  # CLIC ecdf way: B' way
  a=CLIC_Models_stdz
  f=ecdf(a)
  wtCLIC_B[[1]] = f(a)[1]
  wtCLIC_B[[2]] = f(a)[2]
  wtCLIC_B[[3]] = f(a)[3]
  wtCLIC_B[[4]] = f(a)[4]
  wtCLIC_B[[5]] = f(a)[5]
  wtCLIC_B[[6]] = f(a)[6]
  # 4.10.2015, end KJF add.
  
  OmegaVec   <- wtCLIC
  OmegaVec_B <- wtCLIC_B
  
  A<- computeWeightMatrixAVE_by6()   # to re-write for 20x120 case KF 25.6.2015, todo:
  
  #Omega <- matrix(0,P*nbetas,P*nbetas)
  
  #A2 <- matrix(0, Q*(nbetas)/2, length(thetas))
  A2   <- matrix(0, 20, 120)
  A2_B <- matrix(0, 20, 120)
  ##A2 <- new("mpfrMatrix", mpfr(rep(0,20*120),80),Dim = c(20L, 120L))
  ##validObject(A2)
  
  for (i in seq_len(Q*nbetas/4)) {
    #ii <- inter == levels(inter)[i]  just take all 6 which are intcpt, so odd
    ###ii <- rep(c(TRUE,FALSE,FALSE),choose(length(times),2))
    ii <- A[i,]==(1/6)  #adjacency matrix for global parameter i according to A = computeWeightMatrixAVE_by6()
    
    i2 <- which(ii==T)
    
    jj <- c(intervalBIG(i2[1]),intervalBIG(i2[2]),intervalBIG(i2[3]),intervalBIG(i2[4]),intervalBIG(i2[5]),intervalBIG(i2[6]))
    
    #weights <- OmegaVec[jj]
    
    weights   <- c(OmegaVec[[jj[1]]],  OmegaVec[[jj[2]]],  OmegaVec[[jj[3]]],  OmegaVec[[jj[4]]],  OmegaVec[[jj[5]]],  OmegaVec[[jj[6]]])
    weights_B <- c(OmegaVec_B[[jj[1]]],OmegaVec_B[[jj[2]]],OmegaVec_B[[jj[3]]],OmegaVec_B[[jj[4]]],OmegaVec_B[[jj[5]]],OmegaVec_B[[jj[6]]])
    
    A2[i, ii]   <- weights / sum(weights)
    A2_B[i, ii] <- weights_B / sum(weights_B)
    
  }  
  
  #ind.outcome <- c(apply(pairs, 2, rep, length.out = nbetas))
  #ind.param <- rep(rep(seq_len(nbetas/2), each = 2), length.out = length(thetas))
  #inter <- interaction(ind.outcome, ind.param)
  #inter <- factor(inter, levels = sort(levels(inter)))
  #A <- matrix(0, Q*nbetas/2, length(thetas))
  #for (i in seq_len(Q*nbetas/2)) {
  #    A[i, inter == levels(inter)[i]] <- 1 / sum(pairs == 1)
  #}
  
  #Compute A to be 16 x 48 ad hoc (by pencil)
  
  A<- computeWeightMatrixAVE_by6()    
  write.table(A,paste(m,"_A.txt",sep=""),row.names=F,col.names=F)
  
  #A2 <- asNumeric(A2)
  
  #A2 <- matrix(A2,nrow=20,ncol=120)
  write.table(A2,paste(m,"_A2.txt",sep=""),row.names=F,col.names=F)
  write.table(A2_B,paste(m,"_A2_B.txt",sep=""),row.names=F,col.names=F)
  
  
  # pairwise average betas
  ##ave.betas <- c(A %*% thetas)
  # standrd errors for betas
  ##se.betas <- sqrt(diag(A %*% solve(J, K) %*% solve(J) %*% t(A)))
  # CLIC heuristic
  ave.betas2 <- c(A2 %*% thetas)
  # standrd errors for betas
  se.betas2 <- sqrt(diag(A2 %*% solve(J, K) %*% solve(J) %*% t(A2)))
  # 4.10.2015. KJF add, B' way for CLIC heuristic with ecdf()
  # CLIC heuristic
  ave.betas5 <- c(A2_B %*% thetas)
  # standrd errors for betas
  se.betas5 <- sqrt(diag(A2_B %*% solve(J, K) %*% solve(J) %*% t(A2_B)))
  # 4.10.2015. end KJF add, B' way for CLIC heuristic with ecdf()  
  
  #now regular methods
  #JJ,KK
  #KK
  
  KKlis <- vector("list", n)
  for (i in 1:n) {
    cat("Individual No...:", i, "for sscores computation\n")
    SScores    <- vector("list", P)
    ic <-  vector("list",Q)
    sl <-  vector("list",Q)
    sigz <- vector("list",Q)
    sigw <- vector("list",Q)
    ###rh <-  vector("list",P)
    
    for (p in 1:P) {
      prs <- pairs[, p]
      yi <- Data[[paste("y", prs[1], sep = "")]]
      yj <- Data[[paste("y", prs[2], sep = "")]]
      DD$y <- c(yi, yj)
      DD.i <- DD[(ind.i <- DD$id == i), ]
      D <- nearPD(VarCorr(Models[[p]])$id)
      #Scores[[p]] <- score.bin(fixef(Models[[p]]), id = DD.i$id, y = DD.i$y, 
      #                         X = model.matrix(Models[[p]])[ind.i, ], 
      #                         Z = model.matrix(Models[[p]])[ind.i, seq_len(2*ncz)], 
      #                         GHk = GHk)            
      ###Scores[[p]] <- score.bin(c(fixef(Models[[p]]),SigmaYi[[p]],SigmaYj[[p]],rhoYiYj[[p]]), id = DD.i$id, y = DD.i$y, 
      ###                         X = model.matrix(Models[[p]])[ind.i, ], 
      ###                         Z = cbind(model.matrix(Models[[p]])[ind.i, seq_len(2*ncz)]), 
      ###                         GHk = GHk)            
      SScores[[p]] <- score.bin(c(fixef(Models[[p]]),Sigma2YiRI[[p]],Sigma2YiRS[[p]],Sigma2YiRIRS[[p]],Sigma2YjRI[[p]],Sigma2YjRS[[p]],Sigma2YjRIRS[[p]]),
                                #score.bin.BIG(c(fixef(Models[[p]]),Sigma2YiRI[[p]],Sigma2YiRS[[p]],Sigma2YiRIRS[[p]],Sigma2YjRI[[p]],Sigma2YjRS[[p]],Sigma2YjRIRS[[p]]),
                                #score.bin.BIG(c(fixef(Models[[p]]),Sigma2YiRI[[p]],Sigma2YiRS[[p]],Sigma2YiRIRS[[p]],Sigma2YjRI[[p]],Sigma2YjRS[[p]],Sigma2YjRIRS[[p]]),                                   
                                id = DD.i$id, y = DD.i$y, 
                                X = model.matrix(Models[[p]])[ind.i, ], 
                                Z = cbind(model.matrix(Models[[p]])[ind.i, ]), 
                                GHk = GHk,
                                extraParam = extraParam[[p]])            
      
      #thetas <- c(fixef(Models[[p]]),SigmaYi[[p]],SigmaYj[[p]],rhoYiYj[[p]])
      #ii<- prs[1]
      #jj<- prs[2]
      #ic[ii]  <- thetas[1]
      #ic[jj]  <- thetas[2]
      #sl[ii]  <- thetas[3]
      #sl[jj]  <- thetas[4]
      #sig[ii] <- thetas[5]
      #sig[jj] <- thetas[6]
      #rh[p]   <- thetas[7]
      
      
    }
    #kkk <- P * length( c( fixef(Models[[1]]) , SigmaYi[[1]], SigmaYj[[1]], rhoYiYj[[1]] )) 
    kkkk <- P * length( c( fixef(Models[[1]]) , Sigma2YiRI[[1]], Sigma2YiRS[[1]], Sigma2YiRIRS[[1]], Sigma2YjRI[[1]], Sigma2YjRS[[1]], Sigma2YjRIRS[[1]]  )) 
    KK <- matrix(0, kkkk, kkkk)
    ee <- expand.grid(1:P, 1:P)
    ss <- sapply(SScores, length)
    ss2 <- cumsum(ss)
    ss1 <- c(1, ss2[-P] + 1)
    for (ii in 1:nrow(ee)) {
      k <- ee$Var1[ii]
      j <- ee$Var2[ii]
      row.ind <- seq(ss1[k], ss2[k])
      col.ind <- seq(ss1[j], ss2[j])
      KK[row.ind, col.ind] <- SScores[[k]] %o% SScores[[j]]
      KK[row.ind, col.ind] <- SScores[[k]] %o% SScores[[j]]
    }
    KKlis[[i]] <- KK
  }
  KK <- Reduce("+", KKlis)
  # Extract J matrix
  HHessians<- vector("list", P) 
  for (p in 1:P) {
    cat("HHessian No: ...",p,"\n")
    prs <- pairs[, p]
    yi <- Data[[paste("y", prs[1], sep = "")]]
    yj <- Data[[paste("y", prs[2], sep = "")]]
    DD$y <- c(yi, yj)
    #DD.i <- DD[(ind.i <- DD$id == i), ]  #take all observations now, from all i=1:n
    D <- nearPD(VarCorr(Models[[p]])$id)
    
    #Hessians[[p]] <- hess.bin(c( fixef(Models[[p]]),Sigma2YiRI[[p]],Sigma2YiRS[[p]],Sigma2YiRIRS[[p]],Sigma2YjRI[[p]],Sigma2YjRS[[p]],Sigma2YjRIRS[[p]] ),
    #                          id = DD$id, y = DD$y, 
    #                          X = model.matrix(Models[[p]])[, ], 
    #                          Z = cbind(model.matrix(Models[[p]])[, seq_len(2*ncz)]), 
    #                          GHk = GHk,
    #                          extraParam = extraParam[[p]])  
    HHessians[[p]] <- hess.bin.Vect(c( fixef(Models[[p]]),Sigma2YiRI[[p]],Sigma2YiRS[[p]],Sigma2YiRIRS[[p]],Sigma2YjRI[[p]],Sigma2YjRS[[p]],Sigma2YjRIRS[[p]] ),
                                    #hess.bin.BIG.Vect(c( fixef(Models[[p]]),Sigma2YiRI[[p]],Sigma2YiRS[[p]],Sigma2YiRIRS[[p]],Sigma2YjRI[[p]],Sigma2YjRS[[p]],Sigma2YjRIRS[[p]] ),               
                                    id = DD$id, y = DD$y,
                                    X = model.matrix(Models[[p]])[, ],
                                    Z = cbind(model.matrix(Models[[p]])[, ]),
                                    GHk = GHk,
                                    extraParam = extraParam[[p]])  
    
  }
  
  JJ <- do.call(bdiag, HHessians)
  
  # Compute A matrix
  nbetasS <- length( c(fixef(Models[[1]]) , Sigma2YiRI[[1]], Sigma2YiRS[[1]], Sigma2YiRIRS[[1]], Sigma2YjRI[[1]], Sigma2YjRS[[1]], Sigma2YjRIRS[[1]]) )
  #thetas <- c( sapply(Models, fixef) )
  thetasS <- numeric(0)
  for (ii in 1:10*P) {
    thetasS[ii]=0
  }
  for (ii in 1:P) {
    thetasS[(ii-1)*10+1]=fixef(Models[[ii]])[1]
    thetasS[(ii-1)*10+2]=fixef(Models[[ii]])[2]     
    thetasS[(ii-1)*10+3]=fixef(Models[[ii]])[3]     
    thetasS[(ii-1)*10+4]=fixef(Models[[ii]])[4]             
    thetasS[(ii-1)*10+5]=Sigma2YiRI[[ii]]   
    thetasS[(ii-1)*10+6]=Sigma2YiRS[[ii]]   
    thetasS[(ii-1)*10+7]=Sigma2YiRIRS[[ii]]   
    thetasS[(ii-1)*10+8]=Sigma2YjRI[[ii]]         
    thetasS[(ii-1)*10+9]=Sigma2YjRS[[ii]]         
    thetasS[(ii-1)*10+10]=Sigma2YjRIRS[[ii]]               
    
  }
  
  
  AA <- computeWeightMatrixAVE()
  #Compute JJ,KK,thetasS,AA READY
  # pairwise average betas
  ave.betas <- c(AA %*% thetasS)
  # standrd errors for betas
  se.betas <- sqrt(diag(AA %*% solve(JJ, KK) %*% solve(JJ) %*% t(AA)))
  
  # DWAVE algorithm is method no. 3:
  #Step 1: Vi computation
  Sigma <- solve(JJ, KK) %*% solve(JJ)
  Vi    <- Sigma
  #Step 2: Omega computation with filter (mask)
  Omega<-matrix(0,P*nbetasS,P*nbetasS)
  crow=0
  ccol=0
  for (ii in 1:P) {
    crow= (ii-1)*nbetasS
    for (jj in 1:P) {
      ccol = (jj-1)*nbetasS
      for (kk in 1:nbetasS) {
        Omega[crow+kk, ccol+kk] <- Vi[crow+kk, ccol+kk]
      }
    }
  }
  #Step 3: Omega inversion in new Omega
  if ( abs(det(Omega)) > 10^(-10) ) 
  {
    Omega <- solve(Omega)
  }
  if ( abs(det(Omega)) < 10^(-10) ) 
  {
    Omega <- ginv(Omega)
  }
  #Step 4: Formulas for A
  A<- matrix(0,nbetasS,P*nbetasS)
  crow=0
  ccol=0
  crow2=0
  ccol2=0
  for (kk in 1:nbetasS) {
    
    #compute denom, which does not depend on r
    denom <-0
    
    for (ii in 1:P) {
      for (jj in 1:P) {
        crow2 <- (ii-1)*nbetasS
        ccol2 <- (jj-1)*nbetasS
        denom <- denom + Omega[crow2+kk,ccol2+kk]
      }
    }
    
    #compute nom, which depends on r
    for (r in 1:P) {
      nom <-0
      for (l in 1:P) {
        crow <- (r-1)*nbetasS
        ccol <- (l-1)*nbetasS
        
        nom  <- nom + Omega[crow+kk,ccol+kk]
      }
      # nom is ready      
      # for everyy kk (outer loop) and r (inner loop) compute nom/denom       
      A[kk,(r-1)*nbetasS+kk] <- nom/denom
    }
    
  }
  # pairwise average betas
  
  A5<-A  #for now
  
  
  A3 <- computeWeightMatrixWAVE(A5)
  
  ave.betas3 <- c(A3 %*% thetasS)
  # standrd errors for betas
  se.betas3 <- sqrt(diag(A3 %*% solve(JJ, KK) %*% solve(JJ) %*% t(A3)))    
  
  
  # WAVE algorithm is method no. 4:
  ###Omega <- solve(Sigma * JJ)
  if (det(Sigma) > 10^(-10)) {
    #Sigma <- solve(J, K) %*% solve(J)
    Omega <- solve(Sigma * JJ)
  }
  if (det(Sigma) < 10^(-10)) {
    #  Sigma <- solve(J, K) %*% solve(J)
    Omega <- ginv(Sigma * JJ)   # use ginv() function from MASS, generalized inverse, seems robust to det -->0
  }
  npairs <- sum(pairs == 1)
  A4 <- matrix(0, Q*nbetasS/2, length(thetasS))
  for (i in seq_len(Q*nbetasS/2)) {
    #ii <- inter == levels(inter)[i]
    ii <- AA[i,]==(1/3)  #adjacency matrix for global parameter i according to AA = computeWeightMatrixAVE()    
    weights <- diag(Omega[ii, ii])
    A4[i, ii] <- weights / sum(weights)
  }
  ave.betas4 <- c(A4 %*% thetasS)
  # standrd errors for betas
  ###se.betas2 <- sqrt(diag(A2 %*% solve(J, K) %*% solve(J) %*% t(A2)))
  ###se.betas3 <- sqrt(diag(A3 %*% solve(J, K) %*% solve(J) %*% t(A3)))    
  se.betas4 <- sqrt(diag(A4 %*% solve(JJ, KK) %*% solve(JJ) %*% t(A4)))    
  
  write.table(AA,paste(m,"_AA.txt",sep=""),row.names=F,col.names=F)
  write.table(A2,paste(m,"_A2.txt",sep=""),row.names=F,col.names=F)
  write.table(A3,paste(m,"_A3.txt",sep=""),row.names=F,col.names=F)              
  write.table(A4,paste(m,"_A4.txt",sep=""),row.names=F,col.names=F)              
  write.table(A2_B,paste(m,"_A2_B.txt",sep=""),row.names=F,col.names=F)                
  
  write.table(ave.betas,paste(m,"_ave.betas.txt",sep=""),row.names=F,col.names=F)
  write.table(ave.betas2,paste(m,"_ave.betas2.txt",sep=""),row.names=F,col.names=F)
  write.table(ave.betas3,paste(m,"_ave.betas3.txt",sep=""),row.names=F,col.names=F)
  write.table(ave.betas4,paste(m,"_ave.betas4.txt",sep=""),row.names=F,col.names=F)
  write.table(ave.betas5,paste(m,"_ave.betas5.txt",sep=""),row.names=F,col.names=F)
  
  write.table(se.betas,paste(m,"_se.betas.txt",sep=""),row.names=F,col.names=F)
  write.table(se.betas2,paste(m,"_se.betas2.txt",sep=""),row.names=F,col.names=F)
  write.table(se.betas3,paste(m,"_se.betas3.txt",sep=""),row.names=F,col.names=F)
  write.table(se.betas4,paste(m,"_se.betas4.txt",sep=""),row.names=F,col.names=F)
  write.table(se.betas5,paste(m,"_se.betas5.txt",sep=""),row.names=F,col.names=F)
  
  # results
  cbind("Value(ave)" = ave.betas, "SE(ave)" = se.betas, 
        "Value(clic-h-exp)" = ave.betas2, "SE(clic-h-exp)" = se.betas2,
        "Value(dwave)" = ave.betas3, "SE(dwave)" = se.betas3,
        "Value(wave)" = ave.betas4, "SE(wave)" = se.betas4,
        "Value(clic-h-ecdf)" = ave.betas5, "SE(clic-h-ecdf)" = se.betas5)
  ###cbind("Value(ave)" = ave.betas.Florios) #, "SE(ave)" = se.betas, 
  #    "Value(wave)" = ave.betas2, "SE(wave)" = se.betas2)
  
}

interval <- function(x) {
  if ( x %% 10 != 0) {
    res <- x %/% 10 + 1
  }
  else {
    x <- x-1  
    res <- x %/% 10 + 1
  }
  return(res)
}

intervalBIG <- function(x) {
  if ( x %% 20 != 0) {
    res <- x %/% 20 + 1
  }
  else {
    x <- x-1  
    res <- x %/% 20 + 1
  }
  return(res)
}

computeWeightMatrixAVE_by6 <- function() {
  
  #A is 18 x 42 so that se are also computed
  
  A <- matrix(0,20,120)
  
  #A <- matrix(0,Q*5,P*(5+5))
  
  A[1,c(1,21,41,71,91,111)]=1
  
  A[2,c(3,23,43,72,92,112)]=1
  
  A[3,c(5,25,45,73,93,113)]=1
  
  A[4,c(6,26,46,74,94,114)]=1
  
  A[5,c(7,27,47,75,95,115)]=1
  
  A[6,c(2,31,51,61,81,116)]=1
  
  A[7,c(4,32,52,63,83,117)]=1
  
  A[8,c(8,33,53,65,85,118)]=1
  
  A[9,c(9,34,54,66,86,119)]=1
  
  A[10,c(10,35,55,67,87,120)]=1
  
  A[11,c(11,22,56,62,96,101)]=1
  
  A[12,c(12,24,57,64,97,103)]=1
  
  A[13,c(13,28,58,68,98,105)]=1
  
  A[14,c(14,29,59,69,99,106)]=1
  
  A[15,c(15,30,60,70,100,107)]=1
  
  A[16,c(16,36,42,76,82,102)]=1
  
  A[17,c(17,37,44,77,84,104)]=1
  
  A[18,c(18,38,48,78,88,108)]=1
  
  A[19,c(19,39,49,79,89,109)]=1
  
  A[20,c(20,40,50,80,90,110)]=1
  
  ## now also for rho's
  ## skip rho's
  
  ##finalize A
  
  A<- (1/6)*A
  
  A
}



computeWeightMatrixAVE <- function() {
  
  #A is 18 x 42 so that se are also computed
  
  A <- matrix(0,20,60)
  
  #A <- matrix(0,Q*5,P*(5+5))
  
  A[1,1]=1
  A[1,11]=1
  A[1,21]=1
  
  A[2,3]=1
  A[2,13]=1
  A[2,23]=1
  
  A[3,5]=1
  A[3,15]=1
  A[3,25]=1
  
  A[4,6]=1
  A[4,16]=1
  A[4,26]=1
  
  A[5,7]=1
  A[5,17]=1
  A[5,27]=1
  
  A[6,2]=1
  A[6,31]=1
  A[6,41]=1
  
  A[7,4]=1
  A[7,33]=1
  A[7,43]=1
  
  A[8,8]=1
  A[8,35]=1
  A[8,45]=1
  
  A[9,9]=1
  A[9,36]=1
  A[9,46]=1
  
  A[10,10]=1
  A[10,37]=1
  A[10,47]=1
  
  A[11,12]=1
  A[11,32]=1
  A[11,51]=1
  
  A[12,14]=1
  A[12,34]=1
  A[12,53]=1
  
  A[13,18]=1
  A[13,38]=1
  A[13,55]=1
  
  A[14,19]=1
  A[14,39]=1
  A[14,56]=1
  
  A[15,20]=1
  A[15,40]=1
  A[15,57]=1
  
  A[16,22]=1
  A[16,42]=1
  A[16,52]=1
  
  A[17,24]=1
  A[17,44]=1
  A[17,54]=1
  
  A[18,28]=1
  A[18,48]=1
  A[18,58]=1
  
  A[19,29]=1
  A[19,49]=1
  A[19,59]=1
  
  A[20,30]=1
  A[20,50]=1
  A[20,60]=1
  
  ## now also for rho's
  ## skip rho's
  
  ##finalize A
  
  A<- (1/3)*A
  
  A
}


computeWeightMatrixWAVE <- function(A2) {
  
  #A2 is 16 x 48 so that se are also computed
  A3 <- matrix(0,20,60)
  #A2 <- matrix(0,Q*5,P*(5+5))
  
  
  nzPattern <- computeWeightMatrixAVE()
  
  #the logic is to loop on the columns of A2 and fill in elements (1 for each column), in the spot of nzPattern with diagOmega elements
  for (j in 1:60) {
    for (i in 1:20) {
      if (nzPattern[i,j]!=0) { 
        A3[i,j] <- colSums(A2)[j]  # an easy way to assign the unique element of each column in A2 to A3 suitable position
      }
    }    
  }
  
  ##finalize A3, so that rows add up to 1
  
  valueDenom <- rowSums(A3)
  for (i in 1:20) {
    A3[i,]<- A3[i,] / valueDenom[i]
    
  }  
  
  
  A3
  
}


hess.bin.Vect <- function (thetas, id, y, X, Z, GHk = GHk, extraParam) {
  #cddVect(thetas, logLik.bin, id = id, y = y, X = X, Z = Z, GHk = GHk)
  res <- matrix(nrow=length(thetas),ncol=length(thetas))
  
  #for (i in 1:length(thetas)) {
  #for (j in 1:length(thetas)) {
  #thetasIJ <- thetas
  res <- cddVect(thetas, logLik.bin, id = id, y = y, X = X, Z = Z, GHk = GHk, extraParam)
  #}
  #}
  
  res
}



cddVect <- function (x0, f, ..., eps = 0.0005) {
  
  # Translate Matlab to R from URL:
  # http://grizzly.la.psu.edu/~suj14/programs.html
  # Matlab code: SUNG Jae Jun, PhD
  # R code: Kostas Florios, PhD
  
  # Matlab comments
  
  #%   Compute the Hessian of a real-valued function numerically
  #%   This is a translation of the Gauss command, hessp(fun,x0), considering only
  #%   the real arguments.
  #%   f: real-valued function (1 by 1)
  #%   x0: k by 1, real vector
  #%   varargin: various passing arguments
  #%   H: k by k, Hessian of f at x0, symmetric matrix
  
  #initializations
  xx0 <- x0
  k <- length(x0)
  x0 <- matrix(0,k,1)
  x0[,1] <- xx0
  dax0 <- matrix(0,k,1)
  hessian <- matrix(0,nrow=k,ncol=k)
  grdd <- matrix(0,nrow=k,ncol=1)
  #eps <- 6.0554544523933429e-6
  eps <- 0.0005
  H <- matrix(0,nrow=k,ncol=k)
  
  # Computaion of stepsize (dh)
  ax0=abs(x0)
  for (i in 1:k) {
    if (x0[i,1] != 0) {
      dax0[i,1] <- x0[i,1]/ax0[i,1]
    }
    else {
      dax0[i,1] <- 1
    }
  }
  
  #dh <- eps*max(ax0, (1e-2)*matrix(1,k,1))*dax0
  dh <- eps*dax0
  xdh=x0+dh
  dh=xdh-x0;  # This increases precision slightly
  ee <- matrix(0,nrow=k,ncol=k)
  I <- diag(1,k)
  for (i in 1:k) {
    ee[,i] <- I[,i]*dh
  }
  
  # Computation of f0=f(x0)
  f0 <- f(x0, ...)
  
  # Compute forward step
  for (i in 1:k) {
    grdd[i,1] <- f(x0+ee[,i], ...)
  }
  
  # Compute 'double' forward step
  for (i in 1:k) {
    cat("Computing Row No...:", i, "of hessian\n")
    for (j in i:k) {
      hessian[i,j] <- f(x0+(ee[,i]+ee[,j]), ...)
      if ( i!=j) {
        hessian[j,i] <- hessian[i,j]
      }
    }
  }
  
  l <- t(matrix(1,k,1))
  grdd <- kronecker(l,grdd)
  
  #H <- (((hessian - grdd) - t(grdd)) + f0[1,1]*matrix(1,nrow=k,ncol=k) ) / kronecker(dh,t(dh))
  H <- (((hessian - grdd) - t(grdd)) +  f0*matrix(1,nrow=k,ncol=k) ) / kronecker(dh,t(dh))
  
  return(H)
}


score.bin.BIG <- function (thetas, id, y, X, Z, GHk = 5, extraParam, Data=Data, r=r, i=i) {
  fd (thetas, logLik.bin.BIG, id = id, y = y, X = X, Z = Z, GHk = GHk, extraParam, Data=Data, r=r, i=i)
}

hess.bin.BIG.Vect <- function (thetas, id, y, X, Z, GHk = GHk, extraParam, Data, r, i) {
  #cddVect(thetas, logLik.bin, id = id, y = y, X = X, Z = Z, GHk = GHk)
  res <- matrix(nrow=length(thetas),ncol=length(thetas))
  
  #for (i in 1:length(thetas)) {
  #for (j in 1:length(thetas)) {
  #thetasIJ <- thetas
  res <- cddVect(thetas, logLik.bin.BIG, id = id, y = y, X = X, Z = Z, GHk = GHk, extraParam, Data, r, i)
  #}
  #}
  
  res
}


logLik.bin.BIG <- function (thetas, id, y, X, Z, GHk = 5, extraParam, Data, r, i) {
  
  #environment(aveThetas)<- environment(logLik.bin) <- environment(score.bin) <- environment()
  #environment(aveThetas)<-  environment()
  environment(aveThetas2)<-  environment()
  
  
  res<-0
  cp<-0
  ##thetasP <- vector("list", P) 
  #1: pair 1-2
  ###thetasP[[1]] = c( thetas[c(1,2,5,6,9,13)],  thetas[c(17)], thetas[c(10,14,18)] )
  #2: pair 1-3
  ###thetasP[[2]] = c( thetas[c(1,3,5,7,9,13)],  thetas[c(17)], thetas[c(11,15,19)] )
  #3: pair 1-4
  ###thetasP[[3]] = c( thetas[c(1,4,5,8,9,13)],  thetas[c(17)], thetas[c(12,16,20)] )
  #4: pair 2-3
  ###thetasP[[4]] = c( thetas[c(2,3,6,7,10,14)], thetas[c(18)], thetas[c(11,15,19)] )
  #5: pair 2-4
  ###thetasP[[5]] = c( thetas[c(2,4,6,8,10,14)], thetas[c(18)], thetas[c(12,16,20)] )
  #6: pair 3-4
  ###thetasP[[6]] = c( thetas[c(3,4,7,8,11,15)], thetas[c(19)], thetas[c(12,16,20)] )
  
  ##thetasP <- thetas
  ##thetasP <- as.vector(thetasP)
  thetasP <- thetas[1:10]
  thetasP <- as.vector(thetasP)
  thQ     <- thetas[11:20]
  thQ     <- as.vector(thQ)
  
  p <- r
  
  #DD <- do.call(rbind, list(Data[1:3], Data[1:3]))
  #DD$outcome <- gl(2, nrow(Data))
  
  p<-r
  prs <- pairs[, p]
  
  DD <- Data
  
  res <-       logLik.bin(thetasP,
                          id = DD$id, y = DD$y, 
                          X = rbind(model.matrix(Models[[p]])[1:(dim(Data)[1]/2),],
                                    model.matrix(Models[[p]])[(n*length(times)+1):(n*length(times)+(dim(Data)[1])/2),]),
                          Z = rbind(model.matrix(Models[[p]])[1:(dim(Data)[1]/2),],
                                    model.matrix(Models[[p]])[(n*length(times)+1):(n*length(times)+(dim(Data)[1]/2)),]),                  
                          GHk = GHk,
                          extraParam = extraParam)  
  
  indices <- 1:Q
  indic   <- indices[c(-prs[1],-prs[2])]
  
  ##thetasQ <- vector("list", Q) 
  #1: item 1
  ##thetasQ[[1]] = c(fixef(ModelsOne[[1]]),Sigma2YiRIone[[1]],Sigma2YiRSone[[1]],Sigma2YiRIRSone[[1]])
  #2: item 2
  ##thetasQ[[2]] = c(fixef(ModelsOne[[2]]),Sigma2YiRIone[[2]],Sigma2YiRSone[[2]],Sigma2YiRIRSone[[2]])
  #3: item 3
  ##thetasQ[[3]] = c(fixef(ModelsOne[[3]]),Sigma2YiRIone[[3]],Sigma2YiRSone[[3]],Sigma2YiRIRSone[[3]])
  #4: item 4
  ##thetasQ[[4]] = c(fixef(ModelsOne[[4]]),Sigma2YiRIone[[4]],Sigma2YiRSone[[4]],Sigma2YiRIRSone[[4]])
  thetasQ <- vector("list", 2) 
  thetasQ[[1]] <- thetas[11:15]
  thetasQ[[2]] <- thetas[16:20]
  
  DDD <- do.call(rbind, list(Data[1:3]))
  DDD$outcome <- gl(1, nrow(Data))
  
  
  for (k in 1:length(indic) ) {
    ii <- indic[k]
    if (i <= n) {
      ind.i <- DataRaw$id == i
    }
    else {
      ind.i <- rep(T,n*length(times))
    }
    yi <- DataRaw[[paste("y", indic[k], sep = "")]][ind.i]
    DDD <- DataRaw[ind.i,]
    DDD$y <- c(yi)
    res <- res + logLik.bin.One(thetasQ[[k]],
                                id = DDD$id, y=DDD$y,
                                X = model.matrix(ModelsOne[[ii]])[1:(dim(Data)[1]/2),],
                                Z = model.matrix(ModelsOne[[ii]])[1:(dim(Data)[1]/2),],
                                GHk = GHk,
                                extraParamOne = extraParamOne[[ii]])                            
    
  }
  
  res  
}

logLik.bin.One <- function (thetas, id, y, X, Z, GHk = 5, extraParamOne) {
  #thetas <- relist(thetas, lis.thetas)
  #betas <- thetas$betas
  #ncz <- ncol(Z)
  #D <- matrix(0, ncz, ncz) 
  #D[lower.tri(D, TRUE)] <- thetas$D
  #D <- D + t(D)
  #diag(D) <- diag(D) / 2
  #
  betas<-thetas[1:2]
  Sigma2YiRI<-thetas[3]
  Sigma2YiRS<-thetas[4]
  Sigma2YiRIRS<-thetas[5]
  
  
  Dold <- matrix(0,ncol=2,nrow=2)
  Dnew <- matrix(0,ncol=2,nrow=2)
  
  Dold <- extraParamOne   # extraParam is VarCorr(Models[[i]])$id
  #re-align columns, rows
  #old: glmer, 1,2,3,4
  #new: Florios, 1,3,2,4
  
  #now plug-in the thetas, in order to perform numerical derivatives wrt theta (scores, hessians)
  Dnew[1,1] = Sigma2YiRI
  Dnew[2,2] = Sigma2YiRS
  Dnew[2,1] = Sigma2YiRIRS
  Dnew[1,2] = Dnew[2,1]
  
  #simplify notation: D=Dnew, and proceed
  
  D <-nearPD(Dnew)
  ncz <- ncol(Z)
  GH <- gauher(GHk)
  b <- as.matrix(expand.grid(rep(list(GH$x), ncz)))
  dimnames(b) <- NULL
  k <- nrow(b)
  wGH <- as.matrix(expand.grid(rep(list(GH$w), ncz)))
  wGH <- 2^(ncz/2) * apply(wGH, 1, prod) * exp(rowSums(b * b))
  b <- sqrt(2) * b
  Ztb <- Z %*% t(b)
  #
  mu.y <- plogis(as.vector(X %*% betas) + Ztb)
  logBinom <- dbinom(y, 1, mu.y, TRUE)
  log.p.yb <- rowsum(logBinom, id)
  log.p.b <- dmvnorm(b, rep(0, ncol(Z)), D, TRUE)
  p.yb <- exp(log.p.yb + rep(log.p.b, each = nrow(log.p.yb)))
  p.y <- c(p.yb %*% wGH)
  #-sum(log(p.y), na.rm = TRUE)  #logLik as min, original
  sum(log(p.y), na.rm = TRUE)  #logLik as max, CLIC heuristic
}



gc()
