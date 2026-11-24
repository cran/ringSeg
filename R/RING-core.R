# single change-point
Uscan.point <- function(n, Ebynode, sumR, sumRi){
  # "n" is the total number of nodes.
  # Ebynode[[i]][['node']] is the vector of nodes that are connect to i by an edge.
  # Ebynode[[i]][['weight']] is the vector of weights that are connect to i by an edge.
  # The nodes are numbered by their order in the sequence.
  g = c(0,rep(1,n-1))
  U0 = U1 = U2 = rep(0,n-1)
  U0[1] = 2*sum(Ebynode[[1]][['weight']])
  U2[1] = sumR - U0[1]
  for(i in 2:(n-1)){
    g[i] = 0  # update g
    links = Ebynode[[i]][['node']]
    if(length(links)>0){
      add = sum(Ebynode[[i]][['weight']][which(rep(g[i],length(links)) != g[links])])
      subtract = sumRi[i] - add
      #print(c(i,add,subtract))
      U0[i] = U0[i-1] + 2*(add - subtract)
      U1[i] = U1[i-1] + 2*subtract
    }else{
      U0[i] = U0[i-1]
      U1[i] = U1[i-1]
    }
    U2[i] = sumR - U0[i] - U1[i]
  }
  tt = 1:(n-1)
  Uw = (n-tt-1)/(n-2)*U1 + (tt-1)/(n-2)*U2
  Udiff = U1 - U2
  return(list('Uw'=Uw,'Udiff'=Udiff))
}

# changed interval
Uscan.interval <- function(n, Ebynode, sumR, sumRi){
  # "n" is the total number of nodes.
  # Ebynode[[i]][['node']] is the vector of nodes that are connect to i by an edge.
  # Ebynode[[i]][['weight']] is the vector of weights that are connect to i by an edge.
  # The nodes are numbered by their order in the sequence.
  U0 = U1 = U2 = Uw = matrix(0,n,n)
  for (i in 1:(n-1)){
    g = rep(0,n)
    j = i+1
    g[j] = 1 # update g
    U0[i,j] = 2*sumRi[j]
    U2[i,j] = sumR - 2*sumRi[j]
    
    for (j in seq.int(i+1L,n)[-1L]){
      g[j] = 1 # update g

      links = Ebynode[[j]][['node']]
      if (length(links)>0){
        #add = sum(Ebynode[[j]][['weight']][which(rep(g[j],length(links)) != g[links])])
        add = sum(Ebynode[[j]][['weight']][rep(g[j],length(links)) != g[links]])
        subtract = sumRi[j] - add
        U0[i,j] = U0[i,j-1] + 2*(add-subtract)
        U1[i,j] = U1[i,j-1] + 2*subtract
      }else{
        U0[i,j] = U0[i,j-1]
        U1[i,j] = U1[i,j-1]
      }
      
      U2[i,j] = sumR - U0[i,j] - U1[i,j]
      Uw[i,j] = ((n-j+i-1)/(n-2))*U1[i,j] + (j-i-1)/(n-2)*U2[i,j]
    }
  }
  Udiff = U1 - U2
  return(list('Uw'=Uw[upper.tri(Uw)],'Udiff'=Udiff[upper.tri(Udiff)]))
}

#alternative='single' or 'interval'.
#If alternative='interval', l0=n0 and l1=n1.
#B is the number of permutations; if B=0, don't perform permutation.
rcpd<-function(R, n0=NULL, n1=NULL, pval.appr=TRUE, skew.corr=TRUE, B=0, alternative=c('single','interval')){
  alternative = match.arg(alternative)
  if(!is.logical(pval.appr) || length(pval.appr) != 1L || is.na(pval.appr)){
    stop("'pval.appr' must be a single TRUE or FALSE.")
  }
  if(!is.logical(skew.corr) || length(skew.corr) != 1L || is.na(skew.corr)){
    stop("'skew.corr' must be a single TRUE or FALSE.")
  }
  if(length(B) != 1L || !is.numeric(B) || !is.finite(B) || B < 0 || B != floor(B)){
    stop("'B' must be a single non-negative integer.")
  }
  R = as.matrix(R)
  if(nrow(R) != ncol(R)){
    stop("'R' must be a square n x n rank/weight matrix.")
  }
  if(!is.numeric(R) || any(!is.finite(R))){
    stop("'R' must be a finite numeric rank/weight matrix (no NA/NaN/Inf).")
  }
  if(any(R < 0)){
    stop("'R' must have non-negative weights.")
  }
  if(any(diag(R) != 0)){
    stop("'R' must have a zero diagonal.")
  }
  if(any(R != t(R))){
    if(!isTRUE(all.equal(R, t(R)))){
      warning("'R' is not symmetric; symmetrized as (R + t(R))/2.")
    }
    R = R/2 + t(R)/2
  }
  scale = max(R)
  if(scale == 0){
    stop("'R' is empty (all entries zero); no edges to scan. Check the graph / k.")
  }
  n = nrow(R)
  if(n < 4){
    stop("RING change-point detection needs at least 4 observations (got n = ", n, ").")
  }
  # Adding a constant to all off-diagonal weights adds a deterministic
  # term at each split size. Remove a large common offset before summing.
  offset = min(R[upper.tri(R)])
  spread = scale-offset
  if(spread > 0 && spread <= 1e-3*scale){
    diag(R) = offset
    # Binary scaling preserves the represented differences and exact ties.
    R = (R-offset)/2^floor(log2(spread))
  }else if(scale < 1e-50 || scale > 1e50){
    # Protect products of up to three weights at extreme magnitudes.
    R = R/scale
  }
  if(!is.null(n0)){
    if(length(n0) != 1L || !is.numeric(n0) || !is.finite(n0))
      stop("'n0' must be NULL or a single finite number.")
    if(n0 < 0 || n0 > n)
      stop("'n0' must be between 0 and n (", n, ").")
    if(n0 != round(n0)){
      warning("'n0' rounded to the nearest integer (", round(n0), ").")
      n0 = round(n0)
    }
  }
  if(!is.null(n1)){
    if(length(n1) != 1L || !is.numeric(n1) || !is.finite(n1))
      stop("'n1' must be NULL or a single finite number.")
    if(n1 < 0 || n1 > n)
      stop("'n1' must be between 0 and n (", n, ").")
    if(n1 != round(n1)){
      warning("'n1' rounded to the nearest integer (", round(n1), ").")
      n1 = round(n1)
    }
  }
  if(is.null(n0)){
    n0 = 0.05*n
    n0 = ceiling(n0)
  }
  if(is.null(n1)){
    n1 = 0.95*n
    n1 = floor(n1)
  }
  if(n0<2){
    message("  Note: Starting index has been set to n0 = 2 as the weighted edge-count test statistic is not well-defined for t<2.")
    n0=2
  }
  if(n1>(n-2)){
    message("  Note: Ending index has been set to n1 = ", n-2, " as the weighted edge-count test statistic is not well-defined for t> ", n-2, ".")
    n1=n-2
  }
  if(n0>n1){
    stop("Invalid scan range: n0 (", n0, ") must be <= n1 (", n1, ").")
  }

  Ebynode = lapply(1:n, function(i){
    ei = which(R[i,] > 0)
    list('node'=ei,'weight'=R[i,ei])
  })
  
  res = list() #save result

  R2 = R*R; sumR = sum(R); sumR2 = sum(R2)
  sumRi = rowSums(R)
  
  r0 = sumR/(n*(n-1)); rd = sqrt(sumR2/(n*(n-1)))
  r1 = sqrt( mean(sumRi^2)/(n-1)^2)
  
  var.w = (n-2)*rd^2 - 2*(n-1)*r1^2 + n*r0^2
  var.diff = r1^2 - r0^2
  # Recompute close to zero using sums of squares to avoid cancellation.
  eps = 32*.Machine$double.eps
  if(var.diff <= eps*r1^2){
    var.diff = mean((sumRi - mean(sumRi))^2)/(n-1)^2
  }
  if(var.w <= eps*n*rd^2){
    row.effect = (sumRi - mean(sumRi))/(n-2)
    residual = R - r0 - outer(row.effect,row.effect,"+")
    diag(residual) = 0
    residual[abs(residual) <= eps*max(R)] = 0
    var.w = (n-2)*sum(residual^2)/(n*(n-1))
  }

  tt = 1:(n-1)
  muw = n*( tt-1 )*( n - tt-1 )*r0/(n-2)
  mudiff = (n-1)*( 2*tt - n)*r0
  sd_w = sqrt( 2*tt*( tt-1 )*(n - tt)*( n - tt-1 )*var.w/( (n-2)^2*(n-3)) )
  sd_diff = sqrt( 4*(n-1)*tt*(n - tt)*var.diff )
  score = function(U, indices, sizes){
    Zw = if(var.w == 0) rep(0,length(indices)) else (U$Uw[indices]-muw[sizes])/sd_w[sizes]
    Zdiff = if(var.diff == 0) rep(0,length(indices)) else (U$Udiff[indices]-mudiff[sizes])/sd_diff[sizes]
    list(TR=Zdiff^2 + Zw^2, MR=pmax(abs(Zdiff),Zw), WR=Zw)
  }

  permutation_p = function(permuted, observed){
    # Different addition orders can move an exact tie by a few ulps.
    tolerance = 32*.Machine$double.eps*max(1,abs(observed))
    (sum(permuted >= observed-tolerance)+1)/(B+1)
  }

  analytic = function(scanZ){
    out = list(TR=NA_real_, MR=NA_real_, WR=NA_real_)
    if(skew.corr){
      out$MR.cor = NA_real_
      out$WR.cor = NA_real_
    }
    if(n0 == n1 || var.w == 0 || var.diff == 0){
      warning("Analytic p-values are unavailable for a fixed scan bound or zero null variance; use B > 0 for permutation p-values.")
      return(out)
    }

    # The leading factors of the large-threshold approximations peak at
    # these values. Use only the upper side of those factors.
    minimum = if(alternative == "single") c(2,1,1) else c(4,sqrt(3),sqrt(3))
    observed = c(scanZ$TR$TR.max,scanZ$MR$MR.max,scanZ$WR$WR.max)
    eligible = observed >= minimum
    if(any(!eligible)){
      omitted = paste(names(out)[seq_along(eligible)][!eligible],collapse=", ")
      warning("Scan maximum below the analytic tail range for ",omitted,
              "; returning NA. Use B > 0 for permutation p-values.")
    }
    if(any(eligible)){
      statistics = c("generalized","max","weighted")[eligible]
      correct = skew.corr && n>=6 && any(eligible[-1L])
      if(alternative == "single"){
        values = pASingle(R,scanZ,Ebynode,statistics=statistics,skew.corr=correct,n0=n0,n1=n1)
      }else{
        values = pAInterval(R,scanZ,Ebynode,statistics=statistics,skew.corr=correct,l0=n0,l1=n1)
      }
      out[names(values)] = values
    }
    if(skew.corr && n<6){
      warning("Skewness-corrected p-values require at least six observations; use B > 0 for permutation p-values.")
    }
    out
  }
  scanZ = list()
  
  if(alternative=='single'){
    U = Uscan.point(n, Ebynode, sumR,sumRi)
    valid.sizes = 2:(n-2)
    scores = score(U,valid.sizes,valid.sizes)
    TR = c(NA_real_,scores$TR,NA_real_)
    MR = c(NA_real_,scores$MR,NA_real_)
    WR = c(NA_real_,scores$WR,NA_real_)
    
    #Time interval
    T0 = n0:n1
    id0 = which.max(TR[T0])
    scanZ$TR$tau = T0[id0]
    scanZ$TR$TR.max = TR[T0][id0]
    scanZ$TR$TR = TR
    
    id0 = which.max(MR[T0])
    scanZ$MR$tau = T0[id0]
    scanZ$MR$MR.max = MR[T0][id0]
    scanZ$MR$MR = MR
    
    id0 = which.max(WR[T0])
    scanZ$WR$tau = T0[id0]
    scanZ$WR$WR.max = WR[T0][id0]
    scanZ$WR$WR = WR
    
    if(var.w == 0) scanZ$WR$tau = NA_integer_
    if(var.w == 0 && var.diff == 0){
      scanZ$TR$tau = scanZ$MR$tau = NA_integer_
    }
    res$scanZ = scanZ
    
    if(pval.appr==TRUE){
      res$pval.appr = analytic(scanZ)
    }
    if(B>0){
      permmatch = rep(0,n)
      Z.B = sapply(1:B, function(b){
        perm = sample(n)
        permmatch[perm] = 1:n
        Ebnstar =  lapply(1:n, function(i){
          oldlinks = Ebynode[[permmatch[i]]][['node']]
          list('node'= perm[oldlinks],'weight'=Ebynode[[permmatch[i]]][['weight']])
        })
        Ub = Uscan.point(n,Ebnstar,sumR,sumRi[permmatch])
        unname(vapply(score(Ub,T0,T0),max,numeric(1)))
      })
      
      res$pval.perm = list('TR'=permutation_p(Z.B[1,],scanZ$TR$TR.max),
                           'MR'=permutation_p(Z.B[2,],scanZ$MR$MR.max),
                           'WR'=permutation_p(Z.B[3,],scanZ$WR$WR.max))
      res$perm.dat = list('TR'=Z.B[1,],
                          'MR'=Z.B[2,],
                          'WR'=Z.B[3,])
    }
    return(res)
  }
  if(alternative=='interval'){
    l0 = n0; l1 = n1
    U = Uscan.interval(n, Ebynode, sumR, sumRi)
    candidates = which(upper.tri(R),arr.ind=TRUE)
    dif = candidates[,2]-candidates[,1]
    ids = which(dif>=l0 & dif<=l1)

    scores = score(U,ids,dif[ids])
    scanZ = lapply(names(scores),function(statistic){
      values = scores[[statistic]]
      best = which.max(values)
      result = list(tau=as.numeric(candidates[ids[best],]))
      result[[paste0(statistic,".max")]] = values[best]
      result
    })
    names(scanZ) = names(scores)
    if(var.w == 0) scanZ$WR$tau = rep(NA_integer_,2L)
    if(var.w == 0 && var.diff == 0){
      scanZ$TR$tau = scanZ$MR$tau = rep(NA_integer_,2L)
    }
    res$scanZ = scanZ

    if(pval.appr==TRUE){
      res$pval.appr = analytic(scanZ)
    }
    if(B>0){
      permmatch = rep(0,n)
      Z.B = sapply(1:B, function(b){
        perm = sample(n)
        permmatch[perm] = 1:n
        Ebnstar =  lapply(1:n, function(i){
          oldlinks = Ebynode[[permmatch[i]]][['node']]
          list('node'= perm[oldlinks],'weight'=Ebynode[[permmatch[i]]][['weight']])
        })
        Ub = Uscan.interval(n,Ebnstar,sumR,sumRi[permmatch])
        unname(vapply(score(Ub,ids,dif[ids]),max,numeric(1)))
      })
      
      res$pval.perm = list('TR'=permutation_p(Z.B[1,],scanZ$TR$TR.max),
                           'MR'=permutation_p(Z.B[2,],scanZ$MR$MR.max),
                           'WR'=permutation_p(Z.B[3,],scanZ$WR$WR.max))
      res$perm.dat = list('TR'=Z.B[1,],
                          'MR'=Z.B[2,],
                          'WR'=Z.B[3,])
    }
    return(res)
  }
}

# the Nu function
Nu = function(x){
  y = x/2
  (1/y)*(pnorm(y)-0.5)/(y*pnorm(y) + dnorm(y))
}

# hw(n,t) = hw(n,t/n)
rho_one_Rw = function(n, t){
  (n-1)*(2*t^2 - 2*n*t + n)/( 2*t*(n-t)*(t-1)*(t+1-n) )
}

# The p-value of TR, WR and MR: are independent of the graph, and only depend on n, n0, n1 and the statistics b.
# Single Change point 
pTR1 = function(b, n, n0=ceiling(0.05*n), n1=floor(0.95*n)){
  if(b>0){
    integrandG = function(t,w){
      x1 = n/(2*t*(n - t))
      x2 = rho_one_Rw(n,t)
      (x1*cos(w)^2+x2*sin(w)^2)*b*Nu(sqrt(2*b*(x1*cos(w)^2+x2*sin(w)^2)))/(pi)
    }
    integrand0 = function(t) {integrate(integrandG,0,2*pi,t=t,subdivisions=3000, stop.on.error=FALSE)$value}
    pval.generalized = dchisq(b,2)*integrate(Vectorize(integrand0),n0,n1,subdivisions=3000, stop.on.error=FALSE)$value
    return(min(pval.generalized,1))
  }else{
    return(1)
  }
}
# Changed Interval 
pTR2 = function(b, n, l0=ceiling(0.05*n), l1=floor(0.95*n)){
  if (b>0){
    integrandG = function(t,w){
      x1 = n/(2*t*(n - t))
      x2 = rho_one_Rw(n,t)
      (n-t)*(2*(x1*cos(w)^2+x2*sin(w)^2) *b*Nu(sqrt(2*b*(x1*cos(w)^2+x2*sin(w)^2))))^2/(2*pi)
    }
    integrand0 = function(t) {integrate(integrandG,0,2*pi,t=t,subdivisions=3000, stop.on.error=FALSE)$value}
    pval.generalized = dchisq(b,2)*integrate(Vectorize(integrand0), l0, l1, subdivisions=3000, stop.on.error=FALSE)$value
    return(min(pval.generalized,1))
  }else{
    return(1)
  }
}

# Single Change point 
pWR1 = function(b, n, n0=ceiling(0.05*n), n1=floor(0.95*n)){
  if (b>0){
    integrandW = function(t){
      x = rho_one_Rw(n,t)
      x*Nu(sqrt(2*b^2*x))
    }
    pval.weighted = dnorm(b)*b*integrate(integrandW, n0, n1, subdivisions=3000, stop.on.error=FALSE)$value
    return(min(pval.weighted,1))
  }else{
    return(1)
  }
}
# Changed Interval 
pWR2 = function(b, n, l0=ceiling(0.05*n), l1=floor(0.95*n)){
  if (b>0){
    integrandW = function(t){
      x = rho_one_Rw(n,t)
      (b^2*x*Nu(sqrt(2*b^2*x)))^2*(n-t)
    }
    pval.weighted = try(dnorm(b)/b*integrate(integrandW, l0, l1, subdivisions=3000, stop.on.error=FALSE)$value, silent=T)
    return(min(pval.weighted,1))
  }else{
    return(1)
  }
}

# Single Change point 
pDR1 = function(b, n, n0=ceiling(0.05*n), n1=floor(0.95*n)){
  if (b>0){
    integrand1 = function(t){
      x1 = n/(2*t*(n - t))
      x1*Nu(sqrt(2*b^2*x1))
    }
    pval.diff = 2*dnorm(b)*b*integrate(integrand1, n0, n1, subdivisions=3000, stop.on.error=FALSE)$value
    return(min(pval.diff,1))
  }else{
    return(1)
  }
}
# Changed Interval 
pDR2 = function(b, n, l0=ceiling(0.05*n), l1=floor(0.95*n)){
  if (b>0){
    integrand1 = function(t){
      x1 = n/(2*t*(n - t))
      (b^2*x1*Nu(sqrt(2*b^2*x1)))^2*(n-t)
    }
    pval.diff = try(2*dnorm(b)/b*integrate(integrand1, l0, l1, subdivisions=3000, stop.on.error=FALSE)$value,silent=T)
    return(min(pval.diff,1))
  }else{
    return(1)
  }
}

# Single Change point 
pMR1 = function(b, n, n0=ceiling(0.05*n), n1=floor(0.95*n)){
  if (b>0){
    pval.w = pWR1(b,n,n0,n1)
    pval.diff = pDR1(b,n,n0,n1)
    pval.max.type = pval.w + (1-pval.w)*pval.diff
    return(pval.max.type)
  }else{
    return(1)
  }
}
# Changed Interval 
pMR2 = function(b, n, l0=ceiling(0.05*n), l1=floor(0.95*n)){
  if (b>0){
    pval.w = pWR2(b,n,l0,l1)
    pval.diff = pDR2(b,n,l0,l1)
    pval.max.type = pval.w + (1-pval.w)*pval.diff
    return(pval.max.type)
  }else{
    return(1)
  }
}

# Test: the critical values should be the same as those in Lynna's paper
# n=1000; n0=25; n1=n-n0; l0=50; l1 = n-l0 
# c(pTR1(14.11,n,n0,n1),pWR1(3.14,n,n0,n1),pMR1(3.38,n,n0,n1))
# c(pTR2(23.96,n,l0,l1),pWR2(4.22,n,l0,l1),pMR2(4.34,n,l0,l1))

# p-value approximation for single change-point
pASingle <- function(R, scanZ, Ebynode, statistics="all", skew.corr=FALSE, n0=NULL, n1=NULL){
  n = nrow(R)
  output = list()
  if (length(which(!is.na(match(c("g","generalized","all"), statistics))))>0){
    b = scanZ$TR$TR.max
    output$TR = pTR1(b,n,n0,n1)
  }
  #  if (skew.corr==FALSE){
  if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
    b = scanZ$MR$MR.max
    output$MR = pMR1(b,n,n0,n1)
  }
  if (length(which(!is.na(match(c("w","weighted","all"), statistics))))>0){
    b = scanZ$WR$WR.max
    output$WR = pWR1(b,n,n0,n1)
  }
  if (skew.corr){
    ### Calculate the third moment ###
    thirdlist = ThirdMoment(R,Ebynode)
    EZw3 = thirdlist$EZw3
    EZdiff3 = thirdlist$EZdiff3
    TI = 1:(n-1)

    ### Calculate the corrected p-value for Zw ###
    x.w = rho_one_Rw(n,TI)
    x.w[is.infinite(x.w)] = 0
    if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
      # Use the observed max-type threshold for both component tails.
      b = scanZ$MR$MR.max
      p.w.corr = pS_correct_w(n,b,EZw3,x.w,n0,n1)
      x.diff = n/(2*TI*(n - TI))
      p.diff.corr = pS_correct_diff(n,b,EZdiff3,x.diff,n0,n1)
      output$MR.cor = p.w.corr + (1-p.w.corr)*p.diff.corr
    }
    if (length(which(!is.na(match(c("w","weighted","all"), statistics))))>0){
      b = scanZ$WR$WR.max
      output$WR.cor = pS_correct_w(n,b,EZw3,x.w,n0,n1)
    }
  }
  return(output)
}

pAInterval <- function(R, scanZ, Ebynode, statistics="all", skew.corr=FALSE, l0=NULL, l1=NULL){
  n = nrow(R)
  output = list()
  if (length(which(!is.na(match(c("g","generalized","all"), statistics))))>0){
    b = scanZ$TR$TR.max
    output$TR = pTR2(b,n,l0,l1)
  }
  #  if (skew.corr==FALSE){
  if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
    b = scanZ$MR$MR.max
    output$MR = pMR2(b,n,l0,l1)
  }
  if (length(which(!is.na(match(c("w","weighted","all"), statistics))))>0){
    b = scanZ$WR$WR.max
    output$WR = pWR2(b,n,l0,l1)
  }
  if (skew.corr){
    ### Calculate the third moment ###
    thirdlist = ThirdMoment(R,Ebynode)
    EZw3 = thirdlist$EZw3
    EZdiff3 = thirdlist$EZdiff3
    TI = 1:(n-1)

    ### Calculate the corrected p-value for Zw ###
    x.w = rho_one_Rw(n,TI)
    x.w[is.infinite(x.w)] = 0
    if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
      # Use the observed max-type threshold for both component tails.
      b = scanZ$MR$MR.max
      p.w.corr = pI_correct_w(n,b,EZw3,x.w,l0,l1)
      x.diff = n/(2*TI*(n - TI))
      p.diff.corr = pI_correct_diff(n,b,EZdiff3,x.diff,l0,l1)
      output$MR.cor = p.w.corr + (1-p.w.corr)*p.diff.corr
    }
    if (length(which(!is.na(match(c("w","weighted","all"), statistics))))>0){
      b = scanZ$WR$WR.max
      output$WR.cor = pI_correct_w(n,b,EZw3,x.w,l0,l1)
    }
  }
  return(output)
}


##Calculate the corrected p-value for Zdiff
pS_correct_diff <- function(n,b,r,x,n0,n1){
  if (b<=0){
    return(1)
  }
  fallback = function(){
    message("Zdiff: p-value approximation without skewness correction is reported.")
    pDR1(b,n,n0,n1)
  }
  theta_b = rep(NA_real_,n-1)
  pos = which(1+2*r*b>0)
  # Rationalized root remains well-defined and accurate as skewness tends to zero.
  theta_b[pos] = 2*b/(1+sqrt((1+2*r*b)[pos]))
  ratio = exp((b-theta_b)^2/2 + r*theta_b^3/6)/sqrt(1+r*theta_b)
  a = x*Nu(sqrt(2*b^2*x)) * ratio
  
  nn.l = ceiling(n/2)-length(which(1+2*r[1:ceiling(n/2)]*b>0))
  nn.r = ceiling(n/2)-length(which(1+2*r[ceiling(n/2):(n-1)]*b>0))
  if (nn.l>0.35*n || nn.r>0.35*n){
    message("Zdiff: nn.l>0.35*n || nn.r>0.35*n.")
    message("Zdiff: p-value approximation without skewness correction is reported.")
    return(pDR1(b,n,n0,n1))
  }
  if (nn.l>=n0){
    message(paste0('Zdiff: nn.l>=n0, b:',b,'.'))
    neg = which(1+2*r[1:ceiling(n/2)]*b<=0)
    dif = c(diff(neg),n/2-nn.l)
    id1 = which.max(dif)
    id2 = id1 + ceiling(0.03*n)
    id3 = id2 + ceiling(0.09*n)
    if(id2+1L > length(a) || id3 > length(a) ||
       any(!is.finite(a[c(id2,id2+1L,id3)]))){
      return(fallback())
    }
    inc = (a[id3]-a[id2])/(id3-id2)
    a[id2:1] = a[id2+1]-inc*(1:id2)
  }
  if (nn.r>=(n-n1)){
    message(paste0('Zdiff: nn.r>=(n-n1), b:',b,'.'))
    neg = which(1+2*r[ceiling(n/2):(n-1)]*b<=0 )
    id1 = min(neg+ceiling(n/2)-1,ceiling(n/2)-1)
    id2 = id1 - ceiling(0.03*n)
    id3 = id2 - ceiling(0.09*n)
    if(id3 < 1L || id2 < 2L || id2 > length(ratio) ||
       any(!is.finite(ratio[c(id3,id2-1L,id2)]))){
      return(fallback())
    }
    inc = (ratio[id3]-ratio[id2])/(id3-id2)
    ratio[id2:(n-1)] = ratio[id2-1]+inc*((id2:(n-1))-id2)
    ratio[ratio<0]=0
    a[(n/2):(n-1)] = (x*Nu(sqrt(2*b^2*x)) * ratio)[(n/2):(n-1)] # update a after extrapolation 
  }
  neg2 = which(a<0)
  a[neg2] = 0
  if(any(!is.finite(a[n0:n1]))){
    return(fallback())
  }
  integrand = function(s){
    a[s]
  }
  result = tryCatch(2*dnorm(b)*b*integrate(integrand, n0, n1, subdivisions=3000, stop.on.error=FALSE)$value,
                    error=function(e) NA_real_)
  if(!is.finite(result) || result < 0 || result > 1){
    return(fallback())
  }
  return(result)
}

##Calculate the corrected p-value for Zw
pS_correct_w <- function(n,b,r,x,n0,n1){
  if (b<=0){
    return(1)
  }
  fallback = function(){
    message("Zw: p-value approximation without skewness correction is reported.")
    pWR1(b,n,n0,n1)
  }
  theta_b = rep(NA_real_,n-1)
  pos = which(1+2*r*b>0)
  theta_b[pos] = 2*b/(1+sqrt((1+2*r*b)[pos]))
  ratio = exp((b-theta_b)^2/2 + r*theta_b^3/6)/sqrt(1+r*theta_b)
  a = x*Nu(sqrt(2*b^2*x)) * ratio
  a[which(x==0)] = 0
  nn = n-1-length(pos)
  if (nn>0.75*n){
    message("Extrapolation of Zw for skewness-corrected p-value approximation could not be performed.")
    return(pWR1(b, n, n0, n1))
  }
  if (nn>=(n0-1)+(n-n1)){
    message(paste0('Zw: nn>=(n0-1)+(n-n1), b:',b,'.'))
    neg = which(1+2*r*b<=0)
    if(length(neg) < 2L){
      return(fallback())
    }
    dif = diff(neg)
    id1 = which.max(dif)
    id2 = id1 + ceiling(0.03*n)
    id3 = id2 + ceiling(0.09*n)
    if(id2+1L > length(a) || id3 > length(a) ||
       any(!is.finite(a[c(id2,id2+1L,id3)]))){
      return(fallback())
    }
    inc = (a[id3]-a[id2])/(id3-id2)
    a[id2:1] = a[id2+1]-inc*(1:id2)
    a[(n/2+1):n] = a[(n/2):1]
    neg2 = which(a<0)
    a[neg2] = 0
  }
  if(any(!is.finite(a[n0:n1]))){
    return(fallback())
  }
  integrand = function(s){
    a[s]
  }
  result = tryCatch(dnorm(b)*b*integrate(integrand, n0, n1, subdivisions=3000, stop.on.error=FALSE)$value,
                    error=function(e) NA_real_)
  if(!is.finite(result) || result < 0 || result > 1){
    return(fallback())
  }
  return(result)
}

##Calculate the corrected p-value for Zdiff
pI_correct_diff <- function(n,b,r,x,l0,l1){
  if (b<=0){
    return(1)
  }
  fallback = function(){
    message("Zdiff: p-value approximation without skewness correction is reported.")
    pDR2(b,n,l0,l1)
  }
  theta_b = rep(NA_real_,n-1)
  pos = which(1+2*r*b>0)
  theta_b[pos] = 2*b/(1+sqrt((1+2*r*b)[pos]))
  ratio = exp((b-theta_b)^2/2 + r*theta_b^3/6)/sqrt(1+r*theta_b)
  a = (b^2*x*Nu(sqrt(2*b^2*x)))^2 * ratio  ## this is different from single change-point alternative
  
  nn.l = ceiling(n/2)-length(which(1+2*r[1:ceiling(n/2)]*b>0))
  nn.r = ceiling(n/2)-length(which(1+2*r[ceiling(n/2):(n-1)]*b>0))
  if (nn.l>0.35*n || nn.r>0.35*n){
    message("Zdiff: nn.l>0.35*n || nn.r>0.35*n.")
    message("Zdiff: p-value approximation without skewness correction is reported.")
    return(pDR2(b,n,l0,l1))
  }
  neg = which(1+2*r[1:ceiling(n/2)]*b<=0)
  if (nn.l>=l0){
    message(paste0('Zdiff: nn.l>=l0, b:',b,'.'))
    dif = c(diff(neg),n/2-nn.l)
    id1 = which.max(dif)
    id2 = id1 + ceiling(0.03*n)
    id3 = id2 + ceiling(0.09*n)
    if(id2+1L > length(a) || id3 > length(a) ||
       any(!is.finite(a[c(id2,id2+1L,id3)]))){
      return(fallback())
    }
    inc = (a[id3]-a[id2])/(id3-id2)
    a[id2:1] = a[id2+1]-inc*(1:id2)
  }
  neg = which(1+2*r[ceiling(n/2):(n-1)]*b<=0 )
  if (nn.r>=(n-l1)){
    message(paste0('Zdiff: nn.r>=(n-l1), b:',b,'.'))
    id1 = min(neg+ceiling(n/2)-1,ceiling(n/2)-1)
    id2 = id1 - ceiling(0.03*n)
    id3 = id2 - ceiling(0.09*n)
    if(id3 < 1L || id2 < 2L || id2 > length(ratio) ||
       any(!is.finite(ratio[c(id3,id2-1L,id2)]))){
      return(fallback())
    }
    inc = (ratio[id3]-ratio[id2])/(id3-id2)
    ratio[id2:(n-1)] = ratio[id2-1]+inc*((id2:(n-1))-id2)
    ratio[ratio<0]=0
    right = seq.int(ceiling(n/2),n-1)
    a[right] = ((b^2*x*Nu(sqrt(2*b^2*x)))^2 * ratio)[right]
  }
  neg2 = which(a<0)
  a[neg2] = 0
  if(any(!is.finite(a[l0:l1]))){
    return(fallback())
  }
  integrand = function(s){
    a[s]*(n-s)   ## this is different from single change-point alternative
  }
  result = tryCatch(2*dnorm(b)/b*integrate(integrand, l0, l1, subdivisions=3000, stop.on.error=FALSE)$value,
                    error=function(e) NA_real_)
  if(!is.finite(result) || result < 0 || result > 1){
    return(fallback())
  }
  return(result)
}
  
##Calculate the corrected p-value for Zw
pI_correct_w <- function(n,b,r,x,l0,l1){
  if (b<=0){
    return(1)
  }
  fallback = function(){
    message("Zw: p-value approximation without skewness correction is reported.")
    pWR2(b,n,l0,l1)
  }
  theta_b = rep(NA_real_,n-1)
  pos = which(1+2*r*b>0)
  theta_b[pos] = 2*b/(1+sqrt((1+2*r*b)[pos]))
  ratio = exp((b-theta_b)^2/2 + r*theta_b^3/6)/sqrt(1+r*theta_b)
  a = (b^2*x*Nu(sqrt(2*b^2*x)))^2 * ratio  ## this is different from single change-point
  
  a[which(x==0)] = 0
  nn = n-1-length(pos)
  if (nn>0.75*n){
    message("Extrapolation of Zw for skewness-corrected p-value approximation could not be performed.")
    return(pWR2(b, n, l0, l1))
  }
  if (nn>=(l0-1)+(n-l1)){
    message(paste0('Zw: nn>=(l0-1)+(n-l1), b:',b,'.'))
    neg = which(1+2*r*b<=0)
    if(length(neg) < 2L){
      return(fallback())
    }
    dif = diff(neg)
    id1 = which.max(dif)
    id2 = id1 + ceiling(0.03*n)
    id3 = id2 + ceiling(0.09*n)
    if(id2+1L > length(a) || id3 > length(a) ||
       any(!is.finite(a[c(id2,id2+1L,id3)]))){
      return(fallback())
    }
    inc = (a[id3]-a[id2])/(id3-id2)
    a[id2:1] = a[id2+1]-inc*(1:id2)
    a[(n/2+1):n] = a[(n/2):1]
    neg2 = which(a<0)
    a[neg2] = 0
  }
  if(any(!is.finite(a[l0:l1]))){
    return(fallback())
  }
  integrand = function(s){
    a[s]*(n-s)
  }
  result = tryCatch(dnorm(b)/b*integrate(integrand, l0, l1, subdivisions=3000, stop.on.error=FALSE)$value,
                    error=function(e) NA_real_)
  if(!is.finite(result) || result < 0 || result > 1){
    return(fallback())
  }
  return(result)
}

#The rank function to calculate rank of elements of a matrix.
#Two possible methods: the overall rank and the row-wise rank.
Rise_Rank<-function(S,method=c('overall','row')){
  #S: a N X N similarity matrix
  #R: the output of N X N rank matrix
  method = match.arg(method)
  R = diag(0,nrow(S))
  if(method=='overall'){
    R[upper.tri(R)] = rank(S[upper.tri(S)])
    R = R + t(R)
  }
  if(method=='row'){
    dimnames(R) = dimnames(S)
    for(i in seq_len(nrow(S))){
      R[i,-i] = rank(S[i,-i],ties.method="average")
    }
  }
  return(R)
}

### Calculate the third moment ###
###return a list: list('EZw3','EZdiff3')
ThirdMoment<-function(R,Ebynode){
  n = nrow(R)
  R2 = R*R; R3 = R2*R; Ri = colSums(R); Ri2 = colSums(R2); Rdd = sum(Ri)
  S = C = rep(0,8)
  S[1] = sum(R3); S[2] = sum(Ri2*Ri); 
  
  R.up = matrix(0,n,n)
  R.up[upper.tri(R.up)] = R[upper.tri(R.up)]
  E = which(R.up>0,arr.ind = TRUE)
  E = cbind(E,R[E])
  x0 = 0
  for (i in 1:nrow(E)){
    j = E[i,1]
    k = E[i,2]
    tmp = match(Ebynode[[j]][['node']], Ebynode[[k]][['node']])
    i3 = which(!is.na(tmp))
    x0 = x0 + sum(Ebynode[[j]][['weight']][i3]*Ebynode[[k]][['weight']][tmp[i3]])*E[i,3]
  }
  S[3] = 2*x0
  
  S[4] = sum(R2)*Rdd; S[5] = sum(Ri^3); S[6] = sum(t(R*Ri)*Ri)
  S[7] = sum(Ri^2)*Rdd; S[8] = Rdd^3
  
  C[1] = 4*S[1]; C[2] = 24*(S[2]-S[1]); C[3] = 8*S[3]; 
  C[4] = 6*(S[4]-4*S[2]+2*S[1])
  C[5] = 8*(S[5]-3*S[2]+2*S[1])
  C[6] = 24*(S[6]-S[3]-2*S[2]+S[1])
  C[7] = 12*(S[7]-4*S[6]-2*S[5]-S[4]+2*S[3]+10*S[2]-4*S[1])
  C[8] = S[8] - 12*S[7] + 24*S[6] + 16*S[5] + 6*S[4] - 8*S[3] - 48*S[2] + 16*S[1]    
  
  t = 1:(n-1)
  p1 = t*(t-1)/(n*(n-1))
  p2 = p1*(t-2)/(n-2)
  p3 = p2*(t-3)/(n-3)
  p4 = p3*(t-4)/(n-4)
  p5 = p4*(t-5)/(n-5)
  p6 = p1*(n-t)*(n-t-1)/((n-2)*(n-3))
  p7 = p2*(n-t)*(n-t-1)/((n-3)*(n-4))
  p8 = p3*(n-t)*(n-t-1)/((n-4)*(n-5))
  
  t = (n-1):(1)
  q1 = t*(t-1)/(n*(n-1))
  q2 = q1*(t-2)/(n-2)
  q3 = q2*(t-3)/(n-3)
  q4 = q3*(t-4)/(n-4)
  q5 = q4*(t-5)/(n-5)
  q6 = q1*(n-t)*(n-t-1)/((n-2)*(n-3))
  q7 = q2*(n-t)*(n-t-1)/((n-3)*(n-4))
  q8 = q3*(n-t)*(n-t-1)/((n-4)*(n-5))
  
  #E(U1^3), E(U2^3), E(U1^2 U2), E(U1 U2^2) 
  EU1.3 = p1*C[1] + p2*(C[2]+C[3]) + p3*(C[4]+C[5]+C[6]) + p4*C[7] + p5*C[8]
  EU2.3 = q1*C[1] + q2*(C[2]+C[3]) + q3*(C[4]+C[5]+C[6]) + q4*C[7] + q5*C[8]
  EU1.2U2.1 = p6*C[4]/3 + p7*C[7]/3 + p8*C[8]
  EU1.1U2.2 = q6*C[4]/3 + q7*C[7]/3 + q8*C[8]
  
  t = 1:(n-1)
  pt = (t-1)/(n-2); qt = 1- pt
  EUw3 = qt^3*EU1.3 + 3*qt*pt*(qt*EU1.2U2.1 + pt*EU1.1U2.2) + pt^3*EU2.3
  pt = -1; qt = 1
  EUdiff3 = qt^3*EU1.3 + 3*qt*pt*(qt*EU1.2U2.1 + pt*EU1.1U2.2) + pt^3*EU2.3
  
  r0 = Rdd/(n*(n-1))
  rd = sqrt(sum(R2)/(n*(n-1)))
  r1 = sqrt(mean(Ri^2)/(n-1)^2)
  TI = 1:(n-1)
  muw = n*( TI-1 )*( n - TI-1 )*r0/(n-2)
  mudiff = (n-1)*( 2*TI - n)*r0
  sd_w = sqrt( 2*TI*( TI-1 )*(n - TI)*( n - TI-1 )*(  (n-2)*rd^2 - 2*(n-1)*r1^2 + n*r0^2)/( (n-2)^2*(n-3)) )
  sd_diff = sqrt( 4*(n-1)*TI*(n - TI)*(r1^2 - r0^2) )
  
  EZw3 = (EUw3 - 3*muw*sd_w^2 - muw^3)/sd_w^3
  EZdiff3 = (EUdiff3 - 3*mudiff*sd_diff^2 - mudiff^3)/sd_diff^3 #r
  EZw3[is.na(EZw3)] = 0
  return(list('EZw3'=EZw3,'EZdiff3'=EZdiff3))
}