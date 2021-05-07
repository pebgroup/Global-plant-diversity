lavSpatialCorrect <- function(obj, xvar, yvar, alpha=0.05){
  require(lavaan)
  require(ape)
  
  #first, get the residuals from the model
  if(length(lavNames(obj, type="lv"))!=0){
    resids <- as.data.frame(residuals(obj, "casewise"))
  }else{
    resids <- as.data.frame(residuals_lavaan(obj))
  }
  
  #get only endogenous variables
  resids <- resids[,which(apply(resids, 2, function(x) length(unique(x))) !=1)]
  

  #make a distance matrix
  distMat <- as.matrix(dist(cbind(xvar, yvar)))
  
  #invert this matrix for weights
  distsInv <- 1/distMat
  diag(distsInv) <- 0
  
  morans_i <- lapply(resids,  function(x){
    mi <- Moran.I(c(x), distsInv)
    if(mi$p.value>alpha){
      mi$n.eff <- nrow(resids) #don't correct sample size
    }else{
      #large sample size approximation
      mi$n.eff <- nrow(resids)*(1-mi$observed)/(1+mi$observed)
    }
    
    #return the list
    mi
    
    })
  
  
  #get the vcov matrix
  v <- diag(vcov(obj))
  n <- nrow(resids)
  #using new sample sizes, for each variable, calculate new Z-scores
  params <- lapply(names(morans_i), function(acol){
    idx <- grep(paste0(acol, "~"),names(v))  #regression or covariances
    idx <- c(idx, grep(paste0("=~",acol),names(v)))  #latent variable definitions
    v_idx <- v[idx]*n/morans_i[[acol]]$n.eff
    ret <-  data.frame(Parameter = names(v)[idx], Estimate=coef(obj)[idx], 
                       n.eff = morans_i[[acol]]$n.eff, Std.err = sqrt(v_idx))
    ret[["Z-value"]] <- ret$Estimate/ret$Std.err
    ret[["P(>|z|)"]] <- 2*pnorm(abs(ret[["Z-value"]]), lower.tail=F)
  
    ret
    
  })
  names(params) <- names(morans_i)

  mi <- lapply(morans_i, function(m) {
         data.frame(observed=m$observed, expected=m$expected, 
                    sd=m$sd, p.value = m$p.value, n.eff = m$n.eff)
  })

 return(list(Morans_I = mi, parameters=params)) 
  
}

lavSpatialCorrect_correct <- function(obj, xvar, yvar, alpha=0.05){
  require(lavaan)
  require(ape)
  
  # # test data
  obj= fm.p
  xvar=dat_no.na$x
  yvar=dat_no.na$y
  alpha=0.05
  
  #first, get the residuals from the model
  if(length(lavNames(obj, type="lv"))!=0){
    resids <- as.data.frame(residuals(obj, "casewise"))
  }else{
    resids <- as.data.frame(my.residuals_lavaan(obj))
  }
  
  #get only endogenous variables
  resids <- resids[,which(apply(resids, 2, function(x) length(unique(x))) !=1)]
  
  
  #make a distance matrix
  distMat <- as.matrix(dist(cbind(xvar, yvar)))
  
  #invert this matrix for weights
  distsInv <- 1/distMat
  diag(distsInv) <- 0
  
  morans_i <- lapply(resids,  function(x){
    mi <- Moran.I(c(x), distsInv)
    if(mi$p.value>alpha){
      mi$n.eff <- nrow(resids) #don't correct sample size
    }else{
      #large sample size approximation
      mi$n.eff <- nrow(resids)*(1-mi$observed)/(1+mi$observed)
    }
    
    #return the list
    mi
    
  })
  
  
  #get the vcov matrix
  v <- diag(vcov(obj))
  n <- nrow(resids)
  #using new sample sizes, for each variable, calculate new Z-scores
  params <- lapply(names(morans_i), function(acol){
    idx <- grep(paste0(acol, "~"),names(v))  #regression or covariances
    idx <- c(idx, grep(paste0("=~",acol),names(v)))  #latent variable definitions
    v_idx <- v[idx]*n/morans_i[[acol]]$n.eff
    ret <-  data.frame(Parameter = names(v)[idx], Estimate=coef(obj)[idx], 
                       n.eff = morans_i[[acol]]$n.eff, Std.err = sqrt(v_idx))
    ret[["Z-value"]] <- ret$Estimate/ret$Std.err
    ret[["P(>|z|)"]] <- 2*pnorm(abs(ret[["Z-value"]]), lower.tail=F)
    
    ret
    
  })
  names(params) <- names(morans_i)
  
  mi <- lapply(morans_i, function(m) {
    data.frame(observed=m$observed, expected=m$expected, 
               sd=m$sd, p.value = m$p.value, n.eff = m$n.eff)
  })
  
  return(list(Morans_I = mi, parameters=params)) 
  
}

# lavSpatialCorrect_correct.std <- function(obj, xvar, yvar, alpha=0.05){
#   require(lavaan)
#   require(ape)
#   
#   #first, get the residuals from the model
#   if(length(lavNames(obj, type="lv"))!=0){
#     resids <- as.data.frame(residuals(obj, "casewise"))
#   }else{
#     resids <- as.data.frame(my.residuals_lavaan(obj))
#   }
#   
#   #get only endogenous variables
#   resids <- resids[,which(apply(resids, 2, function(x) length(unique(x))) !=1)]
#   
#   
#   #make a distance matrix
#   distMat <- as.matrix(dist(cbind(xvar, yvar)))
#   
#   #invert this matrix for weights
#   distsInv <- 1/distMat
#   diag(distsInv) <- 0
#   
#   morans_i <- lapply(resids,  function(x){
#     mi <- Moran.I(c(x), distsInv)
#     if(mi$p.value>alpha){
#       mi$n.eff <- nrow(resids) #don't correct sample size
#     }else{
#       #large sample size approximation
#       mi$n.eff <- nrow(resids)*(1-mi$observed)/(1+mi$observed)
#     }
#     
#     #return the list
#     mi
#     
#   })
#   
#   
#   #get the vcov matrix
#   v <- diag(vcov(obj))
#   n <- nrow(resids)
#   #using new sample sizes, for each variable, calculate new Z-scores
#   params <- lapply(names(morans_i), function(acol){
#     idx <- grep(paste0(acol, "~"),names(v))  #regression or covariances
#     idx <- c(idx, grep(paste0("=~",acol),names(v)))  #latent variable definitions
#     
#     v_idx <- v[idx]*n/morans_i[[acol]]$n.eff
#     stand_solution <- standardizedSolution(obj)
#     ret <-  data.frame(Parameter = names(v)[idx], 
#                        Estimate=coef(obj)[idx], 
#                        n.eff = morans_i[[acol]]$n.eff,
#                        Std.err = sqrt(v_idx)),
#                       
#     ret[["Z-value"]] <- ret$Estimate/ret$Std.err
#     ret[["P(>|z|)"]] <- 2*pnorm(abs(ret[["Z-value"]]), lower.tail=F)
#     
#     ret
#     
#   })
#   names(params) <- names(morans_i)
#   
#   mi <- lapply(morans_i, function(m) {
#     data.frame(observed=m$observed, expected=m$expected, 
#                sd=m$sd, p.value = m$p.value, n.eff = m$n.eff)
#   })
#   
#   return(list(Morans_I = mi, parameters=params)) 
#   
# }

#lavSpatialCorrect(obj, xvar, yvar)

lavSpatialCorrect_moran.mc <- function(obj, xvar, yvar, alpha=0.05){
  require(lavaan)
  require(ape)

  #  TEST DATA
  obj=c2.fit.p
  xvar=dat_no.na$x
  yvar=dat_no.na$y
  distance=1600
  alpha=0.05
  distance=c(500,2000)
  
  #first, get the residuals from the model
  if(length(lavNames(obj, type="lv"))!=0){
    resids <- as.data.frame(residuals(obj, "casewise"))
  }else{
    resids <- as.data.frame(my.residuals_lavaan(obj))
  }
  
  #get only endogenous variables
  resids <- resids[,which(apply(resids, 2, function(x) length(unique(x))) !=1)]
  
  
  #make a distance matrix
  distMat <- as.matrix(dist(cbind(xvar, yvar)))
  
  #invert this matrix for weights
  distsInv <- 1/distMat
  diag(distsInv) <- 0
  
  morans_i <- lapply(resids,  function(x){
    mi <- moran.mc(c(x), mat2listw(distsInv, style = "B"), nsim=599, zero.policy=T)
    #mi <- Moran.I(c(x), distsInv)
    if(mi$p.value>alpha){
      mi$n.eff <- nrow(resids) #don't correct sample size
    }else{
      #large sample size approximation
      mi$n.eff <- nrow(resids)*(1-mi$statistic)/(1+mi$statistic)
    }
    
    #return the list
    mi
    
  })
  
  #get the vcov matrix
  v <- diag(vcov(obj))
  n <- nrow(resids)
  #using new sample sizes, for each variable, calculate new Z-scores
  params <- lapply(names(morans_i), function(acol){
    idx <- grep(paste0(acol, "~"),names(v))  #regression or covariances
    idx <- c(idx, grep(paste0("=~",acol),names(v)))  #latent variable definitions
    v_idx <- v[idx]*n/morans_i[[acol]]$n.eff
    ret <-  data.frame(Parameter = names(v)[idx], Estimate=coef(obj)[idx], 
                       n.eff = morans_i[[acol]]$n.eff, Std.err = sqrt(v_idx))
    ret[["Z-value"]] <- ret$Estimate/ret$Std.err
    ret[["P(>|z|)"]] <- 2*pnorm(abs(ret[["Z-value"]]), lower.tail=F)
    
    ret
    
  })
  names(params) <- names(morans_i)
  
  mi <- lapply(morans_i, function(m) {
    data.frame(observed=m$observed, expected=m$expected, 
               sd=m$sd, p.value = m$p.value, n.eff = m$n.eff)
  })
  
  return(list(Morans_I = mi, parameters=params)) 
  
}
# does not work yet





my.lavSpatialCorrect <- function(obj, xvar, yvar, distance, alpha=0.05){
  # include changes to calc morans I on the basis of distance, not globally.
  # change to use moran.mc for distance and also Monte Carlo p value estimation
  
#  TEST DATA
  obj=c2.fit.p
  xvar=dat_no.na$x
  yvar=dat_no.na$y
  distance=1600
  alpha=0.05
  distance=c(500,2000)
  
  require(lavaan)
  require(ape)
  
  #first, get the residuals from the model
  if(length(lavNames(obj, type="lv"))!=0){
    resids <- as.data.frame(residuals(obj, "casewise"))
  }else{
    resids <- as.data.frame(my.residuals_lavaan(obj))
  }
  
  #get only endogenous variables
  resids <- resids[,which(apply(resids, 2, function(x) length(unique(x))) !=1)]
  
  
  # #make a distance matrix
  # distMat <- as.matrix(dist(cbind(xvar, yvar)))
  # 
  # #invert this matrix for weights
  # distsInv <- 1/distMat
  # diag(distsInv) <- 0
  # distsInv[is.infinite(distsInv)] <- 0 # fix infinite cases
  
  # get neighbors
  coo <- cbind(yvar, xvar)
  S.dist  <-  dnearneigh(coo, 0, distance, longlat = TRUE)
  lw <- nb2listw(S.dist, style="W", zero.policy=T) 
  # mi <- moran.mc(resids[,4], lw, nsim=599, zero.policy=T)   
  # mi$p.value
  # mi$observed
  # #plot(dat_no.na$sr_trans, resids[,4])
    
  morans_i <- lapply(resids,  function(x){
    mi <- moran.mc(c(x), lw, nsim=599, zero.policy=T)
    #mi <- Moran.I(c(x), distsInv)
    if(mi$p.value>alpha){
      mi$n.eff <- nrow(resids) #don't correct sample size
    }else{
      #large sample size approximation
      mi$n.eff <- nrow(resids)*(1-mi$statistic)/(1+mi$statistic)
    }
    
    #return the list
    mi
    
  })
  
  
  #get the vcov matrix
  v <- diag(vcov(obj))
  n <- nrow(resids)
  #using new sample sizes, for each variable, calculate new Z-scores
  params <- lapply(names(morans_i), function(acol){
    idx <- grep(paste0(acol, "~"),names(v))  #regression or covariances
    idx <- c(idx, grep(paste0("=~",acol),names(v)))  #latent variable definitions
    v_idx <- v[idx]*n/morans_i[[acol]]$n.eff
    ret <-  data.frame(Parameter = names(v)[idx], Estimate=coef(obj)[idx], 
                       n.eff = morans_i[[acol]]$n.eff, Std.err = sqrt(v_idx))
    ret[["Z-value"]] <- ret$Estimate/ret$Std.err
    ret[["P(>|z|)"]] <- 2*pnorm(abs(ret[["Z-value"]]), lower.tail=F)
    
    ret
    
  })
  names(params) <- names(morans_i)
  
  mi <- lapply(morans_i, function(m) {
    data.frame(observed=m$statistic, expected=mean(m$res), 
               sd=sd(m$res), p.value = m$p.value, n.eff = m$n.eff)
  })
  
  return(list(Morans_I = mi, parameters=params)) 
  
}

#lavSpatialCorrect(obj, xvar, yvar)
