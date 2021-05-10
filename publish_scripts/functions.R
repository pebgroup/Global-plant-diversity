# Function for resolving cases in which multiple tips in the tree point to the same accepted wcsp name 
# Used in match_data.R
resolve_multiple <- function(MATCHES, wcsp, phylo, phylogb){
  # remove multiple linkages, i.e. when multiple tips in the tree are assigned to the same accepted name, using following criteria
  # 1) preferably keep tips that have molecular data behind them
  # 2) preferably keep tips that have the same genus name as the species they link to in WCSP - important for tips that have been added based on taxonomy
  # 3) randomly thereafter
  for(n in names(table(MATCHES)[table(MATCHES)>1])){
    tips <- which(MATCHES == n)
    erase <- rep(FALSE,length(tips))
    #first criterion: if any of the tips comes from GenBank, erase all that don't
    if(sum(phylo$tip.label[tips] %in% phylogb$tip.label)>0){
      erase <- !phylo$tip.label[tips] %in% phylogb$tip.label
      #break
    } else {#second criterion (only relevant if choosing among non-genbank tips): if at least one of the tips has the same genus as in wcsp, erase the rest
      #retrieve WCSP genus
      gen <- as.vector(wcsp[n,"genus"])
      #retrieve tip genera
      GEN <- vector("character", length(tips))
      for(tip in 1:length(tips)){
        GEN[tip] <- strsplit(phylo$tip.label[tips],"_")[[tip]][1]
      }
      #if at least on e of the tips has the same genus as in wcsp, erase the rest
      if(gen %in% GEN){
        erase <- !GEN == gen
      }
    }
    if(sum(!erase)>1){#if there are still n>1 items to be kept (i.e. not erased), chose n-1 randomly for erasing
      erase[sample(which(erase==FALSE), size=(sum(!erase)-1))] <- TRUE
    } 
    
    #erase tips 
    MATCHES[tips[erase]] <- NA
  }
  return(MATCHES)
}


# A modified version of bind.tip from phytools, in which the "untangling" of the tree is silenced
# This saves time when mulitple tips need to be bound to the tree. However, this means the tree neeeds to be 
# "untangled" manually after adding the last tip using untangle()
# Used in add_species.R
# Requires phytools
bind.tip2 <- function (tree, tip.label, edge.length = NULL, where = NULL, 
                       position = 0, interactive = FALSE, ...) 
{
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  use.edge.length <- if (is.null(tree$edge.length)) 
    FALSE
  else TRUE
  if (use.edge.length == FALSE) 
    tree <- compute.brlen(tree)
  if (interactive == TRUE) {
    plotTree(tree, ...)
    cat(paste("Click where you would like to bind the tip \"", 
              tip.label, "\"\n", sep = ""))
    flush.console()
    obj <- get.treepos(message = FALSE)
    where <- obj$where
    position <- obj$pos
  }
  else if (is.null(where)) 
    where <- Ntip(tree) + 1
  if (where <= Ntip(tree) && position == 0) {
    pp <- 1e-12
    if (tree$edge.length[which(tree$edge[, 2] == where)] <= 
        1e-12) {
      tree$edge.length[which(tree$edge[, 2] == where)] <- 2e-12
      ff <- TRUE
    }
    else ff <- FALSE
  }
  else pp <- position
  if (is.null(edge.length) && is.ultrametric(tree)) {
    H <- nodeHeights(tree)
    if (where == (Ntip(tree) + 1)) 
      edge.length <- max(H)
    else edge.length <- max(H) - H[tree$edge[, 2] == where, 
                                   2] + position
  }
  tip <- list(edge = matrix(c(2, 1), 1, 2), tip.label = tip.label, 
              edge.length = edge.length, Nnode = 1)
  class(tip) <- "phylo"
  obj <- bind.tree(tree, tip, where = where, position = pp)
  if (where <= Ntip(tree) && position == 0) {
    nn <- obj$edge[which(obj$edge[, 2] == which(obj$tip.label == 
                                                  tip$tip.label)), 1]
    obj$edge.length[which(obj$edge[, 2] == nn)] <- obj$edge.length[which(obj$edge[, 
                                                                                  2] == nn)] + 1e-12
    obj$edge.length[which(obj$edge[, 2] == which(obj$tip.label == 
                                                   tip$tip.label))] <- 0
    obj$edge.length[which(obj$edge[, 2] == which(obj$tip.label == 
                                                   tree$tip.label[where]))] <- 0
  }
  #root.time <- if (!is.null(obj$root.time)) 
  #  obj$root.time
  #else NULL
  #obj <- untangle(obj, "read.tree")
  #if (!is.null(root.time)) 
  #  obj$root.time <- root.time
  #if (interactive) 
  #  plotTree(obj, ...)
  if (!use.edge.length) 
    obj$edge.length <- NULL
  obj
}


# function that estimates the root distance (RD) for all tips in a phylo object
# polytomies are accounted for assuming a birth-death process and using an empirically estimated function (-2/3 + 2*log(N))
# requires ape
# used in phylostruct.R
root.distance <- function(tree, mc.cores=8){
  # score polytomies in tree and estimate how many nodes need to be added for each one
  polys <- table(tree$edge[,1]) #calculates for each edge how many descendents it has
  polys <- -1.677857 + log(polys) * 2.00427 #transforms this into expected RD
  polys[polys <0] = 0
  # calculate root distances for all tips in parallel
  RDs <- mclapply(1:Ntip(tree), countnodes, tree=tree, polys=polys, mc.cores = mc.cores)
  return(RDs)
}

# function that modifies corrplot from ggcorrplot package
# to allow for visual highlighting of values above a certain thereshold
# adjusting corrplot function to highlight correlations higher lower than threshold
my_corrplot <- function (corr, method = c("square", "circle"), type = c("full", "lower", "upper"), 
                         ggtheme = ggplot2::theme_minimal, title = "", 
                         show.legend = TRUE, legend.title = "rho", show.diag = FALSE, 
                         colors = c("blue", "white", "red"), outline.color = "gray", 
                         hc.order = FALSE, hc.method = "complete", lab = FALSE, lab_col = "black", 
                         lab_size = 4, p.mat = NULL, sig.level = 0.05, insig = c("pch", "blank"), 
                         pch = 4, pch.col = "black", pch.cex = 5, tl.cex = 12, 
                         tl.col = "black", tl.srt = 45, digits = 2, highlight=TRUE, highthreshold=0.7) 
{
  type <- match.arg(type)
  method <- match.arg(method)
  insig <- match.arg(insig)
  if (inherits(corr, "cor_mat")) {
    cor.mat <- corr
    corr <- .tibble_to_matrix(cor.mat)
    p.mat <- .tibble_to_matrix(attr(cor.mat, "pvalue"))
  }
  if (!is.matrix(corr) & !is.data.frame(corr)) {
    stop("Need a matrix or data frame!")
  }
  corr <- as.matrix(corr)
  corr <- base::round(x = corr, digits = digits)
  if (hc.order) {
    ord <- .hc_cormat_order(corr)
    corr <- corr[ord, ord]
    if (!is.null(p.mat)) {
      p.mat <- p.mat[ord, ord]
      p.mat <- base::round(x = p.mat, digits = digits)
    }
  }
  if (type == "lower") {
    corr <- .get_lower_tri(corr, show.diag)
    p.mat <- .get_lower_tri(p.mat, show.diag)
  }
  else if (type == "upper") {
    corr <- .get_upper_tri(corr, show.diag)
    p.mat <- .get_upper_tri(p.mat, show.diag)
  }
  corr <- reshape2::melt(corr, na.rm = TRUE)
  colnames(corr) <- c("Var1", "Var2", "value")
  corr$pvalue <- rep(NA, nrow(corr))
  corr$signif <- rep(NA, nrow(corr))
  if (!is.null(p.mat)) {
    p.mat <- reshape2::melt(p.mat, na.rm = TRUE)
    corr$coef <- corr$value
    corr$pvalue <- p.mat$value
    corr$signif <- as.numeric(p.mat$value <= sig.level)
    p.mat <- subset(p.mat, p.mat$value > sig.level)
    if (insig == "blank") {
      corr$value <- corr$value * corr$signif
    }
  }
  corr$abs_corr <- abs(corr$value) * 10
  p <- ggplot2::ggplot(data = corr, mapping = ggplot2::aes_string(x = "Var1", 
                                                                  y = "Var2", fill = "value"))
  if (method == "square") {
    p <- p + ggplot2::geom_tile(color = outline.color)
  }
  else if (method == "circle") {
    p <- p + ggplot2::geom_point(color = outline.color, shape = 21, 
                                 ggplot2::aes_string(size = "abs_corr")) + ggplot2::scale_size(range = c(4, 
                                                                                                         10)) + ggplot2::guides(size = FALSE)
  }
  p <- p + ggplot2::scale_fill_gradient2(low = colors[1], high = colors[3], 
                                         mid = colors[2], midpoint = 0, limit = c(-1, 1), space = "Lab", 
                                         name = legend.title)
  if (class(ggtheme)[[1]] == "function") {
    p <- p + ggtheme()
  }
  else if (class(ggtheme)[[1]] == "theme") {
    p <- p + ggtheme
  }
  p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = tl.srt, 
                                                              vjust = 1, size = tl.cex, hjust = 1), axis.text.y = ggplot2::element_text(size = tl.cex)) + 
    ggplot2::coord_fixed()
  label <- round(x = corr[, "value"], digits = digits)
  if (!is.null(p.mat) & insig == "blank") {
    ns <- corr$pvalue > sig.level
    if (sum(ns) > 0) 
      label[ns] <- " "
  }
  if (lab) {
    if(highlight==TRUE) {
      nh <- abs(corr$value) >= highthreshold
      typo <- c(1, 2)[as.numeric(nh)+1]
      #   p <- p + ggplot2::aes_string(face = typo)
      p <- p + ggplot2::geom_text(mapping = ggplot2::aes_string(x = "Var1", 
                                                                y = "Var2", fontface=typo), label = label, color = lab_col, size = lab_size)
    }else
      p <- p + ggplot2::geom_text(mapping = ggplot2::aes_string(x = "Var1", 
                                                                y = "Var2"), label = label, color = lab_col, size = lab_size)
  }
  
  if (!is.null(p.mat) & insig == "pch") {
    p <- p + ggplot2::geom_point(data = p.mat, mapping = ggplot2::aes_string(x = "Var1", 
                                                                             y = "Var2"), shape = pch, size = pch.cex, color = pch.col)
  }
  if (title != "") {
    p <- p + ggplot2::ggtitle(title)
  }
  if (!show.legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  p <- p + .no_panel()
  p
}
environment(my_corrplot) <- asNamespace('ggcorrplot')



# function for z-transformation to avoid variance issues
ztrans <- function(x){(x - mean(x)) / sd(x)} 


# scale data between values of 0 and 1
range01 <- function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))}




# function to relabel edges in SEM plots to include thresholds of significance
sem_sig_labels <- function(sem_fit_object){
  table2<-parameterEstimates(sem_fit_object,standardized=TRUE)[!is.na(parameterEstimates(sem_fit_object)$pvalue) & 
                                                                 parameterEstimates(sem_fit_object)$op!=":=",]
  sig <- rep(" ",nrow(table2))
  sig[table2$pvalue<=0.05] <- "*"
  sig[table2$pvalue<=0.01] <- "**"
  sig[table2$pvalue<=0.001] <- "***"
  b <- paste0(round(table2$std.all,2), sig)
  return(b)
}



##### lavaan code ##############################################
# Code by Dr. Jarrett E.K. Byrnes; Last Modified 12/06/2016
# created in regards to Issue #44 on lavaan
# https://github.com/yrosseel/lavaan/issues/44
# added personal modification to 

my.predict_lavaan <- function(fit, newdata = NULL){
  # test this, uncomment when done
  # fit <-  c2.fit.p
  # newdata = NULL
  
  stopifnot(inherits(fit, "lavaan"))
  
  #Make sure we can use this
  if(!inspect(fit, "meanstructure")) stop("Need to supply meanstructure = TRUE in fit\n")
  if(is.null(newdata)){
    newdata <- data.frame(inspect(fit, "data"))
    names(newdata) <- lavNames(fit)
  }
  
  if(length(lavNames(fit, type="lv"))!=0) stop("Does not currently work with latent variables\n")
  
  #check for new data
  if(sum(!(lavNames(fit, type="ov.x") %in% names(newdata)))>0) stop("Not all exogenous variables supplied!")
  
  #Add some new columns to newdata
  newdata$Intercept <- 1
  
  #newdata[lavNames(fit, "ov.nox")] <- 0 # this sets non-exogenous observed variables to zero 
  
  
  
  mod_df <- data.frame(lhs = fit@ParTable$lhs,
                       op = fit@ParTable$op,
                       rhs = fit@ParTable$rhs,
                       exo = fit@ParTable$exo,
                       est = fit@ParTable$est,
                       se = fit@ParTable$se,
                       stringsAsFactors=FALSE)
  
  #Drop covariances
  mod_df <- mod_df[-which(mod_df$op=="~~"),]
  mod_df[which(mod_df$op=="~1"),]$rhs <- "Intercept"
  
  #get rid of exogenous on lhs
  mod_df <- mod_df[-which(mod_df$exo==1),]
  
  #Order by lhs
  mod_df <- mod_df[sort(mod_df$lhs, index.return=TRUE)$ix,]
  
  #let us know which variables on the rhs are exogenous
  mod_df$ord <- 0
  mod_df[which(!(mod_df$rhs %in% mod_df$lhs)),]$ord <- 1
  
  #Make a "order"
  ord_current <- 1
  while(sum(mod_df$ord==0)>0){
    for(r in unique(mod_df$lhs)){
      val <-  sum(mod_df[which(mod_df$lhs==r),]$ord==0)
      if(val==0) {
        mod_df[which(mod_df$lhs==r),]$ord <- ord_current
        
        if(sum(mod_df$rhs==r)>0)
          mod_df[which(mod_df$rhs==r),]$ord <- ord_current+1
      }
    }
    ord_current <- ord_current +1
  }
  
  #correct for ragged ordering
  for(r in unique(mod_df$lhs)){
    mod_df[which(mod_df$lhs==r),]$ord <- max(mod_df[which(mod_df$lhs==r),]$ord)
  }
  
  #sort by order 
  mod_df <- mod_df[sort(mod_df$ord, index.return=TRUE)$ix,]
  
  #now do the fitting in order
  fit_df <- data.frame(base = rep(1, nrow(newdata)))
  
  for(r in unique(mod_df$lhs)){
    #r="sr_trans"
    subdf <- subset(mod_df, mod_df$lhs==r)
    #make a formula
    rhs <- paste0(subdf$rhs, collapse=" + ")
    form <- as.formula(paste0(r, " ~ ", rhs))
    
    #use formula to get right part of the data in right format
    mod_mat <- model.matrix(form, newdata)[,-1] 
    new_val = mod_mat %*% subdf$est # this is the matrix multiplication
    fit_df[[r]] <- new_val
    
    # TEST
    #plot(fit_df[["sr_trans"]],lm1$fitted.values) # works with commented line, does not work with uncommented line
    
    
    #newdata[[r]] <- new_val # this replaces the base data while working with it.... nope
  }
  
  return(fit_df[,-1]) # minus 1 removes the base column
  
}

my.fitted_lavaan <- function(fit){
  my.predict_lavaan(fit)
}

my.residuals_lavaan <- function(fit){
  fitted_vals <- my.fitted_lavaan(fit)
  
  rawdata <- data.frame(inspect(fit, "data"))
  names(rawdata) <- lavNames(fit)
  
  res <- data.frame(base = rep(1, nrow(rawdata)))
  for(vals in names(fitted_vals)){
    res[[vals]] <- rawdata[[vals]] - fitted_vals[[vals]] 
  }
  
  return(res[,-1])
}


lavSpatialCorrect_correct <- function(obj, xvar, yvar, alpha=0.05){
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



